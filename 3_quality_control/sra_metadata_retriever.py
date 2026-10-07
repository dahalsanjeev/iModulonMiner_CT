#!/usr/bin/env python3
"""
sra_metadata_retriever.py
=========================
Retrieves biological growth-condition metadata for a list of SRA experiment
accession IDs (SRX / ERX / DRX) from NCBI SRA and BioSample databases.

Target fields
-------------
  strain, media, carbon_source, nitrogen_source, temperature, pH, OD,
  growth_phase, growth_condition, culture_type, aeration, agitation,
  time_point, organism, genotype

Pipeline
--------
  1. esearch  → SRA UID
  2. efetch   → SRA experiment XML  (title, design, sample attributes)
  3. elink    → linked BioSample ID
  4. efetch   → BioSample XML       (all structured attributes)
  5. Synonym mapping → normalised target columns
  6. Regex extraction → attempt to pull values from free-text description
  7. Completeness score + missing-field flag for manual review
  8. Progressive CSV save

Dependencies
------------
  pip install requests pandas

Optional (for bulk project fetches):
  pip install pysradb

Usage
-----
  python sra_metadata_retriever.py [METADATA_FILE]   # use the metadata file from interim that passed QC1 steps
"""

import re
import sys
import time
import logging
import argparse
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Optional

import requests
import pandas as pd

# ── Logging ───────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────────────────
#  CONFIGURATION  (edit here)
# ─────────────────────────────────────────────────────────────────────────────

# Your e-mail address (REQUIRED by NCBI policy)
NCBI_EMAIL = [EMAIL_ADDRESS]

# Optional NCBI API key → raises rate limit from 3 to 10 req/s
# Register free at: https://www.ncbi.nlm.nih.gov/account/
NCBI_API_KEY = [API_KEY]

# Seconds between requests (0.34 without API key, 0.11 with key)
REQUEST_DELAY = 0.50

# Output CSV path
OUTPUT_CSV = [OUTPUT_FILE]

# ─────────────────────────────────────────────────────────────────────────────
#  FIELD SYNONYMS
#  Maps a canonical target field → list of known BioSample attribute names
#  (case-insensitive; spaces/dashes/underscores treated as equivalent)
# ─────────────────────────────────────────────────────────────────────────────

FIELD_SYNONYMS: dict = {
    "strain": [
        "strain", "strain_clone", "isolate", "host_strain", "bacterial_strain",
        "strain_background", "genotype_strain",
    ],
    "media": [
        "media", "medium", "growth_medium", "culture_medium", "growth_media",
        "cultivation_medium", "defined_medium", "fermentation_medium",
    ],
    "carbon_source": [
        "carbon_source", "c_source", "carbon", "carb_source", "sole_carbon_source",
        "primary_carbon_source",
    ],
    "nitrogen_source": [
        "nitrogen_source", "n_source", "nitrogen", "nitro_source",
        "sole_nitrogen_source",
    ],
    "temperature": [
        "temperature", "temp", "growth_temperature", "incubation_temperature",
        "culture_temperature", "growth_temp", "fermentation_temperature",
    ],
    "pH": [
        "ph", "growth_ph", "culture_ph", "initial_ph", "starting_ph",
        "ph_value", "medium_ph",
    ],
    "OD": [
        "od", "od600", "od_600", "optical_density", "od_at_harvest",
        "starting_od", "od_at_inoculation", "final_od",
    ],
    "growth_phase": [
        "growth_phase", "culture_phase", "physiological_state", "cell_phase",
        "growth_stage", "culture_stage", "harvest_phase",
    ],
    "growth_condition": [
        "growth_condition", "condition", "growth_conditions", "culture_condition",
        "experimental_condition", "treatment", "stress_condition",
        "perturbation", "experimental_treatment",
    ],
    "culture_type": [
        "culture_type", "fermentation_type", "bioreactor_type", "cultivation_type",
        "culture_format", "flask_type", "vessel_type",
    ],
    "aeration": [
        "aeration", "oxygen", "dissolved_oxygen", "do", "aerobic_anaerobic",
        "oxygen_tension", "oxygen_level", "anaerobic", "aerobic",
    ],
    "agitation": [
        "agitation", "stirring", "rpm", "shaking_speed", "agitation_speed",
        "mixing_speed", "rotation",
    ],
    "time_point": [
        "time_point", "sampling_time", "harvest_time", "induction_time",
        "time_of_sampling", "sample_time", "collection_time_point",
    ],
    "organism": [
        "organism", "species", "scientific_name", "organism_name",
    ],
    "genotype": [
        "genotype", "genetic_background", "genotype_variation", "mutation",
        "gene_knockout", "deletion",
    ],
    "replicate": [
        "replicate", "biological_replicate", "technical_replicate",
        "rep", "sample_number",
    ],
}

# ─────────────────────────────────────────────────────────────────────────────
#  REGEX PATTERNS for free-text extraction
#  Applied to: experiment_title, design_description, study_abstract
# ─────────────────────────────────────────────────────────────────────────────

FREE_TEXT_PATTERNS: dict = {
    "temperature": [
        r"(\d{1,3})\s*°C",
        r"(\d{1,3})\s*degrees?\s*C(?:elsius)?",
        r"temp(?:erature)?[:\s=]+(\d{1,3})\s*°?C",
    ],
    "pH": [
        r"pH\s*[=:~]?\s*(\d+\.?\d*)",
        r"pH\s+of\s+(\d+\.?\d*)",
    ],
    "OD": [
        r"OD(?:600)?\s*[=:~]?\s*(\d+\.?\d*)",
        r"optical\s+density\s*(?:of|=|:)?\s*(\d+\.?\d*)",
    ],
    "growth_phase": [
        (r"(exponential|log|stationary|lag|mid.?log|early.?log"
         r"|late.?log|mid.?exponential|early.?stationary)\s+phase"),
    ],
    "carbon_source": [
        (r"(?:grown|cultured)\s+(?:on|with|in)\s+"
         r"(glucose|glycerol|acetate|succinate|lactate|xylose|arabinose"
         r"|fructose|maltose|mannose|ethanol|pyruvate|citrate|galactose"
         r"|sorbitol|cellobiose)"),
        r"carbon\s+source[:\s]+([A-Za-z0-9\-]+)",
    ],
    "media": [
        (r"\b(LB|M9|MOPS|RPMI|DMEM|BHI|SOC|TB|2xYT|YPD|YNB|EMM"
         r"|rich\s+medium|minimal\s+medium|defined\s+medium"
         r"|complex\s+medium)\b"),
    ],
    "aeration": [
        r"(aerobic|anaerobic|microaerophilic|anoxic|oxic|semi.?aerobic)",
    ],
    "strain": [
        (r"\b(WT|wild.?type|MG1655|BW25113|DH5a|BL21|W3110|K-12|REL606"
         r"|BY4741|BY4742|W303|CEN\.PK|S288c)\b"),
        r"strain\s+([A-Za-z0-9_\-]+)",
    ],
}

# ─────────────────────────────────────────────────────────────────────────────
#  NCBI E-utilities helpers
# ─────────────────────────────────────────────────────────────────────────────
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _ncbi_params(extra: dict) -> dict:
    p = {"email": NCBI_EMAIL}
    if NCBI_API_KEY:
        p["api_key"] = NCBI_API_KEY
    p.update(extra)
    return p


def ncbi_get(endpoint: str, params: dict, retries: int = 3) -> requests.Response:
    """GET an NCBI E-utilities endpoint with retry + rate-limit backoff."""
    url = f"{NCBI_BASE}/{endpoint}"
    for attempt in range(retries):
        try:
            resp = requests.get(url, params=_ncbi_params(params), timeout=30)
            resp.raise_for_status()
            time.sleep(REQUEST_DELAY)
            return resp
        except requests.exceptions.RequestException as exc:
            wait = 2 ** attempt
            log.warning(f"    Attempt {attempt+1}/{retries} failed ({exc}); "
                        f"retrying in {wait}s ...")
            time.sleep(wait)
    raise RuntimeError(f"NCBI request failed after {retries} retries: {url}")


def accession_to_sra_uid(accession: str) -> Optional[str]:
    """Convert an SRA experiment accession to its internal NCBI UID."""
    resp = ncbi_get("esearch.fcgi", {
        "db": "sra", "term": accession, "retmax": 1, "retmode": "xml",
    })
    root = ET.fromstring(resp.text)
    ids = root.findall(".//Id")
    return ids[0].text if ids else None


def fetch_sra_xml(uid: str) -> Optional[ET.Element]:
    """Fetch full SRA experiment XML for a given NCBI UID."""
    resp = ncbi_get("efetch.fcgi", {
        "db": "sra", "id": uid, "rettype": "full", "retmode": "xml",
    })
    try:
        return ET.fromstring(resp.text)
    except ET.ParseError:
        log.error("    Could not parse SRA XML")
        return None


def sra_uid_to_biosample_id(sra_uid: str) -> Optional[str]:
    """Use elink to find the BioSample linked to an SRA UID."""
    resp = ncbi_get("elink.fcgi", {
        "dbfrom": "sra", "db": "biosample", "id": sra_uid, "retmode": "xml",
    })
    root = ET.fromstring(resp.text)
    ids = root.findall(".//LinkSetDb/Link/Id")
    return ids[0].text if ids else None


def fetch_biosample_xml(biosample_id: str) -> Optional[ET.Element]:
    """Fetch BioSample XML by internal BioSample UID."""
    resp = ncbi_get("efetch.fcgi", {
        "db": "biosample", "id": biosample_id,
        "rettype": "xml", "retmode": "xml",
    })
    try:
        root = ET.fromstring(resp.text)
        # efetch wraps result in <BioSampleSet>; return the inner <BioSample>
        bs = root.find(".//BioSample")
        return bs if bs is not None else root
    except ET.ParseError:
        log.error("    Could not parse BioSample XML")
        return None



# ─────────────────────────────────────────────────────────────────────────────
#  PUBLICATION RETRIEVAL
#
#  Three complementary strategies (tried in order, results merged):
#
#  1. Direct elink:  sra  ──► pubmed
#     Works when submitters explicitly linked a paper at submission time.
#
#  2. BioProject path:  sra ──► bioproject ──► pubmed
#     BioProject records often carry publication links even when the SRA
#     experiment entry does not.
#
#  3. PubMed text search:  esearch db=pubmed term="SRX123456"
#     Fallback – catches papers that cite the accession in their abstract/
#     methods but were never formally linked in NCBI.
#
#  For all PMIDs found, efetch returns structured PubMed XML from which
#  title, first-author, full author list, journal, year, volume/pages,
#  DOI and abstract are extracted.
# ─────────────────────────────────────────────────────────────────────────────

def sra_uid_to_bioproject_uid(sra_uid: str) -> Optional[str]:
    """elink: sra → bioproject.  Returns BioProject NCBI UID (not accession)."""
    resp = ncbi_get("elink.fcgi", {
        "dbfrom": "sra", "db": "bioproject", "id": sra_uid, "retmode": "xml",
    })
    root = ET.fromstring(resp.text)
    ids = root.findall(".//LinkSetDb/Link/Id")
    return ids[0].text if ids else None


def bioproject_uid_to_pmids(bp_uid: str) -> list:
    """elink: bioproject → pubmed.  Returns list of PMID strings."""
    resp = ncbi_get("elink.fcgi", {
        "dbfrom": "bioproject", "db": "pubmed", "id": bp_uid, "retmode": "xml",
    })
    root = ET.fromstring(resp.text)
    return [el.text for el in root.findall(".//LinkSetDb/Link/Id") if el.text]


def sra_uid_to_pmids_direct(sra_uid: str) -> list:
    """elink: sra → pubmed (direct link).  Returns list of PMID strings."""
    resp = ncbi_get("elink.fcgi", {
        "dbfrom": "sra", "db": "pubmed", "id": sra_uid, "retmode": "xml",
    })
    root = ET.fromstring(resp.text)
    return [el.text for el in root.findall(".//LinkSetDb/Link/Id") if el.text]


def pmids_from_text_search(accession: str) -> list:
    """
    esearch in PubMed for the accession string.
    Catches papers that mention the accession in abstract/methods but were
    never formally linked in NCBI.  Capped at 5 results to avoid false hits.
    """
    resp = ncbi_get("esearch.fcgi", {
        "db": "pubmed",
        "term": f'"{accession}"',
        "retmax": 5,
        "retmode": "xml",
    })
    root = ET.fromstring(resp.text)
    return [el.text for el in root.findall(".//Id") if el.text]


def fetch_pubmed_details(pmids: list) -> list:
    """
    efetch PubMed XML for a list of PMIDs.
    Returns a list of dicts, one per PMID, with fields:
      pmid, title, authors, first_author, journal, year,
      volume, pages, doi, abstract
    """
    if not pmids:
        return []

    resp = ncbi_get("efetch.fcgi", {
        "db": "pubmed",
        "id": ",".join(pmids),
        "rettype": "xml",
        "retmode": "xml",
    })

    try:
        root = ET.fromstring(resp.text)
    except ET.ParseError:
        log.error("    Could not parse PubMed XML")
        return []

    records = []
    for article in root.findall(".//PubmedArticle"):
        rec: dict = {}

        # PMID
        pmid_el = article.find(".//PMID")
        rec["pmid"] = _safe_strip(pmid_el.text) if pmid_el is not None else None

        # Title
        title_el = article.find(".//ArticleTitle")
        rec["title"] = _safe_strip(title_el.text) if title_el is not None else None

        # Journal
        journal_el = article.find(".//Journal/Title")                   or article.find(".//ISOAbbreviation")
        rec["journal"] = _safe_strip(journal_el.text) if journal_el is not None else None

        # Year  (PubDate may have Year or MedlineDate)
        year_el = article.find(".//PubDate/Year")                or article.find(".//PubDate/MedlineDate")
        if year_el is not None:
            # MedlineDate is like "2021 Jan-Feb"; grab first 4 chars
            rec["year"] = _safe_strip(year_el.text)[:4] if year_el.text else None

        # Volume / pages
        rec["volume"] = _safe_strip(_elem_text(article, ".//Volume"))
        rec["pages"]  = _safe_strip(_elem_text(article, ".//MedlinePgn"))

        # Authors  – collect LastName + Initials for each author
        authors = []
        for author in article.findall(".//AuthorList/Author"):
            last  = _safe_strip(_elem_text(author, "LastName"))
            init  = _safe_strip(_elem_text(author, "Initials"))
            cname = _safe_strip(_elem_text(author, "CollectiveName"))
            if last:
                authors.append(f"{last} {init}".strip() if init else last)
            elif cname:
                authors.append(cname)
        rec["first_author"] = authors[0] if authors else None
        rec["authors"]      = "; ".join(authors) if authors else None

        # DOI
        doi_el = article.find(".//ArticleIdList/ArticleId[@IdType='doi']")
        rec["doi"] = _safe_strip(doi_el.text) if doi_el is not None else None

        # Abstract (first 500 chars to keep CSV compact)
        abstract_parts = [
            _safe_strip(el.text)
            for el in article.findall(".//AbstractText")
            if el.text
        ]
        abstract = " ".join(p for p in abstract_parts if p)
        rec["abstract"] = abstract[:500] if abstract else None

        records.append({k: v for k, v in rec.items() if v is not None})

    return records


def gather_publications(sra_uid: str, accession: str) -> list:
    """
    Run all three publication-finding strategies for a single SRA UID.
    Returns a deduplicated list of PubMed detail dicts.
    """
    pmids: set = set()

    # Strategy 1: direct sra → pubmed elink
    try:
        pmids.update(sra_uid_to_pmids_direct(sra_uid))
    except Exception as e:
        log.warning(f"   sra→pubmed elink failed: {e}")

    # Strategy 2: sra → bioproject → pubmed
    try:
        bp_uid = sra_uid_to_bioproject_uid(sra_uid)
        if bp_uid:
            pmids.update(bioproject_uid_to_pmids(bp_uid))
    except Exception as e:
        log.warning(f"   bioproject→pubmed elink failed: {e}")

    # Strategy 3: PubMed text search (fallback for unlinked accessions)
    try:
        text_hits = pmids_from_text_search(accession)
        pmids.update(text_hits)
    except Exception as e:
        log.warning(f"   PubMed text search failed: {e}")

    if not pmids:
        log.info(f"   No publications found for {accession}")
        return []

    log.info(f"   Found {len(pmids)} PMID(s): {', '.join(sorted(pmids))}")
    return fetch_pubmed_details(list(pmids))

# ─────────────────────────────────────────────────────────────────────────────
#  SAFE HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def _safe_strip(v) -> Optional[str]:
    """Strip a string; return None if the value is None or blank."""
    if v is None:
        return None
    s = str(v).strip()
    return s if s else None


def _normalize_key(s) -> str:
    """
    Lowercase and collapse whitespace/dashes/dots/slashes to underscores.
    Returns empty string for None or blank input (never raises).
    """
    if not s:
        return ""
    return re.sub(r"[\s\-./]+", "_", str(s).strip().lower())


def _elem_text(element, xpath: str) -> Optional[str]:
    """findtext on an element, returning None if missing or blank."""
    if element is None:
        return None
    return _safe_strip(element.findtext(xpath))


# ─────────────────────────────────────────────────────────────────────────────
#  PARSERS
# ─────────────────────────────────────────────────────────────────────────────

def parse_sra_xml(root: ET.Element) -> dict:
    """
    Extract key fields from an SRA experiment XML element.
    Returns a flat dict of raw values prefixed with 'sra__' or 'sra_attr__'.
    """
    d: dict = {}

    d["sra__experiment_title"]   = _elem_text(root, ".//EXPERIMENT/TITLE") \
                                   or _elem_text(root, ".//Title")
    d["sra__design_description"] = _elem_text(root, ".//DESIGN/DESIGN_DESCRIPTION")
    d["sra__library_strategy"]   = _elem_text(root, ".//LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY")
    d["sra__library_source"]     = _elem_text(root, ".//LIBRARY_DESCRIPTOR/LIBRARY_SOURCE")
    d["sra__library_selection"]  = _elem_text(root, ".//LIBRARY_DESCRIPTOR/LIBRARY_SELECTION")

    if root.find(".//LIBRARY_LAYOUT/PAIRED") is not None:
        d["sra__library_layout"] = "PAIRED"
    elif root.find(".//LIBRARY_LAYOUT/SINGLE") is not None:
        d["sra__library_layout"] = "SINGLE"

    # PLATFORM: structure is <PLATFORM><ILLUMINA><INSTRUMENT_MODEL>...</INSTRUMENT_MODEL>
    # The immediate child of <PLATFORM> is the platform-type element (e.g. ILLUMINA, OXFORD_NANOPORE)
    platform_child = root.find(".//PLATFORM/*")
    if platform_child is not None:
        d["sra__platform"] = _safe_strip(platform_child.tag)
        # Instrument model lives one level deeper
        instr = platform_child.find("INSTRUMENT_MODEL")
        if instr is not None:
            d["sra__instrument_model"] = _safe_strip(instr.text)
        else:
            # Some submissions put a flat text in the platform child (rare)
            d["sra__instrument_model"] = _safe_strip(platform_child.text)

    study_elem = root.find(".//STUDY")
    if study_elem is not None:
        d["sra__study_accession"] = _safe_strip(study_elem.get("accession"))
    d["sra__study_title"]    = _elem_text(root, ".//STUDY/DESCRIPTOR/STUDY_TITLE")
    raw_abstract             = _elem_text(root, ".//STUDY/DESCRIPTOR/STUDY_ABSTRACT")
    d["sra__study_abstract"] = raw_abstract[:600] if raw_abstract else None

    sample_elem = root.find(".//SAMPLE")
    if sample_elem is not None:
        d["sra__sample_accession"] = _safe_strip(sample_elem.get("accession"))
    d["sra__organism"] = _elem_text(root, ".//SAMPLE_NAME/SCIENTIFIC_NAME")

    # Sample attributes embedded inside the SRA XML
    for attr in root.findall(".//SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"):
        tag = _safe_strip(attr.findtext("TAG"))
        val = _safe_strip(attr.findtext("VALUE"))
        if tag and val:
            key = "sra_attr__" + _normalize_key(tag)
            d[key] = val

    # Run statistics (first run only)
    run = root.find(".//RUN")
    if run is not None:
        d["sra__run_accession"] = _safe_strip(run.get("accession"))
        d["sra__total_spots"]   = _safe_strip(run.get("total_spots"))
        d["sra__total_bases"]   = _safe_strip(run.get("total_bases"))

    return {k: v for k, v in d.items() if v is not None}


def parse_biosample_xml(bs_elem: ET.Element) -> dict:
    """
    Extract all attributes from a <BioSample> (or <BioSampleSet>) XML element.
    Returns a flat dict with 'bs__' prefix.
    """
    d: dict = {}

    if bs_elem is None:
        return d

    acc = _safe_strip(bs_elem.get("accession"))
    if acc:
        d["bs__accession"] = acc

    org = _elem_text(bs_elem, ".//Organism/OrganismName")
    if org:
        d["bs__organism"] = org

    title = _elem_text(bs_elem, ".//Description/Title")
    if title:
        d["bs__title"] = title

    # All structured attributes (the primary source of growth metadata)
    for attr in bs_elem.findall(".//Attributes/Attribute"):
        name = _safe_strip(attr.get("attribute_name") or attr.get("harmonized_name"))
        val  = _safe_strip(attr.text)
        if name and val:
            key = "bs__" + _normalize_key(name)
            d[key] = val

    return d


# ─────────────────────────────────────────────────────────────────────────────
#  NORMALISATION
# ─────────────────────────────────────────────────────────────────────────────

def _key_matches_synonym(raw_key: str, synonyms: list) -> bool:
    """Check whether a stripped raw key matches any synonym."""
    stripped = re.sub(r"^(bs__|sra_attr__|sra__)", "", raw_key)
    for syn in synonyms:
        if _normalize_key(syn) == stripped:
            return True
    return False


def map_to_target_fields(all_raw: dict) -> dict:
    """
    Scan raw metadata keys for synonym matches.
    BioSample attributes (bs__) take priority over SRA-embedded ones.
    """
    bs_keys   = {k: v for k, v in all_raw.items() if k.startswith("bs__")}
    sra_attrs = {k: v for k, v in all_raw.items() if k.startswith("sra_attr__")}

    normalised: dict = {}
    for target, synonyms in FIELD_SYNONYMS.items():
        for pool in (bs_keys, sra_attrs, all_raw):
            for key, val in pool.items():
                if _key_matches_synonym(key, synonyms):
                    normalised[target] = val
                    break
            if target in normalised:
                break
    return normalised


def extract_from_free_text(all_raw: dict) -> dict:
    """
    Apply regex patterns to free-text fields for fields not yet populated
    by structured synonym matching.
    Returns {field: extracted_value} tagged with '[regex]' prefix.
    """
    text_parts = [
        all_raw.get("sra__experiment_title",   ""),
        all_raw.get("sra__design_description", ""),
        all_raw.get("sra__study_abstract",     ""),
        all_raw.get("bs__title",               ""),
    ]
    combined = " ".join(t for t in text_parts if t)

    extracted: dict = {}
    for field, patterns in FREE_TEXT_PATTERNS.items():
        for pat in patterns:
            m = re.search(pat, combined, re.IGNORECASE)
            if m:
                extracted[field] = m.group(1).strip()
                break
    return extracted


# ─────────────────────────────────────────────────────────────────────────────
#  MAIN RETRIEVAL FUNCTION
# ─────────────────────────────────────────────────────────────────────────────

TARGET_FIELDS = list(FIELD_SYNONYMS.keys())


def fetch_metadata(accession: str) -> dict:
    """
    Full pipeline for one SRA experiment accession.
    Returns a flat dict ready for a DataFrame row.
    """
    row: dict = {"accession": accession}
    log.info(f"Processing  {accession}")

    # ── 1. SRA UID ─────────────────────────────────────────────────────────────
    uid = accession_to_sra_uid(accession)
    if not uid:
        log.warning(f"   No SRA UID found for {accession}")
        row["error"] = "SRA UID not found"
        return row
    row["sra_uid"] = uid

    # ── 2. SRA experiment XML ───────────────────────────────────────────────────
    sra_xml = fetch_sra_xml(uid)
    sra_raw: dict = {}
    if sra_xml is not None:
        sra_raw = parse_sra_xml(sra_xml)
        row.update(sra_raw)
    else:
        log.warning(f"   SRA XML parse failed for {accession}")

    # ── 3. BioSample ────────────────────────────────────────────────────────────
    bs_raw: dict = {}
    biosample_id = sra_uid_to_biosample_id(uid)
    if biosample_id:
        row["biosample_id"] = biosample_id
        bs_xml = fetch_biosample_xml(biosample_id)
        if bs_xml is not None:
            bs_raw = parse_biosample_xml(bs_xml)
            row.update(bs_raw)
    else:
        log.warning(f"   No BioSample linked for {accession}")

    # ── 4. Synonym mapping ───────────────────────────────────────────────────────
    all_raw = {**sra_raw, **bs_raw}
    mapped = map_to_target_fields(all_raw)

    # ── 5. Free-text regex (fill missing fields only) ────────────────────────────
    free_text_hits = extract_from_free_text(all_raw)
    for field, val in free_text_hits.items():
        if field not in mapped:
            mapped[field] = f"[regex] {val}"

    # ── 6. Write target field columns ────────────────────────────────────────────
    for field in TARGET_FIELDS:
        row[f"META__{field}"] = mapped.get(field, "")

    # ── 7. Publications ──────────────────────────────────────────────────────────
    pubs = gather_publications(uid, accession)
    if pubs:
        # Primary publication = first result (by PMID ascending for determinism)
        pubs_sorted = sorted(pubs, key=lambda p: int(p.get("pmid", 0)))
        primary = pubs_sorted[0]
        row["pub__pmid"]         = primary.get("pmid", "")
        row["pub__title"]        = primary.get("title", "")
        row["pub__first_author"] = primary.get("first_author", "")
        row["pub__authors"]      = primary.get("authors", "")
        row["pub__journal"]      = primary.get("journal", "")
        row["pub__year"]         = primary.get("year", "")
        row["pub__doi"]          = primary.get("doi", "")
        row["pub__abstract"]     = primary.get("abstract", "")
        # All PMIDs (semicolon-separated) in case multiple papers are linked
        row["pub__all_pmids"]    = "; ".join(
            p["pmid"] for p in pubs_sorted if p.get("pmid")
        )
        if len(pubs_sorted) > 1:
            row["pub__additional_refs"] = "; ".join(
                f"{p.get('first_author','')} ({p.get('year','')}) {p.get('title','')}"
                for p in pubs_sorted[1:]
            )
    else:
        for col in ["pub__pmid", "pub__title", "pub__first_author", "pub__authors",
                    "pub__journal", "pub__year", "pub__doi", "pub__abstract",
                    "pub__all_pmids"]:
            row[col] = ""

    # ── 8. Completeness score ─────────────────────────────────────────────────
    key_fields = ["strain", "media", "carbon_source", "temperature",
                  "growth_phase", "growth_condition"]
    found = sum(1 for f in key_fields if mapped.get(f))
    row["completeness_score"] = f"{found}/{len(key_fields)}"
    row["needs_manual_review"] = "YES" if found < len(key_fields) // 2 else "maybe"

    return row


# ─────────────────────────────────────────────────────────────────────────────
#  BATCH RUNNER
# ─────────────────────────────────────────────────────────────────────────────

def run_batch(accessions: list, output_csv: str = OUTPUT_CSV) -> pd.DataFrame:
    """
    Fetch metadata for every accession and save progressively to CSV.
    Returns the final DataFrame.
    """
    records: list = []
    total = len(accessions)

    for i, acc in enumerate(accessions):
        acc = acc.strip()
        if not acc or acc.startswith("#"):
            continue
        log.info(f"[{i+1}/{total}]  {acc}")
        try:
            record = fetch_metadata(acc)
        except Exception as exc:
            log.error(f"   FAILED: {exc}")
            record = {"accession": acc, "error": str(exc)}
        records.append(record)
        _save(records, output_csv)

    df = pd.DataFrame(records)
    _save(records, output_csv)
    log.info(f"\nDone. Results saved to {output_csv}")
    _print_summary(df)
    return df


def _save(records: list, path: str) -> None:
    df = pd.DataFrame(records)
    meta_cols = [c for c in df.columns if c.startswith("META__")]
    pub_cols  = [c for c in df.columns if c.startswith("pub__")]
    rest_cols = [c for c in df.columns if c not in meta_cols and c not in pub_cols]
    ordered   = rest_cols[:1] + meta_cols + pub_cols + rest_cols[1:]
    df[ordered].to_csv(path, index=False)


def _print_summary(df: pd.DataFrame) -> None:
    meta_cols = [c for c in df.columns if c.startswith("META__")]
    print("\n── Target field coverage ──────────────────────────────────────────")
    for col in meta_cols:
        filled = df[col].astype(bool).sum()
        field  = col.replace("META__", "")
        print(f"  {field:<22} {filled}/{len(df)} samples filled")
    if "pub__pmid" in df.columns:
        n_pubs = df["pub__pmid"].astype(bool).sum()
        print(f"\n── Publications found: {n_pubs}/{len(df)} accessions have a linked paper")
    print()


# ─────────────────────────────────────────────────────────────────────────────
#  ENTRY POINT
# ─────────────────────────────────────────────────────────────────────────────

# =============================================================================
#  INPUT: Edit ACCESSIONS or pass a text file as argv[1]
# =============================================================================


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Retrieve SRA growth-condition metadata for experiment accessions."
    )
    parser.add_argument(
        "input_file",
        nargs="?",
        default=None,
        help="Text file with one SRA accession per line (optional).",
    )
    parser.add_argument(
        "-o", "--output",
        default=OUTPUT_CSV,
        help=f"Output CSV path (default: {OUTPUT_CSV})",
    )
    parser.add_argument(
        "--email",
        default=NCBI_EMAIL,
        help="Your e-mail for NCBI (required by NCBI policy).",
    )
    parser.add_argument(
        "--api-key",
        default=NCBI_API_KEY,
        help="NCBI API key (optional, raises rate limit to 10 req/s).",
    )
    args = parser.parse_args()

    NCBI_EMAIL   = args.email
    NCBI_API_KEY = args.api_key
    REQUEST_DELAY = 0.11 if NCBI_API_KEY else 0.50

    # Load accessions
    if args.input_file:
        df_metadata = pd.read_csv(args.input_file, sep='\t', index_col=0)
        accessions = df_metadata.index.unique()
        accessions = [k.strip() for k in df_metadata.index.to_list()]
        assert df_metadata.shape[0] == len(accessions)

    if not accessions:
        log.error("No accessions provided. Edit DEFAULT_ACCESSIONS or pass a file.")
        sys.exit(1)

    run_batch(accessions, output_csv=args.output)
