import os
import pandas as pd
from mpi4py import MPI
import argparse
import math
import re
import numpy as np
import sys

def check_and_delete_small_files(tmp_dir):
    """
    Checks the size of all CSV files in the tmp directory and deletes those that are less than 1B in size,
    along with their corresponding S or A files.
    """
    all_files = os.listdir(tmp_dir)
    for file in all_files:
        file_path = os.path.join(tmp_dir, file)
        # Check if the file is a CSV file and its size is less than 1B
        if file.endswith('.csv') and os.path.getsize(file_path) < 1:
            os.remove(file_path)  # Delete the file
            
            # Determine the corresponding file name
            base_name = file.rsplit('_', 1)[0]  # Remove the last part after '_'
            counterpart_suffix = 'A.csv' if '_S.csv' in file else 'S.csv'
            counterpart_file = f"{base_name}_{counterpart_suffix}"
            counterpart_path = os.path.join(tmp_dir, counterpart_file)
            if os.path.exists(counterpart_path):
                os.remove(counterpart_path)  # Delete the corresponding file if it exists
     
def rename_files_sequentially(tmp_dir):
    """
    Renames proc_i_A.csv and proc_i_S.csv files to a sequential order in the given directory.
    """
    a_files = sorted([f for f in os.listdir(tmp_dir) if f.endswith('_A.csv')], key=lambda x: int(x.split('_')[1]))
    s_files = sorted([f for f in os.listdir(tmp_dir) if f.endswith('_S.csv')], key=lambda x: int(x.split('_')[1]))
    
    # Rename A files
    for i, file_name in enumerate(a_files):
        new_name = f"proc_{i}_A.csv"
        os.rename(os.path.join(tmp_dir, file_name), os.path.join(tmp_dir, new_name))
        
    # Rename S files
    for i, file_name in enumerate(s_files):
        new_name = f"proc_{i}_S.csv"
        os.rename(os.path.join(tmp_dir, file_name), os.path.join(tmp_dir, new_name))

def main():
    parser = argparse.ArgumentParser(description='Adjust S and A matrix files across MPI processes.')
    parser.add_argument('-o', '--output_folder', required=True, help='Output folder name')
    parser.add_argument('-n', '--num_cores', type=int, required=True, help='Number of cores')
    args = parser.parse_args()

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    if size != args.num_cores:
        if rank == 0:
            print("Error: The number of MPI processes does not match the number of specified cores.")
        sys.exit()

    tmp_dir = os.path.join(args.output_folder, "tmp")
    
    if rank == 0:
        print("\nChecking and deleting small files.")
        check_and_delete_small_files(tmp_dir)

    # Wait for all processes to finish adjusting matrices
    comm.Barrier()

    # Sequential renaming performed by rank 0
    if rank == 0:
        rename_files_sequentially(tmp_dir)
        print("\nRenaming to sequential order completed.")

    # Ensure the renaming print statement completes before any process exits
    comm.Barrier()

    if rank == 0:
        print("\nAdjustment completed")

if __name__ == "__main__":
    main()