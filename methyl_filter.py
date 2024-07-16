
import gzip
import os
import pandas as pd
import numpy as np
import sys
import argparse
from tabulate import tabulate
from typing import Tuple, List
from cov_read_module import Cov_Read
import time
import cudf


# cov_reader = cov_read()
# df = cov_reader.build_df("GSM2465633_M00018724_27wk_Liver.cov.txt.gz")

# setting_file = sys.argv[1]
# sample_file = sys.argv[2]

def read_setting_file(setting_file: str):
    """
    Reads the Genome Feature chromosome, start location, and end location from setting_file 

    Returns:
    - cudf.DataFrame: Genome Feature locations
    """

    df = cudf.read_csv(setting_file)
    return df

def read_sample_file(sample_file_path: str):
    """
    Reads the 'Bismark coverage' .cov.gz file from the mouse tissue sample

    Returns:
    - df (cudf.DataFrame): methylation data for mouse tissue sample
    """
    sample_files = []
    for filename in os.listdir(sample_file_path):    
            sample_files.append(filename)
    return sample_files

    
    # reader = Cov_Read()
    # return reader.build_df(sample_file_path)


def filter_dataframe(df: cudf.DataFrame, name_str: str, chromosome: str, start_loc: int, end_loc: int) -> cudf.DataFrame:
    """
    Returns a filtered DataFrame based on the chromosome and the inclusive range between start_loc and end_loc.
    """

    # Work with a filtered copy of the DataFrame
    filtered_df = df[(df['chromosome'] == chromosome) &
                     (df['s_loc'] >= start_loc) &
                     (df['s_loc'] <= end_loc)].copy()

    # Modify the copy
    filtered_df['methyl rate'] = filtered_df['methyl rate'].round(2)
    filtered_df['probe'] = name_str

    column_order = ['probe', 'chromosome', 's_loc', 'e_loc', 'methyl rate', 'methylated reads', 'unmethylated reads']
    return filtered_df[column_order]


def main():
    # Create argument parser for Setting and Sample file
    parser = argparse.ArgumentParser(
                    prog='methyl_filter',
                    description='returns methylation level of locations in a specified range and chromosome')
    parser.add_argument('SampleDirectoryPath')
    parser.add_argument('SettingFile')
    args = parser.parse_args()
    

    samples_directory = read_sample_file(args.SampleDirectoryPath)
    settings_df = read_setting_file(args.SettingFile)

    reader = Cov_Read()

    
    '''
    For each region in the settings_file, 
    generate an output for each sample in the sample directory
    '''
    start_time = time.perf_counter()
    for idx, sample in enumerate(samples_directory):
        final_path = f"{args.SampleDirectoryPath}/{sample}"
        settings_df["chromosome"] = settings_df["chromosome"].astype(str)
        settings_df["start_loc"] = settings_df["start_loc"].astype(int)
        settings_df["end_loc"] = settings_df["end_loc"].astype(int)

        # Initialize an empty list to store filtered dataframes
        filtered_dfs = []


        cov_df = reader.build_df(final_path)
        settings_length = len(settings_df)
        counter = 0
        for row in settings_df.to_pandas().itertuples(index=True, name='Pandas'):
            
            probe = row.name
            filtered_df = filter_dataframe(cov_df, probe, row.chromosome, row.start_loc, row.end_loc) # .to_csv(f'{sample}_out_{idx+1}.csv', index=False, header=True)
            filtered_dfs.append(filtered_df)
            counter += 1
            print(f"iteration: {counter}/{settings_length}")

        final_df = cudf.concat(filtered_dfs, ignore_index=True)
        final_df.to_csv(f'{sample}_{idx+1}.csv', index=False, header=True)

        # Additional code for the new CSV
        summary_df = final_df.groupby(['probe', 'chromosome']).agg({
            's_loc': ['min', 'max'],
            'methyl rate': 'mean'
        }).reset_index()
        print(summary_df.tail())

        # Flatten the MultiIndex columns
        summary_df.columns = ['probe', 'chromosome', 's_loc', 'e_loc', 'average_methylation_rate']
        summary_df["average_methylation_rate"] = summary_df["average_methylation_rate"].round(2)
        # Save the summary dataframe to a new CSV file
        summary_df.to_csv(f'{sample}_summary_{idx+1}.csv', index=False, header=True)
        # print(f"Time elapsed: {(time.perf_counter()-start_time):.3f}s")

if __name__ == "__main__":
    main()