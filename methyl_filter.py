import gzip
import os
import pandas as pd
import numpy as np
import sys
import argparse
from tabulate import tabulate
from typing import Tuple, List
from cov_read_module import Cov_Read


# cov_reader = cov_read()
# df = cov_reader.build_df("GSM2465633_M00018724_27wk_Liver.cov.txt.gz")

# setting_file = sys.argv[1]
# sample_file = sys.argv[2]

def read_setting_file(setting_file: str):
    """
    Reads the Genome Feature chromosome, start location, and end location from setting_file 

    Returns:
    - pd.DataFrame: Genome Feature locations
    """

    df = pd.read_csv(setting_file)
    return df

def read_sample_file(sample_file_path: str):
    """
    Reads the 'Bismark coverage' .cov.gz file from the mouse tissue sample

    Returns:
    - df (pd.DataFrame): methylation data for mouse tissue sample
    """
    sample_files = []
    for filename in os.listdir(sample_file_path):    
            sample_files.append(filename)
    return sample_files

    
    # reader = Cov_Read()
    # return reader.build_df(sample_file_path)


def filter_dataframe(df: pd.DataFrame, name_str: str, chromosome: str, start_loc: int, end_loc: int) -> pd.DataFrame:
    """
    Returns a filtered DataFrame based on the chromosome and the inclusive range between start_loc and end_loc.

    Parameters:
    - df (pd.DataFrame): The DataFrame to filter.
    - chromosome (int): The chromosome number to filter by.
    - start_loc (int): The inclusive starting location to filter by.
    - end_loc (int): The inclusive ending location to filter by.

    Returns:
    - pd.DataFrame: Filtered DataFrame.
    """
    # Apply the filter conditions, ensuring the range is inclusive
    filtered_df = df[(df['chromosome'] == chromosome) &
                     (df['s_loc'] >= start_loc) &
                     (df['s_loc'] <= end_loc)].copy() # Use copy() to defined filtered_df as copy, not view and avoid SettingWithCopyWarning
    
    # Trim CG_Sites column from file dump
    filtered_df.drop(labels= 'CG site', axis=1, inplace= True)

    # TODO: change cov_read_module to simply not read past the first 2 digits rather than round here
    filtered_df['methyl rate'] = filtered_df['methyl rate'].round(2)

    # set name column equal to single value
    filtered_df['name'] = name_str
 
    column_order = ['name', 'chromosome', 's_loc', 'e_loc',
                    'methyl rate', 'methylated reads', 'unmethylated reads']
    
    filtered_df = filtered_df[column_order]

    return filtered_df

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
    for row in settings_df.itertuples(index=True, name='Pandas'):
        name = row.name
        chromosome = str(row.chromosome)
        start_loc = int(row.start_loc)
        end_loc = int(row.end_loc)

        for idx, sample in enumerate(samples_directory):
            final_path = f"{args.SampleDirectoryPath}/{sample}"
            cov_df = reader.build_df(final_path)
            filter_dataframe(cov_df, name, chromosome, start_loc, end_loc).to_csv(f'{name}_out_{idx+1}.csv', index=False, header=True)

        # with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        #     print(filter_dataframe(cov_df, name, chromosome, start_loc, end_loc).to_csv(index=False, header=False))


if __name__ == "__main__":
    main()