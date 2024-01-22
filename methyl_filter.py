import gzip
import os
import pandas as pd
import numpy as np
import sys
import argparse
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

def read_sample_file(sample_file: str, sample_directory: str):
    """
    Reads the 'Bismark coverage' .cov.gz file from the mouse tissue sample

    Returns:
    - pd.DataFrame: methylation data for mouse tissue sample
    """

    reader = Cov_Read()
    df = reader.build_df(sample_file)
    return df


def filter_dataframe(df: pd.DataFrame, chromosome: str, start_loc: int, end_loc: int) -> pd.DataFrame:
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
    filtered_df = df[(df['chromosome'] == str(chromosome)) &
                     (df['s_loc'] >= start_loc) &
                     (df['s_loc'] <= end_loc)]
    return filtered_df

def main():
    # Create argument parser for Setting and Sample file
    parser = argparse.ArgumentParser(
                    prog='methyl_filter',
                    description='returns methylation level of locations in a specified range and chromosome')
    parser.add_argument('SampleFile')
    parser.add_argument('SampleDirectory')
    parser.add_argument('SettingFile')
    args = parser.parse_args()
    
    cov_df = read_sample_file(args.SampleFile, args.SampleDirectory)
    settings_df = read_setting_file(args.SettingFile)
    
    for row in settings_df.itertuples(index=True, name='Pandas'):
        print(row)
        # 'row' is a named tuple
        chromosome = str(row.chromosome)
        start_loc = int(row.start_loc)
        end_loc = int(row.end_loc)

        print(filter_dataframe(cov_df.head(20), chromosome, start_loc, end_loc).head(10))



if __name__ == "__main__":
    main()