import gzip
import os
import pandas as pd
import numpy as np
import sys
from typing import Tuple, List
from tabulate import tabulate


class Cov_Read:
    """
    A class to read, process, and analyze DNA methylation data from .cov.gz files.
    """
    
    pd.set_option('display.max_rows', None)
   

    sample_list = ["Lung", "Heart", "Liver", "Cortex"]

    def build_df(self, file_path, print_path = False) -> pd.DataFrame:
        '''                
        build_df(filename: str) -> pd.DataFrame:
        Reads a .cov.gz file and returns a DataFrame with methylation data.
        
        Parameters:
            filename (str): The name of the .cov.gz file to be processed.
        
        Returns:
            pd.DataFrame: A DataFrame with columns ['chromosome', 's_loc', 'e_loc',
            'methyl rate', 'methylated_reads', 'unmethylated_reads', 'CG site'].
        '''
        if print_path:
            print(file_path)
            
        # Create empty DataFrame to be filled
        df = pd.DataFrame()
        columns_names = ['chromosome', 's_loc', 'e_loc', 'methyl rate', 'methylated reads', 'unmethylated reads']

        # Open the gzipped file in text mode
        with gzip.open(file_path, 'rt') as file:  
            df = pd.read_csv(file_path, sep="\t", header=None, names=columns_names, low_memory=False)

            # Add Column for CG-Sites
            loc_list = df['s_loc'].tolist()
            cg_list = []
            for idx in range(len(loc_list)):
                # check if index isn't the last value, and check if next methylated read is sequential to current index
                if idx < len(loc_list)-1 and (loc_list[idx] == loc_list[idx+1]-1):
                    cg_list.append(True)
                # check if index isn't the first value, and check if the previous methylated read is sequential to current index
                elif idx > 0 and (loc_list[idx] == loc_list[idx-1]+1):
                    cg_list.append(True)
                # Append false otherwise
                else:
                    cg_list.append(False)
            df['CG site'] = cg_list

        return df
    
    def build_result(self, test = False):
        '''
        Processes all .cov.gz files in the directory and compiles a DataFrame
        with summary results for each file.
        
        Parameters:
            test (bool, optional): If True, process a test file. Defaults to False.
        '''
        column_names = ['id', 'age', 'tissue', 'num_sites', 'ave depth', 'ave methylation', 'ave methylation > 2 depth', 'ave methylation > 5 depth']
        result = pd.DataFrame(columns=column_names)
        
        if not test:
            for filename in os.listdir(self.directory):
                self._build_result_helper(filename, result)
        else: 
            filename = "GSM2465667_M04NB_1wk_Liver.cov.txt.gz"
            self._build_result_helper(filename, result)
            
        return result

    def _build_result_helper(self, filename, result):
        df = self.build_df(filename)
        length = len(df)
        depth = (df['Methylated Reads'] + df['Unmethylated Reads']).mean()
        mean = df.loc[df['cg_site'] == True, 'methyl rate'].mean()
        
        # Get the average methylation for rows with read_count >= 2 & read_count <= 100
        df = df.drop(df[df['unmethylated reads'] + df['methylated reads'] < 2].index)
        df = df.drop(df[df['unmethylated reads'] + df['methylated reads'] > 100].index)
        two_depth_mean = df.loc[df['cg_site'] == True, 'methyl rate'].mean()

        # Get the average methylation for rows with read_count >= 5 & read_count <= 100
        df = df.drop(df[df['Unmethylated Reads'] + df['Methylated Reads'] < 5].index)
        five_depth_mean = df.loc[df['cg_site'] == True, 'methyl rate'].mean()

        # Split the filename into id and age_sample
        splits = filename.split("_")
        id = splits[0]+"_"+splits[1]
        age = splits[2]
        tissue = splits[3].split(".")[0]
        result.loc[len(result)] = [id, age, tissue, length, depth, mean, two_depth_mean, five_depth_mean]
