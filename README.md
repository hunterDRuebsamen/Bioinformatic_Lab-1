# Bioinformatic_Lab-1

Sample files must be stored in directory in order to be passed to methyl_filter.py

#### Methyl_filter Tutorial:
  Parameters: SampleFolderPath, SettingsFile.csv

  `Python3 methyl_filter.py data_directory/mousefile.cov.gz input.csv`

#### SettingsFile format (*Currently only accepts .csv files*):

gene_name | chromosome_number | start_location | end_location

"cd3g" | 7 | 165722407 | 165724371 

**NOTE:** chromosome_number, start_location, and end_location are the only required columns in the .csv file. The rest are optional and will be ignored by the program
