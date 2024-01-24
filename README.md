# Bioinformatic_Lab-1

#### Methyl_filter Tutorial:
  Parameters: SampleFilePath, SettingFile.csv

  `Python3 methyl_filter data_directory/mousefile.cov.gz input.csv`

#### SettingFile format (*Currently only accepts .csv files*):

gene_name | chromosome_number | start_location | end_location

"cd3g" | 7 | 165722407 | 165724371 

**NOTE:** chromosome_number, start_location, and end_location are the only required columns in the .csv file. The rest are optional and will be ignored by the program
