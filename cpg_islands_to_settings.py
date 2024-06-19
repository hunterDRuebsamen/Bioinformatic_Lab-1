import pandas as pd
from io import StringIO


# Reading the data into a DataFrame
df = pd.read_csv("cpgIsland-GRCm38.txt", sep='\t')

df = df.drop(columns=["name"])

# Renaming columns
df = df.rename(columns={
    "#bin": "#probe",
    "chrom": "chromosome",
    "chromStart": "start_loc",
    "chromEnd": "end_loc"
})

# Dropping unnecessary columns
df = df[['#probe', 'chromosome', 'start_loc', 'end_loc']]

# Removing the "chr" from the start of every entry in the chromosome column
df['chromosome'] = df['chromosome'].str.replace('chr', '')

df['#probe'] = df.apply(lambda row: f"{row.name}-{row['#probe']}", axis=1)


df.to_csv('cpgIsland-GRCm38-settings.txt', index=False)

print("Saved CpgIsland-GRCm38-settings.txt to file")
