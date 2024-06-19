import pandas as pd
from io import StringIO


# Reading the data into a DataFrame
df = pd.read_csv("CpgIsland-GRCm39.txt", sep='\t')

df = df.drop(columns=["name"])

# Renaming columns
df = df.rename(columns={
    "#bin": "name",
    "chrom": "chromosome",
    "chromStart": "start_loc",
    "chromEnd": "end_loc"
})

# Dropping unnecessary columns
df = df[['name', 'chromosome', 'start_loc', 'end_loc']]

# Removing the "chr" from the start of every entry in the chromosome column
df['chromosome'] = df['chromosome'].str.replace('chr', '')

df['name'] = df.apply(lambda row: f"{row.name}-{row['name']}", axis=1)


df.to_csv('CpgIsland-GRCm39-settings.txt', index=False)

print("Saved CpgIsland-GRCm39-settings.txt to file")
