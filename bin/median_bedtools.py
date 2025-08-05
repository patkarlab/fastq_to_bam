import pandas as pd
import sys
import os

tsv = sys.argv[1]

df = pd.read_csv(tsv, delim_whitespace=True, header=None)

median = df[4].median()
sample = os.path.basename(tsv).replace('.counts.bed', '')

print(f"{sample}\t{median}")
