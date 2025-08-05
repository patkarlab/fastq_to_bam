import pandas as pd
import sys
import os

input_file = sys.argv[1]
df = pd.read_csv(input_file)

median = df.iloc[:, 6].median()  
sample = os.path.basename(input_file).replace('.coverview_regions.csv', '')

print(f"{sample}\t{median}")
