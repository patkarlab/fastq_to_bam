import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

fastq_file = sys.argv[1]      # File for FASTQ_TO_BAM
dragen_file = sys.argv[2]     # File for Dragen_BAM

df_fastq = pd.read_csv(fastq_file, sep=r"\s+", header=None, names=["Sample", "FASTQ_TO_BAM"])
df_dragen = pd.read_csv(dragen_file, sep=r"\s+", header=None, names=["Sample", "Dragen_BAM"])

df = pd.merge(df_fastq, df_dragen, on="Sample")

plt.figure(figsize=(8, 6))
plt.scatter(df['Dragen_BAM'], df['FASTQ_TO_BAM'], color='blue', edgecolor='k', s=80, label='Samples')

max_val = max(df['Dragen_BAM'].max(), df['FASTQ_TO_BAM'].max())
plt.plot([0, max_val], [0, max_val], 'r--', label='y = x')

slope, intercept = np.polyfit(df['Dragen_BAM'], df['FASTQ_TO_BAM'], 1)
fit_line = slope * df['Dragen_BAM'] + intercept
plt.plot(df['Dragen_BAM'], fit_line, color='green', label=f'Fit: y = {slope:.2f}x + {intercept:.2f}')

corr = df['Dragen_BAM'].corr(df['FASTQ_TO_BAM'])

plt.xlabel('Dragen_BAM Median Coverage')
plt.ylabel('FASTQ_TO_BAM Median Coverage')
plt.title(f'Scatter Plot - Coverview\nPearson r = {corr:.3f}')
plt.legend()

plt.tight_layout()
plt.savefig("scatter_plot_coverview.png", dpi=300)
print("Scatter plot saved to scatter_plot_coverview.png")

