import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load the TSV file into a pandas DataFrame
file_path = 'domain_comparison_results_human.tsv'
df = pd.read_csv(file_path, sep='\t')

# Transform the e-value to logarithmic scale
df['log_hmm_evalue'] = np.log10(df['hmm_evalue'] + 1e-300)  # Adding a small constant to avoid log(0)

# Plot the density plot
plt.figure(figsize=(10, 6))
sns.kdeplot(
    x=df['log_hmm_evalue'], 
    y=df['overlap_percentage'], 
    fill=True, 
    cmap='viridis',
    bw_adjust=0.5,  # Adjust the bandwidth to see if it helps
    levels=100     # Increase the number of contour levels for more detail
)

plt.xlabel('Log10(Evalue)')
plt.ylabel('Overlap Percentage')
plt.title('Density Plot of Overlap vs Log10(Evalue)')

# Save the plot as a PNG file with 300 dpi
output_file = 'density_human.png'
plt.savefig(output_file, dpi=300)

plt.show()