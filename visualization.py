import pysam
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Set non-interactive backend

import matplotlib.pyplot as plt
# Your plotting code here

# Load VCF file using pysam
vcf_file = 'filtered_variants.vcf'
vcf = pysam.VariantFile(vcf_file, 'r')

# Extract variant quality scores
quality_scores = [record.qual for record in vcf if record.qual is not None]

# Create a DataFrame for better plotting and analysis
df = pd.DataFrame({
    'Quality Score': quality_scores
})

# Plot the distribution of quality scores
plt.figure(figsize=(10, 6))
plt.hist(df['Quality Score'], bins=30, color='skyblue', edgecolor='black')
plt.title('Distribution of Variant Quality Scores')
plt.xlabel('Quality Score')
plt.ylabel('Frequency')
plt.grid(axis='y', alpha=0.75)

# Save the plot as a PNG file
plt.savefig('variant_quality_distribution.png', dpi=300, bbox_inches='tight')

# Show the plot (optional)
plt.show()
