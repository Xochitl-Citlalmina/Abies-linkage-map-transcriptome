import pandas as pd
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

# Load files
expression_data = pd.read_csv('allgenesorder.txt', sep='\t')  # Gene expression data
map_data = pd.read_csv('combined_map_data.txt', sep='\t')     # Linkage map

# Separate samples into damaged (DC, DS) and healthy (SC, SS)
damaged_samples = expression_data[expression_data['muestra'].str.startswith(('DC', 'DS'))]
healthy_samples = expression_data[expression_data['muestra'].str.startswith(('SC', 'SS'))]

# Function to compute correlations by cM interval
def correlation_by_interval(samples, linkage_map, expression_threshold=500, cm_interval=1):
    # Filter genes above expression threshold
    filtered = samples[samples['expresion'] > expression_threshold]
    
    # Merge with linkage map to get positions
    genes_with_pos = pd.merge(filtered, linkage_map[['GENE', 'POS_avg', 'CHR']], on='GENE')

    # Compute distances and average expression between gene pairs in same chromosome
    distances = []
    for (gene1, chr1, pos1, expr1), (gene2, chr2, pos2, expr2) in combinations(
        zip(genes_with_pos['GENE'], genes_with_pos['CHR'], genes_with_pos['POS_avg'], genes_with_pos['expresion']), 2):
        if chr1 == chr2:
            dist = abs(pos1 - pos2)
            avg_expr = (expr1 + expr2) / 2
            distances.append((chr1, dist, avg_expr))
    
    distances_df = pd.DataFrame(distances, columns=['LinkageGroup', 'Distance', 'AverageExpression'])
    distances_df['Interval'] = np.floor(distances_df['Distance'] / cm_interval) * cm_interval

    # Compute Pearson correlation in each interval and group
    correlations = []
    for chr_id, group_chr in distances_df.groupby('LinkageGroup'):
        for interval, group in group_chr.groupby('Interval'):
            if len(group) > 1:
                corr, _ = pearsonr(group['Distance'], group['AverageExpression'])
                correlations.append((chr_id, interval, corr))
    
    return pd.DataFrame(correlations, columns=['LinkageGroup', 'Interval', 'Correlation'])

# Loop through desired resolutions
for cm in [1, 5]:
    # Calculate correlations
    corr_damaged = correlation_by_interval(damaged_samples, map_data, cm_interval=cm)
    corr_healthy = correlation_by_interval(healthy_samples, map_data, cm_interval=cm)
    
    # Get linkage groups
    linkage_groups = sorted(corr_damaged['LinkageGroup'].unique())

    # Plot per linkage group
    for lg in linkage_groups:
        plt.figure(figsize=(8, 5))
        damaged = corr_damaged[corr_damaged['LinkageGroup'] == lg]
        healthy = corr_healthy[corr_healthy['LinkageGroup'] == lg]

        plt.plot(damaged['Interval'], damaged['Correlation'], label='Damaged', marker='o', color='red')
        plt.plot(healthy['Interval'], healthy['Correlation'], label='Healthy', marker='o', color='blue')
        
        plt.title(f'Correlation in Linkage Group {lg} ({cm} cM bins)')
        plt.xlabel('Distance (cM)')
        plt.ylabel('Pearson Correlation')
        plt.axhline(y=0.05, color='gray', linestyle='--', label='Threshold 0.05')
        plt.legend()
        plt.tight_layout()
        
        plt.savefig(f'correlation_lg{lg}_{cm}cM.png', dpi=300)
        plt.show()
