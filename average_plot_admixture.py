import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
from scipy.cluster.hierarchy import linkage, leaves_list

def plot_admixture(admixture_file, ethnicities_file):
    # Load ADMIXTURE data
    admixture_data = np.loadtxt(admixture_file)
    
    # Load ethnicities data
    ethnicities_df = pd.read_csv(ethnicities_file, delim_whitespace=True, header=None)
    
    # Extract ethnicities (assuming ethnicities are in the 2nd column)
    ethnicities = ethnicities_df.iloc[:, 1].values
    
    # Ensure that the number of ethnicities matches the number of rows in admixture data
    if len(ethnicities) != admixture_data.shape[0]:
        raise ValueError("Number of ethnicities does not match number of rows in ADMIXTURE data")
    
    # Create a DataFrame for easy manipulation
    df = pd.DataFrame(admixture_data, columns=[f'Ancestry {i+1}' for i in range(admixture_data.shape[1])])
    df['Ethnicity'] = ethnicities

    # Group by ethnicity and calculate the mean for each ancestry
    grouped_df = df.groupby('Ethnicity').mean()
    
    # Perform hierarchical clustering to find similarity between rows
    ancestry_data = grouped_df.values
    linkage_matrix = linkage(ancestry_data, method='ward')  # Use Ward's method for clustering
    ordered_indices = leaves_list(linkage_matrix)  # Get the order of indices after clustering
    
    # Reorder the grouped_df based on the clustering results
    grouped_df = grouped_df.iloc[ordered_indices]
    
    # Plotting
    fig, ax = plt.subplots(figsize=(12, 8))
    bar_width = 0.9  # Control the width of the bars to remove white space
    bottom = np.zeros(grouped_df.shape[0])
    
    for column in grouped_df.columns:  # Loop through ancestry columns
        ax.bar(grouped_df.index, grouped_df[column], bottom=bottom, width=bar_width, label=column)
        bottom += grouped_df[column].values

    plt.xlabel('Ethnicity')
    plt.ylabel('Average Ancestry Proportion')
    plt.title('Average ADMIXTURE Results by Ethnicity (Sorted by Similarity)')
    plt.xticks(rotation=90, ha='center')
    plt.legend(title='Ancestry', bbox_to_anchor=(1.05, 1))
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python plot_admixture.py <admixture_file> <ethnicities_file>")
        sys.exit(1)

    admixture_file = sys.argv[1]
    ethnicities_file = sys.argv[2]
    plot_admixture(admixture_file, ethnicities_file)
