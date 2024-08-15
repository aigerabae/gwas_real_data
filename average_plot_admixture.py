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
    
    # Extract populations and labels (assuming populations are in the 2nd column and grouping in the 3rd column)
    ethnicities = ethnicities_df.iloc[:, 1].values
    populations = ethnicities_df.iloc[:, 2].values
    
    # Ensure that the number of populations matches the number of rows in admixture data
    if len(populations) != admixture_data.shape[0]:
        raise ValueError("Number of populations does not match number of rows in ADMIXTURE data")
    
    # Create a DataFrame for easy manipulation
    df = pd.DataFrame(admixture_data, columns=[f'Ancestry {i+1}' for i in range(admixture_data.shape[1])])
    df['Ethnicity'] = ethnicities
    df['Population'] = populations

    # Group by population and calculate the mean for each ancestry
    grouped_df = df.groupby(['Population', 'Ethnicity']).mean()
    
    # Sort by population
    grouped_df = grouped_df.sort_index(level=0)
    
    # Perform hierarchical clustering to find similarity between rows
    ancestry_data = grouped_df.values
    linkage_matrix = linkage(ancestry_data, method='ward')  # Use Ward's method for clustering
    ordered_indices = leaves_list(linkage_matrix)  # Get the order of indices after clustering
    
    # Reorder the grouped_df based on the clustering results
    grouped_df = grouped_df.iloc[ordered_indices]
    
    # Flatten the multi-index for plotting
    grouped_df.reset_index(inplace=True)
    
    # Create a combined label for the x-axis
    grouped_df['Labels'] = grouped_df['Population'] + ' - ' + grouped_df['Ethnicity']
    
    # Plotting
    fig, ax = plt.subplots(figsize=(12, 8))
    bar_width = 0.9  # Control the width of the bars to remove white space
    bottom = np.zeros(grouped_df.shape[0])
    
    for column in grouped_df.columns[2:-1]:  # Exclude 'Population', 'Ethnicity', and 'Labels' columns
        ax.bar(grouped_df['Labels'], grouped_df[column], bottom=bottom, width=bar_width, label=column)
        bottom += grouped_df[column].values

    plt.xlabel('Ethnicity (Grouped by Population)')
    plt.ylabel('Average Ancestry Proportion')
    plt.title('Average ADMIXTURE Results by Ethnicity and Population (Sorted by Similarity)')
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
