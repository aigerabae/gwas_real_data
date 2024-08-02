#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sys

def main(file_path):
    # Read the Fst data into a pandas DataFrame
    df = pd.read_csv(file_path, sep='\t', header=0)

    # Ensure that both Pop1 and Pop2 are treated as categorical variables
    df['POP1'] = df['POP1'].astype('category')
    df['POP2'] = df['POP2'].astype('category')

    # Get the list of all unique populations
    populations = sorted(set(df['POP1']).union(set(df['POP2'])))

    # Pivot the DataFrame to create a matrix suitable for a heatmap
    heatmap_data = df.pivot(index='POP1', columns='POP2', values='HUDSON_FST')

    # Ensure that all populations are included in the heatmap matrix
    heatmap_data = heatmap_data.reindex(index=populations, columns=populations, fill_value=np.nan)

    # Fill in missing values by copying the upper triangle to the lower triangle
    heatmap_data = heatmap_data.fillna(heatmap_data.T)

    # Optional: Fill NaN values with 0 for visualization if necessary
    heatmap_data = heatmap_data.fillna(0)

    # Sort the values to ensure increasing order in the heatmap
    sorted_columns = heatmap_data.mean().sort_values().index
    heatmap_data = heatmap_data[sorted_columns].reindex(sorted_columns)

    # Create the heatmap using seaborn with three decimal places and smaller font size
    plt.figure(figsize=(10, 8))
    sns.heatmap(
        heatmap_data,
        annot=True,
        fmt=".3f",
        cmap="YlGnBu",
        cbar=True,
        square=True,
        annot_kws={"size": 8}  # Adjust the font size here
    )
    plt.title('Pairwise Fst Heatmap')
    plt.xlabel('POP2')
    plt.ylabel('POP1')

    # Save the plot to a file
    plt.savefig('fst_heatmap_sorted.png')
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python plot_fst_heatmap.py <fst_file>")
        sys.exit(1)
    
    file_path = sys.argv[1]
    main(file_path)
