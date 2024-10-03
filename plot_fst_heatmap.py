#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sys

def main(fst_file_path, sorting_order_file_path):
    # Read the Fst data into a pandas DataFrame
    df = pd.read_csv(fst_file_path, sep='\t', header=0)

    # Read the sorting order from the sorting_order_fst.txt file
    with open(sorting_order_file_path, 'r') as f:
        sorting_order = [line.strip() for line in f.readlines()]

    # Ensure that both Pop1 and Pop2 are treated as categorical variables
    df['POP1'] = df['#POP1'].astype('category')
    df['POP2'] = df['POP2'].astype('category')

    # Pivot the DataFrame to create a matrix suitable for a heatmap
    heatmap_data = df.pivot(index='POP1', columns='POP2', values='HUDSON_FST')

    # Ensure that all populations are included in the heatmap matrix
    heatmap_data = heatmap_data.reindex(index=sorting_order, columns=sorting_order, fill_value=np.nan)

    # Fill in missing values by copying the upper triangle to the lower triangle
    heatmap_data = heatmap_data.fillna(heatmap_data.T)

    # Optional: Fill NaN values with 0 for visualization if necessary
    heatmap_data = heatmap_data.fillna(0)

    # Create the heatmap using seaborn with three decimal places and smaller font size
    plt.figure(figsize=(10, 8))
    ax = sns.heatmap(
        heatmap_data,
        annot=True,
        fmt=".3f",
        cmap="YlGnBu",
        cbar=True,
        square=True,
        annot_kws={"size": 2}  # Adjust the font size here
    )
    plt.title('Pairwise Fst Heatmap')
    plt.xlabel('POP2')
    plt.ylabel('POP1')

    # Manually set the tick positions and labels
    ax.set_xticks(np.arange(len(sorting_order)) + 0.5)
    ax.set_yticks(np.arange(len(sorting_order)) + 0.5)
    ax.set_xticklabels(sorting_order, fontsize=4, rotation=90)
    ax.set_yticklabels(sorting_order, fontsize=4, rotation=0)

    # Display x-axis labels on both the top and bottom
    ax.xaxis.tick_top()
    ax.tick_params(axis='x', which='both', top=True, bottom=True, labeltop=True, labelbottom=True)

    # Save the plot to PNG and SVG formats with high resolution
    plt.savefig('fst_heatmap_sorted.png', dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.savefig('fst_heatmap_sorted.svg', bbox_inches='tight', pad_inches=0.1)

    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python plot_fst_heatmap.py <fst_file> <sorting_order_file>")
        sys.exit(1)
    
    fst_file_path = sys.argv[1]
    sorting_order_file_path = sys.argv[2]
    main(fst_file_path, sorting_order_file_path)
