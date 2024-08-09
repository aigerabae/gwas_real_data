import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
from scipy.cluster.hierarchy import linkage, leaves_list

def main(Q_file, ethnicities_file):
    # Load the ancestry proportions data from the .Q file
    data = np.loadtxt(Q_file)

    # Load the sample information from the ethnicities file
    ethnicities = pd.read_csv(ethnicities_file, sep='\t', header=None, names=['ID', 'Ethnicity', 'Population'])

    # Number of individuals and populations
    num_individuals, num_populations = data.shape

    # Merge the data and ethnicities based on their order
    ethnicities['Ancestry'] = list(data)

    # Create a DataFrame where each ancestry component is a separate column
    ancestry_df = pd.DataFrame(data, columns=[f'Population {i+1}' for i in range(num_populations)])
    ethnicities = pd.concat([ethnicities, ancestry_df], axis=1)

    # Calculate the mean ancestry proportions for each population
    mean_proportions = ethnicities.groupby('Population').mean().iloc[:, 3:]  # exclude ID, Ethnicity, Population columns

    # Perform hierarchical clustering to determine the order of populations based on similarity
    linkage_matrix = linkage(mean_proportions, method='average')
    ordered_population_indices = leaves_list(linkage_matrix)
    ordered_populations = mean_proportions.index[ordered_population_indices]

    # Reorder the ethnicities DataFrame based on the ordered populations
    ethnicities['Population'] = pd.Categorical(ethnicities['Population'], categories=ordered_populations, ordered=True)
    ethnicities = ethnicities.sort_values(by=['Population', 'Ethnicity'])

    # Split data into two parts for two subplots
    split_index = num_individuals // 2
    part1 = ethnicities.iloc[:split_index]
    part2 = ethnicities.iloc[split_index:]

    # Plot the data in two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(35, 6))  # Adjust the figsize as needed for width and height

    def plot_admixture(ax, data, title):
        bottom = np.zeros(len(data))
        grouped_x = []
        x = np.arange(len(data))
        
        for population_name, population_group in data.groupby('Population'):
            for ethnicity_name, ethnicity_group in population_group.groupby('Ethnicity'):
                grouped_x.extend(x[:len(ethnicity_group)])
                x = x[len(ethnicity_group):]
        
        for i in range(num_populations):
            ancestry_component = data[f'Population {i+1}']
            ax.bar(grouped_x, ancestry_component, width=1.0, bottom=bottom, color=colors[i % len(colors)], label=f'Population {i+1}')
            bottom += ancestry_component
        
        ax.set_xlabel('Individuals grouped by Ethnicity and Population')
        ax.set_ylabel('Ancestry Proportion')
        ax.set_title(title)
        
        ethnicities_labels = []
        current_pos = 0
        ethnicity_pos = 0

        for population_name, population_group in data.groupby('Population'):
            for ethnicity_name, ethnicity_group in population_group.groupby('Ethnicity'):
                label_pos = ethnicity_pos + len(ethnicity_group) / 2
                ethnicities_labels.append((label_pos, f'{ethnicity_name}'))
                current_pos += len(ethnicity_group)
                ethnicity_pos += len(ethnicity_group)
                ax.axvline(x=current_pos, color='black', linewidth=1.5)

        ax.set_xticks([pos for pos, label in ethnicities_labels])
        ax.set_xticklabels([label for pos, label in ethnicities_labels], rotation=90, fontsize=8)
    
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'orange', 'purple', 'pink', 'brown', 'gray', 'black']
    
    # Plot the first half
    plot_admixture(ax1, part1, 'Ancestry Proportions - Part 1')
    
    # Plot the second half
    plot_admixture(ax2, part2, 'Ancestry Proportions - Part 2')
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python plot_admixture.py <Q_file> <ethnicities_file>")
        sys.exit(1)

    Q_file = sys.argv[1]
    ethnicities_file = sys.argv[2]
    
    main(Q_file, ethnicities_file)
