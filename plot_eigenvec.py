import sys
import pandas as pd
import plotly.express as px
import numpy as np

# Check if the correct number of arguments are provided
if len(sys.argv) != 3:
    print("Usage: python plot_eigenvec.py <eigenvec_file> <ethnicities_file>")
    sys.exit(1)

# Read the eigenvec file, skipping the first line
eigenvec_file = sys.argv[1]
with open(eigenvec_file, 'r') as f:
    lines = f.readlines()[1:]

# Read the ethnicities file
ethnicities_file = sys.argv[2]
ethnicities_df = pd.read_csv(ethnicities_file, delim_whitespace=True, header=None, names=['sample_id', 'ethnicity', 'region'])
ethnicities = ethnicities_df.set_index('sample_id').to_dict()['ethnicity']

# Extract data from the eigenvec file
data = {
    'PC1': [-float(line.split()[2]) for line in lines],
    'PC2': [-float(line.split()[3]) for line in lines],
    'sample_id': [line.split()[0] for line in lines],
    'ethnicity': [ethnicities[line.split()[0]] for line in lines],
    'region': [ethnicities_df[ethnicities_df['sample_id'] == line.split()[0]]['region'].values[0] for line in lines]
}

df = pd.DataFrame(data)

# Add a new column for marker size
df['marker_size'] = np.where(df['ethnicity'] == 'Kazakh', 14, 6)  # 14 for 'Kazakh', 6 for others

# Sort the dataframe by region and then by ethnicity
df = df.sort_values(by=['region', 'ethnicity'], ascending=True)

# Determine unique regions and ethnicities
unique_regions = df['region'].unique()
unique_ethnicities = df['ethnicity'].unique()

# Define a set of shape symbols (ensure enough for various regions)
shape_symbols = ['circle', 'square', 'triangle-up', 'triangle-down', 'diamond', 'cross', 
                 'x', 'star', 'hexagram', 'triangle-left', 'triangle-right', 'pentagon', 
                 'hourglass', 'bowtie', 'circle-open']

# Map regions to shape symbols
shape_map = dict(zip(unique_regions, shape_symbols[:len(unique_regions)]))

# Use color palette to assign colors dynamically
color_palette = px.colors.qualitative.Plotly
color_map = dict(zip(unique_ethnicities, color_palette[:len(unique_ethnicities)]))

# Create interactive scatter plot with dynamic marker size, shape, and color
fig = px.scatter(df, x='PC1', y='PC2', color='ethnicity', symbol='region', 
                 hover_data={'sample_id': True, 'ethnicity': True},
                 size='marker_size',  # Use the marker_size column to adjust sizes
                 size_max=14,  # Maximum size for Kazakh
                 symbol_map=shape_map,  # Map regions to shapes
                 color_discrete_map=color_map,  # Map ethnicities to colors
                 category_orders={'region': list(unique_regions), 'ethnicity': list(unique_ethnicities)})  # Ensure legend order

# Update layout for larger font sizes and remove the grid
fig.update_layout(
    title='PCA Plot',
    title_font_size=24,
    xaxis_title='PC1',
    xaxis_title_font_size=20,
    yaxis_title='PC2',
    yaxis_title_font_size=20,
    legend_title_font_size=18,
    font=dict(size=16),
    xaxis=dict(showgrid=False, zeroline=False),
    yaxis=dict(showgrid=False, zeroline=False),
    plot_bgcolor='rgba(0,0,0,0)'  # Make background transparent
)

# Save the plot as an HTML file
fig.write_html('interactive_plot.html')

# Show the plot in a browser
fig.show()
