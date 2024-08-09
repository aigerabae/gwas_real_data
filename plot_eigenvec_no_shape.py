import sys
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
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
    'PC2': [float(line.split()[3]) for line in lines],
    'sample_id': [line.split()[0] for line in lines],
    'ethnicity': [ethnicities[line.split()[0]] for line in lines],
    'region': [ethnicities_df[ethnicities_df['sample_id'] == line.split()[0]]['region'].values[0] for line in lines]
}

df = pd.DataFrame(data)

# Determine unique areas and create dynamic mappings
unique_regions = df['region'].unique()
unique_ethnicities = df['ethnicity'].unique()

# Define a set of shape symbols (more than enough to cover various regions)
shape_symbols = ['circle', 'square', 'triangle-up', 'triangle-down', 'square-open', 'diamond',
                 'cross', 'x', 'star', 'hexagram', 'triangle-left', 'triangle-right', 
                 'pentagon', 'hourglass', 'bowtie', 'circle-open']

# Use color palette to assign colors dynamically
color_palette = px.colors.qualitative.Plotly

# Create mappings for regions to shapes and ethnicities to colors
shape_map = dict(zip(unique_regions, shape_symbols[:len(unique_regions)]))
color_map = dict(zip(unique_ethnicities, color_palette[:len(unique_ethnicities)]))

# Debug print statements
print("Shape Map:", shape_map)
print("Color Map:", color_map)

# Create interactive scatter plot
fig = px.scatter(df, x='PC1', y='PC2', color='ethnicity', symbol='region', 
                 hover_data={'sample_id': True, 'ethnicity': True})  # Display ethnicity on hover

# Update marker shapes and colors based on dynamic mapping
for region, shape in shape_map.items():
    fig.update_traces(marker=dict(symbol=shape), selector=dict(mode='markers+text', legendgroup=region))
    
for ethnicity in unique_ethnicities:
    color = color_map.get(ethnicity)
    fig.update_traces(marker=dict(color=color), selector=dict(name=ethnicity))

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
    xaxis=dict(showgrid=False, zeroline=False),  # Reverse x-axis
    yaxis=dict(showgrid=False, zeroline=False),  # Reverse y-axis
    plot_bgcolor='rgba(0,0,0,0)'  # Make background transparent
)

# Save the plot as an HTML file
fig.write_html('interactive_plot.html')

# Show the plot in a browser
fig.show()