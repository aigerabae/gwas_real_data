import sys
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy.spatial import ConvexHull
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
with open(ethnicities_file, 'r') as f:
    ethnicities = {line.split()[0]: line.split()[1] for line in f.readlines()}

# Extract data from the file
data = {
    'PC1': [float(line.split()[2]) for line in lines],
    'PC2': [float(line.split()[3]) for line in lines],
    'sample_id': [line.split()[0] for line in lines],
    'ethnicity': [ethnicities[line.split()[0]] for line in lines]
}

df = pd.DataFrame(data)

# Define marker shapes and colors based on ethnicity
shape_map = {
    'Uygur': 'circle', 'Kazakh': 'circle', 'Hazara': 'circle',
    'Russian': 'square', 'French': 'square', 'Basque': 'square', 'Bergamo': 'square',
    'Pathan': 'triangle-up', 'Sindhi': 'triangle-up', 'Kalash': 'triangle-up',
    'Adygei': 'diamond',
    'Bedouin': 'square-open', 'Mozabite': 'square-open',
    'Japanese': 'triangle-down', 'Northern': 'triangle-down', 'Mongolian': 'triangle-down', 'Yakut': 'triangle-down', 'Han': 'triangle-down'
}

color_map = {
    'circle': 'rgba(150,0,100,0.6)',   # Darker Pink
    'square': 'rgba(100,100,0,0.6)',        # Darker Blue
    'triangle-up': 'rgba(100,0,150,0.6)',    # Darker Green
    'diamond': 'rgba(0,150,50,0.6)',       # Darker Yellow
    'square-open': 'rgba(50,200,100,0.6)',   # Darker Beige
    'triangle-down': 'rgba(200,200,100,0.6)'  # Darker Khaki
}

# Create interactive scatter plot
fig = px.scatter(df, x='PC1', y='PC2', color='ethnicity', symbol='ethnicity', hover_data=['sample_id', 'ethnicity'])

# Update marker shapes and colors based on ethnicity
for ethnicity, shape in shape_map.items():
    fig.update_traces(marker=dict(symbol=shape), selector=dict(name=ethnicity))

# Function to add convex hulls with custom names
def add_convex_hull(df, shape, color, legend_name):
    shape_df = df[df['ethnicity'].map(lambda x: shape_map[x] == shape)]
    if shape_df.shape[0] < 3:  # Need at least 3 points to form a convex hull
        return
    
    points = shape_df[['PC1', 'PC2']].values
    try:
        hull = ConvexHull(points)
        # Add convex hull trace
        hull_points = np.append(points[hull.vertices], [points[hull.vertices][0]], axis=0)  # close the hull
        fig.add_trace(go.Scatter(
            x=hull_points[:, 0],
            y=hull_points[:, 1],
            mode='lines',
            line=dict(color=color, width=1),  # Thinner lines
            name=legend_name  # Set custom legend name
        ))
    except:
        print(f"Warning: Convex hull computation failed for shape: {shape}")

# Add convex hulls for each shape type with custom legend names
shape_legends = {
    'circle': 'Central Asia',
    'square': 'Europe',
    'triangle-up': 'South Asia',
    'triangle-down': 'East Asia',
    'square-open': 'Middle East'
}

for shape, legend_name in shape_legends.items():
    add_convex_hull(df, shape, color_map[shape], legend_name)

# Update layout for larger font sizes
fig.update_layout(
    title='PCA Plot',
    title_font_size=24,
    xaxis_title='PC1',
    xaxis_title_font_size=20,
    yaxis_title='PC2',
    yaxis_title_font_size=20,
    legend_title_font_size=18,
    font=dict(size=16)
)

# Save the plot as an HTML file
fig.write_html('interactive_plot.html')

# Show the plot in a browser
fig.show()
