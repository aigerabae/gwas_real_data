import numpy as np
import pandas as pd
import sys
import os

def create_percentages_table(admixture_file, ethnicities_file):
    # Load ADMIXTURE data
    admixture_data = np.loadtxt(admixture_file)
    
    # Load ethnicities data (assuming ethnicities are in the 2nd column)
    ethnicities_df = pd.read_csv(ethnicities_file, delim_whitespace=True, header=None)
    
    # Extract ethnicities (from the 2nd column)
    ethnicities = ethnicities_df.iloc[:, 1].values
    
    # Ensure that the number of ethnicities matches the number of rows in admixture data
    if len(ethnicities) != admixture_data.shape[0]:
        raise ValueError("Number of ethnicities does not match number of rows in ADMIXTURE data")
    
    # Create a DataFrame for the ADMIXTURE data
    k = admixture_data.shape[1]  # Number of ancestry components (K)
    admixture_df = pd.DataFrame(admixture_data, columns=[f'K{i+1}' for i in range(k)])
    admixture_df['Ethnicity'] = ethnicities
    
    # Group by ethnicity and calculate the mean (in percentage) for each ancestry component
    averages_df = admixture_df.groupby('Ethnicity').mean() * 100  # Convert to percentage
    
    # Reset index so that 'Ethnicity' becomes a column
    averages_df.reset_index(inplace=True)
    
    # Rename the columns to 'Ethnicity', '1', '2', ..., 'k'
    new_column_names = ['Ethnicity'] + [f'{i+1}' for i in range(k)]
    averages_df.columns = new_column_names
    
    # Create output filename based on the input .Q file
    output_file = f"{os.path.splitext(admixture_file)[0]}_percentages.tsv"
    
    # Save the DataFrame to a TSV file
    averages_df.to_csv(output_file, sep='\t', index=False)
    print(f"Output saved to: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python create_percentages_table.py <admixture_file> <ethnicities_file>")
        sys.exit(1)

    admixture_file = sys.argv[1]
    ethnicities_file = sys.argv[2]
    create_percentages_table(admixture_file, ethnicities_file)

