#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import os

# Read command-line arguments
ibd_file = sys.argv[1]
ethnicity_file = sys.argv[2]

# Load the IBD file with specified dtypes to avoid mixed type warnings
ibd_df = pd.read_csv(ibd_file, delim_whitespace=True, dtype={'FID1': str, 'IID1': str, 'FID2': str, 'IID2': str})

# Load the ethnicity file
ethnicity_df = pd.read_csv(ethnicity_file, sep='\t', header=None, names=['IID', 'Ethnicity', 'Region'])

# Merge the IBD file with the ethnicity file for both IID1 and IID2
merged_df = ibd_df.merge(ethnicity_df[['IID', 'Ethnicity']], left_on='IID1', right_on='IID', how='left') \
                  .merge(ethnicity_df[['IID', 'Ethnicity']], left_on='IID2', right_on='IID', suffixes=('_1', '_2'))

# Ensure symmetry by sorting Ethnicity_1 and Ethnicity_2 columns
merged_df[['Ethnicity_1', 'Ethnicity_2']] = pd.DataFrame(
    np.sort(merged_df[['Ethnicity_1', 'Ethnicity_2']], axis=1), 
    index=merged_df.index
)

# Calculate average PI_HAT between ethnic groups
avg_ibd = merged_df.groupby(['Ethnicity_1', 'Ethnicity_2'])['PI_HAT'].mean().reset_index()

# Create the output file name
output_file = os.path.splitext(ibd_file)[0] + '_average.txt'

# Save the result to the output file
avg_ibd.to_csv(output_file, sep='\t', index=False)

print(f"Average IBD values saved to {output_file}")
