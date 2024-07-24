import matplotlib.pyplot as plt
import pandas as pd

# Read the sexcheck file into a DataFrame
df = pd.read_csv('kaz4.sexcheck', sep='\s+')

# Plot histogram of F values
plt.hist(df['F'], bins=20, edgecolor='black')
plt.xlabel('F Value')
plt.ylabel('Frequency')
plt.title('Distribution of F Values')
plt.grid(True)
plt.savefig('f_values_histogram.pdf')
plt.show()

