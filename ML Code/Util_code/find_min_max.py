import os
import csv
import os
import pandas as pd

# Define the directory path where the CSV files are located
directory = './Data/Test'



# Initialize a dictionary to store the max and min values for each column
column_values = {}

# Iterate over all files in the directory
for filename in os.listdir(directory):
    if filename.endswith('.csv'):
        file_path = os.path.join(directory, filename)
        
        # Read the CSV file using pandas
        df = pd.read_csv(file_path)
        
        # Iterate over each column
        for column in df.columns:
            # Update the max and min values for each column
            if column not in column_values:
                column_values[column] = {'max': float('-inf'), 'min': float('inf')}
            column_values[column]['max'] = max(column_values[column]['max'], df[column].max())
            column_values[column]['min'] = min(column_values[column]['min'], df[column].min())

# Print the max and min values for each column
for column, values in column_values.items():
    print(f"Column: {column}")
    print(f"Max Value: {values['max']}")
    print(f"Min Value: {values['min']}")
    print()

