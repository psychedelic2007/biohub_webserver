import os
import pandas as pd

# Path to the folder containing the CSV files. This folder is the one generated using "combine_f1_f3.py"
folder_path = 'path/to/combine/folder/'

# Path to the main CSV file ("A")
main_file_path = 'path/to/entropy/csv/file'

# Read the main CSV file ("entropy.csv")
main_df = pd.read_csv(main_file_path)

# Iterate over the files in the folder
for filename in os.listdir(folder_path):
    if filename.endswith('.csv'):
        file_path = os.path.join(folder_path, filename)

        # Read each file in the folder
        df = pd.read_csv(file_path)

        # Get the length of the second column of the current file
        column_length = len(df.iloc[:, 1])

        # Add a new column with the values from the second column of "entropy.csv" up to the length of the current file's column
        df['F4'] = main_df.iloc[:column_length, 1].values

        # Write the updated DataFrame back to the file
        df.to_csv(file_path, index=False)
