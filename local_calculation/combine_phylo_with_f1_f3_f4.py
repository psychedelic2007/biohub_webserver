import os
import pandas as pd

folder1 = "path/to/folder/with/combine/f1_f3_f4"
folder2 = "path/to/phylogenetics/folder"
output_folder = "path/to/output/folder"

# Get a list of csv files in folder 1
csv_files_folder1 = [file for file in os.listdir(folder1) if file.endswith(".csv")]

# combine csv files
for file in csv_files_folder1:
    #write the name (accession_id.csv) of the file that was the last csv file while comparing the phylogenetics.
    #basically the csv file you are left with that you can't compare at the last with any other csv
    if (file == "name_of_csv_that_was_main_ancestor_in_phylogenetics"):
        continue
    else:
        file_path_folder1 = os.path.join(folder1, file)
        file_path_folder2 = os.path.join(folder2, file)
        output_path = os.path.join(output_folder, file)

        # read csv file
        df_folder1 = pd.read_csv(file_path_folder1)
        df_folder2 = pd.read_csv(file_path_folder2)

        # combine data
        combined_df = pd.concat([df_folder1.iloc[:, :4], df_folder2.iloc[:, 1]], axis=1)

        # save combined csv
        combined_df.to_csv(output_path, index=False)