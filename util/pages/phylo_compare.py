import streamlit as st
import pandas as pd
import os
import zipfile
from io import StringIO

def phylo_compare():

    def compare_sequences(seq1, seq2):
        comparison_results = []
        for i in range(len(seq1)):
            if(seq1[i] == seq2[i]):
                comparison_results.append([i, 0])
            else:
                comparison_results.append([i, 1])

        df = pd.DataFrame(comparison_results, columns=["Residue Position", "Target"])
        return df

    def process_fasta_file(file):
        sequences = {}
        with open(file, "r") as f:
            lines = f.readlines()
            current_seq = ""
            current_id = ""
            for line in lines:
                if(line.startswith(">")):
                    if(current_seq != ""):
                        sequences[current_id] = current_seq
                        current_seq = ""
                    current_id = line.strip()[1:]
                else:
                    current_seq += line.strip()
            sequences[current_id] = current_seq

        csv_files = []
        for seq_id, seq in sequences.items():
            df = compare_sequences(seq, st.session_state[seq_id])
            csv_filename = seq_id + ".csv"
            df.to_csv(csv_filename, index=False)
            csv_filename.append(csv_filename)

        with zipfile.ZipFile("reuslts.zip", "w") as zipf:
            for csv_file in csv_files:
                zipf.write(csv_file)

        for csv_file in csf_files:
            os.remove(csv_file)

        return "reuslts.zip"

    st.title("Seq Compare")
    fasta_file = st.file_uploader("upload file", type=["fasta","fa"])

    if(fasta_file is not None):
        stringio = StringIO(fasta_file.getvalue().decode("utf-8"))
        records = list(SeqIO.parse(stringio, "fasta"))
        st.write("Enter Accession ID:")
        for i in range(len(records)):
            accession_id = st.text_input(f"Sequence {i+1} Accession ID")
            st.session_state[f"seq{i+1}"] = accession_id
            
        if st.button("Submit"):
            zip_filename = process_fasta_file(fasta_file.name)
            st.write(f"Process complete download")
            st.download_button("Download", data=zip_filename)
