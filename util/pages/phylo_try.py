import streamlit as st
from ete3 import Tree
import pandas as pd
import os
import zipfile
from io import StringIO

def phylo_try():
    def read_newick_file(file):
        with open(file, 'r') as f:
            stringio = StringIO.read(f.getvalue().decode("utf-8"))
            newick = stringio.read()
        tree = Tree(newick)
        return tree

    def compare_sequences(seq1, seq2):
        result = ''
        for i in range(len(seq1)):
            if seq1[i] == seq2[i]:
                result += '0'
            else:
                result += '1'
        return result

    def process_data(file):
        # Read the newick file
        tree = read_newick_file(file)

        # Initialize an empty dictionary to store the results
        results = {}

        # Iterate through each node in the tree
        for node in tree.traverse():
            # Check if the node has a sequence attribute
            if 'sequence' in node.features:
                # Check if the node has a parent
                if node.up:
                    # Get the sequence of the parent node
                    parent_seq = node.up.sequence
                    # Get the sequence of the current node
                    node_seq = node.sequence
                    # Compare the sequences
                    comparison = compare_sequences(parent_seq, node_seq)
                    # Add the comparison to the results dictionary
                    results[node.name] = comparison

        # Convert the results dictionary to a DataFrame
        df = pd.DataFrame.from_dict(results, orient='index', columns=['Comparison'])

        return df

    file = st.file_uploader('Upload Newick File', type=['nh'])

    # Check if a file has been uploaded
    if file:
        # Process the data
        df = process_data(file)

        # Add a download button for the CSV file
        csv = df.to_csv(index_label='Node')
        button = st.download_button(
            label='Download CSV',
            data=csv,
            file_name='tree_comparison.csv',
            mime='text/csv'
        )

        # Display the results
        st.dataframe(df)