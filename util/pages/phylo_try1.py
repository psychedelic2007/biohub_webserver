import streamlit as st
from Bio import Phylo
from Bio.Seq import Seq
from io import StringIO

def phylo_try1():
    # Define the function to compare the descendant sequence with its ancestor
    def compare_sequences(tree):
        # Initialize an empty list to store the comparison results for all descendants
        comparison_list = []

        # Iterate over each descendant clade in the tree
        for descendant in tree.get_terminals():
            # Get the ancestor node for the current descendant
            ancestor = tree.common_ancestor(descendant)

            # Get the sequences of the descendant and ancestor nodes
            descendant_seq = Seq(descendant.name)
            ancestor_seq = Seq(ancestor.name)

            # Compare the sequences residue-wise
            comparison = [descendant_seq[i] != ancestor_seq[i] for i in range(len(descendant_seq))]

            # Append the comparison result to the list
            comparison_list.append(comparison)

        return comparison_list

    # Set the app title
    st.title("Phylogenetic Sequence Comparison")

    # Allow the user to upload a newick file
    uploaded_file = st.file_uploader("Choose a newick file", type="nh")

    if uploaded_file:
        # Parse the newick file into a Phylo object
        tree = Phylo.read(StringIO(uploaded_file.getvalue().decode()), "newick")

        # Compare the sequences of the descendant node and its ancestor
        comparison = compare_sequences(tree)

        # Display the comparison result as a table
        # result_table = []
        # for i, val in enumerate(comparison):
        #     result_table.append((i + 1, "Different" if val else "Same"))
        #
        # st.write("Comparison result:")
        # st.table(result_table)
        #
        # # Provide a download button for the comparison result as a CSV file
        # csv = "position,difference\n"
        # for i, val in enumerate(comparison):
        #     csv += f"{i + 1},{1 if val else 0}\n"
        # st.download_button(
        #     label="Download comparison result as CSV",
        #     data=csv.encode(),
        #     file_name="comparison_result.csv",
        #     mime="text/csv"
        # )

        # Create a CSV string containing the residue-wise comparison results
        csv = "accession_id,target\n"
        for descendant in tree.get_terminals():
            # Get the ancestor node for the current descendant
            ancestor = tree.common_ancestor(descendant)

            # Get the accession ID of the descendant node
            accession_id = "unknown"
            if "|" in descendant.name:
                accession_id = descendant.name.split("|")[1]

            # Get the sequences of the descendant and ancestor nodes
            descendant_seq = Seq(descendant.name)
            ancestor_seq = Seq(ancestor.name)

            # Compare the sequences residue-wise and append the comparison result to the CSV string
            for i in range(len(descendant_seq)):
                target = "1" if descendant_seq[i] != ancestor_seq[i] else "0"
                csv += f"{accession_id},{target}\n"

        # Display the comparison result as a download link to a CSV file
        st.download_button(
            label="Download comparison result as CSV",
            data=csv.encode(),
            file_name="comparison_result.csv",
            mime="text/csv"
        )

