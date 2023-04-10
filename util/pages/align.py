import streamlit as st
from Bio import AlignIO
from Bio import SeqIO
from io import StringIO

def align():
    st.title("Multiple Sequence Alignment Viewer")

    # Display instructions
    st.write("Please upload a multiple sequence alignment file in Clustal format (.aln)")

    # Create file uploader
    file = st.file_uploader("Upload file", type=["aln"])

    # Display file contents if file is uploaded
    if file is not None:
        stringio = StringIO(file.getvalue().decode("utf-8"))
        with stringio:
            records = AlignIO.read(stringio, "clustal")
            st.write(records)
