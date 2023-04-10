import streamlit as st
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import AlignIO
from io import StringIO
import warnings
warnings.filterwarnings("ignore")

def mrf():
    st.markdown(
        """
        <style>
        [data-testid="stSidebar"][aria-expanded="true"]{
            background-color: #48a2cf;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )

    st.title("MRF")

    # First ask user how many files they would like to upload
    num_files = st.selectbox("How many files would you like to upload?", list(range(1, 13)))

    uploaded_files = []

    # Initialize session state variables
    for i in range(num_files):
        st.session_state[f"month_{i}"] = ""
        st.session_state[f"year_{i}"] = ""

    # def upload_align_file(month, year):
    #     st.write(f"Uploaded file processed for {month} {year}")

    st.write("""***""")

    for i in range(num_files):
        col1, col2 = st.columns([1,1])
        with col1:
            file = st.file_uploader(f"Upload file {i+1}", type="aln")
            if file:
                uploaded_files.append(file)
        with col2:
            month = st.text_input(f"Month for file {i+1}")
            year = st.text_input(f"Year for file {i+1}")

    if st.button("Submit"):
        st.write("""***""")
        with st.spinner("Have some Coffee while SAMOSA compiles your dataset"):
            st.write("Calculating Entropy.......")
            if file is not None:
                for i in range(num_files):
                    stringio = StringIO(file.getvalue().decode("utf-8"))
                    with stringio:
                        align_clustal = AlignIO.read(stringio, "clustal")
                        #st.write(align_clustal)

                        def shannon_entropy(list_input):
                            unique_aa = set(list_input)
                            M = len(list_input)
                            entropy_list = []

                            # number of residues in a column
                            for aa in unique_aa:
                                n_i = list_input.count(aa)
                                P_i = n_i / float(M)
                                entropy_i = P_i * (math.log(P_i, 2))
                                entropy_list.append(entropy_i)
                            sh_entropy = -(sum(entropy_list))
                            return sh_entropy

                        def shannon_entropy_list_msa(alignment_file):
                            shannon_entropy_list = []
                            for col_no in range(len(list(alignment_file[0]))):
                                list_input = list(alignment_file[:, col_no])
                                shannon_entropy_list.append(shannon_entropy(list_input))
                            return shannon_entropy_list

                        clustal_omega = shannon_entropy_list_msa(align_clustal)

                        month = st.session_state[f"month_{i+1}"]
                        year = st.session_state[f"year_{i+1}"]

                        plt.figure(figsize=(8,6))
                        plt.plot(clustal_omega)
                        plt.title(f"Entropy Plot for {month} {year}")
                        plt.xlabel("Residue Position")
                        plt.ylabel("Entropy")
                        plt.legend([f"{month} {year}"])
                        plt.show()
