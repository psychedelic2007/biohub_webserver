import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os.path
import pathlib
import csv
import math
import zipfile
from Bio import AlignIO
from Bio import SeqIO
from io import StringIO
import time
import io

def feature4():
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

    st.title("Welcome to the Entropy Calculation Page!")
    st.subheader("Upload the sequence in CLUSTAL format (.aln)")
    uploaded_file = st.file_uploader("Upload", type=["aln"])
    show_submit = st.button("Submit")

    #Initialise the current index
    if "current_index" not in st.session_state:
        st.session_state.current_index = 0

    st.write("""***""")

    if show_submit:
        with st.spinner("Have some Coffee while SAMOSA compiles your dataset"):
            st.write("Calculating Entropy.......")
            if uploaded_file is not None:
                stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
                with stringio:
                    align_clustal = AlignIO.read(stringio, "clustal")
                    st.write(align_clustal)

                    def shannon_entropy(list_input):
                        unique_aa = set(list_input)
                        M = len(list_input)
                        entropy_list = []

                        #number of residues in a column
                        for aa in unique_aa:
                            n_i = list_input.count(aa)
                            P_i = n_i/float(M)
                            entropy_i = P_i*(math.log(P_i,2))
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

                    fig, ax = plt.subplots(figsize=(18, 10))
                    ax.plot(clustal_omega, color='green')
                    ax.set_xlabel("Residue Position", fontsize=20)
                    ax.set_ylabel("Shannon Entropy", fontsize=20)
                    ax.tick_params(axis='both', which='major', labelsize=15)

                    buffer = io.BytesIO()
                    fig.savefig(buffer, format="pdf", dpi=300)
                    st.pyplot(fig)
                    st.download_button("Download Image",data=buffer.getvalue(), file_name="entropy_plot.pdf", mime="image/pdf")

                stringio_new = StringIO(uploaded_file.getvalue().decode("utf-8"))
                align = SeqIO.parse(stringio_new,"clustal")
                seq1 = next(align).seq
                csv_filename = "Entropy.csv"
                with open(csv_filename,"w") as csv_file:
                    writer = csv.writer(csv_file, delimiter=',')
                    writer.writerow(['Position','Entropy'])
                    for i in range(len(seq1)):
                        writer.writerow([seq1[i],clustal_omega[i]])
            time.sleep(1)
            st.success("Processing Complete!")
            with open(csv_filename,"rb") as file:
                st.download_button("Download Entropy Data",data=file,file_name="entropy.csv",mime="text/csv")
