import streamlit as st
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import AlignIO
from io import StringIO
import base64
import warnings
import zipfile
import os
import io
import random
warnings.filterwarnings("ignore")
st.set_option('deprecation.showPyplotGlobalUse', False)

def entropy():
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

    with open("quotes.txt", "r", encoding="utf-8") as f:
        file_text = f.read()
        quotes = file_text.split("\n\n")

    # First ask user how many files they would like to upload
    num_files = st.selectbox("How many files would you like to upload?", list(range(1, 13)))
    #plt.rcParams["figure.figsize"] = [8, 2]
    plt.rcParams["figure.autolayout"] = True

    uploaded_files = []
    months = []
    years = []

    # Initialize session state variables
    for i in range(num_files):
        st.session_state[f"month_{i}"] = ""
        st.session_state[f"year_{i}"] = ""

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
            months.append(month)
            years.append(year)

    def create_zip_file(df):
        zip_file = io.BytesIO()
        with zipfile.ZipFile(zip_file, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for i in range(num_files):
                entropy_file = pd.DataFrame({"Residue Position": range(1, len(clustal_omega) + 1), "Shannon Entropy": clustal_omega})
                entropy_file.to_csv(f'entropy_data_{months[i]}_{years[i]}.csv', index=False)
                zipf.write(f'entropy_data_{months[i]}_{years[i]}.csv')
        return zip_file.getvalue()

    if st.button("Submit"):
        with st.spinner("Have some Coffee while SAMOSA compiles your dataset"):
            for i, file in enumerate(uploaded_files):
                file_contents = file.getvalue().decode("utf-8")
                stringio = StringIO(file_contents)
                with stringio:
                    align_clustal = AlignIO.read(stringio, "clustal")
                    st.write(align_clustal)

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

                    if num_files == 1:
                        plt.plot(clustal_omega, label=f"{months[i]} {years[i]}")
                        plt.xlabel("Residue Position")
                        plt.ylabel("Shannon Entropy")
                        plt.legend()
                        st.pyplot()
                        plt.savefig(f"entropy_{months[i]}_{years[i]}.png", dpi=300)
                        buffer = io.BytesIO()
                        plt.savefig(buffer, format="png", dpi=300)
                        st.download_button("Download Image", data=buffer.getvalue(), file_name=f"entropy_plot_{months[i]}_{years[i]}.png",mime="image/png")
                    else:
                        plt.plot(clustal_omega, label=f"{months[i]} {years[i]}")
                        plt.xlabel("Residue Position")
                        plt.ylabel("Shannon Entropy")
                        plt.legend()
                        st.pyplot()
                        plt.savefig(f"entropy_{months[i]}_{years[i]}.png", dpi=300)
                        buffer = io.BytesIO()
                        plt.savefig(buffer, format="png", dpi=300)
                        st.download_button("Download Image", data=buffer.getvalue(),
                                           file_name=f"entropy_plot_{months[i]}_{years[i]}.png", mime="image/png", key=i)

            st.pyplot()

            for i in range(num_files):
                entropy_file = pd.DataFrame({"Residue Position": range(1, len(clustal_omega) + 1), "Shannon Entropy": clustal_omega})
                zip_file = create_zip_file(entropy_file)

            with st.spinner("Please wait while SAMOSA compiles your model training data......"):
                st.success("Data compiled!")
                st.download_button(label="Download Entropy Data", data=zip_file, file_name="entropy_output.zip", mime="application/zip")

        st.write("")
        st.write("""***""")
        quote = random.choice(quotes)
        st.write(quote)
        st.write("""***""")

