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
from io import StringIO
import time

def upload():
    st.title("Welcome to the Entropy Calculation Page!")
    st.subheader("Upload the Alignment file in CLUSTAL format (.aln)")
    uploaded_file = st.file_uploader("Upload", type=["clustal"])
    show_submit = st.button("Submit")

    # Initialize the current index
    if "current_index" not in st.session_state:
        st.session_state.current_index = 0

    st.write("""***""")

    if show_submit:
        with st.spinner("Have some Coffee while SAMOSA compiles your dataset"):


