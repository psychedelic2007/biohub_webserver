import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os.path
import pathlib
from collections import Counter
import csv
import math
import zipfile
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
from Bio.SeqRecord import SeqRecord
from io import StringIO
import time
import random

def feature3():
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

    with open("quotes.txt", "r", encoding="utf-8") as f:
        file_text = f.read()
        quotes = file_text.split("\n\n")

    st.title("Welcome to the Feature Calculation Page!")
    st.subheader("Upload the FASTA file")
    uploaded_file = st.file_uploader("Upload", type=["fasta"])
    show_submit = st.button("Submit")

    # Initialize the current index
    if "current_index" not in st.session_state:
        st.session_state.current_index = 0

    st.write("""***""")

    if show_submit:
        with st.spinner("Have some Coffee while SAMOSA compiles your dataset"):
            aa_list = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
            fcaa = []
            with zipfile.ZipFile("output.zip", "w") as output_zip:
                if(uploaded_file is not None):
                    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
                    sequence_file = list(SeqIO.parse(stringio, "fasta"))

                    for seq in sequence_file:
                        sequence = seq.seq
                        A = sequence.count("A")
                        C = sequence.count("C")
                        D = sequence.count("D")
                        E = sequence.count("E")
                        F = sequence.count("F")
                        G = sequence.count("G")
                        H = sequence.count("H")
                        I = sequence.count("I")
                        K = sequence.count("K")
                        L = sequence.count("L")
                        M = sequence.count("M")
                        N = sequence.count("N")
                        P = sequence.count("P")
                        Q = sequence.count("Q")
                        R = sequence.count("R")
                        S = sequence.count("S")
                        T = sequence.count("T")
                        V = sequence.count("V")
                        W = sequence.count("W")
                        Y = sequence.count("Y")

                        A_cur = A / len(sequence)
                        C_cur = C / len(sequence)
                        D_cur = D / len(sequence)
                        E_cur = E / len(sequence)
                        F_cur = F / len(sequence)
                        G_cur = G / len(sequence)
                        H_cur = H / len(sequence)
                        I_cur = I / len(sequence)
                        K_cur = K / len(sequence)
                        L_cur = L / len(sequence)
                        M_cur = M / len(sequence)
                        N_cur = N / len(sequence)
                        P_cur = P / len(sequence)
                        Q_cur = Q / len(sequence)
                        R_cur = R / len(sequence)
                        S_cur = S / len(sequence)
                        T_cur = T / len(sequence)
                        V_cur = V / len(sequence)
                        W_cur = W / len(sequence)
                        Y_cur = Y / len(sequence)

                        A_fut = (((12 * A) + (2 * D) + (2 * E) + (4 * G) + (4 * P) + (4 * S) + (4 * T) + (
                                        4 * V)) / (36 * len(sequence)))
                        C_fut = (((2 * R) + (2 * C) + (2 * G) + (2 * F) + (4 * S) + (2 * W) + (2 * Y)) / (
                                        18 * len(sequence)))
                        D_fut = (((2 * A) + (2 * N) + (2 * D) + (4 * E) + (2 * G) + (2 * H) + (2 * Y) + (2 * V)) / (
                                        18 * len(sequence)))
                        E_fut = (((2 * A) + (4 * D) + (2 * E) + (2 * Q) + (2 * G) + (2 * K) + (2 * V)) / (
                                        18 * len(sequence)))
                        F_fut = (((2 * C) + (2 * I) + (6 * L) + (2 * F) + (2 * S) + (2 * Y) + (2 * V)) / (
                                        18 * len(sequence)))
                        G_fut = (((4 * A) + (6 * R) + (2 * D) + (2 * C) + (2 * E) + (12 * G) + (2 * S) + (1 * W) + (
                                        4 * V)) / (36 * len(sequence)))
                        H_fut = (((2 * R) + (2 * N) + (2 * D) + (4 * Q) + (2 * H) + (2 * L) + (2 * P) + (2 * Y)) / (
                                        18 * len(sequence)))
                        I_fut = (((1 * R) + (2 * N) + (6 * I) + (4 * L) + (1 * K) + (3 * M) + (2 * F) + (2 * S) + (
                                        3 * T) + (3 * V)) / (27 * len(sequence)))
                        K_fut = (((2 * R) + (4 * N) + (2 * E) + (2 * Q) + (1 * I) + (2 * K) + (1 * M) + (2 * T)) / (
                                        18 * len(sequence)))
                        L_fut = (((4 * R) + (2 * Q) + (2 * H) + (4 * I) + (18 * L) + (2 * M) + (6 * F) + (4 * P) + (
                                        2 * S) + (1 * W) + (6 * V)) / (54 * len(sequence)))
                        M_fut = (((1 * R) + (3 * I) + (2 * L) + (1 * K) + (1 * T) + (1 * V)) / (9 * len(sequence)))
                        N_fut = (((2 * N) + (2 * D) + (2 * H) + (2 * I) + (4 * K) + (2 * S) + (2 * T) + (2 * Y)) / (
                                        18 * len(sequence)))
                        P_fut = (((4 * A) + (4 * D) + (2 * Q) + (2 * H) + (4 * L) + (12 * P) + (4 * S) + (
                                        4 * T)) / (36 * len(sequence)))
                        Q_fut = (((2 * R) + (2 * E) + (2 * Q) + (4 * H) + (2 * L) + (2 * K) + (2 * P)) / (
                                        18 * len(sequence)))
                        R_fut = (((18 * R) + (2 * C) + (2 * Q) + (6 * G) + (2 * H) + (1 * I) + (4 * L) + (2 * K) + (
                                        1 * M) + (4 * P) + (6 * S) + (2 * T) + (2 * W)) / (54 * len(sequence)))
                        S_fut = (((4 * A) + (6 * R) + (2 * N) + (4 * C) + (2 * G) + (2 * I) + (2 * L) + (2 * F) + (
                                        4 * P) + (14 * S) + (6 * T) + (1 * W) + (2 * Y)) / (54 * len(sequence)))
                        T_fut = (((4 * A) + (2 * R) + (2 * N) + (3 * I) + (2 * K) + (1 * M) + (4 * P) + (6 * S) + (
                                        12 * T)) / (36 * len(sequence)))
                        V_fut = (((4 * A) + (2 * D) + (2 * E) + (4 * G) + (3 * I) + (6 * L) + (1 * M) + (2 * F) + (
                                        12 * V)) / (36 * len(sequence)))
                        W_fut = (((2 * R) + (2 * C) + (1 * G) + (1 * L) + (1 * S)) / (9 * len(sequence)))
                        Y_fut = (((2 * N) + (2 * D) + (2 * C) + (2 * H) + (2 * F) + (2 * S) + (2 * Y)) / (
                                        18 * len(sequence)))

                        aa_curr_count_list = [A_cur, C_cur, D_cur, E_cur, F_cur, G_cur, H_cur, I_cur, K_cur, L_cur,
                                                  M_cur, N_cur, P_cur, Q_cur, R_cur, S_cur, T_cur, V_cur, W_cur, Y_cur]
                        aa_fut_count_list = [A_fut, C_fut, D_fut, E_fut, F_fut, G_fut, H_fut, I_fut, K_fut, L_fut,
                                                 M_fut, N_fut, P_fut, Q_fut, R_fut, S_fut, T_fut, V_fut, W_fut, Y_fut]

                        fcaa = []
                        for i in range(len(aa_curr_count_list)):
                            if (aa_curr_count_list[i] == 0):
                                ratio = 0
                                fcaa.append(ratio)
                            else:
                                ratio = aa_fut_count_list[i] / aa_curr_count_list[i]
                                fcaa.append(ratio)

                    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
                    sequence_file = list(SeqIO.parse(stringio, "fasta"))
                    for seq in sequence_file:
                        aa_dict2 = dict(zip(aa_list,fcaa))
                        csv_filename = seq.id + ".csv"
                        with open(csv_filename, "w") as csv_file:
                            writer = csv.writer(csv_file, delimiter=',')
                            writer.writerow(['Sequence','F3'])
                            for residue in seq:
                                writer.writerow([residue,aa_dict2.get(residue)])
                        output_zip.write(csv_filename)
                        time.sleep(1)
            output_zip.close()
            st.success("Processing Complete!")
            with open("output.zip","rb") as fp:
                btn = st.download_button(label="Download", data=fp, file_name="output.zip", mime="application/zip")

        st.write("")
        st.write("""***""")
        quote = random.choice(quotes)
        st.write(quote)
        st.write("""***""")