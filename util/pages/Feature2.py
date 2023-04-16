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

def feature2():
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

    def count_partitions(s, amino):
        #print("Total length of sequence is ::", len(s))
        #print("The number of A in sequence is ::", s.count(amino))

        partition = round((len(s) / s.count(amino)))
        g = []
        chunks, chunk_size = len(s), partition

        for i in range(0, chunks, chunk_size):
            x = s[i:i + chunk_size]
            g.append(x)

        count_0 = 0
        count_1 = 0
        count_2 = 0
        count_3 = 0
        count_4 = 0
        count_5 = 0
        c0 = []
        c1 = []
        c2 = []
        c3 = []
        c4 = []
        c5 = []

        for z in range(0, len(g)):
            b = g[z].count(amino)
            if (b == 0):
                count_0 += 1
                c0.append(count_0)
            elif (b == 1):
                count_1 += 1
                c1.append(count_1)
            elif (b == 2):
                count_2 += 1
                c2.append(count_2)
            elif (b == 3):
                count_3 += 1
                c3.append(count_3)
            elif (b == 4):
                count_4 += 1
                c4.append(count_4)
            elif (b == 5):
                count_5 += 1
                c5.append(count_5)

        #print("The partition containing 0", amino, "is ::", len(c0))
        #print("The partition containing 1", amino, "is ::", len(c1))
        #print("The partition containing 2", amino, "is ::", len(c2))
        #print("The partition containing 3", amino, "is ::", len(c3))
        #print("The partition containing 4", amino, "is ::", len(c4))
        #print("The partition containing 5", amino, "is ::", len(c5))

        def factorial(n):
            if (n == 1 or n == 0):
                return 1
            else:
                return n * factorial(n - 1)

        r = s.count(amino)
        n = len(g)
        first = factorial(r) / (
                    factorial(len(c0)) * factorial(len(c1)) * factorial(len(c2)) * factorial(len(c3)) * factorial(
                len(c4)))
        second_den = []
        for i in range(0, len(g)):
            fact = factorial(g[i].count(amino))
            second_den.append(fact)

        product = 1
        for item in second_den:
            product = product * item
        second = factorial(r) / product
        third = n ** (-r)
        final = first * second * third
        return final

    def calculate_max_proba(s, amino):
        #print("Total length of sequence is ::", len(s))
        #print("The number of A in sequence is ::", s.count(amino))
        aa_count = []
        new = s.count(amino)
        aa_count.append(new)

        def factorial(n):
            if (n == 1 or n == 0):
                return 1
            else:
                return n * factorial(n - 1)

        def distribution_probability(R, r, q, n):
            first_part = factorial(R)
            for qi in q:
                first_part = first_part / factorial(qi)

            second_part = factorial(R)
            for ri in r:
                second_part = second_part / factorial(ri)

            third_part = n ** (-R)

            return (first_part * second_part * third_part)

        def integer_partition(n):
            partitions = []

            def generate_partitions(n, max_val, current_partition):
                if (n == 0):
                    partitions.append(current_partition)
                else:
                    for i in range(1, min(n, max_val) + 1):
                        generate_partitions(n - i, i, current_partition + [i])

            generate_partitions(n, n, [])
            return partitions

        max_pdp = []

        for i in range(len(aa_count)):
            R = aa_count[i]
            n = aa_count[i]

            partitions = integer_partition(n)
            for p in partitions:
                while (len(p) < n):
                    p.append(0)

            q_list = []
            for p in partitions:
                q_temp = []
                for i in range(len(p) + 1):
                    qu = p.count(i)
                    q_temp.append(qu)
                q_list.append(q_temp)

            final_dp = []
            for p, a in zip(partitions, q_list):
                r = [i for i in p]
                q = [i for i in a]
                prob = distribution_probability(R, r, q, n)
                final_dp.append(prob)
            max_pdp.append(max(final_dp))
        return max_pdp

    with open("quotes.txt", "r", encoding="utf-8") as f:
        file_text = f.read()
        quotes = file_text.split("\n\n")

    st.subheader("Upload the FASTA file")
    uploaded_file = st.file_uploader("Upload", type=["fasta"])
    show_submit = st.button("Submit")

    # Initialize the current index
    if "current_index" not in st.session_state:
        st.session_state.current_index = 0

    st.write("""***""")

    if show_submit:
        with st.spinner("Have some Coffee while SAMOSA compiles your dataset"):
            with zipfile.ZipFile("output.zip", "w") as output_zip:
                if(uploaded_file is not None):
                    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
                    sequence_file = list(SeqIO.parse(stringio, "fasta"))
                    aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S',
                                       'T', 'V', 'W', 'Y']

                    for seq in sequence_file:
                        final_op = []
                        final_dp = []
                        for aa in aa_list:
                            count_op = count_partitions(seq,aa)
                            count_dp = calculate_max_proba(seq,aa)
                            final_op.append(count_op)
                            final_dp.append(count_dp)

                        final_value = []
                        for i in range(len(final_op)):
                            x = np.log(final_dp[i][0]/final_op[i])
                            final_value.append(x)

                        aa_dict1 = dict(zip(aa_list,final_value))
                        csv_filename = seq.id + ".csv"
                        with open(csv_filename, "w") as csv_file:
                            writer = csv.writer(csv_file, delimiter=',')
                            writer.writerow(['Sequence','F2'])
                            for residue in seq:
                                writer.writerow([residue,aa_dict1.get(residue)])

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
