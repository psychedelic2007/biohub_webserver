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
import requests

def feature1():
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
            max_pdp = []
            fcaa = []
            with zipfile.ZipFile("output.zip", "w") as output_zip:
                if(uploaded_file is not None):
                    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
                    sequence_file = list(SeqIO.parse(stringio, "fasta"))

                    aa_pair = ['AA', 'AC', 'AD', 'AE', 'AF', 'AG', 'AH', 'AI', 'AK', 'AL', 'AM', 'AN', 'AP','AQ', 'AR', 'AS', 'AT', 'AV', 'AW', 'AY',
                               'CA', 'CC', 'CD', 'CE', 'CF', 'CG', 'CH', 'CI', 'CK', 'CL', 'CM', 'CN', 'CP','CQ', 'CR', 'CS', 'CT', 'CV', 'CW', 'CY',
                               'DA', 'DC', 'DD', 'DE', 'DF', 'DG', 'DH', 'DI', 'DK', 'DL', 'DM', 'DN', 'DP','DQ', 'DR', 'DS', 'DT', 'DV', 'DW', 'DY',
                               'EA', 'EC', 'ED', 'EE', 'EF', 'EG', 'EH', 'EI', 'EK', 'EL', 'EM', 'EN', 'EP','EQ', 'ER', 'ES', 'ET', 'EV', 'EW', 'EY',
                               'FA', 'FC', 'FD', 'FE', 'FF', 'FG', 'FH', 'FI', 'FK', 'FL', 'FM', 'FN', 'FP','FQ', 'FR', 'FS', 'FT', 'FV', 'FW', 'FY',
                               'GA', 'GC', 'GD', 'GE', 'GF', 'GG', 'GH', 'GI', 'GK', 'GL', 'GM', 'GN', 'GP','GQ', 'GR', 'GS', 'GT', 'GV', 'GW', 'GY',
                               'HA', 'HC', 'HD', 'HE', 'HF', 'HG', 'HH', 'HI', 'HK', 'HL', 'HM', 'HN', 'HP','HQ', 'HR', 'HS', 'HT', 'HV', 'HW', 'HY',
                               'IA', 'IC', 'ID', 'IE', 'IF', 'IG', 'IH', 'II', 'IK', 'IL', 'IM', 'IN', 'IP','IQ', 'IR', 'IS', 'IT', 'IV', 'IW', 'IY',
                               'KA', 'KC', 'KD', 'KE', 'KF', 'KG', 'KH', 'KI', 'KK', 'KL', 'KM', 'KN', 'KP','KQ', 'KR', 'KS', 'KT', 'KV', 'KW', 'KY',
                               'LA', 'LC', 'LD', 'LE', 'LF', 'LG', 'LH', 'LI', 'LK', 'LL', 'LM', 'LN', 'LP','LQ', 'LR', 'LS', 'LT', 'LV', 'LW', 'LY',
                               'MA', 'MC', 'MD', 'ME', 'MF', 'MG', 'MH', 'MI', 'MK', 'ML', 'MM', 'MN', 'MP','MQ', 'MR', 'MS', 'MT', 'MV', 'MW', 'MY',
                               'NA', 'NC', 'ND', 'NE', 'NF', 'NG', 'NH', 'NI', 'NK', 'NL', 'NM', 'NN', 'NP','NQ', 'NR', 'NS', 'NT', 'NV', 'NW', 'NY',
                               'PA', 'PC', 'PD', 'PE', 'PF', 'PG', 'PH', 'PI', 'PK', 'PL', 'PM', 'PN', 'PP','PQ', 'PR', 'PS', 'PT', 'PV', 'PW', 'PY',
                               'QA', 'QC', 'QD', 'QE', 'QF', 'QG', 'QH', 'QI', 'QK', 'QL', 'QM', 'QN', 'QP','QQ', 'QR', 'QS', 'QT', 'QV', 'QW', 'QY',
                                       'RA', 'RC', 'RD', 'RE', 'RF', 'RG', 'RH', 'RI', 'RK', 'RL', 'RM', 'RN', 'RP',
                                       'RQ', 'RR', 'RS', 'RT', 'RV', 'RW', 'RY',
                                       'SA', 'SC', 'SD', 'SE', 'SF', 'SG', 'SH', 'SI', 'SK', 'SL', 'SM', 'SN', 'SP',
                                       'SQ', 'SR', 'SS', 'ST', 'SV', 'SW', 'SY',
                                       'TA', 'TC', 'TD', 'TE', 'TF', 'TG', 'TH', 'TI', 'TK', 'TL', 'TM', 'TN', 'TP',
                                       'TQ', 'TR', 'TS', 'TT', 'TV', 'TW', 'TY',
                                       'VA', 'VC', 'VD', 'VE', 'VF', 'VG', 'VH', 'VI', 'VK', 'VL', 'VM', 'VN', 'VP',
                                       'VQ', 'VR', 'VS', 'VT', 'VV', 'VW', 'VY',
                                       'WA', 'WC', 'WD', 'WE', 'WF', 'WG', 'WH', 'WI', 'WK', 'WL', 'WM', 'WN', 'WP',
                                       'WQ', 'WR', 'WS', 'WT', 'WV', 'WW', 'WY',
                                       'YA', 'YC', 'YD', 'YE', 'YF', 'YG', 'YH', 'YI', 'YK', 'YL', 'YM', 'YN', 'YP',
                                       'YQ', 'YR', 'YS', 'YT', 'YV', 'YW', 'YY']

                    for seq in sequence_file:
                        aa = seq.seq
                        aa_counts = []
                        for pairs in aa_pair:
                            count = aa.count(pairs)
                            aa_counts.append(count)
                        first = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W','Y']
                        second = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W','Y']
                        final = []

                        for i in first:
                            for j in second:
                                x = aa.count(i) * aa.count(j) / len(aa)
                                final.append(x)

                        diff = [x1 - x2 for (x1, x2) in zip(aa_counts, final)]

                        rounded = []
                        for a in diff:
                            r = round(a)
                            rounded.append(r)

                        f1 = []
                        for i in range(0, len(aa) - 2):
                            a = aa[i]
                            b = aa[i + 1]
                            c = aa[i + 2]
                            d = a + b
                            e = b + c
                            if (i == 0):
                                x1 = (aa.count(d) - round((aa.count(a) * aa.count(b)) / len(aa)))
                                f1.append(x1)
                            if (d, e in aa_pair):
                                x = (aa.count(d) - round((aa.count(a) * aa.count(b)) / len(aa))) + (aa.count(e) - round((aa.count(b) * aa.count(c)) / len(aa)))
                                f1.append(x)

                        list1 = np.arange(0, len(aa) - 1)
                        list2 = f1[0:len(aa) - 1]

                    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
                    sequence_file = list(SeqIO.parse(stringio, "fasta"))
                    for seq in sequence_file:
                        aa_dict1 = dict(zip(aa_list, max_pdp))
                        aa_dict2 = dict(zip(aa_list, fcaa))
                        csv_filename = seq.id + ".csv"
                        with open(csv_filename, "w") as csv_file:
                            writer = csv.writer(csv_file, delimiter=',')
                            writer.writerow(['Sequence', 'F1'])
                            for i in range(len(list1)):
                                writer.writerow([aa[i], list2[i]])

                        output_zip.write(csv_filename)
                        time.sleep(1)
            output_zip.close()
            st.success("Processing Complete!")
            with open("output.zip", "rb") as fp:
                btn = st.download_button(label="Download", data=fp, file_name="output.zip",mime="application/zip")
        st.write("")
        st.write("""***""")
        quote = random.choice(quotes)
        st.write(quote)
        st.write("""***""")
