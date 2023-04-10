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
from io import StringIO

def feature_calculation():
    st.title("Welcome to the Feature Calculation Page!")
    st.subheader("Enter the sequence in FASTA format")

    sequence_input = ""

    sequence = st.text_area("Sequence Input", sequence_input, height=250)
    sequence = sequence.splitlines()
    sequence = sequence[1:]
    sequence = ''.join(sequence)

    st.subheader("Or,")
    st.subheader("Upload the FASTA file")

    uploaded_file = st.file_uploader("Upload your FASTA file", type=["fasta"])
    operations = st.multiselect("Please select the relevant features that you would like to calculate",("Feature 1", "Feature 2", "Feature 3", "Feature 4"))
    show_submit = st.button("Submit")

    # Initialize the current index
    if "current_index" not in st.session_state:
        st.session_state.current_index = 0

    st.write("""***""")

    if show_submit:
        if uploaded_file is not None:
            stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
            ids = []
            for record in SeqIO.parse(stringio, "fasta"):
                x = record.id
                ids.append(x)
            st.success(f"The FASTA file has {len(ids)} sequences")
        else:
            st.header("Input")
            st.write(sequence)
            st.write("The total length of the sequence is ", len(sequence))

            st.write('''***''')

            X = Counter(sequence)
            f1_xticks = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
            f1_xlabel = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

            fig, ax = plt.subplots(figsize=(18, 10))
            ax.bar(X.keys(), X.values(), color='green')
            ax.set_xlabel("Amino Acid", fontsize=20)
            ax.set_ylabel("Frequency", fontsize=20)
            ax.set_xticks(f1_xticks,f1_xlabel)
            ax.tick_params(labelsize=15)
            fig.savefig("amino_acid_count.pdf", dpi=300)
            st.pyplot(fig)

        results = []

        for operation in operations:
            if(operation=="Feature 1"):
                #Calculating Feature 1

                aa = sequence
                aa_pair = ['AA', 'AC', 'AD', 'AE', 'AF', 'AG', 'AH', 'AI', 'AK', 'AL', 'AM', 'AN', 'AP', 'AQ', 'AR', 'AS', 'AT',
                           'AV', 'AW', 'AY',
                            'CA', 'CC', 'CD', 'CE', 'CF', 'CG', 'CH', 'CI', 'CK', 'CL', 'CM', 'CN', 'CP', 'CQ', 'CR', 'CS', 'CT',
                            'CV', 'CW', 'CY',
                            'DA', 'DC', 'DD', 'DE', 'DF', 'DG', 'DH', 'DI', 'DK', 'DL', 'DM', 'DN', 'DP', 'DQ', 'DR', 'DS', 'DT',
                            'DV', 'DW', 'DY',
                            'EA', 'EC', 'ED', 'EE', 'EF', 'EG', 'EH', 'EI', 'EK', 'EL', 'EM', 'EN', 'EP', 'EQ', 'ER', 'ES', 'ET',
                            'EV', 'EW', 'EY',
                            'FA', 'FC', 'FD', 'FE', 'FF', 'FG', 'FH', 'FI', 'FK', 'FL', 'FM', 'FN', 'FP', 'FQ', 'FR', 'FS', 'FT',
                            'FV', 'FW', 'FY',
                            'GA', 'GC', 'GD', 'GE', 'GF', 'GG', 'GH', 'GI', 'GK', 'GL', 'GM', 'GN', 'GP', 'GQ', 'GR', 'GS', 'GT',
                            'GV', 'GW', 'GY',
                            'HA', 'HC', 'HD', 'HE', 'HF', 'HG', 'HH', 'HI', 'HK', 'HL', 'HM', 'HN', 'HP', 'HQ', 'HR', 'HS', 'HT',
                            'HV', 'HW', 'HY',
                            'IA', 'IC', 'ID', 'IE', 'IF', 'IG', 'IH', 'II', 'IK', 'IL', 'IM', 'IN', 'IP', 'IQ', 'IR', 'IS', 'IT',
                            'IV', 'IW', 'IY',
                            'KA', 'KC', 'KD', 'KE', 'KF', 'KG', 'KH', 'KI', 'KK', 'KL', 'KM', 'KN', 'KP', 'KQ', 'KR', 'KS', 'KT',
                            'KV', 'KW', 'KY',
                            'LA', 'LC', 'LD', 'LE', 'LF', 'LG', 'LH', 'LI', 'LK', 'LL', 'LM', 'LN', 'LP', 'LQ', 'LR', 'LS', 'LT',
                            'LV', 'LW', 'LY',
                            'MA', 'MC', 'MD', 'ME', 'MF', 'MG', 'MH', 'MI', 'MK', 'ML', 'MM', 'MN', 'MP', 'MQ', 'MR', 'MS', 'MT',
                            'MV', 'MW', 'MY',
                            'NA', 'NC', 'ND', 'NE', 'NF', 'NG', 'NH', 'NI', 'NK', 'NL', 'NM', 'NN', 'NP', 'NQ', 'NR', 'NS', 'NT',
                            'NV', 'NW', 'NY',
                            'PA', 'PC', 'PD', 'PE', 'PF', 'PG', 'PH', 'PI', 'PK', 'PL', 'PM', 'PN', 'PP', 'PQ', 'PR', 'PS', 'PT',
                            'PV', 'PW', 'PY',
                            'QA', 'QC', 'QD', 'QE', 'QF', 'QG', 'QH', 'QI', 'QK', 'QL', 'QM', 'QN', 'QP', 'QQ', 'QR', 'QS', 'QT',
                            'QV', 'QW', 'QY',
                            'RA', 'RC', 'RD', 'RE', 'RF', 'RG', 'RH', 'RI', 'RK', 'RL', 'RM', 'RN', 'RP', 'RQ', 'RR', 'RS', 'RT',
                            'RV', 'RW', 'RY',
                            'SA', 'SC', 'SD', 'SE', 'SF', 'SG', 'SH', 'SI', 'SK', 'SL', 'SM', 'SN', 'SP', 'SQ', 'SR', 'SS', 'ST',
                            'SV', 'SW', 'SY',
                            'TA', 'TC', 'TD', 'TE', 'TF', 'TG', 'TH', 'TI', 'TK', 'TL', 'TM', 'TN', 'TP', 'TQ', 'TR', 'TS', 'TT',
                            'TV', 'TW', 'TY',
                            'VA', 'VC', 'VD', 'VE', 'VF', 'VG', 'VH', 'VI', 'VK', 'VL', 'VM', 'VN', 'VP', 'VQ', 'VR', 'VS', 'VT',
                            'VV', 'VW', 'VY',
                            'WA', 'WC', 'WD', 'WE', 'WF', 'WG', 'WH', 'WI', 'WK', 'WL', 'WM', 'WN', 'WP', 'WQ', 'WR', 'WS', 'WT',
                            'WV', 'WW', 'WY',
                            'YA', 'YC', 'YD', 'YE', 'YF', 'YG', 'YH', 'YI', 'YK', 'YL', 'YM', 'YN', 'YP', 'YQ', 'YR', 'YS', 'YT',
                            'YV', 'YW', 'YY']

                aa_counts = []

                for pairs in aa_pair:
                    count = aa.count(pairs)
                    aa_counts.append(count)

                first = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
                second = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
                final = []

                for i in first:
                    for j in second:
                        x = aa.count(i) * aa.count(j) / len(sequence)
                        final.append(x)

                diff = [x1 - x2 for (x1, x2) in zip(aa_counts, final)]

                rounded = []
                for a in diff:
                    r = round(a)
                    rounded.append(r)

                f1 = []
                for i in range(0, len(sequence) - 2):
                    a = sequence[i]
                    b = sequence[i + 1]
                    c = sequence[i + 2]
                    d = a + b
                    e = b + c

                if (i == 0):
                    x1 = (aa.count(d) - round((aa.count(a) * aa.count(b)) / len(sequence)))
                    f1.append(x1)
                if (d, e in aa_pair):
                    x = (aa.count(d) - round((aa.count(a) * aa.count(b)) / len(sequence))) + (aa.count(e) - round((aa.count(b) * aa.count(c)) / len(sequence)))
                    f1.append(x)

                list1 = np.arange(0, len(sequence) - 1)
                list2 = f1[0:len(sequence) - 1]

            elif(operation=="Featuer 2"):
                #Calculating the second feature
                aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
                aa_count = []

                for i in aa_list:
                    new = sequence.count(i)
                    aa_count.append(new)

                def factorial(n):
                    if (n == 0 or n == 1):
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

                x = np.arange(0, 20)
                fig, ax = plt.subplots(figsize=(18, 10))
                ax.bar(x, max_pdp, color='green')
                ax.set_xlabel("Nucleotide")
                ax.set_ylabel("Frequency")
                fig.savefig(".pdf", dpi=300)
                # st.pyplot(fig)

            elif(operation=="Feature 3"):
                #calculating third feature
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

                A_fut = (((12 * A) + (2 * D) + (2 * E) + (4 * G) + (4 * P) + (4 * S) + (4 * T) + (4 * V)) / (36 * len(sequence)))
                C_fut = (((2 * R) + (2 * C) + (2 * G) + (2 * F) + (4 * S) + (2 * W) + (2 * Y)) / (18 * len(sequence)))
                D_fut = (((2 * A) + (2 * N) + (2 * D) + (4 * E) + (2 * G) + (2 * H) + (2 * Y) + (2 * V)) / (18 * len(sequence)))
                E_fut = (((2 * A) + (4 * D) + (2 * E) + (2 * Q) + (2 * G) + (2 * K) + (2 * V)) / (18 * len(sequence)))
                F_fut = (((2 * C) + (2 * I) + (6 * L) + (2 * F) + (2 * S) + (2 * Y) + (2 * V)) / (18 * len(sequence)))
                G_fut = (((4 * A) + (6 * R) + (2 * D) + (2 * C) + (2 * E) + (12 * G) + (2 * S) + (1 * W) + (4 * V)) / (36 * len(sequence)))
                H_fut = (((2 * R) + (2 * N) + (2 * D) + (4 * Q) + (2 * H) + (2 * L) + (2 * P) + (2 * Y)) / (18 * len(sequence)))
                I_fut = (((1 * R) + (2 * N) + (6 * I) + (4 * L) + (1 * K) + (3 * M) + (2 * F) + (2 * S) + (3 * T) + (3 * V)) / (27 * len(sequence)))
                K_fut = (((2 * R) + (4 * N) + (2 * E) + (2 * Q) + (1 * I) + (2 * K) + (1 * M) + (2 * T)) / (18 * len(sequence)))
                L_fut = (((4 * R) + (2 * Q) + (2 * H) + (4 * I) + (18 * L) + (2 * M) + (6 * F) + (4 * P) + (2 * S) + (1 * W) + (6 * V)) / (54 * len(sequence)))
                M_fut = (((1 * R) + (3 * I) + (2 * L) + (1 * K) + (1 * T) + (1 * V)) / (9 * len(sequence)))
                N_fut = (((2 * N) + (2 * D) + (2 * H) + (2 * I) + (4 * K) + (2 * S) + (2 * T) + (2 * Y)) / (18 * len(sequence)))
                P_fut = (((4 * A) + (4 * D) + (2 * Q) + (2 * H) + (4 * L) + (12 * P) + (4 * S) + (4 * T)) / (36 * len(sequence)))
                Q_fut = (((2 * R) + (2 * E) + (2 * Q) + (4 * H) + (2 * L) + (2 * K) + (2 * P)) / (18 * len(sequence)))
                R_fut = (((18 * R) + (2 * C) + (2 * Q) + (6 * G) + (2 * H) + (1 * I) + (4 * L) + (2 * K) + (1 * M) + (4 * P) + (6 * S) + (2 * T) + (2 * W)) / (54 * len(sequence)))
                S_fut = (((4 * A) + (6 * R) + (2 * N) + (4 * C) + (2 * G) + (2 * I) + (2 * L) + (2 * F) + (4 * P) + (14 * S) + (6 * T) + (1 * W) + (2 * Y)) / (54 * len(sequence)))
                T_fut = (((4 * A) + (2 * R) + (2 * N) + (3 * I) + (2 * K) + (1 * M) + (4 * P) + (6 * S) + (12 * T)) / (36 * len(sequence)))
                V_fut = (((4 * A) + (2 * D) + (2 * E) + (4 * G) + (3 * I) + (6 * L) + (1 * M) + (2 * F) + (12 * V)) / (36 * len(sequence)))
                W_fut = (((2 * R) + (2 * C) + (1 * G) + (1 * L) + (1 * S)) / (9 * len(sequence)))
                Y_fut = (((2 * N) + (2 * D) + (2 * C) + (2 * H) + (2 * F) + (2 * S) + (2 * Y)) / (18 * len(sequence)))

                aa_curr_count_list = [A_cur, C_cur, D_cur, E_cur, F_cur, G_cur, H_cur, I_cur, K_cur, L_cur, M_cur, N_cur, P_cur,Q_cur, R_cur, S_cur, T_cur, V_cur, W_cur, Y_cur]
                aa_fut_count_list = [A_fut, C_fut, D_fut, E_fut, F_fut, G_fut, H_fut, I_fut, K_fut, L_fut, M_fut, N_fut, P_fut,Q_fut, R_fut, S_fut, T_fut, V_fut, W_fut, Y_fut]

                fcaa = []
                for i in range(len(aa_curr_count_list)):
                    if (aa_curr_count_list[i] == 0):
                        ratio = 0
                        fcaa.append(ratio)
                    else:
                        ratio = aa_fut_count_list[i] / aa_curr_count_list[i]
                        fcaa.append(ratio)

                x = np.arange(0, 20)
                xticks = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
                xlabel = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
                fig, ax = plt.subplots(figsize=(18, 10))
                ax.bar(x, aa_curr_count_list, label="Current Count")
                ax.bar(x, aa_fut_count_list, width=0.5 * 0.8, label="Future Count")
                ax.set_xlabel("Amino Acids", fontsize=20)
                ax.set_ylabel("3rd Frequency", fontsize=20)
                ax.legend(fontsize=15)
                ax.tick_params(labelsize=15)
                ax.set_xticks(xticks, xlabel)
                st.pyplot(fig)

            else:
                #calculating feature 4
                aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
                length = len(sequence)

                aa_prob = []
                for i in aa_list:
                    new = sequence.count(i)
                    prob = float(new / length)
                    aa_prob.append(prob)

                final = []
                for i in aa_prob:
                    if (i == 0):
                        ent = 0
                        final.append(ent)
                    else:
                        ent = (-i) * (math.log(i, 2))
                        final.append(ent)

                x = np.arange(0, 20)
                fig, ax = plt.subplots(figsize=(18, 10))
                ax.bar(x, final, color='green')
                ax.set_xlabel("Nucleotide")
                ax.set_ylabel("Frequency")
                st.pyplot(fig)

                st.write("""***""")
                st.write("### Your Data file Looks Something Like This")

                df = pd.read_csv(data)
                s = df.shape
                st.dataframe(df)
