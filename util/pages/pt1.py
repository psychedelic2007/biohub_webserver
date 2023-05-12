import streamlit as st
import networkx as nx
from Bio import Phylo, Entrez, SeqIO
from io import StringIO
import plotly.graph_objs as go
import random
import csv
import zipfile

# Streamlit app
def pt1():
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

    def compare_sequences(seq1, seq2):
        comparison = []
        for i in range(len(seq1)):
            if seq1[i] == seq2[i]:
                comparison.append(0)
            else:
                comparison.append(1)
        return comparison
    
    def process_tree(uploaded_file):
        stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
        tree = Phylo.read(stringio, "newick")
        return tree

    def process_tree(uploaded_file, entrez_email):
        stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
        tree = Phylo.read(stringio, "newick")
        G = Phylo.to_networkx(tree)
        pos = nx.kamada_kawai_layout(G)  # Layout algorithm
        Xn = [pos[k][0] for k in G.nodes()]
        Yn = [pos[k][1] for k in G.nodes()]
        Xe = []
        Ye = []
        for e in G.edges():
            Xe += [pos[e[0]][0], pos[e[1]][0], None]
            Ye += [pos[e[0]][1], pos[e[1]][1], None]
        trace1 = go.Scatter(x=Xe, y=Ye, mode='lines', line=dict(color='rgb(210,210,210)', width=1), hoverinfo='none')
        trace2 = go.Scatter(x=Xn, y=Yn, mode='markers', name='nodes', marker=dict(symbol='circle-dot', size=18, color='#6175c1', line=dict(color='rgb(50,50,50)', width=1)), text=[node for node in G.nodes()], hoverinfo='text')
        axis = dict(showline=False, zeroline=False, showgrid=False, showticklabels=False)
        layout = go.Layout(title='Phylogenetic Tree', width=800, height=800, showlegend=False, xaxis=axis, yaxis=axis, margin=dict(l=40, r=40, b=85, t=100), hovermode='closest')
        fig = go.Figure(layout=layout)
        fig.add_trace(trace1)
        fig.add_trace(trace2)
        st.plotly_chart(fig)

    st.title("Welcome to the Preprocessing Page!")
    st.subheader("Upload your Newick Tree file")
    uploaded_file = st.file_uploader("Upload", type=["nh", "dnd", "nwk"])
    entrez_email = st.text_input("Enter your email for Entrez login: ")

    if st.button("Submit"):
        if uploaded_file is not None:
            with st.spinner("Processing..."):
                tree = process_tree(uploaded_file)
                process_tree(tree)
                Entrez.email = entrez_email
                seq_list = []

                num_clades = len(list(tree.find_clades()))

                progress_bar = st.progress(0)

                for i, clade in enumerate(tree.find_clades()):
                    if clade.name:
                        name_parts = clade.name.split('_')
                        st.write(f'Branch: {name_parts[0]} | Sequence: {name_parts[1]}')
                        handle = Entrez.efetch(db="protein", id=name_parts[1], rettype="fasta", retmode="text")
                        record = SeqIO.read(handle, 'fasta')
                        seq_list.append(record)
                    else:
                        st.write(f'Node: {clade.confidence}')

            progress_bar.progress((i + 1) / num_clades)

        st.success("Written the sequences in the list")

        zip_filename = "sequences.zip"
        with zipfile.ZipFile(zip_filename, "w") as zip_file:
            for record in seq_list:
                seq_id = record.id
                seq_branch = record.seq
                seq_nodal = seq_list[seq_list.index(record) + 1].seq if seq_list.index(record) + 1 < len(seq_list) else ""
                comparison = compare_sequences(seq_branch, seq_nodal)
                csv_filename = f"{seq_id}.csv"

                with open(csv_filename, "w") as f:
                    writer = csv.writer(f)
                    writer.writerow(["Residue Number", "Comparison"])
                    for residue_num, comparison_value in enumerate(comparison):
                        writer.writerow([residue_num + 1, comparison_value])

            zip_file.write(csv_filename)

        st.download_button(label="Download ZIP", data=zip_filename, file_name=zip_filename, mime="application/zip")


    st.write("")
    st.write("""***""")
    quote = random.choice(quotes)
    st.write(quote)
    st.write("""***""")
