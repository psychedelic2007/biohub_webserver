import streamlit as st
import networkx as nx
from Bio import Phylo
from io import StringIO
import plotly.graph_objs as go
from Bio import Entrez, SeqIO
import zipfile
import random
import csv

def phylo():
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

    def compare_sequences(seq1,seq2):
        comparison = []
        for i in range(len(seq1)):
            if(seq1[i]==seq2[i]):
                comparison.append(0)
            else:
                comparison.append(1)
        return comparison

    st.title("Welcome to the Preprocessing Page!")
    st.subheader("Upload your FASTA file")
    uploaded_file = st.file_uploader("Upload", type=["nh","dnd","nwk"])
    entrez_email = st.text_input("Enter your email for Entrez login: ")
    seq1_input = st.text_input("Enter the accession ID of Descendent")
    seq2_input = st.text_input("Enter the accession ID of Ancestor")
    seq_list = []

    # Initialize the current index
    if "current_index" not in st.session_state:
        st.session_state.current_index = 0

    st.write("""***""")

    if st.button("Submit"):
        with st.spinner("Processing..."):
            if uploaded_file is not None:
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

                #tree = Phylo.read(uploaded_file.decode("utf-8"),"newick")
                Entrez.email = entrez_email
                #seq_list = []

                num_clades = len(list(tree.find_clades()))

                progress_bar = st.progress(0)

                # Print the name of the branch and node sequences
                for i,clade in enumerate(tree.find_clades()):
                    if(clade.name):
                        # Assume that the node name includes sequence information separated by an underscore
                        name_parts = clade.name.split('_')
                        st.write(f'Branch: {name_parts[0]} | Sequence: {name_parts[1]}')
                        handle = Entrez.efetch(db="protein", id=name_parts[1], rettype="fasta", retmode="text")
                        record = SeqIO.read(handle, 'fasta')
                        seq_list.append(record)
                    else:
                        st.write(f'Node: {clade.confidence}')

                    #update the progress bar
                    progress_bar.progress((i+1)/num_clades)

                #st.write(seq_list)
                st.success("Written the sequences in the list")

                seq1_seq = []
                seq2_seq = []
                for record in seq_list:
                    if(seq1_input in record.id):
                        seq1 = record.seq
                        seq1_seq.append(seq1)
                for record in seq_list:
                    if(seq2_input in record.id):
                        seq2 = record.seq
                        seq2_seq.append(seq2)

                #st.write(seq_list[0].seq)
                #st.write(seq1_seq[0])
                #st.write(len(seq1_seq))
                #st.write(seq2_seq)
                #st.write(len(seq2_seq))
                #final = []
                #for i in range(len())
                final_compare = compare_sequences(seq1_seq[0], seq2_seq[0])
                st.write(final_compare)

                final_csv = seq1_input + ".csv"
                list1 = seq1_seq[0]
                list2 = final_compare

                with open(final_csv, "w") as f:
                    writer = csv.writer(f, delimiter=",")
                    writer.writerow(['Sequence','Target'])
                    for i in range(len(list1)):
                        writer.writerow([list1[i], list2[i]])

                    st.download_button(label="Download csv", data=final_csv, file_name='final_csv.csv', mime='text/csv')

        st.write("")
        st.write("""***""")
        quote = random.choice(quotes)
        st.write(quote)
        st.write("""***""")
