import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
import streamlit as st
from streamlit.components.v1 import html
from PIL import Image
from multiapp import MultiApp
from util.pages.home_page import home_page
from util.pages.Preprocessing import preprocessing
from util.pages.entropy_try import entropy
from util.pages.Feature1 import feature1
from util.pages.Feature2 import feature2
from util.pages.Feature3 import feature3
from util.pages.new_model2 import new_model2
from util.pages.Testing import testing
from util.pages.Tutorial import tutorial
from util.pages.About import about
from util.pages.final_ml_model_training import training
from util.pages.entropy_calculation import mrf
from util.pages.phylo import phylo
from util.pages.phylo_try import phylo_try
from util.pages.phylo_try1 import phylo_try1
from util.pages.Theory_corner import theory

app = MultiApp()

app.add_app("Home Page", home_page)
app.add_app("Preprocessing", preprocessing)
app.add_app("Entropy Calculation", entropy)
app.add_app("Feature 1 Calculation", feature1)
app.add_app("Feature 2 Calculation", feature2)
app.add_app("Feature 3 Calculation", feature3)
app.add_app("Model Training", new_model2)
app.add_app("Testing", testing)
app.add_app("Tutorial", tutorial)
app.add_app("Theory Corner", theory)
app.add_app("About", about)
app.add_app("Training", training)
app.add_app("MRF", mrf)
app.add_app("Phylo", phylo)
app.add_app("Phylo Try", phylo_try)
app.add_app("Phylo Try1", phylo_try1)
app.run()
