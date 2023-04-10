from __future__ import absolute_import, division, print_function
import streamlit as st
import warnings
warnings.filterwarnings("ignore")
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tensorflow as tf
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import confusion_matrix, precision_recall_curve, roc_auc_score, roc_curve, accuracy_score
from sklearn.ensemble import RandomForestClassifier
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.models import load_model
from keras import optimizers
import seaborn as sns

## Import Keras objects for Deep Learning
from keras.models  import Sequential
from keras.layers import Input, Dense, Flatten, Dropout, BatchNormalization
from keras.optimizers import Adam, SGD, RMSprop

def model():
    uploaded_file = st.file_uploader("Upload your data", type="csv")
    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file)

        x = df.drop('Target', axis=1)
        y = df['Target']

        x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=1111)

        model_1 = Sequential([
            Dense(4, input_shape=(4,), activation='tanh'),
            Dense(8, activation='tanh'),
            Dense(8, activation='tanh'),
            Dense(1, activation='sigmoid')])

        st.write("The model summary for the developed model is :: ", model_1.summary)

        #model_1.save("path_to_save_model/model.h5")

        adam = optimizers.Adam(lr=0.001)
        model_1.compile(optimizer=adam, loss='binary_crossentropy', metrics=['accuracy'])
        run_1 = model_1.fit(x_train, y_train, validation_data=(x_test, y_test), epochs=500)


        # evaluating the model on train and test set
        st.write("Train split:")
        loss, accuracy = model_1.evaluate(x_train, y_train, verbose=1)
        st.write("Accuracy : {:5.2f}".format(accuracy))

        st.write("Test split:")
        loss, accuracy = model_1.evaluate(x_test, y_test, verbose=1)
        st.write("Accuracy : {:5.2f}".format(accuracy))

        #predict_class = model_1.predict_classes(x_test)
        #st.write("Predicted Classes are :: ", predict_class)

        predict_proba = model_1.predict(x_test)
        st.write("Predicted Probability is :: ", predict_proba)