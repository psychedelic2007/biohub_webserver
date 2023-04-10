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

def model_updated():
    uploaded_file = st.file_uploader("Upload your data", type="csv")
    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file)

        learning_rate = st.slider('Select learning rate', 0.0001, 0.01, step=0.001)
        num_hidden_layers = st.slider('Select number of hidden layers', 1, 5, step=1)
        num_neurons = []
        activation_funcs = []
        for i in range(num_hidden_layers):
            num_neurons.append(st.slider(f'Select number of neurons in layer {i + 1}', 1, 50, step=1))
            activation_funcs.append(st.selectbox(f'Select activation function for layer {i + 1}',
                                                         ['relu', 'sigmoid', 'tanh']))
        test_size = st.slider('Select test size', 0.1, 0.5, step=0.05)
        validation_size = st.slider('Select validation size', 0.1, 0.5, step=0.05)

        if st.button("Start Model Training"):
            x = df.drop('Target',axis=1)
            y = df['Target']
            x_train_val, x_test, y_train_val, y_test = train_test_split(x, y, test_size=test_size, random_state=1111)
            x_train, x_val, y_train, y_val = train_test_split(x_train_val, y_train_val, test_size=validation_size,
                                                          random_state=1111)

            model = Sequential()
            model.add(Dense(num_neurons[0], input_shape=(x_train.shape[1],), activation=activation_funcs[0]))
            for i in range(1, num_hidden_layers):
                model.add(Dense(num_neurons[i], activation=activation_funcs[i]))
            model.add(Dense(1, activation='sigmoid'))

            st.write("The model summary for the developed model is :: ")
            st.write(model.summary())

            adam = optimizers.Adam(lr=learning_rate)
            model.compile(optimizer=adam, loss='binary_crossentropy', metrics=['accuracy'])
            model.build(input_shape=x_train.shape)

            run = model.fit(x_train, y_train, validation_data=(x_val, y_val), epochs=500, verbose=0)

            # evaluating the model on train, validation and test set
            st.write("Train split:")
            loss, accuracy = model.evaluate(x_train, y_train, verbose=1)
            st.write("Accuracy : {:5.2f}".format(accuracy))

            st.write("Validation split:")
            loss, accuracy = model.evaluate(x_val, y_val, verbose=1)
            st.write("Accuracy : {:5.2f}".format(accuracy))

            st.write("Test split:")
            loss, accuracy = model.evaluate(x_test, y_test, verbose=1)
            st.write("Accuracy : {:5.2f}".format(accuracy))

            predict_proba = model.predict(x_test)
            st.write("Predicted Probability is :: ", predict_proba)