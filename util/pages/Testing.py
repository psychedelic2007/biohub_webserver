import streamlit as st
import pandas as pd
import numpy as np
import zipfile
import os
from sklearn.model_selection import train_test_split
import tensorflow as tf
from tensorflow import keras

def testing():
    # Define a function to extract files from zip archive
    def extract_zipfile(zip_file):
        with zipfile.ZipFile(zip_file, 'r') as zip_ref:
            zip_ref.extractall('')


    # Define a function to load and preprocess data
    def load_data(file_path):
        data = pd.read_csv(file_path)
        # Do your data preprocessing here
        return data


    # Define your model architecture and training function here
    def train_model(data):
        # Split the data into training and validation sets
        X_train, X_val, y_train, y_val = train_test_split(data.iloc[:, :-1], data.iloc[:, -1], test_size=0.2, random_state=42)

        # Define your model architecture here
        model = keras.Sequential([
            keras.layers.Dense(64, input_shape=(X_train.shape[1],), activation='relu'),
            keras.layers.Dense(64, activation='relu'),
            keras.layers.Dense(1, activation='sigmoid')
        ])

        # Compile the model
        model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

        # Train the model
        model.fit(X_train, y_train, epochs=10, batch_size=32, validation_data=(X_val, y_val))

        # Return the trained model
        return model


    # Define a function to save the model in h5 format
    def save_model(model, file_path):
        model.save(file_path)


    st.title("Training model on uploaded data")
    st.write("Please upload your training data in zip format")

    uploaded_file = st.file_uploader("Choose a file")
    if uploaded_file is not None:
        extract_zipfile(uploaded_file)
        st.write("Files extracted successfully")

        # Initialize a list to store the trained models
        trained_models = []

        # Loop through extracted files and train the model on each file
        csv_files = [f for f in os.listdir('./') if f.endswith('.csv')]
        num_csv_files = len(csv_files)
        train_files = int(0.8*num_csv_files)
        test_files = num_csv_files - train_files
        st.write(f"Number of CSV files found: {num_csv_files}")
        st.write(f"Using {train_files} files for training and {test_files} files for testing.")

        for i, file_name in enumerate(csv_files):
            if i < train_files:
                st.write(f"Training model on {file_name}...")
                data = load_data(file_name)
                model = train_model(data)
                trained_models.append(model)
                model_name = file_name.split('.')[0] + '.h5'
                save_model(model, model_name)
                st.write(f"Model trained on {file_name} and saved as {model_name}")
            else:
                st.write(f"Skipping file {file_name} for testing.")

        st.write("All models trained successfully!")

        # Add a download button for each trained model
        for i, model in enumerate(trained_models):
            model_name = csv_files[i].split('.')[0] + '.h5'
            with open(model_name, 'rb') as f:
                bytes_data = f.read()
                st.download_button(label=f"Download trained model {i+1}", data=bytes_data, file_name=model_name)

