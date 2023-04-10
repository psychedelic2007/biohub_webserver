import streamlit as st
import pandas as pd
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
import zipfile
import glob


def training():
    # Define function to train model on each CSV
    def train_model(csv_file, hidden_layer_sizes, activation, train_size):
        # Read CSV file
        df = pd.read_csv(csv_file)

        # Separate features and target
        X = df.iloc[:, :-1]
        y = df.iloc[:, -1]

        # Split data into training and testing sets
        X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=train_size)

        # Create MLPClassifier object with specified parameters
        model = MLPClassifier(hidden_layer_sizes=hidden_layer_sizes, activation=activation)

        # Fit the model to the training data
        model.fit(X_train, y_train)

        # Predict the classes of the testing data
        y_pred = model.predict(X_test)

        # Calculate the accuracy of the model
        accuracy = accuracy_score(y_test, y_pred)

        # Return the trained model and its accuracy
        return model, accuracy


    # Define function to create zip file of trained models
    def create_zip_file(models):
        # Create a new ZipFile object
        zip_file = zipfile.ZipFile('trained_models.zip', 'w')

        # Add each trained model to the zip file
        for i, model in enumerate(models):
            zip_file.writestr(f'model_{i}.pkl', pickle.dumps(model))

        # Close the zip file
        zip_file.close()

    # Set page header
    st.title('MLPClassifier Trainer')

    # Allow user to upload zip file
    zip_file = st.file_uploader('Upload zip file of CSVs', type='zip')

    # If zip file is uploaded
    if zip_file:

        # Extract zip file
        with zipfile.ZipFile(zip_file) as zf:
            zf.extractall('data')

        # Get list of CSV files
        csv_files = sorted(glob.glob('data/*.csv'))

        # Allow user to select model parameters
        hidden_layer_sizes = st.slider('Select number of neurons for each hidden layer', 1, 10, 1)
        activation = st.multiselect('Select activation function for each hidden layer',
                                    ['logistic', 'tanh', 'identity', 'relu'])
        train_size = st.slider('Select percentage of data to use for training', min_value=10, max_value=90, step=10)

        # Train model on each CSV file with specified parameters
        models = []
        accuracies = []
        for i, csv_file in enumerate(csv_files):
            model, accuracy = train_model(csv_file, hidden_layer_sizes=hidden_layer_sizes, activation=activation,
                                          train_size=train_size / 100)
            models.append(model)
            accuracies.append(accuracy)
            st.write(f'Model {i + 1} trained with accuracy {accuracy:.2f}')

        # Create zip file of trained models
        create_zip_file(models)

        # Allow user to download zip file of trained models
        st.download_button('Download trained models', data=open('trained_models.zip', 'rb').read(),
                           file_name='trained_models.zip', mime='application/zip')
