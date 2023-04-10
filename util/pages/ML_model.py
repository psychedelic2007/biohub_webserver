#importing the relevant packages
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import joblib
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.neural_network._multilayer_perceptron import MLPClassifier
from sklearn.metrics import classification_report, roc_curve, roc_auc_score, confusion_matrix
from tensorflow.keras.callbacks import EarlyStopping
from zipfile import ZipFile
import io

def model_training():
    st.title("Welcome to the Model Training Page!")
    st.subheader("Upload your Training Data in .csv format")

    uploaded_file = st.file_uploader("Upload your csv file", type="csv")

    def filedownload(df, filename):
        csv = df.get_value().encode()
        b64 = base64.b64encode(csv).decode()
        href = f'<a href="data:file/zip;base64,{b64}" download="{filename}">Download {filename} File</a>'
        return href

    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file)

        learning_rate_options = [0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0]
        num_hidden_layer_options = [1,2,3,4]
        num_neurons_options = [4,8,16,32,64]
        optimizer_options = ['adam','sgd']
        activation_options = ['relu','tanh','sigmoid']
        train_test_split_options = ['60/20/20', '70/15/15', '80/10/10']

        learning_rate = st.selectbox("Select Learning Rate", learning_rate_options)
        num_hidden_layers = st.selectbox("Select number of Hidden Layers", num_hidden_layer_options)
        num_neurons = st.selectbox("Select number of Neurons in a Hidden Layer", num_neurons_options)
        optimizer = st.selectbox("Select Optimizer", optimizer_options)
        activation = st.selectbox("Select Activation Function", activation_options)
        train_test_split_choice = st.selectbox("Select percentage of Training Data, Testing Data and Validation Data", train_test_split_options)

        if(train_test_split_choice=='60/20/20'):
            train_test_split_percentages = [0.6,0.2,0.2]
        elif(train_test_split_choice=='70/15/15'):
            train_test_split_percentages = [0.7,0.15,0.15]
        elif(train_test_split_choice=='80/10/10'):
            train_test_split_percentages = [0.8,0.1,0.1]

        x = df.drop('Target', axis=1)
        y = df['Target']

        x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=train_test_split_percentages[1]+train_test_split_percentages[2], random_state=42)
        x_val, x_test, y_val, y_test = train_test_split(x_test, y_test, test_size=train_test_split_percentages[2]/(train_test_split_percentages[1]+train_test_split_percentages[2]), random_state=42)

        model = MLPClassifier(max_iter=5000, early_stopping=True)
        # param_grid = {'learning_rate_init': [learning_rate],
        #               'hidden_layer_sizes': [(num_neurons,)* num_hidden_layers],
        #               'solver': [optimizer],
        #               'activation': [activation]}

        #grid_search = GridSearchCV(model, param_grid)
        with st.spinner("Have some Dosa while SAMOSA trains your model......"):
            model.fit(x_train,y_train)

        #best_model = grid_search.best_estimator_
        #train_score = best_model.score(x_train, y_train)
        #test_score = best_model.score(x_test, y_test)
        #val_score = best_model.score(x_val, y_val)

        #y_pred = best_model.predict(x_test)
        cm = confusion_matrix(y_test, y_pred)
        cr = classification_report(y_test, y_pred)

        #y_pred_proba = best_model.predict_proba(x_test)[:,1]
        #fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
        #auc = roc_auc_score(y_test, y_pred_proba)

        # fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(18,10))
        #
        # ax[0].plot(best_model.loss_curve_)
        # ax[0].set_xlabel("Epoch", fontsize=20)
        # ax[0].set_ylabel("Loss", fontsize=20)
        # ax[0].set_title("Training and Validation Loss", fontsize=20)
        # ax[0].tick_params(axis='both', which='major', labelsize=15)
        #
        # ax[1].plot(best_model.score(x_train, y_train), label="Training Accuracy")
        # ax[1].plot(best_model.score(x_test, y_test), label="Testing Accuracy")
        # ax[1].plot(best_model.score(x_val, y_val), label="Validation Accuracy")
        # ax[1].set_xlabel("Epoch", fontsize=20)
        # ax[1].set_ylabel("Training, Testing and Validation Accuracy", fontsize=20)
        # ax[1].tick_params(axis='both', which='major', labelsize=15)
        # ax[1].legend()
        #
        # fig2, ax2 = plt.subplots(figsize=(18,10))
        # ax2.plot(fpr, tpr, label=f'AUC = {auc:.3f}')
        # ax2.plot([0,1],[0,1], 'k--')
        # ax2.set_xlabel("False Positive Rate", fontsize=20)
        # ax2.set_ylabel("True Positive Rate", fontsize=20)
        # ax2.set_title('ROC Curve', fontsize=20)
        # ax2.tick_params(axis='both', which='major', labelsize=15)
        # ax2.legend()

        st.subheader("Here are the best model performance metrics")
        #st.write(f"Best parameters: {grid_search.best_params_}")
        st.write(f"Training Score: {train_score:.3f}")
        st.write(f"Testing Score: {test_score:.3f}")
        st.write(f"Validation Score: {val_score:.3f}")
        st.write(f"Confusion Matrix:\n{cm}")
        st.write(f"Classification report:\n{cr}")
        st.pyplot(fig)
        st.pyplot(fig2)

        with st.spinner("Compiling your model data......."):
            buffer = io.BytesIO()
            with ZipFile(buffer, "w") as zip_file:
                zip_file.writestr('model.pkl', joblib.dump(best_model, 'model.pkl'))
                zip_file.writestr('Data.csv', data.to_csv(index=False))
                zip_file.writestr('Training_data.csv', pd.concat([x_train, y_train], axis=1).to_csv(index=False))
                zip_file.writestr('Testing_data.csv', pd.concat([x_test, y_test], axis=1).to_csv(index=False))
                zip_file.writestr('Validation_data.csv', pd.concat([x_val,y_val], axis=1).to_csv(index=False))
                fig.savefig('Accuracy_loss.pdf', dpi=300)
                fig2.savefig('ROC.pdf', dpi=300)
                zip_file.write('Accuracy_loss.pdf')
                zip_file.write('ROC.pdf')

        st.markdown(filedownload(buffer, f"{model_name}_results.zip"), unsafe_allow_html=True)
        st.success("Download Complete!")
