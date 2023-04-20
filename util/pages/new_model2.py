import pickle
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import tensorflow as tf
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, classification_report, roc_curve, roc_auc_score
from zipfile import ZipFile
import zipfile
import io
import joblib
import random

def new_model2():
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

    st.title('Model Training Web App')

    uploaded_file = st.file_uploader('Upload the data file', type='zip')

    def read_zipfile(file):
        with zipfile.ZipFile(file) as zip:
            filenames = [f for f in zip.namelist() if f.endswith('.csv')]
            dfs = []
            for filename in filenames:
                with zip.open(filename) as f:
                    dfs.append(pd.read_csv(io.StringIO(f.read().decode())))
        return dfs

    def create_zip_file(df, df_train, df_test, acc, loss, confmat, roc, model):
        zip_file = io.BytesIO()

        with zipfile.ZipFile(zip_file, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for i, df in enumerate(dfs):
                df.to_csv(f'original_data_model_{i}.csv', index=False)
                df_train.to_csv(f'training_data_model_{i}.csv', index=False)
                df_test.to_csv(f'testing_data_model_{i}.csv', index=False)
                zipf.write(f'original_data_model_{i}.csv')
                zipf.write(f'training_data_model_{i}.csv')
                zipf.write(f'testing_data_model_{i}.csv')

            for i, his in enumerate(hist):
                fig2, ax2 = plt.subplots(figsize=(18, 10))
                ax2.plot(his.history['accuracy'], linewidth=3)
                ax2.plot(his.history['val_accuracy'], linewidth=3)
                ax2.set_title('Model accuracy', fontsize=20)
                ax2.set_ylabel('Accuracy', fontsize=20)
                ax2.set_xlabel('Epoch', fontsize=20)
                ax2.tick_params(axis='both', which='major', labelsize=15)
                ax2.legend(['Training Accuracy', 'Validation Accuracy'], fontsize=20)
                plt.savefig(f"model_{i}_accuracy.pdf", dpi=300)
                zipf.write(f'model_{i}_accuracy.pdf')

                fig1, ax1 = plt.subplots(figsize=(18, 10))
                ax1.plot(his.history['loss'], linewidth=3)
                ax1.plot(his.history['val_loss'], linewidth=3)
                ax1.set_title('Model Loss', fontsize=20)
                ax1.set_ylabel('Loss', fontsize=20)
                ax1.set_xlabel('Epoch', fontsize=20)
                ax1.tick_params(axis='both', which='major', labelsize=15)
                ax1.legend(['Training Loss', 'Validation Loss'], fontsize=20)
                plt.savefig(f"model_{i}_loss.pdf", dpi=300)
                zipf.write(f'model_{i}_loss.pdf')

                fig, ax = plt.subplots(figsize=(18, 10))
                sns.heatmap(cm, annot=True, cmap='Blues', fmt='g', ax=ax)
                ax.set_xticks(np.arange(len(cm)), fontsize=15)
                ax.set_yticks(np.arange(len(cm)), fontsize=15)
                ax.set_xticklabels(['Negative', 'Positive'], rotation=45, ha='center', fontsize=20)
                ax.set_yticklabels(['Negative', 'Positive'], va='center', fontsize=20)
                ax.set_xlabel("Predicted", fontsize=20)
                ax.set_ylabel("Actual", fontsize=20)
                cbar = ax.collections[0].colorbar
                cbar.ax.tick_params(labelsize=15)
                plt.savefig(f"confusion_matrix_model_{i}.pdf", dpi=300)
                zipf.write(f'confusion_matrix_model_{i}.pdf')

                fig3, ax3 = plt.subplots(figsize=(18, 10))
                ax3.plot(fpr, tpr, label=f'AUC = {auc:.3f}')
                ax3.plot([0, 1], [0, 1], 'k--')
                ax3.set_xlabel("False Positive Rate", fontsize=20)
                ax3.set_ylabel("True Positive Rate", fontsize=20)
                ax3.set_title('ROC Curve', fontsize=20)
                ax3.tick_params(axis='both', which='major', labelsize=15)
                ax3.legend(fontsize=20)
                plt.savefig(f"roc_curve_model{i}.pdf", dpi=300)
                zipf.write(f'roc_curve_model_{i}.pdf')

            for i, model in enumerate(models):
                 model.save(f"model_{i}.h5")

        return zip_file.getvalue()

    if uploaded_file is not None:
        st.success('Data uploaded successfully!')

        dfs = read_zipfile(uploaded_file)
        st.success('Files Extracted Successfully!')

        # User selects model parameters
        learning_rate = st.slider('Select learning rate', 0.01, 0.5, 0.1)
        hidden_layers = st.slider('Select number of hidden layers', 1, 10, 3)
        neurons = st.slider('Select number of neurons in each hidden layer', 1, 100, 10)
        activation_func = st.selectbox('Select activation function', ['relu', 'sigmoid', 'tanh'])
        epochs = st.slider('Select number of epochs', 10, 1000, 50)
        train_size = st.slider('Select percentage of data for Training', 10,90,70)
        test_size = 100 - train_size

        # User clicks submit button to start training
        if st.button('Submit'):
            st.write('''***''')
            st.write('Training started...')

            # Define the model architecture
            models = []
            hist = []
            for i, df in enumerate(dfs):

                # Split the data into training and testing sets
                x = df.drop('Target', axis=1)
                x = np.asarray(X).astype(np.float32)
                y = df['Target']
                x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=test_size / 100, random_state=42)

                st.write(f'Training the model for CSV file_{i}')
                model = tf.keras.Sequential()
                model.add(tf.keras.layers.Dense(neurons, input_shape=(x_train.shape[1],), activation=activation_func))
                for i in range(hidden_layers - 1):
                    model.add(tf.keras.layers.Dense(neurons, activation=activation_func))

                model.add(tf.keras.layers.Dense(1, activation='sigmoid'))

                # Compile the model
                model.compile(loss='binary_crossentropy', optimizer=tf.keras.optimizers.Adam(lr=learning_rate),metrics=['accuracy'])

                # Train the model
                history = model.fit(x_train, y_train, epochs=epochs, validation_data=(x_test, y_test))

                # Display training progress using a progress bar
                progress_bar = st.progress(0)
                status_text = st.empty()
                for i in range(epochs):
                    model.fit(x_train, y_train, epochs=1, batch_size=32, verbose=0)
                    status_text.text(f'Epoch {i + 1}/{epochs}')
                    progress_bar.progress((i + 1) / epochs)
                status_text.text(f'Training completed')
                models.append(model)
                hist.append(history)

                st.write('''***''')

            # evaluating the model on train and test set
            st.header("Your Model Evaluation")
            st.subheader("Train split:")
            i = np.arange(1,len(models))
            for i, model in enumerate(models):
                loss, accuracy = model.evaluate(x_train, y_train, verbose=1)
                st.write("Accuracy is: {:5.2f}".format(accuracy)+ f" for model {i}")
                #st.write(f"for model{i}")

            st.subheader("Test split:")
            for i, model in enumerate(models):
                loss, accuracy = model.evaluate(x_test, y_test, verbose=1)
                st.write("Accuracy is: {:5.2f}".format(accuracy)+ f" for model {i}")

            for i, model in enumerate(models):
                predict_proba = model.predict(x_test)
                st.write(f"Predicted Probability according to model {i} is:: ", predict_proba)

            st.write('''***''')

            st.subheader('Training and validation Accuracy:')
            for i, his in enumerate(hist):
                fig2, ax2 = plt.subplots(figsize=(18,10))
                ax2.plot(his.history['accuracy'], linewidth=3)
                ax2.plot(his.history['val_accuracy'], linewidth=3)
                ax2.set_title('Model accuracy', fontsize=20)
                ax2.set_ylabel('Accuracy', fontsize=20)
                ax2.set_xlabel('Epoch', fontsize=20)
                ax2.tick_params(axis='both', which='major', labelsize=15)
                ax2.legend(['Training Accuracy', 'Validation Accuracy'], fontsize=20)
                plt.savefig(f"model_{i}_accuracy.pdf", dpi=300)
                st.subheader(f"Accuracy Plot for model {i}")
                st.pyplot(fig2)

            st.subheader('Training and validation Loss:')
            for i, his in enumerate(hist):
                fig1, ax1 = plt.subplots(figsize=(18, 10))
                ax1.plot(his.history['loss'], linewidth=3)
                ax1.plot(his.history['val_loss'], linewidth=3)
                ax1.set_title('Model Loss', fontsize=20)
                ax1.set_ylabel('Loss', fontsize=20)
                ax1.set_xlabel('Epoch', fontsize=20)
                ax1.tick_params(axis='both', which='major', labelsize=15)
                ax1.legend(['Training Loss', 'Validation Loss'], fontsize=20)
                plt.savefig(f"model_{i}_loss.pdf", dpi=300)
                st.subheader(f"Loss Plot for model {i}")
                st.pyplot(fig1)

            st.write('''***''')

            #model predictions for test set
            for i, model in enumerate(models):
                y_pred = model.predict(x_test)
                y_pred = (y_pred > 0.5)

                #confusion matrix
                st.subheader(f'Confusion Matrix for model {i}')
                cm = confusion_matrix(y_test, y_pred)
                fig, ax = plt.subplots(figsize=(18,10))
                sns.heatmap(cm, annot=True, cmap='Blues', fmt='g', ax=ax)
                ax.set_xticks(np.arange(len(cm)), fontsize=15)
                ax.set_yticks(np.arange(len(cm)), fontsize=15)
                ax.set_xticklabels(['Negative','Positive'], rotation=45, ha='center', fontsize=20)
                ax.set_yticklabels(['Negative','Positive'], va='center', fontsize=20)
                ax.set_xlabel("Predicted", fontsize=20)
                ax.set_ylabel("Actual", fontsize=20)
                cbar = ax.collections[0].colorbar
                cbar.ax.tick_params(labelsize=15)
                plt.savefig(f"confusion_matrix_model_{i}.pdf", dpi=300)
                st.pyplot(fig)

                #classification report
                st.subheader(f"Classification Report for model {i}")
                cr = classification_report(y_test, y_pred, output_dict=True)
                cr_df = pd.DataFrame(cr).transpose()
                st.write(cr_df)

            #plot ROC curve
            for i, model in enumerate(models):
                y_pred_proba = model.predict(x_test)[:, 0]
                fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
                auc = roc_auc_score(y_test, y_pred_proba)
                st.subheader(f"Receiver Operating Curve for model {i}")
                fig3, ax3 = plt.subplots(figsize=(18, 10))
                ax3.plot(fpr, tpr, label=f'AUC = {auc:.3f}')
                ax3.plot([0,1],[0,1], 'k--')
                ax3.set_xlabel("False Positive Rate", fontsize=20)
                ax3.set_ylabel("True Positive Rate", fontsize=20)
                ax3.set_title('ROC Curve', fontsize=20)
                ax3.tick_params(axis='both', which='major', labelsize=15)
                ax3.legend(fontsize=20)
                plt.savefig(f"roc_curve_model_{i}.pdf", dpi=300)
                st.pyplot(fig3)

            for i, df in enumerate(dfs):
                x = df.drop('Target', axis=1)
                y = df['Target']
                x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=test_size / 100, random_state=42)
                train_data = pd.concat([x_train, y_train], axis=1)
                test_data = pd.concat([x_test, y_test], axis=1)

            for i, model in enumerate(models):
                m = model.save(f"model_{i}.h5")
                zip_file = create_zip_file(df, train_data, test_data, fig2, fig1, fig, fig3, m)

            with st.spinner("Please wait while SAMOSA compiles your model training data......"):
                st.success("Data compiled!")
                st.subheader("Download the model output")
                st.download_button(label="Download", data=zip_file, file_name="model_output.zip", mime="application/zip")

    quote = random.choice(quotes)
    st.write(quote)
