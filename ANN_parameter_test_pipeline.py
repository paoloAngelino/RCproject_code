import rpy2.robjects as ro
import rpy2.robjects.packages as rpackages
import anndata as ad
import scanpy as sc
import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers, models
import pandas as pd
import umap
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from keras.callbacks import EarlyStopping
import seaborn as sb
from scipy.spatial.distance import pdist
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import MinMaxScaler
from scipy.stats import pearsonr
import time
from keras import layers, models
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc, roc_auc_score
from sklearn.decomposition import PCA
import math
import random
import rpy2.robjects as ro
from rpy2 import robjects
from rpy2.robjects import r
from rpy2.robjects import pandas2ri
import pandas as pd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.metrics import roc_curve, auc
from sklearn.ensemble import RandomForestClassifier
import os

rds_path = '/tmp/work/RCproject/gene_lists.rds'
rds_data = r.readRDS(rds_path)

# Extract the names of the lists and their contents
gene_lists = {}
for name, item in zip(rds_data.names, rds_data):
    # Each 'item' is a list associated with the 'name'
    inner_list = list(item)  # Convert the inner R list to a Python list
    gene_lists[name] = inner_list

# Now `python_data` is a dictionary with names as keys and lists as values
gene_list_names = gene_lists.keys()

print('test')

dropout_rates = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]

def test_params(current_value):

    print('lets go')
    
    epoch_count = 0
    pandas2ri.activate()
    test_set_size = 0.1
    dropout_rate = current_value #0.1
    balance = True
    l2_reg = 0.1
    batch_size = 16  #determines how many samples are processed per batch, each epoch will process multiple batches
    learning_rate = 0.00001
    num_epochs = 3000
    report_frequency = 5
    accuracy_threshold = 0.95
    clipnorm = 2.0
    simplifly_categories = True
    holdout_size = 0.5
    use_gene_list = True
    current_gene_list = 'de_intersect_plus_bulk_genes'
    PCA_reduce = False
    n_comp_PCA = 16

    # Load the RDS file
    rds_path = '/tmp/work/RCproject/gene_lists.rds'
    rds_data = r.readRDS(rds_path)
    
    # Extract the names of the lists and their contents
    gene_lists = {}
    for name, item in zip(rds_data.names, rds_data):
        # Each 'item' is a list associated with the 'name'
        inner_list = list(item)  # Convert the inner R list to a Python list
        gene_lists[name] = inner_list

    current_genes = gene_lists[current_gene_list]

    # Function to read RDS file and extract counts and metadata
    def read_rds_to_matrix_and_metadata(file_path):
        # Load the RDS file in R
        ro.r(f"sce <- readRDS('{file_path}')")
    
        # Extract count data (assumed to be stored in assays)
        counts = ro.r('assay(sce, "scalelogcounts")')
        # Extract row (gene) and column (cell) names
        gene_names = ro.r('rownames(sce)')
        cell_names = ro.r('colnames(sce)')
        
        # Convert to a NumPy array
        counts_np = ro.conversion.rpy2py(counts)
    
        # Convert the counts matrix to a pandas DataFrame
        counts_df = pd.DataFrame(counts_np, index=gene_names, columns=cell_names)
    
        # Extract metadata from colData and convert to a pandas DataFrame directly
        metadata = ro.r('as.data.frame(colData(sce))')  # Get the colData as an R data frame
        metadata_df = pd.DataFrame(metadata)  # Convert R data frame to pandas DataFrame directly
    
        return counts_df, metadata_df
    
    # Usage example
    file_path = '/tmp/work/RCproject/GEO_singlecellexperiment.rds'
    counts_df, metadata_df = read_rds_to_matrix_and_metadata(file_path)
    
    if use_gene_list:
        counts_df = counts_df.loc[current_genes]

    filtered_metadata_df = metadata_df[~metadata_df['Response'].isin(['partial'])].copy()

    row_names = filtered_metadata_df.index.tolist()  # Convert index to a list
    
    filtered_counts_df = counts_df[row_names].copy()

    categories_technology = filtered_metadata_df['batch']

    if simplifly_categories:
        category_map = {'GSE133057': 'micro', 'GSE145037': 'micro', 'GSE150082': 'micro','GSE190826':'seq','GSE209746':'seq',
                        'GSE45404_GPL1': 'micro', 'GSE45404_GPL2': 'micro', 'GSE93375': 'micro','GSE94104': 'micro'}
        categories_technology = np.vectorize(category_map.get)(categories_technology)
        
    # Create a LabelEncoder instance
    label_encoder = LabelEncoder()
    # Fit and transform the categories to integers
    numerical_categories_technology = label_encoder.fit_transform(categories_technology)        
    # categories_outcome = adata.obs['Response']
    categories_outcome = filtered_metadata_df['Response']
    numerical_categories_outcome = label_encoder.fit_transform(categories_outcome)


    frequency_counts = pd.Series(numerical_categories_technology).value_counts()
    frequency_counts = pd.Series(numerical_categories_outcome).value_counts()    
    unique_combinations_array = (numerical_categories_outcome + (numerical_categories_technology+1)*2)-2
    np.unique(unique_combinations_array)

    # Assuming filtered_metadata_df is your filtered DataFrame
    filtered_metadata_df.loc[:, 'numerical_categories_technology'] = numerical_categories_technology
    filtered_metadata_df.loc[:, 'numerical_categories_outcome'] = numerical_categories_outcome
    filtered_metadata_df.loc[:, 'combination_tech_outcome'] = unique_combinations_array

    #normalizaiton
    # gene_expression_data = adata.layers['scalelogcounts']
    
    gene_expression_data = filtered_counts_df.T
    
    #Min-max normalization
    
    scaler = MinMaxScaler()
    gene_expression_data = scaler.fit_transform(gene_expression_data)
    
    number_genes = gene_expression_data.shape[1]
    input_dim = number_genes
    
    if PCA_reduce:
    # Initialize PCA and fit it to X_train
        n_components = n_comp_PCA  # You can adjust this based on your data
        pca = PCA(n_components=n_components)
        gene_expression_data = pca.fit_transform(gene_expression_data)


    # setup the test and train datasets
    X_train, X_test, y_train, y_test = train_test_split(gene_expression_data, filtered_metadata_df, test_size=test_set_size, random_state=1)
    
    y_train_outcome = y_train['numerical_categories_outcome']
    y_test_outcome = y_test['numerical_categories_outcome']
    
    y_train_tech = y_train['numerical_categories_technology']
    y_test_tech = y_test['numerical_categories_technology']
    
    y_train_comb = y_train['combination_tech_outcome']
    y_test_comb = y_test['combination_tech_outcome']

    # # Define the input shape
    # # input_shape = (gene_expression_data.shape[1],)[0]  # Number of genes

    input_shape = (X_train.shape[1],)[0]  # Number of genes
    
    def build_outcome_classifier():
        model = keras.Sequential()
        model.add(layers.Input(shape=(input_shape,)))  # Input shape matches your data
        
        model.add(layers.Dense((512),kernel_regularizer=tf.keras.regularizers.l2(l2_reg),kernel_initializer='he_normal'))
        model.add(layers.LeakyReLU(alpha=0.1))  # Leaky ReLU helps with vanishing gradients
        model.add(layers.BatchNormalization())
        model.add(layers.Dropout(dropout_rate))
        
        model.add(layers.Dense((256),kernel_regularizer=tf.keras.regularizers.l2(l2_reg),kernel_initializer='he_normal'))
        model.add(layers.LeakyReLU(alpha=0.1))
        model.add(layers.BatchNormalization())
        model.add(layers.Dropout(dropout_rate))
    
        model.add(layers.Dense((128),kernel_regularizer=tf.keras.regularizers.l2(l2_reg),kernel_initializer='he_normal'))
        model.add(layers.LeakyReLU(alpha=0.1))
        model.add(layers.BatchNormalization())
        model.add(layers.Dropout(dropout_rate))
        
        model.add(layers.Dense((64),kernel_regularizer=tf.keras.regularizers.l2(l2_reg),kernel_initializer='he_normal'))
        model.add(layers.LeakyReLU(alpha=0.1))
        model.add(layers.BatchNormalization())
        model.add(layers.Dropout(dropout_rate))
    
        model.add(layers.Dense((32),kernel_initializer='he_normal'))
        model.add(layers.LeakyReLU(alpha=0.1))
        model.add(layers.BatchNormalization())
        
        # Output layer for binary classification with sigmoid activation
        model.add(layers.Dense(1, activation='sigmoid'))
        
        return model

    # Define the mode
    outcome_classifier = build_outcome_classifier()
    
    # set up the optimizer
    optimizer = tf.keras.optimizers.Adam(learning_rate=learning_rate, clipnorm = clipnorm)
    
    # determine the sample and batch size
    num_samples = math.floor(X_train.shape[0]* (1-holdout_size))  # number of samples used in each training epoch
    
    # batch_size = adata.shape[0]
    
    # Calculate the number of steps per epoch
    num_steps_per_epoch = num_samples // batch_size
    
    # Compile the outcome discriminator
    outcome_classifier.compile(optimizer=optimizer, loss='binary_crossentropy', metrics=['accuracy'])


    test_accuracy_list = []
    train_accuracy_list = []
    
    num_outcomes = len(np.unique(y_test_outcome))
    num_conditions = len(np.unique(unique_combinations_array))
    
    # Define the training step for only the outcome classifier
    def train_step(data, outcome_labels):
        with tf.GradientTape() as tape:
            # Forward pass through the outcome classifier
            outcome_predictions = outcome_classifier(data)
    
            # Compute the biological discriminator loss
            outcome_loss = tf.keras.losses.binary_crossentropy(outcome_labels, outcome_predictions)
            outcome_loss = tf.reduce_mean(outcome_loss)  # Average over the batch
    
        # Compute gradients for the outcome classifier
        classifier_grads = tape.gradient(outcome_loss, outcome_classifier.trainable_variables)
        
        # Calculate accuracy for the outcome classifier
        predicted_outcome_labels = tf.cast(outcome_predictions > 0.5, tf.float32)  # Threshold at 0.5
        outcome_labels_float = tf.cast(outcome_labels, tf.float32)
    
        # Calculate accuracy
        accuracy = tf.reduce_mean(tf.cast(tf.equal(predicted_outcome_labels, outcome_labels_float), tf.float32))
    
        return outcome_loss, accuracy, classifier_grads
    
    # Training loop
    for epoch in range(num_epochs):
        total_loss = 0.0  # To accumulate losses
        total_accuracy = 0.0  # To accumulate accuracy
        accumulated_grads = [tf.zeros_like(var) for var in outcome_classifier.trainable_variables]  # Initialize gradient accumulator
    
        # Split train data randomly, holding out a portion for generalization
        X_train_temp, X_test_temp, y_train_temp, y_test_temp = train_test_split(X_train, y_train, test_size=holdout_size, random_state=None)
        y_train_comb_temp = y_train_temp['combination_tech_outcome']
        y_train_temp = y_train_temp['numerical_categories_outcome']
        
        # Mini-batch training loop
        for step in range(num_steps_per_epoch):
            # Balance batches if necessary
            batch_indices = []
            if balance:
                for condition in range(num_conditions):
                    condition_indices = np.where(y_train_comb_temp == condition)[0]
                    condition_batch_indices = np.random.choice(condition_indices, size=batch_size // num_conditions, replace=True)
                    batch_indices.append(condition_batch_indices)
            else:
                all_indices = np.arange(len(X_train_temp))
                random_indices = np.random.choice(all_indices, size=batch_size, replace=True)
                batch_indices.append(random_indices)
            X_batch = X_train_temp[np.concatenate(batch_indices)]
            y_batch = y_train_temp[np.concatenate(batch_indices)]
            y_batch = tf.expand_dims(y_batch, axis=-1)  # Adjust labels shape for binary_crossentropy
                    
            # Perform the training step and collect gradients
            outcome_loss, accuracy, classifier_grads = train_step(X_batch, y_batch)
            
            # Accumulate gradients and losses
            total_loss += outcome_loss.numpy()
            total_accuracy += accuracy.numpy()
            accumulated_grads = [acc_grad + grad for acc_grad, grad in zip(accumulated_grads, classifier_grads)]
    
        # Average the accumulated gradients
        averaged_grads = [grad / num_steps_per_epoch for grad in accumulated_grads]
    
        # Apply averaged gradients to update model weights
        optimizer.apply_gradients(zip(averaged_grads, outcome_classifier.trainable_variables))
    
        # Calculate average loss and accuracy for the epoch
        avg_loss = total_loss / num_steps_per_epoch
        avg_accuracy = total_accuracy / num_steps_per_epoch
    
        # Print average accuracy at the end of each epoch and calculate the accuracy for the test set
        if epoch % report_frequency == 0:
            # Evaluate on test data
            outcome_predictions = outcome_classifier(X_test)
            predicted_outcome_labels = tf.cast(outcome_predictions > 0.5, tf.float32)  # Threshold at 0.5
            outcome_labels = tf.expand_dims(y_test_outcome, axis=-1)  # Reshape to match logits shape
            outcome_labels_float = tf.cast(outcome_labels, tf.float32)
    
            # Calculate accuracy
            test_accuracy = tf.reduce_mean(tf.cast(tf.equal(predicted_outcome_labels, outcome_labels_float), tf.float32))
            
            # Store and print metrics
            train_accuracy_list.append(avg_accuracy)
            test_accuracy_list.append(test_accuracy)
            print(f'Epoch {epoch}, Average Outcome Loss: {avg_loss}, Average Accuracy: {avg_accuracy}, Test Accuracy: {test_accuracy}')
    
            # Early stopping condition for accuracy
            if test_accuracy > accuracy_threshold:
                print('Early stopping: test set performance high enough')
                break

    current_val_str = str(current_value)
    title = 'dropout_is_' + current_val_str
    
    # Create the PDF file name based on a variable
    output_dir = "tmp_plots"
    pdf_filename = os.path.join(output_dir, f"{title}.pdf")
    
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)  # Ensure the directory is created
    
    
    # Create a new figure to hold all plots
    with PdfPages(pdf_filename) as pdf:
        fig, axs = plt.subplots(2, 2, figsize=(12, 10))  # 2 rows, 2 columns
    
        # Training accuracy plot
        x_values = np.arange(1, len(train_accuracy_list) + 1) * report_frequency
        frequency_counts = pd.Series(y_train_outcome).value_counts()
        train_chance_level = frequency_counts[0] / len(y_train_outcome)
    
        axs[0, 0].plot(x_values, train_accuracy_list, label='Training Accuracy', color='blue')
        axs[0, 0].axhline(train_chance_level, color='black', linestyle='--')
        axs[0, 0].set_title('Training Set Accuracy Over Epochs')
        axs[0, 0].set_xlabel('Epoch')
        axs[0, 0].set_ylabel('Training Accuracy')
        axs[0, 0].grid()
        axs[0, 0].legend()
    
        # Testing accuracy plot
        frequency_counts = pd.Series(y_test_outcome).value_counts()
        test_chance_level = frequency_counts[0] / len(y_test_outcome)
    
        axs[0, 1].plot(x_values, test_accuracy_list, label='Test Accuracy', color='orange')
        axs[0, 1].axhline(test_chance_level, color='black', linestyle='--')
        axs[0, 1].set_title('Test Set Accuracy Over Epochs')
        axs[0, 1].set_xlabel('Epoch')
        axs[0, 1].set_ylabel('Test Accuracy')
        axs[0, 1].grid()
        axs[0, 1].legend()
    
        # ROC curve for ANN
        fpr, tpr, thresholds = roc_curve(y_test_outcome, predicted_outcome_labels)
        roc_auc = auc(fpr, tpr)
    
        axs[1, 0].plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (area = {roc_auc:.2f})')
        axs[1, 0].plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
        axs[1, 0].set_xlim([0.0, 1.0])
        axs[1, 0].set_ylim([0.0, 1.05])
        axs[1, 0].set_xlabel('False Positive Rate')
        axs[1, 0].set_ylabel('True Positive Rate')
        axs[1, 0].set_title('ROC Curve for ANN')
        axs[1, 0].legend(loc="lower right")
    
        # ROC curve for Random Forest
        clf = RandomForestClassifier(n_estimators=100, random_state=42)
        clf.fit(X_train, y_train_outcome)
        y_pred_prob = clf.predict_proba(X_test)[:, 1]
        fpr, tpr, thresholds = roc_curve(y_test_outcome, y_pred_prob)
        roc_auc = auc(fpr, tpr)
    
        axs[1, 1].plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (area = {roc_auc:.2f})')
        axs[1, 1].plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
        axs[1, 1].set_xlim([0.0, 1.0])
        axs[1, 1].set_ylim([0.0, 1.05])
        axs[1, 1].set_xlabel('False Positive Rate')
        axs[1, 1].set_ylabel('True Positive Rate')
        axs[1, 1].set_title('ROC Curve for Random Forest')
        axs[1, 1].legend(loc="lower right")
    
        plt.tight_layout()
        pdf.savefig(fig)  # Save the figure to the PDF
        plt.close(fig)     # Close the figure

for current_value in dropout_rates:
    test_params(current_value)