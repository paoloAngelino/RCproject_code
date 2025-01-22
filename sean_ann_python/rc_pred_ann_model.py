import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers, models

class PredAnnModel:
    def __init__(
        self,
        current_genes,
        dropout_rate=0.3,
        balance=True,
        l2_reg=-0.2,
        batch_size=16,
        num_epochs=5000,
        report_frequency=1,
        auc_threshold=0.9,
        clipnorm=2.0,
        simplify_categories=True,
        holdout_size=0.5,
        multiplier=3,
        auc_thresholds = [0.6, 0.7, 0.8, 0.85, 0.88,0.89,0.90,0.91,0.92],
        lr_dict = lr_dict = {
                    0.6:  0.005,
                    0.7:  0.001,
                    0.8:  0.0005,
                    0.85: 0.0001,
                    0.88: 0.00005,
                    0.89: 0.00001,
                    0.9:  0.000005,
                    0.91: 0.000001,
                    0.92: 0.0000005
}
    ):
        """
        Initializes the PredAnnModel with specified hyperparameters and configuration.

        Parameters:
        - current_genes (list): A non-empty list of genes to be used as model features.
        - dropout_rate (float): Dropout rate to prevent overfitting (default: 0.3).
        - balance (bool): Whether to balance technology and outcome variables during training (default: True).
        - l2_reg (float): Strength of L2 regularization (default: -0.2).
        - batch_size (int): Batch size for training (default: 16).
        - num_epochs (int): Total number of training epochs (default: 5000).
        - report_frequency (int): Frequency of reporting model metrics (AUC and Accuracy) during training (default: 1).
        - auc_threshold (float): AUC threshold for early stopping (default: 0.9).
        - clipnorm (float): Gradient clipping norm to prevent exploding gradients (default: 2.0).
        - simplify_categories (bool): Whether to simplify categories in the dataset (default: True).
        - holdout_size (float): Proportion of samples withheld during training (default: 0.5).
        - multiplier (int): Scales the number of nodes in most network layers (default: 3).
        - auc_thresold (list): auc values for the test set for which the learning rate should be adjusted
        - lr_dict (dict): Dictionary defining the learning rate based on the measured test set accuracy.

        Raises:
        - ValueError: If `current_genes` is not a non-empty list.
        """
        if not isinstance(current_genes, list) or not current_genes:
            raise ValueError("The 'current_genes' parameter must be a non-empty list of genes.")

        self.current_genes = current_genes  # List of genes provided by the user to define model features.
        self.dropout_rate = dropout_rate  # Dropout rate for regularization.
        self.balance = balance  # Balance technology and outcomes during training.
        self.l2_reg = l2_reg  # Degree of L2 regularization.
        self.batch_size = batch_size  # Batch size for training.
        self.num_epochs = num_epochs  # Total number of training epochs.
        self.report_frequency = report_frequency  # Frequency for collecting metrics during training.
        self.auc_threshold = auc_threshold  # AUC threshold for early stopping.
        self.clipnorm = clipnorm  # Gradient clipping value to prevent exploding gradients.
        self.simplify_categories = simplify_categories  # Whether to reduce data categories (e.g., microarray vs. sequencing).
        self.holdout_size = holdout_size  # Proportion of samples withheld during training.
        self.multiplier = multiplier  # Scales the number of nodes in most layers of the network.
        self.auc_thresholds = auc_thresholds  # AUC values at which the learning rate should be adjusted
        self.lr_dict =  lr_dict  # Dynamically adjusts the learning rate based on the test set accuracy
