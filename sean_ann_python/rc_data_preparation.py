import rpy2.robjects as ro
import rpy2.robjects.packages as rpackages
from rpy2.robjects import pandas2ri
import pandas as pd
import numpy as np
import os
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler

# Enable pandas-to-R conversion
pandas2ri.activate()

class RcDataPreparation:
    def __init__(
        self,
        data_dir='/tmp/work/RCproject/GEO_singlecellexperiment.rds',
        test_set_size=0.1,
        simplify_categories=True,
        random_seed=1
    ):
        """
        Initializes the RcDataPreparation class with specified parameters.

        Parameters:
        - data_dir (str): Path to the RDS file used for building the model.
        - test_set_size (float): Fraction of data to be used as a test set (default: 0.1).
        - simplify_categories (bool): Whether to simplify categories in the dataset (default: True).
        - random_seed (int): Seed for reproducible dataset splits (default: 1).
        """
        self.data_dir = data_dir
        self.test_set_size = test_set_size
        self.simplify_categories = simplify_categories
        self.random_seed = random_seed
        
        # Initialize placeholders for class attributes
        self.counts_df = None
        self.metadata_df = None
        self.categories_technology = None
        self.categories_outcome = None
        self.numerical_categories_technology = None
        self.numerical_categories_outcome = None
        self.unique_combinations_array = None
        self.x_train = None
        self.x_test = None
        self.y_train = None
        self.y_test = None

        # Automatically execute the data preparation pipeline
        self.retrieve_data()
        self.filter_data()
        self.simplify_the_categories()
        self.encode_labels()
        self.add_meta_data()
        self.retrieve_all_genes()
        self.scale_data()
        self.establish_test_train()

    def retrieve_data(self):
        """
        Loads the RDS file and extracts counts and metadata into Pandas DataFrames.

        Raises:
        - FileNotFoundError: If the specified RDS file does not exist.
        - ValueError: If an error occurs during data extraction or conversion.
        """
        if not os.path.exists(self.data_dir):
            raise FileNotFoundError(f"RDS file not found at path: {self.data_dir}")

        try:
            # Load the RDS file in R
            ro.r(f"sce <- readRDS('{self.data_dir}')")
            
            # Extract count data (assumed to be stored in assays)
            counts = ro.r('assay(sce, "scalelogcounts")')
            # Extract row (gene) and column (cell) names
            gene_names = list(ro.r('rownames(sce)'))
            cell_names = list(ro.r('colnames(sce)'))
            
            # Convert to a NumPy array and store as a DataFrame
            counts_np = ro.conversion.rpy2py(counts)
            self.counts_df = pd.DataFrame(counts_np, index=gene_names, columns=cell_names)
            
            # Extract metadata and convert to a Pandas DataFrame
            metadata = ro.r('as.data.frame(colData(sce))')
            self.metadata_df = ro.conversion.rpy2py(metadata)

            print("Data successfully loaded.")
        
        except Exception as e:
            raise ValueError(f"Failed to retrieve data from RDS file: {e}")

    def filter_data(self):
        """
        Filters metadata to exclude unwanted categories and synchronizes counts with metadata.
        """
        # Remove rows with 'partial' responses
        self.metadata_df = self.metadata_df[~self.metadata_df['Response'].isin(['partial'])].copy()
        # Filter counts DataFrame to match metadata indices
        self.counts_df = self.counts_df.loc[:, self.metadata_df.index.tolist()]

    def simplify_the_categories(self):
        """
        Simplifies technology categories by mapping original labels to simplified ones.
        """
        if self.simplify_categories:
            category_map = {
                'GSE133057': 'micro', 'GSE145037': 'micro', 'GSE150082': 'micro',
                'GSE190826': 'seq', 'GSE209746': 'seq', 'GSE45404_GPL1': 'micro',
                'GSE45404_GPL2': 'micro', 'GSE93375': 'micro', 'GSE94104': 'micro'
            }
            self.categories_technology = self.metadata_df['batch'].map(category_map)
        else:
            self.categories_technology = self.metadata_df['batch']

    def encode_labels(self):
        """
        Encodes categorical labels (technology and outcome) into numerical values.
        """
        label_encoder = LabelEncoder()

        # Encode technology categories
        self.numerical_categories_technology = label_encoder.fit_transform(self.categories_technology)

        # Encode response categories
        self.categories_outcome = self.metadata_df['Response']
        self.numerical_categories_outcome = label_encoder.fit_transform(self.categories_outcome)

    def add_meta_data(self):
        """
        Adds numerical labels and combinations of technology and outcomes to metadata.
        """
        self.unique_combinations_array = (
            self.numerical_categories_outcome + 
            (self.numerical_categories_technology + 1) * 2
        ) - 2
        
        self.metadata_df['numerical_categories_technology'] = self.numerical_categories_technology
        self.metadata_df['numerical_categories_outcome'] = self.numerical_categories_outcome
        self.metadata_df['combination_tech_outcome'] = self.unique_combinations_array

    def retrieve_all_genes(self):
        self.genes_list = self.counts_df.index.to_list()

    def scale_data(self):
        scaler = MinMaxScaler()
        self.counts_df = scaler.fit_transform(self.counts_df)

    def establish_test_train(self):
        """
        Splits the dataset into training and testing sets and exports training sample IDs.
        """
        self.x_train, self.x_test, self.y_train, self.y_test = train_test_split(
            self.counts_df.T, self.metadata_df, 
            test_size=self.test_set_size, random_state=self.random_seed
        )

        # Extract target labels for training and testing sets
        self.y_train_outcome = self.y_train['numerical_categories_outcome']
        self.y_test_outcome = self.y_test['numerical_categories_outcome']
        self.y_train_tech = self.y_train['numerical_categories_technology']
        self.y_test_tech = self.y_test['numerical_categories_technology']
        self.y_train_comb = self.y_train['combination_tech_outcome']
        self.y_test_comb = self.y_test['combination_tech_outcome']

        # Export training sample indices to a file
        with open('train_samples.txt', 'w') as f:
            for name in self.y_train.index.tolist():
                f.write(f"{name}\n")
