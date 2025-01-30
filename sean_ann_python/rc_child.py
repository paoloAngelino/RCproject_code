class RcChild:
    """ 
    A class for instantiating a child for the purpose of a genetic algorithm for feature selection. Can either start from scratch or will depend on cross breeding and mutation,

    Attributes:
        XXX (type): XXX
        XXX (type): XXX

    """

    def __init__(self, folds_object, current_fold=0):
        """
        Initializes the RcFoldForANN object with the data from the specified fold.
        
        Args:
            folds_object (object): The object containing the folds (e.g., rcFolds).
            current_fold (int, optional): The index of the current fold to be used for training and testing. Default is 0.
        """
        self.folds_object = folds_object  # The folds object containing the training and test sets
        self.current_fold = current_fold  # The current fold index
        
        # Extracting relevant data from the folds object
        self.unique_combinations_array = self.folds_object.y_train_folds[0]['combination_tech_outcome']  # Target variable
        self.x_train = self.folds_object.x_train_folds[current_fold]  # Gene expression data for training
        self.genes_list = self.folds_object.genes_list  # List of genes for subsetting
        self.x_test = self.folds_object.x_test_folds[current_fold]  # Gene expression data for testing
        self.y_train = self.folds_object.y_train_folds[current_fold]  # Target values for training set
        self.y_test_outcome = self.folds_object.y_test_folds[current_fold]['numerical_categories_outcome']  # Numerical outcomes for testing
