input_shape = (X_train.shape[1],)[0]  # Number of genes

def build_outcome_classifier():
    model = keras.Sequential()
    model.add(layers.Input(shape=(input_shape,)))  # Input shape matches your data
    
    model.add(layers.Dense((512),kernel_regularizer=tf.keras.regularizers.l2(l2_reg),kernel_initializer='he_normal'))
    model.add(layers.BatchNormalization())
    model.add(layers.LeakyReLU(alpha=0.1))  # Leaky ReLU helps with vanishing gradients
    model.add(layers.Dropout(dropout_rate))
    
    model.add(layers.Dense((256),kernel_regularizer=tf.keras.regularizers.l2(l2_reg),kernel_initializer='he_normal'))
    model.add(layers.BatchNormalization())
    model.add(layers.LeakyReLU(alpha=0.1))
    model.add(layers.Dropout(dropout_rate))

    model.add(layers.Dense((128),kernel_regularizer=tf.keras.regularizers.l2(l2_reg),kernel_initializer='he_normal'))
    model.add(layers.BatchNormalization())
    model.add(layers.LeakyReLU(alpha=0.1))
    model.add(layers.Dropout(dropout_rate))
    
    model.add(layers.Dense((64),kernel_regularizer=tf.keras.regularizers.l2(l2_reg),kernel_initializer='he_normal'))
    model.add(layers.BatchNormalization())
    model.add(layers.LeakyReLU(alpha=0.1))
    model.add(layers.Dropout(dropout_rate))

    model.add(layers.Dense((32),kernel_initializer='he_normal'))
    model.add(layers.BatchNormalization())
    model.add(layers.LeakyReLU(alpha=0.1))
    
    # Output layer for binary classification with sigmoid activation
    model.add(layers.Dense(1, activation='sigmoid'))
    
    return model