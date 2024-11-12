# Define the callback-like function to adjust learning rate based on AUC
def adjust_learning_rate_by_auc(epoch, model, X_test, y_test_outcome, lr_dict, auc_thresholds,test_auc):
    # Get model's current learning rate
    current_lr = tf.keras.backend.get_value(model.optimizer.lr)

    # Adjust learning rate based on AUC thresholds
    new_lr = current_lr
    for threshold in auc_thresholds:
        if test_auc >= threshold:
            new_lr = lr_dict[threshold]
    
    # Set new learning rate
    if new_lr != current_lr:
        tf.keras.backend.set_value(model.optimizer.lr, new_lr)
        print(f"Epoch {epoch + 1}: Adjusting learning rate to {new_lr:.6f}")
    
    return test_auc
