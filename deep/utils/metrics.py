from sklearn.metrics import accuracy_score, recall_score, matthews_corrcoef, f1_score, roc_auc_score

def evaluate(y_true, y_pred):
    return {
        'accuracy': accuracy_score(y_true, y_pred > 0.5),
        'recall': recall_score(y_true, y_pred > 0.5),
        'mcc': matthews_corrcoef(y_true, y_pred > 0.5),
        'f1': f1_score(y_true, y_pred > 0.5),
        'auc': roc_auc_score(y_true, y_pred)
    }
