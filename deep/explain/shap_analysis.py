import shap
import numpy as np
import torch
from config.params import Config
import matplotlib.pyplot as plt
from models.sramp import SRAMP
from utils.data_loader import load_data

feature_names = [
    'RBP_Num', 'SplicingSite_Num', 'ORF', 'ORF_Count', 'm6a_in_ORF',
    'm6a_in_IRES', 'HNRNPC', 'YTHDC1', 'YTHDF1', 'YTHDF2', 'METTL3',
    'METTL14', 'METTL16', 'WTAP', 'ALKBH5', 'FTO', 'miRNA',
    'Dist to ORF', 'Dist to IRES'
]

for i in range(1, 11):
    print(f"🔁 Processing round {i}...")
    x_train, x_test, emb_train, emb_test, y_train, y_test = load_data(i)
    model = SRAMP(mode='genomeonly')
    model.eval()
    model.load_state_dict(torch.load(f'model_run{i}.pt'))  # 需保存模型

    def model_predict_only_emb(emb_array_np):
        emb_tensor = torch.from_numpy(emb_array_np).float()
        emb_tensor = emb_tensor.unsqueeze(1).repeat(1, 41, 1)
        dummy_x = torch.zeros(emb_tensor.shape[0], 41, 41)
        with torch.no_grad():
            preds = model(dummy_x, emb_tensor)
        return preds.numpy()

    background = emb_train[:100, Config.half_length, :].numpy()
    test = emb_test[:300, Config.half_length, :].numpy()
    masker = shap.maskers.Independent(background)
    explainer = shap.Explainer(model_predict_only_emb, masker)
    shap_values = explainer(test)
    importance = np.mean(np.abs(shap_values.values), axis=0)
    top10_idx = np.argsort(importance)[::-1][:10]

    # Summary plot
    plt.figure()
    shap.summary_plot(
        shap_values.values[:, top10_idx],
        features=test[:, top10_idx],
        feature_names=np.array(feature_names)[top10_idx],
        show=False
    )
    plt.tight_layout()
    plt.savefig(f"shap1_run{i}.pdf")
    plt.close()

    np.save(f'shap_values_run{i}.npy', shap_values.values)

    # Bar plot
    plt.figure()
    plt.barh(range(10), importance[top10_idx], color='skyblue')
    plt.yticks(range(10), np.array(feature_names)[top10_idx])
    plt.xlabel("Mean |SHAP value|")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(f"shap_bar1_run{i}.pdf")
    plt.close()