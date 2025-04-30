import numpy as np
import torch
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from torch.utils.data import TensorDataset, DataLoader

def load_data(index):
    x = torch.from_numpy(np.load(f'../score_matrices_{index}.npy')).float()
    emb = np.load(f'../phenotype_features_{index}.npy')
    labels = torch.from_numpy(np.load(f'../labels_{index}.npy')).float().unsqueeze(-1)

    emb = StandardScaler().fit_transform(emb)
    emb = np.repeat(emb[:, np.newaxis, :], 41, axis=1)
    emb = torch.from_numpy(emb).float()

    return train_test_split(x, emb, labels, test_size=0.2, random_state=42)

def get_loader(x_train, emb_train, y_train, batch_size):
    dataset = TensorDataset(x_train, emb_train, y_train)
    return DataLoader(dataset, batch_size=batch_size, shuffle=True)