import torch
from torch import nn
from config.params import Config

class Attention(nn.Module):
    def __init__(self, hidden_size):
        super().__init__()
        self.attn = nn.Linear(hidden_size, 1)

    def forward(self, x):
        weights = torch.softmax(self.attn(x), dim=1)
        output = (x * weights).sum(dim=1)
        return output.unsqueeze(1)

class Permute(nn.Module):
    def __init__(self, *order):
        super().__init__()
        self.order = order

    def forward(self, x):
        return x.permute(self.order)

class SRAMP(nn.Module):
    def __init__(self, mode='genomeonly', halfseqlen=20):
        super().__init__()
        self.mode = mode
        self.halfseqlen = halfseqlen
        self.conv_model = nn.Sequential(
            Permute(0, 2, 1),
            nn.Conv1d(41, Config.conv_len, 3, padding=1),
            nn.Conv1d(Config.conv_len, Config.conv_len, 5, padding=2),
            Permute(0, 2, 1),
        )
        self.token_emb = nn.Linear(41, Config.latent - Config.conv_len - Config.emb_len)
        self.genometrans = nn.Linear(Config.emb_len, Config.emb_len)
        self.trans = nn.TransformerEncoderLayer(Config.latent, 4, 4 * Config.latent, dropout=0.2, activation='gelu', batch_first=True, norm_first=True)
        self.lstm = nn.GRU(Config.latent, Config.latent // 2, num_layers=2, bidirectional=True, batch_first=True, dropout=0.2)
        self.attn = Attention(Config.latent)
        self.tail = nn.Sequential(
            nn.Linear(Config.latent, Config.latent // 4),
            nn.GELU(),
            nn.Linear(Config.latent // 4, 2)
        )

    def forward(self, x, emb):
        if Config.half_length > self.halfseqlen:
            x = x[:, Config.half_length - self.halfseqlen: Config.half_length + self.halfseqlen + 1]
            emb = emb[:, Config.half_length - self.halfseqlen: Config.half_length + self.halfseqlen + 1]
        te = self.token_emb(x)
        me = self.conv_model(x)
        x = torch.cat([te, me], axis=-1)
        if self.mode == 'seqonly':
            emb = emb.zero_()
        elif self.mode == 'genomeonly':
            x = x.zero_()
        y = self.genometrans(emb)
        x = torch.cat([x, y], axis=-1)
        x = self.trans(x)
        x, _ = self.lstm(x)
        x = self.attn(x)
        x = self.tail(x)
        return torch.sigmoid(x[:, :, -1])