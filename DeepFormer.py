import torch
import torch.nn as nn


class DeepFormer(nn.Module):
    def __init__(self, seq_length, d_model, n_heads, n_layers):
        super(DeepFormer, self).__init__()
        self.embedding = nn.Embedding(seq_length, d_model)
        encoder_layer = nn.TransformerEncoderLayer(d_model=d_model, nhead=n_heads)
        self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=n_layers)
        self.fc = nn.Linear(d_model, seq_length)

    def forward(self, x):
        x = self.embedding(x)
        x = self.transformer(x)
        x = self.fc(x)
        return x


# 模型初始化
seq_length = 1024
model = DeepFormer(seq_length=seq_length, d_model=512, n_heads=8, n_layers=6)
input_seq = torch.randint(0, seq_length, (32, seq_length))
output = model(input_seq)
print(output.shape)
