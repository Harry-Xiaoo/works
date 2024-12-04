import torch
import torch.nn as nn
from torch.nn.functional import cosine_similarity


class SimCLR(nn.Module):
    def __init__(self, encoder, projection_dim):
        super(SimCLR, self).__init__()
        self.encoder = encoder
        self.projection = nn.Linear(encoder.output_dim, projection_dim)

    def forward(self, x):
        h = self.encoder(x)
        z = self.projection(h)
        return z


# Contrastive Loss
def contrastive_loss(z_i, z_j, tau=0.5):
    sim = cosine_similarity(z_i, z_j)
    exp_sim = torch.exp(sim / tau)
    return -torch.log(exp_sim / exp_sim.sum())


# 示例
encoder = DeepFormer(seq_length=128, d_model=256, n_heads=4, n_layers=4)
model = SimCLR(encoder, projection_dim=128)
