import torch
import torch.nn.functional as F
from torch.nn import Sequential, Linear, ReLU, BatchNorm1d
from torch_geometric.nn import GINConv, global_add_pool

class GINEncoder(torch.nn.Module):
    def __init__(self, num_node_features, dim):
        super(GINEncoder, self).__init__()
        nn1 = Sequential(Linear(num_node_features, dim), ReLU(), Linear(dim, dim))
        self.conv1 = GINConv(nn1)
        self.bn1 = BatchNorm1d(dim)

        nn2 = Sequential(Linear(dim, dim), ReLU(), Linear(dim, dim))
        self.conv2 = GINConv(nn2)
        self.bn2 = BatchNorm1d(dim)

        nn3 = Sequential(Linear(dim, dim), ReLU(), Linear(dim, dim))
        self.conv3 = GINConv(nn3)
        self.bn3 = BatchNorm1d(dim)

    def forward(self, x, edge_index, batch):
        x = F.relu(self.conv1(x, edge_index))
        x = self.bn1(x)
        x = F.relu(self.conv2(x, edge_index))
        x = self.bn2(x)
        x = F.relu(self.conv3(x, edge_index))
        x = self.bn3(x)
        x = global_add_pool(x, batch)
        return x

class SiameseGNN(torch.nn.Module):
    def __init__(self, num_node_features, dim=64):
        super(SiameseGNN, self).__init__()
        self.encoder = GINEncoder(num_node_features, dim)

    def forward(self, graph_control, graph_screen):
        emb_control = self.encoder(graph_control.x, graph_control.edge_index, graph_control.batch)
        emb_screen = self.encoder(graph_screen.x, graph_screen.edge_index, graph_screen.batch)
        similarity = F.cosine_similarity(emb_control, emb_screen)
        mimicry_score = (similarity + 1.0) / 2.0
        return mimicry_score
