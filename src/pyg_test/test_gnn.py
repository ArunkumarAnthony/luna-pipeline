import torch
from torch_geometric.data import Data, Batch
from siamese_gnn import SiameseGNN

def create_dummy_graph():
    """Creates a fake molecular interaction graph (e.g., 5 atoms, 4 bonds)"""
    # 5 nodes (atoms/features), each with an 8-dimensional feature vector
    x = torch.rand((5, 8)) 
    
    # 4 edges (bonds/interactions) connecting the 5 nodes
    edge_index = torch.tensor([
        [0, 1, 1, 2, 3, 4], 
        [1, 0, 2, 1, 4, 3]
    ], dtype=torch.long)
    
    return Data(x=x, edge_index=edge_index)

def main():
    print("Initializing Siamese GNN Engine...")
    
    # Check if GPU is available (otherwise fallback to your CPU)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Compute Device: {device}")

    # Initialize the model (8 input features, 64-dimensional embedding)
    model = SiameseGNN(num_node_features=8, dim=64).to(device)
    model.eval() # Set to evaluation mode

    # Create two fake LUNA graphs
    graph_control = create_dummy_graph()
    graph_hit = create_dummy_graph()

    # Batch them (PyG requires batching even for single graphs)
    batch_control = Batch.from_data_list([graph_control]).to(device)
    batch_hit = Batch.from_data_list([graph_hit]).to(device)

    print("\n Feeding graphs into the Siamese GNN...")
    
    # Run the forward pass!
    with torch.no_grad():
        mimicry_score = model(batch_control, batch_hit)

    print(f"Success! Generated Mimicry Score: {mimicry_score.item():.4f}")
    print("(Score is between 0.0 and 1.0. 1.0 means perfect thermodynamic match)")

if __name__ == "__main__":
    main()