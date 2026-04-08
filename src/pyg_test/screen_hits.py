import torch
from torch_geometric.data import Data, Batch
from rdkit import Chem
import numpy as np
from siamese_gnn import SiameseGNN
import time

def featurize_mol(mol):
    """Converts an RDKit molecule into a PyTorch Geometric Graph."""
    if mol is None:
        return None
    
    # 8-dimensional atom features (Atomic Num, Degree, Charge, Aromatic, etc.)
    node_features = []
    for atom in mol.GetAtoms():
        features = [
            atom.GetAtomicNum(),
            atom.GetDegree(),
            atom.GetFormalCharge(),
            int(atom.GetIsAromatic()),
            atom.GetImplicitValence(),
            atom.GetMass(),
            atom.GetExplicitValence(),
            atom.GetTotalNumHs()
        ]
        node_features.append(features)
        
    x = torch.tensor(node_features, dtype=torch.float)
    
    # Bond connectivity (Edge Index)
    edge_indices = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        edge_indices += [[i, j], [j, i]] # Undirected graph
        
    if len(edge_indices) > 0:
        edge_index = torch.tensor(edge_indices, dtype=torch.long).t().contiguous()
    else:
        edge_index = torch.empty((2, 0), dtype=torch.long)
        
    return Data(x=x, edge_index=edge_index)

def main():
    print("🚀 Firing up the Siamese GPU Engine...")
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    try:
        model = SiameseGNN(num_node_features=8, dim=64).to(device)
        model.eval()
    except NameError:
        print("❌ ERROR: Could not find SiameseGNN. Make sure siamese_gnn.py is in this folder!")
        return

    # UPDATED: The exact, verified path to your Abaucin SDF file!
    control_path = "../../data/a_baumannii/loldf/ligand/abaucin.sdf" 
    hits_path = "../../data/drug_repurposing/repo_hits.sdf"

    print(f"🧬 Loading Abaucin Control from {control_path}...")
    
    # Using the SDMolSupplier because we know it's a valid SDF file now
    control_suppl = Chem.SDMolSupplier(control_path)
    
    try:
        control_mol = next(control_suppl)
    except StopIteration:
        print("❌ ERROR: The Abaucin SDF file is empty.")
        return
        
    if control_mol is None:
        print("❌ ERROR: Could not read the Abaucin SDF file. Please check the path!")
        return
        
    graph_control = featurize_mol(control_mol)
    batch_control = Batch.from_data_list([graph_control]).to(device)

    print(f"📂 Loading the 39,099 Screening Hits...")
    hits_suppl = Chem.SDMolSupplier(hits_path)
    
    results = []
    start_time = time.time()
    
    print("\n🔬 Screening in progress. GPU is accelerating...")
    with torch.no_grad():
        for idx, mol in enumerate(hits_suppl):
            if mol is None:
                continue
            
            # Extract molecule name or assign ID
            mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"Hit_{idx}"
            
            # Convert to Graph and Batch
            graph_hit = featurize_mol(mol)
            batch_hit = Batch.from_data_list([graph_hit]).to(device)
            
            # Calculate Mimicry!
            score = model(batch_control, batch_hit).item()
            results.append((mol_name, score))
            
            # Print progress every 5000 molecules
            if (idx + 1) % 5000 == 0:
                print(f"   ... Processed {idx + 1} candidates.")

    # Sort results from highest match to lowest
    results.sort(key=lambda x: x[1], reverse=True)
    
    # Save the Top 100
    print("\n🏆 Screening Complete! Time: {:.2f} seconds".format(time.time() - start_time))
    print("Saving top candidates to 'top_gnn_hits.csv'...")
    
    with open("top_gnn_hits.csv", "w") as f:
        f.write("Molecule_Name,Mimicry_Score\n")
        for name, score in results[:100]:
            f.write(f"{name},{score:.4f}\n")
            
    print("✅ Pipeline Complete. Check top_gnn_hits.csv for your best leads!")

if __name__ == "__main__":
    main()