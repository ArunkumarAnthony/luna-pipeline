#!/bin/bash

# ==============================================================================
# LUNA AUTOMATION SCRIPT: DUDE-Z MAPK2 BENCHMARK (V3 - Header-Safe)
# ==============================================================================
set -e

# --- 1. CONFIGURATION ---
PROTEIN_IN="inputs/MAPK2/rec.crg.pdb"
ACTIVES_IN="inputs/MAPK2/dudez_1pt0LD_ligand_poses.mol2"
DECOYS_IN="inputs/MAPK2/dudez_1pt0LD_decoy_poses.mol2"

RESULTS_DIR="./results"
TOY_MOL2="toy_library_poses.mol2"
ENTRIES_FILE="entries.txt"

NUM_ACTIVES=100
NUM_DECOYS=1900

PROTEIN_ID=$(basename "$PROTEIN_IN" .crg.pdb)

echo "Starting Benchmark Setup for $PROTEIN_ID..."

# --- 2. CLEANUP ---
rm -rf "$RESULTS_DIR" "$TOY_MOL2" "$ENTRIES_FILE"
mkdir -p "$RESULTS_DIR/pdbs"
mkdir -p "$RESULTS_DIR/configs"

# --- 3. PROTEIN PREP ---
echo "Fixing PDB atom naming (ILE CD -> CD1)..."
cp "$PROTEIN_IN" "$RESULTS_DIR/pdbs/${PROTEIN_ID}.pdb"
sed -i 's/ CD  ILE/ CD1 ILE/g' "$RESULTS_DIR/pdbs/${PROTEIN_ID}.pdb"

# --- 4. HEADER-SAFE SAMPLING ---
echo "Sampling actives and decoys (skipping headers)..."
python3 <<EOF
import random
import os

def get_mols(filename, count):
    if not os.path.exists(filename):
        return []
    with open(filename, 'r') as f:
        content = f.read()
    
    # Split by molecule tag
    parts = content.split("@<TRIPOS>MOLECULE")
    
    # The first part (parts[0]) is always the file header/comments. Skip it.
    mols = ["@<TRIPOS>MOLECULE" + p for p in parts[1:] if p.strip()]
    
    print(f"Found {len(mols)} valid molecules in {filename}")
    return random.sample(mols, min(count, len(mols)))

# Collect sampled records
sampled_actives = get_mols("$ACTIVES_IN", $NUM_ACTIVES)
sampled_decoys = get_mols("$DECOYS_IN", $NUM_DECOYS)
sampled_data = sampled_actives + sampled_decoys

if not sampled_data:
    print("ERROR: No molecules found! Check your input paths.")
    exit(1)

# Save combined Mol2 and build entries.txt
with open("$TOY_MOL2", "w") as f_mol, open("$ENTRIES_FILE", "w") as f_ent:
    for mol_block in sampled_data:
        f_mol.write(mol_block)
        
        # Extract name: Skip the tag and find the first line that isn't a comment
        lines = mol_block.split('\n')
        for line in lines:
            line = line.strip()
            if not line or line.startswith("@") or line.startswith("#"):
                continue
            # The first non-comment, non-tag line is our Name
            f_ent.write(f"{line}\n")
            break

print(f"Created $ENTRIES_FILE and $TOY_MOL2 with {len(sampled_data)} compounds.")
EOF

# --- 5. EXECUTE LUNA ---
echo "🧬 Running LUNA calculation..."
python ./luna/run.py \
    -p "$PROTEIN_ID" \
    -l "$TOY_MOL2" \
    -e "$ENTRIES_FILE" \
    -w "$RESULTS_DIR" \
    --pse \
    --ifp \
    --nproc 11

echo "=================================================="
echo "BENCHMARK COMPLETE!"
echo "Results Folder: $RESULTS_DIR"
echo "IFP File:      $RESULTS_DIR/results/fingerprints/ifp.csv"
echo "=================================================="
