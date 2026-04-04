import pandas as pd

ACTIVES_FILE = 'inputs/MAPK2/dudez_1pt0LD_ligand_poses.mol2'
IFP_CSV = 'results/results/fingerprints/ifp.csv'
OUTPUT_CSV = 'labeled_ml_dataset.csv'

print("🔍 Scanning original files for Active IDs...")
active_ids = set()

# 1. Memorize all the Active IDs
with open(ACTIVES_FILE, 'r') as f:
    lines = f.readlines()
    for i, line in enumerate(lines):
        if "@<TRIPOS>MOLECULE" in line:
            for j in range(i+1, len(lines)):
                candidate = lines[j].strip()
                if not candidate.startswith('#') and not candidate.startswith('@'):
                    # Strip invisible spaces just in case
                    active_ids.add(candidate.strip())
                    break

print(f"📂 Loading LUNA fingerprints from {IFP_CSV}...")
df = pd.read_csv(IFP_CSV)
id_col = df.columns[0]

print("🏷️ Applying labels...")
# BULLETPROOF FIX: Remove 'rec:' and strip invisible spaces before checking
df['Label'] = df[id_col].apply(lambda x: 1 if str(x).replace("rec:", "").strip() in active_ids else 0)

# Move Label to the front
cols = ['Label'] + [c for c in df.columns if c != 'Label']
df = df[cols]

df.to_csv(OUTPUT_CSV, index=False)

print("\n========================================")
print("✅ ML DATASET READY!")
print(f"Total compounds processed: {len(df)}")
print(f"🎯 Actives (Class 1): {df['Label'].sum()}")
print(f"🛑 Decoys (Class 0): {len(df) - df['Label'].sum()}")
print("========================================")
