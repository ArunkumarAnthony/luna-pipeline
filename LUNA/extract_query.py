import pandas as pd

INPUT_CSV = 'labeled_ml_dataset.csv'
QUERY_OUT = 'master_query.csv'
TRAIN_OUT = 'training_dataset.csv'

print(f"📂 Loading {INPUT_CSV}...")
df = pd.read_csv(INPUT_CSV)

# Isolate the Actives (Label == 1)
actives = df[df['Label'] == 1]

if len(actives) == 0:
    print("❌ Error: No actives found in the dataset.")
    exit()

# Take the very first Active to be our One-Shot Control
query_df = actives.iloc[[0]]
# Column 0 is Label, Column 1 is the Molecule ID (LUNA puts it as rec:ID)
query_name = query_df.iloc[0, 1] 

print(f"🎯 Selected '{query_name}' as the Master Control Query!")

# Save the Query to its own file
query_df.to_csv(QUERY_OUT, index=False)

# Remove the Query from the main dataset so the AI doesn't cheat by training on it
remaining_df = df.drop(query_df.index)
remaining_df.to_csv(TRAIN_OUT, index=False)

print(f"💾 Saved Master Query to: {QUERY_OUT}")
print(f"💾 Saved remaining {len(remaining_df)} compounds to: {TRAIN_OUT}")
print("\n✅ READY FOR PHASE 4 (PyTorch)!")

