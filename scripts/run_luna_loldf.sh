# 1. Create a quick entries file for your ligand
echo "abaucin" > data/a_baumannii/loldf/entries.txt

# 2. Run the LUNA interaction extraction and fingerprint generation
python LUNA/luna/run.py \
  -p loldf-abaumannii \
  --pdbdir data/a_baumannii/loldf/protein \
  -l data/a_baumannii/loldf/ligand/abaucin.mol2 \
  -e data/a_baumannii/loldf/entries.txt \
  -w data/a_baumannii/loldf/luna_results \
  --ifp \
  -T EIFP \
  --ifp-out data/a_baumannii/loldf/luna_results/abaucin_eifp.csv \
  --pse \
  --nproc 1 \
  --overwrite
