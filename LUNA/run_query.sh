#!/bin/bash
set -e

echo "🚀 Generating Master Query Fingerprint..."

# 1. Bypass the LUNA path bug
mkdir -p ./results_query/pdbs
cp inputs/MAPK2/rec.crg.pdb ./results_query/pdbs/rec.pdb

# 2. Apply the chemistry fix (ILE CD -> CD1)
sed -i 's/ CD  ILE/ CD1 ILE/g' ./results_query/pdbs/rec.pdb

# 3. Create the entry file
echo "xtal-lig" > xtal_entry.txt

# 4. Run LUNA (Letting bash handle the pathing naturally)
python ./luna/run.py \
    -p rec \
    -l inputs/MAPK2/xtal-lig.pdb \
    -e xtal_entry.txt \
    -w ./results_query \
    --ifp

echo "=================================================="
echo "✅ Master Query generated successfully!"
echo "=================================================="
