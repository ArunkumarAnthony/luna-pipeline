import sys
import CDPL.Chem as Chem
import CDPL.Biomol as Biomol
import CDPL.Pharm as Pharm

def main():
    receptor_file = "data/a_baumannii/loldf/protein/loldf-abaumannii.pdb"
    ligand_file = "data/a_baumannii/loldf/ligand/abaucin.sdf"
    
    # We use .pml (LigandScout format), which is the standard 3D format for PharmacoMatch
    output_file = "data/a_baumannii/loldf/control_pharm.pml"

    print("Loading structures...")
    
    # 1. Load Receptor
    reader_rec = Chem.MoleculeReader(receptor_file)
    rec_mol = Chem.BasicMolecule()
    if not reader_rec.read(rec_mol):
        sys.exit("Error: Failed to read receptor.")
        
    # 2. Load Ligand
    reader_lig = Chem.MoleculeReader(ligand_file)
    lig_mol = Chem.BasicMolecule()
    if not reader_lig.read(lig_mol):
        sys.exit("Error: Failed to read ligand.")
        
    print("Preparing molecules...")
    # Prepare Receptor (Perceive bonds, rings, and aromaticity)
    Chem.calcImplicitHydrogenCounts(rec_mol, False)
    Chem.perceiveHybridizationStates(rec_mol, False)
    Chem.perceiveSSSR(rec_mol, False)
    Chem.setRingFlags(rec_mol, False)
    Chem.setAromaticityFlags(rec_mol, False)
    Biomol.perceiveResidues(rec_mol, False)
    
    # Prepare Ligand
    Pharm.prepareForPharmacophoreGeneration(lig_mol)
    
    print("Extracting interaction pharmacophore...")
    ph4_gen = Pharm.InteractionPharmacophoreGenerator()
    ph4 = Pharm.BasicPharmacophore()
    
    # Generate the pharmacophore from the complex
    ph4_gen.generate(rec_mol, lig_mol, ph4)
    
    print(f"Saving to {output_file}...")
    writer = Pharm.FeatureContainerWriter(output_file)
    if not writer.write(ph4):
        sys.exit("Error: Failed to write pharmacophore.")
    writer.close()
    
    print("Success! Pharmacophore template generated.")

if __name__ == "__main__":
    main()