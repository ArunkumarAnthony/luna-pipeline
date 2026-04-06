import sys
import argparse
import CDPL.Chem as Chem
import CDPL.Biomol as Biomol
import CDPL.Pharm as Pharm
import CDPL.MolProp as MolProp

def readAndPrepareReceptorStructure(args: argparse.Namespace) -> Chem.Molecule:
    reader = Chem.MoleculeReader(args.receptor_file) 
    
    sup_fmts = [ Chem.DataFormat.MOL2, Biomol.DataFormat.PDB, Biomol.DataFormat.MMTF, Biomol.DataFormat.MMCIF ]
                        
    if reader.getDataFormat() not in sup_fmts:
        sys.exit('Error: receptor input file format not supported')

    rec_mol = Chem.BasicMolecule()
    try:
        if not reader.read(rec_mol):
            sys.exit('Error: reading receptor structure failed')
    except Exception as e:
        sys.exit('Error: reading receptor structure failed:\n' + str(e))            

    try:
        rem_atoms = False

        if args.strip_res_list:            
            atoms_to_rem = Chem.Fragment() 
            res_to_strip = { tlc.upper() for tlc in args.strip_res_list }
        
            for atom in rec_mol.atoms:     
                if Biomol.getResidueCode(atom).upper() in res_to_strip:
                    atoms_to_rem.addAtom(atom)

            if atoms_to_rem.numAtoms > 0:
                rec_mol -= atoms_to_rem    
                rem_atoms = True
                if not args.quiet:
                    print('! Removed %s atoms from the receptor structure' % str(atoms_to_rem.numAtoms))

        Chem.perceiveSSSR(rec_mol, rem_atoms)
        Chem.setRingFlags(rec_mol, rem_atoms)
        Chem.calcImplicitHydrogenCounts(rec_mol, rem_atoms)
        Chem.perceiveHybridizationStates(rec_mol, rem_atoms)
        Chem.setAromaticityFlags(rec_mol, rem_atoms)

        if Chem.makeHydrogenComplete(rec_mol):                    
            Chem.calcHydrogen3DCoordinates(rec_mol)               
            Biomol.setHydrogenResidueSequenceInfo(rec_mol, False) 

        MolProp.calcAtomHydrophobicities(rec_mol, False)          

    except Exception as e:
        sys.exit('Error: processing of receptor structure failed: ' + str(e))            

    return rec_mol

def parseArgs() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Generates pharmacophores describing the interactions between a given receptor structure and a set of ligand molecules.')
    parser.add_argument('-r', dest='receptor_file', required=True, metavar='<file>', help='Receptor structure input file (*.mol2, *.pdb, *.mmtf, *.cif, *.mmcif)')
    parser.add_argument('-l', dest='ligands_file', required=True, metavar='<file>', help='Ligand structure input file (*.sdf, *.mol2, *.cdf)')
    parser.add_argument('-o', dest='out_file', required=True, metavar='<file>', help='Pharmacophore output file (*.pml, *.cdf)')
    parser.add_argument('-s', dest='strip_res_list', required=False, metavar='<three-letter code>', nargs='+', help='Residues to remove')
    parser.add_argument('-q', dest='quiet', required=False, action='store_true', default=False, help='Disable progress output')
    parser.add_argument('-x', dest='gen_x_vols', required=False, action='store_true', default=False, help='Generate exclusion volumes')
    return parser.parse_args()

def main() -> None:
    args = parseArgs()

    rec_mol = readAndPrepareReceptorStructure(args)          
    lig_reader = Chem.MoleculeReader(args.ligands_file)      
    ph4_writer = Pharm.FeatureContainerWriter(args.out_file) 
 
    lig_mol = Chem.BasicMolecule()          
    ia_ph4 = Pharm.BasicPharmacophore()     
    ph4_gen = Pharm.InteractionPharmacophoreGenerator() 
    ph4_gen.addExclusionVolumes(args.gen_x_vols)        
                                                        
    try:
        i = 1
        while lig_reader.read(lig_mol):
            mol_id = Chem.getName(lig_mol).strip() 
            if mol_id == '':
                mol_id = '#' + str(i)  
            else:
                mol_id = '\'%s\' (#%s)' % (mol_id, str(i))

            if not args.quiet:
                print('- Generating interaction pharmacophore of molecule %s...' % mol_id)

            try:
                Pharm.prepareForPharmacophoreGeneration(lig_mol) 
                ph4_gen.generate(lig_mol, rec_mol, ia_ph4, True) 

                if not args.quiet:
                     print(' -> Generated %s features: %s' % (str(ia_ph4.numFeatures), Pharm.generateFeatureTypeHistogramString(ia_ph4)))
                
                try:
                    if not ph4_writer.write(ia_ph4): 
                        sys.exit('Error: writing interaction pharmacophore of molecule %s failed' % mol_id)
                except Exception as e:               
                    sys.exit('Error: writing interaction pharmacophore of molecule %s failed: %s' % (mol_id, str(e)))
                
            except Exception as e:                   
                sys.exit('Error: interaction pharmacophore generation for molecule %s failed: %s' % (mol_id, str(e)))
            i += 1
            
    except Exception as e:                           
        sys.exit('Error: reading molecule %s failed: %s' % (str(i), str(e)))

    if not args.quiet:
        print('Done!')

    ph4_writer.close()
    sys.exit(0)
        
if __name__ == '__main__':
    main()
