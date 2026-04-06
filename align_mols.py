import sys
import argparse
import CDPL.Chem as Chem
import CDPL.Pharm as Pharm
import CDPL.Math as Math

def readRefPharmacophore(filename: str) -> Pharm.Pharmacophore:
    reader = Pharm.PharmacophoreReader(filename)
    ph4 = Pharm.BasicPharmacophore()
    if not reader.read(ph4):
        sys.exit('Error: reading reference pharmacophore failed')
    return ph4

def genPharmacophore(mol: Chem.Molecule) -> Pharm.Pharmacophore:
    Pharm.prepareForPharmacophoreGeneration(mol)       
    ph4_gen = Pharm.DefaultPharmacophoreGenerator()    
    ph4 = Pharm.BasicPharmacophore()                   
    ph4_gen.generate(mol, ph4)                         
    return ph4

def clearFeatureOrientations(ph4: Pharm.BasicPharmacophore) -> None:
    for ftr in ph4:
        Pharm.clearOrientation(ftr)
        Pharm.setGeometry(ftr, Pharm.FeatureGeometry.SPHERE)
        
def parseArgs() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Aligns molecules onto a reference pharmacophore.')
    parser.add_argument('-r', dest='ref_ph4_file', required=True, metavar='<file>', help='Reference pharmacophore (*.pml)')
    parser.add_argument('-i', dest='in_file', required=True, metavar='<file>', help='Molecule input file')
    parser.add_argument('-o', dest='out_file', required=True, metavar='<file>', help='Aligned molecule output file')
    parser.add_argument('-n', dest='num_out_almnts', required=False, default=1, type=int, help='Max output solutions')
    parser.add_argument('-x', dest='exhaustive', required=False, action='store_true', default=False, help='Exhaustive search')
    parser.add_argument('-d', dest='min_pose_rmsd', required=False, default=0.0, type=float, help='Min RMSD')
    parser.add_argument('-q', dest='quiet', required=False, action='store_true', default=False, help='Disable progress')
    parser.add_argument('-p', dest='pos_only', required=False, action='store_true', default=False, help='Position matching only')
    return parser.parse_args()

def main() -> None:
    args = parseArgs()
    ref_ph4 = readRefPharmacophore(args.ref_ph4_file) 
    mol_reader = Chem.MoleculeReader(args.in_file) 
    
    # Crucial: Treat each conformer in the SDF as an individual molecule to screen
    Chem.setMultiConfImportParameter(mol_reader, False) 
    
    mol_writer = Chem.MolecularGraphWriter(args.out_file) 
    mol = Chem.BasicMolecule()
    almnt = Pharm.PharmacophoreAlignment(True) 

    if args.pos_only:                          
        clearFeatureOrientations(ref_ph4)
    
    almnt.addFeatures(ref_ph4, True)               
    almnt.performExhaustiveSearch(args.exhaustive) 
    almnt_score = Pharm.PharmacophoreFitScore()
    
    try:
        i = 1
        while mol_reader.read(mol):
            mol_id = Chem.getName(mol).strip() 
            if mol_id == '': mol_id = '#' + str(i)
            else: mol_id = "'%s' (#%s)" % (mol_id, str(i))

            if not args.quiet:
                print('- Aligning molecule conformer %s...' % mol_id)

            try:
                mol_ph4 = genPharmacophore(mol)    
                if args.pos_only:                  
                    clearFeatureOrientations(mol_ph4)

                almnt.clearEntities(False)         
                almnt.addFeatures(mol_ph4, False)  

                almnt_solutions = []               
                while almnt.nextAlignment():                                     
                    score = almnt_score(ref_ph4, mol_ph4, almnt.getTransform())  
                    xform = Math.Matrix4D(almnt.getTransform())                  
                    almnt_solutions.append((score, xform))

                if not args.quiet:
                    print(' -> Found %s alignment solutions' % str(len(almnt_solutions)))
                
                saved_coords = Math.Vector3DArray()      
                Chem.get3DCoordinates(mol, saved_coords) 

                if Chem.hasStructureData(mol):           
                    struct_data = Chem.getStructureData(mol)
                else:                                    
                    struct_data = Chem.StringDataBlock()
                    Chem.setStructureData(mol, struct_data)

                struct_data.addEntry('<PharmFitScore>', '') 
                output_cnt = 0
                last_pose = None
                almnt_solutions = sorted(almnt_solutions, key=lambda entry: entry[0], reverse=True)

                for solution in almnt_solutions:
                    if output_cnt == args.num_out_almnts: break
                    curr_pose = Math.Vector3DArray(saved_coords)
                    Math.transform(curr_pose, solution[1])  

                    if args.min_pose_rmsd > 0.0 and last_pose and Math.calcRMSD(last_pose, curr_pose) < args.min_pose_rmsd:
                        continue

                    Chem.set3DCoordinates(mol, curr_pose)  
                    struct_data[len(struct_data) - 1].setData(format(solution[0], '.4f'))     
                    
                    if not mol_writer.write(mol):
                        sys.exit('Error: writing alignment pose failed')

                    last_pose = curr_pose
                    output_cnt += 1

                if not args.quiet:
                    print(' -> %s alignment poses saved' % str(output_cnt))

            except Exception as e:
                print('Error: pharmacophore alignment failed: %s' % str(e))

            i += 1
                
    except Exception as e: 
        sys.exit('Error: reading input molecule failed: ' + str(e))

    mol_writer.close()
    sys.exit(0)

if __name__ == '__main__':
    main()
