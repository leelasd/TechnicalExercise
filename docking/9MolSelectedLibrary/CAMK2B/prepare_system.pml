load model.pdb,protein 
load lig_3bhh_CAMK2B.sdf, ligand 
h_add elem O or elem N
save native.pdb, ligand
save receptor.pdb, protein
