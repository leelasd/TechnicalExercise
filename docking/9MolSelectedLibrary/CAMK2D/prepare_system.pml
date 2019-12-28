load ali-model.pdb,protein 
load native.pdb, ligand 
h_add elem O or elem N
save native.pdb, ligand
save receptor.pdb, protein
