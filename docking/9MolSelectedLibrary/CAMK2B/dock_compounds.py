import glob
import os
inp_files = glob.glob('../PDBQT/*.pdbqt')
for inp_file in inp_files:
    label = inp_file.split('/')[-1].split('.')[0]
    print(label)
    command = 'smina -r receptor.pdbqt -l %s --autobox_ligand native.pdbqt --autobox_add 8 --accurate_line --exhaustiveness 10 -o docked_%s.sdf'%(inp_file,label)
    os.system(command)
#flex_command = 'smina -r receptor.pdbqt -l %s --autobox_ligand native.pdbqt --flexdist_ligand native.pdbqt --flexdist 4 --autobox_add 8 --accurate_line --exhaustiveness 10 -o docked_%s.sdf'%(inp_file,label)
