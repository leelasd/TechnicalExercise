import rdkit
from rdkit import Chem
from rdkit.Chem import PandasTools
import glob
sdflist = glob.glob('*.sdf') 
for sdf in sdflist:
    print(sdf)
    db = [m for m in Chem.SDMolSupplier(sdf,removeHs=False)]
    label = db[0].GetProp('_Name')
    for i,mol in enumerate(db): 
        Chem.MolToPDBFile(mol,label+'_%d.pdb'%i)
