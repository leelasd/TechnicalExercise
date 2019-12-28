# coding: utf-8
import oddt
from oddt.fingerprints import SPLIF, similarity_SPLIF
from oddt.interactions import hbonds, hydrophobic_contacts, pi_stacking, salt_bridges, halogenbonds
import numpy as np
from rdkit import Chem
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
import matplotlib.pyplot as plt
#protein = next(oddt.toolkit.readfile('pdb', 'ali-model.pdb'))
protein = next(oddt.toolkit.readfile('pdb', '../CAMK2G/receptor.pdbqt'))
import glob
files= glob.glob('../CAMK2G/docked_*.sdf')
for ifile in files:
    print(ifile)
    mols = list(oddt.toolkit.readfile('sdf', ifile))
    print('HphobCont,Hbonds')
    for mol in mols: 
        conts = hydrophobic_contacts(protein,mol)
        res = hbonds(protein,mol)
        print(ifile, len(conts[0]),len(res[0]))
