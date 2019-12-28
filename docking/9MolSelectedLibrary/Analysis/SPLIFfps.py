# coding: utf-8
import oddt
from oddt.fingerprints import SPLIF, similarity_SPLIF
import numpy as np
from rdkit import Chem
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
import matplotlib.pyplot as plt
protein = next(oddt.toolkit.readfile('pdb', 'ali-model.pdb'))
mols = list(oddt.toolkit.readfile('sdf', 'docked_DB12964_0.sdf'))
all_fps = [ SPLIF(mols[i],protein,depth=2, size=4096, distance_cutoff=4.5) for i in range(len(mols))]
X = np.array(all_fps)
dists = np.zeros((len(all_fps),len(all_fps)))
for i in range(len(all_fps)):
    for j in range(i,len(all_fps)): 
        dist = similarity_SPLIF(all_fps[i],all_fps[j])
        dists[i][j] = dist 
        dists[j][i] = dist 

model_com = AgglomerativeClustering(affinity='precomputed', n_clusters=4, linkage='complete').fit(dists)
model_avg = AgglomerativeClustering(affinity='precomputed', n_clusters=4, linkage='average').fit(dists)
mol_names = ['Name' for i in range(len(all_fps))]
df = pd.DataFrame({'names': mol_names,'cid_avg': model_avg.labels_,'cid_com': model_com.labels_})
df.to_csv('res.csv',index=False)
