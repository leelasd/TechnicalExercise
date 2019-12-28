"""
Topological fingerprints.
"""
from __future__ import division
from __future__ import unicode_literals

__author__ = "Leela Dodda"
__copyright__ = "Copyright 2014, Silicon Therapeutics"
__license__ = "MIT"

from deepchem.feat import Featurizer
import numpy as np

class Pharm2D(Featurizer):
  """
  Circular (Morgan) fingerprints.

  Parameters
  ----------
  radius : int, optional (default 2)
      Fingerprint radius.
  size : int, optional (default 2048)
      Length of generated bit vector.
  chiral : bool, optional (default False)
      Whether to consider chirality in fingerprint generation.
  bonds : bool, optional (default True)
      Whether to consider bond order in fingerprint generation.
  features : bool, optional (default False)
      Whether to use feature information instead of atom information; see
      RDKit docs for more info.
  sparse : bool, optional (default False)
      Whether to return a dict for each molecule containing the sparse
      fingerprint.
  smiles : bool, optional (default False)
      Whether to calculate SMILES strings for fragment IDs (only applicable
      when calculating sparse fingerprints).
  """
  name = 'circular'
  from rdkit.Chem.Pharm2D import Gobbi_Pharm2D,Generate

  def __init__(self,
               factory=Gobbi_Pharm2D.factory):
    self.factory = factory

  def _featurize(self, mol):
    """
    Calculate circular fingerprint.

    Parameters
    ----------
    mol : RDKit Mol
        Molecule.
    """
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors, DataStructs
    from rdkit.Chem.Pharm2D import Gobbi_Pharm2D,Generate
    fp = Generate.Gen2DFingerprint(mol,Gobbi_Pharm2D.factory)
    exp_fp = DataStructs.cDataStructs.ConvertToExplicit(fp)
    arr = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(exp_fp, arr)
    return arr

