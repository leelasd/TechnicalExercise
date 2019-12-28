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
  factory : RDKit Pharmacophore feat factory 
      Fingerprint radius.
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

class MACCS(Featurizer):
  """
  Circular (Morgan) fingerprints.

  Parameters
  ----------
  factory : RDKit Pharmacophore feat factory 
      Fingerprint radius.
  """
  name = 'MACCS'
#  from rdkit.Chem.Pharm2D import Gobbi_Pharm2D,Generate
  def __init__(self):
    pass

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
    from rdkit.Chem.rdMolDescriptors import GetMACCSKeysFingerprint
    fp = rdMolDescriptors.GetMACCSKeysFingerprint(mol)
    arr = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

