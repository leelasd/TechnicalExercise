import deepchem as dc
from sklearn.ensemble import RandomForestRegressor
from deepchem.utils.evaluate import Evaluator
from sklearn.externals import joblib
from deepchem.utils.save import save_dataset_to_disk
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error, r2_score
from math import sqrt

def  GetLable(y,yp):
     #pearson_r2 = pearsonr(y, yp)[0] 
     if r2_score(y,yp) > 0:
         pearson_r2 = sqrt(r2_score(y,yp))
     else: pearson_r2 =0.0
     mse = mean_squared_error(y,yp)
     label = 'RMSE:%.2f Pearson:%.2f'%(mse, pearson_r2)
     return(label)

print('''
Creates RF Regression Model for Prediction of cMet/EGFR Binding pIC50s
#splitter = splitters['Butina']
#Only for Butina Clustering train_dataset, valid_dataset, test_dataset = splitter.train_valid_test_split(dataset,log_every_n=1000,cutoff=0.2)
'''
)

def rf_model_builder(model_params, model_dir):
    sklearn_model = RandomForestRegressor(**model_params)
    return dc.models.SklearnModel(sklearn_model, model_dir)
dataset_file = '../../../data/pki/CAMK2G_pki.csv'
df = pd.read_csv(dataset_file)
print(df.describe())
featurizer = dc.feat.CircularFingerprint(radius=2,size=1024)
loader = dc.data.CSVLoader(
      tasks=["pchembl_value"], smiles_field="canonical_smiles",
      featurizer=featurizer)
dataset = loader.featurize(dataset_file)
#### Different Splitters 
splitters = {
		'Scaffold': dc.splits.ScaffoldSplitter(dataset),
		'RandomStratified':dc.splits.RandomStratifiedSplitter(dataset),
		'MaxMin': dc.splits.MaxMinSplitter(dataset),
		'FingerPrint': dc.splits.FingerprintSplitter(dataset),
                'Butina': dc.splits.ButinaSplitter(dataset),
             }
type_split = 'Scaffold' #'RandomStratified'
splitter = splitters[type_split]
train_dataset, valid_dataset, test_dataset = splitter.train_valid_test_split(dataset,frac_train=.60,frac_valid=.20,frac_test=.20,seed=42)
transformers = [
    dc.trans.NormalizationTransformer(transform_y=True, dataset=train_dataset)
    ]

for dataset in [train_dataset, valid_dataset, test_dataset]:
    for transformer in transformers:
        dataset = transformer.transform(dataset)

