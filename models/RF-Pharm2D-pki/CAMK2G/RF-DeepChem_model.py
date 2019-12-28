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
from utils import *
from PharmFinger import Pharm2D 

def FitModel(dataset_file,splitter_type):
    #featurizer = dc.feat.CircularFingerprint(radius=2,size=1024)
    featurizer = Pharm2D() 
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
    #for type_split in  ['Scaffold','RandomStratified','MaxMin','Butina']: #list(splitters.keys())
    type_split = splitter_type 
    print(type_split)
    metric = dc.metrics.Metric(dc.metrics.r2_score)
    splitter = splitters[type_split]
    train_dataset, valid_dataset, test_dataset = splitter.train_valid_test_split(dataset,frac_train=.60,frac_valid=.20,frac_test=.20,seed=42)
    transformers = [
        dc.trans.NormalizationTransformer(transform_X=True, dataset=train_dataset),
        ]
    
    for dataset in [train_dataset, valid_dataset, test_dataset]:
        for transformer in transformers:
            dataset = transformer.transform(dataset)
    
    params_dict = {
        'bootstrap': [True],
        'max_depth': [5,10,20,25,50],
        'max_features': [3],
        'min_samples_leaf': [3],
        'min_samples_split': [10],
        'n_estimators': [10,20,25,50,100],
        'max_features': ["auto", "sqrt", "log2", None], 
    }
    
    optimizer = dc.hyper.HyperparamOpt(rf_model_builder,verbose=True)
    best_rf, best_rf_hyperparams, all_rf_results = optimizer.hyperparam_search(
        params_dict, train_dataset, valid_dataset, transformers,
        metric=metric,logdir='opt_model_CAMK2G_%s'%type_split)
    filename = 'finalized_model_CAMK2G_%s.sav'%type_split
    print(type_split)
    PlotModel(best_rf,type_split,train_dataset,test_dataset,valid_dataset)
    joblib.dump(best_rf.model_instance, filename)
    return None

if __name__ == "__main__": 
    dataset_file = '../../../data/pki/CAMK2G_pki.csv'
    df = pd.read_csv(dataset_file)
    print(df.describe())
    FitModel(dataset_file,splitter_type='RandomStratified')
    FitModel(dataset_file,splitter_type='Scaffold')
    FitModel(dataset_file,splitter_type='MaxMin')
    FitModel(dataset_file,splitter_type='FingerPrint')
