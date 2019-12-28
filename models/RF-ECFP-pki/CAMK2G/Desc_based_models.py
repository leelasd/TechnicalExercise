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
#featurizer = dc.feat.CircularFingerprint(radius=2,size=1024)
featurizer = dc.feat.RDKitDescriptors() 
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
type_split = 'RandomStratified' #'RandomStratified'
splitter = splitters[type_split]
train_dataset, valid_dataset, test_dataset = splitter.train_valid_test_split(dataset,frac_train=.60,frac_valid=.20,frac_test=.20,seed=42)
transformers = [
    dc.trans.NormalizationTransformer(transform_y=True, dataset=train_dataset)
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

metric = dc.metrics.Metric(dc.metrics.r2_score)
optimizer = dc.hyper.HyperparamOpt(rf_model_builder)
best_rf, best_rf_hyperparams, all_rf_results = optimizer.hyperparam_search(
    params_dict, train_dataset, valid_dataset, transformers,
    metric=metric)
print(best_rf_hyperparams)
filename = 'finalized_model_CAMK2G_%s.sav'%type_split
joblib.dump(best_rf.model_instance, filename)
joblib.load(filename)
#dataset.save_to_disk()
#save_dataset_to_disk('./',train_dataset, valid_dataset, test_dataset, transformer) 
train_y_pred = best_rf.predict(train_dataset)
test_y_pred = best_rf.predict(test_dataset)
valid_y_pred = best_rf.predict(valid_dataset)
### Plotting Dataset
fig, (ax1,ax2,ax3) = plt.subplots(1, 3,figsize=(15,5))
ax1.scatter(train_dataset.y,train_y_pred,c='r',label=GetLable(train_dataset.y,train_y_pred))
ax2.scatter(test_dataset.y,test_y_pred,c='b',  label=GetLable(test_dataset.y,test_y_pred))
ax3.scatter(valid_dataset.y,valid_y_pred,c='g',label=GetLable(valid_dataset.y,valid_y_pred))
ax1.legend()
ax2.legend()
ax3.legend()
ax1.set_title('Training Data')
ax2.set_title('Test Data')
ax3.set_title('Validation Data')
ax1.set_ylabel('Predicted')
ax2.set_ylabel('Predicted')
ax3.set_ylabel('Predicted')
ax1.set_xlabel('Experimental')
ax2.set_xlabel('Experimental')
ax3.set_xlabel('Experimental')
plt.tight_layout()
plt.savefig('ModelPerformance_CAMK2G_%s.png'%type_split,dpi=300)
plt.show()
