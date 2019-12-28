import deepchem as dc
from sklearn.ensemble import RandomForestRegressor
from matplotlib import pyplot as plt
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error, r2_score
from math import sqrt

def GetLable(y,yp):
     if r2_score(y,yp) > 0:
         pearson_r2 = sqrt(r2_score(y,yp))
     else: pearson_r2 =0.0
     mse = mean_squared_error(y,yp)
     label = 'RMSE:%.2f Pearson:%.2f'%(mse, pearson_r2)
     return(label)

def PlotModel(best_rf,type_split,train_dataset,test_dataset,valid_dataset):
    train_y_pred = best_rf.predict(train_dataset)
    test_y_pred = best_rf.predict(test_dataset)
    valid_y_pred = best_rf.predict(valid_dataset)
    ## Plotting Dataset
    train_label = GetLable(train_dataset.y,train_y_pred) 
    test_label  = GetLable(test_dataset.y,test_y_pred)
    valid_label = GetLable(valid_dataset.y,valid_y_pred) 
    fig, (ax1,ax2,ax3) = plt.subplots(1, 3,figsize=(15,5))
    ax1.scatter(train_dataset.y,train_y_pred,c='r',label=train_label)
    ax2.scatter(test_dataset.y,test_y_pred,c='b',  label=test_label)
    ax3.scatter(valid_dataset.y,valid_y_pred,c='g',label=valid_label)
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
    #plt.show()
    print('Training   %s'%train_label)
    print('Testing    %s'%test_label)
    print('Validation %s'%valid_label)
    return None

def rf_model_builder(model_params, model_dir):
    sklearn_model = RandomForestRegressor(**model_params)
    return dc.models.SklearnModel(sklearn_model, model_dir)

