#This code perform TF to PRF1 regressoin without removing cells that have values of PRF1 <= 0
#Updated the code later including TF ZNF354C in the tf_list
# import dependencies
import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
import seaborn as sns
import matplotlib.pyplot as plt
import math
from sklearn.preprocessing import StandardScaler
from sklearn import ensemble
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
import pickle


original_data_dir = "/dccstor/bmfm-targets/data/omics/transcriptome/scRNA/finetune/TIL/GSE156728/h5ad/btrmv2/combat with normalize and log1p.qc.all.h5ad"#combat with normalize and log1p.h5ad"
adata_1 = sc.read_h5ad(original_data_dir)

# Reading TF files
f_tfs = "allTFs_hg38.txt"
tf_list = pd.read_csv(f_tfs, index_col=False, sep=',')
## adding "ZNF354C" in TF list as it was read as header
tf_list.loc[len(tf_list.index)] = ["ZNF354C"]

#retaining the TFs that avaiable in our dataset
refined_tf = []
cnt =0
for i in range(0,len(tf_list)):
    tf = tf_list["ZNF354C"][i]
    if tf in adata_1.var_names:
        #print(tf)
        cnt +=1
        refined_tf.append(tf)

TF_adata_PRF1 = adata_1[:,refined_tf]
y_adata = adata_1[:,'PRF1']


X = TF_adata_PRF1.X
y = y_adata.X


#training GradientBoosting model and testing
X_train, X_test, y_train, y_test = train_test_split(X,y, test_size=0.3, random_state=44)
params = {
    "n_estimators": 500,
    "max_depth": 4,
    "min_samples_split": 5,
    "learning_rate": 0.01,
    "loss": "squared_error",
}

reg_PRF1 = ensemble.GradientBoostingRegressor(**params)
reg_PRF1.fit(X_train, y_train)

print("training complete")

predictions = reg_PRF1.predict(X_test)
mse = mean_squared_error(y_test.astype('float64'), predictions)
error = math.sqrt(mse)
print('Mean Srquared Error:', mse)
print("root mean square error {}".format(error))

error_list = [['mse',mse], ['error', error], ['len of TF', len(refined_tf)]]
df = pd.DataFrame(error_list, columns=['name', 'value'])
df.to_csv( "/dccstor/bmfm-targets/users/tanwi/research_data/TIL_DataAnalysis/TFtoPRF1_GBR_allCells_ZNF354C.txt", sep='\t')

# Saving model trained with grid search

saved_model_GBR = "/dccstor/bmfm-targets/users/tanwi/research_data/TIL_DataAnalysis/GBR_TFtoPRF1_allCells_ZNF354C.sav"
pickle.dump(reg_PRF1, open(saved_model_GBR, 'wb'))