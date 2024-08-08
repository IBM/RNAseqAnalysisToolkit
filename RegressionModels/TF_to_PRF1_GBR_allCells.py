
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


original_data_dir = "combat with normalize and log1p.qc.all.h5ad"#combat with normalize and log1p.h5ad"
adata_1 = sc.read_h5ad(original_data_dir)



X = TF_adata_PRF1.X
y = y_adata.X


#training GradientBoosting model and testing
X_train, X_test, y_train, y_test = train_test_split(X,y, test_size=0.3, random_state=44)
params = {
    "n_estimators": 100,
    "max_depth": 3,
    "learning_rate": 0.01,
    "loss": "squared_error",
}


print("training complete")

predictions = reg_PRF1.predict(X_test)
mse = mean_squared_error(y_test.astype('float64'), predictions)
error = math.sqrt(mse)
print("root mean square error {}".format(error))

