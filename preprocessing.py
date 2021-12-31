# -*- coding: utf-8 -*-
"""
Created on Thu Dec 23 11:21:41 2021

@author: idscmm
"""

import pandas as pd
import copy
import numpy as np
from sklearn.impute import KNNImputer

allvar = 0

def f_get_Normalization(X, norm_mode):    
    num_Patient, num_Feature = np.shape(X)
    
    if norm_mode == 'standard': #zero mean unit variance
        for j in range(num_Feature):
            if np.nanstd(X[:,j]) != 0:
                X[:,j] = (X[:,j] - np.nanmean(X[:, j]))/np.nanstd(X[:,j])
            else:
                X[:,j] = (X[:,j] - np.nanmean(X[:, j]))
    elif norm_mode == 'normal': #min-max normalization
        for j in range(num_Feature):
            X[:,j] = (X[:,j] - np.nanmin(X[:,j]))/(np.nanmax(X[:,j]) - np.nanmin(X[:,j]))
    else:
        print("INPUT MODE ERROR!")
    
    return X

Patient_characteristics=pd.read_csv(r"C:\Users\idscmm\Downloads\RWDDrugEfficacy\Biostats - RWD drug efficacy\Patient_characteristics.csv")
Event_duration = pd.read_csv(r"C:\Users\idscmm\Downloads\RWDDrugEfficacy\Biostats - RWD drug efficacy\Event_duration.csv")
data = copy.deepcopy(Patient_characteristics)

cat_col = ['treatment_variable', 'sex','other_drugs_1','other_drugs_2', 'other_drugs_3', 'other_drugs_4', 'other_drugs_5',
           'other_drugs_6', 'other_drugs_7', 'other_drugs_8', 'diagnosis_1',
           'diagnosis_2', 'diagnosis_3', 'diagnosis_4', 'diagnosis_5',
           'diagnosis_6', 'diagnosis_7', 'diagnosis_8', 'diagnosis_9',
           'diagnosis_10', 'diagnosis_11', 'diagnosis_12', 'diagnosis_13',
           'diagnosis_14', 'diagnosis_15','Diag_Score_1', 'Diag_Score_2']
cont_col = ['age','lab_1', 'lab_6', 'lab_7', 'lab_8','lab_5','lab_2', 'lab_4', 'lab_3']

for i in cat_col:
    locals()[i] = pd.get_dummies(Patient_characteristics[i])
    
Diag_Score_1.columns = ['Diag_Score_1_' + str(x) for x in list(Diag_Score_1.columns)]
Diag_Score_2.columns = ['Diag_Score_2_' + str(x) for x in list(Diag_Score_2.columns)]

data = data.drop(columns = ['Diag_Score_1','Diag_Score_2'])  
data = pd.concat([data, Diag_Score_1.iloc[:,:-1], Diag_Score_2.iloc[:,:-1]],axis=1)  
for i in cat_col[:-2]:
    if i == 'treatment_variable':
        data[i] = treatment_variable['Drug_A']
    elif i == 'sex':
        data[i] = sex[1]
    else:
        data[i] = locals()[i]['Yes']

cat_col = cat_col[:-2] + list(Diag_Score_1.columns)[:-1] + list(Diag_Score_2.columns)[:-1]

# Scaling
#data.loc[:,cont_col] = f_get_Normalization(np.array(data[cont_col]),'standard')


    
imputer = KNNImputer(n_neighbors=3)
if allvar == 1:
    data_complete = imputer.fit_transform(data[cat_col + cont_col])
    data.loc[:,cat_col + cont_col] = data_complete
    data.to_csv(r'C:\Users\idscmm\Downloads\RWDDrugEfficacy\Processed_Data\data_knn_allvariables.csv',index=False)
else:
    data_complete = imputer.fit_transform(data[cat_col + cont_col[:6]])
    data.loc[:,cat_col + cont_col[:6]] = data_complete
    data.to_csv(r'C:\Users\idscmm\Downloads\RWDDrugEfficacy\Processed_Data\data_knn_4variables.csv',index=False)    
    
    



