# ------------------------------------------------------------------------------------
# Generate Files As Input For DBN Analysis
# Derive Expression data for genes of interest and named DEGs-time-Activation*.csv
# Derive Expression data for TFs of interest and named TFs-DEGs-time-Activation*.csv

import pandas as pd
import os
List_directory = ['Data','Data/TDEGs','Data/TTFs']
for directory in List_directory:
	if not os.path.exists(directory):
    		os.makedirs(directory)

db = pd.read_csv('myDB_Validation.csv',index_col='id')
db_DETarList = db[db['Chen_isTar'] == 1]

DETarList = list(db_DETarList.index)

# Generate DEGs-time-Activation*.csv
for atp in range(1,7):
	db_tmp = db_DETarList.loc[db['ATP']==atp,'kat_1_1':'kat_3_7']
	db_tmp.to_csv('Data/TDEGs/DEGs-time-Activation' + str(atp) + '.csv')

# Generate TFs-DEGs-time-Activation*.csv
db_TFList = db[db['TF']==1]
for atp in range(1,7):
	db_tmp = db_TFList.loc[db['ATP']==atp,'kat_1_1':'kat_3_7']
	db_tmp.to_csv('Data/TTFs/TFs-DEGs-time-Activation' + str(atp) + '.csv')