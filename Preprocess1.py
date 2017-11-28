# For each target generate a summary of corresponding TFs and ATP;
# myDB_ValidationSet.csv: Contains All DEGs information including ATP, Expr, isTF,
# Functional Categories, if in Validation set; 

import pandas as pd
import os
directory = 'Data/Target_Summary/'
if not os.path.exists(directory):
    	os.makedirs(directory)

db = pd.read_csv('myDB_Validation.csv',index_col='id')
db_DETarList = db[db['Chen_isTar'] == 1]

DETarList = list(db_DETarList.index)

vSet = pd.read_csv('ValidationSet.csv',index_col='Transcription factor')

# DETarList = ['AT2G37090']
for x in DETarList:
	db_tmp = vSet[vSet.Promoter==x]

	# db_tmp = vSet[vSet.Promoter==x].copy()
	for y in list(db_tmp.index):
		db_tmp.loc[y,'Reg_ATP'] = db.loc[y,'ATP']
		db_tmp.loc[y,'Tar_ATP'] = db.loc[x,'ATP']
		db_tmp.to_csv(directory+x+'.csv')




