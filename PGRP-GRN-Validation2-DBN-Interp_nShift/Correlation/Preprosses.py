# -------------------------------------------
# Haonan Tong
# -------------------------------------------
import pandas as pd
import numpy as np
import os

if not os.path.exists('data'):
	os.makedirs('data')
	
db = pd.read_csv('myDB_Validation.csv',index_col = 'id')



# Expression for TFs 
for atp in range(1,7):
	db_tmp = db.loc[(db['ATP']==atp) & (db['TF']==1), 't0':'t6']
	db_tmp = np.log2(db_tmp)
	db_tmp.to_csv('data/log2Expr_Val_TF_AT'+str(atp)+'.csv')

# Expression for nTFs 
for atp in range(1,7):
	db_tmp = db.loc[(db['Chen_isTar']==1) & (db['ATP']==atp), 't0':'t6']
	db_tmp = np.log2(db_tmp)
 	db_tmp.to_csv('data/log2Expr_Val_nTF_AT'+str(atp)+'.csv')