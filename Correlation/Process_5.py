# -------------------------------------------
# Haonan Tong
# -------------------------------------------
import pandas as pd 
import numpy as np 
from scipy.stats import pearsonr



def f_getPosition(db, bdeu):
	l = sum(db.abs() > np.absolute(bdeu))
	u = sum(db.abs() >=np.absolute(bdeu))
	return l,u


# Use Pearson correlation analyze correlation of target activated at a time point
# and potential regulator activated at any time
# Output Pearson correlation matrix would look like this
# 	G1	G2	G3 ...
# TF1
# TF2
# TF3
# ...
# For nTF activated at 1 - 5

Directory_nTF = 'Data/log2Expr_Val_nTF_AT'
Directory_TF = 'Data/log2Expr_Val_TF_AT'
Directory_output = 'Data/PearsonCorr_Val'


# DataFrame for all Target	
frames = [];
for atp_tar in range(2,7):
	# DataFrame for Regulator
	df_reg_tmp = pd.read_csv(Directory_nTF+str(atp_tar)+'.csv', index_col = 'id')
	frames.append(df_reg_tmp)
	df_tar = pd.concat(frames)
 

# DataFrame for Potential Regulator
frames = [];
for atp_reg in range(1,7):
	# DataFrame for Regulator
	df_reg_tmp = pd.read_csv(Directory_TF+str(atp_reg)+'.csv', index_col = 'id')
	frames.append(df_reg_tmp)
	df_reg = pd.concat(frames)

# Derive Pearson Correltion Matrix
ntar = len(df_tar)
nreg = len(df_reg)
Corr_Matrx = np.zeros((nreg,ntar))
for tar in range(0,ntar):
	for reg in range(0,nreg):
		star = df_tar.iloc[tar,:]
		sreg = df_reg.iloc[reg,:]

		# print star.corr( sreg, method='pearson')
		Corr_Matrx[reg, tar] = star.corr( sreg, method='pearson')

# print Corr_Matrx
df_Corr_Mtrx = pd.DataFrame(Corr_Matrx, index=df_reg.index.values, columns=df_tar.index.values )
df_Corr_Mtrx.index.name = 'id'
df_Corr_Mtrx.to_csv(Directory_output+'.csv')

# Analysis 
Dir_Analysis = 'Target_Analysis/'

# DETarList
db = pd.read_csv('../myDB_Validation.csv',index_col='id')
db_DETarList = db[db['Chen_isTar'] == 1]
DETarList = list(db_DETarList.index)

# print df_Corr_Mtrx
# print DETarList
for gtar in DETarList:
	df_score_tmp = df_Corr_Mtrx.loc[:,gtar]
	dir_tar = Dir_Analysis + gtar + '.csv'
	df_target_tmp = pd.read_csv(dir_tar, index_col='Transcription factor')

 	lower_index = [];
	upper_index = [];
	for i in range(0,len(df_score_tmp)):
		bdeu = df_score_tmp.iloc[i]
		[l, u] = f_getPosition(df_score_tmp, bdeu)
		lower_index.append(l)
		upper_index.append(u)

	df_score_tmp = df_score_tmp.to_frame()

	df_score_tmp['lower_index'] = lower_index
	df_score_tmp['upper_index'] = upper_index
	result = pd.concat([df_target_tmp, df_score_tmp], axis=1, join_axes=[df_target_tmp.index])
	result.to_csv(Dir_Analysis+gtar+'_report.csv')


















