import pandas as pd 
import numpy as np 
from scipy.stats import pearsonr
import seaborn as sns
import matplotlib.pyplot as plt

# Script for Pearson Correlation Analysis Without Interpolation
# Use Pearson correlation analyze correlation of target activated at a time point
# and potential regulator activated by (at and previous to) the time
# Output Pearson correlation matrix would look like this
# 	G1	G2	G3 ...
# TF1
# TF2
# TF3
# ...

######################################################################
######################################################################
# Generate Correlation matrix for each target in validation set vs TF 
db = pd.read_csv('myDB_Validation.csv',index_col = 'id')

df_tar = db.loc[db.Chen_isTar==1,'t0':'t6']
df_reg = db.loc[db.TF==1,'t0':'t6']

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
# print df_reg.index.values
# print df_Corr_Mtrx


df_Corr_Mtrx.to_csv('Correlation.csv')


######################################################################
######################################################################
Tarlist = list(db[db.Chen_isTar==1].index)
for tar in Tarlist:
	df_tmp = pd.read_csv(tar+'.csv',index_col = 'Transcription factor')
	Rlist = list(df_tmp.index)
	Pearsonr = []
	for reg in Rlist:
		Pearsonr.append(df_Corr_Mtrx.loc[reg, tar])
	df_tmp['Pearsonr'] = Pearsonr
	# print df_tmp
	df_tmp.to_csv(tar+'_validation.csv')



######################################################################
######################################################################
cmap = sns.diverging_palette(250, 10, sep=20, l=55,  as_cmap=True)
plt.figure()
g_hp = sns.heatmap(df_Corr_Mtrx, 
	robust=True, yticklabels=False, center=0,vmin = -1, vmax = 1, cmap=cmap, annot=False)
# g_hp.set_yticklabels(g_hp.get_yticklabels(), rotation = 0, fontsize = 8)
# g_hp.set_ylabel('TF', fontsize = 8)
# g_hp.set_xlabel('nTF', fontsize = 8)
# g_hp.set_title('Correlation For Targets Of Interest With Potential Regulators', fontsize = 35 )

fig = plt.gcf()
fig.set_size_inches(18.5,10.5)
fig.savefig('Heatmap.pdf',bbox_inches='tight')




















