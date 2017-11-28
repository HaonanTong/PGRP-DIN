# -------------------------------------------
# Haonan Tong
# -------------------------------------------
import pandas as pd 
import numpy as np 
from scipy.stats import pearsonr
from scipy.interpolate import interp1d

# Linearly interpolate time course data and then
# Use Pearson correlation analyze correlation of target activated at a time point
# and potential regulator activated by (at and previous to) the time
# Output Pearson correlation matrix would look like this
# 	G1	G2	G3 ...
# TF1
# TF2
# TF3
# ...
# For nTF activated at 1 - 5
def f_interp( data ):
	t = [0, 0.25, 0.5, 1, 4, 12, 24]
	f = interp1d(t,data)
	x_new = np.arange(0, 24.25, 0.25)
	y_interp = f(x_new)
	return pd.Series(y_interp)

for atp_tar in range(2,6):
	# DataFrame for Target
	df_tar = pd.read_csv('Data/log2Expr_C2H4_nTF_AT'+str(atp_tar)+'.csv', index_col = 'id')

	# DataFrame for Potential Regulator
	frames = [];
	for atp_reg in range(1,(atp_tar+1)):
		# DataFrame for Regulator
		df_reg_tmp = pd.read_csv('Data/log2Expr_C2H4_TF_AT'+str(atp_reg)+'.csv', index_col = 'id')
		frames.append(df_reg_tmp)
	df_reg = pd.concat(frames)
	
	# Derive Pearson Correaltion Matrix
	ntar = len(df_tar)
	nreg = len(df_reg)
	Corr_Matrx = np.zeros((nreg,ntar))
	for tar in range(0,ntar):
		for reg in range(0,nreg):
			star = f_interp( df_tar.iloc[tar,:] )
			sreg = f_interp( df_reg.iloc[reg,:] )

			# print star.corr( sreg, method='pearson')
			Corr_Matrx[reg, tar] = star.corr( sreg, method='pearson')
	# print Corr_Matrx
	df_Corr_Mtrx = pd.DataFrame(Corr_Matrx, index=df_reg.index.values, columns=df_tar.index.values )
	df_Corr_Mtrx.index.name = 'id'
	df_Corr_Mtrx.to_csv('Data/PearsonCorr_C2H4_Interp_Tar_ATP_'+str(atp_tar)+'.csv')


# Target AT6
atp_tar = 6
df_tar = pd.read_csv('Data/log2Expr_C2H4_nTF_AT'+str(atp_tar)+'.csv', index_col = 'id')
atp_reg = 5
df_reg = pd.read_csv('Data/log2Expr_C2H4_TF_AT'+str(atp_reg)+'.csv', index_col = 'id')
# Derive Pearson Correaltion Matrix
ntar = len(df_tar)
nreg = len(df_reg)
Corr_Matrx = np.zeros((nreg,ntar))
for tar in range(0,ntar):
	for reg in range(0,nreg):
		star = f_interp( df_tar.iloc[tar,:] )
		sreg = f_interp( df_reg.iloc[reg,:] )

		# print star.corr( sreg, method='pearson')
		Corr_Matrx[reg, tar] = star.corr( sreg, method='pearson')
# print Corr_Matrx
df_Corr_Mtrx = pd.DataFrame(Corr_Matrx, index=df_reg.index.values, columns=df_tar.index.values )
df_Corr_Mtrx.index.name = 'id'
df_Corr_Mtrx.to_csv('Data/PearsonCorr_C2H4_Interp_Tar_ATP_'+str(atp_tar)+'.csv')