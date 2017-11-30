import pandas as pd
import os

# Create Summary directories for different levels of discretization
directory_BDeu = 'Target_Analysis/Target_Analysis_nlevels' # +nlevels
directory_Desp = 'Target_Analysis/Target_Desp/'
directory = 'Target_Analysis/Summary/'


# List_n_levels = [2,3,5,7,10,15]
List_n_levels = [2,3,4,5]
for n_levels in List_n_levels:
	try:
	    os.stat(directory+"n_levels_"+str(n_levels))
	except:
	    os.mkdir(directory+"n_levels_"+str(n_levels))   



# recursively walk a directory and fnmatch.filter
# to match against a simple expression:
import fnmatch
import shutil
# Remove Files
# for root, dirnames, filenames in os.walk(directory):
#     for filename in fnmatch.filter(filenames, '*_1Parent.csv'):
#         os.remove(os.path.join(root, filename))

# Copy files *_1Parent to Target_Analysis/Summary/
# matches = []
# for root, dirnames, filenames in os.walk('.'):
#     for filename in fnmatch.filter(filenames, '*_1Parent.csv'):
#         matches.append(os.path.join(root, filename))

# for file in matches:
# 	try:
# 		shutil.copy2(file,directory)
# 	except IOError:
# 		os.chmod(file, 777)
# 		shutil.move(file,directory)

def f_getPosition(db, bdeu):
	l = sum(db.BDeu>bdeu)
	u = sum(db.BDeu>=bdeu)
	return l,u



db = pd.read_csv('../myDB_Validation.csv',index_col='id')
db_DETarList = db[db['Chen_isTar'] == 1]

DETarList = list(db_DETarList.index)

for n_levels in List_n_levels:
	for x in DETarList:
		df_target_tmp = pd.read_csv(directory_Desp + x + '.csv', index_col='Transcription factor')
		df_score_tmp = pd.read_csv(directory_BDeu+str(n_levels)+'/'+x+'_1Parent.csv',index_col='Source')
		BDeu_lower_index = [];
		BDeu_upper_index = [];
		for i in range(0,len(df_score_tmp)):
			bdeu = df_score_tmp.BDeu[i]
			[l, u] = f_getPosition(df_score_tmp, bdeu)
			BDeu_lower_index.append(l)
			BDeu_upper_index.append(u)

		df_score_tmp['BDeu_lower_index'] = BDeu_lower_index
		df_score_tmp['BDeu_upper_index'] = BDeu_upper_index

		df_score_tmp = df_score_tmp.loc[:,'BDeu':'BDeu_upper_index']

		# df_score_tmp = df_score_tmp.BDeu
		result = pd.concat([df_target_tmp, df_score_tmp], axis=1, join_axes=[df_target_tmp.index])
		# print result
		result.to_csv(directory+"n_levels_"+str(n_levels)+'/'+x+'_report.csv')
