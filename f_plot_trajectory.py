import pandas as pd
import sys
import matplotlib.pyplot as plt


g_Reg = str(sys.argv[1])
g_Tar = str(sys.argv[2])

db = pd.read_csv('myDataFrame.csv',index_col='id')

dfs = db.loc[[g_Reg,g_Tar],'t0':'t6']
dfs = dfs.T
dfs = ( dfs - dfs.mean() ) / dfs.std()



gh = dfs.plot(
	title='Trajectories'
	, lw=2 , legend=True )#, style = 'lightgrey')
# data_nTF.AT1G04310.plot(lw=4,ax=gh, legend = True)
# data_TF.AT5G13910.plot(lw=4,ax=gh, legend = True)

gh.set_xlabel("Time (hour)")
gh.set_ylabel("Normalized Expression Level")
gh.set_xticklabels(['0','0.25','0.5','1','4','12','24'])
# gh.legend(loc = 'best')

# Output Figures
fig = plt.gcf()
fig.set_size_inches(18.5,10.5)
fig.savefig('Img/Trajectories_'+g_Reg+'_'+g_Tar+'.pdf',bbox_inches='tight')