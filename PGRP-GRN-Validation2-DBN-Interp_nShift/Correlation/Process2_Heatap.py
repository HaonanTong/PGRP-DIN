import seaborn as sns
import pandas as pd
import numpy as np;
import matplotlib.pyplot as plt
# plt.plot([1,2,3,4])
# plt.ylabel('some numbers')
# plt.show()

# uniform_data = np.random.rand(10, 12)
# g_hp = sns.heatmap(uniform_data)
# plt.show()

for tar_atp in range(2,7):
	db = pd.read_csv('Data/PearsonCorr_C2H4_Tar_ATP_'+str(tar_atp)+'.csv',index_col = 'id')
	# print db
	cmap = sns.diverging_palette(250, 10, sep=20, l=55,  as_cmap=True)
	plt.figure()
	g_hp = sns.heatmap(db, robust=True, yticklabels=True, center=0,vmin = -1, vmax = 1, cmap=cmap, annot=True)
	g_hp.set_yticklabels(g_hp.get_yticklabels(), rotation = 0, fontsize = 8)
	g_hp.set_ylabel('TF', fontsize = 8)
	g_hp.set_xlabel('nTF', fontsize = 8)
	g_hp.set_title('Correlation For nTF Actived At Time Point ' + str(tar_atp) + '\nWith Potential Regulators', fontsize = 35 )

	fig = plt.gcf()
	fig.set_size_inches(18.5,10.5)
	fig.savefig('Img/Heatmap_atp_'+str(tar_atp)+'.pdf',bbox_inches='tight')
	# plt.show()
# 	raw_input('Press to continue...')