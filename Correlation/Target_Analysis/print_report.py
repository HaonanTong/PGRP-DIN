# Copy files *_1Parent to Target_Analysis/Summary/
import pandas as pd
import os
import fnmatch
import sys


matches = []
for root, dirnames, filenames in os.walk('.'):
    for filename in fnmatch.filter(filenames,'AT*_report.csv'):
    	print filename
        matches.append(os.path.join(root, filename))

# pd.set_option('display.width', pd.util.terminal.get_terminal_size()[0])


for file in matches:
	db = pd.read_csv(file,index_col='Transcription factor')
	# pd.set_option('display.precision', 2)
	print '------------------------------------------------'
	print
	print db.to_string()
	print
	raw_input("Press to Continue...")