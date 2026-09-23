# Add a replica of the original data to the end of the old
# datafile but with realization # changed from 0 to 1

import subprocess 

fin_name = 'permanent_test_positions.txt'
  
subprocess.run(['cp ' + fin_name + ' test_temp.txt'], shell=True, check=True)

with open(fin_name, 'r') as fin:
	old_file = fin.readlines()

with open(fin_name, 'a') as fout:
	for il, line in enumerate(old_file):
		temp = line.split()
		temp[-1] = '1'
		if il < len(old_file) - 1:
			temp.append('\n')
		fout.write((' ').join(temp)) 
