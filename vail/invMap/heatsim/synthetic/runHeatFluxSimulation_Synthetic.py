
import numpy as np
import os

for eq in range(59,64):

	if eq <= 9:
		eqstr = '0' + str(eq)
	else:
		eqstr = str(eq)

	os.mkdir('eq' + eqstr)

	os.system('cp run.sbatch ' + 'eq' + eqstr + '/run.sbatch')
	os.system('cp HeatFluxSimulation_Synthetic.m ' + 'eq' + eqstr + '/HeatFluxSimulation_Synthetic.m')

	f = open('eq' + eqstr + '/option.txt', 'w')
	f.write(eqstr)
	f.close

	os.system('cd eq' + eqstr + '; sbatch run.sbatch')
