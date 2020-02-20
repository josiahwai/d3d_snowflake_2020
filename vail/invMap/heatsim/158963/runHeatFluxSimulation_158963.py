
import numpy as np
import os

for time in np.arange(3600, 6000, 100):

	os.mkdir(str(time))

	os.system('cp run.sbatch ' + str(time) + '/run.sbatch')
	os.system('cp HeatFluxSimulation_158963.m ' + str(time) + '/HeatFluxSimulation_158963.m')

	f = open(str(time) + '/time.txt', 'w')
	f.write(str(time))
	f.close

	os.system('cd ' + str(time) + '; sbatch run.sbatch')
