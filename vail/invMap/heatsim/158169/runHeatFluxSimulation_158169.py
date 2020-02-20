
import numpy as np
import os

for time in np.arange(3800, 6000, 100):

	os.mkdir(str(time))

	os.system('cp run.sbatch ' + str(time) + '/run.sbatch')
	os.system('cp HeatFluxSimulation_158169.m ' + str(time) + '/HeatFluxSimulation_158169.m')

	f = open(str(time) + '/time.txt', 'w')
	f.write(str(time))
	f.close

	os.system('cd ' + str(time) + '; sbatch run.sbatch')
