
import numpy as np
import os

times = np.arange(3500, 4600, 100)

for time in times:

	os.mkdir(str(time))

	os.system('cp run.sbatch ' + str(time) + '/run.sbatch')
	os.system('cp HeatFluxSimulation_158960.m ' + str(time) + '/HeatFluxSimulation_158960.m')

	f = open(str(time) + '/time.txt', 'w')
	f.write(str(time))
	f.close

	os.system('cd ' + str(time) + '; sbatch run.sbatch')
