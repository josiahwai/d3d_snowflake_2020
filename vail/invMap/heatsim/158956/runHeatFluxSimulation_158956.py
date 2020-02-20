
import numpy as np
import os

times = np.arange(3700, 6000, 100)
times = np.delete(times, 19) # Remove 5600 ms

for time in times:

	os.mkdir(str(time))

	os.system('cp run.sbatch ' + str(time) + '/run.sbatch')
	os.system('cp HeatFluxSimulation_158956.m ' + str(time) + '/HeatFluxSimulation_158956.m')

	f = open(str(time) + '/time.txt', 'w')
	f.write(str(time))
	f.close

	os.system('cd ' + str(time) + '; sbatch run.sbatch')
