import os
import shutil



dirpath = '/u/jwai/ITER/Ip/scans4/001'         #MODIFY dirpath AS APPROPRIATE

sbatch_path = '/u/jwai/ITER/python_scripts/'
os.system('dos2unix submit.sbatch')

if os.path.exists(dirpath + '/output'):
  shutil.rmtree(dirpath + '/output')
os.mkdir(dirpath + '/output')
  

os.chdir(dirpath)
os.system('sbatch ' + sbatch_path + 'submit.sbatch')

