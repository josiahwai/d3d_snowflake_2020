import os
import shutil
import numpy as np

#MODIFY dirpath AS APPROPRIATE
shot = 155355
times = np.arange(2000,2200,100)
fn = '/design_eq_batch.m'

shot_dir = '/u/jwai/d3d_snowflake_2020/current/batch_jobs/' + str(shot)
if os.path.exists(shot_dir):
    shutil.rmtree(shot_dir)
os.mkdir(shot_dir)

scripts_dir = '/u/jwai/d3d_snowflake_2020/current/batch_jobs/scripts/'
os.system('dos2unix submit.sbatch')


for time in times:
    os.chdir(shot_dir)
    time_dir = shot_dir + '/' + str(time)

    # clean up folder
    if os.path.exists(time_dir):
        shutil.rmtree(time_dir)
    os.mkdir(time_dir)
    
    # copy job to local dir
    shutil.copy2(scripts_dir + fn, time_dir + fn)
    
    # pass arguments to text file    
    os.chdir(time_dir)
    f = open('args.txt', 'w')
    txt = str(shot) +  ',' + str(time)
    f.write(txt)
    f.close()
    
    # run job
    os.system('sbatch ' + scripts_dir + 'submit.sbatch')
    







