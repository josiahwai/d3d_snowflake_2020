import os
import shutil
import numpy as np


shot = 155334
shot_dir = '/u/jwai/d3d_snowflake_2020/current/batch_jobs/' + str(shot)

if os.path.exists(time_dir):
        shutil.rmtree(time_dir)
os.mkdir(time_dir)