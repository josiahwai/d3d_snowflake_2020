#!/bin/bash -l

#SBATCH -n 1
#SBATCH -t 00:60:00
#SBATCH -J snowSYN
#SBATCH -o job_snowf-%j.out
#SBATCH -e job_snowf-%j.err
#SBATCH --export=ALL

module purge
module load matlab/r2016b

ls -lan
pwd
echo '<<>> BEGINNING <<>>'
date

matlab -nodisplay < HeatFluxSimulation_Synthetic.m

ls -lan *.mat

echo '<<>> END <<>>'
date
