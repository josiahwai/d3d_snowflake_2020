#!/bin/bash

#SBATCH -n 1
#SBATCH -p general
#SBATCH -t 01:00:00
#SBATCH -J heatsim
#SBATCH --mem=1000
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --export=ALL


echo '<<>> BEGINNING <<>>'
date

matlab -nodisplay < sim_gridsearch.m


date
echo '<<>> END <<>>'