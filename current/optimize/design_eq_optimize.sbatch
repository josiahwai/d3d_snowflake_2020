#!/bin/bash

#SBATCH -n 1
#SBATCH -p general
#SBATCH -t 01:00:00
#SBATCH -J heatsim
#SBATCH --mem=4000
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --export=ALL


echo '<<>> BEGINNING <<>>'
date

matlab -nodisplay < design_eq_optimize.m


date
echo '<<>> END <<>>'