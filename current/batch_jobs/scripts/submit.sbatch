#!/bin/bash

#SBATCH -n 1
#SBATCH -p general
#SBATCH -t 00:30:00
#SBATCH -J sfd_sim
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --export=ALL


echo '<<>> BEGINNING <<>>'
date

matlab -nodisplay < design_eq_batch.m


date
echo '<<>> END <<>>'