#!/bin/bash

#SBATCH -n 1
#SBATCH -p general
#SBATCH -t 00:30:00
#SBATCH -J heatsim
#SBATCH --mem=2000
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --export=ALL


echo '<<>> BEGINNING <<>>'
date

matlab -nodisplay < sfd_costfun.m


date
echo '<<>> END <<>>'