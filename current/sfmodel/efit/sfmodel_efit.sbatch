#!/bin/bash

#SBATCH -n 1
#SBATCH -p general
#SBATCH -t 3:00:00
#SBATCH -J sfmodel
#SBATCH --mem=8000
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --export=ALL


echo '<<>> BEGINNING <<>>'
date

matlab -nodisplay < sfmodel_efit.m


date
echo '<<>> END <<>>'