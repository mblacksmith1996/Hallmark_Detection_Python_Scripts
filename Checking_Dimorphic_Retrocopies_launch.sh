#!/bin/bash
#SBATCH --mail-user=blacksmi@umich.edu
#SBATCH --mail-type=FAIL,ARRAY_TASKS
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=18G
#SBATCH --time=4:00:00
#SBATCH --account=jmkidd1
#SBATCH --partition=debug
#SBATCH --output=logs/%x-%A_%a.out.log
#SBATCH --export=ALL

python  Checking_Dimorphic_Retrocopies.py > CDR.log
