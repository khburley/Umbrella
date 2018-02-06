#!/bin/bash
#
# Author: Gaetano Calabro, UCI gcalabro@uci.edu
# Amber pmemd single gpu execution
# Job name:
#----------------
#SBATCH -J "Umbrella"
#----------------

#----------------
#SBATCH -p mf_titanx
#----------------

#----------------
#SBATCH -o out.tex
#----------------

#----------------
#SBATCH -e err.tex
#----------------

# Specifying resources needed for run:
#
#--------------
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH --distribution=block:cyclic
#SBATCH --gres=gpu:1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-core=2
#SBATCH --threads-per-core=2

#--------------

# Informational output
echo "=================================== SLURM JOB ==================================="
echo
echo "The job will be started on the following node(s):"
echo $SLURM_JOB_NODELIST
echo
echo "Slurm User:         $SLURM_JOB_USER"
echo "Run Directory:      $(pwd)"
echo "Job ID:             $SLURM_JOB_ID"
echo "Job Name:           $SLURM_JOB_NAME"
echo "Partition:          $SLURM_JOB_PARTITION"
echo "Number of nodes:    $SLURM_JOB_NUM_NODES"
echo "Number of tasks:    $SLURM_NTASKS"
echo "Submitted From:     $SLURM_SUBMIT_HOST"
echo "Submit directory:   $SLURM_SUBMIT_DIR"
echo "=================================== SLURM JOB ==================================="
echo

export MODULEPATH=/home/gcalabro/local/modulefiles:$MODULEPATH
source activate py35

cd $SLURM_SUBMIT_DIR

echo 'Working Directory:'
pwd

 
date

python restrOpenMM_GPU.py -i mol2files/ -f gaff -g mol2files/ > minGAFF.out

date

