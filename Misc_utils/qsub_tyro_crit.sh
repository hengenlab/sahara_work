#!/bin/bash
#PBS -S /bin/bash
#PBS -N spkint_job_1
#PBS -V
#PBS -l nodes=1:ppn=4:gpus=1:exclusive_process,walltime=24:00:00


python3 big_crit.py '' '' True True