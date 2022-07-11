#!/bin/bash
#SBATCH --job-name="worker"
#SBATCH --partition="secondary-fdr"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH -t 4:00:00
#SBATCH --output=scratch/logs.txt
huey_consumer.py tasks.remote_tasks.huey -n -w 4 -k process