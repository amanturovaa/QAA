#!/bin/bash
#script to run python quality score dist 

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH --cpus-per-task=8                 #optional: number of cpus, default is 1
#SBATCH --mem=16GB                        #optional: amount of memory, default is 4GB per cpu
#SBATCH --time=10:00:00 
#SBATCH --job-name=qscore            #optional: job name
#SBATCH --output=qscore%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=qscore%j.err        #optional: file to store stderr from job, %j adds the assigned jobID


/usr/bin/time -v python qscore_dist.py -f SRR25630297_1.fastq -l SRR25630297_R1

/usr/bin/time -v python qscore_dist.py -f SRR25630297_2.fastq -l SRR25630297_R2

/usr/bin/time -v python qscore_dist.py -f SRR25630381_1.fastq -l SRR25630381_R1

/usr/bin/time -v python qscore_dist.py -f SRR25630381_2.fastq -l SRR25630381_R2