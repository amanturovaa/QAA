#!/bin/bash
#gzip script

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH --cpus-per-task=8                 #optional: number of cpus, default is 1
#SBATCH --time=10:00:00 
#SBATCH --job-name=gzip            #optional: job name
#SBATCH --output=gzip%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=gzip%j.err        #optional: file to store stderr from job, %j adds the assigned jobID



gzip CcoxCrh_comrhy60_E0_6cm_1_read1.fastq
gzip CcoxCrh_comrhy60_E0_6cm_1_read2.fastq
gzip CcoxCrh_comrhy112_EO_adult_2_read1.fastq
gzip CcoxCrh_comrhy112_EO_adult_2_read2.fastq



