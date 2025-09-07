#!/bin/bash
#mapped_unmapped script

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH --cpus-per-task=8                 #optional: number of cpus, default is 1
#SBATCH --time=10:00:00 
#SBATCH --mem=16GB                        #optional: amount of memory, default is 4GB per cpu 
#SBATCH --job-name=mapped_unmapped            #optional: job name
#SBATCH --output=mapped_unmapped%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=mapped_unmapped%j.err        #optional: file to store stderr from job, %j adds the assigned jobID



/usr/bin/time -v python mapped_unmapped.py -f campylo_SRR25630381_picard_deduplicated.sam

/usr/bin/time -v python mapped_unmapped.py -f campylo_SRR25630297_picard_deduplicated.sam


