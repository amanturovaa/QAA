#!/bin/bash
#trimmomatic script

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH --cpus-per-task=8                 #optional: number of cpus, default is 1
#SBATCH --mem=16GB                        #optional: amount of memory, default is 4GB per cpu
#SBATCH --time=10:00:00 
#SBATCH --job-name=trimmomatic            #optional: job name
#SBATCH --output=trimmomatic%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=trimmomatic%j.err        #optional: file to store stderr from job, %j adds the assigned jobID


/usr/bin/time -v trimmomatic PE \
    SRR25630297_1_trimmed.fastq SRR25630297_2_trimmed.fastq \
    SRR25630297_1_paired.fastq.gz SRR25630297_1_unpaired.fastq.gz \
    SRR25630297_2_paired.fastq.gz SRR25630297_2_unpaired.fastq.gz \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:35


/usr/bin/time -v trimmomatic PE \
    SRR25630381_1_trimmed.fastq SRR25630381_2_trimmed.fastq \
    SRR25630381_1_paired.fastq.gz SRR25630381_1_unpaired.fastq.gz \
    SRR25630381_2_paired.fastq.gz SRR25630381_2_unpaired.fastq.gz \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:35