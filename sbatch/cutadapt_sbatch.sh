#!/bin/bash
#cutadapt script

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH --cpus-per-task=8                 #optional: number of cpus, default is 1
#SBATCH --mem=16GB                        #optional: amount of memory, default is 4GB per cpu
#SBATCH --time=10:00:00 
#SBATCH --job-name=cutadapt            #optional: job name
#SBATCH --output=cutadapt%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=cutadapt%j.err        #optional: file to store stderr from job, %j adds the assigned jobID


/usr/bin/time -v cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o SRR25630297_1_trimmed.fastq -p SRR25630297_2_trimmed.fastq SRR25630297_1.fastq SRR25630297_2.fastq

/usr/bin/time -v cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o SRR25630381_1_trimmed.fastq -p SRR25630381_2_trimmed.fastq SRR25630381_1.fastq SRR25630381_2.fastq



