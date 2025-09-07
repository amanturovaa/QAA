#!/bin/bash
#picard script

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH --cpus-per-task=8                 #optional: number of cpus, default is 1
#SBATCH --time=10:00:00 
#SBATCH --job-name=picard            #optional: job name
#SBATCH --output=picard%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=picard%j.err        #optional: file to store stderr from job, %j adds the assigned jobID


/usr/bin/time -v picard MarkDuplicates \
    INPUT=campylo_SRR25630297_Aligned_sorted.bam \
    OUTPUT=campylo_SRR25630297_Aligned_sorted_picard.sam \
    METRICS_FILE=campylo_SRR25630297_sorted_picard.metrics \
    REMOVE_DUPLICATES=TRUE \
    VALIDATION_STRINGENCY=LENIENT





/usr/bin/time -v picard MarkDuplicates \
    INPUT=campylo_SRR25630381_Aligned_sorted.bam \
    OUTPUT=campylo_SRR25630381_Aligned_sorted_picard.sam \
    METRICS_FILE=campylo_SRR25630381_sorted_picard.metrics \
    REMOVE_DUPLICATES=TRUE \
    VALIDATION_STRINGENCY=LENIENT





