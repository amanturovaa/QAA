#!/bin/bash
#STAR script SRR25630297

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH --cpus-per-task=8                 #optional: number of cpus, default is 1
#SBATCH --time=10:00:00
#SBATCH --mem=32GB                        #optional: amount of memory, default is 4GB per cpu 
#SBATCH --job-name=STAR_SRR25630297           #optional: job name
#SBATCH --output=STAR_SRR25630297%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=STAR_SRR25630297%j.err        #optional: file to store stderr from job, %j adds the assigned jobID


#align reads
/usr/bin/time -v STAR \
    --runThreadN 8 \
    --runMode alignReads \
    --outFilterMultimapNmax 3 \
    --outSAMunmapped Within KeepPairs \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --readFilesCommand zcat \
    --readFilesIn SRR25630297_1_paired.fastq.gz SRR25630297_2_paired.fastq.gz \
    --genomeDir campylomormyrus_STAR \
    --outFileNamePrefix campylo_SRR25630297_

#SAM to BAM, sort, index
/usr/bin/time -v samtools view -bS campylo_SRR25630297_Aligned.out.sam > campylo_SRR25630297_Aligned.bam

/usr/bin/time -v samtools sort campylo_SRR25630297_Aligned.bam -o campylo_SRR25630297_Aligned_sorted.bam

/usr/bin/time -v samtools index campylo_SRR25630297_Aligned_sorted.bam