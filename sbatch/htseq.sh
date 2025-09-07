#!/bin/bash
#htseq for SRR25630381 and SRR25630297

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH --cpus-per-task=8                 #optional: number of cpus, default is 1
#SBATCH --time=10:00:00 
#SBATCH --mem=32GB                        #optional: amount of memory, default is 4GB per cpu 
#SBATCH --job-name=htseq          #optional: job name
#SBATCH --output=htseq_SRR25630381%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=htseq_SRR25630381%j.err        #optional: file to store stderr from job, %j adds the assigned jobID



/usr/bin/time -v htseq-count \
    --format sam \
    --stranded=yes \
    --idattr gene_id \
    campylo_SRR25630381_picard_deduplicated.sam \
    campylomormyrus.gtf > SRR25630381_htseq_stranded_yes_counts.txt

/usr/bin/time -v htseq-count \
    --format sam \
    --stranded=reverse \
    --idattr gene_id \
    campylo_SRR25630381_picard_deduplicated.sam \
    campylomormyrus.gtf > SRR25630381_htseq_stranded_reverse_counts.txt

/usr/bin/time -v htseq-count \
    --format sam \
    --stranded=yes \
    --idattr gene_id \
    campylo_SRR25630297_picard_deduplicated.sam \
    campylomormyrus.gtf > SRR25630297_htseq_stranded_yes_counts.txt

/usr/bin/time -v htseq-count \
    --format sam \
    --stranded=reverse \
    --idattr gene_id \
    campylo_SRR25630297_picard_deduplicated.sam \
    campylomormyrus.gtf > SRR25630297_htseq_stranded_reverse_counts.txt