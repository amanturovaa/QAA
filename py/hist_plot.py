#!/usr/bin/env python
import gzip
import matplotlib.pyplot as plt

def get_read_lengths(filename):
    lengths = []
    with gzip.open(filename, 'rt') as f:
        for i, line in enumerate(f):
            if i % 4 == 1:  
                lengths.append(len(line.strip()))
            if len(lengths) >= 1000000: 
                break
    return lengths

#reading files
r1_297 = get_read_lengths("SRR25630297_1_paired.fastq.gz")
r2_297 = get_read_lengths("SRR25630297_2_paired.fastq.gz")
r1_381 = get_read_lengths("SRR25630381_1_paired.fastq.gz")
r2_381 = get_read_lengths("SRR25630381_2_paired.fastq.gz")


#plot SRR25630297
plt.figure()
plt.hist(r1_297, bins=30, alpha=0.6, color='blue', label='R1')
plt.hist(r2_297, bins=30, alpha=0.6, color='pink', label='R2')
plt.title('SRR25630297')
plt.ylabel('Frequency')
plt.xlabel('Read Length')
plt.legend()
plt.savefig('SRR25630297_read_lengths.png')
plt.show()

#plot SRR25630381
plt.figure()
plt.hist(r1_381, bins=30, alpha=0.6, color='blue', label='R1')
plt.hist(r2_381, bins=30, alpha=0.6, color='pink', label='R2')
plt.title('SRR25630381')
plt.ylabel('Frequency')
plt.xlabel('Read Length')
plt.legend()
plt.savefig('SRR25630381_read_lengths.png')
plt.show()

