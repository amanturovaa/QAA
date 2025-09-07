#!/usr/bin/env python
import argparse 
import matplotlib.pyplot as plt 

def convert_phred(letter): #convert ascii character to phred score 33
    return ord(letter) - 33

def process_fastq(file): #process fastq file;calc mean quality scores
    qsums = [0] * 151 #initialize running sums for each position up to 150 bp
    qcounts = [0] * 151 #initialize counts for each position to calc means
    
    with open(file, "r") as fh: #open fastq file for reading
        for i, line in enumerate(fh): #loop through each line in file
            if i % 4 == 3: #quality lines are every 4th line starting at 3
                line = line.strip() #remove newline characters
                for pos, letter in enumerate(line[:150]): #iterate through each quality character up to position 150
                    q = convert_phred(letter) #convert quality character to phred score
                    qsums[pos] += q #add quality score to running sum for this position
                    qcounts[pos] += 1 #increment count for this position
    
    means = [sum_val / count if count != 0 else 0 for sum_val, count in zip(qsums, qcounts)] #calculate means using running sums
    return means

def plot_quality(means, label): #function to create quality score plot
    plt.figure(figsize=(12, 6)) 
    plt.plot(range(len(means)), means, linewidth=2) #plot mean quality scores
    plt.xlabel("Base Position")
    plt.ylabel("Mean Quality Score")
    plt.title(f"Mean Quality Scores - {label}") 
    plt.xlim(0, 150)
    plt.ylim(0, 42) 
    plt.savefig(f"{label}_quality_distribution.png", dpi=300, bbox_inches='tight') #save plot as png
    plt.show() 


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate mean quality score plots per base using running sum strategy.")
    parser.add_argument("-f", "--file", required=True)
    parser.add_argument("-l", "--label", required=True) 
    args = parser.parse_args() 
    
    print(f"{args.file}") 
    means = process_fastq(args.file) 
    plot_quality(means, args.label) 