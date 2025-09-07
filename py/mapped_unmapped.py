#!/usr/bin/env python
import argparse 

def parse_sam(sam_file): 
    mapped_count = 0 #initialize the mapped and unmapped reads
    unmapped_count = 0
    
    with open(sam_file, 'r') as file: #opening the sam file and loop through file
        for line in file: #starting my iteration
            if line.startswith('@'): #select line that start with the @ and to then skip them since those are headers
                continue
            
            no_newline_tab = line.strip().split('\t')
            flag = int(no_newline_tab[1]) #using 1 here to grab the FLAG out of sam
            
            if (flag & 256) != 256: #not counting the secondary reads
                if (flag & 4) != 4: #checking the flags
                    mapped_count += 1 #if the flag value is not 4, then it is mapped and value is +1
                else: #if the flag value is 4, then it is unmapped and value is +1
                    unmapped_count += 1
    
    print(f"Mapped reads: {mapped_count}") #printing out the mapped or unmapped reads 
    print(f"Unmapped reads: {unmapped_count}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--filename", required=True)
    args = parser.parse_args()
    parse_sam(args.filename) 