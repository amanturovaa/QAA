Lab notebook PS2

---
title: "PS2_template"
output: html_document
date: "2025-04-25"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# RNA-seq Quality Assessment Assignment - Bi 623 (Summer 2025 Assignment/PS 2)

## Overall assignment:

In this assignment, you will process electric organ and/or skeletal muscle RNA-seq reads for a future differential gene expression analysis. We will be completing the differential gene expression analysis in our last bioinformatics assignment of this class. You will learn how to use existing tools for quality assessment and read trimming, compare quality assessments to those created by your own software, and how to align and count reads. Additionally, you will learn how to summarize important information in a high-level report. You should create a cohesive, well written report for your "PI" about what you've learned about/from your data.

**This template is provided as reference for instructions. Files with specific naming conventions are requested to be turned in at the end of this problem set. You can use this template to gather notes while completing this assignment.** Be sure to upload all relevant materials by the deadline and **double check** to be sure that your offline repository is up-to-date with your online repository. Answers to questions should be included in your final, high-level, report as a `pdf`. This pdf should be generated using Rmarkdown and submitted to Canvas as well as GitHub. Be sure to keep a well-organized, detailed lab notebook!


### Dataset:

Each of you will be working with 2 RNA-seq files from two different electric fish studies (PRJNA1005245 and PRJNA1005244). The methods for the PRJNA1005244 dataset are [published](https://doi.org/10.1093/molbev/msae021) and the methods for the PRJNA1005245 dataset are written in the third chapter of a [thesis](https://canvas.uoregon.edu/courses/266187/files/22059308?module_item_id=5380118). For all steps below, process the two libraries separately. SRR assignments are here: ```/projects/bgmp/shared/Bi623/PS2/QAA_data_Assignments.txt```. If you have time, consider claiming and processing additional RNA-seq raw sequencing files via this [google doc](https://docs.google.com/document/d/1vEmVEzUaTjbDF4JyNsWH-wFpi8dm4wkcvWgSoYZzoCY/edit?usp=sharing). Although this is not extra credit, it will make our downstream RNA-seq analysis more interesting and your classmates will appreciate your efforts.

You are responsible for downloading this data from NCBI SRA, dumping into FASTQ files, and zipping those files (check ICA1 for a refresher). We are processing this data for use in a future assignment, so please keep your files well organized. Finally, rename the files to the convention
Species_sample_tissue_age/size_sample#_readnumber.fastq.gz.

SpeciesSmaller,Sample,Age/Size,Tissue,Sample#,FastqFileName,YourSubmittedHTSeqCountsFileName

CcoxCrh_comrhy112_EO_adult_2_read1.fastq.gz #381
CcoxCrh_comrhy112_EO_adult_2_read2.fastq.gz #381

CcoxCrh_comrhy60_E0_6cm_1_read1.fastq.gz #297
CcoxCrh_comrhy60_E0_6cm_1_read2.fastq.gz #297


**Reminder: This template file IS not your final product; however, it gives you a space to record all of the necessary information for your final report.**

```{bash, eval=FALSE}
## Download your data

srun -A bgmp -p bgmp --time=0-8:00:00 --cpus-per-task=6 --pty bash

conda install bioconda::sra-tools

prefetch SRR25630297
prefetch SRR25630381

fasterq-dump SRR25630297
fasterq-dump SRR25630381
```


SRR25630297     Anna                                                                                
SRR25630381     Anna 


fasterq-dump SRR25630297/
spots read      : 58,398,524
reads read      : 116,797,048
reads written   : 116,797,048


fasterq-dump SRR25630381
spots read      : 6,776,454
reads read      : 13,552,908
reads written   : 13,552,908

## Part 1 – Read quality score distributions

1. Create a new conda environment called `QAA` and install `FastQC`, `cutadapt`, and `Trimmomatic`. Google around if you need a refresher on how to create conda environments. Recommend doing this in an interactive session, not the login node! Record details of how you created this environment in your lab notebook! Make sure you check your installation with:
   - `fastqc --version` (should be 0.12.1)  

[Record details on how you made the conda environment]

```{bash, eval=FALSE}
conda create --name QAA
conda activate QAA
conda install fastqc
conda install bioconda::cutadapt
conda install trimmomatic
fastqc --version
conda install python=3.10
conda install bioconda::cutadapt=5.0
```

2. Using `FastQC` via the command line on Talapas, produce plots of the per-base quality score distributions for R1 and R2 reads. Also, produce plots of the per-base N content, and comment on whether or not they are consistent with the quality score plots.

[Include FastQC commands, plots of per-base N content, comments on consistency with quality score plots]

```{r}
fastqc SRR25630297_1.fastq SRR25630297_2.fastq
fastqc SRR25630381_1.fastq SRR25630381_2.fastq


knitr::include_graphics("SRR25630297_1.fastq.png")
knitr::include_graphics("SRR25630297_2.fastq.png")
knitr::include_graphics("SRR25630381_1.fastq.png")
knitr::include_graphics("SRR25630381_2.fastq.png")

```
The per sequence quality scores look the same between the 4 fastq files, with the average quality per read being being a Phred score of 36. 
The per base sequence quality looks similar between SRR25630297_1.fastq.png, SRR25630297_2.fastq.png, SRR25630381_1.fastq.png, but SRR25630381_2.fastq.png but the last third of the reads are in yellow and last position in the reads is in the red.




3. Run your quality score plotting script from your Demultiplexing assignment in Bi622. (Make sure you're using the "running sum" strategy!!) Describe how the `FastQC` quality score distribution plots compare to your own. If different, propose an explanation. Also, does the run time differ? Mem/CPU usage? If so, why?

sbatch qscore_dist.sh
The mean quality scores look the same between the fastqc and my script from the demultiplexing assignment, but the fastqc ran much quicker than the python demultiplexing script. The larger of the fastq files took 20 min to run with the demultiplexing script while the fastqc took less than 20 min to run all 4 files. 


	Command being timed: "python qscore_dist.py -f SRR25630297_1.fastq -l SRR25630297_R1"
	User time (seconds): 1237.05
	System time (seconds): 9.90
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 20:55.60
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 109312
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 23384
	Voluntary context switches: 2957
	Involuntary context switches: 335
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0

[Include quality score distribution plot, a comparison to FastQC, comments on any differences between your quality score plotting and FastQC]

```{r}
knitr::include_graphics("SRR25630297_1.fastq_py.png")
knitr::include_graphics("SRR25630297_2.fastq_py.png")
knitr::include_graphics("SRR25630381_1.fastq_py.png")
knitr::include_graphics("SRR25630381_2.fastq_py.png")

```





4. Comment on the overall data quality of your two libraries. Go beyond per-base qscore distributions. Examine the `FastQC` [documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/) for guidance on interpreting results and planning next steps. Make and justify a recommendation on whether these data are of high enough quality to use for further analysis. 

[Include comments on data quality and recommendation on whether this can be used for further analysis]
Looking at the per tile sequence quality for all the plots, which shows if there is a loss of quality associated with a part of the flow cell. If the plot is completely blue, that indicates a good quality plot, if there are warmer colors, that indicates poorer quality. SRR25630297 overall had a "good" looking plot. Read1 had a couple lines of orange, and Read2 was completely blue. SRR25630381 had streaks of warmer colors around the 125 to 150 bp range at 2109, 2160, and 2607. Read2 had more of the yellows as opposed to greens like Read2 indicating poorer quality from the flow cell in Read2. This could be caused by inhibitors from the library, incorrect insert size, or incorrect library quality or quantity. SRR25630381 Read2 shows some adapter content starting from reads 105 to 135, so there are Illumina Universal Adapters in 5-10% of the reads. The first SRR25630297 is good for further analysis. 


## Part 2 – Adaptor trimming comparison

5.  If you haven't already in your QAA environment, install `Cutadapt` and `Trimmomatic`. Check your installations with:
    - `cutadapt --version` (should be 5.0) sbatch this
    - `trimmomatic -version` (should be 0.39)

[Record details on install (if happened here) and/or version checking]

conda install bioconda::cutadapt=5.0
cutadapt --version
5.0

conda install trimmomatic=0.39
trimmomatic -version
0.39



6. Using `Cutadapt`, properly trim adapter sequences from your assigned files. Be sure to read how to use `Cutadapt`. Use default settings. What proportion of reads (both R1 and R2) were trimmed?

    <details>
    <summary>Try to determine what the adapters are on your own. If you cannot (or if you do, and want to confirm), click here to see the actual adapter sequences used.</summary>
  
    R1: `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`
    
    R2: `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`
    </details>

    - *Sanity check*: Use your Unix skills to search for the adapter sequences in your datasets and confirm the expected sequence orientations. Report the commands you used, the reasoning behind them, and how you confirmed the adapter sequences.
  
SRR25630297
=== Summary ===
  Total read pairs processed:         58,398,524
  Read 1 with adapter:               3,372,644 (5.8%)
  Read 2 with adapter:               3,758,353 (6.4%)
Pairs written (passing filters):    58,398,524 (100.0%)

Total basepairs processed: 17,519,557,200 bp
  Read 1: 8,759,778,600 bp
  Read 2: 8,759,778,600 bp
Total written (filtered):  17,420,247,099 bp (99.4%)
  Read 1: 8,710,933,998 bp
  Read 2: 8,709,313,101 bp


=== First read: Adapter 1 ===

Sequence: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA; Type: regular 3'; Length: 33; Trimmed: 3372644 times

Minimum overlap: 3
No. of allowed errors:
1-9 bp: 0; 10-19 bp: 1; 20-29 bp: 2; 30-33 bp: 3

Bases preceding removed adapters:
  A: 12.8%
  C: 34.0%
  G: 35.4%
  T: 17.7%
  none/other: 0.0%





=== Second read: Adapter 2 ===

Sequence: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT; Type: regular 3'; Length: 33; Trimmed: 3758353 times

Minimum overlap: 3
No. of allowed errors:
1-9 bp: 0; 10-19 bp: 1; 20-29 bp: 2; 30-33 bp: 3

Bases preceding removed adapters:
  A: 14.9%
  C: 33.0%
  G: 38.3%
  T: 13.7%
  none/other: 0.0%




SRR25630381
=== Summary === 

Total read pairs processed:          6,776,454
  Read 1 with adapter:                 971,693 (14.3%)
  Read 2 with adapter:                 969,365 (14.3%)
Pairs written (passing filters):     6,776,454 (100.0%)

Total basepairs processed: 2,032,936,200 bp
  Read 1: 1,016,468,100 bp
  Read 2: 1,016,468,100 bp
Total written (filtered):  1,992,831,463 bp (98.0%)
  Read 1:   995,754,215 bp
  Read 2:   997,077,248 bp

=== First read: Adapter 1 ===

Sequence: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA; Type: regular 3'; Length: 33; Trimmed: 971693 times

Minimum overlap: 3
No. of allowed errors:
1-9 bp: 0; 10-19 bp: 1; 20-29 bp: 2; 30-33 bp: 3

Bases preceding removed adapters:
  A: 10.9%
  C: 27.3%
  G: 45.5%
  T: 16.3%
  none/other: 0.0%



=== Second read: Adapter 2 ===

Sequence: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT; Type: regular 3'; Length: 33; Trimmed: 969365 times

Minimum overlap: 3
No. of allowed errors:
1-9 bp: 0; 10-19 bp: 1; 20-29 bp: 2; 30-33 bp: 3

Bases preceding removed adapters:
  A: 12.7%
  C: 28.5%
  G: 46.2%
  T: 12.6%
  none/other: 0.0%


[Include commands and report out the proportion of reads trimmed]    

```{bash, eval=FALSE}

sbatch cutadapt_sbatch.sh

grep "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" SRR25630297_1_trimmed.fastq | wc
grep "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" SRR25630297_2_trimmed.fastq | wc
grep "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" SRR25630381_1_trimmed.fastq | wc
```


grep commands showed 0 instances of the adaptor sequences appearing, which means they were successfully trimmed. 




7. Use `Trimmomatic` to quality trim your reads. Specify the following, **in this order**:
    - LEADING: quality of 3
    - TRAILING: quality of 3
    - SLIDING WINDOW: window size of 5 and required quality of 15
    - MINLENGTH: 35 bases

    Be sure to output compressed files and clear out all intermediate files.
    
    
 chmod 755 trimmomatic_sbatch.sh
```{bash, eval=FALSE}
sbatch trimmomatic_sbatch.sh
zcat SRR25630381_2_paired.fastq.gz | head

```
@SRR25630381.1 A00821:326:H5KKHDSXY:3:1101:15646:1000 length=150
TGTAAGTCCAGGTAAATGGTTTAAGAGGTCAAAAAAGACATAGCATGTGAGAAATCTGTATTTCTTTTCCAGCATAACATACCATAAAGTGATACAGAGATAAAGAAAAAAAATTGGATCGATATGGGAAGCCATCCCAGCAAAGTA
+
FFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFF:FFFFF,FFFFFFFFFFFFFFF:FFFFF
@SRR25630381.2 A00821:326:H5KKHDSXY:3:1101:27850:1000 length=150
ATTAGAGCCAACCCGTCTCTGTGGCAAAAGAGTGGGAAGATCTTCGAGTAGAGGTGATAAACCTACCGAACCTAGTGATAGCTGGTTGCTTAGGAAATGGATATTAGTTCAGCTTACTGCCATTCTCAGATCAAAATAATAAGGACCAAC
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:,FFFFFFFFFFFFFF::FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFF:FFFFFFFFF:FFFFFFFFFFFF:FFFFFF:FFFFF
@SRR25630381.3 A00821:326:H5KKHDSXY:3:1101:30273:1000 length=150
TGGCATGATGCCAATTCAGAAGAAACCAAAGTCGTCCAACATGCTGTTGGTCTATCAAACCGATTCCAGCAACCAAAGGAGCATTCACCTCGCAGTACCCTGTCTCCTTCAGATTCTACTGCCACAAAGGAACCACACGGTTATTTCAGG





8. Plot the trimmed read length distributions for both paired R1 and paired R2 reads (on the same plot - yes, you will have to use Python or R to plot this. See ICA4 from Bi621). You can produce 2 different plots for your 2 different RNA-seq samples. There are a number of ways you could possibly do this. One useful thing your plot should show, for example, is whether R1s are trimmed more extensively than R2s, or vice versa. Comment on whether you expect R1s and R2s to be adapter-trimmed at different rates and why.

[Include your plot and comment on R1/R2 adapter trimming]
The R1/R2 for both files show minimal differences on the adapter trimming. Read2 shows slightly more, but in general, R2 is expected to have a higher rate of adapter trimming compared to R1. Differences can be attributed to R2 being read after R1 on illumina, allowing the sample to degrade slighly and R2 reads can have lower quality scores towards the 3' ends.



```{R, eval=TRUE}
knitr::include_graphics("SRR25630297.png")
knitr::include_graphics("SRR25630381.png")
```

9. Bonus - Run `FastQC` on your trimmed data. Comment on differences you observe between the trimmed and untrimmed data. Include any figures needed to support your conclusions.

[Include command, comments on differences, and plot/s]

```{bash, eval=FALSE}

```

## Part 3 – Alignment and strand-specificity
10. Install additional software for alignment and counting of RNA-seq reads. In your QAA environment, use conda to install:
    - Star
    - Picard
    - Samtools
    - NumPy
    - Matplotlib
    - HTSeq

[Record details on how you installed these packages]

```{bash, eval=FALSE}
conda install star
conda install picard=2.18
conda install samtools
conda install numpy
mamba install -c bioconda htseq
```

11. Download the publicly available *Campylomormyrus compressirostris* genome fasta and gff file from [Dryad](https://datadryad.org/dataset/doi:10.5061/dryad.c59zw3rcj) and generate an alignment database from it. If the download fails, the files are available `/projects/bgmp/shared/Bi623/PS2/campylomormyrus.fasta`, `/projects/bgmp/shared/Bi623/PS2/campylomormyrus.gff`. Align the reads to your *C. compressirostris* database using a splice-aware aligner. Use the settings specified in PS8 from Bi621. 

  > [!IMPORTANT]
  > You will need to use gene models to perform splice-aware alignment, see PS8 from Bi621. You may need to convert the gff file into a gtf file for this to work successfully.

[Record details on how you downloaded the genome, prepared the dataset for alignment, and commands for generating the alignment database and aligning reads]
    
    
```{bash, eval=FALSE}
conda install bioconda::gffread
chmod 755 STAR_sbatch.sh
sbatch STAR_sbatch.sh
```

12. Remove PCR duplicates using [Picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates). You may need to sort your reads with `samtools` before running Picard. 
- Use the following for running picard: picard MarkDuplicates INPUT=[FILE] OUTPUT=[FILE] METRICS_FILE=[FILENAME].metrics REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT
 

```{bash, eval=FALSE}
chmod 755 picard_sbatch.sh 
sbatch picard_sbatch.sh 
wc campylo_Aligned_sorted_picard.bam 
samtools view -h campylo_SRR25630297_Aligned_sorted_picard.bam > campylo_SRR25630297_picard_deduplicated.sam

```
  1517044   8831425 432406554 campylo_Aligned_sorted_picard.bam


13. Using your script from PS8 in Bi621, report the number of mapped and unmapped reads from each of your 2 SAM files post deduplication with picard. Make sure that your script is looking at the bitwise flag to determine if reads are primary or secondary mapping (update/fix your script if necessary).


python mapped_unmapped.py -f campylo_SRR25630297_picard_deduplicated.sam
[Include the number of mapped and unmapped reads from both same files]

SRR25630381
Mapped reads: 7580940
Unmapped reads: 751006

SRR25630297
Mapped reads: 36643440
Unmapped reads: 7255102

SRR25630297: 83.6% mapped (36,643,440/43,898,542)
SRR25630381: 91.0% mapped (7,580,940/8,331,946)



```{r}
sample <- c("SRR25630381", "SRR25630297")
mapped <- c(7580940, 36643440)
unmapped <- c(751006, 7255102)
total <- mapped + unmapped

mapped_unmapped_results <- data.frame(
  Sample = sample,
  Mapped_Reads = mapped,
  Unmapped_Reads = unmapped,
  Total_Reads = total)

```


14. Count deduplicated reads that map to features using `htseq-count`. You should run htseq-count twice: once with `--stranded=yes` and again with `--stranded=reverse`. Use default parameters otherwise. You may need to use the `-i` parameter for this run.

```{bash, eval=FALSE}
sbatch htseq.sh
```

15. Demonstrate convincingly whether or not the data are from "strand-specific" RNA-Seq libraries **and** which `stranded=` parameter should you use for counting your reads for a future differential gene expression analyses. Include any commands/scripts used. Briefly describe your evidence, using quantitative statements (e.g. "I propose that these data are/are not strand-specific, because X% of the reads are y, as opposed to z."). This [kit](https://www.revvity.com/product/nex-rapid-dir-rna-seq-kit-2-0-8rxn-nova-5198-01) was used during library preparation. This [paper](https://academic.oup.com/bfg/article/19/5-6/339/5837822) may provide helpful information.

  > [!TIP]
  > Recall ICA4 from Bi621.

tail -10 SRR25630381_htseq_stranded_yes_counts.txt
snap_masked-ptg003120l-processed-gene-0.9       0
snap_masked-ptg003128l-processed-gene-0.5       0
snap_masked-ptg003215l-processed-gene-0.5       0
snap_masked-ptg003293l-processed-gene-0.3       0
snap_masked-ptg003325l-processed-gene-0.4       0
__no_feature    6616210
__ambiguous     313
__too_low_aQual 0
__not_aligned   391801
__alignment_not_unique  235423

tail -10 SRR25630381_htseq_stranded_reverse_counts.txt
snap_masked-ptg003120l-processed-gene-0.9       0
snap_masked-ptg003128l-processed-gene-0.5       61
snap_masked-ptg003215l-processed-gene-0.5       15
snap_masked-ptg003293l-processed-gene-0.3       0
snap_masked-ptg003325l-processed-gene-0.4       12
__no_feature    2911264
__ambiguous     22964
__too_low_aQual 0
__not_aligned   391801
__alignment_not_unique  235423

tail -10 SRR25630297_htseq_stranded_yes_counts.txt
snap_masked-ptg003120l-processed-gene-0.9       0
snap_masked-ptg003128l-processed-gene-0.5       0
snap_masked-ptg003215l-processed-gene-0.5       0
snap_masked-ptg003293l-processed-gene-0.3       0
snap_masked-ptg003325l-processed-gene-0.4       0
__no_feature    33667781
__ambiguous     1603
__too_low_aQual 0
__not_aligned   3854905
__alignment_not_unique  1167098

tail -10 SRR25630297_htseq_stranded_reverse_counts.txt
snap_masked-ptg003120l-processed-gene-0.9       7
snap_masked-ptg003128l-processed-gene-0.5       840
snap_masked-ptg003215l-processed-gene-0.5       56
snap_masked-ptg003293l-processed-gene-0.3       0
snap_masked-ptg003325l-processed-gene-0.4       95
__no_feature    12501778
__ambiguous     146152
__too_low_aQual 0
__not_aligned   3854905
__alignment_not_unique  1167098
  
  
  

[Describe whether your reads are "string-specific", why you think they are, any evidence, and which stranded parameter is appropriate and why]
The reads are strand specific since the reverse strand specified reads are showing there are fewer unassigned reads compared to the forward strand parameter. The no_feature for the "yes" paramerter has 33667781 reads that are not assigned to any gene feautre but the "reverse" parameter shows 12501778 reads that are unassigned. This means the "reverse" is assigning reads to genes successfully.




16. BONUS - Turn your commands from part 1 and 2 into a script with a loop going through your two SRA files

## Bonus (optional!)

Review the [publication](https://doi.org/10.1093/molbev/msae021) from PRJNA1005244 or the third chapter of the [thesis](https://canvas.uoregon.edu/courses/266187/files/22059308?module_item_id=5380118) for the PRJNA1005245 dataset. See if this information leads to any additional insight of your analysis.

[Add insights to the dataset]

## Upload your:
- [ ] lab notebook
- [ ] Talapas batch script/code
- [ ] FastQC plots
- [ ] counts files generated from htseq-count (in a folder would be nice; **only include the counts files that would be used in a future differential RNA-seq analysis: use the format Species_sample_tissue_age/size_sample#_readnumber_htseqcounts_[revORyes]stranded.txt**)
- [ ] pdf report (see below; turn into both Github AND Canvas)
- [ ] and any additional plots, code, or code output

to GitHub.
    
### Pdf report details
You should create a pdf file (using Rmarkdown) with a high-level report including:

- [ ] all requested plots
- [ ] answers to questions
- [ ] mapped/unmapped read counts from PS8 script (in a nicely formatted table)
- [ ] It should be named `QAA_report.pdf`
- [ ] Include at the top level of your repo
- [ ] ALSO, submit it to Canvas.

> [!TIP]
> You may need to install LaTeX to knit your rmarkdown into a pdf file. Run `tinytex::install_tinytex()` to install it on R.
   
The three parts of the assignment should be clearly labeled. Be sure to title and write a descriptive figure caption for each image/graph/table you present. 

> [!TIP]
> Think about figure captions you've read and discussed in Journal Club. Find some good examples to model your captions on.




















































