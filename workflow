### PREPARE FILES FOR CHEERS ###

# 1) Download .tagalign files from consolidated epigenomes on H3K27ac from 
https://egg2.wustl.edu/roadmap/data/byFileType/alignments/consolidated/

#2) Download .narrowPeak files from consolidated epigenomes on H3K27ac from 
https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/

#3) Merge all peaks across all cell-types:
ml bedtools/2.27.1
cat *Peak > H3K27ac_ENCODE_peaks.bed
sort -k1,1 -k2,2n H3K27ac_ENCODE_peaks.bed > tmp && mv tmp H3K27ac_ENCODE_peaks.bed
bedtools merge -i H3K27ac_ENCODE_peaks.bed > H3K27ac_ENCODE_peaks_merged.bed
wc -l H3K27ac_ENCODE_peaks_merged.bed
1216858

#4) Convert tagAlign files to BAM
for file in *tagAlign.gz; do bedToBam -i $file -g human.hg19.genome > $file\.bam ; done


#5) Run multicov to count reads within each peak
bedtools multicov -bams reads/*.bam -bed peaks/H3K27ac_ENCODE_peaks_merged.bed > H3K27ac_ENCODE_counts.txt

#6) Remove bottom 10% peaks with lower count
#Include a column in the end of H3K27ac_ENCODE_counts.txt with the sum of counts across a row
awk

# Sort descending based on that sum value of last column

# Calculate how much is 10% of wc -l

# Select wc -l - 10% value with head




Link: 


create_LD_blocks.py

Description:
For a list of SNPs, it finds all the SNPs in LD (r2>0.8).

Usage:
Create a directory named after the trait. Then move into the directory and run the command.

mkdir trait_name
cd trait_name
python create_LD_blocks.py SNP_LIST OUTPUT_DIR LD_DIR

Output:
Parent directory is the name of the trait; within that directory there are subdirectories for all the chromosomes (chr1, chr2...) and within each of these there are .txt files named after the lead SNP (for example results_ld_rs9989735.txt). Each file contains all the SNPs in the LD with the lead (python indexing: snp name is at position 3, chr at position 0 and base pair information is at position 4)

****
python /hpc/users/hemerd01/hemerd01-old-home/daiane/projects/BP/CHEERS/create_LD_blocks.py /hpc/users/hemerd01/hemerd01-old-home/daiane/projects/BP/gwas_for_input/WHR_CHEERS.txt /hpc/users/hemerd01/hemerd01-old-home/daiane/projects/BP/CHEERS/CHEERS_output /hpc/users/hemerd01/hemerd01-old-home/daiane/projects/BP/CHEERS/CHEERS_output/ld
