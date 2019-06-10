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
