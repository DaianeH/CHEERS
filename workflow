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

bedToBam -i E116-H3K27ac.tagAlign -g human.hg19.genome > E116-H3K27ac.tagAlign.bam


#5) Run multicov to count reads within each peak
#Generate bai index from bam
ml samtools/1.9
for file in *bam; do samtools index $file ; done
bedtools multicov -bams reads/*.bam -bed peaks/H3K27ac_ENCODE_peaks_merged.bed > H3K27ac_ENCODE_counts.txt

bedtools multicov -bams reads/E003-H3K27ac.tagAlign.gz.bam  reads/E034-H3K27ac.tagAlign.gz.bam  reads/E062-H3K27ac.tagAlign.gz.bam  reads/E090-H3K27ac.tagAlign.gz.bam  reads/E112-H3K27ac.tagAlign.gz.bam reads/E004-H3K27ac.tagAlign.gz.bam  reads/E037-H3K27ac.tagAlign.gz.bam  reads/E063-H3K27ac.tagAlign.gz.bam  reads/E091-H3K27ac.tagAlign.gz.bam  reads/E113-H3K27ac.tagAlign.gz.bam reads/E005-H3K27ac.tagAlign.gz.bam  reads/E038-H3K27ac.tagAlign.gz.bam  reads/E065-H3K27ac.tagAlign.gz.bam  reads/E092-H3K27ac.tagAlign.gz.bam  reads/E114-H3K27ac.tagAlign.gz.bam reads/E006-H3K27ac.tagAlign.gz.bam  reads/E039-H3K27ac.tagAlign.gz.bam  reads/E066-H3K27ac.tagAlign.gz.bam  reads/E093-H3K27ac.tagAlign.gz.bam  reads/E115-H3K27ac.tagAlign.gz.bam reads/E007-H3K27ac.tagAlign.gz.bam  reads/E040-H3K27ac.tagAlign.gz.bam  reads/E067-H3K27ac.tagAlign.gz.bam  reads/E094-H3K27ac.tagAlign.gz.bam  reads/E116-H3K27ac.tagAlign.bam reads/E008-H3K27ac.tagAlign.gz.bam  reads/E041-H3K27ac.tagAlign.gz.bam  reads/E068-H3K27ac.tagAlign.gz.bam  reads/E095-H3K27ac.tagAlign.gz.bam  reads/E117-H3K27ac.tagAlign.gz.bam reads/E011-H3K27ac.tagAlign.gz.bam  reads/E042-H3K27ac.tagAlign.gz.bam  reads/E069-H3K27ac.tagAlign.gz.bam  reads/E096-H3K27ac.tagAlign.gz.bam  reads/E118-H3K27ac.tagAlign.gz.bam reads/E012-H3K27ac.tagAlign.gz.bam  reads/E043-H3K27ac.tagAlign.gz.bam  reads/E071-H3K27ac.tagAlign.gz.bam  reads/E097-H3K27ac.tagAlign.gz.bam  reads/E119-H3K27ac.tagAlign.gz.bam reads/E013-H3K27ac.tagAlign.gz.bam  reads/E044-H3K27ac.tagAlign.gz.bam  reads/E072-H3K27ac.tagAlign.gz.bam  reads/E098-H3K27ac.tagAlign.gz.bam  reads/E120-H3K27ac.tagAlign.gz.bam reads/E014-H3K27ac.tagAlign.gz.bam  reads/E045-H3K27ac.tagAlign.gz.bam  reads/E073-H3K27ac.tagAlign.gz.bam  reads/E099-H3K27ac.tagAlign.gz.bam  reads/E121-H3K27ac.tagAlign.gz.bam reads/E015-H3K27ac.tagAlign.gz.bam  reads/E046-H3K27ac.tagAlign.gz.bam  reads/E074-H3K27ac.tagAlign.gz.bam  reads/E100-H3K27ac.tagAlign.gz.bam  reads/E122-H3K27ac.tagAlign.gz.bam reads/E016-H3K27ac.tagAlign.gz.bam  reads/E047-H3K27ac.tagAlign.gz.bam  reads/E075-H3K27ac.tagAlign.gz.bam  reads/E101-H3K27ac.tagAlign.gz.bam  reads/E123-H3K27ac.tagAlign.gz.bam reads/E017-H3K27ac.tagAlign.gz.bam  reads/E048-H3K27ac.tagAlign.gz.bam  reads/E076-H3K27ac.tagAlign.gz.bam  reads/E102-H3K27ac.tagAlign.gz.bam  reads/E124-H3K27ac.tagAlign.gz.bam reads/E019-H3K27ac.tagAlign.gz.bam  reads/E049-H3K27ac.tagAlign.gz.bam  reads/E078-H3K27ac.tagAlign.gz.bam  reads/E103-H3K27ac.tagAlign.gz.bam  reads/E125-H3K27ac.tagAlign.gz.bam reads/E020-H3K27ac.tagAlign.gz.bam  reads/E050-H3K27ac.tagAlign.gz.bam  reads/E079-H3K27ac.tagAlign.gz.bam  reads/E104-H3K27ac.tagAlign.gz.bam  reads/E126-H3K27ac.tagAlign.gz.bam reads/E021-H3K27ac.tagAlign.gz.bam  reads/E055-H3K27ac.tagAlign.gz.bam  reads/E080-H3K27ac.tagAlign.gz.bam  reads/E105-H3K27ac.tagAlign.gz.bam  reads/E127-H3K27ac.tagAlign.gz.bam reads/E022-H3K27ac.tagAlign.gz.bam  reads/E056-H3K27ac.tagAlign.gz.bam  reads/E084-H3K27ac.tagAlign.gz.bam  reads/E106-H3K27ac.tagAlign.gz.bam  reads/E128-H3K27ac.tagAlign.gz.bam reads/E026-H3K27ac.tagAlign.gz.bam  reads/E058-H3K27ac.tagAlign.gz.bam  reads/E085-H3K27ac.tagAlign.gz.bam  reads/E108-H3K27ac.tagAlign.gz.bam  reads/E129-H3K27ac.tagAlign.gz.bam reads/E029-H3K27ac.tagAlign.gz.bam  reads/E059-H3K27ac.tagAlign.gz.bam  reads/E087-H3K27ac.tagAlign.gz.bam  reads/E109-H3K27ac.tagAlign.gz.bam reads/E032-H3K27ac.tagAlign.gz.bam  reads/E061-H3K27ac.tagAlign.gz.bam  reads/E089-H3K27ac.tagAlign.gz.bam  reads/E111-H3K27ac.tagAlign.gz.bam -bed peaks/H3K27ac_ENCODE_peaks_merged.bed > H3K27ac_ENCODE_counts.txt



#6) Remove bottom 10% peaks with lower count
#Include a column in the end of H3K27ac_ENCODE_counts.txt with the sum of counts across a row
awk '{sum=0; for(i=4; i<=101; i++) sum += $i; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,$41,$42,$43,$44,$45,$46,$47,$48,$49,$50,$51,$52,$53,$54,$55,$56,$57,$58,$59,$60,$61,$62,$63,$64,$65,$66,$67,$68,$69,$70,$71,$72,$73,$74,$75,$76,$77,$78,$79,$80,$81,$82,$83,$84,$85,$86,$87,$88,$89,$90,$91,$92,$93,$94,$95,$96,$97,$98,$99,$100,$101, sum}' H3K27ac_ENCODE_counts.txt | sort -k 102,102  > H3K27ac_ENCODE_counts_sums_complete.txt

# 121686 bottom 10% lower counts
tail -1095172 H3K27ac_ENCODE_counts_sums_complete.txt > H3K27ac_ENCODE_counts_sums_complete_nodecile.txt

# Remove sum column
awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,$41,$42,$43,$44,$45,$46,$47,$48,$49,$50,$51,$52,$53,$54,$55,$56,$57,$58,$59,$60,$61,$62,$63,$64,$65,$66,$67,$68,$69,$70,$71,$72,$73,$74,$75,$76,$77,$78,$79,$80,$81,$82,$83,$84,$85,$86,$87,$88,$89,$90,$91,$92,$93,$94,$95,$96,$97,$98,$99,$100,$101}' H3K27ac_ENCODE_counts_sums_complete_nodecile.txt > H3K27ac_ENCODE_counts_sums_nodecile.txt

# Sed rep space with tab
sed -i "s/ /\t/g" H3K27ac_ENCODE_counts_sums_nodecile.txt

# Insert header:
sed  -i '1i chr start end H1_Cell_Line  CD3_Primary_Cells_Peripheral_UW Peripheral_Blood_Mononuclear_Primary_Cells  Fetal_Muscle_Leg  Thymus  H1_BMP4_Derived_Mesendoderm_Cultured_Cells  CD4_Memory_Primary_Cells  Adipose_Nuclei  Fetal_Placenta  Spleen  H1_BMP4_Derived_Trophoblast_Cultured_Cells  CD4_Naive_Primary_Cells Aorta	Fetal_Stomach A549_EtOH_0.02pct_Lung_Carcinoma  H1_Derived_Mesenchymal_Stem_Cells 	Adult_Liver Fetal_Thymus  Dnd41_TCell_Leukemia  H1_Derived_Neuronal_Progenitor_Cultured_Cells CD4+_CD25-_CD45RO+_Memory_Primary_Cells Brain_Angular_Gyrus Gastric GM12878_Lymphoblastoid  H9_Cell_Line  CD4+_CD25-_IL17-_PMA-Ionomycin_stimulated_MACS_purified_Th_Primary_Cells  Brain_Anterior_Caudate  Left_Ventricle  HeLa-S3_Cervical_Carcinoma  hESC_Derived_CD184+_Endoderm_Cultured_Cells CD4+_CD25-_IL17+_PMA-Ionomcyin_stimulated_Th17_Primary_Cells  Brain_Cingulate_Gyrus Lung  HepG2_Hepatocellular_Carcinoma  hESC_Derived_CD56+_Ectoderm_Cultured_Cells  CD4+_CD25-_Th_Primary_Cells Brain_Hippocampus_Middle  Ovary HMEC_Mammary_Epithelial hESC_Derived_CD56+_Mesoderm_Cultured_Cells  CD4+_CD25+_CD127-_Treg_Primary_Cells  Brain_Inferior_Temporal_Lobe  Pancreas  HSMM_Skeletal_Muscle_Myoblasts  HUES48_Cell_Line  CD4+_CD25int_CD127+_Tmem_Primary_Cells  Brain_Mid_Frontal_Lobe  Placenta_Amnion HSMMtube_Skeletal_Muscle_Myotubes_Derived_from_HSMM HUES6_Cell_Line CD56_Primary_Cells  Brain_Substantia_Nigra  Psoas_Muscle  HUVEC_Umbilical_Vein_Endothelial_Cells  HUES64_Cell_Line  CD8_Naive_Primary_Cells Colonic_Mucosa  Rectal_Mucosa.Donor_29  K562_Leukemia IMR90_Cell_Line CD8_Memory_Primary_Cells  Colon_Smooth_Muscle	Rectal_Mucosa.Donor_31  Monocytes-CD14+_RO01746 iPS-18_Cell_Line  Chondrocytes_from_Bone_Marrow_Derived_Mesenchymal_Stem_Cell_Cultured_Cells  Duodenum_Smooth_Muscle  Rectal_Smooth_Muscle  NH-A_Astrocytes iPS-20b_Cell_Line	Mobilized_CD34_Primary_Cells_Female Esophagus Right_Atrium  NHDF-Ad_Adult_Dermal_Fibroblasts  iPS_DF_6.9_Cell_Line  Penis_Foreskin_Fibroblast_Primary_Cells_skin01  Fetal_Adrenal_Gland Right_Ventricle NHEK-Epidermal_Keratinocytes  iPS_DF_19.11_Cell_Line  Penis_Foreskin_Fibroblast_Primary_Cells_skin02  Fetal_Intestine_Large Sigmoid_Colon NHLF_Lung_Fibroblasts Bone_Marrow_Derived_Mesenchymal_Stem_Cell_Cultured_Cells  Penis_Foreskin_Keratinocyte_Primary_Cells_skin03  Fetal_Intestine_Small Skeletal_Muscle_Female  Osteoblasts CD14_Primary_Cells  Penis_Foreskin_Melanocyte_Primary_Cells_skin01  Pancreatic_Islets Small_Intestine CD19_Primary_Cells_Peripheral_UW  Penis_Foreskin_Melanocyte_Primary_Cells_skin03  Fetal_Muscle_Trunk  Stomach_Smooth_Muscle' H3K27ac_ENCODE_counts_sums_nodecile_TEST.txt




Link: 


create_LD_blocks.py

python /hpc/users/hemerd01/daiane/projects/BP/CHEERS/create_LD_blocks.py /hpc/users/hemerd01/daiane/projects/BP/CHEERS/gwas_for_input/BP_SNPs_CRCh37.txt /hpc/users/hemerd01/daiane/projects/BP/CHEERS/BP /hpc/users/hemerd01/daiane/projects/BP/CHEERS/ld


Output:
Parent directory is the name of the trait; within that directory there are subdirectories for all the chromosomes (chr1, chr2...) and within each of these there are .txt files named after the lead SNP (for example results_ld_rs9989735.txt). Each file contains all the SNPs in the LD with the lead (python indexing: snp name is at position 3, chr at position 0 and base pair information is at position 4)

****
python /hpc/users/hemerd01/hemerd01-old-home/daiane/projects/BP/CHEERS/create_LD_blocks.py /hpc/users/hemerd01/hemerd01-old-home/daiane/projects/BP/gwas_for_input/WHR_CHEERS.txt /hpc/users/hemerd01/hemerd01-old-home/daiane/projects/BP/CHEERS/CHEERS_output /hpc/users/hemerd01/hemerd01-old-home/daiane/projects/BP/CHEERS/CHEERS_output/ld
