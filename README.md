Rejoin-seq used custom scripts and terminal commands to pull the necessary information to to quantify the fraction of template switched using template_switch_quantify.py.  CRISPResso2 was used for indel analysis of rejoined outcomes. CIGAR strings were extracted using cigar_mapping_reading.py to resolve deletion and insertion sizes along the amplicon used.  Finally, repair outcomes were collated by barcode (gene targeted) using ddr_data_cluster.py to generate heatmap clustering of common repair patterns.  

Details below highlight script usage, dependencies, and intermediate terminal commands for arranging data.  Example files are provided along with expected output files to confirm successful execution of scripts.

Requirements:  
seqtk - 1.3-r106 or 1.3-r116-dirty  
bedtools - 2.27 or 2.26  
python -  3.11.5 or 3.10  
pandas - 2.1.3 or 2.0.3  
numpy - 1.20.3, 1.21.5, 1.22.3, 1.24.4 or 1.26.2  
seaborn - 0.13.2
biopython
openpyxl 3.1.5

Usage:  
template_switch_quantify.py
python template_switch_quantify1.py -r1 U6DDR_Cycing_CLV_1_R1.csv -r2 U6DDR_Cycing_CLV_1_R2.csv -ref human_DDR_minipool_ref.csv
(input R1 and R2 .csv files --> ID, barcode, seq quality + ref file)
as 'w_score.csv' file
requires: seaborn-0.13.2, biopython

CRISPResso2 --> https://github.com/pinellolab/CRISPResso2 
CRISPResso --fastq_r1 U6DDR_Cycing_CLV_1_S16_L001_R1_001_10k.fastq.gz --fastq_r2 U6DDR_Cycing_CLV_1_S16_L001_R2_001_10k.fastq.gz --amplicon_seq GCAAACCTGGACAAGATGCTGGCATCGCCATCCAGCAGAGCGGTACCgatccgacgcgccatctctaggcccgcgccggccccctcgcacggacttgtgggagaagctcggctactcccctgccccggttaatttgcatataatatttcctagtaactatagaggcttaatgtgcgataaaagacagataatctgttctttttaatactagctacattttacatgataggcttggatttctataacttcgtatagcatacattatacgaagttataaacagcacaaaaggaaactcaccctaactgtaaagtaattgtgtgttttgagactataaGtatcccttggagaaCCAcctTGTTGG --bam_output -n U6DDR_Cycing_CLV_1sm

cigar_mapping_reading.py:
1. CRISPResso2 to get BAM file
2. extract ID and CIGAR using: samtools view your_file.bam | cut -f 1,6 > output.txt
3. add an ID and Cigar header
4. cigar_mapping_reading.py
python cigar_mapping_reading2.py -i U6DDR/

ddr_data_cluster.py:
python ddr_data_cluster1.py -i_excel analyzed_U6DDR_Cycing_CLV_1_w_annotation.xlsx
requires: openpyxl 3.1.5