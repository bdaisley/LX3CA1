# CALI-LX3-field-study

This code is related to a manuscript currently under review. Public access will be granted upon manuscript acceptance for publication. 

# Figures
Code used to generated each figure is provided in a separate file for clarity. If running in a singular script format, ensure scripts are run in order (e.g. Fig1, Fig2, Fig3, etc) to prevent errors.

# Custom scripts for assigning denovo taxonomy to sequences derived from unclassified microbial dark matter
"custom_tax_script.sh" and "custom_tax_script_Rscript.R" should both be saved in the same directory, along with two necessary files: 1) ASV count table generated from DADA2 pipeline, and 2) Fasta file of ASV sequences (names should match rownames of ASV count table). See example files for formatting reference.

Ensure permissions are granted to both files (e.g. chmod +x custom_tax_script.sh | custom_tax_script_Rscript.R) and that you set the database locations in the "custom_tax_script.sh" file before starting.

Dependencies:

usearch (v11.0.667_i86linux32 used for testing)
vsearch (vsearch v2.7.0_linux_x86_64 used for testing)

To run the full pipeline assuming a standard 16S rRNA sequencing dataset, simply open a new terminal in the directory containing the 2 scripts and 2 input files described above and run:

Basic command line format:

./custom_tax_script.sh file1 file2

---------- Example ----------
  
Lets consider a case where your files are named as follows:
  
1) ASV_seqs.fasta
2) ASV_counts.txt

./custom_tax_script.sh ASV_seqs.fasta ASV_counts.txt

