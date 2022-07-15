## California field study on patty- and spray-based LX3 suppelementation in honey bees

This supplementary code is related to a manuscript currently under review. Public access will be granted upon manuscript acceptance for publication. 

### Figures
Code used to generated each figure is provided in a separate file for clarity. If running in a singular script format, ensure scripts are run in order (e.g. Fig1, Fig2, Fig3, etc) to prevent errors.

### Custom scripts for assigning denovo taxonomy to sequences derived from unclassified microbial dark matter
"custom_tax_script.sh" and "custom_tax_script_Rscript.R" should both be saved in the same directory, along with two necessary files: 1) a fasta file of the relevant ASV sequences, and 2) ASV count table generated from DADA2 pipeline. Header names in fasta file should match rownames of ASV count table. See example files for formatting reference.

Ensure permissions are granted to both files (e.g. chmod +x custom_tax_script.sh | chmod +x custom_tax_script_Rscript.R, on a linux-based system) and that you set the database locations in the "custom_tax_script.sh" file before starting.

Dependencies:
usearch (v11.0.667_i86linux32 used for testing) # Available from https://www.drive5.com/usearch/download.html

To run the full pipeline assuming a standard 16S rRNA sequencing dataset, simply open a new terminal in the directory containing the 2 scripts and 2 input files described above and run:


<pre class="r"><code>#Basic command line format:
./custom_tax_script.sh [input_file_1] [input_file_2]

# input_file_1 being the fasta file, and input_file2 being the ASV counts table</code></pre>

---------- Example ----------</br>
Lets consider a case where your files are named as follows:
  
1) ASV_seqs.fasta
2) ASV_counts.txt

<pre class="r"><code>./custom_tax_script.sh ASV_seqs.fasta ASV_counts.txt</code></pre>


