#!/bin/bash

export db_BEEXACT="16S_dbs/BEExact/seeds_typestrains_and_bxid___mar2022.fasta"
export db_NBCI="16S_dbs/16S_ribosomal_RNA___NCBI/BLASTdb_vsearch_R-adjusted.fasta"

seqs=$1
counts=$2

cat $counts > ASV_bac_counts.txt

filename1=$(basename "$1")
fname1="${filename1%.*}"

filename2=$(basename "$2")
fname2="${filename2%.*}"

mkdir bac
mkdir bac/clustering

#Set permissions
#chmod +x custom_tax_script_Rscript.R


# Denovo custom tax script

##########################################################
########################################################## BACTERIA
##########################################################

#::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 1 - Assignment of taxonomy
#::::::::::::::::::::::::::::::::::::::::::::::::::
#BEExact seedsonly
usearch --usearch_global $seqs --db $db_BEEXACT \
            --strand plus --id 0.51 --maxaccepts 0 --maxrejects 0 --maxhits 1 --threads 18 --uc bac/UCout_BEExact_seedsonly_usearch.uc

    #NCBI R-adjusted
    usearch --usearch_global $seqs --db $db_NBCI \
                  --strand plus --id 0.51 --maxaccepts 0 --maxrejects 0 --maxhits 1 --threads 18 --uc bac/UCout_NCBI_usearch_R-adjusted.uc

#::::::::::::::::::::::::::::::::
# Step 2 - Denovo clustering
#::::::::::::::::::::::::::::::::
#eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
#eval "$(conda shell.bash hook)"
#conda activate qiime2-2021.4
#vsearch --sortbylength $seqs --output ${fname1}_sorted.fasta
#conda deactivate

#Using usearch for consistency - updated july/2022
usearch -sortbylength $seqs -fastaout ${fname1}_sorted.fasta -minseqlength 1

usearch -cluster_smallmem ${fname1}_sorted.fasta -id 0.990 -maxhits 1 -maxaccepts 0 -maxrejects 0 -uc bac/clustering/ASV_SPECIES.uc
usearch -cluster_smallmem ${fname1}_sorted.fasta -id 0.958 -maxhits 1 -maxaccepts 0 -maxrejects 0 -uc bac/clustering/SPECIES_GENUS.uc
usearch -cluster_smallmem ${fname1}_sorted.fasta -id 0.896 -maxhits 1 -maxaccepts 0 -maxrejects 0 -uc bac/clustering/GENUS_FAMILY.uc
usearch -cluster_smallmem ${fname1}_sorted.fasta-id 0.862 -maxhits 1 -maxaccepts 0 -maxrejects 0 -uc bac/clustering/FAMILY_ORDER.uc
usearch -cluster_smallmem ${fname1}_sorted.fasta -id 0.835 -maxhits 1 -maxaccepts 0 -maxrejects 0 -uc bac/clustering/ORDER_CLASS.uc
usearch -cluster_smallmem ${fname1}_sorted.fasta -id 0.808 -maxhits 1 -maxaccepts 0 -maxrejects 0 -uc bac/clustering/CLASS_PHYLUM.uc

#::::::::::::::::::::::::::::::::
# Step 3 - Run R script
#::::::::::::::::::::::::::::::::

Rscript custom_tax_script_Rscript.R

cp ASV_bac_counts_annotated.txt ${fname2}_annotated.txt
rm ASV_bac_counts.txt
rm ASV_bac_counts_annotated.txt
rm ${fname1}_sorted.fasta
