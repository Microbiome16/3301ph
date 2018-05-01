#!/bin/bash --login
#PBS -l walltime=01:00:00
#PBS -l select=2:ncpus=4
#PBS -N join_reads
#PBS -A d411-training

cd $PBS_O_WORKDIR

#load modules and qiime
module load miniconda/python2

#loading virtualenv
echo "loading virtual env"

source activate qiime1

#setting temporary directory
mkdir -p ~/qiime_tmp
export TMPDIR=~/qiime_tmp

#join_paired_ends.py -f Read1.fastq.gz -r Read2.fastq.gz -b Index.fastq.gz -m SeqPrep -o joined1

#echo "splitting"
#time split_libraries_fastq.py -m ~/2018_02_smb/map_incomplete.tsv -i ~/2018_02_smb/joined1/seqprep_assembled.fastq.gz -b ~/2018_02_smb/joined1/seqprep_assembled.fastq_barcodes.fastq  -o dem3 -q 19 --store_qual_scores --rev_comp_barcode --rev_comp_mapping_barcodes

#count_seqs.py -i ~/2018_02_smb/dem2/seqs.fna -o countseqs

#pick_closed_reference_otus.py -i ~/2018_02_smb/dem2/seqs.fna -o ~/2018_02_smb/otus -a -o 8

#biom summarize-table -i ./merged_otu_table.biom

#core_diversity_analysis.py --recover_from_faliure -o cdout/ -i merged_otu_table.biom -m map.tsv -t 97_otus.tree -e 371590 --recover_from_faliure

#make_2d_plots.py -i unweighted_unifrac_pc.txt -m map.tsv_corrected.txt

#make_2d_plots.py -i weighted_unifrac_pc.txt -m map.tsv_corrected.txt

#compare_categories.py --method anosim -i ~/2018_02_smb/baron/cdout/bdiv_even30000/unweighted_unifrac_dm.txt -m ~/2018_02_smb/map.tsv -c SamplePh -o anosim_out -n 999

source deactivate
