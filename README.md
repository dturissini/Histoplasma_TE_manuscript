## Histoplasma TE analysis
This git repo contain scripts used for the analyses of Transposable Elements (TEs) for the manuscript "The dynamics of transposable element content in the genome of the human pathogen Histoplasma"

The scripts in this git repo make extensive use of a sqlite database that was created for this project. Note that multiple naming schemes had been used for isolates and file names prior to the creation of scripts in ths git repo, as a result multiple ids may be used for the same same isolate: sample_id, tree_name, sample_short_name. The database table histoplasma tree names was created to convert between the different naming schemes.

The specific commands used to run dnaPipeTE and process results can be found in dnaPipeTE/dnaPipeTE_steps.txt The specific commands used for all subsequent steps can be found in te_analysis_manuscript_steps.txt.

### Running dnaPipeTE

1. load_histoplasma_fastq.py 
    * Creates a database table linking fastq files to sample ids.

2.  make_run_dnaPipeTE_scripts.py
    * Creates scripts for running dnaPipeTE with a RepeatModeler (RM) library using slurm.

3.  make_run_dnaPipeTE_EDTA_lib_scripts.py
    * Creates scripts for running dnaPipeTE with a isolate-specific EDTA libraries using slurm.

4.  process_dnaPipeTE_results.py
    * Processes dnaPipeTE RM lib results and stores them in a database table.


5.  process_dnaPipeTE_results_edta_lib.py
    * Processes dnaPipeTE EDTA lib results and stores them in a database table.


### Process metadata and assembly stats

1. python3 load_histoplasma_manuscript_samples.py
    * Creates a database table with just the isolates used in the manuscript.
    
2.  load_histoplasma_assembly_stats.py
    * Creates a database table with information about each isolate's genome assembly.
    
3.  python3 make_histoplasma_tree_names.py
    * Creates a database table containing the different names used for each isolate across multiple preexisting files. It also contains the species for each sample.
    
4.  python3 get_histoplasma_scaffold_lengths.py
    * Creates database tables containing the length of each scaffold and the whole genome for each isolate.

### Process EDTA and EDTA-RM results

1. python3 process_results_edta_rm.py
    * Processes EDTA-RM results and stores them in a database table.

2. python3 process_results_edta_only.py
    * Processes EDTA-only results and stores them in a database table.


### Steps run to confirm that there is no overlap between DNA/DTT and TIR/TC1_Mariner
EDTA separately reported DNA/DTT and TIR/TC1_Mariner TEs despite them being considered synonymous. Before renaming TIR/TC1_Mariner to DNA/DTT, we confirmed that the TEs were fully independent and never overlapped to prevent double-counting.
After no overlap was confirmed, process_results_edta_rm.py and process_results_edta_only.py were modified to convert TIR/TC1_Mariner to DNA/DTT and were rerun.


1. python3 get_dna_dtt_tir_tc1_mariner_overlaps.py
    * Identify any possibly cases where DNA/DTT and TIR/TC1_Mariner TEs overlap in the genomes of any sample

2. dna_dtt_tir_tc1_mariner_queries.sql
    * queries to confirm no overlap 



###Process results from all four TE approaches
1. python3 make_histoplasma_te_genome_pers.py
    * Creates database tables containing the genomic percentages for each TE for each isolate.

2. Rscript get_te_per_genome_len_corrs.R
    * Th\is R script does multiple regressions for the genomic percentage for TE classes versus genome length and stores the results in a database table.


###Generate figures and tables for manuscript
1. Rscript histoplasma_te_genomic_per_comps.R
    * Creates plots in pdf files comparing the genomic TE percentages between the different TE approaches.
    
2. Rscript histoplasma_te_genomic_per_heatmaps.R
    * Creates heatmaps of genomic TE percentages in pdf files for the different TE approaches.

3. Rscript histoplasma_te_genomic_per_length_corr.R
    * Creates plots in pdf files of the genomic TE percentages versus genome lengths.

4. Rscript histoplasma_TE_hit_props.R
    * Creates plots in pdf files showing the proportions of the TEs from the libraries that were found in the isolate genomes for the different TE approaches.

5. histoplasma_te_manuscript_tables.sql
    * Contains queries used to make tables and supplemental tables for the manuscript. Each query creates a .csv file.

