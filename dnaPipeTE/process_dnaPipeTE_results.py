import sqlite3
import os
import sys
import tempfile
import glob
import re

#python3 process_dnaPipeTE_results.py test



def main():  
  base_dir = '/work/users/d/t/dturissi/histoplasma/dnaPipeTE'
  output_base = 'output'
  longleaf_output_dir = os.path.join(base_dir, output_base)
   
  db_file = '/work/users/d/t/dturissi/histoplasma/te_analysis/histoplasma_tes.db'
  
  conn = sqlite3.connect(db_file)    
  conn.execute(f"""drop table if exists histoplasma_dnaPipeTE""")
  conn.execute(f"""create table histoplasma_dnaPipeTE
                   (hd_id int primary key,
                    sample_id varchar(100), 
                    genome_coverage decimnal(5,2),
                    sample_number int,
                    RM_t varchar(3,2),
                    RM_anno varchar(50),
                    dnaPipeTE_contig varchar(100),
                    RM_subclass varchar(50),
                    RM_class varchar(50),
                    contig_len int,
                    lib_start int,
                    lib_end int,
                    lib_hit_prop float,
                    genome_per float)""")

  conn.execute(f"create index idx_hd_sample_id on histoplasma_dnaPipeTE(sample_id)")  

  
  sample_query = conn.execute("""select sample_id 
                                 from histoplasma_manuscript_samples m, histoplasma_tree_names t
                                 where m.tree_name = t.tree_name""")
  
  manuscript_samples = []
  for sample_id, in sample_query:
    manuscript_samples.append(sample_id)

  conn.close()
       
  hd_id = 0                             
  with tempfile.NamedTemporaryFile(mode='w') as t: 
    for results_dir in glob.glob('output/*'):
      #get values from name of result file
      sample_id, genome_coverage, sample_number, RM_t = re.match(r'output/(.+)_([\.0-9]+)_(\d+)_([\.0-9]+)', results_dir).groups()    
      
      if sample_id in manuscript_samples: 
        reads_per_comp_file = os.path.join(results_dir, 'reads_per_component_and_annotation')  
        counts_file = os.path.join(results_dir, 'Counts.txt')  
        sample_results_file = os.path.join(results_dir, 'Annotation', 'one_RM_hit_per_Trinity_contigs')  
        
        if os.path.exists(sample_results_file):
          count_total = 0
          with open(counts_file, 'r') as c:
            for line in c:
              values = line.strip().split()
              if len(values) > 0:
                label, count = values
                if label == 'Total':
                  count_total = int(count)
          
          aligned_nts = {}
          with open(reads_per_comp_file, 'r') as a:
            for line in a:
              values = line.strip().split()  
              if len(values) == 7:  
                num_reads, aligned_nt, contig, RM_len, RM_anno, RM_subclass, hit_prop = values 
                aligned_nts[contig] = int(aligned_nt)
          
          
          with open(sample_results_file, 'r') as r:  
            for line in r:
              values = line.strip().split()  
              contig, RM_len, hit_prop, RM_anno, RM_subclass, lib_len, lib_coord, lib_hit_prop = values
              
              if contig in aligned_nts:
                #convert TIR/Tc1_Mariner to DNA/DTT since they are synonymous but EDTA reports them separately
                #it was previous confirmed that there is no overlap between hits so this doesn't result in double counting
                #TIR/PiggyBac and LTR/Gypsy are simply renamed
                if RM_subclass == 'TIR/Tc1_Mariner':
                  RM_subclass = 'DNA/DTT'
                elif RM_subclass == 'TIR/PiggyBac':
                  RM_subclass = 'DNA/DTB'
                elif RM_subclass == 'LTR/Gypsy':
                  RM_subclass = 'LTR/Ty3'
                              
                RM_class = RM_subclass.split('/')[0]
                
                lib_start, lib_end = lib_coord.replace('[', '').replace(']', '').split('-')
                
                if RM_class in ['LINE', 'LTR', 'SINE', 'PLE', 'Retroposon']:
                  RM_class = 'RNA TE'
                
                if RM_class in ['DNA', 'RC', 'TIR']:
                  RM_class = 'DNA TE'
                
                genome_per = 100 * aligned_nts[contig] / count_total
                hd_id += 1
                t.write(f"""{hd_id}\t{sample_id}\t{genome_coverage}\t{sample_number}\t{RM_t}\t{RM_anno}\t{contig}\t{RM_subclass}\t{RM_class}\t{aligned_nt}\t{lib_start}\t{lib_end}\t{lib_hit_prop}\t{genome_per}\n""")

    t.flush() 
    os.system(f"""sqlite3 {db_file} ".mode tabs" ".import {t.name} histoplasma_dnaPipeTE" """)
                       



       
if __name__ == '__main__':
  main()
