import sqlite3
import os
import sys
import tempfile

#python3 load_histoplasma_genome_sizes.py



def main():  
  db_file = '/work/users/d/t/dturissi/histoplasma/te_analysis/histoplasma_tes.db'
  fasta_dir = '/proj/matutelb/projects/fungus/histoplasma/Assembly/genomes'
  
  conn = sqlite3.connect(db_file)    
  conn.execute(f"""drop table if exists histoplasma_scaffold_lens""")
  conn.execute(f"""create table histoplasma_scaffold_lens
                   (tree_name varchar(100), 
                    sample_id varchar(100), 
                    scaffold varchar(100),
                    scaffold_len int)""")

  conn.execute(f"create index idx_hsl_sample_id on histoplasma_scaffold_lens(sample_id)")  
  conn.execute(f"create index idx_hsl_tree_name on histoplasma_scaffold_lens(tree_name)")  

  conn.execute(f"""drop table if exists histoplasma_genome_lens""")
  conn.execute(f"""create table histoplasma_genome_lens
                   (sample_id varchar(100) primary key, 
                    genome_len int)""")


  nonstandard_fasta_files = {'NACVFR_Histo_HC1070058_2_S17_L002': 'NACVFR_Histo_HC1070058.sorted.fasta',
                             '104_P_19_S5': '104_P_19.sorted.fasta',
                             'HISSP-FGBON2001-xx-CL-GUF-xxxx-036-BB_S21': 'HISSP-FGBON2001.sorted.fasta',
                             'HISSP-FGAMA2041-uu-CL-GUF-xxxx-036-BB_S12_L001': 'HISSP-FGAMA2041.sorted.fasta',
                             'HISSP-FGPOE2043-uu-CL-GUF-xxxx-036-BB_S4_L001': 'HISSP-FGPOE2043.sorted.fasta',
                             'NACVFR_Histo_HC4137_S92_L008': 'NACVFR_Histo_HC4137.sorted.fasta',
                             'NACVFR_Histo_HC970591_S81_L007': 'NACVFR_Histo_HC970591.sorted.fasta',
                             'HISSP-SDOR2042-xx-CL-SUR-xxxx-036-BB_S14': 'HISSP-SDOR2042.sorted.fasta'}
  
  tree_name_query = conn.execute(f"""select sample_id, t.tree_name
                                     from histoplasma_manuscript_samples m, histoplasma_tree_names t
                                     where m.tree_name = t.tree_name""")
     
                              
  with tempfile.NamedTemporaryFile(mode='w') as t: 
    for sample_id, tree_name in tree_name_query:
      trimmed_sample_id = sample_id.replace('Histoplasma_capsulatum_', '').replace('Histoplasma_mississippiense_I_', '').replace('Histoplasma_mississippiense_II_', '').replace('Histoplasma_ohiense_', '')
      
      tree_name_fastq_file = os.path.join(fasta_dir, tree_name + '.sorted.fasta')
      sample_id_fastq_file = os.path.join(fasta_dir, sample_id + '.sorted.fasta')
      sample_id_trimmed_fastq_file = os.path.join(fasta_dir, trimmed_sample_id + '.sorted.fasta')
      
      fastq_file = ''
      if tree_name in nonstandard_fasta_files:
        fastq_file = os.path.join(fasta_dir, nonstandard_fasta_files[tree_name])
      elif os.path.isfile(tree_name_fastq_file):
        fastq_file = tree_name_fastq_file
      elif os.path.isfile(sample_id_fastq_file):
        fastq_file = sample_id_fastq_file
      elif os.path.isfile(sample_id_trimmed_fastq_file):
        fastq_file = sample_id_trimmed_fastq_file
      else:
        print(f"{tree_name} {sample_id} fasta does not exist")
      
      if fastq_file != '':
        with open(fastq_file, 'r') as f:  
          scaffold_len = 0
          for line in f:
            line = line.strip()
            if line[0] == '>':
              if scaffold_len > 0:
                t.write(f"""{tree_name}\t{sample_id}\t{scaffold}\t{scaffold_len}\n""")
              
              scaffold = line[1:]
              scaffold_len = 0
            else:
              scaffold_len += len(line)
            
          t.write(f"""{tree_name}\t{sample_id}\t{scaffold}\t{scaffold_len}\n""")
        
    conn.close()
    t.flush() 
    os.system(f"""sqlite3 {db_file} ".mode tabs" ".import {t.name} histoplasma_scaffold_lens" """)
                       

  conn = sqlite3.connect(db_file)    
  conn.execute(f"""insert into histoplasma_genome_lens
                   select sample_id, sum(scaffold_len)
                   from histoplasma_scaffold_lens
                   group by sample_id""")
  
  conn.execute(f"""insert into histoplasma_genome_lens
                   select 'Ep_130_s_7', 27837035""")

  conn.commit()             
  conn.close()             
       
if __name__ == '__main__':
  main()
