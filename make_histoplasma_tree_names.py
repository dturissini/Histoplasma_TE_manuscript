import sqlite3
import os
import sys
import tempfile
import re


#python3 make_histoplasma_tree_names.py


def main():  
  db_file = '/work/users/d/t/dturissi/histoplasma/te_analysis/histoplasma_tes.db'
  
  conn = sqlite3.connect(db_file)    
  conn.execute(f"""drop table if exists histoplasma_tree_names""")
  conn.execute(f"""create table histoplasma_tree_names
                   (sample_id varchar(100) primary key, 
                    tree_name varchar(100),
                    sample_short_name varchar(100),
                    species varchar(50))""")

  conn.execute(f"create index idx_htn_tree_name on histoplasma_tree_names(tree_name)")  
  conn.execute(f"create index idx_htn_sample_short_name on histoplasma_tree_names(sample_short_name)")  

  conn.close()
  

  species_convert = {'nam2': 'H. ohiense', 
                     'nam1': 'H. mississippiense', 
                     'suramericamun': 'H. suramericanum', 
                     'Amazon_III': 'Amazon-III', 
                     'Africa': 'Africa',
                     'capsu': 'H. capsulatum', 
                     'mz5-like': 'mz5-like', 
                     'india': 'India',
                     'odd_clinical': 'B05821',
                     'outgroup': 'B. parva'}
  
  species_names = {}
  sample_short_names = {}
  with open('shorten_names_isolates.csv') as s:
    next(s)
    for line in s:
      species, tree_name, short_name = line.strip().split(',')
      species_names[tree_name] = species_convert[species]
      sample_short_names[tree_name] = short_name
            
                              
  with tempfile.NamedTemporaryFile(mode='w') as t_tree: 
    with open('histoplasma_te_fastqs.txt', 'r') as r:  
      for line in r:
        sample_id, fastq_file = line.strip().split("\t") 
         
        tree_name = fastq_file.split('/')[-1].replace('_R1_001.fastq.gz', '').replace('_1.fastq.gz', '')
        
        if tree_name[:4] == 'ES2_':
          tree_name = tree_name[:7]
                      
        if tree_name == 'HISSP-FGBON2001-xx-CL-GUF-xxxx-036-BB_S21_L001':
          tree_name = 'HISSP-FGBON2001-xx-CL-GUF-xxxx-036-BB_S21'
       
        if tree_name == 'SECH_101-Nam2_G184A_ATCACG_L006':
          tree_name = 'caps_G184A_ATCACG_L006'
        
        sample_short_name = ''
        if tree_name in sample_short_names:
          sample_short_name = sample_short_names[tree_name]
        
        species = ''
        if tree_name in species_names:
          species = species_names[tree_name]
          
        t_tree.write(f"""{sample_id}\t{tree_name}\t{sample_short_name}\t{species}\n""")
          
          
    t_tree.flush() 
    os.system(f"""sqlite3 {db_file} ".mode tabs" ".import {t_tree.name} histoplasma_tree_names" """)

                       
             

       
if __name__ == '__main__':
  main()
