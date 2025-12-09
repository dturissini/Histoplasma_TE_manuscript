import sqlite3
import os
import sys
import tempfile

#python3 load_histoplasma_fastq.py



def main():  
  db_file = '/work/users/d/t/dturissi/histoplasma/te_analysis/histoplasma_tes.db'
  
  conn = sqlite3.connect(db_file)    
  conn.execute(f"""drop table if exists histoplasma_manuscript_samples""")
  conn.execute(f"""create table histoplasma_manuscript_samples
                   (tree_name varchar(100) primary key)""")

  conn.close()     
  
  #exclude samples with problematic genome assemblies
  exclude_samples = ['HISSP-FGPOE2043-uu-CL-GUF-xxxx-036-BB_S4_L001',
                     'HISSP-SDOR2042-xx-CL-SUR-xxxx-036-BB_S14',
                     'NACVFR_Histo_HC4137_S92_L008',
                     'NACVFR_Histo_HC970591_S81_L007']
                                
  with tempfile.NamedTemporaryFile(mode='w') as t:
    with open('trimmed_histo_dat.csv', 'r') as r:  
      next(r)
      for line in r:
        values = line.strip().split(",") 
        
        if values[1] not in exclude_samples:
          t.write(f"""{values[1]}\n""")
                    
    t.flush() 
    os.system(f"""sqlite3 {db_file} ".mode tabs" ".import {t.name} histoplasma_manuscript_samples" """)
                       
             

       
if __name__ == '__main__':
  main()
