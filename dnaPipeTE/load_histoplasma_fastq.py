import sqlite3
import os
import sys
import tempfile

#python3 load_histoplasma_fastq.py



def main():  
  db_file = '/work/users/d/t/dturissi/histoplasma/te_analysis/histoplasma_tes.db'
  
  conn = sqlite3.connect(db_file)    
  conn.execute(f"""drop table if exists histoplasma_fastq""")
  conn.execute(f"""create table histoplasma_fastq
                   (sample_id varchar(100) primary key, 
                    fastq_file varchar(255))""")

  conn.close()
  
                              
  with tempfile.NamedTemporaryFile(mode='w') as t_fastq: 
    with open('histoplasma_te_fastqs.txt', 'r') as r:  
      for line in r:
        sample_id, fastq_file = line.strip().split("\t")          
        t_fastq.write(f"""{sample_id}\t{fastq_file}\n""")
          
          
    t_fastq.flush() 
    os.system(f"""sqlite3 {db_file} ".mode tabs" ".import {t_fastq.name} histoplasma_fastq" """)

                       
             

       
if __name__ == '__main__':
  main()
