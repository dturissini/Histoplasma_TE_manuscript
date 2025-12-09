import sqlite3
import os
import sys
import tempfile

#python3 load_histoplasma_assembly_stats.py



def main():  
  db_file = '/work/users/d/t/dturissi/histoplasma/te_analysis/histoplasma_tes.db'
  
  conn = sqlite3.connect(db_file)    
  conn.execute(f"""drop table if exists histoplasma_assembly_stats""")
  conn.execute(f"""create table histoplasma_assembly_stats
                   (tree_name varchar(100) primary key, 
                    L50 float,
                    N50 float,
                    busco_complete float)""")

  conn.close()

                              
  with tempfile.NamedTemporaryFile(mode='w') as t: 
    with open('assembly_stats_tk.csv', 'r') as r:  
      next(r)
      for line in r:
        species, tree_name, genome_len, L50, N50, busco_complete, short_name = line.strip().split(",")  
             
        t.write(f"""{tree_name}\t{L50}\t{N50}\t{busco_complete}\n""")
          
    t.flush() 
    os.system(f"""sqlite3 {db_file} ".mode tabs" ".import {t.name} histoplasma_assembly_stats" """)
                       
             

       
if __name__ == '__main__':
  main()
