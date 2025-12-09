import sqlite3
import os
import sys
import tempfile

#python3 get_edta_te_overlaps.py


def get_edta_overlaps(conn, te_approach, te_table, te_table_id):
  print(te_approach, 'A')
  conn.execute(f"""insert into edta_dna_dtt_tir_tc1_mariner_overlaps
                   select 'A', '{te_approach}', d.{te_table_id}, t.{te_table_id}, d.sample_id, d.scaffold, t.scaffold_end - d.scaffold_start + 1 overlap_len
                   from {te_table} d, {te_table} t
                   where d.subclass = 'DNA/DTT'
                   and t.subclass = 'TIR/Tc1_Mariner'
                   and d.scaffold = t.scaffold
                   and d.scaffold_start between t.scaffold_start and t.scaffold_end
                   and d.scaffold_end > t.scaffold_end
                   and d.sample_id = t.sample_id""")

  print(te_approach, 'B')
  conn.execute(f"""insert into edta_dna_dtt_tir_tc1_mariner_overlaps
                   select 'B', '{te_approach}', d.{te_table_id}, t.{te_table_id}, d.sample_id, d.scaffold, d.scaffold_end - t.scaffold_start + 1 overlap_len
                   from {te_table} d, {te_table} t
                   where d.subclass = 'DNA/DTT'
                   and t.subclass = 'TIR/Tc1_Mariner'
                   and d.scaffold = t.scaffold
                   and d.scaffold_end between t.scaffold_start and t.scaffold_end
                   and d.scaffold_start < t.scaffold_start
                   and d.sample_id = t.sample_id""")

  print(te_approach, 'C')
  conn.execute(f"""insert into edta_dna_dtt_tir_tc1_mariner_overlaps
                   select 'C', '{te_approach}', d.{te_table_id}, t.{te_table_id}, d.sample_id, d.scaffold, d.scaffold_end - d.scaffold_start + 1 overlap_len
                   from {te_table} d, {te_table} t
                   where d.subclass = 'DNA/DTT'
                   and t.subclass = 'TIR/Tc1_Mariner'
                   and d.scaffold = t.scaffold
                   and t.scaffold_start < d.scaffold_start
                   and t.scaffold_end > d.scaffold_end
                   and d.sample_id = t.sample_id""")

  print(te_approach, 'D')
  conn.execute(f"""insert into edta_dna_dtt_tir_tc1_mariner_overlaps
                   select 'D', '{te_approach}', d.{te_table_id}, t.{te_table_id}, d.sample_id, d.scaffold, t.scaffold_end - t.scaffold_start + 1 overlap_len
                   from {te_table} d, {te_table} t
                   where d.subclass = 'DNA/DTT'
                   and t.subclass = 'TIR/Tc1_Mariner'
                   and d.scaffold = t.scaffold
                   and d.scaffold_start <= t.scaffold_start
                   and d.scaffold_end >= t.scaffold_end
                   and d.sample_id = t.sample_id
                   and not exists (select 'x'
                                   from edta_dna_dtt_tir_tc1_mariner_overlaps eo
                                   where dd_id = d.{te_table_id}
                                   and ttm_id = t.{te_table_id})""")



def main():  
  db_file = '/work/users/d/t/dturissi/histoplasma/te_analysis/histoplasma_tes.db'
  
  conn = sqlite3.connect(db_file)    
  conn.execute(f"""drop table if exists edta_dna_dtt_tir_tc1_mariner_overlaps""")
  conn.execute(f"""create table edta_dna_dtt_tir_tc1_mariner_overlaps
                   (overlap_query varchar(1),
                   te_approach varchar(20),
                   dd_id int,
                   ttm_id int,
                   sample_id varchar(50),
                   scaffold varchar(50),
                   overlap_len int)""")

  conn.execute(f"""create index idx_ddttm_dd_id on edta_dna_dtt_tir_tc1_mariner_overlaps(dd_id)""")
  conn.execute(f"""create index idx_ddttm_ttm_id on edta_dna_dtt_tir_tc1_mariner_overlaps(ttm_id)""")
  conn.execute(f"""create index idx_ddttm_sample_id on edta_dna_dtt_tir_tc1_mariner_overlaps(sample_id)""")
  
  get_edta_overlaps(conn, 'EDTA-only', 'histoplasma_edta_tes_edta_only', 'heeo_id')
  get_edta_overlaps(conn, 'EDTA-RM', 'histoplasma_edta_tes', 'he_id')
  
  conn.commit()
  conn.close()

        
                              
             

       
if __name__ == '__main__':
  main()




  
