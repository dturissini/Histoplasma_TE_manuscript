import sqlite3
import os
import sys

#python3 make_histoplasma_te_genome_pers.py



def main():  
  db_file = '/work/users/d/t/dturissi/histoplasma/te_analysis/histoplasma_tes.db'
  
  conn = sqlite3.connect(db_file)    
  conn.execute(f"""drop table if exists histoplasma_te_genome_pers_class""")
  conn.execute(f"""create table histoplasma_te_genome_pers_class
                   (te_approach varchar(20),
                    class varchar(50),
                    sample_id varchar(100), 
                    genome_per float)""")
  
  conn.execute(f"create index idx_htgpc_te_approach on histoplasma_te_genome_pers_class(te_approach)")  
  conn.execute(f"create index idx_htgpc_sample_id on histoplasma_te_genome_pers_class(sample_id)")  
  conn.execute(f"create index idx_htgpc_class on histoplasma_te_genome_pers_class(class)")  
                    

  conn.execute(f"""drop table if exists histoplasma_te_genome_pers_subclass""")
  conn.execute(f"""create table histoplasma_te_genome_pers_subclass
                   (te_approach varchar(20),
                    subclass varchar(50),
                    sample_id varchar(100), 
                    genome_per float)""")
  
  conn.execute(f"create index idx_htgps_te_approach on histoplasma_te_genome_pers_subclass(te_approach)")  
  conn.execute(f"create index idx_htgps_sample_id on histoplasma_te_genome_pers_subclass(sample_id)")  
  conn.execute(f"create index idx_htgps_class on histoplasma_te_genome_pers_subclass(subclass)")  
                    

  get_genome_pers(conn)
  
  #get te_approaches and sample_ids
  te_approach_query = conn.execute(f"""select distinct te_approach from histoplasma_te_genome_pers_subclass""")
  
  te_approaches = []
  for te_approach, in te_approach_query:
    te_approaches.append(te_approach)

  sample_id_query = conn.execute(f"""select distinct sample_id from histoplasma_te_genome_pers_subclass""")

  sample_ids = []
  for sample_id, in sample_id_query:
    sample_ids.append(sample_id)
  
  #add zeros for missing te classes and subclasses
  add_zeros(conn, 'class', te_approaches, sample_ids)     
  add_zeros(conn, 'subclass', te_approaches, sample_ids)     
                       
             
  conn.commit()
  conn.close()


def get_genome_pers(conn):
  #EDTA-RM
  print('EDTA-RM')
  conn.execute(f"""insert into histoplasma_te_genome_pers_class
                   select 'EDTA-RM', class, e.sample_id, sum(100.0 * (scaffold_end - scaffold_start + 1) / genome_len)
                   from histoplasma_edta_tes e, edta_te_classes c, histoplasma_genome_lens l
                   where e.sample_id = l.sample_id
                   and e.subclass = c.subclass
                   group by 'EDTA-RM', class, e.sample_id""")

  conn.execute(f"""insert into histoplasma_te_genome_pers_subclass
                   select 'EDTA-RM', subclass, e.sample_id, sum(100.0 * (scaffold_end - scaffold_start + 1) / genome_len)
                   from histoplasma_edta_tes e, histoplasma_genome_lens l
                   where e.sample_id = l.sample_id
                   group by 'EDTA-RM', subclass, e.sample_id""")


  #EDTA-only
  print('EDTA-only')
  conn.execute(f"""insert into histoplasma_te_genome_pers_class
                   select 'EDTA-only', class, e.sample_id, sum(100.0 * (scaffold_end - scaffold_start + 1) / genome_len)
                   from histoplasma_edta_tes_EDTA_only e, edta_te_classes_EDTA_only c, histoplasma_genome_lens l
                   where e.sample_id = l.sample_id
                   and e.subclass = c.subclass
                   group by 'EDTA-only', class, e.sample_id""")

  conn.execute(f"""insert into histoplasma_te_genome_pers_subclass
                   select 'EDTA-only', subclass, e.sample_id, sum(100.0 * (scaffold_end - scaffold_start + 1) / genome_len)
                   from histoplasma_edta_tes_EDTA_only e, histoplasma_genome_lens l
                   where e.sample_id = l.sample_id
                   group by 'EDTA-only', subclass, e.sample_id""")
  
  
  #dnaPipeTE_RM_lib
  print('dnaPipeTE_RM_lib')
  conn.execute(f"""insert into histoplasma_te_genome_pers_class
                   select 'dnaPipeTE RM lib', RM_class, e.sample_id, sum(genome_per)
                   from histoplasma_dnaPipeTE e
                   where genome_coverage == 1
                   and sample_number == 5 
                   and RM_t == 0.2
                   group by 'dnaPipeTE RM lib', e.sample_id, RM_class""")

  conn.execute(f"""insert into histoplasma_te_genome_pers_subclass
                   select 'dnaPipeTE RM lib', RM_subclass subclass, e.sample_id, sum(genome_per)
                   from histoplasma_dnaPipeTE e
                   where genome_coverage == 1
                   and sample_number == 5 
                   and RM_t == 0.2
                   group by 'dnaPipeTE RM lib', e.sample_id, RM_subclass""")


  #dnaPipeTE_EDTA_lib
  print('dnaPipeTE_EDTA_lib')
  conn.execute(f"""insert into histoplasma_te_genome_pers_class
                   select 'dnaPipeTE EDTA lib', RM_class, e.sample_id, sum(genome_per)
                   from histoplasma_dnaPipeTE_edta_lib e
                   where genome_coverage == 1
                   and sample_number == 5 
                   and RM_t == 0.2
                   group by 'dnaPipeTE EDTA lib', e.sample_id, RM_class""")
  
  conn.execute(f"""insert into histoplasma_te_genome_pers_subclass
                   select 'dnaPipeTE EDTA lib', RM_subclass subclass, e.sample_id, sum(genome_per)
                   from histoplasma_dnaPipeTE_edta_lib e
                   where genome_coverage == 1
                   and sample_number == 5 
                   and RM_t == 0.2
                   group by 'dnaPipeTE EDTA lib', e.sample_id, RM_subclass""")
  



def add_zeros(conn, class_type, te_approaches, sample_ids):
  class_query = conn.execute(f"""select {class_type} 
                                 from histoplasma_te_genome_pers_{class_type}
                                 group by {class_type}""")
  
  for te_class, in class_query:
    print(class_type, te_class)
    for te_approach in te_approaches:
      for sample_id in sample_ids:
        te_check = conn.execute(f"""select count(*) 
                                    from histoplasma_te_genome_pers_{class_type}
                                    where {class_type} = '{te_class}'
                                    and te_approach = '{te_approach}'
                                    and sample_id = '{sample_id}'""").fetchone()[0]
        
        if te_check == 0:      
          conn.execute(f"""insert into histoplasma_te_genome_pers_{class_type}
                           values
                           ('{te_approach}', '{te_class}', '{sample_id}', 0)""")
  
       
if __name__ == '__main__':
  main()
