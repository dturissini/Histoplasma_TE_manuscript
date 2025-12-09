import sqlite3
import os
import tempfile
import pandas as pd



#python3 process_results_edta_rm.py 



def main():  
  db_file = '/work/users/d/t/dturissi/histoplasma/te_analysis/histoplasma_tes.db'
  conn = sqlite3.connect(db_file)
    
  conn.execute("drop table if exists histoplasma_edta_tes")
  conn.execute("""create table histoplasma_edta_tes
                 (he_id int primary key,
                  sample_id varchar(50),
                  score float,
                  per_div float,
                  per_del float,
                  per_ins float,
                  scaffold varchar(50),
                  scaffold_start int,
                  scaffold_end int,
                  scaffold_left varchar(50),
                  scaffold_strand varchar(50),
                  te_family varchar(50),
                  subclass varchar(50),
                  te_start varchar(50),
                  te_end varchar(50),
                  te_left varchar(50),
                  te_id varchar(50))""")
                  
  conn.execute("create index idx_pt_sample_id on histoplasma_edta_tes(sample_id)")
  conn.execute("create index idx_pt_scaffold on histoplasma_edta_tes(scaffold)")

  conn.execute("drop table if exists edta_te_classes")
  conn.execute("""create table edta_te_classes
                 (subclass varchar(50) primary key,
                  class varchar(50))""")
                  
  conn.execute("create index idx_tc_class on edta_te_classes(class)")


  rm_edta_dir = '/proj/matutelb/projects/fungus/histoplasma/TE_Analysis/RM_EDTA'
  
  sample_query = conn.execute("""select sample_id 
                                 from histoplasma_manuscript_samples m, histoplasma_tree_names t
                                 where m.tree_name = t.tree_name""")
  
  subclasses = set()
  he_id_start = 1
  for sample_id, in sample_query:
    if sample_id == '104_P_19_S5':
      sample_id_edta = 'Histoplasma_capsulatum_104_P_19'
    elif sample_id == 'Dr_Anuradha_Fungal_WGS_S14':
      sample_id_edta = 'Histoplasma_capsulatum_Dr_Anuradha_Fungal_WGS_S11'
    elif sample_id == 'HISSP-FGBON2001-xx-CL-GUF-xxxx-036-BB_S21':
      sample_id_edta = 'Histoplasma_capsulatum_HISSP-FGBON2001'
    elif sample_id == 'HISSP-FGAMA2041-uu-CL-GUF-xxxx-036-BB_S12_L001':
      sample_id_edta = 'Histoplasma_capsulatum_HISSP-FGAMA2041'
    elif sample_id == 'HISSP-FGPOE2043-uu-CL-GUF-xxxx-036-BB_S4_L001':
      sample_id_edta = 'Histoplasma_capsulatum_HISSP-FGPOE2043'
    elif sample_id == 'NACVFR_Histo_HC4137_S92_L008':
      sample_id_edta = 'Histoplasma_capsulatum_NACVFR_Histo_HC4137'
    elif sample_id == 'NACVFR_Histo_HC970591_S81_L007':
      sample_id_edta = 'Histoplasma_capsulatum_NACVFR_Histo_HC970591'
    elif sample_id == 'HISSP-SDOR2042-xx-CL-SUR-xxxx-036-BB_S14':
      sample_id_edta = 'Histoplasma_capsulatum_HISSP-SDOR2042'
    elif sample_id == 'Ep_130_s_7':
      sample_id_edta = 'Blastomyces_parvus_UAMH130'
    else:
      sample_id_edta = sample_id  
    
    results_file = os.path.join(rm_edta_dir, sample_id_edta + '.RM', sample_id_edta + '.scaffolds.fa.out')
    if os.path.exists(results_file):
      te_df = pd.read_csv(results_file, sep='\s+', header=None, skiprows=3, usecols=list(range(15)))  
      te_df.columns = ['score', 'per_div', 'per_del', 'per_ins', 'scaffold', 'scaffold_start', 'scaffold_end', 'scaffold_left', 'scaffold_strand', 'te_family', 'subclass', 'te_start', 'te_end', 'te_left', 'te_id']
      
      #remove parentheses
      te_df['te_start'] = te_df['te_start'].astype(str)
      te_df['te_start'] = te_df['te_start'].str.replace('(', '')
      te_df['te_start'] = te_df['te_start'].str.replace(')', '')
      te_df['te_start'] = te_df['te_start'].astype(int)

      te_df['te_end'] = te_df['te_end'].astype(str)
      te_df['te_end'] = te_df['te_end'].str.replace('(', '')
      te_df['te_end'] = te_df['te_end'].str.replace(')', '')
      te_df['te_end'] = te_df['te_end'].astype(int)

      te_df['te_left'] = te_df['te_left'].astype(str)
      te_df['te_left'] = te_df['te_left'].str.replace('(', '')
      te_df['te_left'] = te_df['te_left'].str.replace(')', '')
      te_df['te_left'] = te_df['te_left'].astype(int)
      
      #add he_id unique key
      he_id_end = len(te_df.index) + he_id_start
      he_ids = list(range(he_id_start, he_id_end))
      he_id_start = he_id_end    
      te_df.insert(1, "he_id", he_ids, True)
      
      #add sample_id    
      te_df.insert(2, "sample_id", [sample_id] * len(te_df.index), True)
      
      #convert TIR/Tc1_Mariner to DNA/DTT since they are synonymous but EDTA reports them separately
      #it was previous confirmed that there is no overlap between hits so this doesn't result in double counting
      #TIR/PiggyBac and LTR/Gypsy are simply renamed
      te_df.loc[te_df["subclass"] == 'TIR/Tc1_Mariner', "subclass"] = 'DNA/DTT'
      te_df.loc[te_df["subclass"] == 'TIR/PiggyBac', "subclass"] = 'DNA/DTB'
      te_df.loc[te_df["subclass"] == 'LTR/Gypsy', "subclass"] = 'LTR/Ty3'
          
      #load dataframe into db table
      te_df.to_sql('histoplasma_edta_tes', if_exists = 'append', index=False, con=conn)
      
      #get subclasses
      sample_subclasses = set(te_df['subclass'])
      subclasses = subclasses.union(sample_subclasses)
    else:
      print(f"{sample_id} does not exist")
  
  for subclass in subclasses:
    class_prefix = subclass.split('/')[0]
    if class_prefix in ['DNA', 'DNAauto', 'DNAnona', 'TIR', 'RC']:
      te_class = 'DNA TE'
    elif class_prefix in ['LINE', 'LTR', 'rRNA', 'SINE']:
      te_class = 'RNA TE'
    else:
      te_class = class_prefix
    
    conn.execute(f"""insert into edta_te_classes
                     values
                     ('{subclass}', '{te_class}')""")
  
  conn.commit()


if __name__ == '__main__':
  main()
