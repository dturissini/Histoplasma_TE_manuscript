import sqlite3
import os
import tempfile




#python3 process_results_edta_only.py 



def main():  
  db_file = '/work/users/d/t/dturissi/histoplasma/te_analysis/histoplasma_tes.db'
  conn = sqlite3.connect(db_file)
    
  conn.execute("drop table if exists histoplasma_edta_tes_edta_only")
  conn.execute("""create table histoplasma_edta_tes_edta_only
                 (heeo_id int primary key,
                  sample_id varchar(50),
                  score float,
                  scaffold varchar(50),
                  scaffold_start int,
                  scaffold_end int,
                  scaffold_strand varchar(50),
                  subclass varchar(50),
                  full_intact_te varchar(1))""")
                  
  conn.execute("create index idx_heeo_sample_id on histoplasma_edta_tes_edta_only(sample_id)")
  conn.execute("create index idx_heeo_scaffold on histoplasma_edta_tes_edta_only(scaffold)")
  conn.execute("create index idx_heeo_subclass on histoplasma_edta_tes_edta_only(subclass)")

  conn.execute("drop table if exists edta_te_classes_edta_only")
  conn.execute("""create table edta_te_classes_edta_only
                 (subclass varchar(50) primary key,
                  class varchar(50))""")
                  
  conn.execute("create index idx_tceo_class on edta_te_classes_edta_only(class)")


  edta_dir = '/proj/matutelb/projects/fungus/histoplasma/TE_Analysis/EDTA_output'
  
  sample_query = conn.execute("""select sample_id 
                                 from histoplasma_manuscript_samples m, histoplasma_tree_names t
                                 where m.tree_name = t.tree_name""")
  
  subclasses = set()
  heeo_id = 1
  with tempfile.NamedTemporaryFile(mode='w') as t:  
    for sample_id, in sample_query:
      print(sample_id)
      
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
      
      gff_file = os.path.join(edta_dir, sample_id_edta, sample_id_edta + '.scaffolds.fa.mod.EDTA.TEanno.gff3')
      
      
      #get edta results
      with open(gff_file, 'r') as gff: 
        for line in gff:
          if line[0] != "#":
            line = line.strip()
            gff_vals = line.split('\t')
            
            if len(gff_vals) == 9:                 
              (chrom, source, feature_type, start, end, score, strand, phase, attributes) = gff_vals
              
              attribute_dict = process_attributes(attributes)
              te_id = attribute_dict['Name']
              
              
              subclass = attribute_dict['Classification']
              
              #convert TIR/Tc1_Mariner to DNA/DTT since they are synonymous but EDTA reports them separately
              #it was previous confirmed that there is no overlap between hits so this doesn't result in double counting
              #TIR/PiggyBac and LTR/Gypsy are simply renamed
              if subclass == 'TIR/Tc1_Mariner':
                subclass = 'DNA/DTT'
              elif subclass == 'TIR/PiggyBac':
                subclass = 'DNA/DTB'
              elif subclass == 'LTR/Gypsy':
                subclass = 'LTR/Ty3'
              
              subclasses.add(subclass)
              
              te_method = attribute_dict['Method']
              
              heeo_id += 1
              full_intact_te = 'N'
              if te_method == 'structural': 
                full_intact_te = 'Y'
                
              
              t.write(f"{heeo_id}\t{sample_id}\t{score}\t{chrom}\t{start}\t{end}\t{strand}\t{subclass}\t{full_intact_te}\n")

    t.flush() 
    os.system(f"""sqlite3 {db_file} ".mode tabs" ".import {t.name} histoplasma_edta_tes_edta_only" """)
  
    
  for subclass in subclasses:
    class_prefix = subclass.split('/')[0]
    if class_prefix in ['DNA', 'DNAauto', 'DNAnona', 'TIR', 'RC']:
      te_class = 'DNA TE'
    elif class_prefix in ['LINE', 'LTR', 'rRNA', 'SINE']:
      te_class = 'RNA TE'
    else:
      te_class = class_prefix
    
    conn.execute(f"""insert into edta_te_classes_edta_only
                     values
                     ('{subclass}', '{te_class}')""")

  conn.commit()
  conn.close()


def process_attributes(attributes):
  attributes_dict = {}
  for attribute in attributes.split(';'):
    if attribute != '' and '=' in attribute:
      (key, val) = attribute.split('=')
      attributes_dict[key] = val

  return attributes_dict  
  


if __name__ == '__main__':
  main()
