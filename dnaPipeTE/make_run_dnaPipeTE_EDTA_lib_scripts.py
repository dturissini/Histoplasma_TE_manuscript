import sqlite3
import os
import sys


#python3 make_run_dnaPipeTE_scripts.py all_1_5_02



def main():  
  run_prefix = sys.argv[1]
  
  base_dir = '/work/users/d/t/dturissi/histoplasma/dnaPipeTE'
  output_base = 'output_edta_lib'
  mnt_output_dir = os.path.join('/mnt', output_base)
  longleaf_output_dir = os.path.join(base_dir, output_base)
  dnaPipeTE_sh_dir = os.path.join(base_dir, 'dnaPipeTE_sh_scripts_edta_libs')
  master_sh_dir = os.path.join(base_dir, 'master_sh_scripts_edta_libs')
  log_dir = os.path.join(base_dir, 'logs')
   
  os.system(f"mkdir -p {longleaf_output_dir}")
  os.system(f"mkdir -p {dnaPipeTE_sh_dir}")
  os.system(f"mkdir -p {master_sh_dir}")
  os.system(f"mkdir -p {log_dir}")
  
  db_file = '/work/users/d/t/dturissi/histoplasma/te_analysis/histoplasma_tes.db'
  master_sh_file = os.path.join(master_sh_dir, run_prefix + '_run_dnaPipeTE_commands_edta_libs.sh')
  sbatch_file = os.path.join(base_dir, run_prefix + '_run_dnaPipeTE_edta_libs.sbatch')
  edta_lib_base_dir = '/proj/matutelb/projects/fungus/histoplasma/TE_Analysis/EDTA_output'
  
  threads = 16
  image_file = '/work/users/d/t/dturissi/software/dnaPipeTE/dnapipete.img'

  
  conn = sqlite3.connect(db_file)    
  sample_query = conn.execute(f"""select f.sample_id, genome_len, fastq_file
                                  from histoplasma_fastq f, histoplasma_genome_lens l
                                  where f.sample_id = l.sample_id
                                  and f.sample_id in (select sample_id 
                                                      from histoplasma_manuscript_samples m, histoplasma_tree_names t
                                                      where m.tree_name = t.tree_name)""")
      
  genome_coverages = ['0.1', '0.5', '2']  
  sample_numbers = ['5']  
  RM_ts = ['0.2']
                              
  total_cmds = 0
  with open(master_sh_file, 'w') as dm:
    for sample_id, genome_len, fastq_file in sample_query:
      sample_id_edta = sample_id
      
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
      
      edta_lib_file = os.path.join(edta_lib_base_dir, sample_id_edta, sample_id_edta + '.scaffolds.fa.mod.EDTA.TElib.fa')
      
      
      for genome_coverage in genome_coverages:
        for sample_number in sample_numbers: 
          for RM_t in RM_ts:
            total_cmds += 1
            
            sample_sh_dir = os.path.join(dnaPipeTE_sh_dir, sample_id)
            sample_longleaf_output_dir = os.path.join(longleaf_output_dir, sample_id + '_' + genome_coverage + '_' + sample_number + '_' + RM_t)
            sample_mnt_output_dir = os.path.join(mnt_output_dir, sample_id + '_' + genome_coverage + '_' + sample_number + '_' + RM_t)

            os.system(f"mkdir -p {sample_sh_dir}")
            os.system(f"mkdir -p {sample_longleaf_output_dir}")
            
            
            sample_run_dnaPipeTE_sh_file = os.path.join(sample_sh_dir, sample_id + '_run_dnaPipeTE_' + '_' + genome_coverage + '_' + sample_number + '_' + RM_t + '.sh')                      
            sample_apptainer_sh_file = os.path.join(sample_sh_dir, sample_id + '_apptainer_exec_' + '_' + genome_coverage + '_' + sample_number + '_' + RM_t + '.sh')        
                                   
            with open(sample_run_dnaPipeTE_sh_file, 'w') as r:    
              r.write(f"#!/bin/bash\n")
              r.write(f"cd /opt/dnaPipeTE\n")
              r.write(f"""python3 dnaPipeTE.py -input {fastq_file} -output {sample_mnt_output_dir} -RM_lib {edta_lib_file} -genome_size {genome_len} -genome_coverage {genome_coverage} -sample_number {sample_number} -RM_t {RM_t} -cpu {threads}\n""")
            
            dm.write(f"{sample_apptainer_sh_file}\n")        
            with open(sample_apptainer_sh_file, 'w') as a:    
              a.write(f"#!/bin/bash\n")
              a.write(f"""apptainer exec --bind {base_dir}:/mnt {image_file} {sample_run_dnaPipeTE_sh_file}\n""")    
            
                
            os.system(f"chmod +x {sample_run_dnaPipeTE_sh_file}")
            os.system(f"chmod +x {sample_apptainer_sh_file}")



  with open(sbatch_file, 'w') as o:
    o.write(f"#!/bin/bash\n")
    o.write(f"#SBATCH -p general\n")
    o.write(f"#sbatch -J dnaPipeTE\n")
    o.write(f"#SBATCH -t 240:00:00\n")
    o.write(f"#SBATCH --mem=64g\n")
    o.write(f"#SBATCH --cpus-per-task={threads}\n")
    o.write(f"#SBATCH --array=1-{total_cmds}%100\n")
    o.write(f"#SBATCH -o {log_dir}/run_dnaPipeTE_edta_lib_%A_%a.out\n")
    o.write(f"#SBATCH -e {log_dir}/run_dnaPipeTE_edta_lib_%A_%a.err\n\n")
    
    o.write(f"ml apptainer\n")
    o.write(f"ml repeatmasker/4.1.5\n")

    o.write(f"apptainer_sh=$(sed -n ${{SLURM_ARRAY_TASK_ID}}p {master_sh_file})\n")
    o.write(f"""eval "${{apptainer_sh}}"\n""")




       
if __name__ == '__main__':
  main()
