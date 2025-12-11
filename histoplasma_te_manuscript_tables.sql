sqlite3 /work/users/d/t/dturissi/histoplasma/te_analysis/histoplasma_tes.db
.separator ","
.header on

.output '/work/users/d/t/dturissi/histoplasma/te_analysis/tables/table_s2_isolate_genome_data.csv'

select sample_short_name Isolate, species Species, genome_len "Genome length", L50, N50, busco_complete "Busco completeness"
from histoplasma_genome_lens l, histoplasma_tree_names t, histoplasma_assembly_stats a
where l.sample_id = t.sample_id
and t.tree_name = a.tree_name
and t.tree_name in (select tree_name
                    from histoplasma_manuscript_samples)
order by sample_short_name;

.output



.output '/work/users/d/t/dturissi/histoplasma/te_analysis/tables/table_s3_dnaPipeTE_parameter_sweep.csv'

select "dnaPipeTE EDTA lib" "TE approach", genome_coverage, sample_number, RM_t, round(avg(genome_per), 1) "Average genomic percentage", round(avg(num_subclasses), 1) "Average number of superfamilies"
from (select e.sample_id, genome_coverage, sample_number, RM_t, sum(genome_per) genome_per, count(distinct RM_subclass) num_subclasses
      from histoplasma_dnaPipeTE_edta_lib e, histoplasma_genome_lens l
      where e.sample_id = l.sample_id
      group by e.sample_id, genome_coverage, sample_number, RM_t) x
group by genome_coverage, sample_number, RM_t
union all
select "dnaPipeTE RM lib" "TE approach", genome_coverage, sample_number, RM_t, round(avg(genome_per), 1) "Average genomic percentage", round(avg(num_subclasses), 1) "Average number of superfamilies"
from (select e.sample_id, genome_coverage, sample_number, RM_t, sum(genome_per) genome_per, count(distinct RM_subclass) num_subclasses
      from histoplasma_dnaPipeTE e, histoplasma_genome_lens l
      where e.sample_id = l.sample_id
      group by e.sample_id, genome_coverage, sample_number, RM_t) x
group by genome_coverage, sample_number, RM_t;

.output


.output '/work/users/d/t/dturissi/histoplasma/te_analysis/tables/table_4_te_completion_per_edta_rm.csv'

select class Class, c.subclass Superfamily, 
round(100.0 * sum(case when 1.0 * (te_end - te_start + 1) / (te_end + te_left) == 1 then 1 else 0 end) / count(*), 1) per_100,
round(100.0 * sum(case when 1.0 * (te_end - te_start + 1) / (te_end + te_left) >= .95 then 1 else 0 end) / count(*), 1) per_gt_95,
round(100.0 * sum(case when 1.0 * (te_end - te_start + 1) / (te_end + te_left) >= .9 then 1 else 0 end) / count(*), 1) per_gt_90,
round(100.0 * sum(case when 1.0 * (te_end - te_start + 1) / (te_end + te_left) >= .5 then 1 else 0 end) / count(*), 1) per_gt_50                   
from histoplasma_edta_tes d, edta_te_classes c
where d.subclass = c.subclass
and class in ('DNA TE', 'RNA TE')
group by class, c.subclass;

.output



.output '/work/users/d/t/dturissi/histoplasma/te_analysis/tables/table_s8_te_completion_per_edta_only.csv'

select class Class, c.subclass Superfamily, 
round(100.0 * sum(case when full_intact_te = 'Y' then 1 else 0 end) / count(*), 1) per_100                 
from histoplasma_edta_tes_edta_only d, edta_te_classes_edta_only c
where d.subclass = c.subclass
and class in ('DNA TE', 'RNA TE')
group by class, c.subclass;

.output



.output '/work/users/d/t/dturissi/histoplasma/te_analysis/tables/table_s9_te_completion_per_dnapipete_rm_lib.csv'

select RM_class Class, RM_subclass Superfamily, 
round(100.0 * sum(case when lib_hit_prop == 1 then 1 else 0 end) / count(*), 1) per_100,
round(100.0 * sum(case when lib_hit_prop >= .95 then 1 else 0 end) / count(*), 1) per_gt_95,
round(100.0 * sum(case when lib_hit_prop >= .9 then 1 else 0 end) / count(*), 1) per_gt_90,
round(100.0 * sum(case when lib_hit_prop >= .5 then 1 else 0 end) / count(*), 1) per_gt_50                   
from histoplasma_dnaPipeTE
where RM_class in ('DNA TE', 'RNA TE')
and genome_coverage = 1
and sample_number = 5 
and RM_t = 0.2
group by RM_class, RM_subclass;

.output



.output '/work/users/d/t/dturissi/histoplasma/te_analysis/tables/table_s10_te_completion_per_dnapipete_edta_lib.csv'

select RM_class Class, RM_subclass Superfamily, 
round(100.0 * sum(case when lib_hit_prop == 1 then 1 else 0 end) / count(*), 1) per_100,
round(100.0 * sum(case when lib_hit_prop >= .95 then 1 else 0 end) / count(*), 1) per_gt_95,
round(100.0 * sum(case when lib_hit_prop >= .9 then 1 else 0 end) / count(*), 1) per_gt_90,
round(100.0 * sum(case when lib_hit_prop >= .5 then 1 else 0 end) / count(*), 1) per_gt_50                   
from histoplasma_dnaPipeTE_edta_lib
where RM_class in ('DNA TE', 'RNA TE')
and genome_coverage = 1
and sample_number = 5 
and RM_t = 0.2
group by RM_class, RM_subclass;

.output





.output '/work/users/d/t/dturissi/histoplasma/te_analysis/tables/table_s17_te_per_genome_len_corr.csv'

select class "Class", te_approach "TE method", lm_type "Regression", round(r_2, 3) "R^2", round(p_value, 3) "P-value", round(f_statistic, 3) "F statistic", f_statistic_df "DF"
from histo_te_per_genome_len_corrs_classes
where class in ('DNA TE', 'RNA TE')
order by class, te_approach;

.output
