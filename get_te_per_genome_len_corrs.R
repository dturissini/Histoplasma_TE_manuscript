##library(phytools)
library(RSQLite)
library(caper)
library(multcomp)
library(car)


base_dir <- "/work/users/d/t/dturissi/histoplasma/te_analysis"
plot_dir <- paste(base_dir, '/plots/', sep='')


setwd(base_dir)


conn <- dbConnect(dbDriver("SQLite"), 'histoplasma_tes.db')

dbSendQuery(conn, "drop table if exists histo_te_per_genome_len_corrs_classes")
dbSendQuery(conn, "create table histo_te_per_genome_len_corrs_classes
                   (te_approach varchar(20),
                    class varchar(50),
                    lm_type varchar(50),
                    r_2 float,
                    p_value float,
                    f_statistic float,
                    f_statistic_df int)")

dbSendQuery(conn, "create index idx_htpglcc_class on histo_te_per_genome_len_corrs_classes(class)")  

#classes
rm_edta_class_pers <- dbGetQuery(conn, "select t.tree_name, class, round(sum(100.0 * (scaffold_end - scaffold_start + 1) / genome_len), 5) genome_per, genome_len, busco_complete
                                        from histoplasma_edta_tes e, edta_te_classes c, histoplasma_genome_lens l, histoplasma_tree_names t, histoplasma_assembly_stats a
                                        where e.sample_id = l.sample_id
                                        and e.sample_id = t.sample_id
                                        and e.subclass = c.subclass
                                        and t.tree_name = a.tree_name
                                        group by t.tree_name, class, genome_len, busco_complete")


edta_only_class_pers <- dbGetQuery(conn, "select t.tree_name, class, round(sum(100.0 * (scaffold_end - scaffold_start + 1) / genome_len), 5) genome_per, genome_len, busco_complete
                                          from histoplasma_edta_tes_edta_only e, edta_te_classes_edta_only c, histoplasma_genome_lens l, histoplasma_tree_names t, histoplasma_assembly_stats a
                                          where e.sample_id = l.sample_id
                                          and e.sample_id = t.sample_id
                                          and e.subclass = c.subclass
                                          and t.tree_name = a.tree_name
                                          group by t.tree_name, class, genome_len, busco_complete")

dp_edta_lib_class_pers <- dbGetQuery(conn, "select t.tree_name, RM_class class, round(sum(100.0 * contig_len / genome_len), 5) genome_per, genome_len, busco_complete
                                            from histoplasma_dnaPipeTE_edta_lib e, histoplasma_genome_lens l, histoplasma_tree_names t, histoplasma_assembly_stats a
                                            where e.sample_id = l.sample_id
                                            and e.sample_id = t.sample_id
                                            and genome_coverage == 1
                                            and sample_number == 5 
                                            and RM_t == 0.2
                                            and t.tree_name = a.tree_name
                                            group by t.tree_name, RM_class, genome_len, busco_complete")


dp_rm_lib_class_pers <- dbGetQuery(conn, "select t.tree_name, RM_class class, round(sum(100.0 * contig_len / genome_len), 5) genome_per, genome_len, busco_complete
                                          from histoplasma_dnaPipeTE e, histoplasma_genome_lens l, histoplasma_tree_names t, histoplasma_assembly_stats a
                                          where e.sample_id = l.sample_id
                                          and e.sample_id = t.sample_id
                                          and genome_coverage == 1
                                          and sample_number == 5 
                                          and RM_t == 0.2
                                          and t.tree_name = a.tree_name
                                          group by t.tree_name, RM_class, genome_len, busco_complete")



classes <- sort(unique(c(rm_edta_class_pers$class, edta_only_class_pers$class, dp_edta_lib_class_pers$class, dp_rm_lib_class_pers$class)))
                                



get_class_genomic_lm <- function(te_type, class, te_per_df)
  {
  if (nrow(te_per_df) > 1)
    {  
    #lm and lm with busco completeness covariate
    dpt_lm <- lm(genome_per ~ genome_len, subset(te_per_df, tree_name != 'Ep_130_s_7'))
    dpt_lm_r2 <- signif(summary(dpt_lm)$r.squared, 3)
    dpt_lm_p <- signif(coef(summary(dpt_lm))[2,4], 3)  

    dpt_lm_busco_covar <- lm(genome_per ~ genome_len * busco_complete, subset(te_per_df, tree_name != 'Ep_130_s_7'))
    dpt_lm_busco_covar_r2 <- signif(summary(dpt_lm_busco_covar)$r.squared, 3)
    dpt_lm_busco_covar_p <- signif(coef(summary(dpt_lm_busco_covar))[2,4], 3)  
    
    if (is.nan(dpt_lm_p))
      {dpt_lm_p <- 0}
        
    if (is.nan(dpt_lm_busco_covar_p))
      {dpt_lm_busco_covar_p <- 0}
        
    dbSendQuery(conn, paste("insert into histo_te_per_genome_len_corrs_classes
                         values
                         ('", te_type, "', '", class, "', 'lm', ",  dpt_lm_r2, ", ", dpt_lm_p, ", null, null)", sep=''))

    dbSendQuery(conn, paste("insert into histo_te_per_genome_len_corrs_classes
                         values
                         ('", te_type, "', '", class, "', 'lm Busco complete covar', ",  dpt_lm_busco_covar_r2, ", ", dpt_lm_busco_covar_p, ", null, null)", sep=''))
    
    
    #PGLS
    histo.tree <- read.tree("histoplasma_phylo_tree_africa.txt")
  
  
    ## match the data to the tree
    droptips <- setdiff(histo.tree$tip.label, te_per_df$tree_name)
    dropdat <- intersect(te_per_df$tree_name, histo.tree$tip.label)
    smalltree <- drop.tip(histo.tree, droptips)
    trimdat <- te_per_df[te_per_df$tree_name %in% dropdat,]
  

    rooted_small_tree <- root(smalltree, outgroup= 'Ep_130_s_7', resolve.root = TRUE)
    rooted_small_tree$node.label <- "" #removing node labels - these were bootstrap support values with lots of repeats\  
    histo.compdat <- comparative.data(rooted_small_tree, trimdat, names.col = "tree_name", force.root = TRUE, na.omit=F)
    te_per_pgls <- pgls(genome_len ~ genome_per, data=histo.compdat, lambda="ML")  
    
    pgls_r2 <- summary(te_per_pgls)$r.squared
    pgls_p <- summary(te_per_pgls)$coefficients[2,4]
    pgls_f_stat <- summary(te_per_pgls)$fstatistic[1]
    pgls_f_stat_df <- summary(te_per_pgls)$df[2]   
    
    dbSendQuery(conn, paste("insert into histo_te_per_genome_len_corrs_classes
                         values
                         ('", te_type, "', '", class, "', 'PGLS', ",  pgls_r2, ", ", pgls_p, ", ", pgls_f_stat, ", ", pgls_f_stat_df, ")", sep=''))         
    }
  }



  
for (class_i in classes)
  {
  get_class_genomic_lm('EDTA-RM', class_i, subset(rm_edta_class_pers, class == class_i))
  get_class_genomic_lm('EDTA-only', class_i, subset(edta_only_class_pers, class == class_i))   
  get_class_genomic_lm('dnaPipeTE EDTA lib', class_i, subset(dp_edta_lib_class_pers, class == class_i))    
  get_class_genomic_lm('dnaPipeTE RM lib', class_i, subset(dp_rm_lib_class_pers, class == class_i))   
  }
  