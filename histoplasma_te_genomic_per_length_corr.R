##library(phytools)
library(RSQLite)


base_dir <- "/work/users/d/t/dturissi/histoplasma/te_analysis"
plot_dir <- paste(base_dir, '/plots/', sep='')


setwd(base_dir)


conn <- dbConnect(dbDriver("SQLite"), 'histoplasma_tes.db')


#subclasses
rm_edta_subclass_pers <- dbGetQuery(conn, "select tree_name, subclass, round(sum(100.0 * (scaffold_end - scaffold_start + 1) / genome_len), 5) genome_per, genome_len
                                           from histoplasma_edta_tes e, histoplasma_genome_lens l, histoplasma_tree_names t
                                           where e.sample_id = l.sample_id
                                           and e.sample_id = t.sample_id
                                           and l.sample_id != 'Ep_130_s_7'
                                           group by tree_name, subclass, genome_len")


edta_only_subclass_pers <- dbGetQuery(conn, "select tree_name, subclass, round(sum(100.0 * (scaffold_end - scaffold_start + 1) / genome_len), 5) genome_per, genome_len
                                             from histoplasma_edta_tes_edta_only e, histoplasma_genome_lens l, histoplasma_tree_names t
                                             where e.sample_id = l.sample_id
                                             and e.sample_id = t.sample_id
                                             and l.sample_id != 'Ep_130_s_7'
                                             group by tree_name, subclass, genome_len")

dp_edta_lib_subclass_pers <- dbGetQuery(conn, "select tree_name, RM_subclass subclass, round(sum(100.0 * contig_len / genome_len), 5) genome_per, genome_len
                                               from histoplasma_dnaPipeTE_edta_lib e, histoplasma_genome_lens l, histoplasma_tree_names t
                                               where e.sample_id = l.sample_id
                                               and e.sample_id = t.sample_id
                                               and genome_coverage == 1
                                               and sample_number == 5 
                                               and RM_t == 0.2
                                               and l.sample_id != 'Ep_130_s_7'
                                               group by tree_name, RM_subclass, genome_len")


dp_rm_lib_subclass_pers <- dbGetQuery(conn, "select tree_name, RM_subclass subclass, round(sum(100.0 * contig_len / genome_len), 5) genome_per, genome_len
                                             from histoplasma_dnaPipeTE e, histoplasma_genome_lens l, histoplasma_tree_names t
                                             where e.sample_id = l.sample_id
                                             and e.sample_id = t.sample_id
                                             and genome_coverage == 1
                                             and sample_number == 5 
                                             and RM_t == 0.2
                                             and l.sample_id != 'Ep_130_s_7'
                                             group by tree_name, RM_subclass, genome_len")


#classes
rm_edta_class_pers <- dbGetQuery(conn, "select tree_name, class, round(sum(100.0 * (scaffold_end - scaffold_start + 1) / genome_len), 5) genome_per, genome_len
                                        from histoplasma_edta_tes e, edta_te_classes c, histoplasma_genome_lens l, histoplasma_tree_names t
                                        where e.sample_id = l.sample_id
                                        and e.sample_id = t.sample_id
                                        and e.subclass = c.subclass
                                        and l.sample_id != 'Ep_130_s_7'
                                        group by tree_name, class, genome_len")


edta_only_class_pers <- dbGetQuery(conn, "select tree_name, class, round(sum(100.0 * (scaffold_end - scaffold_start + 1) / genome_len), 5) genome_per, genome_len
                                          from histoplasma_edta_tes_edta_only e, edta_te_classes_edta_only c, histoplasma_genome_lens l, histoplasma_tree_names t
                                          where e.sample_id = l.sample_id
                                          and e.sample_id = t.sample_id
                                          and e.subclass = c.subclass
                                          and l.sample_id != 'Ep_130_s_7'
                                          group by tree_name, class, genome_len")

dp_edta_lib_class_pers <- dbGetQuery(conn, "select tree_name, RM_class class, round(sum(100.0 * contig_len / genome_len), 5) genome_per, genome_len
                                            from histoplasma_dnaPipeTE_edta_lib e, histoplasma_genome_lens l, histoplasma_tree_names t
                                            where e.sample_id = l.sample_id
                                            and e.sample_id = t.sample_id
                                            and genome_coverage == 1
                                            and sample_number == 5 
                                            and RM_t == 0.2
                                            and l.sample_id != 'Ep_130_s_7'
                                            group by tree_name, RM_class, genome_len")


dp_rm_lib_class_pers <- dbGetQuery(conn, "select tree_name, RM_class class, round(sum(100.0 * contig_len / genome_len), 5) genome_per, genome_len
                                          from histoplasma_dnaPipeTE e, histoplasma_genome_lens l, histoplasma_tree_names t
                                          where e.sample_id = l.sample_id
                                          and e.sample_id = t.sample_id
                                          and genome_coverage == 1
                                          and sample_number == 5 
                                          and RM_t == 0.2
                                          and l.sample_id != 'Ep_130_s_7'
                                          group by tree_name, RM_class, genome_len")



subclasses <- sort(unique(c(rm_edta_subclass_pers$subclass, edta_only_subclass_pers$subclass, dp_edta_lib_subclass_pers$subclass)))
dp_rm_lib_subclasses <- sort(unique(dp_rm_lib_subclass_pers$subclass))
classes <- sort(unique(c(rm_edta_class_pers$class, edta_only_class_pers$class, dp_edta_lib_class_pers$class, dp_rm_lib_class_pers$class)))
                                
class_genomic_corr_plot <- function(te_type, class, panel_name, te_pers, genome_lens)
  {
  if (length(te_pers) > 0)
    {  
    dpt_lm <- lm(te_pers ~ genome_lens)
    dpt_lm_r2 <- signif(summary(dpt_lm)$r.squared, 3)
  
    if (nrow(coef(summary(dpt_lm))) == 2)
    {
    dpt_lm_p <- signif(coef(summary(dpt_lm))[2,4], 3)  
    
    if (! is.nan(dpt_lm_p) && dpt_lm_p < 0.01)
      {
      dpt_lm_p <- format(dpt_lm_p, scientific = TRUE)
      p_main <- gsub("e.*", "", dpt_lm_p)
      p_exponent <- gsub(".*e", "", dpt_lm_p)
      dpt_lm_p <- bquote(.(p_main) %*% 10^.(as.numeric(p_exponent)))
      }
    } else {
    dpt_lm_p <- NA  
    }
    
    plot(genome_lens, te_pers, pch=20, col='black', xlab='Genome length', ylab=paste(class, ' genomic %', sep=''), main=bquote(atop(.(te_type) ~ ': ' ~ .(class),  R^2 == .(dpt_lm_r2) ~ ', p =' ~ .(dpt_lm_p))))
    mtext(panel_name, side=3, line=1.25, at= min(genome_lens) - 0.1 * (max(genome_lens) - min(genome_lens)), cex=1.5)


    if (length(te_pers) > 1)
      {abline(dpt_lm$coefficients)}
    } else {
    plot(1, type='n', xlab='Genome length', ylab=paste(class, ' genomic %', sep=''), main=c(paste(te_type, ': ', class, sep=''), 'No TEs'))
    } 
  }




pdf(paste(plot_dir, "histoplasma_te_genomic_per_length_corr_class.pdf", sep=''), height=12, width=12)
par(mfrow=c(2,2))
for (class in classes)
  {  
  class_genomic_corr_plot('EDTA-RM', class, '', rm_edta_class_pers$genome_per[rm_edta_class_pers$class == class], rm_edta_class_pers$genome_len[rm_edta_class_pers$class == class])    
  class_genomic_corr_plot('EDTA-only', class, '', edta_only_class_pers$genome_per[edta_only_class_pers$class == class], edta_only_class_pers$genome_len[edta_only_class_pers$class == class])    
  class_genomic_corr_plot('dnaPipeTE EDTA lib', class, '', dp_edta_lib_class_pers$genome_per[dp_edta_lib_class_pers$class == class], dp_edta_lib_class_pers$genome_len[dp_edta_lib_class_pers$class == class])    
  class_genomic_corr_plot('dnaPipeTE RM lib', class, '', dp_rm_lib_class_pers$genome_per[dp_rm_lib_class_pers$class == class], dp_rm_lib_class_pers$genome_len[dp_rm_lib_class_pers$class == class])    
  }
par(mfcol=c(1,1))
dev.off()



pdf(paste(plot_dir, "histoplasma_te_genomic_per_length_corr_subclass.pdf", sep=''), height=12, width=12)
par(mfrow=c(2,2))
for (subclass in subclasses)
  {  
  class_genomic_corr_plot('EDTA-RM', subclass, '', rm_edta_subclass_pers$genome_per[rm_edta_subclass_pers$subclass == subclass], rm_edta_subclass_pers$genome_len[rm_edta_subclass_pers$subclass == subclass])    
  class_genomic_corr_plot('EDTA-only', subclass, '', edta_only_subclass_pers$genome_per[edta_only_subclass_pers$subclass == subclass], edta_only_subclass_pers$genome_len[edta_only_subclass_pers$subclass == subclass])    
  class_genomic_corr_plot('dnaPipeTE EDTA lib', subclass, '', dp_edta_lib_subclass_pers$genome_per[dp_edta_lib_subclass_pers$subclass == subclass], dp_edta_lib_subclass_pers$genome_len[dp_edta_lib_subclass_pers$subclass == subclass])    
  class_genomic_corr_plot('dnaPipeTE RM lib', subclass, '', dp_rm_lib_subclass_pers$genome_per[dp_rm_lib_subclass_pers$subclass == subclass], dp_rm_lib_subclass_pers$genome_len[dp_rm_lib_subclass_pers$subclass == subclass])    
  }
par(mfcol=c(1,1))
dev.off()



pdf(paste(plot_dir, "manuscript_figure_10_te_genomic_per_length_corr_class.pdf", sep=''), height=16, width=8)
par(mfrow=c(4,2))
class_genomic_corr_plot('EDTA-RM', 'RNA TE', 'A)', rm_edta_class_pers$genome_per[rm_edta_class_pers$class == 'RNA TE'], rm_edta_class_pers$genome_len[rm_edta_class_pers$class == 'RNA TE'])    
class_genomic_corr_plot('EDTA-RM', 'DNA TE', 'B)', rm_edta_class_pers$genome_per[rm_edta_class_pers$class == 'DNA TE'], rm_edta_class_pers$genome_len[rm_edta_class_pers$class == 'DNA TE'])    

class_genomic_corr_plot('EDTA-only', 'RNA TE', 'C)', edta_only_class_pers$genome_per[edta_only_class_pers$class == 'RNA TE'], edta_only_class_pers$genome_len[edta_only_class_pers$class == 'RNA TE'])    
class_genomic_corr_plot('EDTA-only', 'DNA TE', 'D)', edta_only_class_pers$genome_per[edta_only_class_pers$class == 'DNA TE'], edta_only_class_pers$genome_len[edta_only_class_pers$class == 'DNA TE'])    

class_genomic_corr_plot('dnaPipeTE RM lib', 'RNA TE', 'E)', dp_rm_lib_class_pers$genome_per[dp_rm_lib_class_pers$class == 'RNA TE'], dp_rm_lib_class_pers$genome_len[dp_rm_lib_class_pers$class == 'RNA TE'])    
class_genomic_corr_plot('dnaPipeTE RM lib', 'DNA TE', 'F)', dp_rm_lib_class_pers$genome_per[dp_rm_lib_class_pers$class == 'DNA TE'], dp_rm_lib_class_pers$genome_len[dp_rm_lib_class_pers$class == 'DNA TE'])    

class_genomic_corr_plot('dnaPipeTE EDTA lib', 'RNA TE', 'G)', dp_edta_lib_class_pers$genome_per[dp_edta_lib_class_pers$class == 'RNA TE'], dp_edta_lib_class_pers$genome_len[dp_edta_lib_class_pers$class == 'RNA TE'])    
class_genomic_corr_plot('dnaPipeTE EDTA lib', 'DNA TE', 'H)', dp_edta_lib_class_pers$genome_per[dp_edta_lib_class_pers$class == 'DNA TE'], dp_edta_lib_class_pers$genome_len[dp_edta_lib_class_pers$class == 'DNA TE'])    
par(mfcol=c(1,1))
dev.off()

