##library(phytools)
library(RSQLite)



base_dir <- "/work/users/d/t/dturissi/histoplasma/te_analysis"
plot_dir <- paste(base_dir, '/plots/', sep='')


setwd(base_dir)


conn <- dbConnect(dbDriver("SQLite"), 'histoplasma_tes.db')


erm_hit_prop <- dbGetQuery(conn, "select distinct sample_id, class, 
                                   case when scaffold_strand = '+' then 1.0 * (te_end - te_start + 1) / (te_end + te_left) else 1.0 * (te_end - te_left + 1) / (te_start + te_end) end hit_prop                     
                                   from histoplasma_edta_tes d, edta_te_classes c
                                   where d.subclass = c.subclass
                                   and sample_id != 'Ep_130_s_7'")

dpe_hit_prop <- dbGetQuery(conn, "select distinct sample_id, RM_class class, lib_hit_prop hit_prop                     
                                  from histoplasma_dnaPipeTE_edta_lib d
                                  where genome_coverage = 1
                                  and sample_number = 5 
                                  and RM_t = 0.2
                                  and RM_class not in ('tRNA', 'rRNA')
                                  and sample_id != 'Ep_130_s_7'")

dprm_hit_prop <- dbGetQuery(conn, "select distinct sample_id, RM_class class, lib_hit_prop hit_prop                     
                                   from histoplasma_dnaPipeTE d
                                   where genome_coverage = 1
                                   and sample_number = 5 
                                   and RM_t = 0.2
                                   and RM_class not in ('tRNA', 'rRNA')
                                   and sample_id != 'Ep_130_s_7'")

e_full_per <- dbGetQuery(conn, "select class, 
                                 round(100.0 * sum(case when full_intact_te = 'Y' then 1 else 0 end) / count(*), 1) per_full                     
                                 from histoplasma_edta_tes_edta_only d, edta_te_classes_edta_only c
                                 where d.subclass = c.subclass
                                 and sample_id != 'Ep_130_s_7'
                                 group by class")

e_full_per_all <- dbGetQuery(conn, "select round(100.0 * sum(case when full_intact_te = 'Y' then 1 else 0 end) / count(*), 1) per_full                     
                                 from histoplasma_edta_tes_edta_only d, edta_te_classes_edta_only c
                                 where d.subclass = c.subclass
                                 and sample_id != 'Ep_130_s_7'")


classes <- sort(unique(c(erm_hit_prop$class, dpe_hit_prop$class, dprm_hit_prop$class)))


pdf(paste(plot_dir, "histoplasma_TE_hit_props.pdf", sep=''), height=12, width=12)
par(mfrow=c(2,2))
hist(erm_hit_prop$hit_prop, breaks=seq(0, 1.02, .01), col='black', xlab='TE hit proportion', ylab='TEs', main=c('EDTA-RM: All TEs', paste(round(100 * sum(erm_hit_prop$hit_prop >= 1) / nrow(erm_hit_prop), 1), '% full TEs', sep='')))
plot(1, type='n', bty='n', xaxt='n', yaxt='n', xlab='', ylab='', main=c('EDTA-only: All TEs', paste(e_full_per_all$per_full, '% full TEs', sep='')))
hist(dpe_hit_prop$hit_prop, breaks=seq(0, 1.02, .01), col='black', xlab='TE hit proportion', ylab='TEs', main=c('dnaPipeTE EDTA lib: All TEs', paste(round(100 * sum(dpe_hit_prop$hit_prop >= 1) / nrow(dpe_hit_prop), 1), '% full TEs', sep='')))
hist(dprm_hit_prop$hit_prop, breaks=seq(0, 1.02, .01), col='black', xlab='TE hit proportion', ylab='TEs', main=c('dnaPipeTE RM lib: All TEs', paste(round(100 * sum(dprm_hit_prop$hit_prop >= 1) / nrow(dprm_hit_prop), 1), '% full TEs', sep='')))
par(mfrow=c(1,1))


for (class in classes)
  {
  par(mfrow=c(2,2))
  hist(erm_hit_prop$hit_prop[erm_hit_prop$class == class], breaks=seq(0, 1.02, .01), col='black', xlab='TE hit proportion', ylab='TEs', main=c(paste('EDTA-RM:', class), paste(round(100 * sum(erm_hit_prop$hit_prop[erm_hit_prop$class == class] >= 1) / sum(erm_hit_prop$class == class), 1), '% full TEs', sep='')))
  plot(1, type='n', bty='n', xaxt='n', yaxt='n', xlab='', ylab='', main=c(paste('EDTA-only:', class), paste(max(c(e_full_per$per_full[e_full_per$class == class], 0), na.rm=T), '% full TEs', sep='')))
  hist(dpe_hit_prop$hit_prop[dpe_hit_prop$class == class], breaks=seq(0, 1.02, .01), col='black', xlab='TE hit proportion', ylab='TEs', main=c(paste('dnaPipeTE EDTA lib:', class), paste(round(100 * sum(dpe_hit_prop$hit_prop[dpe_hit_prop$class == class] >= 1) / sum(dpe_hit_prop$class == class), 1), '% full TEs', sep='')))
  hist(dprm_hit_prop$hit_prop[dprm_hit_prop$class == class], breaks=seq(0, 1.02, .01), col='black', xlab='TE hit proportion', ylab='TEs', main=c(paste('dnaPipeTE RM lib:', class), paste(round(100 * sum(dprm_hit_prop$hit_prop[dprm_hit_prop$class == class] >= 1) / sum(dprm_hit_prop$class == class), 1), '% full TEs', sep='')))
  par(mfrow=c(1,1))
  }
dev.off()
