##library(phytools)
library(RSQLite)
library(colorRamps)



base_dir <- "/work/users/d/t/dturissi/histoplasma/te_analysis"
plot_dir <- paste(base_dir, '/plots/', sep='')
csv_dir <- paste(base_dir, '/genome_per_csvs/', sep='')


setwd(base_dir)


conn <- dbConnect(dbDriver("SQLite"), 'histoplasma_tes.db')

#subclasses
rm_edta_subclass_pers <- dbGetQuery(conn, "select sample_short_name, subclass, round(sum(100.0 * (scaffold_end - scaffold_start + 1) / genome_len), 5) genome_per
                                           from histoplasma_edta_tes e, histoplasma_genome_lens l, histoplasma_tree_names t
                                           where e.sample_id = l.sample_id
                                           and e.sample_id = t.sample_id
                                           group by sample_short_name, subclass")


edta_only_subclass_pers <- dbGetQuery(conn, "select sample_short_name, subclass, round(sum(100.0 * (scaffold_end - scaffold_start + 1) / genome_len), 5) genome_per
                                             from histoplasma_edta_tes_edta_only e, histoplasma_genome_lens l, histoplasma_tree_names t
                                             where e.sample_id = l.sample_id
                                             and e.sample_id = t.sample_id
                                             group by sample_short_name, subclass")

dp_edta_lib_subclass_pers <- dbGetQuery(conn, "select sample_short_name, RM_subclass subclass, round(sum(100.0 * contig_len / genome_len), 5) genome_per
                                               from histoplasma_dnaPipeTE_edta_lib e, histoplasma_genome_lens l, histoplasma_tree_names t
                                               where e.sample_id = l.sample_id
                                               and e.sample_id = t.sample_id
                                               and genome_coverage == 1
                                               and sample_number == 5 
                                               and RM_t == 0.2
                                               group by sample_short_name, RM_subclass")


dp_rm_lib_subclass_pers <- dbGetQuery(conn, "select sample_short_name, RM_subclass subclass, round(sum(100.0 * contig_len / genome_len), 5) genome_per
                                             from histoplasma_dnaPipeTE e, histoplasma_genome_lens l, histoplasma_tree_names t
                                             where e.sample_id = l.sample_id
                                             and e.sample_id = t.sample_id
                                             and genome_coverage == 1
                                             and sample_number == 5 
                                             and RM_t == 0.2
                                             group by sample_short_name, RM_subclass")


#classes
rm_edta_class_pers <- dbGetQuery(conn, "select sample_short_name, class, round(sum(100.0 * (scaffold_end - scaffold_start + 1) / genome_len), 5) genome_per
                                        from histoplasma_edta_tes e, edta_te_classes c, histoplasma_genome_lens l, histoplasma_tree_names t
                                        where e.sample_id = l.sample_id
                                        and e.sample_id = t.sample_id
                                        and e.subclass = c.subclass
                                        group by sample_short_name, class")


edta_only_class_pers <- dbGetQuery(conn, "select sample_short_name, class, round(sum(100.0 * (scaffold_end - scaffold_start + 1) / genome_len), 5) genome_per
                                          from histoplasma_edta_tes_edta_only e, edta_te_classes_edta_only c, histoplasma_genome_lens l, histoplasma_tree_names t
                                          where e.sample_id = l.sample_id
                                          and e.sample_id = t.sample_id
                                          and e.subclass = c.subclass
                                          group by sample_short_name, class")

dp_edta_lib_class_pers <- dbGetQuery(conn, "select sample_short_name, RM_class class, round(sum(100.0 * contig_len / genome_len), 5) genome_per
                                            from histoplasma_dnaPipeTE_edta_lib e, histoplasma_genome_lens l, histoplasma_tree_names t
                                            where e.sample_id = l.sample_id
                                            and e.sample_id = t.sample_id
                                            and genome_coverage == 1
                                            and sample_number == 5 
                                            and RM_t == 0.2
                                            group by sample_short_name, RM_class")


dp_rm_lib_class_pers <- dbGetQuery(conn, "select sample_short_name, RM_class class, round(sum(100.0 * contig_len / genome_len), 5) genome_per
                                          from histoplasma_dnaPipeTE e, histoplasma_genome_lens l, histoplasma_tree_names t
                                          where e.sample_id = l.sample_id
                                          and e.sample_id = t.sample_id
                                          and genome_coverage == 1
                                          and sample_number == 5 
                                          and RM_t == 0.2
                                          group by sample_short_name, RM_class")


samples <- dbGetQuery(conn, "select species, t.sample_short_name
                             from histoplasma_tree_names t, histoplasma_manuscript_samples m
                             where t.tree_name = m.tree_name
                             order by case when species = 'B. parva' then 'AAAAA' else species end desc, t.sample_short_name desc")

species <- dbGetQuery(conn, "select distinct species
                             from histoplasma_tree_names t, histoplasma_manuscript_samples m
                             where t.tree_name = m.tree_name
                             order by case when species = 'B. parva' then 'AAAAA' else species end")$species


subclasses <- sort(unique(c(rm_edta_subclass_pers$subclass, edta_only_subclass_pers$subclass, dp_edta_lib_subclass_pers$subclass)))
dp_rm_lib_subclasses <- sort(unique(dp_rm_lib_subclass_pers$subclass))
classes <- sort(unique(c(rm_edta_class_pers$class, edta_only_class_pers$class, dp_edta_lib_class_pers$class, dp_rm_lib_class_pers$class)))


subclass_max_per <- max(c(rm_edta_subclass_pers$genome_per, edta_only_subclass_pers$genome_per, dp_edta_lib_subclass_pers$genome_per, dp_rm_lib_subclass_pers$genome_per))
subclass_min_per <- min(c(rm_edta_subclass_pers$genome_per, edta_only_subclass_pers$genome_per, dp_edta_lib_subclass_pers$genome_per, dp_rm_lib_subclass_pers$genome_per))

class_max_per <- max(c(rm_edta_class_pers$genome_per, edta_only_class_pers$genome_per, dp_edta_lib_class_pers$genome_per, dp_rm_lib_class_pers$genome_per))
class_min_per <- min(c(rm_edta_class_pers$genome_per, edta_only_class_pers$genome_per, dp_edta_lib_class_pers$genome_per, dp_rm_lib_class_pers$genome_per))



te_per_heatmap <- function(te_type, class_type, genome_pers, classes, samples, species, min_per, max_per, num_cex, bool_write_csv, bool_log_col)
  {
  if (bool_write_csv)
    {genomic_prop_df <- data.frame()}
  
  species_col <- c("black", "deeppink1", "grey45", "chocolate4", "darkgoldenrod1", "cyan3", "blue", "darkorange2", "purple", "red")
  per_cols <- adjustcolor(matlab.like(100), 0.4)


  plot(1, type='n', xlab='', ylab='', xaxt='n', yaxt='n', bty='n', main=te_type, xlim=c(0, length(classes)), ylim=c(0, nrow(samples)))
  for (sample_i in 1:nrow(samples))
    {axis(2, at=sample_i - 0.5, labels=samples$sample_short_name[sample_i], col.axis=species_col[which(species == samples$species[sample_i])], las=2, cex.axis=0.5)} 
  
  if (class_type == 'subclass')
    {
    if (te_type == 'dnaPipeTE RM lib')
      {
      par(xpd=T)
      rect(-10.75, nrow(samples) + 13, -1.25, nrow(samples) + 2)
      for (species_i in 1:length(species))
        {
        segments(-10.75, min(which(samples$species == species[species_i])) - 0.75, -10.75, max(which(samples$species == species[species_i])) - 0.25, col=species_col[species_i])
        if (species[species_i] %in% c("B. parva", "H. capsulatum", "H. mississippiense", "H. ohiense", "H. suramericanum"))
          {
          text(-9.5, nrow(samples) + 13 - 1 * species_i, bquote(italic(.(species[species_i]))), cex=0.7, adj=0, col=species_col[species_i])
          } else {
          text(-9.5, nrow(samples) + 13 - 1 * species_i, species[species_i], cex=0.7, adj=0, col=species_col[species_i])
          }
        }
      par(xpd=F)  
      } else {
      par(xpd=T)
      rect(-2.75, nrow(samples) + 13, -0.25, nrow(samples) + 2)
      for (species_i in 1:length(species))
        {
        segments(-2.75, min(which(samples$species == species[species_i])) - 0.75, -2.75, max(which(samples$species == species[species_i])) - 0.25, col=species_col[species_i])
        if (species[species_i] %in% c("B. parva", "H. capsulatum", "H. mississippiense", "H. ohiense", "H. suramericanum"))
          {
          text(-2.5, nrow(samples) + 13 - 1 * species_i, bquote(italic(.(species[species_i]))), cex=0.7, adj=0, col=species_col[species_i])
          } else {
          text(-2.5, nrow(samples) + 13 - 1 * species_i, species[species_i], cex=0.7, adj=0, col=species_col[species_i])
          }
        }
      par(xpd=F)
      }
    }

  axis(3, at=1:length(classes) - 0.5, labels=classes, las=2, cex.axis=0.5)

  for (sample_i in 1:nrow(samples))
    {
    if (bool_write_csv)
      {csv_line <- c(samples$sample_short_name[sample_i])}
    
    for (class_i in 1:length(classes))
      {
      genome_per_raw <- genome_pers$genome_per[genome_pers$sample_short_name == samples$sample_short_name[sample_i] & genome_pers[, 2] == classes[class_i]]

      
      if (length(genome_per_raw) == 0)
        {
        genome_per_raw <- 0
        genome_per <- 0
        per_col <- adjustcolor('darkgrey', 0.4)
        } else {
        if (bool_log_col)
          {
          per_col_i <- max(round(100 * (log(genome_per_raw, 10) - log(min_per, 10)) / (log(max_per, 10) - log(min_per, 10))), 1)
          per_col <- per_cols[per_col_i]
          } else {
          per_col <- per_cols[max(round(100 * genome_per_raw / max_per), 1)]  
          }
        
        if (genome_per_raw > 1)
          { 
          genome_per <- round(genome_per_raw, 1)  
          } else {
          genome_per <- gsub("e\\-0", "e\\-", formatC(genome_per_raw, format = "e", digits = 1))
          }
        }
          
      rect(class_i - 1, sample_i - 1, class_i, sample_i, col=per_col, border=NA)
      text(class_i - 0.5, sample_i - 0.5, genome_per, cex=num_cex)      
      
      if (bool_write_csv)
        {csv_line <- c(csv_line, genome_per_raw)}
      }
      
    if (bool_write_csv)
      {genomic_prop_df <- rbind(genomic_prop_df, csv_line)}
    } 
  
  if (bool_write_csv)
    {
    csv_file <- paste(csv_dir, 'histoplasma_te_', class_type, '_genomic_prop_', te_type, '.csv', sep='')
    write.table(genomic_prop_df, csv_file, sep=',', row.names=F, col.names=c('sample_id', classes))
    }
  }



#subclasses
pdf(paste(plot_dir, "histoplasma_te_subclass_genomic_per_heatmaps.pdf", sep=''), height=12, width=12)
par(mar=c(1.1, 5.1, 9.1, 1.1))
te_per_heatmap('EDTA-RM', 'subclass', rm_edta_subclass_pers, subclasses, samples, species, subclass_min_per, subclass_max_per, 0.6, T, F)
te_per_heatmap('EDTA-only', 'subclass', edta_only_subclass_pers, subclasses, samples, species, subclass_min_per, subclass_max_per, 0.6, T, F)
te_per_heatmap('dnaPipeTE EDTA lib', 'subclass', dp_edta_lib_subclass_pers, subclasses, samples, species, subclass_min_per, subclass_max_per, 0.6, T, F)
te_per_heatmap('dnaPipeTE RM lib', 'subclass', dp_rm_lib_subclass_pers, dp_rm_lib_subclasses, samples, species, subclass_min_per, subclass_max_per, 0.2, T, F)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

pdf(paste(plot_dir, "histoplasma_te_subclass_genomic_per_heatmaps_log_color.pdf", sep=''), height=12, width=12)
par(mar=c(1.1, 5.1, 9.1, 1.1))
te_per_heatmap('EDTA-RM', 'subclass', rm_edta_subclass_pers, subclasses, samples, species, subclass_min_per, subclass_max_per, 0.6, F, T)
te_per_heatmap('EDTA-only', 'subclass', edta_only_subclass_pers, subclasses, samples, species, subclass_min_per, subclass_max_per, 0.6, F, T)
te_per_heatmap('dnaPipeTE EDTA lib', 'subclass', dp_edta_lib_subclass_pers, subclasses, samples, species, subclass_min_per, subclass_max_per, 0.6, F, T)
te_per_heatmap('dnaPipeTE RM lib', 'subclass', dp_rm_lib_subclass_pers, dp_rm_lib_subclasses, samples, species, subclass_min_per, subclass_max_per, 0.2, F, T)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()


#classes
pdf(paste(plot_dir, "histoplasma_te_class_genomic_per_heatmaps.pdf", sep=''), height=12, width=12)
par(mar=c(1.1, 5.1, 9.1, 1.1))
te_per_heatmap('EDTA-RM', 'class', rm_edta_class_pers, classes, samples, species, class_min_per, class_max_per, 0.6, T, F)
te_per_heatmap('EDTA-only', 'class', edta_only_class_pers, classes, samples, species, class_min_per, class_max_per, 0.6, T, F)
te_per_heatmap('dnaPipeTE EDTA lib', 'class', dp_edta_lib_class_pers, classes, samples, species, class_min_per, class_max_per, 0.6, T, F)
te_per_heatmap('dnaPipeTE RM lib', 'class', dp_rm_lib_class_pers, classes, samples, species, class_min_per, class_max_per, 0.6, T, F)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

pdf(paste(plot_dir, "histoplasma_te_class_genomic_per_heatmaps_log_color.pdf", sep=''), height=12, width=12)
par(mar=c(1.1, 5.1, 9.1, 1.1))
te_per_heatmap('EDTA-RM', 'class', rm_edta_class_pers, classes, samples, species, class_min_per, class_max_per, 0.6, F, T)
te_per_heatmap('EDTA-only', 'class', edta_only_class_pers, classes, samples, species, class_min_per, class_max_per, 0.6, F, T)
te_per_heatmap('dnaPipeTE EDTA lib', 'class', dp_edta_lib_class_pers, classes, samples, species, class_min_per, class_max_per, 0.6, F, T)
te_per_heatmap('dnaPipeTE RM lib', 'class', dp_rm_lib_class_pers, classes, samples, species, class_min_per, class_max_per, 0.6, F, T)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()


#manuscripts
pdf(paste(plot_dir, "manuscript_figure_5_genomic_per_heatmaps_log_color_EDTA_RM.pdf", sep=''), height=12, width=12)
par(mar=c(1.1, 5.1, 9.1, 1.1))
te_per_heatmap('EDTA-RM', 'subclass', rm_edta_subclass_pers, subclasses, samples, species, subclass_min_per, subclass_max_per, 0.6, F, T)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()



pdf(paste(plot_dir, "manuscript_figure_S1_genomic_per_heatmaps_log_color_EDTA_only.pdf", sep=''), height=12, width=12)
par(mar=c(1.1, 5.1, 9.1, 1.1))
te_per_heatmap('EDTA-only', 'subclass', edta_only_subclass_pers, subclasses, samples, species, subclass_min_per, subclass_max_per, 0.6, F, T)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()


pdf(paste(plot_dir, "manuscript_figure_S2_genomic_per_heatmaps_log_color_dnaPipeTE_RM_lib.pdf", sep=''), height=12, width=12)
par(mar=c(1.1, 5.1, 9.1, 1.1))
te_per_heatmap('dnaPipeTE RM lib', 'subclass', dp_rm_lib_subclass_pers, dp_rm_lib_subclasses, samples, species, subclass_min_per, subclass_max_per, 0.2, F, T)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()


pdf(paste(plot_dir, "manuscript_figure_S3_genomic_per_heatmaps_log_color_dnaPipeTE_EDTA_lib.pdf", sep=''), height=12, width=12)
par(mar=c(1.1, 5.1, 9.1, 1.1))
te_per_heatmap('dnaPipeTE EDTA lib', 'subclass', dp_edta_lib_subclass_pers, subclasses, samples, species, subclass_min_per, subclass_max_per, 0.6, F, T)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()


