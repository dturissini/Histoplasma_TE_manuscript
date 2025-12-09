##library(phytools)
library(RSQLite)



base_dir <- "/work/users/d/t/dturissi/histoplasma/te_analysis"
plot_dir <- paste(base_dir, '/plots/', sep='')


setwd(base_dir)


conn <- dbConnect(dbDriver("SQLite"), 'histoplasma_tes.db')

subclasses <- dbGetQuery(conn, "select distinct subclass
                                from histoplasma_te_genome_pers_subclass
                                where genome_per > 0 
                                and te_approach != 'dnaPipeTE RM lib'
                                and sample_id != 'Ep_130_s_7'
                                order by subclass")$subclass
                                                         
classes <- dbGetQuery(conn, "select distinct class
                             from histoplasma_te_genome_pers_class a
                             where genome_per > 0 
                             and te_approach != 'dnaPipeTE RM lib'
                             and sample_id != 'Ep_130_s_7'
                             order by class")$class

plot_per_comp <- function(comp_table, comp_type, panel_name, class, te_type_a, te_type_b)
  { 
  comp <- dbGetQuery(conn, paste("select a.", comp_type, ", a.sample_id, a.genome_per per_a, b.genome_per per_b
                                  from ", comp_table, " a, ", comp_table, " b
                                  where a.sample_id = b.sample_id
                                  and a.", comp_type, " = '", class, "'
                                  and a.", comp_type, " = b.", comp_type, "
                                  and a.te_approach = '", te_type_a, "'
                                  and b.te_approach = '", te_type_b, "'
                                  and a.sample_id != 'Ep_130_s_7'", sep=''))
  
  if (max(comp$per_a) > 0 || max(comp$per_b) > 0)
    {
    max_per <- max(c(comp$per_a, comp$per_b))
    per_lm <- lm(comp$per_b ~ comp$per_a)
    per_lm_r2 <- signif(summary(per_lm)$r.squared, 3)
    
      
    if (nrow(coef(summary(per_lm))) == 2)
    {
    per_lm_p <- signif(coef(summary(per_lm))[2,4], 3)  
    
    if (! is.nan(per_lm_p) && per_lm_p < 0.01)
      {
      per_lm_p <- format(per_lm_p, scientific = TRUE)
      p_main <- gsub("e.*", "", per_lm_p)
      p_exponent <- gsub(".*e", "", per_lm_p)
      per_lm_p <- bquote(.(p_main) %*% 10^.(as.numeric(p_exponent)))
      }
    
    } else {
    per_lm_p <- NA  
    }
  
    plot(comp$per_a, comp$per_b, pch=20, xlim=c(0, max_per), ylim=c(0, max_per), xlab=paste(te_type_a, ' genomic %', sep=''), ylab=paste(te_type_b, ' genomic %', sep=''), main=bquote(atop(.(class) ~ ': ' ~ R^2 == ~ .(per_lm_r2) ~ ', p = ' ~ .(per_lm_p), .(te_type_a) ~ ' - ' ~ .(te_type_b))))
    mtext(panel_name, side=3, line=1.25, at= -.1 * max(comp$per_a), cex=2.2)
        
    if (! is.na(per_lm$coefficients[2]))
      {
      abline(per_lm$coefficients)
      cat(class, te_type_a, te_type_b, per_lm$coefficients, "\n")
      }   
    } else {
    plot(1, type='n', xlab=paste(te_type_a, ' genomic %', sep=''), ylab=paste(te_type_b, ' genomic %', sep=''), main=c(paste(class, ' not present in both', sep=''), paste(te_type_a, ' - ', te_type_b, sep='')))
    }
  }


pdf(paste(plot_dir, "histoplasma_te_genomic_per_comps_subclass.pdf", sep=''), height=15, width=10)
par(mfrow=c(3,2))
for (subclass in subclasses)
  {   
  plot_per_comp('histoplasma_te_genome_pers_subclass', 'subclass', '', subclass, 'EDTA-RM', 'EDTA-only')
  plot_per_comp('histoplasma_te_genome_pers_subclass', 'subclass', '', subclass, 'EDTA-RM', 'dnaPipeTE EDTA lib')
  plot_per_comp('histoplasma_te_genome_pers_subclass', 'subclass', '', subclass, 'EDTA-RM', 'dnaPipeTE RM lib')
  plot_per_comp('histoplasma_te_genome_pers_subclass', 'subclass', '', subclass, 'EDTA-only', 'dnaPipeTE EDTA lib')
  plot_per_comp('histoplasma_te_genome_pers_subclass', 'subclass', '', subclass, 'EDTA-only', 'dnaPipeTE RM lib')
  plot_per_comp('histoplasma_te_genome_pers_subclass', 'subclass', '', subclass, 'dnaPipeTE EDTA lib', 'dnaPipeTE RM lib')
  }
par(mfrow=c(1,1))
dev.off()



pdf(paste(plot_dir, "histoplasma_te_genomic_per_comps_class.pdf", sep=''), height=15, width=10)
par(mfrow=c(3,2))
for (class in classes)
  { 
  plot_per_comp('histoplasma_te_genome_pers_class', 'class', '', class, 'EDTA-RM', 'EDTA-only')
  plot_per_comp('histoplasma_te_genome_pers_class', 'class', '', class, 'EDTA-RM', 'dnaPipeTE EDTA lib')
  plot_per_comp('histoplasma_te_genome_pers_class', 'class', '', class, 'EDTA-RM', 'dnaPipeTE RM lib')
  plot_per_comp('histoplasma_te_genome_pers_class', 'class', '', class, 'EDTA-only', 'dnaPipeTE EDTA lib')
  plot_per_comp('histoplasma_te_genome_pers_class', 'class', '', class, 'EDTA-only', 'dnaPipeTE RM lib')
  plot_per_comp('histoplasma_te_genome_pers_class', 'class', '', class, 'dnaPipeTE EDTA lib', 'dnaPipeTE RM lib')
  }
par(mfrow=c(1,1))
dev.off()


pdf(paste(plot_dir, "manuscript_figure_4_histoplasma_te_genomic_per_comps_class.pdf", sep=''), height=12, width=15)
par(mfrow=c(3,4))
  plot_per_comp('histoplasma_te_genome_pers_class', 'class', 'A)', 'RNA TE', 'EDTA-only', 'EDTA-RM')
  plot_per_comp('histoplasma_te_genome_pers_class', 'class', 'B)', 'DNA TE', 'EDTA-only', 'EDTA-RM')
  
  plot_per_comp('histoplasma_te_genome_pers_class', 'class', 'C)', 'RNA TE', 'EDTA-only', 'dnaPipeTE EDTA lib')
  plot_per_comp('histoplasma_te_genome_pers_class', 'class', 'D)', 'DNA TE', 'EDTA-only', 'dnaPipeTE EDTA lib')

  plot_per_comp('histoplasma_te_genome_pers_class', 'class', 'E)', 'RNA TE', 'EDTA-only', 'dnaPipeTE RM lib')
  plot_per_comp('histoplasma_te_genome_pers_class', 'class', 'F)', 'DNA TE', 'EDTA-only', 'dnaPipeTE RM lib')

  plot_per_comp('histoplasma_te_genome_pers_class', 'class', 'G)', 'RNA TE', 'EDTA-RM', 'dnaPipeTE EDTA lib')
  plot_per_comp('histoplasma_te_genome_pers_class', 'class', 'H)', 'DNA TE', 'EDTA-RM', 'dnaPipeTE EDTA lib')

  plot_per_comp('histoplasma_te_genome_pers_class', 'class', 'I)', 'RNA TE', 'EDTA-RM', 'dnaPipeTE RM lib')
  plot_per_comp('histoplasma_te_genome_pers_class', 'class', 'J)', 'DNA TE', 'EDTA-RM', 'dnaPipeTE RM lib')

  plot_per_comp('histoplasma_te_genome_pers_class', 'class', 'K)', 'RNA TE', 'dnaPipeTE EDTA lib', 'dnaPipeTE RM lib')
  plot_per_comp('histoplasma_te_genome_pers_class', 'class', 'L)', 'DNA TE', 'dnaPipeTE EDTA lib', 'dnaPipeTE RM lib')
dev.off()


