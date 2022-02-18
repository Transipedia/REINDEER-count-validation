rm(list = ls())

library(tidyr)
library(magrittr)
library(ggplot2)
library(ggpointdensity)

tx2gene.path <- "../data/tx2gene_map.tsv" # ...SSFA/Data/CCLE/tx2gene_map.tsv
fos.path <- "../data/fos.txt" # ...SSFA/Data/CCLE/new-rmAdapt/reindeer-index-k31/fos.txt
gene.count.dir <- "../data/gene-counts/" # ...SSFA/MEMBERS/haoliang.xue/reindeer-appli/CCLE12-gene/gene-counts/

tx2gene <- read.table(tx2gene.path, header = F, sep = "\t")
names(tx2gene) <- c("qname", "qgene")

k <- 31
m <- "counts"

qstats.dir <- paste0("../results/interprets-k", k, "/")

kallisto.tab <- read.table(paste0(gene.count.dir, "/gene-", m, "-tximport.tsv"), header = T)
kallisto.tab$qgene <- rownames(kallisto.tab)
kallisto.tab.long <- pivot_longer(kallisto.tab, cols = -qgene, 
                                  names_to = "qsample", values_to = "kallisto")

stats.tab <- read.table(paste0(qstats.dir, "query_stats_P40.tsv"), header = T, sep = "\t")
stats.tab <- merge(x = stats.tab, y = tx2gene, by = "qname")

# if (g == "abundance") {
#     stats.tab <- merge(x = stats.tab, y = kc4norm.tab, by = "qsample")
#     stats.tab[, c("qmean", "qmedian", "qmode", "qmin.positive", "qmax", "qsum")] <-
#         stats.tab[, c("qmean", "qmedian", "qmode", "qmin.positive", "qmax", "qsum")] * 1E9 / stats.tab$real.kc
# }

cmp.tab <- merge(x = stats.tab, y = kallisto.tab.long, by = c("qgene", "qsample"))
cmp.tab <- cmp.tab[, c("qgene", "qmean", "qmedian", "qmode", "qmax", "qmin.positive", "qsum", "kallisto")]
cmp.tab <- cmp.tab[cmp.tab$kallisto >= 1E-5 & cmp.tab$qmedian >= 1E-5 & cmp.tab$qmode >= 1E-5, ]
cmp.corr.ps <- lapply(list("qmean", "qmedian", "qmode", "qmax", "qmin.positive", "qsum"), 
                      FUN = function(m) 
                          return(data.frame("metrics" = m, 
                                            "cor.ps" = cor(cmp.tab[, m], cmp.tab[, "kallisto"], 
                                                           method = "pearson")))) %>%
    do.call(what = rbind)
cmp.corr.sp <- lapply(list("qmean", "qmedian", "qmode", "qmax", "qmin.positive", "qsum"), 
                      FUN = function(m) 
                          return(data.frame("metrics" = m, 
                                            "cor.sp" = cor(cmp.tab[, m], cmp.tab[, "kallisto"], 
                                                           method = "spearman")))) %>%
    do.call(what = rbind)
cmp.corr <- merge(cmp.corr.ps, cmp.corr.sp, by = "metrics")
cmp.corr$ypos <- c(50000, 9000, 9000, no = 5000, 5000, 1E7)
cmp.tab.long <- pivot_longer(cmp.tab, names_to = "metrics", values_to = "value", 
                             cols = c(-qgene, -kallisto))
cmp.tab.long$value <- as.numeric(cmp.tab.long$value)

cmp.tab.long$metrics <- factor(cmp.tab.long$metrics, 
                               levels = c("qmean", "qmedian", "qmode", "qmin.positive", "qmax", "qsum"))

p <- ggplot() +
    geom_pointdensity(data = cmp.tab.long, aes(x = kallisto, y = value)) +
    geom_text(data = cmp.corr, size = 5, color = "black", hjust = 0, vjust = 1, 
              aes(x = min(cmp.tab.long$kallisto), y = ypos, 
                  label = paste0("Pearson: ", round(cor.ps, 3), 
                                 "\nSpearman: ", round(cor.sp, 3)))) +
    scale_color_viridis_c() +
    facet_wrap(metrics ~ ., nrow = 2, scales = "free_y") +
    scale_x_log10() +
    scale_y_log10() +
    xlab(ifelse(m == "abundance", yes = "kallisto-tximport tpm", no = "kallisto-tximport est-counts")) +
    ylab(paste("reindeer query", ifelse(m == "abundance", yes = " (scaled)", no = " (raw)"))) +
    theme_bw() +
    theme(text = element_text(size = 20, family = "Arial"),
          legend.text = element_text(size = 10, family = "Arial"),
          legend.title = element_text(size = 10, family = "Arial"),
          plot.margin = margin(r = 2, unit = "cm"))
ggsave(paste0(qstats.dir, "query-metrics-", m, ".svg"), plot = p, width = 12, height = 6)
