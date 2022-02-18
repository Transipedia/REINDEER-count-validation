rm(list = ls())

library(Biostrings)
library(stringr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(ggpointdensity)
library(MASS)

tx2gene.path <- "../data/tx2gene_map.tsv" # ...SSFA/Data/CCLE/tx2gene_map.tsv
fos.path <- "../data/fos.txt" # ...SSFA/Data/CCLE/new-rmAdapt/reindeer-index-k31/fos.txt
gene.count.dir <- "../data/gene-counts/" # ...SSFA/MEMBERS/haoliang.xue/reindeer-appli/CCLE12-gene/gene-counts/

tx2gene <- read.table(tx2gene.path, header = F, sep = "\t")
names(tx2gene) <- c("qname", "qgene")

k <- 31
g <- "counts"

qstats.dir <- paste0("../results/interprets-k", k, "/")

kallisto.tab <- read.table(paste0(gene.count.dir, "/gene-", g, "-tximport.tsv"), header = T)
kallisto.tab$qgene <- rownames(kallisto.tab)
kallisto.tab.long <- pivot_longer(kallisto.tab, cols = -qgene, names_to = "qsample", values_to = "kallisto")

for (m in c("qmean", "qsum")) {
    cor.df <- NULL
    to_plot.all <- NULL

    for (p in c(0, 40, 70, 100)) {
        print(paste(m, p))
        
        stats.tab <- read.table(paste0(qstats.dir, "query_stats_P", p, ".tsv"), header = T, sep = "\t")
        stats.tab <- stats.tab[, c("qname", "qsample", m)]
        stats.tab <- merge(x = stats.tab, y = tx2gene, by = "qname")
        
        cmp.tab <- merge(x = stats.tab, y = kallisto.tab.long, by = c("qgene", "qsample"))
        
        to_fit1 <- cmp.tab[cmp.tab[, "kallisto"] > 0 | cmp.tab[, m] > 0, ]
        neither_drop <- sum(to_fit1[, "kallisto"] > 0 & to_fit1[, m] > 0)
        reindeer_drop <- sum(to_fit1[, "kallisto"] > 0 & to_fit1[, m] == 0)
        kallisto_drop <- sum(to_fit1[, "kallisto"] == 0 & to_fit1[, m] > 0)
        drop.ratio <- round(reindeer_drop / (neither_drop + reindeer_drop + kallisto_drop) * 100, 3)
        
        to_fit2 <- cmp.tab[cmp.tab[, "kallisto"] > 1E-5 & cmp.tab[, m] > 1E-5, ]
        to_fit2$perc <- paste0("P=", p)
        cor.ps <- round(cor(to_fit2[, "kallisto"], to_fit2[, m], method = "pearson"), 3)
        cor.sp <- round(cor(to_fit2[, "kallisto"], to_fit2[, m], method = "spearman"), 3)
        to_plot.all <- rbind(to_plot.all, to_fit2[, c("perc", "qgene", "qsample", "kallisto", m)]) 
        
        cor.df <- rbind(cor.df, 
                        data.frame("perc" = paste0("P=", p),
                                   "cor.pearson" = cor.ps,
                                   "cor.spearman" = cor.sp,
                                   "drop.ratio" = paste0(drop.ratio, "%")))
    }
    write.table(cor.df[, c("perc", "cor.pearson", "cor.spearman", "drop.ratio")], 
                file = paste0(qstats.dir, "acrossP-corTable-k", k, "-", g, "-", m, ".tsv"), 
                col.names = T, row.names = F, quote = F, sep = "\t")

    to_plot.all$perc <- factor(to_plot.all$perc, levels = c("P=0", "P=40", "P=70", "P=100"))
    cor.df$perc <- factor(cor.df$perc, levels = c("P=0", "P=40", "P=70", "P=100"))
    names(to_plot.all)[names(to_plot.all) == m] <- "qvalue"
    
    plt <- ggplot() +
        geom_pointdensity(data = to_plot.all, aes(x = kallisto, y = qvalue)) +
        geom_text(data = cor.df, aes(x = 230000, y = ifelse(m == "qmean", yes = 2000, no = 1E6), 
                                     label = paste0("Pearson: ", cor.pearson, "\n",
                                                    "Spearman: ", cor.spearman, "\n",
                                                    "Reindeer dropped: ", drop.ratio, "\n")),
                  size = 5, color = "black", hjust = 1) +
        facet_wrap(. ~ perc, nrow = 2) +
        scale_color_viridis_c() +
        xlab(paste0("Kallisto-tximport ", ifelse(g == "counts", yes = "est-counts", no = "TPM"))) +
        ylab(paste0("Reindeer ", m, " ", ifelse(g == "counts", yes = "(raw)", no = "(normalized)"))) +
        theme_bw() +
        theme(text = element_text(size = 20, family = "Arial"),
              legend.title = element_text(size = 15, family = "Arial"),
              legend.text = element_text(size = 15, family = "Arial"))
    ggsave(paste0(qstats.dir, "acrossP-scatter-k", k, "-", g, "-", m, ".png"), plot = plt, 
           width = 12, height = 6)    
}
