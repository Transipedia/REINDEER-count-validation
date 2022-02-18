rm(list = ls())

library(Biostrings)
library(stringr)
library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpointdensity)

tx2gene.path <- "../data/tx2gene_map.tsv" # ...SSFA/Data/CCLE/tx2gene_map.tsv
fos.path <- "../data/fos.txt" # ...SSFA/Data/CCLE/new-rmAdapt/reindeer-index-k31/fos.txt
gene.count.dir <- "../data/gene-counts/" # ...SSFA/MEMBERS/haoliang.xue/reindeer-appli/CCLE12-gene/gene-counts/
qres.dir <- "../data/reindeer-qres/" # ...SSFA/MEMBERS/haoliang.xue/reindeer-appli/CCLE12-gene/reindeer-qres/

out.dir <- "../results/two-versions-mean/"

tx2gene <- read.table(tx2gene.path, header = F, sep = "\t")
names(tx2gene) <- c("qname", "qgene")

### =========================== ###
### processing Kallisto results ###
### =========================== ###
kallisto.tab <- read.table(paste0(gene.count.dir, "/gene-counts-tximport.tsv"), header = T)
kallisto.tab$qgene <- rownames(kallisto.tab)
kallisto.tab.long <- pivot_longer(kallisto.tab, cols = -qgene, 
                                  names_to = "qsample", values_to = "kallisto")

p <- 40
k <- 31

### =========================== ###
### processing REINDEER results ###
### =========================== ###
qres.path <- paste0(qres.dir, "k", k, "/query_results/out_query_Reindeer_P", p,
                    "_krator-1000-random-genes_contigs-", k, "_0.out") # query result
qres <- read.table(qres.path, header = F)
colnames(qres) <- c("qname", read.table(fos.path, header = F)$V1 %>% as.character())

qres$qname <- str_extract(qres$qname, pattern = "ENST[0-9]+") # principal transcript as query name (1 gene = 1 tx)
qres.long <- pivot_longer(qres, names_to = "qsample", values_to = "qstr", cols = -qname) # make pivot table
qres.agg <- aggregate(qres.long$qstr, by = list(qres.long$qname, qres.long$qsample), FUN = paste)
names(qres.agg) <- c("qname", "qsample", "qstr")

print(nrow(qres.agg) / 12) # how many genes (principal transcrits) are considered

parseQuery <- function(query.res) {
    qsvect <- unlist(query.res["qstr"])
    count.vect.ignS <- NULL # skip 0 version
    count.vect.rep0 <- NULL # average on 0 version
    for (qs in qsvect) {
        qs.parsed <- str_split(qs, pattern = ",") %>% unlist()
        for (s in qs.parsed) {
            if (nchar(s) > 1) { # if the query is not a single "*"
                pos.info <- str_split(s, pattern = "-|:") %>% unlist()
                b <- as.integer(pos.info[1])
                e <- as.integer(pos.info[2])
                q <- pos.info[3]
                ## to skip stars ##
                if(q != "*") {
                    count.vect.ignS <- c(count.vect.ignS, rep(as.integer(q), e - b + 1))
                }
                ## to replace stars by 0 ##
                count.vect.rep0 <- c(count.vect.rep0, 
                                     rep(ifelse(q == "*", yes = 0, no = as.integer(q)), e - b + 1)) 
            }
        }
    }
    if (!is.null(count.vect.ignS)) { # if the query does not only have stars
        return(data.frame("qname" = as.character(query.res["qname"]), 
                          "qsample" = as.character(query.res["qsample"]), 
                          "qstr.parsed" = as.character(query.res["qstr.parsed"]), 
                          "qmean.ignS" = mean(count.vect.ignS), "qmean.rep0" = mean(count.vect.rep0)))
    } else {
        return(data.frame("qname" = as.character(query.res["qname"]), 
                          "qsample" = as.character(query.res["qsample"]), 
                          "qstr.parsed" = as.character(query.res["qstr.parsed"]),
                          "qmean.ignS" = 0, "qmean.rep0" = 0))
    }
}

qres.stats <- apply(qres.agg, MARGIN = 1, FUN = function(r) parseQuery(r)) %>%
    do.call(what = rbind) # calculate for each gene
qres.stats <- merge(x = qres.stats, y = tx2gene, by = "qname") # map principal transcripts to genes

### ============================== ###
### comparing REINDEER to Kallisto ###
### ============================== ###
cmp.tab <- merge(x = qres.stats, y = kallisto.tab.long, by = c("qgene", "qsample"))
cmp.tab <- cmp.tab[, c("qgene", "qmean.ignS", "qmean.rep0", "kallisto")]
cmp.tab <- cmp.tab[cmp.tab$kallisto >= 1E-5 & cmp.tab$qmean.rep0 >= 1E-5 & cmp.tab$qmean.ignS >= 1E-5, ]

# Pearson and Spearman correlations
cmp.corr.ps <- lapply(list("qmean.ignS", "qmean.rep0"), 
                      FUN = function(m)
                          return(data.frame("metrics" = m, 
                                            "cor.ps" = cor(cmp.tab[, m], cmp.tab[, "kallisto"], 
                                                           method = "pearson")))) %>%
    do.call(what = rbind)
cmp.corr.sp <- lapply(list("qmean.ignS", "qmean.rep0"), 
                      FUN = function(m) 
                          return(data.frame("metrics" = m, 
                                            "cor.sp" = cor(cmp.tab[, m], cmp.tab[, "kallisto"], 
                                                           method = "spearman")))) %>%
    do.call(what = rbind)
cmp.corr <- merge(cmp.corr.ps, cmp.corr.sp, by = "metrics")

# Make plot
cmp.tab.long <- pivot_longer(cmp.tab, names_to = "metrics", values_to = "value", cols = c(-qgene, -kallisto))
cmp.tab.long$value <- as.numeric(cmp.tab.long$value)
cmp.tab.long$metrics <- factor(cmp.tab.long$metrics,levels = c("qmean.ignS", "qmean.rep0"))

ggp <- ggplot() +
    geom_pointdensity(data = cmp.tab.long, aes(x = kallisto, y = value)) +
    geom_text(data = cmp.corr, size = 6, color = "black", hjust = 0, vjust = 1, 
              aes(x = 1E3, y = 1E-1, 
                  label = paste0("Pearson: ", round(cor.ps, 3), "\nSpearman: ", round(cor.sp, 3)))) +
    scale_color_viridis_c() +
    facet_wrap(metrics ~ ., ncol = 2) +
    scale_x_log10() +
    scale_y_log10() +
    xlab("kallisto-tximport est-counts") +
    ylab(paste("reindeer query (raw)")) +
    theme_bw() +
    theme(text = element_text(size = 20, family = "Arial"),
          legend.text = element_text(size = 15, family = "Arial"), legend.position = "none",
          plot.margin = margin(r = 2, unit = "cm"))
ggsave(paste0(out.dir, "CCLE-12-lung-k", k, "_P", p, ".svg"), plot = ggp, width = 12, height = 5)
