rm(list = ls())

library(Biostrings)
library(stringr)
library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)

fos.path <- "../data/fos.txt" # ...SSFA/Data/CCLE/new-rmAdapt/reindeer-index-k31/fos.txt
qres.dir <- "../data/reindeer-qres/" # ...SSFA/MEMBERS/haoliang.xue/reindeer-appli/CCLE12-gene/reindeer-qres/

smps <- read.table(fos.path, header = F)$V1 %>% as.character()

k <- 31

out.dir <- paste0("../results/interprets-k", k, "/")

for (p in c(0, 40, 70, 100)) {
    print(p)
    
    qres.path <- paste0(qres.dir, "k", k, "/query_results/out_query_Reindeer_P", p,
                        "_krator-1000-random-genes_contigs-", k, "_0.out") # query result
    qres <- read.table(qres.path, header = F)
    colnames(qres) <- c("qname", smps)
    
    qres$qname <- str_extract(qres$qname, pattern = "ENST[0-9]+") # principal transcript as query name (1 gene = 1 tx)
    qres.long <- pivot_longer(qres, names_to = "qsample", values_to = "qstr", cols = -qname) # make pivot table
    qres.agg <- aggregate(qres.long$qstr, by = list(qres.long$qname, qres.long$qsample), FUN = paste)
    names(qres.agg) <- c("qname", "qsample", "qstr")
    qres.agg$qstr.parsed <- lapply(qres.agg$qstr, FUN = function(qs) paste0(qs, collapse = ";")) %>%
        unlist()
    
    print(nrow(qres.agg) / 12)
    
    parseQuery <- function(query.res) {
        qsvect <- unlist(query.res["qstr"])
        count.vect <- NULL
        for (qs in qsvect) {
            qs.parsed <- str_split(qs, pattern = ",") %>% unlist()
            for (s in qs.parsed) {
                if (nchar(s) > 1) { # if the query is not a single "*"
                    pos.info <- str_split(s, pattern = "-|:") %>% unlist()
                    b <- as.integer(pos.info[1])
                    e <- as.integer(pos.info[2])
                    q <- pos.info[3]
                    count.vect <- c(count.vect, 
                                    rep(ifelse(q != "*", yes = as.integer(q), no = 0), e - b + 1)) 
                }
            }
        }
        if (!is.null(count.vect)) {
            return(data.frame("qname" = as.character(query.res["qname"]), "qsample" = as.character(query.res["qsample"]), 
                              "qstr.parsed" = as.character(query.res["qstr.parsed"]), 
                              "qmean" = mean(count.vect), "qmedian" = median(count.vect), 
                              "qmode" = as.numeric(names(which.max(table(count.vect)))),
                              "qmax" = max(count.vect), 
                              "qmin.positive" = ifelse(any(count.vect > 0), yes = min(count.vect[count.vect > 0]), no = 0),
                              "qsum" = sum(count.vect)))
        } else {
            return(data.frame("qname" = as.character(query.res["qname"]), "qsample" = as.character(query.res["qsample"]), 
                              "qstr.parsed" = as.character(query.res["qstr.parsed"]),
                              "qmean" = 0, "qmedian" = 0, 
                              "qmode" = 0,
                              "qmax" = 0, "qmin.positive" = 0, 
                              "qsum" = 0))
        }
    }
    
    qres.stats <- apply(qres.agg, MARGIN = 1, FUN = function(r) parseQuery(r)) %>%
        do.call(what = rbind) # calculate for each gene
    write.table(qres.stats[, c("qname", "qsample", "qstr.parsed", "qmean", "qmedian", "qmode", 
                               "qmax", "qmin.positive", "qsum")], 
                col.names = T, row.names = F, quote = F, sep = "\t",
                file = paste0(out.dir, "/query_stats", "_P", p, ".tsv"))
}
