rm (list = ls())

library(stringr)
library(tidyr)
library(Biostrings)
library(ggplot2)
library(ggpointdensity)

k <- 31
g <- "counts"

tx2gene.path <- "../data/tx2gene_map.tsv" # ...SSFA/Data/CCLE/tx2gene_map.tsv
fos.path <- "../data/fos.txt" # ...SSFA/Data/CCLE/new-rmAdapt/reindeer-index-k31/fos.txt

tx2gene <- read.table(tx2gene.path, header = F, sep = "\t")
names(tx2gene) <- c("qname", "qgene")

ctg.fa.path <- paste0("../data/krator-1000-random-genes_contigs-", k, ".fa") # ...SSFA/Data/CCLE/kmerator-spec-ctgs/krator-1000-random-genes_contigs-31.fa
kamrat.path <- paste0("../data/kamrat-queryRes-k", k, ".tsv") # ...SSFA/MEMBERS/haoliang.xue/reindeer-appli/CCLE12-gene/kamrat-res/kamrat-queryRes-k31.tsv
qres.dir <- "../data/reindeer-qres/" # ...SSFA/MEMBERS/haoliang.xue/reindeer-appli/CCLE12-gene/reindeer-qres/

ctg.fa <- readDNAStringSet(ctg.fa.path)
name2seq <- data.frame("qname" = names(ctg.fa), "contig" = as.character(ctg.fa))

kamrat.res <- read.table(kamrat.path, header = T)
colnames(kamrat.res)[1] <- "contig"
kamrat.res <- merge(kamrat.res, name2seq)[, -1]
kamrat.res.long <- pivot_longer(kamrat.res, cols = -qname, names_to = "file", values_to = "count")
kamrat.res.long$file <- str_remove(kamrat.res.long$file, pattern = "_[1-2]?.trim.fastq.gz")
kamrat.res.agg <- aggregate(kamrat.res.long$count, 
                            by = list(kamrat.res.long$qname, kamrat.res.long$file), FUN = sum)
names(kamrat.res.agg) <- c("qname", "qsample", "count.kamrat")

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
                          "qmean" = mean(count.vect.ignS)))
    } else {
        return(data.frame("qname" = as.character(query.res["qname"]), 
                          "qsample" = as.character(query.res["qsample"]), 
                          "qmean" = 0))
    }
}

cmp.tab <- NULL
for (p in c(0, 40, 70, 100)) {
    print(p)
    
    qres.path <- paste0(qres.dir, "k", k, "/query_results/out_query_Reindeer_P", p,
                        "_krator-1000-random-genes_contigs-", k, "_0.out") # query result
    qres <- read.table(qres.path, header = F)
    colnames(qres) <- c("qname", read.table(fos.path, header = F)$V1 %>% as.character())
    
    qres$qname <- str_remove(qres$qname, pattern = ">")
    qres.long <- pivot_longer(qres, names_to = "qsample", values_to = "qstr", cols = -qname) # make pivot table
    qres.agg <- aggregate(qres.long$qstr, by = list(qres.long$qname, qres.long$qsample), FUN = paste)
    names(qres.agg) <- c("qname", "qsample", "qstr")
    
    print(nrow(qres.agg) / 12) # how many genes (principal transcrits) are considered
    
    qres.stats <- apply(qres.agg, MARGIN = 1, FUN = function(r) parseQuery(r)) %>%
        do.call(what = rbind) # calculate for each gene
    qres.stats$perc <- paste0("P=", p)
    
    cmp.tab <- rbind(cmp.tab, merge(kamrat.res.agg, qres.stats))
}

cmp.tab$perc <- factor(cmp.tab$perc, levels = c("P=0", "P=40", "P=70", "P=100"))
plt <- ggplot(cmp.tab) +
    geom_pointdensity(aes(x = count.kamrat, y = qmean)) +
    scale_color_viridis_c() +
    facet_wrap(perc ~ ., nrow = 1) +
    xlab("countTags + KaMRaTquery mean") +
    ylab("Reindeer query mean") +
    theme_bw() +
    theme(text = element_text(size = 20, family = "Arial"),
          legend.title = element_text(size = 15, family = "Arial"),
          legend.text = element_text(size = 15, family = "Arial"))
ggsave(plt, filename = paste0("../results/cmp2countTags/qmean_kamratquerymean.png"),
       width = 16, height = 4)
