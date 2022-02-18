rm(list = ls())

library(tximport)
library(magrittr)
library(GenomicFeatures)

cmdArgs <- commandArgs(trailingOnly = T)
gtf.path <- cmdArgs[1]
kallisto.dir <- cmdArgs[2]
out.dir <- cmdArgs[3]
samples <- cmdArgs[4 : length(cmdArgs)]
# gtf.path <- "/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/reindeer_appli/data/Homo_sapiens.GRCh38.99.gtf"
# kallisto.dir <- "/store/EQUIPES/SSFA/Data/CCLE/kallisto-ensembl99"

print(paste("Sample number: ", length(samples)))

abundance.files <- file.path(paste0(kallisto.dir), samples, "abundance.h5")
names(abundance.files) <- samples

txdb <- makeTxDbFromGFF(file = gtf.path)
tx2gene <- select(txdb, keys(txdb, keytype = "TXNAME"), "GENEID", "TXNAME")
    
gene.kallisto <- tximport(files = abundance.files, type = "kallisto", 
                          txIn = TRUE, txOut = FALSE, tx2gene = tx2gene, ignoreTxVersion = T)
gene.counts <- as.data.frame(gene.kallisto$counts)
write.table(gene.counts, file = paste0(out.dir, "/gene-counts-tximport.tsv"), 
            col.names = T, row.names = T, sep = "\t", quote = F)    
gene.abund <- as.data.frame(gene.kallisto$abundance)
write.table(gene.abund, file = paste0(out.dir, "/gene-abundance-tximport.tsv"), 
            col.names = T, row.names = T, sep = "\t", quote = F)    
