args<-commandArgs(T)
dir <- args[1]
setwd(dir)

summary_readcounts <- read.table("summary.txt",header=TRUE,stringsAsFactors = F)
summary_readcounts=summary_readcounts[grepl("ENS",rownames(summary_readcounts)),]
library(GenomicFeatures)
  ## Calculate gene length for human-92.
  txdb <- makeTxDbFromGFF("/home/zxchen/data/annotation/Ensembl_92/Homo_sapiens",format="gtf")
 # exons_gene <- exonsBy(txdb, by = "gene")
 # exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})
  #exons_gene_lens=unlist(exons_gene_lens)
  exons_tx <- exonsBy(txdb, by = "tx", use.names=TRUE)
  exons_tx_lens <- lapply(exons_tx,function(x){sum(width(reduce(x)))})
  exons_tx_lens=unlist(exons_tx_lens)
  total_reads=apply(summary_readcounts,2,sum)
  #exons_gene_len=exons_gene_lens[match(rownames(summary_readcounts),names(exons_gene_lens))]
  exons_tx_lens=exons_tx_lens[match(rownames(summary_readcounts),names(exons_tx_lens))]
  FPKM=(10^9*t(t(summary_readcounts)/total_reads))/(exons_tx_lens)
  write.table(FPKM,"FPKM.txt", col.names=NA,sep="\t",quote=F)