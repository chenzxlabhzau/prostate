args<-commandArgs(T)
dir <- args[1]
setwd(dir)
#files=list.files(pattern="SRR.*?.txt")
files=list.files(pattern="TCGA.*?.txt")
#feature=read.table(files[1],col.names=c("gene","count"),stringsAsFactors = F)
feature=read.table(files[1],col.names=c("transcript","count"),stringsAsFactors = F)
cdata=function(a){
  a=read.table(a,col.names=c("transcript","count"),stringsAsFactors = F)
  a=a[match(feature$transcript,a$transcript),]
  return(a$count)
}
reads_all=as.data.frame(sapply(files,cdata))
rownames(reads_all)=feature$transcript
names(reads_all)=stringr::str_sub(names(reads_all),end=-5)
write.table(reads_all,"summary.txt", col.names=NA,sep="\t",quote=F)



