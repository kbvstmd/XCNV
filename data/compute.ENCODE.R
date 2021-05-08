compute.ENCODE=function(data.path,bedtools,infile,infile.sorted,genome.file,indata)
{
anno.file4=paste(data.path,"hg19-ccREs.bed",sep="")
outfile4=gsub(".bed$",".tmp4.bed",infile)
cmd4=paste(bedtools,"intersect -wao -sorted -g", genome.file, "-a",infile.sorted,"-b",anno.file4,">",outfile4)
system(cmd4)
anno.data4=fread(outfile4,header=F)
#labels4=factor(paste(anno.data4$V1,anno.data4$V2,anno.data4$V3,anno.data4$V4,sep="_"),levels=labels)
cols=c("pELS","CTCF-bound","PLS","dELS", "CTCF-only","DNase-H3K4me3")
feature.mat=t(sapply(strsplit(as.matrix(anno.data4$V9),","),function(x){y=rep(0,6);x=x[x!="."];if(length(x)>0){y[match(x,cols)]=1};y}))
labels=factor(anno.data4$V5,levels=1:nrow(indata))
len=as.matrix(anno.data4$V10)
len[len=="."]=0
feature.score=apply(feature.mat,2,function(x){x*len})
feature.sum.score=tapply(1:length(labels),labels,function(idx){y=apply(matrix(feature.score[idx,],ncol=6),2,sum)})
feature.sum.score=eval(as.call(c(rbind,feature.sum.score)))
len1=indata$V3-indata$V2+1
feature.normalized.score=apply(feature.sum.score,2,function(x){x/len1})
colnames(feature.normalized.score)=cols
tmp=file.remove(outfile4)
feature.normalized.score
}

