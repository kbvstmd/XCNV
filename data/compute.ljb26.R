compute.ljb26=function(data.path,bedtools,infile,infile.sorted,genome.file,indata)
{
anno.file1.sites=paste(data.path,"hg19_ljb26_all_converted_sites.vcf",sep="")
anno.file1.score=paste(data.path,"hg19_ljb26_all_converted_scores.txt",sep="")
outfile1=gsub(".bed$",".tmp1.bed",infile)
cmd1=paste(bedtools,"intersect -wao -sorted -g",genome.file ,"-a",infile.sorted,"-b",anno.file1.sites,"| awk '{print $5,$9}' >",outfile1)
system(cmd1)
anno.data1=fread(outfile1,header=F)
anno.data1$V2[anno.data1$V2=="."]=NA
anno.data1$V2=as.numeric(anno.data1$V2)
anno.score1=fread(anno.file1.score,header=T)
scores1=sapply(colnames(anno.score1)[1:9],function(x){
tapply(anno.score1[[x]][anno.data1$V2],factor(anno.data1$V1,levels=1:nrow(indata)),function(t){
length(which(t>0))/length(t)
})
})
missing.values=c(0,0,-12.30,-11.958,-20.000,0.0003)
names(missing.values)=colnames(anno.score1)[10:15]
scores2=sapply(colnames(anno.score1)[10:15],function(x){
tapply(anno.score1[[x]][anno.data1$V2],factor(anno.data1$V1,levels=1:nrow(indata)),function(t){
t[is.na(t)]=missing.values[x]
t=as.numeric(t)
mean(t)
})
})
file.remove(outfile1)
ljb26.scores=cbind(scores1,scores2)
}




