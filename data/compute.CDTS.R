compute.CDTS=function(data.path,bedtools,infile,infile.sorted,genome.file,indata)
{
anno.file2=paste(data.path,"CDTS_percentile.txt",sep="")
outfile2=gsub(".bed$",".tmp2.bed",infile)
cmd2=paste(bedtools,"intersect -wao -sorted -g", genome.file ,"-a",infile.sorted,"-b",anno.file2,"| awk '{print $5,$9}'>",outfile2)
system(cmd2)
anno.data2=fread(outfile2,header=F)
anno.data2$V2[anno.data2$V2=="."]=100
anno.data2$V2=as.numeric(anno.data2$V2)
CDTS_1st=tapply(anno.data2$V2,factor(anno.data2$V1,levels=1:nrow(indata)),function(x){
x[is.na(x)]=100
sum(x<1)*10
})/cnv.length[unique(anno.data2$V1)]
CDTS_5th=tapply(anno.data2$V2,factor(anno.data2$V1,levels=1:nrow(indata)),function(x){
x[is.na(x)]=100
sum(x<5)*10
})/cnv.length[unique(anno.data2$V1)]
CDTS_1st[CDTS_1st>1]=1;CDTS_5th[CDTS_5th>1]=1
tmp=file.remove(outfile2)
CDTS.score=cbind(CDTS_1st,CDTS_5th)
}

