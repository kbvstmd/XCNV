compute.haploinsufficiency=function(data.path,bedtools,infile,infile.sorted,genome.file,indata)
{
anno.file5=paste(data.path,"gencode_v19_features.bed",sep="")
outfile5=gsub(".bed$",".tmp5.bed",infile)
cmd5=paste(bedtools,"intersect -wao -sorted -g", genome.file, "-a",infile.sorted,"-b",anno.file5,">",outfile5)
system(cmd5)
anno.data5=fread(outfile5,header=F)
cols=c("pLI"  ,"Episcore"      ,"GHIS")
labels=factor(anno.data5$V5,levels=1:nrow(indata))
HI.score=tapply(1:nrow(anno.data5),labels,function(idx){tmp=as.matrix(anno.data5[idx,10:12]);tmp[tmp=="."]=0;apply(tmp,2,function(t){max(as.numeric(t))})})
HI.score=eval(as.call(c(rbind,HI.score)))
colnames(HI.score)=cols
tmp=file.remove(outfile5)
HI.score
}


