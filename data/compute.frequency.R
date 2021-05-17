#### frequency start
compute.frequency=function(data.path,bedtools,infile,infile.sorted,genome.file,indata)
{
anno.file3=paste(data.path,"merged.cnv.sites.bed",sep="")
outfile3=gsub(".bed$",".tmp3.bed",infile)
cmd3=paste(bedtools, "intersect -wao -f 0.5 -F 0.5 -sorted -g", genome.file ,"-a",infile.sorted,"-b",anno.file3,"| awk '{print $5,$9,$10}' >",outfile3)
system(cmd3)
anno.data3=fread(outfile3,header=F)
anno.data3$V3[anno.data3$V3<0]=NA
sample.info=fread(paste(data.path,"/sample.info.txt",sep=""))
n.samples=table(sample.info$ethnicity_abbr)
merged.content=fread(paste(data.path,"merged.cnv.sample.info.txt",sep=""))
merged.content=merged.content[!is.na(match(merged.content$V1,anno.data3$V3)),]

merged.content.list1=tapply(merged.content$V2,merged.content$V1,unique)
merged.content.list2=merged.content.list1[as.character(anno.data3$V3)]
V1=rep(anno.data3$V1,sapply(merged.content.list2,length))
V2=rep(anno.data3$V2,sapply(merged.content.list2,length))
V3=eval(as.call(c(c,merged.content.list2)))
V4=merged.content$V3[match(V3,merged.content$V2)]
gain.samples.count=table(factor(V1,levels=1:nrow(indata))[V2=="gain"],
factor(V4[V2=="gain"],levels=names(n.samples)))
loss.samples.count=table(factor(V1,levels=1:nrow(indata))[V2=="loss"],
factor(V4[V2=="loss"],levels=names(n.samples)))

gain.samples=sapply(colnames(gain.samples.count),function(x){gain.samples.count[,x]/n.samples[x]})
loss.samples=sapply(colnames(loss.samples.count),function(x){loss.samples.count[,x]/n.samples[x]})
colnames(gain.samples)=paste("gain_freq_",colnames(gain.samples),sep="")
colnames(loss.samples)=paste("loss_freq_",colnames(loss.samples),sep="")
gain.samples.overall=table(factor(V1,levels=1:nrow(indata))[V2=="gain"])
loss.samples.overall=table(factor(V1,levels=1:nrow(indata))[V2=="loss"])

gain.freq.overall=gain.samples.overall/nrow(sample.info)
loss.freq.overall=loss.samples.overall/nrow(sample.info)
tmp=file.remove(outfile3)
all.freq=cbind(gain.samples,loss.samples,gain.freq=gain.freq.overall,loss.freq=loss.freq.overall)
all.freq
}

