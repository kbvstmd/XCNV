library(data.table)
library(xgboost)
data.path=paste(script.path,"/data/",sep="")
bedtools=paste(script.path,"/tools/bedtools",sep="")
args=commandArgs(T)
if(length(args)==0)
{
message("\tPlease specify a BED file\n\tUsage:\tXCNV CNV.bed")
quit(save="no")
}else if(args=="--help" | args=="-h")
{
message("\tUsage:\tXCNV CNV.bed")
quit(save="no")
}else if(!file.exists(args[1]))
{
message(paste("\tNo such file:",args[1]))
message("\tUsage:\tXCNV CNV.bed")
quit(save="no")
}


infile=args[1]
indata=fread(infile,header=F)[,1:4]
n.cols=ncol(indata)
if(grepl("chr",indata$V1[1]))
{
indata=fread(infile,header=F)
indata$V1=gsub("chr","",indata$V1)
}
n.cnvs=nrow(indata)
if(n.cnvs==1)
{
indata=rbind(indata,indata)
indata$V3[1]=indata$V2[1]+49
}
indata$V2=as.integer(indata$V2)
indata$V3=as.integer(indata$V3)
check.data=function(mat)
{
test.chr=sum(is.na(match(mat$V1,c(1:22,"X","Y"))))
test.pos=sum(is.na(mat$V2))+sum(is.na(mat$V3))
test.type=sum(is.na(match(mat$V4,c("gain","loss"))))
test.width=sum(mat$V2>mat$V3)
if(test.chr>0 | test.pos>0 | test.type>0 | test.width)
{
message("Your file format is invalid!")
quit(save="no")
}
}

tmp=check.data(indata)

idx=match(indata$V1,c(1:22,"X","Y"))
indata=indata[order(idx,indata$V2,indata$V3),]
indata=cbind(indata,rows=1:nrow(indata))
infile.sorted=gsub(".bed$",".sort.bed",infile)
fwrite(indata,infile.sorted,sep='\t',row.names=F,col.names=F,quote=F)
labels=paste(indata$V1,indata$V2,indata$V3,indata$V4,sep="_")
cnv.length=indata$V3-indata$V2+1
genome.file=paste(data.path,"genome.txt",sep="")

timestamp()
message("### Start to compute X-CNV annotations")
timestamp()
source(paste(script.path,"/data/compute.ljb26.R",sep=""))
message("### Computing ljb26 annotations")
ljb26.scores=compute.ljb26(data.path,bedtools,infile,infile.sorted,genome.file,indata)
timestamp()
source(paste(script.path,"/data/compute.CDTS.R",sep=""))
message("### Computing CDTS annotations")
CDTS.scores=compute.CDTS(data.path,bedtools,infile,infile.sorted,genome.file,indata)
timestamp()
source(paste(script.path,"/data/compute.frequency.R",sep=""))
message("### Computing PAF annotations")
PAF=compute.frequency(data.path,bedtools,infile,infile.sorted,genome.file,indata)
timestamp()
source(paste(script.path,"/data/compute.ENCODE.R",sep=""))
message("### Computing ENCODE annotations")
ENCODE.scores=compute.ENCODE(data.path,bedtools,infile,infile.sorted,genome.file,indata)
timestamp()
source(paste(script.path,"/data/compute.HI.R",sep=""))
message("### Computing haploinsufficiency annotations")
HI.scores=compute.haploinsufficiency(data.path,bedtools,infile,infile.sorted,genome.file,indata)

#####
all.vars=cbind(indata[,1:4],ljb26.scores,CDTS.scores,PAF,ENCODE.scores,HI.scores)
output.file=gsub(".bed$",".output.csv",infile)
colnames(all.vars)[1]="#V1"
if(n.cnvs==1)
{
all.vars=all.vars[-1,]
}
all.vars$Type=ifelse(all.vars$V4=="gain",1,0)
all.vars$Length=all.vars$V3-all.vars$V2+1
load(paste(data.path,"xcnv.model.Rdata",sep=""))
features.in.model=xcnv.model$feature_names
all.vars.for.prediction=data.matrix(all.vars[,..features.in.model])
timestamp()
message("### Computing MVP score")
xcnv.score=predict(xcnv.model,newdata=all.vars.for.prediction)
all.vars.for.output=cbind(all.vars,MVP_score=xcnv.score)
colnames(all.vars.for.output)[1:4]=c("Chr","Start","End","Type")
fwrite(all.vars.for.output,output.file,row.names=F,quote=F)
message("### X-CNV runs successfully!")
