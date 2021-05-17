args=commandArgs(T)
### Rscript CNV_unification.R 1000
resolution=as.numeric(args[1])
if(resolution<100)
{
quit(save="no")
}
load("all.cnvs.Rdata")
library(data.table)
all.cnv.list=lapply(unique(all.cnvs$chr),function(x){
  tmp=all.cnvs[all.cnvs$chr==x,]
  y=lapply(c("gain","loss"),function(t){
    tmp[tmp$type==t,]
  })
})

two_step_clustering=function(tmp,resolution=100)
  {
ClusterByPos=function(vec,min.dis,multi.samples=FALSE)
{
  names(vec)=1:length(vec)
  sort.vec=sort(vec)
  diff.sort.vec=diff(sort.vec)
  split.index=which(diff.sort.vec>min.dis)
  start.index=c(1,split.index+1)
  end.index=c(split.index,length(sort.vec))
  cls.list=lapply(apply(cbind(start.index,end.index),1,function(x){y=list(as.numeric(names(sort.vec)[x[1]:x[2]])) }),unlist)
  cls.list
}
#tmp=merged.cnv.list[[1]][[1]]
cls1=ClusterByPos(tmp$start,resolution)
cls1=cbind(rep(1:length(cls1),sapply(cls1,length)),unlist(cls1))
cls2=ClusterByPos(tmp$end,resolution)
cls2=cbind(rep(1:length(cls2),sapply(cls2,length)),unlist(cls2))
cls.step1=paste(cls1[order(cls1[,2]),1],cls2[order(cls2[,2]),1],sep="_")
tmp.list=lapply(unique(cls.step1),function(x){
  y=tmp[cls.step1==x,]
})
library(igraph)
## step2
seperation_by_largest_cliques=function(step1.mat,resolution)
{
  if(nrow(step1.mat)>1)
  {
    coord=as.matrix(step1.mat[,2:3])
    rownames(coord)=1:nrow(coord)
    pairs=which(as.matrix(dist(coord,method="maximum"))<=resolution,arr.ind=T)
    pairs=matrix(pairs[pairs[,1]<pairs[,2],],ncol=2)
    nodes=unique(c(pairs))
    g=NULL
  while(length(nodes)>1 & nrow(pairs)>0)
  {
    g=c(g,list(largest.cliques(graph_from_edgelist(pairs,directed=F))[[1]]))
    nodes=setdiff(nodes,unlist(g))
    pairs=matrix(pairs[!is.na(match(pairs[,1],nodes)) & !is.na(match(pairs[,2],nodes)),],ncol=2)
  }
  if(length(nodes)>0)
  {
    g=c(g,as.list(nodes))
  }
  if(length(setdiff(1:nrow(step1.mat),unlist(g)))>0)
  {
    g=c(g,as.list(setdiff(1:nrow(step1.mat),unlist(g))))
  }
  
  }else 
  {
    g=list(1)
  }

res=lapply(g,function(x){
  if(length(x)>1)
  {
    y=data.table(chr=unique(step1.mat[x,]$chr),
                 start=round(median(step1.mat[x,]$start)),
                 end=round(median(step1.mat[x,]$end)),
                 type=unique(step1.mat[x,]$type),
                 sample_id=paste(step1.mat[x,]$sample_id,collapse = ","),
                 database=paste(sort(unique(unlist(strsplit(step1.mat[x,]$database,",")))),collapse = ";"))
  }else{
    y=step1.mat[x,]
  }
y  
})
res=eval(as.call(c(rbind,res)))

}
res=lapply(tmp.list,function(x){seperation_by_largest_cliques(x,resolution)})
res=eval(as.call(c(rbind,res)))
}


merged.cnv.list1=lapply(all.cnv.list,function(x){
  lapply(x,two_step_clustering,resolution=resolution)
})
merged.cnv.res=eval(as.call(c(rbind,lapply(merged.cnv.list1,function(x){rbind(x[[1]],x[[2]])}))))
merged.cnv.res=merged.cnv.res[order(merged.cnv.res$chr,merged.cnv.res$start,merged.cnv.res$end,merged.cnv.res$type),]

