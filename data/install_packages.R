pkgs=installed.packages()
required.packages=c("xgboost","data.table")
required.packages=setdiff(required.packages,rownames(pkgs))
if(length(required.packages)>0)
{
chooseCRANmirror(ind=1)
install.packages(required.packages)
}

