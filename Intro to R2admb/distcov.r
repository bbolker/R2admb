fitds=function(obs,width,detfct="hn",scale.formula=~1,exponent.formula=~1)
{
  # create scale design matrix with formula and data
  scale_mat=model.matrix(scale.formula,obs)
  if(detfct=="haz")exponent_mat=model.matrix(exponent.formula,obs)
  # write out data file
  con=file(paste(tplfile,".dat",sep=""),open="wt")
  writeLines(as.character(nrow(obs)),con)
  write(obs$distance,con,ncolumns=1)
  writeLines(as.character(width),con)
  if(detfct=="haz")
  {
     writeLines("2",con)
     writeLines("2",con)
     writeLines("2  2",con)
     writeLines(paste(ncol(scale_mat)," ",ncol(exponent_mat),sep=""),con)
     write(t(scale_mat),con,ncolumns=ncol(scale_mat))
     write(t(exponent_mat),con,ncolumns=ncol(exponent_mat))
  }else
  {
     writeLines("1",con)
     writeLines("1",con)
     writeLines("2",con)
     writeLines(paste(ncol(scale_mat),sep=""),con)
     write(t(scale_mat),con,ncolumns=ncol(scale_mat))
   }   
  close(con)
  run_admb(tplfile)
  results=read_admb(tplfile)
  cnames=paste("scale:",colnames(scale_mat),sep="")
  if(detfct=="haz")
    cnames=c(cnames,paste("exponent:",colnames(exponent_mat),sep=""))
  names(results$coefficients)=cnames
  rownames(results$vcov)=cnames
  colnames(results$vcov)=cnames
  return(results)
}

tplfile="distcov"
# compile tpl file
compile_admb(tplfile)
# get data from golf tee data in mrds
library(mrds)
data(book.tee.data)
obs=book.tee.data$book.tee.dataframe
obs=obs[obs$observer==1,]
obs=obs[obs$detected==1,]
# fit for different models
model1=fitds(obs,4,"haz",~1,~sex)
model1
model2=fitds(obs,4,"haz",~sex,~sex)
model2
model3=fitds(obs,4,"hn",~sex+size)
model3
model4=fitds(obs,4,"hn",~sex+exposure)
model4



