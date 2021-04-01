library("kimisc")
library("magic")
library("abind")
library("grid")
library("tiff")

# path: input folder of tiff images
# output_path: output folder to store prediction results
# max_n: split by this number of images
args=commandArgs(trailingOnly = TRUE)
path=args[1]
output_path=args[2]
max_n=as.numeric(args[3])
#deeplearning_R=args[4]

#deeplearning_R=thisfile() # path to the executing R script
deeplearning_R=gsub("/deeplearning.R","",thisfile())
#setwd(gsub("deeplearning.R","",deeplearning_R))
setwd(deeplearning_R)

if (!file.exists(path)) {stop("Input folder doesn't exist!")}
if (!file.exists(output_path)) {stop("Output folder doesn't exist!")}
cat(paste("Working on Tiff files in",path,"\n"))

files=list.files(path,pattern=".tif",full=T)
n_total=length(files)
if (n_total==0) {stop("No image files!")}

for (j in 1:ceiling(n_total/max_n))
{
  cat(paste("Working on partition",j,"\n"))
  #output_file=paste(output_path,"/",gsub("\\\\","",gsub(" ","",gsub("/","",path))),
  #                  "_",j,".txt",sep="")
  temp_path = strsplit(path, "/")[[1]]
  output_file = paste(output_path, "/", temp_path[length(temp_path)], "_", j, ".txt", sep = "")
  file_range=((j-1)*max_n+1):min(j*max_n,n_total)
  n=length(file_range)
  data_x=array(0,dim=c(80,80,3,n))
  for (i in file_range) {data_x[,,,i-(j-1)*max_n]=readTIFF(files[i])}
  data_x_new=matrix(0,nrow=80*80*3,ncol=n)
  data_x_new[]=data_x
  data_x_new=t(data_x_new)
  data_x_new=data.frame(rep("$",n),data_x_new)
  write.table(data_x_new,file=output_file,quote=F,row.names=F,col.names=F,sep=" ")
  
  system(paste(deeplearning_R,"/test.out ",deeplearning_R,
               "/configure_3category.txt ",
               n," ",output_file,"_predicted ",output_file,sep=""))
  results=read.table(paste(output_file,"_predicted",sep=""))
  results[,1]=files[file_range]
  results=results[,-5]
  write.table(results,file=paste(output_file,"_results.txt",sep=""),quote=F,col.names=F,row.names=F,sep="\t")
  unlink(output_file)
  unlink(paste(output_file,"_predicted",sep=""))
}
