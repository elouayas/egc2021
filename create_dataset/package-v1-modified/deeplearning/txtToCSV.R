library(xlsx)
slide <- "382412-2_1"
setwd(paste("F:/perl/package/examples/image01/deeplearning_results/", slide, sep = ""))

#combine all predicted results
allfile <- list.files(full.names=T,pattern=".txt")
datalist <- NA
for (i in allfile){
  datalist <- rbind(datalist, read.delim(i, header=FALSE, stringsAsFactors = F))
}

#write.xlsx(x = datalist, file = "test.xlsx", row.names = FALSE, col.names = FALSE)

datalist$Patch.Name <- gsub("(.+\\/)", "", datalist[, 1])
#find predicted cell type
prediction <- function(x){
  return (which(x == max(x))[1])
}
datalist$predicted <- apply(datalist[, 2:4], MARGIN = 1, prediction)

#get location info
setwd("../../patch_extraction_results/cellLocationInfo")
location <- read.csv(paste(slide, ".csv", sep = ""))

#construct cell_Info file
m <- merge(datalist, location, by.x = "Patch.Name", by.y = "Patch.Name")
to.write <- data.frame(Name_of_Image_Patch = m$Patch.Name, 
                       Location.X = m$Location.in.X.axis, 
                       Location.Y = m$Location.in.Y.axis, 
                       Category = m$predicted)
