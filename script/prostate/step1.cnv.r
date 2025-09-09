setwd("/local/projects-t3/PTRAN/Projects_starting_Jan2023/dnaseq/analysis/yang_cnv_test/cnv_tempus")
path="/local/projects-t3/PTRAN/Projects_starting_Jan2023/dnaseq/analysis/yang_cnv_test/cnv_tempus"
file.names <- dir(path, pattern =".txt")

out.file<-""
for(i in 1:length(file.names)){
	  
	  
	  #https://www.bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/oncoplots.html
	  
	  file=read.table(file.names[i],sep="\t",header=T)
  file$Barcode=file.names[i]
    file1=file[,c(ncol(file),1:(ncol(file)-1))]
    out.file=file1[,c(1:4,13:14)]
      
      filename <- paste(file.names[i], ".mod.txt", sep="")
      
      
      write.table(out.file, file=filename,sep="\t",quote=FALSE,row.names=F)}


