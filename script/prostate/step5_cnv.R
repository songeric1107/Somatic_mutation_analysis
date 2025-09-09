path="/Volumes/projects-t3/PTRAN/Projects_starting_Jan2024/ysong/gain_loss/spatial/perc100/"
setwd(path)

#path="/Volumes/projects-t3/PTRAN/Projects_starting_Jan2023/dnaseq/analysis/yang_cnv_test/cnv_tempus/mod_version/bed_file/bed_100Perc_match/JH001.ann.bed"
file.names <- dir(path, pattern =".bed")
out.file<-""
for(i in 1:length(file.names)){
  
  file <- read.table(file.names[i],header=F, sep="\t", stringsAsFactors=FALSE)
  if (nrow(file)!=0) {
    
    out.file <- file[,c(1:5,9)]
    out.file$sample=file.names[i]
  
    filename <- paste(file.names[i], ".id.txt", sep="")
    write.table(out.file, file=filename,sep="\t",quote=FALSE,row.names=F
    )}}
 
path="/Volumes/projects-t3/PTRAN/Projects_starting_Jan2024/ysong/gain_loss/spatial/perc100"
setwd(path)

#path="/Volumes/projects-t3/PTRAN/Projects_starting_Jan2023/dnaseq/analysis/yang_cnv_test/cnv_tempus/mod_version/bed_file/bed_100Perc_match/JH001.ann.bed"
file.names <- dir(path, pattern ="id.txt")

datalist = lapply(file.names, function(x){read.table(file=x, header=T, sep="\t", stringsAsFactors=FALSE)})
allt=Reduce(function(x,y) {merge(x,y,all=T)}, datalist)

write.table(allt,"all_combine.cnv.raw.txt",sep="\t",quote=F)



allt1=reshape2::dcast( V9~sample, data = allt, value.var = "V5", toString)
rownames(allt1)=allt1$sample
write.table(allt1,"all_combine.cnv.per100.txt",sep="\t",quote=F,row.names=T)

all.f=allt1[which(allt1$V9%in%c("PTEN","TP53","RB1","MYC")),]

#/local/projects-t3/PTRAN/Projects_starting_Jan2024/ysong/gain_loss/spatial/perc100/


# -------------------------------------------------------------------------



stat=read.table("~/Desktop/raw_list.txt",sep="\t",header=T)
stat[is.na(stat)]<-0

statm=melt(stat,id.vars=c("sample"))

pdf("~/Desktop/stat.pdf",30,10)
ggplot(statm,aes(x=sample,y=value,fill=variable))+geom_bar(position = "dodge", stat="identity")+
  theme(axis.text.x=element_text(angle = 90, size=4,hjust=1))+scale_y_continuous(expand = c(0,0))
dev.off()


allt1=read.table("../all_combine.cnv.matrix.v2.txt",sep="\t",header=T,row.names=1)

#allt1=t(allt1[-1])

col = c("Gain" = "blue", "Loss" = "red", "Neutral" = "yellow","NA"="white")

allt1[is.na(allt1)]<-""

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "white", col = NA))
  },
  "NA" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "white", col = NA))
  },
  # big blue
  Gain = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Gain"], col = NA))
  },
  # big red
  Loss = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Loss"], col = NA))
  },
  # small green
  Neutral = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),  
              gp = gpar(fill = col["Neutral"], col = NA))})
  
allt1=read.table("/Volumes/projects-t3/PTRAN/Projects_starting_Jan2023/dnaseq/analysis/yang_cnv_test/cnv_tempus/mod_version/bed_file/bed_100Perc_match/all_combine.cnv.matrix.v2.txt",sep='\t',header=T,row.names=1)
allt1[is.na(allt1)]<-0
#colnames(allt1)=gsub("[.]","-",colnames(allt1))
meta=read.delim2("/Users/ysong/Desktop/P_Tran_ref/metadata/tempus_meta.txt",sep="\t",header=T)



onco_mats=allt1[order(match(colnames(allt1),meta$Tumor_Sample_Barcode))]



pdf("../oncoprint_cnv.pdf",15,50)
oncoPrint(allt1,alter_fun = alter_fun, col = col)
dev.off()




   
    file.names <- dir(path, pattern ="ann.bed")
    Data_file <- map(file.names,read.delim, stringsAsFactors = FALSE,header=F, check.names = FALSE, row.names = NULL)
    
    Merge_All_Samples <- Data_file %>% reduce(inner_join, by = "V9")
    
    
    write.csv(Merge_All_Samples, "./Merge_All_Samples.csv", row.names = F)
    