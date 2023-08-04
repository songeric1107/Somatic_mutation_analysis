
all.maf=readRDS("tempus_foundation_w_conflicP.rds")

##########################################
#read metachronous of all samples
met.maf=readRDS("metachro.maf.417.rds")

meta=met.maf@clinical.data
meta.met=meta[which(meta$Time..2.synchronous..1.metachronous==1),]

meta.met1=meta.met[which(meta.met$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.%in%c("1","2","3")),]

#meta.syn=meta[which(meta$Time..2.synchronous..1.metachronous.==2),]

#https://shixiangwang.github.io/sigminer/reference/bp.html

library(maftools);library(sigminer)
#met.maf=subsetMaf(maf = all.maf, tsb=  meta.met$id)





mt_tally_SNV <- sig_tally(met.maf,
                          ref_genome = "BSgenome.Hsapiens.UCSC.hg38",
                          useSyn = F)


library('NMF');library(maftools)

all.tnm = trinucleotideMatrix(maf = met.maf,  add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38",useSyn = F)


mt_sig <- sig_extract(mt_tally_SNV$nmf_matrix,
                      n_sig = 2:10,
                      nrun = 50,
                      cores = 4)
                      
                      

bpall=bp_extract_signatures(mt_tally_SNV$nmf_matrix, range = 2:10)


#bp2 <- bp_extract_signatures(all.tnm$nmf_matrix, range = 2:50)

bp2 <- bp_extract_signatures(mt_tally_SNV$nmf_matrix, range = 2:10)

save(bp2,file="bp_sig.met.417.RData")


sig.sy=bp2$exposure$K6$exposure_mean




#save(bp1,file="bp_sig.met.RData")

bp1=bp2
pdf("qc_417.pdf")
bp_show_survey(bp1)
bp_show_survey2(bp1, highlight =6)

bp_show_survey(bp1)
bp_show_survey2(bp1, highlight =9)

dev.off()


s2 <- get_sig_db("SBS")

fit.sig.raw=sig_fit(t(mt_tally_SNV$nmf_matrix), sig_index = "ALL",type="relative",sig_db="SBS")
fit.sig.abs=sig_fit(t(mt_tally_SNV$nmf_matrix), sig_index = "ALL",sig_db="SBS")


t2=data.frame(t(data.frame(fit.sig.abs)))
rownames(t2)=gsub("[.]","-",rownames(t2))


meta.met$Pelvic.Node.fail.location=gsub("^0$","Pelvic.Node(NO,81)",meta.met$Pelvic.Node.fail.location)
meta.met$Pelvic.Node.fail.location=gsub("^1$","Pelvic.Node(Yes,33)",meta.met$Pelvic.Node.fail.location)

meta.met$Bone=gsub("0","Bone(NO,79)",meta.met$Bone)
meta.met$Bone=gsub("^1$","Bone(Yes,35)",meta.met$Bone)

meta.met$Distant.Node=gsub("0","Distant.node(NO,88)",meta.met$Distant.Node)
meta.met$Distant.Node=gsub("^1$","Distant.node(Yes,26)",meta.met$Distant.Node)


meta.met$Visceral=gsub("^0$","Visceral(NO,104)",meta.met$Visceral)
meta.met$Visceral=gsub("^1$","Visceral(Yes,10)",meta.met$Visceral)

meta.met$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.=gsub("^1$","oligoprogressor(57)",meta.met$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.)
meta.met$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.=gsub("^2$","polyprogressor(25)",meta.met$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.)
meta.met$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.=gsub("^3$","no.progression.at.last.fu(32)",meta.met$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.)


t2.ori=t2
tmb=data.frame(tmb(met.maf))

t2.ori$sum=rowSums(t2.ori)


t2.ori.tmb=merge(t2.ori,tmb,by.x=0,by.y="Tumor_Sample_Barcode")


t2.ori.tmb$unknown=round(t2.ori.tmb$total-t2.ori.tmb$sum)
rownames(t2.ori.tmb)=t2.ori.tmb$Row.names

sub=data.frame(t(data.frame(t2.ori.tmb[c(2:73,78)])))


rownames(sub)=gsub("[.]","-",rownames(sub))
t2=data.frame(t(sub))

t2m=melt(as.matrix(t2))
t2m$Var1=gsub("[.]","-",t2m$Var1)


#####count#######################################################

t2.rel.ann=merge(data.frame(meta.met),t2m,by.x="Tumor_Sample_Barcode",by.y="Var1")
#rownames(t2.rel.ann)=t2.rel.ann$Tumor_Sample_Barcode

set.seed(100)
library(randomcoloR)
n <- 72
palette <- c(distinctColorPalette(n),"#000000")

t2.rel.ann=t2.rel.ann[which(t2.rel.ann$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.!="na"),]




genus.m1=t2.rel.ann



sort.class <- genus.m1 %>% 
  dplyr::count(Var2, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(Var2)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(Var2 == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)


color=c("red","green","pink","blue","orange","cyan")
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node[(NO,81)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node(Yes,33)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)




pdf("count.ref.sbs.order.pattern.fail.count.w.mutation.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(Var2, levels = levels(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = Var2))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_grid(~New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.,scale="free_x")
p+ylab("relative")+ scale_fill_manual(values=palette)


dev.off()


t2.rel.ann=t2.rel.ann[which(t2.rel.ann$Pelvic.Node.fail.location!="na"),]




genus.m1=t2.rel.ann



sort.class <- genus.m1 %>% 
  dplyr::count(Var2, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(Var2)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(Var2 == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)


color=c("red","green","pink","blue","orange","cyan")
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node[(NO,81)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node(Yes,33)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)




pdf("count.ref.sbs.order.pelvic.node.count.w.mutation.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(Var2, levels = levels(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = Var2))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_grid(~Pelvic.Node.fail.location,scale="free_x")
p+ylab("relative")+ scale_fill_manual(values=palette)


dev.off()

t2.rel.ann=t2.rel.ann[which(t2.rel.ann$Distant.Node!="na"),]




genus.m1=t2.rel.ann



sort.class <- genus.m1 %>% 
  dplyr::count(Var2, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(Var2)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(Var2 == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)

pdf("count.ref.sbs.order.distant.node.count.w.mutation.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(Var2, levels = levels(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = Var2))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_grid(~Distant.Node,scale="free_x")
p+ylab("relative")+ scale_fill_manual(values=palette)

dev.off()

t2.rel.ann=t2.rel.ann[which(t2.rel.ann$Bone!="na"),]




genus.m1=t2.rel.ann



sort.class <- genus.m1 %>% 
  dplyr::count(Var2, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(Var2)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(Var2 == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)

pdf("count.ref.sbs.order.bone.count.w.mutation.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(Var2, levels = levels(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = Var2))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_grid(~Bone,scale="free_x")
p+ylab("relative")+ scale_fill_manual(values=palette)
dev.off()


t2.rel.ann=t2.rel.ann[which(t2.rel.ann$Visceral!="na"),]




genus.m1=t2.rel.ann



sort.class <- genus.m1 %>% 
  dplyr::count(Var2, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(Var2)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(Var2 == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)

pdf("count.ref.sbs.order.viseral.count.w.mutation.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)


p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(Var2, levels = levels(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = Var2))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_grid(~Visceral,scale="free_x")
p+ylab("relative")+ scale_fill_manual(values=palette)

dev.off()





##############relative########################
t2.rel= sweep(t2[1:73],1, rowSums(t2[1:73]), FUN='/')
row.names(t2.rel)=gsub("[.]","-",rownames(t2.rel))


t2.rel.ann=merge(data.frame(meta.met),t2.rel,by.x="Tumor_Sample_Barcode",by.y=0)
rownames(t2.rel.ann)=t2.rel.ann$Tumor_Sample_Barcode

set.seed(100)
library(randomcoloR)
n <- 72
palette <- c(distinctColorPalette(n),"#000000")


t2.rel.ann.m=melt(t2.rel.ann)

t2.rel.ann.m=t2.rel.ann.m[which(t2.rel.ann.m$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.!="na"),]



genus.m1=t2.rel.ann.m



sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)


color=c("red","green","pink","blue","orange","cyan")
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node[(NO,81)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node(Yes,33)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)


#Pelvic.Node.fail.location.lab=c("Pelvic.Node(NO,48)","Pelvic.Node(Yes,18)")

#names(Pelvic.Node.fail.location.lab)=c("Pelvic.Node(NO,48)","Pelvic.Node(Yes,18)")



pdf("ref.sbs.order.pattern.failure.rel.w.mutation.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = levels(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_grid(~New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.,scale="free_x")
p+ylab("relative")+ scale_fill_manual(values=palette)
dev.off()



t2.rel.ann.m=t2.rel.ann.m[which(t2.rel.ann.m$Pelvic.Node.fail.location!="na"),]



genus.m1=t2.rel.ann.m



sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)


color=c("red","green","pink","blue","orange","cyan")
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node[(NO,81)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node(Yes,33)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)


#Pelvic.Node.fail.location.lab=c("Pelvic.Node(NO,48)","Pelvic.Node(Yes,18)")

#names(Pelvic.Node.fail.location.lab)=c("Pelvic.Node(NO,48)","Pelvic.Node(Yes,18)")



pdf("ref.sbs.order.pelvic.node.count.w.mutation.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = levels(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_grid(~Pelvic.Node.fail.location,scale="free_x")
p+ylab("relative")+ scale_fill_manual(values=palette)
dev.off()

t2.rel.ann.m=t2.rel.ann.m[which(t2.rel.ann.m$Distant.Node!="na"),]



genus.m1=t2.rel.ann.m



sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)


color=c("red","green","pink","blue","orange","cyan")
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node[(NO,81)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node(Yes,33)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)


#Pelvic.Node.fail.location.lab=c("Pelvic.Node(NO,48)","Pelvic.Node(Yes,18)")

#names(Pelvic.Node.fail.location.lab)=c("Pelvic.Node(NO,48)","Pelvic.Node(Yes,18)")



pdf("ref.sbs.order.distant.node.count.w.mutation.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = levels(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_grid(~Distant.Node,scale="free_x")
p+ylab("relative")+ scale_fill_manual(values=palette)
dev.off()

t2.rel.ann.m=t2.rel.ann.m[which(t2.rel.ann.m$Bone!="na"),]



genus.m1=t2.rel.ann.m



sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)


color=c("red","green","pink","blue","orange","cyan")
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node[(NO,81)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node(Yes,33)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)



pdf("ref.sbs.order.bone.count.w.mutation.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = levels(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_grid(~Bone,scale="free_x")
p+ylab("relative")+ scale_fill_manual(values=palette)
dev.off()

t2.rel.ann.m=t2.rel.ann.m[which(t2.rel.ann.m$Visceral!="na"),]



genus.m1=t2.rel.ann.m



sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)


color=c("red","green","pink","blue","orange","cyan")
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node[(NO,81)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node(Yes,33)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)



pdf("ref.sbs.order.visceral.count.w.mutation.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = levels(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_grid(~Visceral,scale="free_x")
p+ylab("relative")+ scale_fill_manual(values=palette)
dev.off()




sub1=t2.rel.ann[which(t2.rel.ann$Distant.Node=="Distant.node(NO,88)"),]
sub2=t2.rel.ann[which(t2.rel.ann$Distant.Node=="Distant.node(Yes,26)"),]


t2.relt1=data.frame(t(sub1[55:127]))
t2.relt2=data.frame(t(sub2[55:127]))

t2.relt1[is.na(t2.relt1)]<-0
t2.relt2[is.na(t2.relt2)]<-0


t2.relt1$avg=rowMeans(t2.relt1[,])
t2.relt2$avg=rowMeans(t2.relt2[,])

t2.relt1$group="Distant.node(NO,88)"

t2.relt2$group="Distant.node(Yes,26)"

t2.relt2$sig=rownames(t2.relt2)

t2.relt1$sig=rownames(t2.relt1)

mutation.com=rbind(t2.relt1[89:91],t2.relt2[27:29])
mutation.com.f1=mutation.com[which(mutation.com$sig%in%c("SBS1","SBS10b","SBS87","SBS15","unknown")),]
mutation.com.f2=mutation.com[which(mutation.com$sig%notin%c("SBS1","SBS10b","SBS87","SBS15","unknown")),]

set.seed(100)
library(randomcoloR)
n <- 72
palette <- c(distinctColorPalette(n),"black")


# Load dplyr
library(dplyr)

# Group by sum using dplyr
agg_tbl <- mutation.com.f2 %>% group_by(group) %>% 
  summarise(avg = sum(avg),
            .groups = 'drop')
agg_tbl

# Convert tibble to df
df2 <- agg_tbl %>% as.data.frame()
df2$sig="others"
df2=df2[,c(3,1,2)]
#rownames(df2)=df2$sig

data.new=rbind(mutation.com.f1,df2)


data.new$sig=factor(data.new$sig,levels=c("SBS1" ,   "SBS10b" , "SBS15" ,  "SBS87","others" ,   "unknown"))


pdf("~/Desktop/avg_mutation.new.distance.node.pdf",10,10)
ggplot(mutation.com,aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=palette)
#ggplot(mutation.com[which(mutation.com$avg>0),],aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=palette)
ggplot(data.new,aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=c("green","pink","blue","orange","grey","black"))

dev.off()

set.seed(100)
library(randomcoloR)
n <- 72
palette <- c(distinctColorPalette(n),"black")



pdf("avg_mutation.distant.pdf")
ggplot(mutation.com,aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=palette)
#ggplot(mutation.com[which(mutation.com$avg>0),],aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=palette)

dev.off()


sub1=t2.rel.ann[which(t2.rel.ann$Pelvic.Node.fail.location=="Pelvic.Node(NO,81)"),]
sub2=t2.rel.ann[which(t2.rel.ann$Pelvic.Node.fail.location=="Pelvic.Node(Yes,33)"),]


t2.relt1=data.frame(t(sub1[55:127]))
t2.relt2=data.frame(t(sub2[55:127]))

t2.relt1[is.na(t2.relt1)]<-0
t2.relt2[is.na(t2.relt2)]<-0


t2.relt1$avg=rowMeans(t2.relt1[,])
t2.relt2$avg=rowMeans(t2.relt2[,])

t2.relt1$group="Pelvic.Node(NO,81)"

t2.relt2$group="Pelvic.Node(Yes,33)"

t2.relt2$sig=rownames(t2.relt2)

t2.relt1$sig=rownames(t2.relt1)

mutation.com=rbind(t2.relt1[82:84],t2.relt2[34:36])

`%notin%` <- Negate(`%in%`)

mutation.com.f1=mutation.com[which(mutation.com$sig%in%c("SBS1","SBS10b","SBS87","SBS15","unknown")),]
mutation.com.f2=mutation.com[which(mutation.com$sig%notin%c("SBS1","SBS10b","SBS87","SBS15","unknown")),]

set.seed(100)
library(randomcoloR)
n <- 72
palette <- c(distinctColorPalette(n),"black")


# Load dplyr
library(dplyr)

# Group by sum using dplyr
agg_tbl <- mutation.com.f2 %>% group_by(group) %>% 
  summarise(avg = sum(avg),
            .groups = 'drop')
agg_tbl

# Convert tibble to df
df2 <- agg_tbl %>% as.data.frame()
df2$sig="others"
df2=df2[,c(3,1,2)]
#rownames(df2)=df2$sig

data.new=rbind(mutation.com.f1,df2)


data.new$sig=factor(data.new$sig,levels=c("SBS1" ,   "SBS10b" , "SBS15" ,  "SBS87","others" ,   "unknown"))


pdf("avg_mutation.new.pelvic.pdf",10,10)
ggplot(mutation.com,aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=palette)
#ggplot(mutation.com[which(mutation.com$avg>0),],aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=palette)
ggplot(data.new,aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=c("green","pink","blue","orange","grey","black"))

dev.off()

sub1=t2.rel.ann[which(t2.rel.ann$Bone=="Bone(NO,79)"),]
sub2=t2.rel.ann[which(t2.rel.ann$Bone=="Bone(Yes,35)"),]


t2.relt1=data.frame(t(sub1[55:127]))
t2.relt2=data.frame(t(sub2[55:127]))

t2.relt1[is.na(t2.relt1)]<-0
t2.relt2[is.na(t2.relt2)]<-0


t2.relt1$avg=rowMeans(t2.relt1[,])
t2.relt2$avg=rowMeans(t2.relt2[,])

t2.relt1$group="Bone(NO,79)"

t2.relt2$group="Bone(Yes,35)"

t2.relt2$sig=rownames(t2.relt2)

t2.relt1$sig=rownames(t2.relt1)

mutation.com=rbind(t2.relt1[80:82],t2.relt2[36:38])


set.seed(100)
library(randomcoloR)
n <- 72
palette <- c(distinctColorPalette(n),"black")


pdf("avg_mutation.bone.pdf")
ggplot(mutation.com,aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=palette)
#ggplot(mutation.com[which(mutation.com$avg>0),],aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=palette)

dev.off()


sub1=t2.rel.ann[which(t2.rel.ann$Bone=="Bone(NO,79)"),]
sub2=t2.rel.ann[which(t2.rel.ann$Bone=="Bone(Yes,35)"),]


t2.relt1=data.frame(t(sub1[55:127]))
t2.relt2=data.frame(t(sub2[55:127]))

t2.relt1[is.na(t2.relt1)]<-0
t2.relt2[is.na(t2.relt2)]<-0


t2.relt1$avg=rowMeans(t2.relt1[,])
t2.relt2$avg=rowMeans(t2.relt2[,])

t2.relt1$group="Bone(NO,79)"

t2.relt2$group="Bone(Yes,35)"

t2.relt2$sig=rownames(t2.relt2)

t2.relt1$sig=rownames(t2.relt1)

mutation.com=rbind(t2.relt1[80:82],t2.relt2[36:38])
`%notin%` <- Negate(`%in%`)

mutation.com.f1=mutation.com[which(mutation.com$sig%in%c("SBS1","SBS10b","SBS87","SBS15","unknown")),]
mutation.com.f2=mutation.com[which(mutation.com$sig%notin%c("SBS1","SBS10b","SBS87","SBS15","unknown")),]

set.seed(100)
library(randomcoloR)
n <- 72
palette <- c(distinctColorPalette(n),"black")


# Load dplyr
library(dplyr)

# Group by sum using dplyr
agg_tbl <- mutation.com.f2 %>% group_by(group) %>% 
  summarise(avg = sum(avg),
            .groups = 'drop')
agg_tbl

# Convert tibble to df
df2 <- agg_tbl %>% as.data.frame()
df2$sig="others"
df2=df2[,c(3,1,2)]
#rownames(df2)=df2$sig

data.new=rbind(mutation.com.f1,df2)


data.new$sig=factor(data.new$sig,levels=c("SBS1" ,   "SBS10b" , "SBS15" ,  "SBS87","others" ,   "unknown"))


pdf("avg_mutation.new.bone.pdf",10,10)
ggplot(mutation.com,aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=palette)
#ggplot(mutation.com[which(mutation.com$avg>0),],aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=palette)
ggplot(data.new,aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=c("green","pink","blue","orange","grey","black"))

dev.off()

set.seed(100)
library(randomcoloR)
n <- 72
palette <- c(distinctColorPalette(n),"black")







pdf("avg_mutation.bone.pdf")
ggplot(mutation.com,aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=palette)
#ggplot(mutation.com[which(mutation.com$avg>0),],aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=palette)

dev.off()
#########################

sub1=t2.rel.ann[which(t2.rel.ann$Visceral=="Visceral(NO,104)"),]
sub2=t2.rel.ann[which(t2.rel.ann$Visceral=="Visceral(Yes,10)"),]


t2.relt1=data.frame(t(sub1[55:127]))
t2.relt2=data.frame(t(sub2[55:127]))

t2.relt1[is.na(t2.relt1)]<-0
t2.relt2[is.na(t2.relt2)]<-0


t2.relt1$avg=rowMeans(t2.relt1[,])
t2.relt2$avg=rowMeans(t2.relt2[,])

t2.relt1$group="Visceral(NO,104)"

t2.relt2$group="Visceral(Yes,10)"

t2.relt2$sig=rownames(t2.relt2)

t2.relt1$sig=rownames(t2.relt1)

mutation.com=rbind(t2.relt1[105:107],t2.relt2[11:13])

`%notin%` <- Negate(`%in%`)

mutation.com.f1=mutation.com[which(mutation.com$sig%in%c("SBS1","SBS10b","SBS87","SBS15","unknown")),]
mutation.com.f2=mutation.com[which(mutation.com$sig%notin%c("SBS1","SBS10b","SBS87","SBS15","unknown")),]

set.seed(100)
library(randomcoloR)
n <- 72
palette <- c(distinctColorPalette(n),"black")


# Load dplyr
library(dplyr)

# Group by sum using dplyr
agg_tbl <- mutation.com.f2 %>% group_by(group) %>% 
  summarise(avg = sum(avg),
            .groups = 'drop')
agg_tbl

# Convert tibble to df
df2 <- agg_tbl %>% as.data.frame()
df2$sig="others"
df2=df2[,c(3,1,2)]
#rownames(df2)=df2$sig

data.new=rbind(mutation.com.f1,df2)


data.new$sig=factor(data.new$sig,levels=c("SBS1" ,   "SBS10b" , "SBS15" ,  "SBS87","others" ,   "unknown"))


pdf("avg_mutation.new.visceral.pdf",10,10)
ggplot(mutation.com,aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=palette)
#ggplot(mutation.com[which(mutation.com$avg>0),],aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=palette)
ggplot(data.new,aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=c("green","pink","blue","orange","grey","black"))

dev.off()

set.seed(100)
library(randomcoloR)
n <- 72
palette <- c(distinctColorPalette(n),"black")


pdf("~/Desktop/avg_mutation.visceral.pdf")
ggplot(mutation.com,aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=palette)
#ggplot(mutation.com[which(mutation.com$avg>0),],aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=palette)

dev.off()


sub1=t2.rel.ann[which(t2.rel.ann$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.=="no.progression.at.last.fu(32)"),]
sub2=t2.rel.ann[which(t2.rel.ann$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.=="oligoprogressor(57)"),]
sub3=t2.rel.ann[which(t2.rel.ann$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.=="polyprogressor(25)"),]

t2.relt1=data.frame(t(sub1[55:127]))
t2.relt2=data.frame(t(sub2[55:127]))
t2.relt3=data.frame(t(sub3[55:127]))

t2.relt1[is.na(t2.relt1)]<-0
t2.relt2[is.na(t2.relt2)]<-0
t2.relt3[is.na(t2.relt3)]<-0

t2.relt1$avg=rowMeans(t2.relt1[,])
t2.relt2$avg=rowMeans(t2.relt2[,])
t2.relt3$avg=rowMeans(t2.relt3[,])

t2.relt1$group="no.progression.at.last.fu(32)"

t2.relt2$group="oligoprogressor(57)"
t2.relt3$group="polyprogressor(25)"


t2.relt2$sig=rownames(t2.relt2)

t2.relt1$sig=rownames(t2.relt1)
t2.relt3$sig=rownames(t2.relt3)

`%notin%` <- Negate(`%in%`)
mutation.com=rbind(t2.relt1[33:35],t2.relt2[58:60],t2.relt3[26:28])

mutation.com.f1=mutation.com[which(mutation.com$sig%in%c("SBS1","SBS10b","SBS87","SBS15","unknown")),]
mutation.com.f2=mutation.com[which(mutation.com$sig%notin%c("SBS1","SBS10b","SBS87","SBS15","unknown")),]

set.seed(100)
library(randomcoloR)
n <- 72
palette <- c(distinctColorPalette(n),"black")


# Load dplyr
library(dplyr)

# Group by sum using dplyr
agg_tbl <- mutation.com.f2 %>% group_by(group) %>% 
  summarise(avg = sum(avg),
            .groups = 'drop')
agg_tbl

# Convert tibble to df
df2 <- agg_tbl %>% as.data.frame()
df2$sig="others"
df2=df2[,c(3,1,2)]
#rownames(df2)=df2$sig

data.new=rbind(mutation.com.f1,df2)
data.new=data.new[order(data.new$sig),]
data.new=data.new[c(4:18,1:3),]

data.new$sig=factor(data.new$sig,levels=c("SBS1" ,   "SBS10b" , "SBS15" ,  "SBS87","others" ,   "unknown"))


pdf("avg_mutation.new.pattern.pdf",10,10)
ggplot(mutation.com,aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=palette)
#ggplot(mutation.com[which(mutation.com$avg>0),],aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=palette)
ggplot(data.new,aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=c("green","pink","blue","orange","grey","black"))

dev.off()



#######################
library(dplyr)
t2.change=t2 %>% 
  mutate_all(funs(replace(., .==0, NA))) %>% 
  transmute(posmin = pmin(!!! rlang::syms(names(.)), na.rm = TRUE)) %>%
  bind_cols(t2, .)

t2.change=t2 %>% 
  mutate_all(funs(replace(., .==0, NA))) %>% 
  transmute(posmax = pmax(!!! rlang::syms(names(.)), na.rm = TRUE)) %>%
  bind_cols(t2, .)


t2.f=t2.change[which(t2.change$posmax>1),]

t2.rel= sweep(t2[1:72],1, rowSums(t2[1:72]), FUN='/')
t2.rel.ann=merge(data.frame(meta.met),t2.rel,by.x="Tumor_Sample_Barcode",by.y=0)
rownames(t2.rel.ann)=t2.rel.ann$Tumor_Sample_Barcode

sub1=t2.rel.ann[which(t2.rel.ann$Pelvic.Node.fail.location=="Pelvic.Node(NO,81)"),]
sub2=t2.rel.ann[which(t2.rel.ann$Pelvic.Node.fail.location=="Pelvic.Node(Yes,33)"),]


t2.relt1=data.frame(t(sub1[55:126]))
t2.relt2=data.frame(t(sub2[55:126]))

t2.relt1[is.na(t2.relt1)]<-0
t2.relt2[is.na(t2.relt2)]<-0


t2.relt1$avg=rowMeans(t2.relt1[,])
t2.relt2$avg=rowMeans(t2.relt2[,])

t2.relt1$group="Pelvic.Node(NO,81)"

t2.relt2$group="Pelvic.Node(Yes,33)"

t2.relt2$sig=rownames(t2.relt2)

t2.relt1$sig=rownames(t2.relt1)

mutation.com=rbind(t2.relt1[82:84],t2.relt2[34:36])


set.seed(100)
library(randomcoloR)
n <- 74
palette <- distinctColorPalette(n)



pdf("avg_mutation.pdf")
ggplot(mutation.com,aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=palette)
ggplot(mutation.com[which(mutation.com$avg>0),],aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=palette)

dev.off()


#####################

t2.f$avg=rowMeans(t2.f[,-ncol(t2.f)])


fit.ann.ref=merge(data.frame(meta.met),t2.f,by.x="Tumor_Sample_Barcode",by.y=0,all.x=T)


fit.ann.ref$Tumor_Sample_Barcode=gsub("-T","",fit.ann.ref$Tumor_Sample_Barcode)


rownames(fit.ann.ref)=fit.ann.ref$Tumor_Sample_Barcode

fit.ann.ref=fit.ann.ref[which(fit.ann.ref$Pelvic.Node.fail.location!="na"),]

fit.ann1=fit.ann.ref[-c(2:46,51:54,127:128)]
fit.ann1=fit.ann1[which(fit.ann1$Pelvic.Node.fail.location!="na"),]








fit.ann1=fit.ann1[which(fit.ann1$posmax!="NA"),]
rownames(fit.ann1)=fit.ann1$Tumor_Sample_Barcode

fit.ann1=fit.ann1[which(rowSums(fit.ann1[c(7:78)])>0),]

input=fit.ann1[c(7:78)]
rownames(input)=fit.ann1$Tumor_Sample_Barcode

input[is.na(input)]<-0

input_matrix=t(input)
input_matrix[is.na(input_matrix)]<-0

taxmat=data.frame(rownames(input_matrix))
rownames(taxmat)=taxmat$rownames.input_matrix.


metadata=fit.ann1[1:6]

######do not use###############
OTUALL = otu_table(input_matrix, taxa_are_rows = TRUE)

TAXALL=tax_table(as.matrix(taxmat))

rownames(metadata)=metadata$Tumor_Sample_Barcode
sampledata = sample_data(metadata)

match(rownames(sampledata),colnames(input_matrix))


ps <- phyloseq(OTUALL, 
               TAXALL,sampledata)

alph.div=estimate_richness(ps,measures=c("Shannon"))
alph.div1=cbind(sample_data(ps),alph.div)

                                              
result=aov(Shannon~Pelvic.Node.fail.location,alph.div1)                                                                           
summary(result)
pdf("~/Desktop/shannon.pdf")
plot_richness(ps,x="Pelvic.Node.fail.location",color="Pelvic.Node.fail.location",measures=c(  "Shannon","Simpson"))+ geom_violin(width=1.0) +
  geom_boxplot(width=0.1, alpha=0.2) 



dev.off()


braycurtis <- vegdist(fit.ann1[c(7:78)])
mat1=as.matrix(braycurtis)[fit.ann1$Pelvic.Node.fail.location=="Pelvic.Node(NO,81)",fit.ann1$Pelvic.Node.fail.location=="Pelvic.Node(Yes,33)"]
mat2=as.dist(as.matrix(braycurtis)[fit.ann1$Pelvic.Node.fail.location=="Pelvic.Node(NO,81)",fit.ann1$Pelvic.Node.fail.location=="Pelvic.Node(Yes,33)"])

mat1m=melt(as.matrix(mat1))


bray_dist = phyloseq::distance(ps, method="bray", weighted=F)

meandist(braycurtis, df$PoolNumber)



sampledf <- data.frame(sample_data(ps ))

adonis2(bray_dist  ~ Pelvic.Node.fail.location, data = sampledf)

ordination = ordinate(ps, method="MDS", distance=bray_dist)
plot_ordination(ps, ordination,color="Pelvic.Node.fail.location") + theme(aspect.ratio=1)


pdf("cordinate.pdf")    
plot_ordination(ps, ordination,color="Pelvic.Node.fail.location") + theme(aspect.ratio=1) + geom_point(size = 5)

dev.off()
#################################################
#fit.ann1=fit.ann1[which(rowSums(fit.ann1[-c(1:6)])>0),]

fit.ann1m=melt(fit.ann1[-ncol(fit.ann1)],id.vars=c(colnames(fit.ann1[,1:6])))



data=t(fit.ann.ref[,-c(1:54)])

data.f=data[which(rowSums(data)>0),]


sub=data.f[,c(73:74,113)]

pca <- prcomp(t(data.f), 
              scale = TRUE)

pdf("pca.pelvic.pdf")
fviz_pca_biplot(pca,habillage =fit.ann.ref$Pelvic.Node.fail.location,labelsize = 2, repel = TRUE)+
  ylim(-5,5)+xlim(-10,10)+
  theme(text = element_text(size = 7.5),
        axis.title = element_text(size = 7.5),
        axis.text = element_text(size = 7.5))

dev.off()


pdf("pca.Distant.Node.pdf")

fviz_pca_biplot(pca,habillage =fit.ann.ref$Distant.Node,labelsize = 3, repel = TRUE)+
  ylim(-10,10)+xlim(-15,15)+
  theme(text = element_text(size = 7.5),
        axis.title = element_text(size = 7.5),
        axis.text = element_text(size = 7.5))
biplot(pca)
dev.off()



library(randomForest)
prob=c(0.7,0.3)





genus.m1=fit.ann1m



sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)


color=c("red","green","pink","blue","orange","cyan")
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node[(NO,81)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node(Yes,33)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)


Pelvic.Node.fail.location.lab=c("Pelvic.Node(NO,48)","Pelvic.Node(Yes,18)")

names(Pelvic.Node.fail.location.lab)=c("Pelvic.Node(NO,48)","Pelvic.Node(Yes,18)")


pdf("~/Desktop/ref.sbs.order.pelvic.node.count.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = rev(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_grid(~Pelvic.Node.fail.location,scale="free_x")
p+ylab("Est_count")+ scale_fill_manual(values=palette)
dev.off()



pdf("ref.sbs.order.pelvic.node.count.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = rev(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~Pelvic.Node.fail.location,scale="free_x")
p+ylab("Est_count")+ scale_fill_manual(values=palette)
dev.off()


fit.ann1=fit.ann.ref[-c(2:46,51:53)]
fit.ann1=fit.ann1[which(fit.ann1$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.!="na"),]
#fit.ann1=fit.ann1[which(rowSums(fit.ann1[-c(1:6)])>0),]
rownames(fit.ann1)=fit.ann1$Tumor_Sample_Barcode

library(circlize)
pdf("~/Desktop/heatmap_signature_contri.pdf",20,15)
set.seed(100)
ha1 = HeatmapAnnotation( pelvic.node=fit.ann1$Pelvic.Node.fail.location,distant.node=fit.ann1$Distant.Node,bone=fit.ann1$Bone,
                         visceral=fit.ann1$Visceral,fail_pattern=fit.ann1$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.,
col = list(fail_pattern=c("no.progression.at.last.fu(32)" = "red","oligoprogressor(57)" = "orange","polyprogressor(25)" = "blue"), pelvic.node = c("Pelvic.Node(NO,81)" = "green", "Pelvic.Node(Yes,33)" = "red"),
           bone= c("Bone(NO,79)" = "green", "Bone(Yes,35)" = "red"),distant.node= c("Distant.node(NO,88)" = "blue", "Distant.node(Yes,26)" = "red"),Visceral=c('Visceral(NO,104)'="grey","Visceral(Yes,10)"="blue")))
#ha2 = HeatmapAnnotation( high_risk=cb.n1$`High_Risk_Mutation_(1=yes,_0=no)`,risk_group=cb.n1$Risk_Group,crpc=cb.n1$`CRPC_from_oligomet_(1=yes,_0=no)`,col = list(risk_group=c("H" = "red","Metastatic" = "orange","Node+" = "blue","UI" = "green","VH" = "pink","FI"="yellow"),high_risk = c("0" = "green", "1" = "red"),crpc = c("no" = "green", "yes" = "red"),annotation_height = unit(c(5, 5, 15), "mm") ),time_to_crpc = anno_points(time_to_CRPC_from_oligomet.n, ylim = c(0, max(time_to_CRPC_from_oligomet.n, na.rm = TRUE)), axis = TRUE))



col2 = colorRamp2(c( 0, 2,10), c( "white", "orange","purple"))
Heatmap(t(fit.ann1[-c(1:6)]),col=col2,top_annotation=ha1)
Heatmap(t(fit.ann1[-c(1:6)]),col=col2,top_annotation=ha1,column_split = fit.ann1$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.)
dev.off()

fit.ann1m=melt(fit.ann1,id.vars=c(colnames(fit.ann1[,1:6])))

genus.m1=fit.ann1m

sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)


color=c("red","green","pink","blue","orange","cyan")


pdf("ref.sbs.order.pattern.failure.v2.count.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = rev(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.,scale="free_x")
p+ylab("Est_count")+ scale_fill_manual(values=palette)
dev.off()


fit.ann1=fit.ann.ref[-c(2:46,51:53)]
fit.ann1=fit.ann1[which(fit.ann1$Bone!="na"),]
#fit.ann1=fit.ann1[which(rowSums(fit.ann1[-c(1:6)])>0),]

fit.ann1m=melt(fit.ann1,id.vars=c(colnames(fit.ann1[,1:6])))

genus.m1=fit.ann1m

sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)


color=c("red","green","pink","blue","orange","cyan")


pdf("~/Desktop/ref.sbs.order.bone.v2.count.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = rev(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~Bone,scale="free_x")
p+ylab("Est_count")+ scale_fill_manual(values=palette)
dev.off()



fit.ann1=fit.ann.ref[-c(2:46,51:53)]
fit.ann1=fit.ann1[which(fit.ann1$Distant.Node!="na"),]
#fit.ann1=fit.ann1[which(rowSums(fit.ann1[-c(1:6)])>0),]

fit.ann1m=melt(fit.ann1,id.vars=c(colnames(fit.ann1[,1:6])))

genus.m1=fit.ann1m

sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)


color=c("red","green","pink","blue","orange","cyan")


pdf("ref.sbs.order.distant.node.v2.count.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = rev(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~Distant.Node,scale="free_x")
p+ylab("Est_count")+ scale_fill_manual(values=palette)
dev.off()


fit.ann1=fit.ann.ref[-c(2:46,51:53)]
fit.ann1=fit.ann1[which(fit.ann1$Visceral!="na"),]
#fit.ann1=fit.ann1[which(rowSums(fit.ann1[-c(1:6)])>0),]

fit.ann1m=melt(fit.ann1,id.vars=c(colnames(fit.ann1[,1:6])))

genus.m1=fit.ann1m

sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)


color=c("red","green","pink","blue","orange","cyan")


pdf("~/Desktop/ref.sbs.order.viseral.v2.count.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = rev(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~Visceral,scale="free_x")
p+ylab("Est_count")+ scale_fill_manual(values=palette)
dev.off()


#########exposure.count.denovo
sig.sy=bp2$exposure$K6$exposure_mean




t2=data.frame(t(data.frame(sig.sy)))
rownames(t2)=gsub("[.]","-",rownames(t2))





t2.ori=t2
tmb=data.frame(tmb(met.maf))

t2.ori$sum=rowSums(t2.ori)


t2.ori.tmb=merge(t2.ori,tmb,by.x=0,by.y="Tumor_Sample_Barcode")


t2.ori.tmb$unknown=round(t2.ori.tmb$total-t2.ori.tmb$sum)
rownames(t2.ori.tmb)=t2.ori.tmb$Row.names

sub=data.frame(t(data.frame(t2.ori.tmb[c(2:7,12)])))


t2=data.frame(t(sub))
rownames(t2)=gsub("[.]","-",rownames(t2))


t2m=melt(as.matrix(t2))
t2m$Var1=gsub("[.]","-",t2m$Var1)


###################################
t2.tmb=merge(tmb,t2,by.x="Tumor_Sample_Barcode",by.y=0)
#t2.tmb$others=t2.tmb$total-rowSums(t2.tmb[-c(1:4)])

fit.ann.ref=merge(data.frame(meta.met),t2.tmb,by="Tumor_Sample_Barcode")


fit.ann.ref$Tumor_Sample_Barcode=gsub("-T","",fit.ann.ref$Tumor_Sample_Barcode)

fit.ann1=fit.ann.ref[-c(2:46,51:53,55:57)]

fit.ann1f=fit.ann1[which(fit.ann1$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.!="na"),]
#fit.ann1=fit.ann1[which(rowSums(fit.ann1[-c(1:6)])>0),]

fit.ann1m=melt(fit.ann1f,id.vars=c(colnames(fit.ann1[,1:6])))

genus.m1=fit.ann1m

sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)


color=c("red","green","pink","blue","orange","cyan","black")

pdf("denovo.sbs.order.pattern.failer.w.unknown.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = levels(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.,scale="free_x")
p+ylab("Est_count")+ scale_fill_manual(values=color)
dev.off()


all.tnm = trinucleotideMatrix(maf = met.maf,  add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38",useSyn = F)



fit.ann1f=fit.ann1[which(fit.ann1$Pelvic.Node.fail.location!="na"),]
#fit.ann1=fit.ann1[which(rowSums(fit.ann1[-c(1:6)])>0),]

fit.ann1m=melt(fit.ann1f,id.vars=c(colnames(fit.ann1[,1:6])))

genus.m1=fit.ann1m

sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)


color=c("red","green","pink","blue","orange","cyan","black")
all.tnm = trinucleotideMatrix(maf = met.maf,  add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38",useSyn = F)


library('NMF')
laml.sign <- estimateSignatures(mat = all.tnm ,
                                nTry =100,
                                pConstant = 0.1,
                                parallel = 1)



pdf("denovo.sbs.order.pelvic.node.w.unknown.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = levels(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~Pelvic.Node.fail.location,scale="free_x")
p+ylab("Est_count")+ scale_fill_manual(values=color)
dev.off()
#################################deno,relative



t2.rel= sweep(t2,1, rowSums(t2), FUN='/')
row.names(t2.rel)=gsub("[.]","-",rownames(t2.rel))







t2.rel.ann=merge(data.frame(meta.met),t2.rel,by.x="Tumor_Sample_Barcode",by.y=0)
rownames(t2.rel.ann)=t2.rel.ann$Tumor_Sample_Barcode

set.seed(100)
library(randomcoloR)
n <- 72
palette <- c(distinctColorPalette(n),"#000000")


t2.rel.ann.m=melt(t2.rel.ann)
t2.rel.ann.m=t2.rel.ann.m[which(t2.rel.ann.m$Pelvic.Node.fail.location!="na"),]


pdf("~/Desktop/deno_rel_mutation_w_unknown.pdf",100,10)
ggplot(t2.rel.ann.m,aes(x=Tumor_Sample_Barcode,y=value,fill=variable))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~Pelvic.Node.fail.location,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=palette)
#ggplot(t2.rel.ann.m[which(t2.rel.ann.m$avg>0),],aes(x=Pelvic.Node.fail.location,y=value,fill=variable))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~Pelvic.Node.fail.location,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=palette)

dev.off()


genus.m1=t2.rel.ann.m



sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)


#color=c("red","green","pink","blue","orange","cyan")
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node[(NO,81)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node(Yes,33)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)


#Pelvic.Node.fail.location.lab=c("Pelvic.Node(NO,48)","Pelvic.Node(Yes,18)")

#names(Pelvic.Node.fail.location.lab)=c("Pelvic.Node(NO,48)","Pelvic.Node(Yes,18)")



pdf("~/Desktop/denovo.rel.sbs.order.pelvic.node.count.w.mutation.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = levels(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_grid(~Pelvic.Node.fail.location,scale="free_x")
p+ylab("propertion")+ scale_fill_manual(values=color)
dev.off()



sub1=t2.rel.ann[which(t2.rel.ann$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.=="no.progression.at.last.fu(32)"),]
sub2=t2.rel.ann[which(t2.rel.ann$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.=="oligoprogressor(57)"),]
sub3=t2.rel.ann[which(t2.rel.ann$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.=="polyprogressor(25)"),]


t2.relt1=data.frame(t(sub1[55:61]))
t2.relt2=data.frame(t(sub2[55:61]))
t2.relt3=data.frame(t(sub3[55:61]))

t2.relt1[is.na(t2.relt1)]<-0
t2.relt2[is.na(t2.relt2)]<-0
t2.relt3[is.na(t2.relt3)]<-0

t2.relt1$avg=rowMeans(t2.relt1[,])
t2.relt2$avg=rowMeans(t2.relt2[,])
t2.relt3$avg=rowMeans(t2.relt3[,])

t2.relt1$group="no.progression.at.last.fu(32) "

t2.relt2$group="oligoprogressor(57)"
t2.relt3$group="polyprogressor(25)"


t2.relt2$sig=rownames(t2.relt2)

t2.relt1$sig=rownames(t2.relt1)
t2.relt3$sig=rownames(t2.relt3)


mutation.com=rbind(t2.relt1[31:33],t2.relt2[57:59],t2.relt3[26:28])


set.seed(100)
library(randomcoloR)
n <- 72
palette <- c(distinctColorPalette(n),"black")



pdf("avg_mutation_pattern_denovo.rel.pdf")
ggplot(mutation.com,aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=color)
#ggplot(mutation.com[which(mutation.com$avg>0),],aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=color)

dev.off()



sub1=t2.rel.ann[which(t2.rel.ann$Pelvic.Node.fail.location=="Pelvic.Node(NO,81)"),]
sub2=t2.rel.ann[which(t2.rel.ann$Pelvic.Node.fail.location=="Pelvic.Node(Yes,33)"),]


t2.relt1=data.frame(t(sub1[55:61]))
t2.relt2=data.frame(t(sub2[55:61]))

t2.relt1[is.na(t2.relt1)]<-0
t2.relt2[is.na(t2.relt2)]<-0


t2.relt1$avg=rowMeans(t2.relt1[,])
t2.relt2$avg=rowMeans(t2.relt2[,])

t2.relt1$group="Pelvic.Node(NO,81)"

t2.relt2$group="Pelvic.Node(Yes,33)"

t2.relt2$sig=rownames(t2.relt2)

t2.relt1$sig=rownames(t2.relt1)

mutation.com=rbind(t2.relt1[79:81],t2.relt2[34:36])


set.seed(100)
library(randomcoloR)
n <- 72
palette <- c(distinctColorPalette(n),"black")



pdf("avg_mutation_denovo.rel.pdf")
ggplot(mutation.com,aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=color)
ggplot(mutation.com[which(mutation.com$avg>0),],aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=color)

dev.off()




###############distant_node

t2.rel.ann.m=t2.rel.ann.m[which(t2.rel.ann.m$Distant.Node!="na"),]





genus.m1=t2.rel.ann.m



sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)


#color=c("red","green","pink","blue","orange","cyan")
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node[(NO,81)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node(Yes,33)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)


#Pelvic.Node.fail.location.lab=c("Pelvic.Node(NO,48)","Pelvic.Node(Yes,18)")

#names(Pelvic.Node.fail.location.lab)=c("Pelvic.Node(NO,48)","Pelvic.Node(Yes,18)")



pdf("denovo.rel.sbs.order.distant.node.count.w.mutation.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = levels(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_grid(~Distant.Node,scale="free_x")
p+ylab("propertion")+ scale_fill_manual(values=color)
dev.off()





sub1=t2.rel.ann[which(t2.rel.ann$Distant.Node=="Distant.node(NO,88)"),]
sub2=t2.rel.ann[which(t2.rel.ann$Distant.Node=="Distant.node(Yes,26)"),]


t2.relt1=data.frame(t(sub1[55:61]))
t2.relt2=data.frame(t(sub2[55:61]))

t2.relt1[is.na(t2.relt1)]<-0
t2.relt2[is.na(t2.relt2)]<-0


t2.relt1$avg=rowMeans(t2.relt1[,])
t2.relt2$avg=rowMeans(t2.relt2[,])

t2.relt1$group="Distant.node(NO,88)"

t2.relt2$group="Distant.node(Yes,26)"

t2.relt2$sig=rownames(t2.relt2)

t2.relt1$sig=rownames(t2.relt1)

mutation.com=rbind(t2.relt1[86:88],t2.relt2[27:29])


set.seed(100)
library(randomcoloR)
n <- 72
palette <- c(distinctColorPalette(n),"black")



pdf("avg_mutation_denovo.distant_node.rel.pdf")
ggplot(mutation.com,aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=color)
ggplot(mutation.com[which(mutation.com$avg>0),],aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=color)

dev.off()


##################bone


t2.rel.ann.m=t2.rel.ann.m[which(t2.rel.ann.m$Bone!="na"),]





genus.m1=t2.rel.ann.m



sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)


#color=c("red","green","pink","blue","orange","cyan")
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node[(NO,81)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node(Yes,33)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)


#Pelvic.Node.fail.location.lab=c("Pelvic.Node(NO,48)","Pelvic.Node(Yes,18)")

#names(Pelvic.Node.fail.location.lab)=c("Pelvic.Node(NO,48)","Pelvic.Node(Yes,18)")



pdf("denovo.rel.sbs.order.bone.count.w.mutation.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = levels(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_grid(~Bone,scale="free_x")
p+ylab("propertion")+ scale_fill_manual(values=color)
dev.off()





sub1=t2.rel.ann[which(t2.rel.ann$Bone=="Bone(NO,79)"),]
sub2=t2.rel.ann[which(t2.rel.ann$Bone=="Bone(Yes,35)"),]


t2.relt1=data.frame(t(sub1[55:61]))
t2.relt2=data.frame(t(sub2[55:61]))

t2.relt1[is.na(t2.relt1)]<-0
t2.relt2[is.na(t2.relt2)]<-0


t2.relt1$avg=rowMeans(t2.relt1[,])
t2.relt2$avg=rowMeans(t2.relt2[,])

t2.relt1$group="Bone(NO,79)"

t2.relt2$group="Bone(Yes,35)"

t2.relt2$sig=rownames(t2.relt2)

t2.relt1$sig=rownames(t2.relt1)

mutation.com=rbind(t2.relt1[78:80],t2.relt2[35:37])


set.seed(100)
library(randomcoloR)
n <- 72
palette <- c(distinctColorPalette(n),"black")



pdf("avg_mutation_denovo.bone.rel.pdf")
ggplot(mutation.com,aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=color)
ggplot(mutation.com[which(mutation.com$avg>0),],aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=color)

dev.off()


########################viseral################################################

t2.rel.ann.m=t2.rel.ann.m[which(t2.rel.ann.m$Visceral!="na"),]





genus.m1=t2.rel.ann.m



sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)


#color=c("red","green","pink","blue","orange","cyan")
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node[(NO,81)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node(Yes,33)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)


#Pelvic.Node.fail.location.lab=c("Pelvic.Node(NO,48)","Pelvic.Node(Yes,18)")

#names(Pelvic.Node.fail.location.lab)=c("Pelvic.Node(NO,48)","Pelvic.Node(Yes,18)")



pdf("denovo.rel.sbs.order.visceral.count.w.mutation.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = levels(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_grid(~Visceral,scale="free_x")
p+ylab("propertion")+ scale_fill_manual(values=color)
dev.off()





sub1=t2.rel.ann[which(t2.rel.ann$Visceral=="Visceral(NO,104)"),]
sub2=t2.rel.ann[which(t2.rel.ann$Visceral=="Visceral(Yes,10)"),]


t2.relt1=data.frame(t(sub1[55:61]))
t2.relt2=data.frame(t(sub2[55:61]))

t2.relt1[is.na(t2.relt1)]<-0
t2.relt2[is.na(t2.relt2)]<-0


t2.relt1$avg=rowMeans(t2.relt1[,])
t2.relt2$avg=rowMeans(t2.relt2[,])

t2.relt1$group="Visceral(NO,104)"

t2.relt2$group="Visceral(Yes,10)"

t2.relt2$sig=rownames(t2.relt2)

t2.relt1$sig=rownames(t2.relt1)

mutation.com=rbind(t2.relt1[102:104],t2.relt2[11:13])


set.seed(100)
library(randomcoloR)
n <- 72
palette <- c(distinctColorPalette(n),"black")



pdf("avg_mutation_denovo.visceral.rel.pdf")
ggplot(mutation.com,aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=color)
ggplot(mutation.com[which(mutation.com$avg>0),],aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=color)

dev.off()





#######################################################
t2.rel.ann.m=t2.rel.ann.m[which(t2.rel.ann.m$Distant.Node!="na"),]





genus.m1=t2.rel.ann.m



sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)


#color=c("red","green","pink","blue","orange","cyan")
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node[(NO,81)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)
#genus.m1$Pelvic.Node.fail.location=gsub("Pelvic.Node(Yes,33)","Pelvic.Node(NO,48)",genus.m1$Pelvic.Node.fail.location)


#Pelvic.Node.fail.location.lab=c("Pelvic.Node(NO,48)","Pelvic.Node(Yes,18)")

#names(Pelvic.Node.fail.location.lab)=c("Pelvic.Node(NO,48)","Pelvic.Node(Yes,18)")



pdf("~/Desktop/denovo.rel.sbs.order.distant.node.count.w.mutation.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = levels(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_grid(~Distant.Node,scale="free_x")
p+ylab("propertion")+ scale_fill_manual(values=color)
dev.off()





sub1=t2.rel.ann[which(t2.rel.ann$Distant.Node=="Distant.node(NO,88)"),]
sub2=t2.rel.ann[which(t2.rel.ann$Distant.Node=="Distant.node(Yes,26)"),]


t2.relt1=data.frame(t(sub1[55:61]))
t2.relt2=data.frame(t(sub2[55:61]))

t2.relt1[is.na(t2.relt1)]<-0
t2.relt2[is.na(t2.relt2)]<-0


t2.relt1$avg=rowMeans(t2.relt1[,])
t2.relt2$avg=rowMeans(t2.relt2[,])

t2.relt1$group="Distant.node(NO,88)"

t2.relt2$group="Distant.node(Yes,26)"

t2.relt2$sig=rownames(t2.relt2)

t2.relt1$sig=rownames(t2.relt1)

mutation.com=rbind(t2.relt1[86:88],t2.relt2[27:29])


set.seed(100)
library(randomcoloR)
n <- 72
palette <- c(distinctColorPalette(n),"black")



pdf("~/Desktop/avg_mutation_denovo.distant_node.rel.pdf")
ggplot(mutation.com,aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=color)
ggplot(mutation.com[which(mutation.com$avg>0),],aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=color)

dev.off()









###################
library(dplyr)
t2.change=t2 %>% 
  mutate_all(funs(replace(., .==0, NA))) %>% 
  transmute(posmin = pmin(!!! rlang::syms(names(.)), na.rm = TRUE)) %>%
  bind_cols(t2, .)

t2.change=t2 %>% 
  mutate_all(funs(replace(., .==0, NA))) %>% 
  transmute(posmax = pmax(!!! rlang::syms(names(.)), na.rm = TRUE)) %>%
  bind_cols(t2, .)


t2.f=t2.change[which(t2.change$posmax>1),]

t2.rel= sweep(t2[1:72],1, rowSums(t2[1:72]), FUN='/')
t2.rel.ann=merge(data.frame(meta.met),t2.rel,by.x="Tumor_Sample_Barcode",by.y=0)
rownames(t2.rel.ann)=t2.rel.ann$Tumor_Sample_Barcode

sub1=t2.rel.ann[which(t2.rel.ann$Pelvic.Node.fail.location=="Pelvic.Node(NO,81)"),]
sub2=t2.rel.ann[which(t2.rel.ann$Pelvic.Node.fail.location=="Pelvic.Node(Yes,33)"),]


t2.relt1=data.frame(t(sub1[55:126]))
t2.relt2=data.frame(t(sub2[55:126]))

t2.relt1[is.na(t2.relt1)]<-0
t2.relt2[is.na(t2.relt2)]<-0


t2.relt1$avg=rowMeans(t2.relt1[,])
t2.relt2$avg=rowMeans(t2.relt2[,])

t2.relt1$group="Pelvic.Node(NO,81)"

t2.relt2$group="Pelvic.Node(Yes,33)"

t2.relt2$sig=rownames(t2.relt2)

t2.relt1$sig=rownames(t2.relt1)

mutation.com=rbind(t2.relt1[82:84],t2.relt2[34:36])


set.seed(100)
library(randomcoloR)
n <- 74
palette <- distinctColorPalette(n)



pdf("avg_mutation.pdf")
ggplot(mutation.com,aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=palette)
ggplot(mutation.com[which(mutation.com$avg>0),],aes(x=group,y=avg,fill=sig))+geom_bar(position = "stack",stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~group,scale="free_x")+ylab("Proportion")+ scale_fill_manual(values=palette)

dev.off()












####################################################


sig.k6=bp2$signature$K6
sig.ext=sig.k6$signature_mean





#########exposure.count

t2=data.frame(t(data.frame(sig.sy)))
rownames(t2)=gsub("[.]","-",rownames(t2))

fit.ann.ref=merge(data.frame(meta.met),t2,by.x="Tumor_Sample_Barcode",by.y=0,all.x=T)


fit.ann.ref$Tumor_Sample_Barcode=gsub("-T","",fit.ann.ref$Tumor_Sample_Barcode)

fit.ann1=fit.ann.ref[-c(2:46,51:53)]
fit.ann1=fit.ann1[which(fit.ann1$Bone!="na"),]
#fit.ann1=fit.ann1[which(rowSums(fit.ann1[-c(1:6)])>0),]

fit.ann1m=melt(fit.ann1,id.vars=c(colnames(fit.ann1[,1:6])))

genus.m1=fit.ann1m

sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)


color=c("red","green","pink","blue","orange","cyan")


pdf("ref.sbs.order.bone.count.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = rev(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~Bone,scale="free_x")
p+ylab("Est_count")+ scale_fill_manual(values=color)
dev.off()




###visceral
t2=data.frame(t(data.frame(sig.sy)))
rownames(t2)=gsub("[.]","-",rownames(t2))

fit.ann.ref=merge(data.frame(meta.met),t2,by.x="Tumor_Sample_Barcode",by.y=0,all.x=T)


fit.ann.ref$Tumor_Sample_Barcode=gsub("-T","",fit.ann.ref$Tumor_Sample_Barcode)

fit.ann1=fit.ann.ref[-c(2:46,51:53)]
fit.ann1=fit.ann1[which(fit.ann1$Visceral!="na"),]
#fit.ann1=fit.ann1[which(rowSums(fit.ann1[-c(1:6)])>0),]

fit.ann1m=melt(fit.ann1,id.vars=c(colnames(fit.ann1[,1:6])))

genus.m1=fit.ann1m

sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)


color=c("red","green","pink","blue","orange","cyan")


pdf("~/Desktop/ref.sbs.order.visceral.count.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = rev(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~Visceral,scale="free_x")
p+ylab("Est_count")+ scale_fill_manual(values=color)
dev.off()

################
t2=data.frame(t(data.frame(sig.sy)))
rownames(t2)=gsub("[.]","-",rownames(t2))

fit.ann.ref=merge(data.frame(meta.met),t2,by.x="Tumor_Sample_Barcode",by.y=0,all.x=T)


fit.ann.ref$Tumor_Sample_Barcode=gsub("-T","",fit.ann.ref$Tumor_Sample_Barcode)

fit.ann1=fit.ann.ref[-c(2:46,51:53)]
fit.ann1=fit.ann1[which(fit.ann1$Bone!="na"),]
#fit.ann1=fit.ann1[which(rowSums(fit.ann1[-c(1:6)])>0),]

fit.ann1m=melt(fit.ann1,id.vars=c(colnames(fit.ann1[,1:6])))

genus.m1=fit.ann1m

sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)


color=c("red","green","pink","blue","orange","cyan")


pdf("ref.sbs.order.bone.count.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = rev(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~Bone,scale="free_x")
p+ylab("Est_count")+ scale_fill_manual(values=color)
dev.off()


######


t2=data.frame(t(data.frame(sig.sy)))
rownames(t2)=gsub("[.]","-",rownames(t2))

fit.ann.ref=merge(data.frame(meta.met),t2,by.x="Tumor_Sample_Barcode",by.y=0,all.x=T)


fit.ann.ref$Tumor_Sample_Barcode=gsub("-T","",fit.ann.ref$Tumor_Sample_Barcode)

fit.ann1=fit.ann.ref[-c(2:46,51:53)]
fit.ann1=fit.ann1[which(fit.ann1$Distant.Node!="na"),]
#fit.ann1=fit.ann1[which(rowSums(fit.ann1[-c(1:6)])>0),]

fit.ann1m=melt(fit.ann1,id.vars=c(colnames(fit.ann1[,1:6])))

genus.m1=fit.ann1m

sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)


color=c("red","green","pink","blue","orange","cyan")


pdf("~/Desktop/ref.sbs.order.distant.node.count.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = rev(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~Distant.Node,scale="free_x")
p+ylab("Est_count")+ scale_fill_manual(values=color)
dev.off()

######

t2=data.frame(t(data.frame(sig.sy)))
rownames(t2)=gsub("[.]","-",rownames(t2))

fit.ann.ref=merge(data.frame(meta.met),t2,by.x="Tumor_Sample_Barcode",by.y=0,all.x=T)


fit.ann.ref$Tumor_Sample_Barcode=gsub("-T","",fit.ann.ref$Tumor_Sample_Barcode)

fit.ann1=fit.ann.ref[-c(2:46,51:53)]
fit.ann1=fit.ann1[which(fit.ann1$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.!="na"),]
#fit.ann1=fit.ann1[which(rowSums(fit.ann1[-c(1:6)])>0),]

fit.ann1m=melt(fit.ann1,id.vars=c(colnames(fit.ann1[,1:6])))

genus.m1=fit.ann1m

sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)


color=c("red","green","pink","blue","orange","cyan")


pdf("ref.sbs.order.pattern.failure.count.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = rev(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.,scale="free_x")
p+ylab("Est_count")+ scale_fill_manual(values=color)
dev.off()



#######




fit.ann.ref=merge(data.frame(meta.met),t2,by.x="Tumor_Sample_Barcode",by.y=0)


fit.ann.ref$Tumor_Sample_Barcode=gsub("-T","",fit.ann.ref$Tumor_Sample_Barcode)

fit.ann1=fit.ann.ref[-c(2:46,51:53)]
fit.ann1=fit.ann1[which(fit.ann1$Pelvic.Node.fail.location!="na"),]
#fit.ann1=fit.ann1[which(rowSums(fit.ann1[-c(1:6)])>0),]

fit.ann1m=melt(fit.ann1,id.vars=c(colnames(fit.ann1[,1:6])))

genus.m1=fit.ann1m

sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)


set.seed(100)
library(randomcoloR)
n <- 74
palette <- distinctColorPalette(n)


pdf("ref.sbs.order.pelvic.node.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = rev(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~Pelvic.Node.fail.location,scale="free_x")
p+ylab("Est_count")+ scale_fill_manual(values=palette)
dev.off()


###bone


fit.ann1=fit.ann.ref[-c(2:46,51:53)]
fit.ann1=fit.ann1[which(fit.ann1$Bone!="na"),]
fit.ann1m=melt(fit.ann1,id.vars=c(colnames(fit.ann1[,1:6])))

genus.m1=fit.ann1m

sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)

library(randomcoloR)
#n <- 74
#palette <- distinctColorPalette(n)


pdf("ref.sbs.order.Bone.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = rev(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~Bone,scale="free_x")
p+ylab("Proportion")+ scale_fill_manual(values=palette)
dev.off()


############



fit.ann1=fit.ann.ref[-c(2:46,51:53)]
fit.ann1=fit.ann1[which(fit.ann1$Visceral!="na"),]
fit.ann1m=melt(fit.ann1,id.vars=c(colnames(fit.ann1[,1:6])))

genus.m1=fit.ann1m

sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)

#library(randomcoloR)
#n <- 74
#palette <- distinctColorPalette(n)


pdf("ref.sbs.order.visceral.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = rev(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~Visceral,scale="free_x")
p+ylab("Proportion")+ scale_fill_manual(values=palette)
dev.off()


###distance.node


fit.ann1=fit.ann.ref[-c(2:46,51:53)]
fit.ann1=fit.ann1[which(fit.ann1$Distant.Node!="na"),]
fit.ann1m=melt(fit.ann1,id.vars=c(colnames(fit.ann1[,1:6])))

genus.m1=fit.ann1m

sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)

library(randomcoloR)
#n <- 74
#palette <- distinctColorPalette(n)


pdf("ref.sbs.order.distance.node.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = rev(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~Distant.Node,scale="free_x")
p+ylab("Proportion")+ scale_fill_manual(values=palette)
dev.off()


####pattern failure

fit.ann1=fit.ann.ref[-c(2:46,51:53)]
fit.ann1=fit.ann1[which(fit.ann1$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.!="na"),]
fit.ann1m=melt(fit.ann1,id.vars=c(colnames(fit.ann1[,1:6])))

genus.m1=fit.ann1m

sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)

library(randomcoloR)
#n <- 74
#palette <- distinctColorPalette(n)


pdf("ref.sbs.order.pattern.failure.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = rev(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.,scale="free_x")
p+ylab("Proportion")+ scale_fill_manual(values=palette)
dev.off()




#############

fit.ann1=fit.ann.ref[-c(2:46,51:53)]
fit.ann1=fit.ann1[which(fit.ann1$Pelvic.Node.fail.location!="na"),]

fit.ann1.dup=fit.ann1

`%notin%` <- Negate(`%in%`)

name=c("SBS10b","SBS1","SBS87","SBS15","SBS32")

mat1=fit.ann1.dup[,which(colnames(fit.ann1.dup)%in%name)]

mat2=fit.ann1.dup[,which(colnames(fit.ann1.dup)%notin%name)]
mat2$other=rowSums(mat2[-c(1:6)])

mat.new=cbind(mat2[c(1:6,ncol(mat2))],mat1)

#fit.ann1.dup$others=rowSums(fit.ann1.dup[-c(1:6),which(colnames(fit.ann1.dup)%notin%name)])
fit.new.m=melt(mat.new,id.vars=c(colnames(mat.new[,1:6])))


genus.m1=fit.new.m

sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)

#library(randomcoloR)
#set.seed(10)
#n <- 74
#palette <- distinctColorPalette(n)


pdf("ref.sbs.match.pelvic_node.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = rev(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~Pelvic.Node.fail.location,scale="free_x")
p+ylab("Proportion")+ scale_fill_manual(values=c("blue","red","green","yellow","purple","grey"))
dev.off()



###

fit.ann1=fit.ann.ref[-c(2:46,51:53)]
fit.ann1=fit.ann1[which(fit.ann1$Visceral!="na"),]

fit.ann1.dup=fit.ann1

`%notin%` <- Negate(`%in%`)

name=c("SBS10b","SBS1","SBS87","SBS15","SBS32")

mat1=fit.ann1.dup[,which(colnames(fit.ann1.dup)%in%name)]

mat2=fit.ann1.dup[,which(colnames(fit.ann1.dup)%notin%name)]
mat2$other=rowSums(mat2[-c(1:6)])

mat.new=cbind(mat2[c(1:6,ncol(mat2))],mat1)

#fit.ann1.dup$others=rowSums(fit.ann1.dup[-c(1:6),which(colnames(fit.ann1.dup)%notin%name)])
fit.new.m=melt(mat.new,id.vars=c(colnames(mat.new[,1:6])))


genus.m1=fit.new.m

sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)

#library(randomcoloR)
#set.seed(10)
#n <- 74
#palette <- distinctColorPalette(n)


pdf("ref.sbs.match.visceral.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = rev(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~Visceral,scale="free_x")
p+ylab("Proportion")+ scale_fill_manual(values=c("blue","red","green","yellow","purple","grey"))
dev.off()


######

fit.ann1=fit.ann.ref[-c(2:46,51:53)]
fit.ann1=fit.ann1[which(fit.ann1$Bone!="na"),]

fit.ann1.dup=fit.ann1

`%notin%` <- Negate(`%in%`)

name=c("SBS10b","SBS1","SBS87","SBS15","SBS32")

mat1=fit.ann1.dup[,which(colnames(fit.ann1.dup)%in%name)]

mat2=fit.ann1.dup[,which(colnames(fit.ann1.dup)%notin%name)]
mat2$other=rowSums(mat2[-c(1:6)])

mat.new=cbind(mat2[c(1:6,ncol(mat2))],mat1)

#fit.ann1.dup$others=rowSums(fit.ann1.dup[-c(1:6),which(colnames(fit.ann1.dup)%notin%name)])
fit.new.m=melt(mat.new,id.vars=c(colnames(mat.new[,1:6])))


genus.m1=fit.new.m

sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)

#library(randomcoloR)
#set.seed(10)
#n <- 74
#palette <- distinctColorPalette(n)


pdf("ref.sbs.match.bone.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = rev(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~Bone,scale="free_x")
p+ylab("Proportion")+ scale_fill_manual(values=c("blue","red","green","yellow","purple","grey"))
dev.off()


###

fit.ann1=fit.ann.ref[-c(2:46,51:53)]
fit.ann1=fit.ann1[which(fit.ann1$Distant.Node!="na"),]

fit.ann1.dup=fit.ann1

`%notin%` <- Negate(`%in%`)

name=c("SBS10b","SBS1","SBS87","SBS15","SBS32")

mat1=fit.ann1.dup[,which(colnames(fit.ann1.dup)%in%name)]

mat2=fit.ann1.dup[,which(colnames(fit.ann1.dup)%notin%name)]
mat2$other=rowSums(mat2[-c(1:6)])

mat.new=cbind(mat2[c(1:6,ncol(mat2))],mat1)

#fit.ann1.dup$others=rowSums(fit.ann1.dup[-c(1:6),which(colnames(fit.ann1.dup)%notin%name)])
fit.new.m=melt(mat.new,id.vars=c(colnames(mat.new[,1:6])))


genus.m1=fit.new.m

sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)

#library(randomcoloR)
#set.seed(10)
#n <- 74
#palette <- distinctColorPalette(n)


pdf("ref.sbs.match.distant.node.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = rev(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~Distant.Node,scale="free_x")
p+ylab("Proportion")+ scale_fill_manual(values=c("blue","red","green","yellow","purple","grey"))
dev.off()



#############################################vis.deno.vo

fit.ann1=fit.ann.ref[-c(2:46,51:53)]
fit.ann1=fit.ann1[which(fit.ann1$Visceral!="na"),]

fit.ann1.dup=fit.ann1

`%notin%` <- Negate(`%in%`)

name=c("SBS10b","SBS1","SBS87","SBS15","SBS32")

mat1=fit.ann1.dup[,which(colnames(fit.ann1.dup)%in%name)]

mat2=fit.ann1.dup[,which(colnames(fit.ann1.dup)%notin%name)]
mat2$other=rowSums(mat2[-c(1:6)])

mat.new=cbind(mat2[c(1:6,ncol(mat2))],mat1)

#fit.ann1.dup$others=rowSums(fit.ann1.dup[-c(1:6),which(colnames(fit.ann1.dup)%notin%name)])
fit.new.m=melt(mat.new,id.vars=c(colnames(mat.new[,1:6])))


genus.m1=fit.new.m

sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)

#library(randomcoloR)
#set.seed(10)
#n <- 74
#palette <- distinctColorPalette(n)


pdf("ref.sbs.match.patter.failure.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = rev(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.,scale="free_x")
p+ylab("Proportion")+ scale_fill_manual(values=c("blue","red","green","yellow","purple","grey"))
dev.off()

####################







fit.ann1=fit.ann.ref[-c(2:46,51:53)]
fit.ann1=fit.ann1[which(fit.ann1$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.!="na"),]

fit.ann1.dup=fit.ann1

`%notin%` <- Negate(`%in%`)

name=c("SBS10b","SBS1","SBS87","SBS15","SBS32")

mat1=fit.ann1.dup[,which(colnames(fit.ann1.dup)%in%name)]

mat2=fit.ann1.dup[,which(colnames(fit.ann1.dup)%notin%name)]
mat2$other=rowSums(mat2[-c(1:6)])

mat.new=cbind(mat2[c(1:6,ncol(mat2))],mat1)

#fit.ann1.dup$others=rowSums(fit.ann1.dup[-c(1:6),which(colnames(fit.ann1.dup)%notin%name)])
fit.new.m=melt(mat.new,id.vars=c(colnames(mat.new[,1:6])))


genus.m1=fit.new.m

sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)

#library(randomcoloR)
#set.seed(10)
#n <- 74
#palette <- distinctColorPalette(n)


pdf("ref.sbs.match.patter.failure.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = rev(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.,scale="free_x")
p+ylab("Proportion")+ scale_fill_manual(values=c("blue","red","green","yellow","purple","grey"))
dev.off()






###################pattern9,facet_wrap

fit.sig.abs=sig_fit(t(mt_tally_SNV$nmf_matrix), sig_index = "ALL",sig_db="SBS")

t2=data.frame(t(data.frame(fit.sig.abs)))
rownames(t2)=gsub("[.]","-",rownames(t2))


fit.ann.ref=merge(data.frame(meta.met),t2,by.x="Tumor_Sample_Barcode",by.y=0)


fit.ann.ref$Tumor_Sample_Barcode=gsub("-T","",fit.ann.ref$Tumor_Sample_Barcode)



fit.ann1=fit.ann.ref[-c(2:46,51:53)]

pdf("ref.sbs.all.count.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=105.9/2.54 * 1,pagecentre=FALSE)
fit.ann1t=melt(fit.ann1,id.vars=colnames(fit.ann1)[1:6])
  p=ggplot(fit.ann1t,aes(x = Tumor_Sample_Barcode, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_grid(variable~New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.,scale="free")
p+ylab("Est_count")+ scale_fill_manual(values=palette)
dev.off()

###################

fit.ann1=fit.ann1[which(fit.ann1$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.!="na"),]

fit.ann1.dup=fit.ann1

`%notin%` <- Negate(`%in%`)

name=c("SBS10b","SBS1","SBS87","SBS15","SBS32","SBS30","SBS5","SBS11")

mat1=fit.ann1.dup[,which(colnames(fit.ann1.dup)%in%name)]

mat2=fit.ann1.dup[,which(colnames(fit.ann1.dup)%notin%name)]
mat2$other=rowSums(mat2[-c(1:6)])

mat.new=cbind(mat2[c(1:6,ncol(mat2))],mat1)

#fit.ann1.dup$others=rowSums(fit.ann1.dup[-c(1:6),which(colnames(fit.ann1.dup)%notin%name)])
fit.new.m=melt(mat.new,id.vars=c(colnames(mat.new[,1:6])))


genus.m1=fit.new.m

sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)

#library(randomcoloR)
#set.seed(10)
#n <- 74
#palette <- distinctColorPalette(n)


pdf("ref.sbs.match.patter.failure.pattern9.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = rev(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.,scale="free_x")
p+ylab("est_count")+ scale_fill_manual(values=c("blue","red","green","yellow","purple","lightslateblue","darkorange","lightcoral","darkslategray"))
dev.off()












###
genus.m1 %>% 
  group_by(Tumor_Sample_Barcode) %>% 
  distinct(Tumor_Sample_Barcode, Pelvic.Node.fail.location,variable) %>%
  ungroup %>%
  summarise(pval = chisq.test(Pelvic.Node.fail.location,variable)$p.value)


#fit.ann1.dup$sum=rowSums(fit.ann1.dup[-c(1:6)] > 0)

####################################

fit2m=fit.ann1m[which(fit.ann1m$variable%in%c("SBS10b","SBS1","SBS87","SBS15","SBS32")),]



###find the top signature

unique(head(fit.ann1m[order(fit.ann1m$value, decreasing= T),], n = 30)$variable)







myColors<-c("deeppink3","chartreuse4","mediumorchid1","darkblue","darkorange","lightcoral","mediumorchid4","darkseagreen","orangered","cyan","pink","blue","grey","red","orange","yellow",
            "darkviolet","hotpink4","magenta","goldenrod4","lightslateblue","grey","aquamarine4","tan3","steelblue4","dodgerblue1","darkslategray","black",
            "lightcoral","mediumorchid4","darkseagreen","orangered","cyan","black","pink","blue","grey","red","orange","yellow")

#######################visceral
fit.ann1=fit.ann.ref[-c(2:46,51:53)]
 fit.ann1=fit.ann1[which(fit.ann1$Visceral!="na"),]
fit.ann1m=melt(fit.ann1,id.vars=c(colnames(fit.ann1[,1:6])))

genus.m1=fit.ann1m

sort.class <- genus.m1 %>% 
  dplyr::count(variable, wt = value ) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::pull(variable)

library(dplyr)
ID.order <- genus.m1 %>%
  dplyr::filter(variable == sort.class[1]) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::pull(Tumor_Sample_Barcode)

library(randomcoloR)
n <- 74
palette <- distinctColorPalette(n)


pdf("ref.sbs.order.visceral.pdf",paper="special",  pointsize=10,  width=55.7/2.54,height=15.9/2.54 * 1,pagecentre=FALSE)
p=genus.m1 %>%
  mutate(sample = factor(Tumor_Sample_Barcode, levels = ID.order)) %>%
  mutate(variable = factor(variable, levels = rev(sort.class))) %>%
  ggplot(aes(x = sample, y = value, fill = variable))+ labs(fill = "variable") +theme(axis.text.x = element_text(angle = 90, size = 10,hjust = 1,vjust=1),axis.text.y = element_text(size=10))+xlab("")+
  geom_bar(stat = "identity", color="black", width = 0.6) +  scale_y_continuous(expand = c(0, 0))+ylab(element_text(size=16))+facet_wrap(~Visceral,scale="free_x")
p+ylab("Proportion")+ scale_fill_manual(values=palette)
dev.off()






################visceral



fit.ann1=fit.ann.ref[-c(2:46,51:53)]
fit.ann1=fit.ann1[which(fit.ann1$Visceral!="na"),]
#fit.ann1=fit.ann1[which(rowSums(fit.ann1[-c(1:6)])>0),]

fit.ann1.dup=fit.ann1

`%notin%` <- Negate(`%in%`)

name=c("SBS10b","SBS1","SBS87","SBS15","SBS32")

mat1=fit.ann1.dup[,which(colnames(fit.ann1.dup)%in%name)]

mat2=fit.ann1.dup[,which(colnames(fit.ann1.dup)%notin%name)]
mat2$other=rowSums(mat2[-c(1:6)])

mat.new=cbind(mat2[c(1:6,ncol(mat2))],mat1)

#fit.ann1.dup$others=rowSums(fit.ann1.dup[-c(1:6),which(colnames(fit.ann1.dup)%notin%name)])
fit.new.m=melt(mat.new,id.vars=c(colnames(mat.new[,1:6])))














##################################

write.table(fit.sig.raw,"raw_signature.txt",sep='\t',quote=F)


sub=bp1$object$K6


fit.sig=sig_fit(t(mt_tally_SNV$nmf_matrix), sig = (sub$Signature),sig_db="SBS")

fit.sig0=sig_fit(t(mt_tally_SNV$nmf_matrix), sig = (sub$Signature),sig_db="SBS",type="relative")

t2=data.frame(t(data.frame(fit.sig0)))
rownames(t2)=gsub("[.]","-",rownames(t2))

fit.ann=merge(data.frame(meta.met),t2,by.x="Tumor_Sample_Barcode",by.y=0)

fit.ann1=fit.ann[-c(2:46,51:53)]

sub1=aggregate( fit.ann1[7:12], by=list( fit.ann1$Pelvic.Node.fail.location), FUN=median)

sub1a=aggregate( fit.ann1[7:12], by=list( fit.ann1$Pelvic.Node.fail.location), FUN=mean)
 sub1a$Group.1=gsub("0","Pelvic.Node(NO)",sub1a$Group.1)
 sub1a$Group.1=gsub("1","Pelvic.Node(Yes)",sub1a$Group.1)



rownames(sub1a)=sub1a$Group.1
 sub1a=sub1a[-1]
 sub1ab=  sweep(sub1a,1, rowSums(sub1a), FUN='/')

sub1am=melt(as.matrix(sub1ab[1:2,]))



pdf("pelvic_feature.pdf")

# Grouped

ggplot(sub1am, aes(fill=Var2, y=value, x=Var1)) + 
  geom_bar(position="stack", stat="identity") +
  #scale_fill_viridis(discrete = T) +
  ggtitle("mutation_signature_pelvic_node") +

  xlab("")+ylab("Proportion")+labs(fill = "Signatures")
dev.off()



####

sub1a=aggregate( fit.ann1[7:12], by=list( fit.ann1$Bone), FUN=mean)
sub1a$Group.1=gsub("0","Bone(NO)",sub1a$Group.1)
sub1a$Group.1=gsub("1","Bone(Yes)",sub1a$Group.1)



rownames(sub1a)=sub1a$Group.1
sub1a=sub1a[-1]
sub1ab=  sweep(sub1a,1, rowSums(sub1a), FUN='/')

sub1am=melt(as.matrix(sub1ab[1:2,]))



pdf("~/Desktop/bone_feature.pdf")

# Grouped

ggplot(sub1am, aes(fill=Var2, y=value, x=Var1)) + 
  geom_bar(position="stack", stat="identity") +
  #scale_fill_viridis(discrete = T) +
  ggtitle("mutation_signature_Bone") +
  
  xlab("")+ylab("Proportion")+labs(fill = "Signatures")
dev.off()

#####
sub1a=aggregate( fit.ann1[7:12], by=list( fit.ann1$Distant.Node), FUN=mean)
sub1a$Group.1=gsub("0","Distant.Node(NO)",sub1a$Group.1)
sub1a$Group.1=gsub("1","Distant.Node(Yes)",sub1a$Group.1)



rownames(sub1a)=sub1a$Group.1
sub1a=sub1a[-1]
sub1ab=  sweep(sub1a,1, rowSums(sub1a), FUN='/')

sub1am=melt(as.matrix(sub1ab[1:2,]))



pdf("distant.node_feature.pdf")

# Grouped

ggplot(sub1am, aes(fill=Var2, y=value, x=Var1)) + 
  geom_bar(position="stack", stat="identity") +
  #scale_fill_viridis(discrete = T) +
  ggtitle("mutation_signature_Distant.node") +
  
  xlab("")+ylab("Proportion")+labs(fill = "Signatures")
dev.off()






###
sub1=aggregate( fit.ann1[7:12], by=list( fit.ann1$Visceral), FUN=median)

sub1a=aggregate( fit.ann1[7:12], by=list( fit.ann1$Visceral), FUN=mean)
sub1a$Group.1=gsub("0","Visceral(NO)",sub1a$Group.1)
sub1a$Group.1=gsub("1","Visceral(Yes)",sub1a$Group.1)



rownames(sub1a)=sub1a$Group.1
sub1a=sub1a[-1]
sub1ab=  sweep(sub1a,1, rowSums(sub1a), FUN='/')

sub1am=melt(as.matrix(sub1ab[1:2,]))



pdf("visceral_feature.pdf")

# Grouped

ggplot(sub1am, aes(fill=Var2, y=value, x=Var1)) + 
  geom_bar(position="stack", stat="identity") +
  #scale_fill_viridis(discrete = T) +
  ggtitle("mutation_signature_Visceral") +
  
  xlab("")+ylab("Proportion")+labs(fill = "Signatures")
dev.off()




#####


sub1a=aggregate( fit.ann1[7:12], by=list( fit.ann1$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.), FUN=mean)
sub1a$Group.1=gsub("1","oligoprogressor",sub1a$Group.1)
sub1a$Group.1=gsub("2","polyprogressor",sub1a$Group.1)
sub1a$Group.1=gsub("3","no.progression.at.last.fu",sub1a$Group.1)


rownames(sub1a)=sub1a$Group.1
sub1a=sub1a[-1]
sub1ab=  sweep(sub1a,1, rowSums(sub1a), FUN='/')

sub1am=melt(as.matrix(sub1ab[,]))



pdf("pattern_failure.pdf")

# Grouped

ggplot(sub1am, aes(fill=Var2, y=value, x=Var1)) + 
  geom_bar(position="stack", stat="identity") +
  #scale_fill_viridis(discrete = T) +
  ggtitle("mutation_signature_pattern_failure") +
  
  xlab("")+ylab("Proportion")+labs(fill = "Signatures")
dev.off()



fit.sig.ref=sig_fit(t(mt_tally_SNV$nmf_matrix), sig = (sub$Signature),sig_db="SBS")



test1=get_groups(sub,method="samples",n_cluster=2)

show_sig_fit(fit.sig0)


sub=bp1$object$K6
#Besides de novo signature discovery shown in previous chapters, another common task is that you have gotten some reference signatures (either from known database like COSMIC or de novo discovery step), you want to know how these signatures contribute (fit) in a sample. That’s the target of sig_fit().


fit.sig=sig_fit(t(mt_tally_SNV$nmf_matrix), sig = (sub$Signature),sig_db="SBS")

show_sig_fit(fit.sig)

## Exposure optimized by sig_fit


H_estimate <- apply(bp2$sig$K6$signature_mean, 2, function(x) x / sum(x)) %*% sig_fit(t(mt_tally_SNV$nmf_matrix), sig = bp2$sig$K6$signature_mean)



H_dt_rel <- sig_fit(t(mt_tally_SNV$nmf_matrix), sig = (sub$Signature), return_class = "data.table", type = "relative")
z <- get_groups(H_dt_rel, method = "k-means")
show_groups(z)


library(sigminer)


obj_suggested <- bp_get_sig_obj(bp2, bp1$suggested)
obj_suggested


obj <- bp_get_sig_obj(bp2, 6)

obj1 <- bp_get_sig_obj(bp2, 9)



p11 <- show_sig_profile(obj, mode = "SBS", y_tr = function(x) x * 100)

pdf("~/Desktop/contri.pdf")
p11
dev.off()


pdf("~/Desktop/comatic.417.sig6.pdf")
show_sig_profile(obj, mode = "SBS", style = "cosmic")

sim <- get_sig_similarity(obj, sig_db = "SBS")

if (require(pheatmap)) {
  pheatmap::pheatmap(sim$similarity)
}

dev.off()

pdf("signature6.417.pdf",20,5)
show_sig_profile(obj, mode = "SBS", style = "cosmic")
obj <- bp_get_sig_obj(bp1, 6)
sim <- get_sig_similarity(obj, sig_db = "SBS")

p=show_sig_profile(obj, mode = "SBS", style = "cosmic")
add_labels(p,x=0.42,y=0.2,y_end=0.85,n_label=6,labels=c(paste(sim$best_match$Sig1$best_match,sim$best_match$Sig1$aetiology,sep="_"),paste(sim$best_match$Sig2$best_match,sim$best_match$Sig2$aetiology,sep="_"),
                                                         paste(sim$best_match$Sig3$best_match,sim$best_match$Sig2$aetiology,sep="_"),paste(sim$best_match$Sig4$best_match,sim$best_match$Sig4$aetiology,sep="_"),
                                                         paste(sim$best_match$Sig5$best_match,sim$best_match$Sig5$aetiology,sep="_"),paste(sim$best_match$Sig6$best_match,sim$best_match$Sig6$aetiology,sep="_")))


#mt_grps <- get_groups(sigs.all7, method = "consensus", match_consensus = TRUE)
dev.off()



pdf("signature6.417.sig9.pdf",20,10)
show_sig_profile(obj1, mode = "SBS", style = "cosmic")
obj <- bp_get_sig_obj(bp1, 9)
sim <- get_sig_similarity(obj, sig_db = "SBS")

p=show_sig_profile(obj, mode = "SBS", style = "cosmic")
add_labels(p,x=0.42,y=0.1,y_end=0.95,n_label=9,labels=c(paste(sim$best_match$Sig1$best_match,sim$best_match$Sig1$aetiology,sep="_"),paste(sim$best_match$Sig2$best_match,sim$best_match$Sig2$aetiology,sep="_"),
                                                        paste(sim$best_match$Sig3$best_match,sim$best_match$Sig2$aetiology,sep="_"),paste(sim$best_match$Sig4$best_match,sim$best_match$Sig4$aetiology,sep="_"),
                                                        paste(sim$best_match$Sig5$best_match,sim$best_match$Sig5$aetiology,sep="_"),paste(sim$best_match$Sig6$best_match,sim$best_match$Sig6$aetiology,sep="_"),
                                                        paste(sim$best_match$Sig7$best_match,sim$best_match$Sig5$aetiology,sep="_"),paste(sim$best_match$Sig8$best_match,sim$best_match$Sig6$aetiology,sep="_"),
                                                        paste(sim$best_match$Sig9$best_match,sim$best_match$Sig5$aetiology,sep="_")))


#mt_grps <- get_groups(sigs.all7, method = "consensus", match_consensus = TRUE)
dev.off()




obj0 <- bp_get_sig_obj(bp1, 6)
sim0 <- get_sig_similarity(obj0, sig_db = "SBS")


obj <- bp_get_sig_obj(bp1, 6)
sim <- get_sig_similarity(obj, sig_db = "SBS")

sub=bp1$object$K6


signature=data.frame(t(sub$Exposure))
rownames(signature)=gsub("[.]","-",rownames(signature))


sig.all.matrixt1=merge(data.frame(meta.met1),signature,by.x="Tumor_Sample_Barcode",by.y=0)

write.table(sig.all.matrixt1,"mutation_signature_6_417.txt",sep="\t",quote=F)

sig.all.matrixt1m=melt(sig.all.matrixt1,id.vars=as.character(colnames(sig.all.matrixt1)[1:54]))

sig.all.matrixt1m$Time..2.synchronous..1.metachronous=gsub("1","metachronous",sig.all.matrixt1m$Time..2.synchronous..1.metachronous)

sig.all.matrixt1m$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.=gsub("1","OLIGO",sig.all.matrixt1m$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.)
sig.all.matrixt1m$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.=gsub("2","POLYPROGRESSOR",sig.all.matrixt1m$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.)

sig.all.matrixt1m$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.=gsub("3","NO.PROGRESSION.AT LAST FU",sig.all.matrixt1m$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.)

library(rstatix)

df=sig.all.matrixt1m

stat.test <- df %>%
  group_by(variable) %>%
  wilcox_test(value ~ New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test
write.table(stat.test,"stat.test.sig6.pattern_failure.417.txt",sep="\t",quote=F)



pdf('signature2_box.all.sig6.pattern_failure.pdf',10,3)

#ggplot(sig.all.matrixt1m,aes(x=variable,y=value,fill=Time..2.synchronous..1.metachronous.))+ geom_boxplot()


#ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$progression.within.1.year!="ex"),],aes(x=variable,y=value,fill=progression.within.1.year))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6))  
#ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$progression.within.2.year!="ex"),],aes(x=variable,y=value,fill=progression.within.2.year))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 
#ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$crpc.within.1y!="ex"),],aes(x=variable,y=value,fill=crpc.within.1y))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 
#ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$crpc.with.2y!="ex"),],aes(x=variable,y=value,fill=crpc.with.2y))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 
#ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$Deek.Paterns.of.Failure..1.LTC..2.olgio..3.poly.0.blank.!="0"),],aes(x=variable,y=value,fill=Deek.Paterns.of.Failure..1.LTC..2.olgio..3.poly.0.blank.))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 
ggplot(sig.all.matrixt1m,aes(x=Tumor_Sample_Barcode,y=value,fill=variable))+ geom_bar(position="stack", stat="identity")+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) +facet_wrap(~New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.)
#ggplot(data, aes(fill=condition, y=value, x=specie)) + geom_bar(position="stack", stat="identity")
dev.off()

exposure=data.frame(get_sig_exposure(sub))


sig.all.matrixt1m$Pelvic.Node.fail.location=gsub("1","Pelvic.node",sig.all.matrixt1m$Pelvic.Node.fail.location)
sig.all.matrixt1m$Distant.Node=gsub("1","Distant.node",sig.all.matrixt1m$Distant.Node)
sig.all.matrixt1m$Bone=gsub("1","Bone",sig.all.matrixt1m$Bone)
sig.all.matrixt1m$Visceral=gsub("1","Visceral",sig.all.matrixt1m$Visceral)


pdf('signature2_failure.location.417.pdf',10,3)
ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$Pelvic.Node.fail.location%in%c(0,"Pelvic.node")),],aes(x=variable,y=value,fill=Pelvic.Node.fail.location))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 
ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$Distant.Node%in%c(0,"Distant.node")),],aes(x=variable,y=value,fill=Distant.Node))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 
ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$Bone%in%c(0,"Bone")),],aes(x=variable,y=value,fill=Bone))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 
ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$Visceral%in%c(0,"Visceral")),],aes(x=variable,y=value,fill=Visceral))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 


#ggplot(sig.all.matrixt1m,aes(x=variable,y=value,fill=Time..2.synchronous..1.metachronous.))+ geom_boxplot()


#ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$progression.within.1.year!="ex"),],aes(x=variable,y=value,fill=progression.within.1.year))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6))  
#ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$progression.within.2.year!="ex"),],aes(x=variable,y=value,fill=progression.within.2.year))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 
#ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$crpc.within.1y!="ex"),],aes(x=variable,y=value,fill=crpc.within.1y))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 
#ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$crpc.with.2y!="ex"),],aes(x=variable,y=value,fill=crpc.with.2y))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 
#ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$Deek.Paterns.of.Failure..1.LTC..2.olgio..3.poly.0.blank.!="0"),],aes(x=variable,y=value,fill=Deek.Paterns.of.Failure..1.LTC..2.olgio..3.poly.0.blank.))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 
#ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$New%in%c(0,"Pelvic.node")),],aes(x=variable,y=value,fill=New.Pattern.o.))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 

dev.off()




pel.uniq=isect$Pelvic.Node


dist.uniq=isect$`Distance Node`

bone.uniq=isect$Bone
vis.uniq=isect$Visceral

uniq.id=data.frame(c(pel.uniq,dist.uniq,bone.uniq,vis.uniq))
names(uniq.id)="id"

meta.uniq=data.frame(meta.metachro[which(meta.metachro$id%in%uniq.id$id),])

meta.uniq$level <- rep(NA, nrow(meta.uniq))

meta.uniq[meta.uniq$Pelvic.Node.fail.location%in%c("1"), ][, "level"] <- "Pelvic.node"

meta.uniq[meta.uniq$Distant.Node%in%c("1"), ][, "level"] <- "Distant.node"
meta.uniq[meta.uniq$Bone%in%c("1"), ][, "level"] <- "Bone"
meta.uniq[meta.uniq$Visceral%in%c("1"), ][, "level"] <- "Visceral"



sub0=merge(meta.uniq[,c(1,55)],sig.all.matrixt1[,c(1,55:60)],by.x="Tumor_Sample_Barcode",by.y="Tumor_Sample_Barcode")

sub0[is.na(sub0)]<-0

rownames(sub0)=sub0$Tumor_Sample_Barcode
sub0=sub0[order(sub0$level),]

sub0t=t(sub0[-c(1:2)])
level=data.frame(sub0$level)

names(level)="group"
ha1=HeatmapAnnotation(location_failure=level$group,col = list(location_failure= c("Bone" = "green", "Distant.node" = "red","Pelvic.node"="orange","Visceral"="yellow")))

scaled_mat = t(scale(t(sub0t)))

pdf("mutation.sig.heatmap.pdf",20,5)
Heatmap(scaled_mat,top_annotation=ha1)
Heatmap(scaled_mat,column_order=colnames(sub0t),top_annotation=ha1)
dev.off()


sub0m=melt(sub0,id.vars=as.character(colnames(sub0)[1:2]))



pdf("~/Desktop/signature.uniq.pdf")
ggplot(sub0m,aes(x=variable,y=value,fill=level))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 
#ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$Distant.Node%in%c(0,"Distant.Node")),],aes(x=variable,y=value,fill=Pelvic.Node.fail.location))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 
#ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$Bone%in%c(0,"Bone")),],aes(x=variable,y=value,fill=Pelvic.Node.fail.location))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 
#ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$Visceral%in%c(0,"Visceral")),],aes(x=variable,y=value,fill=Pelvic.Node.fail.location))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 
dev.off()

#ggplot(sig.all.matrixt1m,aes(x=variable,y=value,fill=Time..2.synchronous..1.metachronous.))+ geom_boxplot()


#ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$progression.within.1.year!="ex"),],aes(x=variable,y=value,fill=progression.within.1.year))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6))  
#ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$progression.within.2.year!="ex"),],aes(x=variable,y=value,fill=progression.within.2.year))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 
#ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$crpc.within.1y!="ex"),],aes(x=variable,y=value,fill=crpc.within.1y))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 
#ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$crpc.with.2y!="ex"),],aes(x=variable,y=value,fill=crpc.with.2y))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 
#ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$Deek.Paterns.of.Failure..1.LTC..2.olgio..3.poly.0.blank.!="0"),],aes(x=variable,y=value,fill=Deek.Paterns.of.Failure..1.LTC..2.olgio..3.poly.0.blank.))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 
ggplot(sig.all.matrixt1m,aes(x=variable,y=value,fill=New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.0.blank.))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 

dev.off()



library("ggpubr");library(rstatix)

pdf('signature2_box.all.pdf',10,3)
sub1=sig.all.matrixt1m[which(sig.all.matrixt1m$New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.0.blank.%in%c(1,2,3)),]

#ggplot(sig.all.matrixt1m,aes(x=variable,y=value,fill=Time..2.synchronous..1.metachronous.))+ geom_boxplot()
# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.0.blank.", dodge = 0.8)
bxp + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
)


#ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$progression.within.1.year!="ex"),],aes(x=variable,y=value,fill=progression.within.1.year))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6))  
#ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$progression.within.2.year!="ex"),],aes(x=variable,y=value,fill=progression.within.2.year))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 
#ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$crpc.within.1y!="ex"),],aes(x=variable,y=value,fill=crpc.within.1y))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 
#ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$crpc.with.2y!="ex"),],aes(x=variable,y=value,fill=crpc.with.2y))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 
#ggplot(sig.all.matrixt1m[which(sig.all.matrixt1m$Deek.Paterns.of.Failure..1.LTC..2.olgio..3.poly.0.blank.!="0"),],aes(x=variable,y=value,fill=Deek.Paterns.of.Failure..1.LTC..2.olgio..3.poly.0.blank.))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 
bxp=ggplot(sig.all.matrixt1m,aes(x=variable,y=value,fill=New.Pattern.of.Failure..1.oligoprogressor..2.polyprogressor..3.no.progression.at.last.fu.0.blank.))+ geom_boxplot()+theme(legend.position = "bottom")+ theme(legend.text=element_text(size=6)) 
bxp + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
)

dev.off()



#############################################

Maf_dt = met.maf@data


Mut_dt = Maf_dt[Variant_Type == "SNP", .(Total = .N), by = Tumor_Sample_Barcode]

obj1 <- bp_get_sig_obj(bp1, 6)

expo=get_sig_exposure(obj1)

#Mut_expo = get_sig_exposure(sigs.all)[, .(sample, Est=Sig1+Sig2)]

gene.exp=get_sig_exposure(obj1)[,.(gene_symbo,Sig1,Sig2,Sig3,Sig4,Sig5,Sig6)]

sample_expo = get_sig_exposure(obj1)[, .(sample,Sig1,Sig2,Sig3,Sig4,Sig5,Sig6)]



sig.mat.ann=merge(meta,sig.mat,by.x="id",by.y="sample")



pdf("gene_contri.sig.pdf")
gene.sig1=gene.exp[,-c(1,2,3,4)]
rownames(gene.sig1)=gene.exp$gene_symbo.Gene.Symbol
Heatmap()
dev.off()





pdf("signature.all.bootrap.consense.330.pdf",20,10)
show_sig_exposure( sub)

sim <- get_sig_similarity(sub, sig_db = "SBS")

p=show_sig_profile(sub, mode = "SBS", style = "cosmic")

add_labels(p,x=0.42,y=0.15,y_end=0.95,n_label=6,labels=c(paste(sim$best_match$Sig1$best_match,sim$best_match$Sig1$aetiology,sep="_"),paste(sim$best_match$Sig2$best_match,sim$best_match$Sig2$aetiology,sep="_"),
                                                         paste(sim$best_match$Sig3$best_match,sim$best_match$Sig2$aetiology,sep="_"),paste(sim$best_match$Sig4$best_match,sim$best_match$Sig4$aetiology,sep="_"),
                                                         paste(sim$best_match$Sig5$best_match,sim$best_match$Sig5$aetiology,sep="_"),paste(sim$best_match$Sig6$best_match,sim$best_match$Sig6$aetiology,sep="_"))
)


#mt_grps <- get_groups(sigs.all7, method = "consensus", match_consensus = TRUE)
dev.off()

