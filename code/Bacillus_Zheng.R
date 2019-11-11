#unique gene----
sub_new <- read.table("all_135/number_of_new_genes.Rtab")
sub_new <- sub_new[,2:135]
sub_unique <- read.table("all_135/number_of_unique_genes.Rtab")
sub_unique <- sub_unique[,2:135]
cluster1_new <- read.table("pangenome1/number_of_new_genes.Rtab")
cluster1_new <- cluster1_new[,2:14]
cluster1_unique <- read.table("pangenome1/number_of_unique_genes.Rtab")
cluster1_unique <- cluster1_unique[,2:14]
cluster2_new <- read.table("pangenome2/number_of_new_genes.Rtab")
cluster2_new <- cluster2_new[,2:67]
cluster2_unique <- read.table("pangenome2/number_of_unique_genes.Rtab")
cluster2_unique <- cluster2_unique[,2:67]
cluster3_new <- read.table("pangenome3/number_of_new_genes.Rtab")
cluster3_new <- cluster3_new[,2:52]
cluster3_unique <- read.table("pangenome3/number_of_unique_genes.Rtab")
cluster3_unique <- cluster3_unique[,2:52]

par(mfrow = c(1, 1),mar=c(6,7,5,3),lwd=3)
boxplot(cluster2_unique,xlim=c(0,65),ylim=c(0,12000),border = mycolor[1],
        ylab=expression(bold('Unique protein count')),lty=1,cex.lab=1.8,font.lab=2,outline=F,axes=F)
par(new=T)
boxplot(cluster3_unique,xlim=c(0,65),ylim=c(0,12000),border = mycolor[2],col = mycolor[2],
        ylab=expression(bold('Unique protein count')),lty=1,cex.lab=1.8,font.lab=2,outline=F,axes=F)
par(new=T)
axis(2,lwd=3,lwd.ticks = 3,cex.axis=1.5,font.axis=2)
axis(1,lwd=3,lwd.ticks = 3,cex.axis=1.3,font.axis=2)
legend(
  x = "topleft",
  legend = c("Lineage 1", "Lineage 2"),
  bty = "n",
  fill = c(mycolor[1],mycolor[2]),
  border = NA,
  text.font = 2
)

library(ggplot2)
library(RColorBrewer)
mycolor<-brewer.pal(9,"Set1")
#distribution of virulence gene----
name <- read.table("name",header = F,as.is = T)
name[which(name$V2=="ff0000"),"V2"] <- 1
name[which(name$V2=="008aff"),"V2"] <- 2
name[which(name$V2=="0db8c8"),"V2"] <- 3
name[which(name$V2=="18e9bb"),"V2"] <- 4
name[which(name$V2=="e98d18"),"V2"] <- 5
name[which(name$V2=="4b00ff"),"V2"] <- 6
name[which(name$V2=="b700ff"),"V2"] <- 7
name[which(name$V2=="ff00bf"),"V2"] <- 8
name <- name[order(name$V2),]
name1 <- name[which(name$V2 !=1 & name$V2 !=2),]
name2 <- name[which(name$V2 %in% c(3,4,5)),]
name3 <- name[which(name$V2 %in% c(6,7,8)),]
sub_content <- read.table("gene_presence_absence.Rtab",header = T)
sub_content_n <- sub_content
rownames(sub_content_n) <- sub_content_n$Gene
sub_content_n <- sub_content_n[,2:116]
sub_content_n["plc",] <- 1
sub_content_n["group_2663",] <- 1
sub_content_n["plcA","Y16G1233DL_spades"]
sub_content_n["plcA","Y16G1234DL_spades"]
sub_content_n["plcA","Y16G1230DL_spades"] <- 1
sub_content_n["plcA","Y15G1082DL_spades"] <- 1
sub_content_n["plcA","Y15G1084DL_spades"]
sub_content_n["ina_1",] <- 1
sub_content_n["ina_2",] <- 1
sub_content_n["ina_3",] <- 1
sub_content_n["plcR",] <- 1
VF <- read.table("VF1",header = F,sep="\t")
sub_content_n_vf <- sub_content_n[which(rownames(sub_content_n) %in% VF$V2),]
sub_content_n_vf <- sub_content_n_vf[order(match(rownames(sub_content_n_vf), VF$V2)),]
sub_content_n_vf <- sub_content_n_vf[,which(colnames(sub_content_n_vf) %in% name1$V1)]
sub_content_n_vf <- sub_content_n_vf[, order(match(colnames(sub_content_n_vf), name1$V1))]
sub_content_n_vf_t <- as.data.frame(t(sub_content_n_vf))
sub_content_n_vf_t$sum <- apply(sub_content_n_vf_t,1,sum)
sub_content_n_vf_t <- sub_content_n_vf_t[order(sub_content_n_vf_t$sum),]

annotation_row <- data.frame(lineage=paste("SC",name1$V2,sep = ""))
rownames(annotation_row) <- name1$V1
annotation_col <- data.frame(VF=VF$V1)
rownames(annotation_col) <- VF$V2
Paired <- brewer.pal(10,"RdYlBu")
mycolor<-brewer.pal(9,"Set1")
color <- list(lineage= c(SC3 = Paired[1],
                         SC4 = Paired[3],
                         SC5 = Paired[5],
                         SC6 = Paired[10],
                         SC7 = Paired[8],
                         SC8 = Paired[6]),
              VF = c(VF_Enzyme = mycolor[1],
                     VF_Immune_evasion = mycolor[2],
                     VF_Immune_inhibitor = mycolor[3],
                     VF_Internalin = mycolor[4],
                     VF_Regulation = mycolor[5],
                     VF_Toxin = mycolor[6]))
pheatmap(as.matrix(sub_content_n_vf_t[,1:33]),
         cluster_rows = F,
         cluster_cols = F,
         show_colnames = F,
         show_rownames = F,
         border_color=NA,
         annotation_row = annotation_row,
         annotation_col = annotation_col,
         legend = T,
         legend_breaks = c(0,1),
         legend_labels = c("loss","gain"),
         color = mycolor[2:1],
         annotation_colors = color)


order_cluster1 <- name1[which(name1$V2 %in% c(3,4,5)),]
sub_content_n_vf_cluster1 <- sub_content_n_vf[,which(colnames(sub_content_n_vf) %in% order_cluster1$V1)]
sub_content_n_vf_cluster1 <- apply(sub_content_n_vf_cluster1,2,sum)
order_cluster2 <- name1[which(name1$V2 %in% c(6,7,8)),]
sub_content_n_vf_cluster2 <- sub_content_n_vf[,which(colnames(sub_content_n_vf) %in% order_cluster2$V1)]
sub_content_n_vf_cluster2 <- apply(sub_content_n_vf_cluster2,2,sum)
boxplot(sub_content_n_vf_cluster1,sub_content_n_vf_cluster2)
pdf(file = "vf_numbers.pdf")
par(mfrow = c(1, 1),
    mar = c(6, 15, 5, 10),
    lwd = 3)
boxplot(
  sub_content_n_vf_cluster1,
  sub_content_n_vf_cluster2,
  ylab = expression(bold('Number of virulence factors')),
  ylim = c(10,26),
  border = c("red", "blue"),
  cex.lab = 2,
  font.lab = 2,
  outline = F,
  axes = F
)
stripchart(
  list(sub_content_n_vf_cluster1,
       sub_content_n_vf_cluster2),
  vertical = TRUE,
  method = "jitter",
  pch = 20,
  cex = 2,
  col = c("red", "blue"),
  bg = "bisque",
  at = c(1, 2),
  add = TRUE
)
axis(
  2,
  lwd = 3,
  lwd.ticks = 3,
  cex.axis = 1.5,
  font.axis = 2
)
dev.off()
a <- wilcox.test(sub_content_n_vf_cluster1, sub_content_n_vf_cluster2)
#enrichment----
up <- replicate(1000, {
  b <- sample(name$V1, 8)
  c <- length(intersect(b, name2$V1))
})
h2 = data.frame(table(up))
hist(up)
mp <- barplot(
  h2$Freq,
  ylab = expression(bold('Frequency')),
  xlab = expression(bold('Gene number')),
  cex.lab = 1.5,
  font.lab = 2,
  col = mycolor[1],
  border = F,
  main = "",
  axes = F
)
axis(
  2,
  lwd = 3,
  lwd.ticks = 3,
  cex.axis = 1.5,
  font.axis = 2
)
axis(
  1,
  at=mp,
  labels = as.character(h2$up),
  lwd = 3,
  lwd.ticks = 3,
  cex.axis = 1.5,
  font.axis = 2
)
chisq.test(matrix(c(6,66,2,135), nrow = 2))
phyper(6,60,135-60,8,lower.tail = F,log.p = F)
0.01423263


lineage_eye <- list(
  aHBG=name2$V1,
  hyper=c("Bacillus_toyonensis_zheng557","Bacillus_toyonensis_zheng7","Bacillus_toyonensis_zheng178",
          "Y16G1233DL_spades","Y16G1234DL_spades","Y16G1230DL_spades","Y15G1082DL_spades","Y15G1084DL_spades"))
names(lineage_eye)<-c("","")
venn.diagram(lineage_eye,filename="lineage_eye.pdf",height=1000,width=1000,fill=mycolor[1:2],cex=0.8,fontface=2,lty=0)
phyper(6,60,135-60,8,lower.tail = F,log.p = F)
0.01423263
phyper(10,37,135-37,15,lower.tail = F,log.p = F)
0.01423263
#VF----
sub_content <- read.table("gene_presence_absence_23.Rtab",header = T)
rownames(sub_content) <- sub_content$Gene
sub_content <- sub_content[,2:24]
sub_content["plcA","GCA_002811445.1_ASM281144v1_genomic"] <- 0
sub_content["plcA","GCA_001909175.1_ASM190917v1_genomic"] <- 0
VF <- read.table("VF",header = F,sep="\t")
sub_content_vf <- sub_content[which(rownames(sub_content) %in% VF$V2),]
sub_content_vf <- sub_content_vf[order(match(rownames(sub_content_vf), VF$V2)),]
ocular <-
  c(
    "Bacillus_toyonensis_zheng557",
    "Bacillus_toyonensis_zheng7",
    "Bacillus_toyonensis_zheng178",
    "Y16G1234DL_spades",
    "Y16G1230DL_spades",
    "Y16G1233DL_spades",
    "Y15G1082DL_spades",
    "Y15G1084DL_spades"
  )
stomach <- c(name[grep("GCA", name$V1), "V1"],"Y16G1232DL_spades","Y15G1085DL_spades")
sub_content_vf_ocular <- sub_content_vf[,which(colnames(sub_content_vf) %in% ocular)]
sub_content_vf_ocular$sum <- apply(sub_content_vf_ocular,1,sum)
sub_content_vf_ocular$pro <- sub_content_vf_ocular$sum/8
sub_content_vf_stomach <- sub_content_vf[,which(colnames(sub_content_vf) %in% stomach)]
sub_content_vf_stomach$sum <- apply(sub_content_vf_stomach,1,sum)
sub_content_vf_stomach$pro <- sub_content_vf_stomach$sum/15
a <- rbind(sub_content_vf_ocular$pro,sub_content_vf_stomach$pro)
barplot(a,beside=T)
barplot(
  a,
  ylab = c(""),
  beside=T,
  col = mycolor[1:2],
  border = mycolor[1:2],
  cex.lab = 1.6,
  font.lab = 2,
  axes = F
)
axis(
  2,
  lwd = 3,
  lwd.ticks = 3,
  cex.axis = 1.5,
  font.axis = 2
)
fisher.test(matrix(c(8,8-8,6,15-6), nrow = 2),alternative = "g")
fisher.test(matrix(c(5,8-5,2,15-2), nrow = 2),alternative = "g")
fisher.test(matrix(c(7,8-7,7,15-7), nrow = 2),alternative = "g")
fisher.test(matrix(c(7,8-7,7,15-7), nrow = 2),alternative = "g")
fisher.test(matrix(c(7,8-7,6,15-6), nrow = 2),alternative = "g")

#compaired ocular and intestine(gene and mutation)----
ocular <-
  c(
    "Bacillus_toyonensis_zheng557",
    "Bacillus_toyonensis_zheng7",
    "Bacillus_toyonensis_zheng178",
    "Y16G1234DL_spades",
    "Y16G1230DL_spades",
    "Y16G1233DL_spades",
    "Y15G1082DL_spades",
    "Y15G1084DL_spades"
  )
stomach <- c(name[grep("GCA", name$V1), "V1"],"Y16G1232DL_spades","Y15G1085DL_spades")

sub_content <- read.table("gene_presence_absence_23.Rtab",header = T)
rownames(sub_content) <- sub_content$Gene
sub_content <- sub_content[,2:24]
sub_content["plcA","GCA_002811445.1_ASM281144v1_genomic"] <- 0
sub_content["plcA","GCA_001909175.1_ASM190917v1_genomic"] <- 0
remove <- read.table("remove_name")
sub_content <- sub_content[which(rownames(sub_content) %in% setdiff(rownames(sub_content),remove$V1)),]
sub_content <-
  sub_content[which(rownames(sub_content) %in% setdiff(
    rownames(sub_content),
    rownames(sub_content[grep("group", rownames(sub_content)), ])
  )), ]
sub_content_n_c1 <- sub_content[,which(colnames(sub_content) %in% ocular)]
sub_content_n_c1$count <- apply(sub_content_n_c1,1,sum)
sub_content_n_c2 <- sub_content[,which(colnames(sub_content) %in% stomach)]
sub_content_n_c2$count <- apply(sub_content_n_c2,1,sum)

sub_content_n_combine <- as.data.frame(cbind(sub_content_n_c1$count,sub_content_n_c2$count))
colnames(sub_content_n_combine) <- c("cluster1","cluster2")
rownames(sub_content_n_combine) <- rownames(sub_content_n_c1)
estimate <- function(x,c1,c2){fisher <-
  fisher.test(matrix(c(
    as.numeric(x[c1]), 8 - as.numeric(x[c1]), as.numeric(x[c2]), 15 - as.numeric(x[c2])
  ), nrow = 2),alternative = "g");fisher$estimate}
sub_content_n_combine$odds1 <- apply(sub_content_n_combine,1,estimate,c1="cluster1",c2="cluster2")
pvalue <- function(x,c1,c2){fisher <-
  fisher.test(matrix(c(
    as.numeric(x[c1]), 8 - as.numeric(x[c1]), as.numeric(x[c2]), 15 - as.numeric(x[c2])
  ), nrow = 2),alternative = "g");fisher$p.value}
sub_content_n_combine$P1 <- apply(sub_content_n_combine,1,pvalue,c1="cluster1",c2="cluster2")
sub_content_n_combine_VF <- sub_content_n_combine[which(rownames(sub_content_n_combine)%in%VF$V2),]
sub_content_n_combine1 <- sub_content_n_combine[which(sub_content_n_combine$P1<=0.05),]
sub_content_n_combine1 <- sub_content_n_combine1[order(sub_content_n_combine1$P1),]
sub_content_n_combine1 <- sub_content_n_combine1[which(rownames(sub_content_n_combine1) %in% setdiff(rownames(sub_content_n_combine1),rownames(sub_content_n_combine1[grep("group",rownames(sub_content_n_combine1)),]))),]

#PGWAS----
PGWAS_ocularcluster <- read.table("PGWAS.csv",header = T,sep=",")

#selection analysis----
omega <- read.table("omega",header = F,sep="\t",col.names = c("gene","omega","dn","ds"))
omega <- omega[order(omega$omega,decreasing = T),]
dir(path = "/home/yuanjian/Project/Bacillus_Z./",
    pattern = "*distance.new",
    full.names = T) -> similarity
dir(path = "/home/yuanjian/Project/Bacillus_Z./",
    pattern = "*distance.new",
    full.names = F) -> similarity_name
lapply(
  similarity,
  read.table,
  header = F,
  quote = "",
  sep = "\t",
  as.is = T,
  col.names = c("pdistance")
) -> similarity_all
par(mar=c(3,3,1,1),mfrow = c(2, 1),lwd=3)
plot(omega$omega,ylim=c(0,0.2),xlab=c(""),ylab=c(""),type = "b",pch=16,col = mycolor[1], axes=F)
axis(2,lwd=3,lwd.ticks = 3,cex.axis=1.4,font.axis=2)
axis(1,at = c(1:27),lwd=3,lwd.ticks = 3,cex.axis=1.4,font.axis=2)
par(new=T)
plot(omega$ds,xlab=c(""),ylab=c(""),type = "b",pch=16,col = mycolor[2], axes=F)
axis(2,lwd=3,lwd.ticks = 3,cex.axis=1.4,font.axis=2)
axis(1,at = c(1:27),lwd=3,lwd.ticks = 3,cex.axis=1.4,font.axis=2)
par(new=T)
plot(omega$dn,xlab=c(""),ylab=c(""),type = "b",pch=16,col = mycolor[3], axes=F)
axis(2,lwd=3,lwd.ticks = 3,cex.axis=1.4,font.axis=2)
axis(1,at = c(1:27),lwd=3,lwd.ticks = 3,cex.axis=1.4,font.axis=2)

all_omega_gene <- read.table("all_omega_gene2",header = F,sep="\t",col.names = c("gene","omega"))
all_omega_gene <- all_omega_gene[1:20,]
plot(all_omega_gene$omega,ylim=c(0,0.5),xlab=c(""),ylab=c(""),type = "b",pch=16,col = mycolor[2], axes=F)
axis(2,lwd=3,lwd.ticks = 3,cex.axis=1.4,font.axis=2)
axis(1,at = c(1:20),lwd=3,lwd.ticks = 3,cex.axis=1.4,font.axis=2)

#modelling----
model <- read.table("basic-example/1gym.profile.txt",sep="\t")
template1 <- read.table("basic-example/ECHKIKFN_01472.B99990004.profile.txt",sep="\t")
template2 <- read.table("advanced-example/multiple_template/ECHKIKFN_01472.B99990002.profile.txt",sep="\t")
plot(model$V3,xlim=c(0,300),ylim=c(-0.06,-0.01),col="red",type="l")
par(new=T)
plot(template1$V3,xlim=c(0,300),ylim=c(-0.06,-0.01),col="blue",type="l")
par(new=T)
plot(template2$V3,xlim=c(0,300),ylim=c(-0.06,-0.01),col="green",type="l")

model <- read.table("basic-example/1aod.profile.txt",sep="\t")
template1 <- read.table("basic-example/ECHKIKFN_00012.B99990004.profile.txt",sep="\t")
template2 <- read.table("advanced-example/multiple_template/ECHKIKFN_00012.B99990002.profile.txt",sep="\t")
plot(model$V3,xlim=c(0,500),ylim=c(-0.06,0),col="red",type="l")
par(new=T)
plot(template1$V3,xlim=c(0,500),ylim=c(-0.06,0),col="blue",type="l")
par(new=T)
plot(template2$V3,xlim=c(0,500),ylim=c(-0.06,0),col="green",type="l")

#ERG----
ERG_A <- read.table("Bacillus/ERG_A.csv",header = F)
rownames(ERG_A) <- ERG_A$V1
ERG_A_zheng557 <- ERG_A[,2:4]
ERG_A_zheng178 <- ERG_A[,5:7]
ERG_A_zheng7 <- ERG_A[,8:10]
ERG_A_zheng557$mean <- apply(ERG_A_zheng557[,1:3],1,mean)
ERG_A_zheng557$sd <- apply(ERG_A_zheng557[,1:3],1,sd)
ERG_A_zheng557$se <- ERG_A_zheng557$sd/sqrt(3)
ERG_A_zheng178$mean <- apply(ERG_A_zheng178[,1:3],1,mean)
ERG_A_zheng178$sd <- apply(ERG_A_zheng178[,1:3],1,sd)
ERG_A_zheng178$se <- ERG_A_zheng178$sd/sqrt(3)
ERG_A_zheng7$mean <- apply(ERG_A_zheng7[,1:3],1,mean)
ERG_A_zheng7$sd <- apply(ERG_A_zheng7[,1:3],1,sd)
ERG_A_zheng7$se <- ERG_A_zheng7$sd/sqrt(3)
ERG_B <- read.table("Bacillus/ERG_B.csv",header = F)
rownames(ERG_B) <- ERG_B$V1
ERG_B_zheng557 <- ERG_B[,2:4]
ERG_B_zheng178 <- ERG_B[,5:7]
ERG_B_zheng7 <- ERG_B[,8:10]
ERG_B_zheng557$mean <- apply(ERG_B_zheng557[,1:3],1,mean)
ERG_B_zheng557$sd <- apply(ERG_B_zheng557[,1:3],1,sd)
ERG_B_zheng557$se <- ERG_B_zheng557$sd/sqrt(3)
ERG_B_zheng178$mean <- apply(ERG_B_zheng178[,1:3],1,mean)
ERG_B_zheng178$sd <- apply(ERG_B_zheng178[,1:3],1,sd)
ERG_B_zheng178$se <- ERG_B_zheng178$sd/sqrt(3)
ERG_B_zheng7$mean <- apply(ERG_B_zheng7[,1:3],1,mean)
ERG_B_zheng7$sd <- apply(ERG_B_zheng7[,1:3],1,sd)
ERG_B_zheng7$se <- ERG_B_zheng7$sd/sqrt(3)


par(mfrow = c(1, 2),mar=c(3,5,1,1),lwd=3)
plot(
  c(1:5),
  ERG_A_zheng557$mean,
  ylim=c(0,100),
  xlab="",
  ylab="",
  pch = 19,
  cex=1.5,
  type="b",
  col = mycolor[1],
  axes = F
)
arrows(
  c(1:5),
  ERG_A_zheng557$mean - ERG_A_zheng557$se,
  c(1:5),
  ERG_A_zheng557$mean + ERG_A_zheng557$se,
  length = 0.05,
  angle = 90,
  code = 3,
  col = mycolor[1]
)
par(new=T)
plot(
  c(1:5),
  ERG_A_zheng178$mean,
  ylim=c(0,100),
  xlab="",
  ylab="",
  pch = 19,
  cex=1.5,
  type="b",
  col = mycolor[2],
  axes = F
)
arrows(
  c(1:5),
  ERG_A_zheng178$mean - ERG_A_zheng178$se,
  c(1:5),
  ERG_A_zheng178$mean + ERG_A_zheng178$se,
  length = 0.05,
  angle = 90,
  code = 3,
  col = mycolor[2]
)
par(new=T)
plot(
  c(1:5),
  ERG_A_zheng7$mean,
  ylim=c(0,100),
  xlab="",
  ylab="",
  pch = 19,
  cex=1.5,
  type="b",
  col = mycolor[3],
  axes = F
)
arrows(
  c(1:5),
  ERG_A_zheng7$mean - ERG_A_zheng7$se,
  c(1:5),
  ERG_A_zheng7$mean + ERG_A_zheng7$se,
  length = 0.05,
  angle = 90,
  code = 3,
  col = mycolor[3]
)
axis(2,lwd=3,lwd.ticks = 3,cex.axis=1,font.axis=2)
axis(1,lwd=3,lwd.ticks = 3,cex.axis=1,font.axis=2)
legend(
  2,100,
  legend = c(
    "",
    "",
    ""),
  bty = "n",
  cex=1.5,
  lty=1,
  col = c(mycolor[1:3]),
  text.font = 2
)



plot(
  c(1:5),
  ERG_B_zheng557$mean,
  ylim=c(0,100),
  xlab="",
  ylab="",
  pch = 19,
  cex=1.5,
  type="b",
  col = mycolor[1],
  axes = F
)
arrows(
  c(1:5),
  ERG_B_zheng557$mean - ERG_B_zheng557$se,
  c(1:5),
  ERG_B_zheng557$mean + ERG_B_zheng557$se,
  length = 0.05,
  angle = 90,
  code = 3,
  col = mycolor[1]
)
par(new=T)
plot(
  c(1:5),
  ERG_B_zheng178$mean,
  ylim=c(0,100),
  xlab="",
  ylab="",
  pch = 19,
  cex=1.5,
  type="b",
  col = mycolor[2],
  axes = F
)
arrows(
  c(1:5),
  ERG_B_zheng178$mean - ERG_B_zheng178$se,
  c(1:5),
  ERG_B_zheng178$mean + ERG_B_zheng178$se,
  length = 0.05,
  angle = 90,
  code = 3,
  col = mycolor[2]
)
par(new=T)
plot(
  c(1:5),
  ERG_B_zheng7$mean,
  ylim=c(0,100),
  xlab="",
  ylab="",
  pch = 19,
  cex=1.5,
  type="b",
  col = mycolor[3],
  axes = F
)
arrows(
  c(1:5),
  ERG_B_zheng7$mean - ERG_B_zheng7$se,
  c(1:5),
  ERG_B_zheng7$mean + ERG_B_zheng7$se,
  length = 0.05,
  angle = 90,
  code = 3,
  col = mycolor[3]
)
axis(2,lwd=3,lwd.ticks = 3,cex.axis=1,font.axis=2)
axis(1,lwd=3,lwd.ticks = 3,cex.axis=1,font.axis=2)
legend(
  2,100,
  legend = c(
    "",
    "",
    ""),
  bty = "n",
  lty=1,
  cex=1.5,
  col = c(mycolor[1:3]),
  text.font = 2
)

a0 <- c(0,100.00,100.00,100.00,100.00,100.00,100.00,100.00,100.00,100.00,100,100,100)
a6 <- c(6,77.38,67.43,86.9,69.38,57.43,46.90,46.23,55.27,35.45,60.77,58.56,55.39)
a <- as.data.frame(rbind(c(100.00,100.00,100.00),c(77.38,80.0,86.9),
           c(69.38,57.43,46.90),c(46.23,55.27,35.45),
           c(60.77,58.56,55.39)))
b0 <- c(0,100.00,100.00,100.00,100.00,100.00,100.00,100.00,100.00,100.00,100,100,100)
b6 <- c(6,71.54,83.35,80.55,41.54,53.35,40.55,37.08,43.92,30.55,57.30,54.33,49.57)
b <- as.data.frame(rbind(c(100.00,100.00,100.00),c(71.54,83.35,80.55),
           c(41.54,53.35,40.55),c(37.08,43.92,30.55),
           c(57.30,54.33,49.57)))

a$mean <- apply(a[,1:3],1,mean)
a$sd <- apply(a[,1:3],1,sd)
a$se <- a$sd/sqrt(3)
b$mean <- apply(b[,1:3],1,mean)
b$sd <- apply(b[,1:3],1,sd)
b$se <- b$sd/sqrt(3)
par(mfrow = c(1, 2),mar=c(3,5,1,1),lwd=3)
x <- barplot(
  a$mean,
  beside = T,
  ylim = c(0,100),
  col = c("grey",mycolor[1:4]),
  cex.lab = 1.6,
  font.lab = 2,
  border = NA,
  axes = F
)
arrows(
  x,
  a$mean - a$se,
  x,
  a$mean + a$se,
  length = 0.05,
  angle = 90,
  code = 3,
  col = c("grey",mycolor[1:4])
)
axis(2,lwd=3,lwd.ticks = 3,cex.axis=1.4,font.axis=2)


x <- barplot(
  b$mean,
  beside = T,
  ylim = c(0,100),
  col = c("grey",mycolor[1:4]),
  cex.lab = 1.6,
  font.lab = 2,
  border = NA,
  axes = F
)
arrows(
  x,
  b$mean - b$se,
  x,
  b$mean + b$se,
  length = 0.05,
  angle = 90,
  code = 3,
  col = c("grey",mycolor[1:4])
)
axis(2,lwd=3,lwd.ticks = 3,cex.axis=1.4,font.axis=2)


t.test(c(100.00,100.00,100.00),c(77.38,80.0,86.9))
t.test(c(77.38,80.0,86.9),c(69.38,57.43,46.90,46.23,55.27,35.45,60.77,58.56,55.39))

t.test(c(100.00,100.00,100.00),c(71.54,83.35,80.55))
t.test(c(77.38,80.0,86.9),c(41.54,53.35,40.55,37.08,43.92,30.55,57.30,54.33,49.57))
#CK----
IL_1β <- read.table("IL-1β.csv",header = F)
rownames(IL_1β) <- IL_1β$V1
IL_1β_PBS <- IL_1β[,2:4]
IL_1β_zheng557 <- IL_1β[,5:7]
IL_1β_zheng178 <- IL_1β[,8:10]
IL_1β_zheng7 <- IL_1β[,11:13]
IL_1β_PBS$mean <- apply(IL_1β_PBS[,1:3],1,mean)
IL_1β_PBS$sd <- apply(IL_1β_PBS[,1:3],1,sd)
IL_1β_PBS$se <- IL_1β_PBS$sd/sqrt(3)
IL_1β_zheng557$mean <- apply(IL_1β_zheng557[,1:3],1,mean)
IL_1β_zheng557$sd <- apply(IL_1β_zheng557[,1:3],1,sd)
IL_1β_zheng557$se <- IL_1β_zheng557$sd/sqrt(3)
IL_1β_zheng178$mean <- apply(IL_1β_zheng178[,1:3],1,mean)
IL_1β_zheng178$sd <- apply(IL_1β_zheng178[,1:3],1,sd)
IL_1β_zheng178$se <- IL_1β_zheng178$sd/sqrt(3)
IL_1β_zheng7$mean <- apply(IL_1β_zheng7[,1:3],1,mean)
IL_1β_zheng7$sd <- apply(IL_1β_zheng7[,1:3],1,sd)
IL_1β_zheng7$se <- IL_1β_zheng7$sd/sqrt(3)
IL_1β_mean <- rbind(IL_1β_PBS$mean,IL_1β_zheng557$mean,IL_1β_zheng178$mean,IL_1β_zheng7$mean)
IL_1β_se <- rbind(IL_1β_PBS$se,IL_1β_zheng557$se,IL_1β_zheng178$se,IL_1β_zheng7$se)

IL_6 <- read.table("IL-6.csv",header = F)
rownames(IL_6) <- IL_6$V1
IL_6_PBS <- IL_6[,2:4]
IL_6_zheng557 <- IL_6[,5:7]
IL_6_zheng178 <- IL_6[,8:10]
IL_6_zheng7 <- IL_6[,11:13]
IL_6_PBS$mean <- apply(IL_6_PBS[,c(1,3)],1,mean)
IL_6_PBS$sd <- apply(IL_6_PBS[,c(1,3)],1,sd)
IL_6_PBS$se <- IL_6_PBS$sd/sqrt(3)
IL_6_zheng557$mean <- apply(IL_6_zheng557[,1:3],1,mean)
IL_6_zheng557$sd <- apply(IL_6_zheng557[,1:3],1,sd)
IL_6_zheng557$se <- IL_6_zheng557$sd/sqrt(3)
IL_6_zheng178$mean <- apply(IL_6_zheng178[,1:3],1,mean)
IL_6_zheng178$sd <- apply(IL_6_zheng178[,1:3],1,sd)
IL_6_zheng178$se <- IL_6_zheng178$sd/sqrt(3)
IL_6_zheng7$mean <- apply(IL_6_zheng7[,1:3],1,mean)
IL_6_zheng7$sd <- apply(IL_6_zheng7[,1:3],1,sd)
IL_6_zheng7$se <- IL_6_zheng7$sd/sqrt(3)
IL_6_mean <- rbind(IL_6_PBS$mean,IL_6_zheng557$mean,IL_6_zheng178$mean,IL_6_zheng7$mean)
IL_6_se <- rbind(IL_6_PBS$se,IL_6_zheng557$se,IL_6_zheng178$se,IL_6_zheng7$se)

TNF <- read.table("TNF-α.csv",header = F)
rownames(TNF) <- TNF$V1
TNF_PBS <- TNF[,2:4]
TNF_zheng557 <- TNF[,5:7]
TNF_zheng178 <- TNF[,8:10]
TNF_zheng7 <- TNF[,11:13]
TNF_PBS$mean <- apply(TNF_PBS[,1:3],1,mean)
TNF_PBS$sd <- apply(TNF_PBS[,1:3],1,sd)
TNF_PBS$se <- TNF_PBS$sd/sqrt(3)
TNF_zheng557$mean <- apply(TNF_zheng557[,1:3],1,mean)
TNF_zheng557$sd <- apply(TNF_zheng557[,1:3],1,sd)
TNF_zheng557$se <- TNF_zheng557$sd/sqrt(3)
TNF_zheng178$mean <- apply(TNF_zheng178[,1:3],1,mean)
TNF_zheng178$sd <- apply(TNF_zheng178[,1:3],1,sd)
TNF_zheng178$se <- TNF_zheng178$sd/sqrt(3)
TNF_zheng7$mean <- apply(TNF_zheng7[,1:3],1,mean)
TNF_zheng7$sd <- apply(TNF_zheng7[,1:3],1,sd)
TNF_zheng7$se <- TNF_zheng7$sd/sqrt(3)
TNF_mean <- rbind(TNF_PBS$mean,TNF_zheng557$mean,TNF_zheng178$mean,TNF_zheng7$mean)
TNF_se <- rbind(TNF_PBS$se,TNF_zheng557$se,TNF_zheng178$se,TNF_zheng7$se)

par(mfrow = c(3, 1),mar=c(3,5,1,1),lwd=3)
x <- barplot(
  IL_1β_mean,
  beside = T,
  ylim = c(0,500),
  col = c("grey",mycolor[1:3]),
  cex.lab = 1.6,
  font.lab = 2,
  border = NA,
  axes = F
)
arrows(
  x,
  IL_1β_mean - IL_1β_se,
  x,
  IL_1β_mean + IL_1β_se,
  length = 0.05,
  angle = 90,
  code = 3,
  col = c("grey",mycolor[1:3])
)
axis(2,lwd=3,lwd.ticks = 3,cex.axis=1.4,font.axis=2)
legend(
  0,580,
  cex = 2,
  legend = c("", "", "",""),
  bty = "n",
  fill = c("grey",mycolor[1:3]),
  border = NA,
  text.font = 2
)

x <- barplot(
  IL_6_mean,
  beside = T,
  ylim = c(0,560),
  col = c("grey",mycolor[1:3]),
  cex.lab = 1.6,
  font.lab = 2,
  border = NA,
  axes = F
)
arrows(
  x,
  IL_6_mean - IL_6_se,
  x,
  IL_6_mean + IL_6_se,
  length = 0.05,
  angle = 90,
  code = 3,
  col = c("grey",mycolor[1:3])
)
axis(2,lwd=3,lwd.ticks = 3,cex.axis=1.4,font.axis=2)

x <- barplot(
  TNF_mean,
  beside = T,
  ylim = c(0,500),
  col = c("grey",mycolor[1:3]),
  cex.lab = 1.6,
  font.lab = 2,
  border = NA,
  axes = F
)
arrows(
  x,
  TNF_mean - TNF_se,
  x,
  TNF_mean + TNF_se,
  length = 0.05,
  angle = 90,
  code = 3,
  col = c("grey",mycolor[1:3])
)
axis(2,lwd=3,lwd.ticks = 3,cex.axis=1.4,font.axis=2)


