### Create RSS graphs for Laurent David staging

#set working directory
setwd("/mnt/nfs/data/Vincent_Lab/Joel_Chappell/SCENIC_petropoulos/RSS")

# read in data
pre_morula<-read.table("RSS.1.Pre-morula.txt", sep = "\t")
morula <- read.table("RSS.2.Morula.txt", sep = "\t")
blast <- read.table("RSS.3.Early blastocyst.txt", sep = "\t")
ICM <- read.table("RSS.4.Inner cell mass.txt", sep = "\t")
E_TE <- read.table("RSS.5.Early trophectoderm.txt", sep = "\t")
Epi <- read.table("RSS.6.Epiblast.txt", sep = "\t")
PE <- read.table("RSS.7.Primitive endoderm.txt", sep = "\t")
TE_Nneg <- read.table("RSS.8.TE.NR2F2-.txt", sep = "\t")
TE_Npos <- read.table("RSS.9.TE.NR2F2+.txt", sep = "\t")

#rename columns
colnames(pre_morula)[1] <- "Regulons"
colnames(pre_morula)[2] <- "Specificity_score"
colnames(morula)[1] <- "Regulons"
colnames(morula)[2] <- "Specificity_score"
colnames(blast)[1] <- "Regulons"
colnames(blast)[2] <- "Specificity_score"
colnames(ICM)[1] <- "Regulons"
colnames(ICM)[2] <- "Specificity_score"
colnames(E_TE)[1] <- "Regulons"
colnames(E_TE)[2] <- "Specificity_score"
colnames(Epi)[1] <- "Regulons"
colnames(Epi)[2] <- "Specificity_score"
colnames(PE)[1] <- "Regulons"
colnames(PE)[2] <- "Specificity_score"
colnames(TE_Nneg)[1] <- "Regulons"
colnames(TE_Nneg)[2] <- "Specificity_score"
colnames(TE_Npos)[1] <- "Regulons"
colnames(TE_Npos)[2] <- "Specificity_score"

#load dependancies 
library(ggplot2)
library(ggrepel)

#make pre_morula graph
ggplot(pre_morula, aes(Regulons, Specificity_score)) +
  ggtitle("Pre-morula") +
  geom_point() +
  geom_point(color = ifelse(pre_morula$Specificity_score > 0.65, "red", "grey"), size = 1) +
  geom_text_repel(data=filter(pre_morula, Specificity_score>0.65), aes(label=Regulons), size = 2, xlim = c(15,NA)) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "cm"),
        plot.title = element_text(hjust=0.5, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_x_discrete(limits= pre_morula$Regulon,
                   breaks = c("FOXN2 (431g)", "PAX4_extended (10g)", "IRF3 (1271g)"),
                   labels = c("FOXN2 (431g)"="1", "PAX4_extended (10g)"="191", "IRF3 (1271g)" ="382"),
                   expand = c(0,5)) +
  scale_y_continuous(limits = c(0.15,0.9),
                     breaks = seq(0.15, 0.9, 0.05))
                  
#make morula graph
ggplot(morula, aes(Regulons, Specificity_score)) +
  ggtitle("Morula") +
  geom_point() +
  geom_point(color = ifelse(morula$Specificity_score > 0.5, "red", "grey"), size = 1) +
  geom_text_repel(data=filter(morula, Specificity_score>0.5), aes(label=Regulons), size = 2, xlim = c(15,NA)) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "cm"),
        plot.title = element_text(hjust=0.5, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_x_discrete(limits= morula$Regulon,
                   breaks = c("TGIF1_extended (55g)", "KLF6 (118g)", "IRF3 (1271g)"),
                   labels = c("TGIF1_extended (55g)"="1", "KLF6 (118g)"="191", "IRF3 (1271g)" ="382"),
                   expand = c(0,5)) +
  scale_y_continuous(limits = c(0.15,0.6),
                     breaks = seq(0.15, 0.6, 0.05))

#make early blastocyst graph
ggplot(blast, aes(Regulons, Specificity_score)) +
  ggtitle("Early blastocyst") +
  geom_point() +
  geom_point(color = ifelse(blast$Specificity_score > 0.41, "red", "grey"), size = 1) +
  geom_text_repel(data=filter(blast, Specificity_score>0.41), aes(label=Regulons), size = 2, xlim = c(15,NA)) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "cm"),
        plot.title = element_text(hjust=0.5, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_x_discrete(limits= blast$Regulon,
                   breaks = c("POU5F1 (64g)", "TWIST2 (54g)", "DLX3_extended (396g)"),
                   labels = c("POU5F1 (64g)"="1", "TWIST2 (54g)"="191", "DLX3_extended (396g)" ="382"),
                   expand = c(0,5)) +
  scale_y_continuous(limits = c(0.15,0.45),
                     breaks = seq(0.15, 0.45, 0.05))

#make ICM graph
ggplot(ICM, aes(Regulons, Specificity_score)) +
  ggtitle("Inner cell mass") +
  geom_point() +
  geom_point(color = ifelse(ICM$Specificity_score > 0.205, "red", "grey"), size = 1) +
  geom_text_repel(data=filter(ICM, Specificity_score>0.205), aes(label=Regulons), size = 2, xlim = c(15,NA)) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "cm"),
        plot.title = element_text(hjust=0.5, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_x_discrete(limits= ICM$Regulon,
                   breaks = c("NFKB2 (53g)", "ZNF384 (13g)", "ZNF143 (512g)"),
                   labels = c("NFKB2 (53g)"="1", "ZNF384 (13g)"="191", "ZNF143 (512g)" ="382"),
                   expand = c(0,5)) +
  scale_y_continuous(limits = c(0.16,0.22),
                     breaks = seq(0.16, 0.22, 0.02))

#make Early TE graph
ggplot(E_TE, aes(Regulons, Specificity_score)) +
  ggtitle("Ealry trophectoderm") +
  geom_point() +
  geom_point(color = ifelse(E_TE$Specificity_score > 0.3, "red", "grey"), size = 1) +
  geom_text_repel(data=filter(E_TE, Specificity_score>0.3), aes(label=Regulons), size = 2, xlim = c(15,NA), segment.size =0.1) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "cm"),
        plot.title = element_text(hjust=0.5, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_x_discrete(limits= E_TE$Regulon,
                   breaks = c("ESRRB_extended (26g)", "HOXD9 (10g)", "FOSB (27g)"),
                   labels = c("ESRRB_extended (26g)"="1", "HOXD9 (10g)"="191", "FOSB (27g)" ="382"),
                   expand = c(0,5)) +
  scale_y_continuous(limits = c(0.16,0.36),
                     breaks = seq(0.16, 0.36, 0.02))
  
#make Epiblast graph
#adapted to include top 10 + selected regulons

A<- c("POU5F1 (64g)", "NANOG (550g)", "SOX2 (284g)")

ggplot(Epi, aes(Regulons, Specificity_score)) +
  ggtitle("Epiblast") +
  geom_point() +
  geom_point(color = ifelse(Epi$Specificity_score > 0.49587, "red", "grey"), size = 1) +
  geom_point(data = subset(Epi, Regulons == "VENTX_extended (58g)"), aes(Regulons, Specificity_score), colour = "blue", size = 1) +
  geom_point(data = subset(Epi, Regulons == "ETV4_extended (303g)"), aes(Regulons, Specificity_score), colour = "blue", size = 1) +
  geom_point(data = subset(Epi, Regulons == "POU5F1 (64g)"), aes(Regulons, Specificity_score), colour = "blue", size = 1) +
  geom_point(data = subset(Epi, Regulons == "NANOG (550g)"), aes(Regulons, Specificity_score), colour = "blue", size = 1) +
  geom_point(data = subset(Epi, Regulons == "SOX2 (284g)"), aes(Regulons, Specificity_score), colour = "blue", size = 1) +
  geom_text_repel(data=filter(Epi, Specificity_score>0.49587), aes(label=Regulons), size = 3, xlim = c(15,NA), segment.size =0.1) +
  geom_text_repel(data=filter(Epi, Regulons %in% A),  aes(label=Regulons), size = 3, xlim = c(30,NA), segment.size =0.1) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "cm"),
        plot.title = element_text(hjust=0.5, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_x_discrete(limits= Epi$Regulon,
                   breaks = c("ETV4_extended (303g)", "FOXP1 (44g)", "TBP (703g)"),
                   labels = c("ETV4_extended (303g)"="1", "FOXP1 (44g)"="191", "TBP (703g)" ="382"),
                   expand = c(0,5)) +
  scale_y_continuous(limits = c(0.15,0.6),
                     breaks = seq(0.15, 0.6, 0.05))

#make PE graph
#adapted to include top 10 + selected regulons

ggplot(PE, aes(Regulons, Specificity_score)) +
  ggtitle("Primitive endoderm") +
  geom_point() +
  geom_point(color = ifelse(PE$Specificity_score > 0.27367, "red", "grey"), size = 1) +
  geom_point(data = subset(PE, Regulons == "FOXA2 (31g)"), aes(Regulons, Specificity_score), colour = "blue", size = 1) +
  geom_point(data = subset(PE, Regulons == "GATA4 (143g)"), aes(Regulons, Specificity_score), colour = "blue", size = 1) +
  geom_point(data = subset(PE, Regulons == "SOX17 (24g)"), aes(Regulons, Specificity_score), colour = "blue", size = 1) +
  geom_text_repel(data=filter(PE, Specificity_score>0.27367), aes(label=Regulons), size = 3, xlim = c(15,NA), segment.size =0.1) +
  geom_text_repel(data=filter(PE, Regulons == "SOX17 (24g)"),  aes(label=Regulons), size = 3, xlim = c(65,NA), segment.size =0.1) +
    theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "cm"),
        plot.title = element_text(hjust=0.5, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_x_discrete(limits= PE$Regulon,
                   breaks = c("HMGN3_extended (37g)", "EZH2 (22g)", "MECOM_extended (27g)"),
                   labels = c("HMGN3_extended (37g)"="1", "EZH2 (22g)"="191", "MECOM_extended (27g)" ="382"),
                   expand = c(0,5)) +
  scale_y_continuous(limits = c(0.16,0.34),
                     breaks = seq(0.16, 0.34, 0.02))

#make TE.NR2F2- graph
B<- c("GATA2 (825g)", "GATA3 (243g)")

ggplot(TE_Nneg, aes(Regulons, Specificity_score)) +
  ggtitle("Trophectoderm (NR2F2_negative)") +
  geom_point() +
  geom_point(color = ifelse(TE_Nneg$Specificity_score > 0.44917, "red", "grey"), size = 1) +
  geom_point(data = subset(TE_Nneg, Regulons == "GATA2 (825g)"), aes(Regulons, Specificity_score), colour = "blue", size = 1) +
  geom_point(data = subset(TE_Nneg, Regulons == "GATA3 (243g)"), aes(Regulons, Specificity_score), colour = "blue", size = 1) +
  geom_text_repel(data=filter(TE_Nneg, Specificity_score>0.44917), aes(label=Regulons), size = 2, xlim = c(15,NA), segment.size =0.1) +
  geom_text_repel(data=filter(TE_Nneg, Regulons %in% B),  aes(label=Regulons), size = 2, xlim = c(40,NA), ylim=c(0.3,0.43), segment.size =0.1) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "cm"),
        plot.title = element_text(hjust=0.5, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_x_discrete(limits= TE_Nneg$Regulon,
                   breaks = c("TFAP2A (99g)", "TP73 (35g)", "CEBPZ_extended (474g)"),
                   labels = c("TFAP2A (99g)"="1", "TP73 (35g)"="191", "CEBPZ_extended (474g)" ="382"),
                   expand = c(0,5)) +
  scale_y_continuous(limits = c(0.15,0.5),
                     breaks = seq(0.15, 0.5, 0.05))

#make TE.NR2F2+ graph
B<- c("GATA2 (825g)", "GATA3 (243g)")

ggplot(TE_Npos, aes(Regulons, Specificity_score)) +
  ggtitle("Trophectoderm (NR2F2+)") +
  geom_point() +
  geom_point(color = ifelse(TE_Npos$Specificity_score >  0.5087, "red", "grey"), size = 1) +
  geom_point(data = subset(TE_Npos, Regulons == "GATA2 (825g)"), aes(Regulons, Specificity_score), colour = "blue", size = 1) +
  geom_point(data = subset(TE_Npos, Regulons == "GATA3 (243g)"), aes(Regulons, Specificity_score), colour = "blue", size = 1) +
  geom_text_repel(data=filter(TE_Npos, Specificity_score> 0.5087), aes(label=Regulons), size = 3, xlim = c(15,NA),segment.size =0.1) +
  geom_text_repel(data=filter(TE_Npos, Regulons %in% B),  aes(label=Regulons), size = 3, xlim = c(30,NA),segment.size =0.1) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "cm"),
        plot.title = element_text(hjust=0.5, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_x_discrete(limits= TE_Npos$Regulon,
                   breaks = c("TEAD3_extended (180g)", "HOXA11_extended (57g)", "CEBPZ_extended (474g)"),
                   labels = c("TEAD3_extended (180g)"="1", "HOXA11_extended (57g)"="191", "CEBPZ_extended (474g)" ="382"),
                   expand = c(0,5)) +
  scale_y_continuous(limits = c(0.15,0.6),
                     breaks = seq(0.15, 0.6, 0.05))



