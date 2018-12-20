### Create RSS graphs for Stirparo staging

#set working directory
setwd("/mnt/nfs/data/Vincent_Lab/Joel_Chappell/SCENIC_petropoulos/RSS")

# read in data
EPI<-read.table("RSS.EPI.txt", sep = "\t")
ICM <- read.table("RSS.ICM.txt", sep = "\t")
inter <- read.table("RSS.intermediate.txt", sep = "\t")
PE <- read.table("RSS.prE.txt", sep = "\t")
TE <- read.table("RSS.TE.txt", sep = "\t")

#rename columns
colnames(EPI)[1] <- "Regulons"
colnames(EPI)[2] <- "Specificity_score"
colnames(ICM)[1] <- "Regulons"
colnames(ICM)[2] <- "Specificity_score"
colnames(inter)[1] <- "Regulons"
colnames(inter)[2] <- "Specificity_score"
colnames(PE)[1] <- "Regulons"
colnames(PE)[2] <- "Specificity_score"
colnames(TE)[1] <- "Regulons"
colnames(TE)[2] <- "Specificity_score"

#load dependancies 
library(ggplot2)
library(ggrepel)

#make EPI graph
#adapted to include top 10 + selected regulons
A<- c("POU5F1 (64g)", "NANOG (550g)", "SOX2 (284g)")

ggplot(EPI, aes(Regulons, Specificity_score)) +
  ggtitle("Epiblast") +
  geom_point() +
  geom_point(color = ifelse(EPI$Specificity_score > 0.32955, "red", "grey"), size = 1) +
  geom_point(data = subset(EPI, Regulons == "VENTX_extended (58g)"), aes(Regulons, Specificity_score), colour = "blue", size = 1) +
  geom_point(data = subset(EPI, Regulons == "ETV4_extended (303g)"), aes(Regulons, Specificity_score), colour = "blue", size = 1) +
  geom_point(data = subset(EPI, Regulons == "POU5F1 (64g)"), aes(Regulons, Specificity_score), colour = "blue", size = 1) +
  geom_point(data = subset(EPI, Regulons == "NANOG (550g)"), aes(Regulons, Specificity_score), colour = "blue", size = 1) +
  geom_point(data = subset(EPI, Regulons == "SOX2 (284g)"), aes(Regulons, Specificity_score), colour = "blue", size = 1) +
  geom_text_repel(data=filter(EPI, Specificity_score>0.32955), aes(label=Regulons), size = 2, xlim = c(15,NA), segment.size =0.1) +
  geom_text_repel(data=filter(EPI, Regulons %in% A),  aes(label=Regulons), size = 2, xlim = c(30,NA), segment.size =0.1) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "cm"),
        plot.title = element_text(hjust=0.5, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_x_discrete(limits= EPI$Regulon,
                   breaks = c("ZSCAN10 (164g)", "FOSL1_extended (24g)", "TBP (703g)"),
                   labels = c("ZSCAN10 (164g)"="1", "FOSL1_extended (24g)"="191", "TBP (703g)" ="382"),
                   expand = c(0,5)) +
  scale_y_continuous(limits = c(0.16,0.43),
                     breaks = seq(0.16, 0.43, 0.02))

#make ICM graph
ggplot(ICM, aes(Regulons, Specificity_score)) +
  ggtitle("Inner cell mass") +
  geom_point() +
  geom_point(color = ifelse(ICM$Specificity_score > 0.35, "red", "grey"), size = 1) +
  geom_text_repel(data=filter(ICM, Specificity_score>0.35), aes(label=Regulons), size = 2, xlim = c(15,NA)) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "cm"),
        plot.title = element_text(hjust=0.5, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_x_discrete(limits= ICM$Regulon,
                   breaks = c("PBX3_extended (35g)", "RAX_extended (44g)", "TCF7L2_extended (384g)"),
                   labels = c("PBX3_extended (35g)"="1", "RAX_extended (44g)"="191", "TCF7L2_extended (384g)" ="382"),
                   expand = c(0,5)) +
  scale_y_continuous(limits = c(0.15,0.45),
                     breaks = seq(0.15, 0.45, 0.05))

#make inter graph
ggplot(inter, aes(Regulons, Specificity_score)) +
  ggtitle("Intermediate") +
  geom_point() +
  geom_point(color = ifelse(inter$Specificity_score > 0.25, "red", "grey"), size = 1) +
  geom_text_repel(data=filter(inter, Specificity_score>0.25), aes(label=Regulons), size = 2, xlim = c(15,NA)) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "cm"),
        plot.title = element_text(hjust=0.5, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_x_discrete(limits= inter$Regulon,
                   breaks = c("JUNB (35g)", "STAT5B_extended (53g)", "TBP (703g)"),
                   labels = c(" JUNB (35g)"="1", "STAT5B_extended (53g)"="191", "TBP (703g)" ="382"),
                   expand = c(0,5)) +
  scale_y_continuous(limits = c(0.15,0.3),
                     breaks = seq(0.15, 0.3, 0.05))

#make PE graph
#adapted to include top 10 + selected regulons
ggplot(PE, aes(Regulons, Specificity_score)) +
  ggtitle("Primitive endoderm") +
  geom_point() +
  geom_point(color = ifelse(PE$Specificity_score > 0.29727, "red", "grey"), size = 1) +
  geom_point(data = subset(PE, Regulons == "FOXA2 (31g)"), aes(Regulons, Specificity_score), colour = "blue", size = 1) +
  geom_point(data = subset(PE, Regulons == "GATA4 (143g)"), aes(Regulons, Specificity_score), colour = "blue", size = 1) +
  geom_point(data = subset(PE, Regulons == "SOX17 (24g)"), aes(Regulons, Specificity_score), colour = "blue", size = 1) +
  geom_text_repel(data=filter(PE, Specificity_score>0.29727), aes(label=Regulons), size = 2, xlim = c(15,NA), segment.size =0.1) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "cm"),
        plot.title = element_text(hjust=0.5, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_x_discrete(limits= PE$Regulon,
                   breaks = c("ISL1_extended (17g)", "TAF1 (13g)", "TBP (703g)"),
                   labels = c("ISL1_extended (17g)"="1", "TAF1 (13g)"="191", "TBP (703g)" ="382"),
                   expand = c(0,5)) +
  scale_y_continuous(limits = c(0.16,0.52),
                     breaks = seq(0.16, 0.52, 0.02))

#make TE graph
#adapted to include top 10 + selected regulons
ggplot(TE, aes(Regulons, Specificity_score)) +
  ggtitle("Trophectoderm") +
  geom_point() +
  geom_point(color = ifelse(TE$Specificity_score > 0.67304, "red", "grey"), size = 1) +
  geom_point(data = subset(TE, Regulons == "GATA2 (825g)"), aes(Regulons, Specificity_score), colour = "blue", size = 1) +
  geom_point(data = subset(TE, Regulons == "GATA3 (243g)"), aes(Regulons, Specificity_score), colour = "blue", size = 1) +
  geom_text_repel(data=filter(TE, Specificity_score>0.67304), aes(label=Regulons), size = 2, xlim = c(30,NA), segment.size =0.1) +
  geom_text_repel(data=filter(TE, Regulons == "GATA2 (825g)"),  aes(label=Regulons), size = 2, xlim = c(40, 40), ylim = c(0.61,0.61), segment.size =0.1) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "cm"),
        plot.title = element_text(hjust=0.5, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_x_discrete(limits= TE$Regulon,
                   breaks = c("GRHL2_extended (61g)", "MIXL1_extended (138g)", "THAP1 (27g)"),
                   labels = c("GRHL2_extended (61g)"="1", "MIXL1_extended (138g)"="191", "THAP1 (27g)" ="382"),
                   expand = c(0,5)) +
  scale_y_continuous(limits = c(0.16,0.75),
                     breaks = seq(0.16, 0.75, 0.05))

