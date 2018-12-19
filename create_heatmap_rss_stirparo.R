#set wd
setwd("/mnt/nfs/data/Vincent_Lab/Joel_Chappell/SCENIC_petropoulos/int")

#read in binary regulon file
nonDupl <- readRDS("4.2_binaryRegulonActivity_nonDupl.Rds")

#convert to data frame
nonDupl1<-as.data.frame(nonDupl)

##load cellInfo from SCENIC if needed

cellInfo[1:5,]

#subset important data based on each lineage 
#view variables under lineage
table(cellInfo$Stirparo.lineage)

sub_EPI<-subset(cellInfo, Stirparo.lineage=="EPI")
sub_ICM<-subset(cellInfo, Stirparo.lineage=="ICM")
sub_prE<-subset(cellInfo, Stirparo.lineage=="prE")
sub_TE<-subset(cellInfo, Stirparo.lineage=="TE")


#subset binary matrix on each lineage subset 
mat_EPI <- which(colnames(nonDupl1) %in% rownames(sub_EPI)) 
mat_EPI <- nonDupl1[, mat_EPI]
mat_ICM1 <- which(colnames(nonDupl1) %in% rownames(sub_ICM)) 
mat_ICM1 <- nonDupl1[, mat_ICM1]
mat_prE <- which(colnames(nonDupl1) %in% rownames(sub_prE)) 
mat_prE <- nonDupl1[, mat_prE]
mat_TE <- which(colnames(nonDupl1) %in% rownames(sub_TE)) 
mat_TE <- nonDupl1[, mat_TE]

#sum rows and divide by number of columns

sum_EPI <- (rowSums(mat_EPI))/length(colnames(mat_EPI))
sum_ICM1 <- (rowSums(mat_ICM1))/length(colnames(mat_ICM1))
sum_prE <- (rowSums(mat_prE))/length(colnames(mat_prE))
sum_TE <- (rowSums(mat_TE))/length(colnames(mat_TE))

#combine data into one column
heat_ma<-cbind(sum_EPI, sum_ICM1, sum_prE, sum_TE)

colnames(heat_ma) <- c("Epiblast", "ICM","Primitive endoderm", "TE")

####
#somthing weird happeining all 0's for HOXC11_extended (27g) which creates an NA for input... Shouldnt be there in the first place though?
#remove this row

#heat_ma<-subset(heat_ma, !rownames(heat_ma) %in% "HOXC11_extended (27g)")
heat_ma = heat_ma[ rowSums(heat_ma)!=0, ] 

#create heatmap for all
input <- t(scale(t(heat_ma)))
input<-na.omit(input)

heatmap.2(input, col=colorRampPalette(rev(brewer.pal(9, "RdBu")) )(255), Rowv=T, scale="none", trace="none", dendrogram = "none")


#subset and create heatmap for selected genes only 

#sel_genes <- c("ETV4_extended (303g)", "VENTX_extended (58g)", "POU5F1 (64g)", "SOX2 (284g)", "GATA3 (243g)", 
#               "GATA2 (825g)", "FOXA2 (31g)", "GATA4 (143g)", "SOX17 (24g)")
#heat_m_sel<-subset(heat_m, rownames(heat_m) %in% sel_genes)

#input_sel <- t(scale(t(heat_m_sel)))
#heatmap.2(input_sel, col=colorRampPalette(rev(brewer.pal(9, "RdBu")) )(255), Rowv=TRUE, scale="none", trace="none", 
          dendrogram = "both", mar=c(10,10))

#subset and create heatmap for selected (+top 5) genes 


sel_genes_1 <- c("ZSCAN10 (164g)", "POU3F1_extended (32g)", "VENTX_extended (58g)", "ARX_extended (111g)", "PBX1_extended (10g)",
                 "PBX3_extended (35g)", "GSC (12g)", "PRDM1 (13g)", "ETV4_extended (303g)", "ISL1_extended (17g)", "GATA5 (20g)",
                 "HNF1B (178g)", "FOXA2 (31g)", "MAF (41g)", "GRHL2_extended (61g)", "JDP2_extended (204g)", "IRF3 (1271g)",
                 "CEBPA_extended (11g)", "NFE2L2_extended (214g)", "GATA2 (825g" , "GATA4 (143g)", "POU5F1 (64g)", "SOX2 (284g)",
                 "GATA3 (243g)", "SOX17 (24g)")

heat_ma_sel1<-subset(heat_ma, rownames(heat_ma) %in% sel_genes_1)

input_sel1 <- t(scale(t(heat_ma_sel1)))
heatmap.2(input_sel1, col=colorRampPalette(rev(brewer.pal(9, "RdBu")) )(255), Rowv=TRUE, Colv=F, scale="none", trace="none", 
          dendrogram = "none", mar=c(10,10))







