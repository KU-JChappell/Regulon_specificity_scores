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
table(cellInfo$Cell.state)

sub_Pm<-subset(cellInfo, Cell.state=="1.Pre-morula")
sub_m<-subset(cellInfo, Cell.state=="2.Morula")
sub_Eb<-subset(cellInfo, Cell.state=="3.Early blastocyst")
sub_ICM<-subset(cellInfo, Cell.state=="4.Inner cell mass")
sub_eTE<-subset(cellInfo, Cell.state=="5.Early trophectoderm")
sub_epi<-subset(cellInfo, Cell.state=="6.Epiblast")
sub_PE<-subset(cellInfo, Cell.state=="7.Primitive endoderm")
sub_TEn<-subset(cellInfo, Cell.state=="8.TE.NR2F2-")
sub_TEp<-subset(cellInfo, Cell.state=="9.TE.NR2F2+")

#subset binary matrix on each lineage subset 
mat_Pm <- which(colnames(nonDupl1) %in% rownames(sub_Pm)) 
mat_Pm <- nonDupl1[, mat_Pm]
mat_m <- which(colnames(nonDupl1) %in% rownames(sub_m)) 
mat_m <- nonDupl1[, mat_m]
mat_Eb <- which(colnames(nonDupl1) %in% rownames(sub_Eb)) 
mat_Eb <- nonDupl1[, mat_Eb]
mat_ICM <- which(colnames(nonDupl1) %in% rownames(sub_ICM)) 
mat_ICM <- nonDupl1[, mat_ICM]
mat_eTE <- which(colnames(nonDupl1) %in% rownames(sub_eTE)) 
mat_eTE <- nonDupl1[, mat_eTE]
mat_PE <- which(colnames(nonDupl1) %in% rownames(sub_PE)) 
mat_PE <- nonDupl1[, mat_PE]
mat_epi <- which(colnames(nonDupl1) %in% rownames(sub_epi)) 
mat_epi <- nonDupl1[, mat_epi]
mat_TEn <- which(colnames(nonDupl1) %in% rownames(sub_TEn)) 
mat_TEn <- nonDupl1[, mat_TEn]
mat_TEp <- which(colnames(nonDupl1) %in% rownames(sub_TEp)) 
mat_TEp <- nonDupl1[, mat_TEp]

#sum rows and divide by number of columns

sum_Pm <- (rowSums(mat_Pm))/length(colnames(mat_Pm))
sum_m <- (rowSums(mat_m))/length(colnames(mat_m))
sum_Eb <- (rowSums(mat_Eb))/length(colnames(mat_Eb))
sum_ICM <- (rowSums(mat_ICM))/length(colnames(mat_ICM))
sum_eTE <- (rowSums(mat_eTE))/length(colnames(mat_eTE))
sum_PE <- (rowSums(mat_PE))/length(colnames(mat_PE))
sum_epi <- (rowSums(mat_epi))/length(colnames(mat_epi))
sum_TEn <- (rowSums(mat_TEn))/length(colnames(mat_TEn))
sum_TEp <- (rowSums(mat_TEp))/length(colnames(mat_TEp))

#combine data into one column
heat_m<-cbind(sum_Pm, sum_m, sum_Eb, sum_ICM, sum_eTE, sum_epi, sum_PE, sum_TEn, sum_TEp)

colnames(heat_m) <- c("Pre-morula", "Morula","Early blastocyst", "ICM", "Early trophectoderm",  "Epiblast",
                      "Primitive endoderm", "TE.NR2F2-", "9.TE.NR2F2+")

####
#somthing weird happeining all 0's for HOXC11_extended (27g) which creates an NA for input... Shouldnt be there in the first place though?
#remove this row

heat_m<-subset(heat_m, !rownames(heat_m) %in% "HOXC11_extended (27g)")

#create heatmap for all
input <- t(scale(t(heat_m)))
heatmap.2(input, col=colorRampPalette(rev(brewer.pal(9, "RdBu")) )(255), Rowv=TRUE, scale="none", trace="none", dendrogram = "both")


#subset and create heatmap for selected genes only 

sel_genes <- c("ETV4_extended (303g)", "VENTX_extended (58g)", "POU5F1 (64g)", "SOX2 (284g)", "GATA3 (243g)", 
               "GATA2 (825g)", "FOXA2 (31g)", "GATA4 (143g)", "SOX17 (24g)")
heat_m_sel<-subset(heat_m, rownames(heat_m) %in% sel_genes)

input_sel <- t(scale(t(heat_m_sel)))
heatmap.2(input_sel, col=colorRampPalette(rev(brewer.pal(9, "RdBu")) )(255), Rowv=TRUE, scale="none", trace="none", 
          dendrogram = "both", mar=c(10,10))

#subset and create heatmap for selected (+top 5) genes 


sel_genes_1 <- c("ETV4_extended (303g)", "VENTX_extended (58g)", "POU5F1 (64g)", "SOX2 (284g)", "GATA3 (243g)", 
               "GATA2 (825g)", "FOXA2 (31g)", "GATA4 (143g)", "SOX17 (24g)", "FOXN2 (431g)", 
               "CREB1_extended (443g)", "IRF4_extended (81g)", "NFYA (488g)	0.78361", "FIGLA (127g)	0.74511",
               "TGIF1_extended (55g)", "KDM5B (1177g)", "SOX15_extended (25g)", "ZNF256 (15g)",
               "CREM_extended (515g)", "GATA6 (345g)", "KLF3_extended (43g)", "KLF8_extended (349g)", 
               "FOXI3_extended (37g)", "NFKB2 (53g)", "ESRRB_extended (26g)", "MECOM_extended (27g)", 
               "TRIM28_extended (1692g)", "TCF7L1 (14g)", "KLF8_extended (349g)", "MYBL2_extended (23g)",
               "PRDM1 (13g)", "STAT2 (24g)", "GSC (12g)", "HMGN3_extended (37g)", "SMAD6 (19g)", "ISL1_extended (17g)",
               "ARX_extended (111g)", "TFAP2A (99g)", "NR2F2_extended (222g)", "CDX1 (19g)", "DLX3_extended (396g)",
               "GRHL1 (354g)", "TEAD3_extended (180g)", "RXRA (12g)", "GCM1 (33g)", "MSX2_extended (237g)",
               "PPARG (29g)")

heat_m_sel1<-subset(heat_m, rownames(heat_m) %in% sel_genes_1)

input_sel1 <- t(scale(t(heat_m_sel1)))
heatmap.2(input_sel1, col=colorRampPalette(rev(brewer.pal(9, "RdBu")) )(255), Rowv=TRUE, Colv=F, scale="none", trace="none", 
          dendrogram = "none", mar=c(10,10))

               





