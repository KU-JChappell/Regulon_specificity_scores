#set wd
setwd("/mnt/nfs/data/Vincent_Lab/Joel_Chappell/SCENIC_petropoulos/int")

#read in binary regulon file
nonDupl <- readRDS("4.2_binaryRegulonActivity_nonDupl.Rds")

#convert to data frame
nonDupl1<-as.data.frame(nonDupl)

#write regulon matrix text file
write.table(nonDupl1, "regulon.activity.matrix.txt", sep = "\t", quote = F,
            row.names = T, col.names = T)

#set wd
setwd("/mnt/nfs/data/Vincent_Lab/Joel_Chappell/SCENIC_petropoulos")

#select required annotation column from cellInfo/format  
anno <- cellInfo[6]
anno <- cbind(Cell_names = rownames(anno), anno)
rownames(anno) <- NULL
colnames(anno)[2] <- "Stirparo_lineage"

#write annotation table - change for each dataset
write.table(anno, "cell.annotation.stirparo.txt", sep = "\t", quote = F,
            row.names = TRUE, col.names = TRUE)

#may have to edit text file in excel if not written appropriatley...