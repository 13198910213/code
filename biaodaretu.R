setwd('D:/first-project/data/')
rt <- read.table('pp-data/hcount0.6(1,0)-order by survival.txt',header=T, sep="\t",row.names= 1,fileEncoding="windows-1252")
cs <- read.table('raw-data/xh/cell-state.txt',header=T, sep="\t",row.names= NULL,fileEncoding="windows-1252")
fpkmfc <- read.table('process-data/fpkm-fc.txt',header=T, sep="\t",row.names= 1,fileEncoding="windows-1252")
View(fpkmfc[1:20,1:20])

sub <- rt[1,]
sub1 <- t(sub)
sub2 <- cbind(id=rownames(sub1),sub1)
rownames(sub2) <- NULL    #样本+按生存排的亚型11223344
sub3 <- data.frame(substring(sub2,1,12))   #保留样本名前12位+亚型11223344



csfc <- merge(cs,fpkmfc,by.x = "GeneName",by.y = "geneid_1") #1515
csfc1 <- csfc[order(csfc$Signature),]

csfc2 <- t(csfc1)
colnames(csfc2) <- csfc2[1,]
css <- data.frame(csfc2[-1,])
csfc3 <- data.frame(cbind(id=substring(rownames(css),1,12),css))
r <- data.frame(csfc3[1,])
r$Subtype <- "Subtype"
r1 <- cbind(id=r[,1],Subtype=r[,1517],r[,2:1516])
csfc4 <- merge(sub3,csfc3,by="id")
csfc44 <- csfc4[!duplicated(csfc4$id),]
csfc5 <- csfc44[order(csfc44$Subtype),]
csfc6 <- data.frame(rbind(r1,csfc5))

cc <- t(csfc6)
colnames(cc) <- cc[1,]
cc1 <- data.frame(cc[-1,])

cc2 <- cc1[-1,-1]
jy <- rownames(cc2)
cc3 <- data.frame(apply(cc2,2,as.numeric))
rownames(cc3) <- jy
boxplot(cc3)
cc3[cc3 == 0] <- 1
cc3[cc3 > 2.5] <- 2.5
cc4 <- cc3
boxplot(cc4)

subtype <- data.frame(t(cc1[1,][-1]))
co <- list(Subtype = c("1" = "#1F78B4","2" = "#B2DF8A","3" = "#33A02C","4" = "#A6CEE3"))

lb <- data.frame(lb=cc1[,1][-1])
rownames(lb) <- jy

co1 <- list(lb = colorRampPalette(c("#4682B4","#E6E6FA","#FF8C00"))(14))


bk <- c(seq(0.2,0.9,by=0.01),seq(0.91,1.29,by=0.01),seq(1.3,2.5,by=0.01))
mycoldown<-colorRampPalette(c("green"))(71)
mycolmedian<-colorRampPalette(c("white"))(39)
mycolup<-colorRampPalette(c("steelblue"))(120)

mycol<-c(mycoldown,mycolmedian,mycolup)
mycol<-colorRampPalette(c("steelblue","white","firebrick3"))(230)

HR <- cc4
library('pheatmap')
result <- pheatmap(HR,
                   scale = "none",                     # scale = "row"的含义是绘图时按行进行均一化。进行均一化可以降低个别特殊样品与其它样品间的差异，这会使得其它样品间的差异在图形中更加显著。一般我们在基于表达量进行聚类分析时，均是常用的参数。
                   #color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
                   color = mycol,
                   main = "23 Tumor-cell Example heatmap",                                             #设置标题
                   #fontsize_row = 4,                                                    #行字体大小
                   #fontsize_col = 1,                                                    #列字体大小
                   annotation_col = subtype,
                   annotation_colors = co,
                   annotation_row = lb,
                   #annotation_colors = co1,
                   show_colnames = F,
                   show_rownames = F,
                   #legend_breaks = c(0,1),
                   #legend_labels=c("0","1"),
                   cluster_rows = F,
                   cluster_cols = F,
                   #cutree_rows = 3,
                   breaks = bk,
                   clustering_distance_rows = "euclidean",
                   clustering_method_rows = "hc"
)

pheatmap(HR,fontsize=9, fontsize_row=6,
         colorRampPalette(c("steelblue","white","firebrick3"))(length(seq(-ceiling(max(HR[!is.na(HR) & !is.infinite(as.matrix(HR))])),
                                                                          ceiling(max(HR[!is.na(HR) & !is.infinite(as.matrix(HR))])),0.1))),
         cluster_cols = F,cluster_rows = F, 
         #cellwidth = 5, cellheight = 3,
         display_numbers = F,
         breaks = seq(-ceiling(max(HR[!is.na(HR) & !is.infinite(as.matrix(HR))])),
                      ceiling(max(HR[!is.na(HR) & !is.infinite(as.matrix(HR))])),0.1),
         legend_break = seq(-ceiling(max(HR[!is.na(HR) & !is.infinite(as.matrix(HR))])),
                            ceiling(max(HR[!is.na(HR) & !is.infinite(as.matrix(HR))])),0.1)
         )
dev.off()
