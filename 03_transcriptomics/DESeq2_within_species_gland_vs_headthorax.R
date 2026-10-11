# NOTE: paths are relative; place the input files listed in README.md next to this script.
# Input expected in the working directory: all_counts.txt
# (featureCounts matrix: rows = genes, columns = 3 gland + 3 head-thorax samples)
library(DESeq2)

df <- read.table("all_counts.txt",sep="\t",header = T,row.names = 1)
any_na <- any(is.na(df))
if (any_na) {
  # 如果有缺失值，可以选择填充缺失值为特定的值，或者删除包含缺失值的行
  # 例如，填充缺失值为 0
  df[is.na(df)] <- 0
  # 或者删除包含缺失值的行
  # df <- na.omit(df)
}

devstage <-factor(c("defensive","chest","defensive","chest","defensive","chest"))
head(devstage)
columnData <-data.frame(row.names = colnames(df), devstage)
ddsFullCountTable <-DESeqDataSetFromMatrix(countData = df, colData = columnData, design= ~ devstage)
dds = DESeq(ddsFullCountTable)  #直接用DESeq函数做normalization
head(dds)
res <- results(dds)
diff_res <- as.data.frame(res)
diff_res$gene_name <- rownames(diff_res)
write.table(diff_res,file = 'DESeq2_diff_results.csv',quote = F,sep = ',',row.names = F,col.names = T)#输出差异结果

####设置阈值筛选差异表达基因####
library(dplyr)

significant_genes2 <- diff_res %>%
  filter(padj < 0.05 & abs(log2FoldChange) >= 1.0)
write.csv(significant_genes2,"PJA_DE_genes.csv")
####筛选上调下调表达基因####

DESeq2_DEG <- na.omit(diff_res)#删除包含缺失值的行
diff <- subset(DESeq2_DEG,padj < 0.05) #先筛选P值
up2 <- subset(diff,log2FoldChange > 1) #上调
down2 <- subset(diff,log2FoldChange < -1) #下调
write.csv(up2,"significant_up_genes2.csv")
write.csv(down2,"significant_down_genes2.csv")

#####另外标准化方法#####
normalized_counts <-counts(dds, normalized=TRUE)  #获得标准化后的数据
normalized_counts_mad <-apply(normalized_counts, 1, mad) #计算每一行的中位绝对偏差（MAD）
normalized_counts <-normalized_counts[order(normalized_counts_mad, decreasing=T), ]  #根据基因在不同的样本中表达变化的差异程度mad值对数据排序，差异越大的基因排位越前
rld <-rlog(dds, blind=FALSE)  # rlog转换后的标准化后的表达结果
rlogMat <-assay(rld) # 返回了一个DGEList对象rld 的数值部分
rlogMat <-rlogMat[order(normalized_counts_mad, decreasing=T), ]
write.table(rlogMat, file="hisat_gene_count_matrix.DESeq2.normalized.rlog.xls", quote=F, sep="\t", row.names=T, col.names=T)  ## rlog转换后的标准化后的表达结果



#####################PCA分析##################
#参考https://mp.weixin.qq.com/s/_ECziRSnb49GwEBkNQ-L5Q
##############################################
#install.packages("FactoMineR")
#install.packages("factoextra")
library(FactoMineR)
library(factoextra)
dat <- as.data.frame(t(counts(dds)))   # 转置dat矩阵
dat_pca <- PCA(dat, graph = FALSE)
group_list <- c("group1","group2","group1","group2","group1","group2")
fviz_pca_ind(dat_pca, col.ind=group_list, mean.point=F, addEllipses = T, legend.title="Groups", ellipse.type="confidence", ellipse.level=0.9, palette = c("#CC3333", "#339999"))+theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))#加个边框

##################样本相关性与热图##############
#install.packages("pheatmap")
library(pheatmap)
expcor <- cor(df, method = "spearman")
head(expcor)
pheatmap::pheatmap(expcor, clustering_method = "average",
                   treeheight_row = 0,treeheight_col = 0,
                   display_numbers = T)


# 筛选显著差异基因（padj < 0.05 且 |log2FoldChange| >= 1）
significant_genes <- subset(diff_res, padj < 0.05 & abs(log2FoldChange) >= 1)
significant_gene_names <- rownames(significant_genes)
# 提取显著差异基因的表达矩阵
heat_matrix1 <- assay(rld)[significant_gene_names, ]
# 设置颜色方案
colors <- colorRampPalette(c("blue", "white", "red"))(50)
# 绘制热图
pheatmap(heat_matrix1,
         cluster_rows = TRUE,  # 对行（基因）进行聚类
         cluster_cols = TRUE,  # 对列（样本）进行聚类
         color = colors,
         annotation_col = columnData, # 在列旁边添加样本的注释信息
         annotation_colors = list(devstage = c("defensive" = "blue", "chest" = "red")), # 自定义注释颜色
         annotation_legend = TRUE, # 显示注释图例
         show_rownames = FALSE) # 如果基因名太长，可以选择不显示












