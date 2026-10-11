# NOTE: paths are relative; place the input files listed in README.md next to this script.

# Requires: DESeq2, ggplot2, pheatmap, EnhancedVolcano, RColorBrewer (Bioconductor/CRAN)

# 加载必要的库
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(RColorBrewer)

# 1. 读取数据
data <- read.csv("RBH_with_expression.tsv", sep = "\t", header = TRUE)

# 2. 提取表达量矩阵和分组信息
expr_data <- data[, c("PJA_expr_1", "PJA_expr_2", "PJA_expr_3", "TCA_expr_1", "TCA_expr_2")]
rownames(expr_data) <- data$PJA_gene  # 使用 PJA_gene 作为行名

# 检查 NA 值并处理
if (any(is.na(expr_data))) {
  cat("NA values found in the count matrix. Replacing NA with 0...\n")
  expr_data[is.na(expr_data)] <- 0  # 替换 NA 为 0
}

# 确保计数数据为整数
expr_data <- round(expr_data)

# 创建分组信息
condition <- factor(c(rep("PJA", 3), rep("TCA", 2)))  # 3个PJA，2个TCA

# 确认 condition 的水平顺序
print("Levels of condition factor:")
print(levels(condition))

# 确保 PJA 为第一水平（对照组）
condition <- relevel(condition, ref = "PJA")
print("Re-leveled condition factor:")
print(levels(condition)) 


# 3. 构建 DESeq2 数据对象
col_data <- data.frame(row.names = colnames(expr_data), condition = condition)
dds <- DESeqDataSetFromMatrix(countData = expr_data, colData = col_data, design = ~ condition)


# 4. 运行 DESeq2 分析
dds <- DESeq(dds)
res <- results(dds)

# 5. 筛选 PJA 和 TCA 高表达基因
res_sorted <- res[order(res$padj), ]
PJA_high <- subset(res_sorted, padj < 0.05 & log2FoldChange < -1)  # PJA高表达
TCA_high <- subset(res_sorted, padj < 0.05 & log2FoldChange > 1)  # TCA高表达


# 保存差异表达结果
write.csv(as.data.frame(res_sorted), "DESeq2_results.csv")
write.csv(as.data.frame(PJA_high), "PJA_high_expression_genes.csv")
write.csv(as.data.frame(TCA_high), "TCA_high_expression_genes.csv")

# 6. 绘制 PCA 图
rld <- rlog(dds, blind = FALSE)
pca_data <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA Plot") +
  theme(panel.border  = element_rect(colour = "black", fill = NA))


# 7. 绘制层次差异基因聚类热图
diff_res <- as.data.frame(res)
diff_res$gene_name <- rownames(diff_res)
significant_genes <- subset(diff_res, padj < 0.05 & abs(log2FoldChange) >= 1)
significant_gene_names <- rownames(significant_genes)
write.csv(as.data.frame(significant_gene_names), "PJA_VS_TCA_DE_genes.csv")
# 提取显著差异基因的表达矩阵
heat_matrix1 <- assay(rld)[significant_gene_names, ]
# 设置颜色方案
colors <- colorRampPalette(c("blue", "white", "red"))(50)
# 绘制热图
pheatmap(heat_matrix1,
         cluster_rows = TRUE,  # 对行（基因）进行聚类
         cluster_cols = TRUE,  # 对列（样本）进行聚类
         color = colors,
         annotation_col = col_data, # 在列旁边添加样本的注释信息
         annotation_colors = list(devstage = c("PJA" = "blue", "TCA" = "red")), # 自定义注释颜色
         annotation_legend = TRUE, # 显示注释图例
         show_rownames = FALSE) # 如果基因名太长，可以选择不显示


# 8. 绘制差异火山图
# 添加一个新的列来标记显著性
res$significance <- "Not Significant"
res$significance[which(res$padj < 0.05 & res$log2FoldChange >= 1)] <- "TCA_upregulated"
res$significance[which(res$padj < 0.05 & res$log2FoldChange <= -1)] <- "PJA_upregulated"

ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  scale_color_manual(values = c("gray", "green", "red")) +
  labs(x = "log2 Fold Change", y = "-log10(padj)", title = "Volcano Plot") +
  theme(panel.background = element_rect(fill = "transparent",color = "black"))
