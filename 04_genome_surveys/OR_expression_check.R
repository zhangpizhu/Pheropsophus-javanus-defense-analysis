# NOTE: paths are relative; place the input files listed in README.md next to this script.
# ============================================================
# #10 假基因检查 · 步骤1：20 个 OR 基因的表达量核查
# 输入: all_counts.txt (featureCounts 原始计数，与主分析同一文件)
# 输出: OR_expression_check.txt（每个基因的标准化表达量汇总）
# ============================================================
# Input expected in the working directory: all_counts.txt (same featureCounts matrix as the main analysis)
library(DESeq2)

df <- read.table("all_counts.txt", sep = "\t", header = TRUE, row.names = 1)
df[is.na(df)] <- 0
devstage <- factor(c("defensive","chest","defensive","chest","defensive","chest"))
dds <- DESeqDataSetFromMatrix(df, data.frame(row.names = colnames(df), devstage), ~ devstage)
dds <- estimateSizeFactors(dds)          # 只做标准化，不做检验，几秒钟完成
nc <- counts(dds, normalized = TRUE)

# ---- 20 个目标基因：6 收缩家族 OR + 14 扩张家族 OR(对照组) ----
or_genes <- c(
  # 收缩组 (OG0000255 x4, OG0000347 x2) —— 检查主角
  "evm.model.Chr15.729","evm.model.Chr7.150","evm.model.Chr7.158",
  "evm.model.Chr7.159","evm.model.Chr10.410","evm.model.Chr7.77",
  # 扩张组 (OG0000197) —— 对照
  "evm.model.Chr6.790","evm.model.Chr4.822","evm.model.Chr4.829",
  "evm.model.Chr14.471","evm.model.Chr14.469","evm.model.Chr6.874",
  "evm.model.Chr6.873","evm.model.Chr17.166","evm.model.Chr1.365",
  "evm.model.Chr1.366","evm.model.Chr1.362","evm.model.Chr1.364",
  "evm.model.Chr6.872","evm.model.Chr6.875")

present <- intersect(or_genes, rownames(nc))
missing <- setdiff(or_genes, rownames(nc))
cat("== 在计数表中找到", length(present), "/ 20 个基因 ==\n")
if (length(missing) > 0) {
  cat("!! 未找到（ID 前缀不符或零计数被过滤）:\n"); print(missing)
  cat(">> 若全部未找到，先执行 head(rownames(nc)) 核对 ID 格式再改 or_genes\n")
}

sub <- nc[present, , drop = FALSE]
gland <- sub[, devstage == "defensive", drop = FALSE]
chest <- sub[, devstage == "chest", drop = FALSE]

out <- data.frame(
  gene        = present,
  group       = ifelse(present %in% or_genes[1:6], "contracted", "expanded"),
  gland_mean  = round(rowMeans(gland), 2),
  gland_max   = round(apply(gland, 1, max), 2),
  chest_mean  = round(rowMeans(chest), 2),
  raw_total   = rowSums(df[present, , drop = FALSE]),
  row.names = NULL)
out <- out[order(out$group, -out$gland_mean), ]

# 判定列：昆虫OR在非嗅觉组织通常沉默；这里只关心"有无转录证据"
out$expressed <- ifelse(out$gland_mean >= 10, "clearly",
                 ifelse(out$gland_mean >= 1,  "low", "near-zero"))
print(out, right = FALSE)
write.table(out, "OR_expression_check.txt", sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n>> 结果已保存: OR_expression_check.txt\n")
cat(">> 判读: gland_mean >= 10 = 明确表达; 1-10 = 低表达; <1 = 近零(腺体沉默,对OR属正常)\n")
