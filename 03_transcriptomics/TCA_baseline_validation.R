# ============================================================
# #11 验证版 第二步：DESeq2 差值比较（difference of differences）
# 输入：
#   1) counts/tca_counts.txt        —— 上一步 featureCounts 输出
#   2) pja_gland_vs_HT_deseq2.csv   —— 您原分析 P. javanus 腺 vs 头胸
#      的 DESeq2 全基因结果表（至少含基因ID + log2FoldChange 两列）
#   3) rbh_orthologs.txt            —— 您原分析的 RBH 直系同源表
#      （两列：PJA基因ID  Tca基因ID）
#   4) core21_genes.txt             —— 21 个交集基因的 P. javanus 基因ID，一行一个
# 输出：
#   TableS30_core21_baseline_validation.csv   —— 21 基因逐一验证表
#   pja_hegs_baseline_summary.txt             —— 全体 PJA_HEGs 的汇总统计
# ============================================================

library(DESeq2)

# ---------- 1. T. castaneum：腺体 vs 腹部对照 ----------
cts <- read.delim("counts/tca_counts.txt", comment.char = "#", row.names = 1)
cts <- cts[, 6:8]                       # 三列计数：腹部对照、雄腺、雌腺
colnames(cts) <- c("SRR1206001", "SRR1206008", "SRR1206010")
coldata <- data.frame(
  row.names = colnames(cts),
  condition = factor(c("abdomen", "gland", "gland"), levels = c("abdomen", "gland"))
)
dds <- DESeqDataSetFromMatrix(cts, coldata, ~ condition)
dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]   # 过滤低表达
# 对照组只有 1 个样本：用局部拟合估计离散度，log2FC 点估计仍可用于效应量比较
dds <- DESeq(dds, fitType = "local")
tca <- as.data.frame(results(dds, contrast = c("condition", "gland", "abdomen")))
tca$Tca_id <- rownames(tca)
# 注意：padj 在 2 vs 1 设计下检验力很低，本分析以 log2FC 效应量为准，
# 不对 T. castaneum 单侧做显著性断言（会写进表注和 Limitations）

# ---------- 2. 读入 P. javanus 腺 vs 头胸的 log2FC ----------
pja <- read.csv("pja_gland_vs_HT_deseq2.csv")     # ← 改成您的实际文件名/列名
# 期望列：PJA_id, log2FC_PJA（腺体相对头胸）

# ---------- 3. RBH 直系同源表合并 ----------
rbh <- read.delim("rbh_orthologs.txt", header = FALSE,
                  col.names = c("PJA_id", "Tca_id"))   # ← 按实际格式调整
m <- merge(rbh, tca[, c("Tca_id", "log2FoldChange")], by = "Tca_id")
colnames(m)[ncol(m)] <- "log2FC_Tca"
m <- merge(m, pja[, c("PJA_id", "log2FC_PJA")], by = "PJA_id")

# ---------- 4. 差值的差值 ----------
m$delta_delta_LFC <- m$log2FC_PJA - m$log2FC_Tca
# 验证通过标准（两层）：
#   pass_strict: PJA 腺体富集(log2FC>=1) 且超出 Tca 基线 >= 1 个 log2
#   pass_lenient: PJA 腺体富集 且 delta_delta > 0
m$pass_strict  <- m$log2FC_PJA >= 1 & m$delta_delta_LFC >= 1
m$pass_lenient <- m$log2FC_PJA >= 1 & m$delta_delta_LFC > 0

# ---------- 5. 输出：21 个交集基因 ----------
core21 <- readLines("core21_genes.txt")           # ← 21 个 PJA 基因ID
tab21 <- m[m$PJA_id %in% core21, ]
tab21 <- tab21[order(-tab21$delta_delta_LFC), ]
write.csv(tab21, "TableS30_core21_baseline_validation.csv", row.names = FALSE)
cat("21 交集基因中成功比对到 T. castaneum 直系同源的：", nrow(tab21), "\n")
cat("pass_strict  (ΔΔLFC >= 1)：", sum(tab21$pass_strict), "\n")
cat("pass_lenient (ΔΔLFC > 0)：", sum(tab21$pass_lenient), "\n")

# ---------- 6. 输出：全体 PJA_HEGs 的汇总（供正文一句话引用） ----------
# pja_hegs.txt：1,688 个 PJA_HEGs 的 PJA 基因ID，一行一个
if (file.exists("pja_hegs.txt")) {
  hegs <- m[m$PJA_id %in% readLines("pja_hegs.txt"), ]
  sink("pja_hegs_baseline_summary.txt")
  cat("PJA_HEGs with RBH ortholog & quantifiable in both species:", nrow(hegs), "\n")
  cat("median log2FC_PJA:", median(hegs$log2FC_PJA, na.rm = TRUE), "\n")
  cat("median log2FC_Tca:", median(hegs$log2FC_Tca, na.rm = TRUE), "\n")
  cat("median delta-delta:", median(hegs$delta_delta_LFC, na.rm = TRUE), "\n")
  cat("fraction pass_lenient:", mean(hegs$pass_lenient, na.rm = TRUE), "\n")
  cat("fraction pass_strict :", mean(hegs$pass_strict,  na.rm = TRUE), "\n")
  sink()
}
cat("完成。\n")
