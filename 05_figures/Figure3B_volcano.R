# NOTE: paths are relative; place the input files listed in README.md next to this script.
# ============================================================================
# Figure 3B — Interspecific comparison: P. javanus vs T. castaneum defensive glands
# 新版图绘脚本 (2026-09): 三色证据分层火山图 + 候选基因标注 + 方向修正
#
# *** 重要修正 (务必阅读) ***
# 旧版脚本以 PJA 为 reference level, results(dds) 输出 log2FC = TCA/PJA,
# 因此 PJA 高表达基因落在负值 (左侧), 而旧图中红色集群在右侧 ——
# 即旧图的 x 轴方向与图注 "1,688 upregulated in P. javanus (red)" 存在不一致。
# 本版使用 results(dds, contrast = c("condition", "PJA", "TCA")),
# 输出 log2FC = PJA/TCA: 正值 (右侧) = PJA 高表达 = 红色, 与图注和正文一致。
# 差异检验本身不变 (同一模型, 仅符号取反), 1,688 / 3,227 / 4,915 数字不受影响。
#
# 输入文件:
#   1) RBH_with_expression.tsv        — RBH 直系同源对表达矩阵 (与旧版一致)
#   2) Fig3_baseline_validation_full.csv — 全基因基线验证 (Table S32 导出)
#   3) Fig3B_gene_annotation.csv      — 标注清单
#   4) [可选] RBH_with_expression_plus_abdomen.tsv — 见第 6 节
#
# 颜色语义 (与 Figure 3A 完全一致):
#   深红 = PJA 上调 且 超过 Tcas 种内基线 (strict ΔΔlog2FC 验证)
#   橙   = PJA 上调 但 未超过基线
#   蓝   = TCA 上调 (PJA 下调)
#   灰   = 不显著
# ============================================================================

# Input files expected in the working directory (see README.md).

library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)

## ---- 用户开关 -------------------------------------------------------------
SHOW_P3 <- FALSE          # TRUE = 连 priority-3 可选基因也标注
USE_LENIENT <- FALSE      # TRUE = 红色层改用 lenient 阈值 (555); FALSE = strict (494)

## ---- 字体 ------------------------------------------------------------------
# 用 "sans" (PDF 内置 Helvetica, 任何机器都能出图)。直接写 "Arial" 会在部分
# 系统触发 "failed to find or load PDF CID font" / "invalid font type"。
# 排版阶段如需 Arial, 在 Illustrator 里统一替换即可。
FONT <- "sans"

COL_VALIDATED <- "#C00000"
COL_UP_NOVAL  <- "#E69F00"
COL_DOWN      <- "#0072B2"
COL_NS        <- "grey82"

# ============================================================================
# 1. 数据读取与 DESeq2 (模型与旧版一致)
# ============================================================================
data <- read.csv("RBH_with_expression.tsv", sep = "\t", header = TRUE)
expr_data <- data[, c("PJA_expr_1", "PJA_expr_2", "PJA_expr_3", "TCA_expr_1", "TCA_expr_2")]
rownames(expr_data) <- data$PJA_gene
if (any(is.na(expr_data))) expr_data[is.na(expr_data)] <- 0
expr_data <- round(expr_data)

condition <- factor(c(rep("PJA", 3), rep("TCA", 2)))
col_data <- data.frame(row.names = colnames(expr_data), condition = condition)
dds <- DESeqDataSetFromMatrix(countData = expr_data, colData = col_data, design = ~ condition)
dds <- DESeq(dds)

# --- 方向修正: 显式对比 PJA vs TCA -> log2FC > 0 即 PJA 高表达 (右侧) ---
res <- results(dds, contrast = c("condition", "PJA", "TCA"))
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# ============================================================================
# 2. 合并基线验证状态 + 证据分层
# ============================================================================
baseline <- read.csv("Fig3_baseline_validation_full.csv", stringsAsFactors = FALSE)

# --- 关键修复 (2026-09-10) ---------------------------------------------------
# 该 CSV 由 Python 导出, 布尔值为 "True"/"False" 文本; 部分 R 环境读入后保持
# 字符型, 使 validated == TRUE 恒为 FALSE → 全部落入 not validated。
# 强制逻辑化, 两种写法都兼容:
baseline$pass_strict  <- toupper(trimws(as.character(baseline$pass_strict)))  == "TRUE"
baseline$pass_lenient <- toupper(trimws(as.character(baseline$pass_lenient))) == "TRUE"

# --- 合并前诊断 (请核对 Console 输出) ---------------------------------------
n_match <- sum(res_df$gene_id %in% baseline$gene_id)
cat(">> 基线表行数:", nrow(baseline),
    "| strict 通过:", sum(baseline$pass_strict), "(预期 494)",
    "| lenient 通过:", sum(baseline$pass_lenient), "(预期 555)\n")
cat(">> gene_id 匹配:", n_match, "/", nrow(res_df),
    "| 匹配率:", round(100 * n_match / nrow(res_df), 1), "%\n")
stopifnot("基线表 gene_id 与 DESeq2 结果零匹配, 请确认使用的是最新 Fig3_baseline_validation_full.csv" = n_match > 0)

baseline$validated <- if (USE_LENIENT) baseline$pass_lenient else baseline$pass_strict

plot_df <- res_df %>%
  filter(!is.na(padj)) %>%
  left_join(baseline[, c("gene_id", "validated")], by = "gene_id") %>%
  mutate(
    is_pja_up = padj < 0.05 & log2FoldChange >=  1,   # PJA 高表达 (右)
    is_tca_up = padj < 0.05 & log2FoldChange <= -1,   # TCA 高表达 (左)
    tier = case_when(
      is_pja_up & !is.na(validated) & validated ~ "Up in PJA, baseline-validated",
      is_pja_up                                 ~ "Up in PJA, not validated",
      is_tca_up                                 ~ "Up in TCA",
      TRUE                                      ~ "Not significant"
    ),
    tier = factor(tier, levels = c("Up in PJA, baseline-validated", "Up in PJA, not validated",
                                   "Up in TCA", "Not significant"))
  )

n_val  <- sum(plot_df$tier == "Up in PJA, baseline-validated")   # 预期 494 (strict)
n_nov  <- sum(plot_df$tier == "Up in PJA, not validated")        # 预期 1,194 (=1,688-494)
n_tca  <- sum(plot_df$tier == "Up in TCA")                       # 预期 3,227

# ============================================================================
# 3. 合并标注清单
# ============================================================================
ann <- read.csv("Fig3B_gene_annotation.csv", stringsAsFactors = FALSE)
if (!SHOW_P3) ann <- ann[ann$priority < 3, ]
# needs_check 已于 2026-09-10 销案: 4 个基因 (ACACA/CYP4G15/Dscam1/yem) 的 evm ID
# 与用户原始 96 基因清单逐一核对一致, 标注表中已全部置为 FALSE。

lab_df <- plot_df %>% inner_join(ann, by = "gene_id")
stopifnot(nrow(lab_df) > 0)

# ============================================================================
# 4. 火山图 (panel a)
# ============================================================================
xmax <- max(abs(plot_df$log2FoldChange), na.rm = TRUE) * 1.05

p_volcano <- ggplot(plot_df, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color = tier), size = 1.1, alpha = 0.55) +
  scale_color_manual(values = c("Up in PJA, baseline-validated" = COL_VALIDATED,
                                "Up in PJA, not validated"      = COL_UP_NOVAL,
                                "Up in TCA"                     = COL_DOWN,
                                "Not significant"               = COL_NS),
                     name = "Evidence tier") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_point(data = lab_df,
             aes(shape = label_class, fill = tier),
             size = 2.6, color = "black", stroke = 0.5) +
  scale_shape_manual(values = c("core" = 24, "candidate" = 21),
                     name = "Candidate status",
                     labels = c("core" = "Core intersection candidates (N = 21)",
                                "candidate" = "PPI-network candidates (dataset-specific)")) +
  scale_fill_manual(values = c("Up in PJA, baseline-validated" = COL_VALIDATED,
                               "Up in PJA, not validated"      = COL_UP_NOVAL,
                               "Up in TCA"                     = COL_DOWN,
                               "Not significant"               = COL_NS),
                    guide = "none") +
  geom_text_repel(data = lab_df[lab_df$label_class == "core", ],
                  aes(label = label), fontface = "bold.italic", size = 3.2,
                  max.overlaps = Inf, box.padding = 0.5, point.padding = 0.3,
                  segment.size = 0.25, seed = 42) +
  geom_text_repel(data = lab_df[lab_df$label_class == "candidate", ],
                  aes(label = label), fontface = "italic", size = 2.6,
                  max.overlaps = Inf, box.padding = 0.4, point.padding = 0.25,
                  segment.size = 0.2, seed = 42) +
  annotate("text", x = xmax * 0.98, y = Inf, vjust = 2.2, hjust = 1, size = 3.2,
           label = paste0("Up in PJA, validated: ", n_val, "\nUp in PJA, not validated: ", n_nov)) +
  annotate("text", x = -xmax * 0.98, y = Inf, vjust = 2.2, hjust = 0, size = 3.2,
           label = paste0("Up in TCA: ", n_tca), color = COL_DOWN) +
  coord_cartesian(xlim = c(-xmax, xmax), clip = "off") +
  labs(x = expression(log[2]~"Fold Change ("*italic("P. javanus")*" / "*italic("T. castaneum")*")"),
       y = expression(-log[10]~(italic(p)[adj]))) +
  theme_classic(base_size = 11, base_family = FONT) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        legend.position = "right",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"))

print(p_volcano)
ggsave("Fig3B_volcano.pdf",  p_volcano, width = 8.5, height = 6)
ggsave("Fig3B_volcano.tiff", p_volcano, width = 8.5, height = 6, dpi = 600, compression = "lzw")

# ============================================================================
# 5. 表达矩阵转换 (rlog) — PCA 与热图共用
# ============================================================================
rld <- rlog(dds, blind = FALSE)

# ============================================================================
# 6. PCA (panel b)
#
# [可选升级 — 回应 R1-4 的可视化证据]
# 旧版 PCA 只有 5 个腺体样本, PC1 占 99% (纯物种效应), PC2≈0%, 信息量低。
# 若将基线验证 (#11) 中定量过的 Tcas 腹部非腺对照样本加入, PCA 可同时呈现
# 组织效应与物种效应。做法:
#   1) 准备 RBH_with_expression_plus_abdomen.tsv: 在原 RBH 表基础上增加腹部
#      样本列, 列名以 "TCA_abd" 开头 (同一批 RBH 直系同源对在 Tcas 腹部样本
#      中的原始计数; 即基线验证使用的 Aabd 腹部对照)。
#      ※ 基线验证实际只用了 1 个腹部对照样本 (Aabd), 因此通常加 1 列即可;
#        脚本会自动识别任意数量 (1 个或 2 个) 的 TCA_abd 列。
#   2) 下方代码检测到该文件即自动启用 6 (或 7) 样本 PCA/热图;
#      未检测到则回退到 5 样本版本 (与旧版一致)。
# 注意: 差异检验 (第 1 节) 始终只用 5 个腺体样本, 1,688/3,227/4,915 不变;
#       腹部样本仅用于 PCA 与热图的展示。
# 单样本 caveat: 腹部组 n=1, 无组内方差, 其在 PCA 中的位置只作方向性参考
#   (illustrate, 非 demonstrate); 若采用, 图注中必须注明 n = 1。
# ============================================================================
EXT_FILE <- "RBH_with_expression_plus_abdomen.tsv"

if (file.exists(EXT_FILE)) {
  data_ext <- read.csv(EXT_FILE, sep = "\t", header = TRUE)
  # 自动识别腹部样本列 (基线验证实际只有 1 个腹部对照 Aabd; 1 或 2 列均可)
  abd_cols <- grep("^TCA_abd", names(data_ext), value = TRUE)
  stopifnot("扩展文件中未找到 TCA_abd* 列" = length(abd_cols) >= 1)
  n_abd <- length(abd_cols)
  cat(sprintf("检测到腹部样本扩展文件 (%d 个腹部样本), 使用 %d 样本 PCA/热图\n",
              n_abd, 5 + n_abd))
  if (n_abd == 1) {
    cat("※ 腹部组 n = 1, 无组内方差: PCA 中该点仅作方向性参考, 图注须注明 n = 1\n")
  }
  ext_cols <- c("PJA_expr_1", "PJA_expr_2", "PJA_expr_3",
                "TCA_expr_1", "TCA_expr_2", abd_cols)
  expr_ext <- round(as.matrix(data_ext[, ext_cols]))
  rownames(expr_ext) <- data_ext$PJA_gene
  if (any(is.na(expr_ext))) expr_ext[is.na(expr_ext)] <- 0
  grp <- factor(c(rep("PJA gland", 3), rep("TCA gland", 2),
                  rep("TCA abdomen (non-gland)", n_abd)),
                levels = c("PJA gland", "TCA gland", "TCA abdomen (non-gland)"))
  col_ext <- data.frame(row.names = ext_cols, group = grp)
  dds_ext <- DESeqDataSetFromMatrix(expr_ext, col_ext, design = ~ group)
  dds_ext <- estimateSizeFactors(dds_ext)
  rld_ext <- rlog(dds_ext, blind = FALSE)
  PCA_GROUPS <- c("PJA gland" = "#C00000", "TCA gland" = "#0072B2",
                  "TCA abdomen (non-gland)" = "#7A7A7A")
  rld_for_plots <- rld_ext; col_for_plots <- col_ext; ann_col_name <- "group"
  pal_for_plots <- PCA_GROUPS
} else {
  cat("未检测到腹部扩展文件, 使用 5 样本版本\n")
  rld_for_plots <- rld; col_for_plots <- col_data; ann_col_name <- "condition"
  pal_for_plots <- c("PJA" = "#C00000", "TCA" = "#0072B2")
}

pca_data <- plotPCA(rld_for_plots, intgroup = ann_col_name, returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = .data[[ann_col_name]])) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = sub("\\.sorted\\.bam$", "", name)),   # 去掉样本名后缀
                  size = 3, show.legend = FALSE, seed = 42) +
  scale_color_manual(values = pal_for_plots, name = "Groups") +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  theme_classic(base_size = 11, base_family = FONT) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8))
print(p_pca)
ggsave("Fig3B_PCA.pdf",  p_pca, width = 6.5, height = 5.5)
ggsave("Fig3B_PCA.tiff", p_pca, width = 6.5, height = 5.5, dpi = 600, compression = "lzw")

# ============================================================================
# 7. 聚类热图 (panel c)
#    与 3A 相同: 新版 scale="row" (行 Z-score), 图注写 row-scaled rlog counts
# ============================================================================
significant_genes <- subset(res_df, padj < 0.05 & abs(log2FoldChange) >= 1)
sig_ids <- intersect(rownames(significant_genes), rownames(rld_for_plots))
heat_matrix <- assay(rld_for_plots)[sig_ids, ]
# scale="row" 除以行标准差, 零方差行会产生 NaN 使聚类崩溃, 先剔除 (仅影响展示)
heat_matrix <- heat_matrix[apply(heat_matrix, 1, sd) > 0, ]
cat(">> 热图基因数 (去除零方差行后):", nrow(heat_matrix), "\n")

colors <- colorRampPalette(c("#0072B2", "white", "#C00000"))(50)
ann_colors <- setNames(list(pal_for_plots), ann_col_name)

pheatmap(heat_matrix,
         scale = "row",
         cluster_rows = TRUE, cluster_cols = TRUE,
         color = colors,
         annotation_col = col_for_plots,
         annotation_colors = ann_colors,
         show_rownames = FALSE,
         treeheight_row = 30, treeheight_col = 30,
         fontsize = 10, fontfamily = FONT,
         filename = "Fig3B_heatmap.pdf", width = 7, height = 6)
pheatmap(heat_matrix, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
         color = colors, annotation_col = col_for_plots, annotation_colors = ann_colors,
         show_rownames = FALSE, treeheight_row = 30, treeheight_col = 30,
         fontsize = 10, fontfamily = FONT,
         filename = "Fig3B_heatmap.tiff", width = 7, height = 6)

# ============================================================================
# 输出核对清单:
#   n_val 应为 494 (strict) 或 555 (lenient); n_val + n_nov = 1,688; n_tca = 3,227
#   lab_df: SHOW_P3=FALSE 时 34 (21 core + 13 candidate), TRUE 时 55
# ============================================================================
cat("\n=== 核对: PJA up validated =", n_val, "| PJA up not validated =", n_nov,
    "| TCA up =", n_tca, "| 标注基因数 =", nrow(lab_df), "===\n")
