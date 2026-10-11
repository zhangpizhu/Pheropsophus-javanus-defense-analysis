# NOTE: paths are relative; place the input files listed in README.md next to this script.
# ============================================================================
# Figure 3A — Intraspecific comparison: P. javanus defensive gland vs head-thorax
# 新版图绘脚本 (2026-09): 三色证据分层火山图 + 候选基因标注 + 色盲安全配色
#
# 输入文件:
#   1) all_counts.txt                 — 原始计数矩阵 (与旧版一致, 勿改)
#   2) Fig3_baseline_validation_full.csv — 全基因 T. castaneum 基线验证结果
#      (gene_id, pass_strict, pass_lenient; 由 Table S32 数据导出)
#   3) Fig3A_gene_annotation.csv      — 标注清单
#      (gene_id, label, label_class[core/candidate], priority[1/2/3],
#       baseline_validated, needs_check)
#
# 颜色语义 (与 Figure 3B 完全一致):
#   深红  = 腺体上调 且 通过跨种基线验证 (最高置信层)
#   橙    = 腺体上调 但 无跨种验证支持
#   蓝    = 下调
#   灰    = 不显著
# 标注语义:
#   ▲ 粗体  = core, 21 个多数据集交集核心候选 (必标)
#   ● 常规  = candidate priority 2, 故事相关基因 (建议标)
#   ● 小号  = candidate priority 3, 可选 (SHOW_P3 开关控制)
# ============================================================================

# Input files expected in the working directory (see README.md).

library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)

## ---- 用户开关 -------------------------------------------------------------
SHOW_P3 <- FALSE   # TRUE = 连 priority-3 可选基因也标注 (标签会更多)

## ---- 配色常量 (Okabe-Ito 色盲安全) ----------------------------------------
COL_VALIDATED <- "#C00000"   # 深红: 上调 + 基线验证
COL_UP_NOVAL  <- "#E69F00"   # 橙:   上调未验证
COL_DOWN      <- "#0072B2"   # 蓝:   下调
COL_NS        <- "grey82"    # 灰:   不显著

## ---- 字体 ------------------------------------------------------------------
# 用 "sans" (PDF 内置 Helvetica, 任何机器都能出图)。
# 不要直接写 "Arial": R 的 pdf() 设备在部分系统 (尤其装有中文字体的 Windows)
# 会把 Arial 解析成 CID 字体而报错 "failed to find or load PDF CID font"。
# 排版阶段如需 Arial, 在 Illustrator 里统一替换即可 (Helvetica 与 Arial 字形几乎一致)。
FONT <- "sans"

# ============================================================================
# 1. 数据读取与 DESeq2 (与旧版一致, 未改动)
# ============================================================================
df <- read.table("all_counts.txt", sep = "\t", header = TRUE, row.names = 1)
if (any(is.na(df))) df[is.na(df)] <- 0

devstage <- factor(c("defensive", "chest", "defensive", "chest", "defensive", "chest"))
columnData <- data.frame(row.names = colnames(df), devstage)
dds <- DESeqDataSetFromMatrix(countData = df, colData = columnData, design = ~ devstage)
dds <- DESeq(dds)

# 显式指定对比方向: defensive (gland) vs chest (head-thorax)
# -> log2FoldChange > 0 表示腺体高表达 (右侧), 与图注一致
res <- results(dds, contrast = c("devstage", "defensive", "chest"))
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# ============================================================================
# 2. 合并基线验证状态 + 证据分层
# ============================================================================
baseline <- read.csv("Fig3_baseline_validation_full.csv", stringsAsFactors = FALSE)

# --- 关键修复 (2026-09-10) ---------------------------------------------------
# 该 CSV 由 Python 导出, 布尔值写作 "True"/"False" 文本。部分 R 环境读入后
# pass_strict 保持字符型, 导致 pass_strict == TRUE 恒为 FALSE → 全部落入
# "not validated" (上次运行 n_val=0 的原因)。这里强制逻辑化, 两种写法都兼容:
baseline$pass_strict  <- toupper(trimws(as.character(baseline$pass_strict)))  == "TRUE"
baseline$pass_lenient <- toupper(trimws(as.character(baseline$pass_lenient))) == "TRUE"

# --- 合并前诊断 (请核对 Console 输出) ---------------------------------------
n_match <- sum(res_df$gene_id %in% baseline$gene_id)
cat(">> 基线表行数:", nrow(baseline),
    "| 其中 strict 通过:", sum(baseline$pass_strict), "(预期 494)\n")
cat(">> gene_id 匹配:", n_match, "/", nrow(res_df),
    "| 匹配率:", round(100 * n_match / nrow(res_df), 1), "%\n")
# 若匹配数为 0, 说明用的不是同一版文件 (如旧版用了 PJA|Chr 格式 ID), 直接报错:
stopifnot("基线表 gene_id 与 DESeq2 结果零匹配, 请确认使用的是最新 Fig3_baseline_validation_full.csv" = n_match > 0)

plot_df <- res_df %>%
  filter(!is.na(padj)) %>%
  left_join(baseline, by = "gene_id") %>%
  mutate(
    is_up   = padj < 0.05 & log2FoldChange >= 1,
    is_down = padj < 0.05 & log2FoldChange <= -1,
    tier = case_when(
      is_up & !is.na(pass_strict) & pass_strict ~ "Up, baseline-validated",
      is_up                                     ~ "Up, not validated",
      is_down                                   ~ "Down",
      TRUE                                      ~ "Not significant"
    ),
    tier = factor(tier, levels = c("Up, baseline-validated", "Up, not validated",
                                   "Down", "Not significant"))
  )

# 各层计数 (图内自动标注, 不硬编码)
n_val  <- sum(plot_df$tier == "Up, baseline-validated")   # 预期 454
n_nov  <- sum(plot_df$tier == "Up, not validated")        # 预期 518
n_down <- sum(plot_df$tier == "Down")

# ============================================================================
# 3. 合并标注清单
# ============================================================================
ann <- read.csv("Fig3A_gene_annotation.csv", stringsAsFactors = FALSE)
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
  scale_color_manual(values = c("Up, baseline-validated" = COL_VALIDATED,
                                "Up, not validated"      = COL_UP_NOVAL,
                                "Down"                   = COL_DOWN,
                                "Not significant"        = COL_NS),
                     name = "Evidence tier") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40", linewidth = 0.4) +
  # 候选基因标记层: 形状编码候选身份, 填充色仍随证据层
  geom_point(data = lab_df,
             aes(shape = label_class, fill = tier),
             size = 2.6, color = "black", stroke = 0.5) +
  scale_shape_manual(values = c("core" = 24, "candidate" = 21),   # 24=实心三角(带边) 21=实心圆(带边)
                     name = "Candidate status",
                     labels = c("core" = "Core intersection candidates (N = 21)",
                                "candidate" = "PPI-network candidates (dataset-specific)")) +
  scale_fill_manual(values = c("Up, baseline-validated" = COL_VALIDATED,
                               "Up, not validated"      = COL_UP_NOVAL,
                               "Down"                   = COL_DOWN,
                               "Not significant"        = COL_NS),
                    guide = "none") +
  # 标签层: core 粗体, candidate 常规
  geom_text_repel(data = lab_df[lab_df$label_class == "core", ],
                  aes(label = label), fontface = "bold.italic", size = 3.2,
                  max.overlaps = Inf, box.padding = 0.5, point.padding = 0.3,
                  segment.size = 0.25, seed = 42) +
  geom_text_repel(data = lab_df[lab_df$label_class == "candidate", ],
                  aes(label = label), fontface = "italic", size = 2.6,
                  max.overlaps = Inf, box.padding = 0.4, point.padding = 0.25,
                  segment.size = 0.2, seed = 42) +
  # 角标计数
  annotate("text", x = xmax * 0.98, y = Inf, vjust = 2.2, hjust = 1, size = 3.2,
           label = paste0("Up, validated: ", n_val, "\nUp, not validated: ", n_nov),
           color = "black") +
  annotate("text", x = -xmax * 0.98, y = Inf, vjust = 2.2, hjust = 0, size = 3.2,
           label = paste0("Down: ", n_down), color = COL_DOWN) +
  coord_cartesian(xlim = c(-xmax, xmax), clip = "off") +
  labs(x = expression(log[2]~"Fold Change (gland vs head-thorax)"),
       y = expression(-log[10]~(italic(p)[adj]))) +
  theme_classic(base_size = 11, base_family = FONT) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        legend.position = "right",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"))

print(p_volcano)
ggsave("Fig3A_volcano.pdf",  p_volcano, width = 8.5, height = 7)
ggsave("Fig3A_volcano.tiff", p_volcano, width = 8.5, height = 7, dpi = 600, compression = "lzw")

# ============================================================================
# 5. PCA (panel b) — 与旧版一致, 仅更新配色与标签防重叠
# ============================================================================
library(FactoMineR)
library(factoextra)

dat <- as.data.frame(t(counts(dds)))
rownames(dat) <- sub("\\.sorted\\.bam$", "", rownames(dat))   # 去掉样本名后缀, 图面更干净
dat_pca <- PCA(dat, graph = FALSE)
# 用有生物学含义的组名 (列顺序与 all_counts.txt 一致: 腺/头胸交替);
# levels 显式固定, 保证 红=腺体、蓝=头胸 的配色不被字母序打乱
group_list <- factor(c("Defensive gland", "Head-thorax", "Defensive gland",
                       "Head-thorax", "Defensive gland", "Head-thorax"),
                     levels = c("Defensive gland", "Head-thorax"))

p_pca <- fviz_pca_ind(dat_pca, col.ind = group_list, mean.point = FALSE,
                      addEllipses = TRUE, ellipse.type = "confidence", ellipse.level = 0.9,
                      legend.title = "Tissue", repel = TRUE, title = "",
                      palette = c("Defensive gland" = "#C00000", "Head-thorax" = "#0072B2"),
                      labelsize = 3.5, pointsize = 2.5) +
  theme_classic(base_size = 11, base_family = FONT) +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 0.8))
print(p_pca)
ggsave("Fig3A_PCA.pdf",  p_pca, width = 6.5, height = 5.5)
ggsave("Fig3A_PCA.tiff", p_pca, width = 6.5, height = 5.5, dpi = 600, compression = "lzw")

# ============================================================================
# 6. 聚类热图 (panel c)
#    注意: 旧版未设 scale 参数 (pheatmap 默认 scale="none", 展示 rlog 绝对值)
#    新版改为 scale="row" (行 Z-score), 展示跨样本的表达模式而非绝对丰度,
#    图注需相应写明 "row-scaled (Z-score) rlog-transformed counts"
# ============================================================================
rld <- rlog(dds, blind = FALSE)

significant_genes <- subset(res_df, padj < 0.05 & abs(log2FoldChange) >= 1)
heat_matrix <- assay(rld)[rownames(significant_genes), ]

# scale="row" 会除以行标准差: 零方差行 (所有样本表达相同) 会产生 NaN 使聚类崩溃,
# 先剔除 (只影响展示, 不影响任何统计结果)
heat_matrix <- heat_matrix[apply(heat_matrix, 1, sd) > 0, ]
cat(">> 热图基因数 (去除零方差行后):", nrow(heat_matrix), "\n")

colors <- colorRampPalette(c("#0072B2", "white", "#C00000"))(50)   # 与火山图同色系
ann_colors <- list(devstage = c(defensive = "#C00000", chest = "#0072B2"))

pheatmap(heat_matrix,
         scale = "row",
         cluster_rows = TRUE, cluster_cols = TRUE,
         color = colors,
         annotation_col = columnData,
         annotation_colors = ann_colors,
         show_rownames = FALSE,
         treeheight_row = 30, treeheight_col = 30,
         fontsize = 10, fontfamily = FONT,
         filename = "Fig3A_heatmap.pdf", width = 7, height = 6)
# 如需 TIFF:
pheatmap(heat_matrix, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
         color = colors, annotation_col = columnData, annotation_colors = ann_colors,
         show_rownames = FALSE, treeheight_row = 30, treeheight_col = 30,
         fontsize = 10, fontfamily = FONT,
         filename = "Fig3A_heatmap.tiff", width = 7, height = 6)

# ============================================================================
# 输出核对清单 (运行后在 Console 检查):
#   n_val 应为 454, n_nov 应为 518
#   lab_df 行数: SHOW_P3=FALSE 时 32 (21 core + 11 candidate), TRUE 时 45
# ============================================================================
cat("\n=== 核对: Up validated =", n_val, "(预期454) | Up not validated =", n_nov,
    "(预期518) | Down =", n_down, "| 标注基因数 =", nrow(lab_df), "===\n")
