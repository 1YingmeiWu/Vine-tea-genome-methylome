# ============================================================
# ============================================================

# 加载所需包
library(ggplot2)
library(reshape2)
library(pheatmap)

# ---------- 1. 用户设置 ----------
file_MYB   <- "MYB-FPKM.xls"
file_bHLH  <- "bHLH-FPKM.xls"
file_DHM   <- "DHM-FPKM.xls"
file_F35H  <- "F35H-FPKM.xls"

cor_method <- "pearson"

# 输出文件前缀
out_prefix <- "TF_vs_synth"

# ---------- 2. 读取并整理单个表达矩阵 ----------
read_exp <- function(file) {
  df <- read.table(file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(df) < 3) stop(paste("文件", file, "列数不足"))
  gene_names <- df[, 2]                     # 基因名列
  if (any(duplicated(gene_names))) {
    warning(paste("文件", file, "中存在重复基因名，将保留第一个"))
    keep <- !duplicated(gene_names)
    df <- df[keep, ]
    gene_names <- gene_names[keep]
  }
  exp_mat <- as.matrix(df[, 3:ncol(df)])
  rownames(exp_mat) <- gene_names
  mode(exp_mat) <- "numeric"
  return(exp_mat)
}

cat("正在读取文件...\n")
myb_exp   <- read_exp(file_MYB)
bhlh_exp  <- read_exp(file_bHLH)
dhm_exp   <- read_exp(file_DHM)
f35h_exp  <- read_exp(file_F35H)

# ---------- 3. 合并转录因子与合成酶基因 ----------
tf_exp <- rbind(myb_exp, bhlh_exp)
cat("转录因子基因数量:", nrow(tf_exp), "\n")
syn_exp <- rbind(dhm_exp, f35h_exp)
cat("合成酶基因数量:", nrow(syn_exp), "\n")

if (!all(colnames(tf_exp) == colnames(syn_exp))) {
  cat("样本列名不完全一致，将按列顺序对齐并保留公共样本...\n")
  common_samples <- intersect(colnames(tf_exp), colnames(syn_exp))
  if (length(common_samples) == 0) stop("没有共同样本，无法计算相关性")
  tf_exp  <- tf_exp[, common_samples, drop = FALSE]
  syn_exp <- syn_exp[, common_samples, drop = FALSE]
}
cat("使用的样本数:", ncol(tf_exp), "\n")

# ---------- 4. 计算相关性矩阵 ----------
cor_mat <- matrix(NA, nrow = nrow(tf_exp), ncol = nrow(syn_exp),
                  dimnames = list(rownames(tf_exp), rownames(syn_exp)))

for (i in 1:nrow(tf_exp)) {
  for (j in 1:nrow(syn_exp)) {
    x <- as.numeric(tf_exp[i, ])
    y <- as.numeric(syn_exp[j, ])
    cor_mat[i, j] <- cor(x, y, method = cor_method, use = "complete.obs")
  }
}
cat("相关性计算完成\n")

# ---------- 5. 绘制热图 ----------
heatmap_file <- paste0(out_prefix, "_heatmap.pdf")
pheatmap(cor_mat,
         main = paste("Correlation (", cor_method, ")", sep = ""),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE,
         fontsize_number = 8,
         filename = heatmap_file,
         width = max(8, nrow(cor_mat)*0.2 + 2),
         height = max(6, ncol(cor_mat)*0.2 + 2))
cat("热图已保存:", heatmap_file, "\n")

# ---------- 6. 绘制散点图 ----------
tf_df <- as.data.frame(t(tf_exp))
tf_df$Sample <- rownames(tf_df)
syn_df <- as.data.frame(t(syn_exp))
syn_df$Sample <- rownames(syn_df)

# 融合为长格式
tf_long <- reshape2::melt(tf_df, id.vars = "Sample", variable.name = "TF", value.name = "TF_exp")
syn_long <- reshape2::melt(syn_df, id.vars = "Sample", variable.name = "Synth", value.name = "Synth_exp")

# 合并（按样本）
merged <- merge(tf_long, syn_long, by = "Sample")

# 添加相关系数标签（从 cor_mat 提取）
cor_df <- as.data.frame(as.table(cor_mat))
colnames(cor_df) <- c("TF", "Synth", "Cor")
merged <- merge(merged, cor_df, by = c("TF", "Synth"))

# 绘制散点图，分面展示
p <- ggplot(merged, aes(x = TF_exp, y = Synth_exp)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 0.8) +
  facet_grid(TF ~ Synth, scales = "free") +
  labs(title = paste("TF vs Synthesis gene expression (", cor_method, " correlation)", sep = ""),
       x = "TF expression (FPKM)", y = "Synthesis gene expression (FPKM)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 8))

# 保存散点图（可能较大，调整宽高）
scatter_file <- paste0(out_prefix, "_scatter.pdf")
ggsave(scatter_file, plot = p, width = min(20, nrow(cor_mat)*2 + 2), 
       height = min(15, ncol(cor_mat)*2 + 2), dpi = 300, limitsize = FALSE)
cat("散点图已保存:", scatter_file, "\n")

# ---------- 7. 输出相关系数表 ----------
cor_file <- paste0(out_prefix, "_correlation.csv")
write.csv(cor_mat, file = cor_file, quote = FALSE)
cat("相关系数表已保存:", cor_file, "\n")

cat("\n分析完成！\n")
