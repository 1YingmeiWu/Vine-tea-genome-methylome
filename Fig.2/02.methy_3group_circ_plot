args <- commandArgs(TRUE)
if (length(args) != 5) {
  print("usage: Rscript $0 <for_bedtools> <chrs.txt> <sample1> <sample2> <sample3>")
  q()
}

# 定义染色体排序函数
sort_chromosomes <- function(chr_names) {
  # 提取数字部分
  chr_numbers <- as.numeric(gsub("chr", "", chr_names))
  # 按数字排序
  sorted_names <- chr_names[order(chr_numbers)]
  return(sorted_names)
}

# 加载所需库
library('tidyr')
library('dplyr')
library('ggsci')
library('circlize')
library('grDevices')
library(ComplexHeatmap)

# 解析输入参数
for_bedtools_file <- args[1]
chrs_file <- args[2]
samples <- args[3:5]
meth_types <- c("CpG", "CHG", "CHH")

# 读取基础数据
input_gene <- read.table('genes_num.txt', header = F, stringsAsFactors = F)
colnames(input_gene) <- c('chr', 'start', 'end', 'value')

if (!file.exists('repeat_num.txt')) {
  stop("repeat_num.txt文件不存在，请检查文件路径")
}
input_repeat <- read.table('repeat_num.txt', header = F, stringsAsFactors = F)
colnames(input_repeat) <- c('chr', 'start', 'end', 'value')

# 读取染色体长度相关数据
for_bedtools <- read.table(for_bedtools_file, header = F, stringsAsFactors = F)
chrs <- read.table(chrs_file, header = F, stringsAsFactors = F)
length_df <- merge(chrs, for_bedtools, by = 'V1')
length_df <- data.frame(V1 = length_df$V1, V2 = 1, V3 = length_df$V2)

# 按染色体数字顺序排序
chr_order <- sort_chromosomes(as.character(length_df$V1))
length_df$V1 <- factor(length_df$V1, levels = chr_order)
length_df <- length_df[order(length_df$V1), ]

# 定义颜色方案
gene_color <- "#555555"
repeat_color <- "#91D1C2"

chr_pals <- c("#8B0000", "#A50000", "#B22222", "#C00000", "#CD0000", "#D70000", "#DC143C", "#E60000", "#EB3644", "#F00000", "#FF0000", "#FF2E2E", "#FF4C4C", "#FF6B6B", "#FF7F7F", "#FF9393", "#FFA7A7", "#FFBBBB", "#FFCFCF", "#FFE3E3")
#chr_pals <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
#              "#8491B4", "#91D1C2", "#EE4C97", "#FFDC91", "#0072B5",
#              "#7E6148", "#B09C85","#976A85","#9BCF5A","#6A5B8C","#C9D45C","#9BA9B9","#D49BB9","#7A8C5B","#C4C4C4","#7A7A7A")
chr_pal <- colorRampPalette(chr_pals)(NROW(length_df))

track_info <- list()

plot_mean_circlize <- function() {
  # 读取所有样本的甲基化数据
  all_meth_data <- list()
  for (sample in samples) {
    sample_meth <- list()
    for (meth in meth_types) {
      meth_file <- paste0(sample, ".", meth, ".methy.mean.txt")
      if (!file.exists(meth_file)) {
        stop(paste("甲基化文件不存在:", meth_file))
      }
      meth_data <- read.table(meth_file, header = F, stringsAsFactors = F)
      colnames(meth_data) <- c('chr', 'start', 'end', 'frac')
      sample_meth[[meth]] <- meth_data[order(meth_data$chr, meth_data$start), ]
    }
    all_meth_data[[sample]] <- sample_meth
  }

  # 计算甲基化数据范围
  meth_values <- unlist(lapply(all_meth_data, function(s)
    lapply(s, function(m) m$frac)))
  meth_min <- min(meth_values)
  meth_mid <- (meth_min + max(meth_values)) / 2
  meth_max <- max(meth_values)

  color_palette <- colorRamp2(c(meth_min, meth_mid, meth_max),
                             c("#4DBBD5","#F5D76E","#E64B35"))

  # 初始化PDF输出
  output_pdf <- paste0(paste(samples, collapse = "_"), '.mean.circlize.pdf')
  pdf(output_pdf, width = 15, height = 15)

  # 完全清除之前的设置
  circos.clear()

  # 设置circos参数
  circos.par(
    gap.degree = 2,
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.005, 0.005),
    start.degree = 90,
    circle.margin = 0.05,
    canvas.xlim = c(-0.9, 0.9),
    canvas.ylim = c(-0.9, 0.9),
    points.overflow.warning = FALSE
  )

  # 初始化染色体
  circos.genomicInitialize(
    length_df,
    plotType = NULL,
    sector.width = length_df$V3 / sum(length_df$V3),
    major.by = NULL
  )

  # ------------------------------
  # 轨道1：染色体标签和长度轨道（最外层，无背景色）- 不标记编号
  # ------------------------------
  circos.track(
    ylim = c(0, 1),
    track.height = 0.1,
    bg.border = NA,
    bg.col = NA,
    panel.fun = function(x, y) {
      xlim = CELL_META$xlim
      ylim = CELL_META$ylim
      seq_ID = CELL_META$sector.index
      chr_length = xlim[2]

      # 只在首尾显示长度标签
      ticks <- c(0, chr_length)

      for(i in seq_along(ticks)) {
        label <- round(ticks[i] / 1e6, 1)

        if(ticks[i] == 0) {
          circos.text(ticks[i] + chr_length * 0.02, 0.3, label, cex = 0.9,
                     facing = "clockwise", niceFacing = TRUE,
                     col = "black", adj = c(0, 0.5))
        } else {
          circos.text(ticks[i] - chr_length * 0.02, 0.3, label, cex = 0.9,
                     facing = "clockwise", niceFacing = TRUE,
                     col = "black", adj = c(1, 0.5))
        }
      }

      # 添加染色体名称
      circos.text(CELL_META$xcenter, 0.6, seq_ID, cex = 1.2,
                 facing = "bending.outside", niceFacing = TRUE,
                 col = "black", adj = c(0.5, 0.5), font = 2)
    },
    track.index = 1
  )
  track_info[["1"]] <- "染色体名称和首尾长度标签（不标记编号）"

  # ------------------------------
  # 轨道2：染色体颜色轨道 - 标记为1
  # ------------------------------
  circos.track(
    ylim = c(0, 1),
    track.height = 0.04,
    bg.border = 'black',
    bg.col = chr_pal,
    panel.fun = function(x, y, ...) {
      if (CELL_META$sector.index == chr_order[1]) {
        circos.text(
          CELL_META$xlim[1], 0.5,
          "1",  # 直接标记为1
          cex = 1.5, adj = c(0, 0.5), facing = "inside", font = 2
        )
      }
    },
    track.index = 2
  )
  track_info[["2"]] <- "染色体轨道（标记为1）"

  # ------------------------------
  # 轨道3：基因数量 - 标记为2
  # ------------------------------
  gene_max <- max(input_gene$value)
  circos.genomicTrack(
    input_gene,
    ylim = c(0, gene_max),
    track.height = 0.04,
    bg.col = 'white',
    bg.border = NA,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(
        region,
        ytop = value$value,
        ybottom = 0,
        col = gene_color,
        border = NA,
        lwd = 0.1
      )
      if (CELL_META$sector.index == chr_order[1]) {
        circos.text(
          CELL_META$xlim[1], gene_max/2,
          "2",  # 标记为2
          cex = 1.5, adj = c(0, 0.5), facing = "inside", font = 2
        )
      }
    },
    track.index = 3
  )
  track_info[["3"]] <- "基因数量（柱形，标记为2）"

  # ------------------------------
  # 轨道4：重复序列数量 - 标记为3
  # ------------------------------
  repeat_max <- max(input_repeat$value)
  circos.genomicTrack(
    input_repeat,
    ylim = c(0, repeat_max),
    track.height = 0.04,
    bg.col = 'white',
    bg.border = NA,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(
        region,
        ytop = value$value,
        ybottom = 0,
        col = repeat_color,
        border = NA,
        lwd = 0.1
      )
      if (CELL_META$sector.index == chr_order[1]) {
        circos.text(
          CELL_META$xlim[1], repeat_max/2,
          "3",  # 标记为3
          cex = 1.5, adj = c(0, 0.5), facing = "inside", font = 2
        )
      }
    },
    track.index = 4
  )
  track_info[["4"]] <- "重复序列数量（柱形，标记为3）"

  # ------------------------------
  # 轨道5-13：甲基化类型 - 标记为4到12
  # ------------------------------
  track_index <- 5
  track_number <- 4  # 从4开始编号
  
  for (sample in samples) {
    for (meth in meth_types) {
      meth_data <- all_meth_data[[sample]][[meth]]
      current_meth_max <- max(meth_data$frac)

      circos.genomicTrack(
        meth_data,
        ylim = c(0, current_meth_max),
        track.height = 0.04,
        bg.col = 'white',
        bg.border = NA,
        panel.fun = function(region, value, ...) {
          circos.genomicRect(
            region,
            value = value,
            col = color_palette(value$frac),
            border = NA,
            lwd = 0.1
          )
          if (CELL_META$sector.index == chr_order[1]) {
            circos.text(
              CELL_META$xlim[1], current_meth_max/2,
              as.character(track_number),  # 使用track_number
              cex = 1.5, adj = c(0, 0.5), facing = "inside", font = 2
            )
          }
        },
        track.index = track_index
      )

      track_info[[as.character(track_index)]] <- paste0("样本", sample, "的", meth, "（标记为", track_number, "）")
      track_index <- track_index + 1
      track_number <- track_number + 1
    }
  }

  # ------------------------------
  # 图例
  # ------------------------------
  # gene_legend <- Legend(
  #   labels = "gene",
  #   labels_gp = gpar(fontsize = 18, font = 2),
  #   grid_height = unit(1.2, "cm"),
  #   grid_width = unit(1.2, "cm"),
  #   background = gene_color,
  #   type = "grid",
  # #  pch = 15,
 
  
  #   border = "black"
  # )

  # repeat_legend <- Legend(
  #   labels = "repeat",
  #   labels_gp = gpar(fontsize = 18, font = 2),
  #   grid_height = unit(1.2, "cm"),
  #   grid_width = unit(1.2, "cm"),
  #   background = repeat_color,
  #   type = "grid",
  #  # pch = 15,
    
  
  #   border = "black"
  # )
  gene_legend <- Legend(
    labels = "gene",
    labels_gp = gpar(fontsize = 18, font = 2),
    grid_height = unit(1, "cm"),
    grid_width = unit(1, "cm"),
    background = gene_color,
    type = "points",
    pch = 15,
    size = unit(0.1, "mm"),  # 非常小的点
    legend_gp = gpar(col = gene_color, fill = gene_color)
  )

  repeat_legend <- Legend(
    labels = "repeat",
    labels_gp = gpar(fontsize = 18, font = 2),
    grid_height = unit(1, "cm"),
    grid_width = unit(1, "cm"),
    background = repeat_color,
    type = "points",
    pch = 15,
    size = unit(0.1, "mm"),  # 非常小的点
    legend_gp = gpar(col = repeat_color, fill = repeat_color)
  )

  # 甲基化类型图例
  gradient_legend <- Legend(
    at = round(seq(meth_min, meth_max, length.out = 3), 0),
    col_fun = color_palette,
    title = "5mC ratio(%)",
    title_gp = gpar(fontsize = 18, font = 2),
    labels_gp = gpar(fontsize = 16),
    direction = "horizontal",
    grid_width = unit(0.8, "npc"),
    grid_height = unit(0.4, "cm"),
    
  )

  # 组合图例
  top_legends <- packLegend(gene_legend, repeat_legend, direction = "vertical", gap = unit(1, "cm"))
  combined_legend <- packLegend(top_legends, gradient_legend, direction = "vertical", gap = unit(1, "cm"))

  # 中心位置放置图例
  pushViewport(viewport(
    x = unit(0.5, "npc"),
    y = unit(0.5, "npc"),
    just = "center",
    width = unit(0.6, "npc"),
    height = unit(0.5, "npc")
  ))
  draw(combined_legend, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = "center")
  upViewport()

  # 关闭PDF
  dev.off()

  cat("=== 环形图输出完成 ===\n")
  cat("文件路径:", output_pdf, "\n")
  cat("PDF尺寸: 15x15 英寸\n\n")
  cat("=== 轨道说明 ===\n")
  cat("轨道1: 染色体名称和首尾长度标签（不标记编号）\n")
  for (i in 2:13) {
    if (!is.null(track_info[[as.character(i)]])) {
      cat(sprintf("轨道%d: %s\n", i, track_info[[as.character(i)]]))
    }
  }
  cat("\n=== 优化说明 ===\n")
  cat("1. 染色体按数字顺序排序（chr1, chr2, chr3...）\n")
  cat("2. 轨道1不标记编号，轨道2开始标记1-12\n")
  cat("3. 顺时针方向排列染色体\n")
}

# 执行绘图函数
plot_mean_circlize()
