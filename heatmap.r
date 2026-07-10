#!/usr/bin/env Rscript
# v1.4.0 2025-07-04
# 重构，并使用ComplexHeatmap包 

#---- 环境设置 ----
suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(htmlwidgets))
suppressMessages(library(plotly))
suppressMessages(library(heatmaply))
suppressMessages(library("optparse"))
suppressMessages(library(logging))
suppressMessages(library(Cairo))

#---- 功能函数 ----

#' @title 解析命令行选项
#' @description 定义并解析脚本运行时使用的所有命令行参数。
#' @return 返回一个包含已解析命令行参数的列表。
parse_options <- function() {
  option_list <- list(
    make_option(c("-i", "--indata"), type = "character", help = "【必须】输入数据文件路径（TSV 格式），第一列为 ID 或行名，其余列为样本值"),
    make_option(c("-g", "--group"), type = "character", help = "【必须】分组信息文件路径（TSV 格式），包含 sample 列及其他分组列，用于列注释"),
    make_option(c("-s", "--show_rownames"), type = "character", default = "auto", help = "是否显示行名。默认：auto（行数≤100 时显示）。格式：T/F/auto。影响：热图行标签可见性"),
    make_option(c("-l", "--show_colnames"), type = "character", default = "auto", help = "是否显示列名。默认：auto（列数≤100 时显示）。格式：T/F/auto。影响：热图列标签可见性及旋转 90 度"),
    make_option(c("--cluster_cols"), type = "character", default = "auto", help = "是否对列进行聚类。默认：auto（列数≤65535 时聚类）。格式：T/F/auto。影响：列 dendrogram 及排序"),
    make_option(c("--cluster_rows"), type = "character", default = "auto", help = "是否对行进行聚类。默认：auto（行数≤65535 时聚类）。格式：T/F/auto。影响：行 dendrogram 及排序，趋势分析需 TRUE"),
    make_option(c("-f", "--idfile"), type = "character", default = "None", help = "ID 筛选文件路径（单列文本）。默认：None（不筛选）。格式：每行一个 ID。影响：仅保留与数据匹配的 ID 行"),
    make_option(c("--label"), type = "character", default = "None", help = "特殊标记基因列表。默认：None（不标记）。格式：逗号分隔的基因名（如\"TP53,EGFR\"）。影响：热图右侧标注特定行"),
    make_option(c("-r", "--rowname"), type = "character", default = "ID", help = "行名构建列名。默认：ID（使用原行名）。格式：数据中的列名。影响：行名格式为\"原 ID(列值)\"，若值为 NA 或\"-\"则保持原 ID"),
    make_option(c("-c", "--trend"), type = "logical", default = FALSE, help = "是否进行趋势分析。默认：FALSE。格式：T/F。影响：生成趋势分析图（需 cluster_rows=TRUE）"),
    make_option(c("--row_anno_name"), type = "character", default = "Class,diffType", help = "行注释列名优先级。默认：Class,diffType。格式：逗号分隔列名。影响：按顺序查找首个存在的列作为行注释"),
    make_option(c("-a", "--scaletf"), type = "character", default = "auto", help = "数据标准化方式。默认：auto（样本>2 时 row，否则 none）。格式：row/none/auto。影响：row 时按行 z-score 标准化，none 时仅 log2 变换"),
    make_option(c("--corr"), type = "logical", default = FALSE, help = "是否绘制相关性热图。默认：FALSE。格式：T/F。影响：生成 Pearson 相关性矩阵及热图（.corr.xls/.png/.pdf）"),
    make_option(c("--heatmap"), type = "logical", default = TRUE, help = "是否绘制主热图。默认：TRUE。格式：T/F。影响：生成 ComplexHeatmap 热图（.png/.pdf）"),
    make_option(c("--circos"), type = "logical", default = FALSE, help = "是否绘制环形热图。默认：FALSE。格式：T/F。影响：生成 circular 热图，独立于主热图"),
    make_option(c("--circos_angle"), type = "numeric", default = 30, help = "环形热图开口角度。默认：30。格式：角度值（度）。影响：环形热图的开口大小"),
    make_option(c("--circos_width"), type = "numeric", default = 8, help = "环形热图宽度。默认：8。格式：英寸。影响：输出图像宽度"),
    make_option(c("--circos_height"), type = "numeric", default = 8, help = "环形热图高度。默认：8。格式：英寸。影响：输出图像高度"),
    make_option(c("--circos_fontsize"), type = "numeric", default = 0.5, help = "环形热图字体大小。默认：0.5。格式：相对单位。影响：标签文字大小"),
    make_option(c("--circos_format"), type = "character", default = "svg", help = "环形热图输出格式。默认：svg。格式：pdf/png/svg。影响：输出文件类型"),
    make_option(c("--cloud"), type = "logical", default = FALSE, help = "是否是云平台"),
    make_option(c("-o", "--prefix"), type = "character", default = "./all", help = "输出文件前缀路径。默认：./all。格式：路径/文件名。影响：所有输出文件的命名基础")
  )

  parser <- OptionParser(option_list = option_list, usage = "Rscript %prog [options]")
  parse_args(parser)
}

#' @title 加载预处理数据
#' @description 负责加载、筛选、格式化输入数据和分组信息。
#' @param opt 包含命令行选项的列表。
#' @return 返回一个列表，包含处理后的数据矩阵和分组信息。
load_and_prepare_data <- function(opt) {
  data_raw <- read_tsv(opt$indata, col_names = TRUE, col_types = cols(), locale = locale(encoding = "UTF-8"), na = c("NA", "N/A", "nan"))
  colnames(data_raw)[1] <- "ID"
  data_raw <- data_raw %>% column_to_rownames(var = "ID")

  if (opt$idfile != "None") {
    tryCatch({
      id_df <- read_tsv(opt$idfile, col_names = FALSE, col_types = cols(), locale = locale(encoding = "UTF-8"), na = c("NA", "N/A", "nan"))
      print(head(id_df[[1]]))
      print(head(rownames(data_raw)))
      id_intersect <- unique(intersect(id_df[[1]], rownames(data_raw)))
      nrow_raw = nrow(data_raw)
      data_raw <- data_raw[id_intersect, ]
      loginfo("已根据idfile筛选数据: 总行数: %s, list数量: %s, 筛选后数量: %s", nrow_raw, nrow(id_df), nrow(data_raw) , logger = script_name)
    }, error = function(e) {
      logerror("处理idfile时出错: %s", e$message, logger = script_name)
    })
  }

  if (opt$rowname %in% colnames(data_raw)) {
    data_raw <- data_raw %>% rownames_to_column("ID_temp") %>%
      mutate(rowname = ifelse(!is.na(.data[[opt$rowname]]) & .data[[opt$rowname]] != "-",
                              paste0(ID_temp, "-", .data[[opt$rowname]]), ID_temp))
      df_name <- data_raw %>% select(any_of(c("ID_temp", "rowname", opt$rowname)))
      colnames(df_name) = c("ID_temp", "rowname", "name")
      data_raw <- data_raw %>% column_to_rownames("rowname") %>% select(-ID_temp)
    loginfo("已使用'%s'列创建行名。", opt$rowname, logger = script_name)
  }

  group_df <- read_tsv(opt$group, col_names = TRUE, col_types = cols(.default = col_character()), locale = locale(encoding = "UTF-8"), na = c("NA", "N/A", "nan", ""))
  out_group_path = paste0(opt$prefix, ".group.xls")
  if (ncol(group_df) == 2 && any(duplicated(group_df[[1]]))) {
      output_path <- "sample_all.txt"
      cmd <- paste("python3",file.path(Bin,"get_pp_sample.py"), opt$group, out_group_path)
      system(cmd)
      loginfo("分组文件 '%s' 只有两列且存在重复样本，已调用 get_pp_sample.py 处理并保存至 %s", opt$group, out_group_path, logger = paste(script_name, "outfile", sep = ":"))
      group_df <- read_tsv(out_group_path, col_names = TRUE, col_types = cols(.default = col_character()), locale = locale(encoding = "UTF-8"), na = c("NA", "N/A", "nan", ""))
  }else{
      write_tsv(group_df, out_group_path)
      loginfo("已保存分组信息至 %s", out_group_path, logger = paste(script_name, "outfile", sep = ":"))
  }

  sample_intersect <- unique(intersect(group_df$sample, colnames(data_raw)))
  if (length(sample_intersect) == 0) {
    logerror("分组信息中没有与数据文件中样本列名匹配的样本。", logger = script_name)
    for (i in c("_count", "_fpkm", "_tpm", "_norm")) {
      logwarn("尝试移除 '%s' 后缀。", i, logger = script_name)
      colnames(data_raw) <- gsub(i, "", colnames(data_raw))
      sample_intersect <- unique(intersect(group_df$sample, colnames(data_raw)))
      if (length(sample_intersect) > 0) {
        loginfo("已根据 '%s' 列筛选数据: 筛选后数量: %s", i, length(sample_intersect), logger = script_name)
        break
      }
    }
  }
  loginfo(paste("共有样本:", length(sample_intersect),":", paste(sample_intersect, collapse = ", ")), logger = script_name)
  data <- data_raw[, sample_intersect]
  group_df <- group_df %>%
    filter(sample %in% sample_intersect) %>%
    distinct(sample, .keep_all = TRUE) %>%
    filter(!if_all(everything(), ~ is.na(.))) %>%
    select(where(~ !all(is.na(.))))

  rows_to_remove <- apply(data, 1, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))
  if (any(rows_to_remove)) {
    loginfo("已移除 %d 个没有方差的行。", sum(rows_to_remove), logger = script_name)
    data <- data[!rows_to_remove, ]
  }

  if (nrow(data) < 2) {
    loginfo("热图至少需要2个或更多的特征。", logger = script_name)
    quit(status = 1)
  }
  if (opt$rowname %in% colnames(data_raw)) {
    out_raw <- data %>%
      rownames_to_column("old_ID") %>%
      left_join(
        df_name %>% select(ID = ID_temp, rowname, name), 
        by = c("old_ID" = "rowname")
      ) %>%
      select(ID, name, everything(), -old_ID)
  }else{
    out_raw <- data %>% rownames_to_column("ID")
  }
  write_tsv(out_raw, file = paste0(opt$prefix, ".plot_data.raw.xls"))
  loginfo("已保存: %s.plot_data.raw.xls", opt$prefix, logger = paste(script_name, "outfile", sep = ":"))


  # class_info <- if (opt$row_anno_name %in% colnames(data_plot)) data_plot[[opt$row_anno_name]] else NULL
  row_anno_name = NULL
  row_anno_name_list = strsplit(opt$row_anno_name, ",")[[1]]
  for (rr in row_anno_name_list) {
    if (rr %in% colnames(data_raw)) {
      row_anno_name = rr
      loginfo("行注释使用 '%s' 列", rr, logger = script_name)
      break
    }
  }
  class_info = if (!is.null(row_anno_name)) data_raw[rownames(data),row_anno_name] else NULL

  return(list(row_anno_name = row_anno_name, class_info = class_info, data = data, groups = group_df))
}

#' @title 生成并保存相关性热图
#' @description 计算样本间的相关性，并绘制热图。
#' @param data_plot 数据矩阵。
#' @param group_df 分组信息的data frame。
#' @param opt 包含命令行选项的列表。
generate_correlation_plot <- function(data_plot, group_df, opt) {
  source(file.path(Bin, "corr.plot.source.r"))

  corr_df <- cor(data_plot, method = 'pearson', use = 'pairwise.complete.obs')
  prefix_corr <- paste0(opt$prefix, ".corr")

  write.table(corr_df, file = paste0(prefix_corr, ".xls"), col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)
  loginfo("相关性矩阵已保存至 %s.xls", prefix_corr, logger = paste(script_name, "outfile", sep = ":"))

  heatmap_p <- plot_corr(corr_df, group_df)

  height <- 4 + 0.5 * ncol(corr_df)
  width <- height

  loginfo("相关性图尺寸 (高x宽): %.2f x %.2f", height, width, logger = script_name)

  png(paste0(prefix_corr, ".png"), height = height, width = width, units = "in", res = 300)
  print(heatmap_p)
  dev.off()

  pdf(paste0(prefix_corr, ".pdf"), height = height, width = width)
  print(heatmap_p)
  dev.off()

  loginfo("相关性图已保存至 %s.png/pdf", prefix_corr, logger = paste(script_name, "outfile", sep = ":"))
}

#' @title 主执行流程
#' @description 整合所有步骤，执行完整的分析流程。
#' @param opt 包含命令行选项的列表。
main <- function(opt) {
  loginfo("脚本开始执行", logger = script_name)

  prepared_data <- load_and_prepare_data(opt)
  row_anno_name <- prepared_data$row_anno_name
  class_info <- prepared_data$class_info
  data_plot <- prepared_data$data
  group_df <- prepared_data$groups

  sample_n <- ncol(data_plot)
  compound_n <- nrow(data_plot)

  if (opt$scaletf == "auto") {
    if (sample_n > 2) {
      opt$scaletf <- "row"
    } else {
      opt$scaletf <- "none"
      data_plot <- t(scale(t(log2(data_plot + 1)), scale = FALSE)) %>% as.data.frame()
      logwarn("样本数等于2，将进行log2变换而非行归一化。", logger = script_name)
    }
  }

  anno_col_list <- get_anno_col(group_df)
  anno_row_list <- get_anno_row(class_info, rownames(data_plot), row_anno_name)

  opt$cluster_rows <- if (opt$cluster_rows == "auto") nrow(data_plot) <= 65535 else as.logical(opt$cluster_rows)
  opt$cluster_cols <- if (opt$cluster_cols == "auto") ncol(data_plot) <= 65535 else as.logical(opt$cluster_cols)

  loginfo("热图维度 (行x列): %d x %d", nrow(data_plot), ncol(data_plot), logger = script_name)
  if (compound_n * sample_n <= 30000 && compound_n <= 10000 && ! opt$cloud) {
    tryCatch({
      loginfo("开始生成交互式HTML热图", logger = script_name)
      heatmap_html(data_plot, anno_col_list, anno_row_list, opt)
    }, error = function(e) {
      logwarn("生成交互式HTML热图失败: %s", e, logger = script_name)
    })
  } else {
    logwarn("数据集过大，或者为云平台数据，跳过生成交互式HTML热图。", logger = script_name)
  }

  opt$show_rownames <- if (opt$show_rownames == "auto") compound_n <= 100 else as.logical(opt$show_rownames)
  opt$show_colnames <- if (opt$show_colnames == "auto") sample_n <= 100 else as.logical(opt$show_colnames)

  if (opt$heatmap) {
    heatmap_p(data_plot, anno_col_list, anno_row_list, opt)
  }

  if (opt$corr) {
    generate_correlation_plot(data_plot, group_df, opt)
  }

  if (opt$circos) {
    loginfo("开始绘制环形热图", logger = script_name)
    heatmap_circos(data_plot, opt)
  }

  loginfo("脚本执行结束", logger = script_name)
}

#---- 脚本执行 ----
argv <- commandArgs(trailingOnly = FALSE)
Bin <- dirname(substring(argv[grep("--file=", argv)], 8))
script_name <- basename(substring(argv[grep("--file=", argv)], 8))

source(file.path(Bin, "heatmap.source.r")) # , encoding = "UTF-8")

opt <- parse_options()
# print(opt)

if (is.null(opt$indata) || is.null(opt$group)) {
  logerror("--indata 和 --group 是必须提供的参数。使用 -h 查看帮助。", logger = script_name)
  quit(status = 1)
}

tryCatch({
  main(opt)
  writeLines("success", paste0(opt$prefix, ".finish"))
}, error = function(e) {
  writeLines(conditionMessage(e), paste0(opt$prefix, ".error"))
})

