##############################################################################
# 脚本名称：single_gene_methylation_plots.R
# 功能描述：读取9个甲基化bed文件，为每个有效基因单独生成甲基化水平柱状图
# 输出文件：每个有效基因对应1个PNG（300DPI）和1个PDF（矢量图），存于当前目录
# 依赖包：tidyverse（包含ggplot2、dplyr等绘图与数据处理包）
# 使用方式：1. 确保gene_list和9个bed文件在工作目录；2. Rscript single_gene_methylation_plots.R
##############################################################################

# ============================================================================
# 第一步：加载依赖包并检查
# ============================================================================
# 检查tidyverse是否已安装，未安装则报错提示
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  stop("错误：需要安装'tidyverse'包！请先运行：install.packages('tidyverse')")
}
# 加载包（suppressMessages隐藏加载时的冲突提示，避免干扰）
suppressMessages(library(tidyverse))


# ============================================================================
# 第二步：配置基础参数（样本、甲基化类型、配色等，便于后续修改）
# ============================================================================
# 1. 样本列表（R/L/S，与bed文件名前缀对应）
samples <- c("R", "L", "S")
# 2. 甲基化类型（CpG/CHG/CHH，与bed文件名后缀对应）
contexts <- c("CpG", "CHG", "CHH")
# 3. 科研级配色（红色/R、蓝色/L、绿色/S，区分度高，适合期刊图）
sample_colors <- c(
  "R" = "#E74C3C",
  "L" = "#3498DB",
  "S" = "#2ECC71"
)
# 4. 图例标签（比样本名更完整，提升图的可读性）
sample_labels <- c(
  "R" = "Sample R",
  "L" = "Sample L",
  "S" = "Sample S"
)
# 5. 图表固定标签（坐标轴、标题前缀）
axis_x <- "Methylation Context"       # X轴标签（甲基化类型）
axis_y <- "Mean Methylation Level"    # Y轴标签（平均甲基化水平）
plot_title_prefix <- "Methylation Level of Gene: "  # 图标题前缀


# ============================================================================
# 第三步：读取基因列表并初始化记录变量
# ============================================================================
# 读取gene_list文件（无表头，第一列即为基因名）
gene_list_path <- "gene_list"
if (!file.exists(gene_list_path)) {
  stop(paste0("错误：基因列表文件'gene_list'不存在！请检查路径：", getwd()))
}
gene_list <- read.table(
  gene_list_path,
  header = FALSE,        # 无表头
  sep = "\t",            # 制表符分隔（兼容纯文本列）
  stringsAsFactors = FALSE  # 基因名按字符处理，不转因子
)$V1  # 取第一列作为基因列表

# 初始化记录变量（跟踪处理结果）
valid_genes <- c()       # 成功绘图的有效基因
skipped_genes <- c()     # 无数据的跳过基因
error_files <- c()       # 缺失的bed文件（用于报错）


# ============================================================================
# 第四步：遍历每个基因，提取数据并生成独立图表
# ============================================================================
cat("=== 开始分析：共", length(gene_list), "个基因 ===\n")

for (gene in gene_list) {
  # 临时变量：存储当前基因的甲基化数据（Context/Sample/Mean_Methylation）
  gene_data <- data.frame(
    Context = character(),
    Sample = character(),
    Mean_Methylation = numeric(),
    stringsAsFactors = FALSE
  )
  # 标记：当前基因是否有有效数据
  has_data <- FALSE
  
  # 嵌套循环：遍历所有9个bed文件（3样本×3甲基化类型）
  for (sample in samples) {
    for (context in contexts) {
      # 构建当前bed文件路径（如：R.CpG.bed、L.CHG.bed）
      bed_path <- paste0(sample, ".", context, ".bed")
      
      # 检查bed文件是否存在：不存在则记录，后续统一报错
      if (!file.exists(bed_path)) {
        error_files <- unique(c(error_files, bed_path))
        next  # 跳过当前文件，继续检查其他文件
      }
      
      # 读取bed文件（按制表符分隔，确保列识别正确）
      bed_data <- read.table(
        bed_path,
        header = FALSE,
        sep = "\t",
        quote = "",          # 忽略文件中的引号，避免列分裂错误
        na.strings = "",     # 空值按NA处理
        stringsAsFactors = FALSE
      )
      
      # 检查bed文件列数：至少7列（需提取第2列基因名、第7列甲基化水平）
      if (ncol(bed_data) < 7) {
        cat("警告：文件", bed_path, "列数不足7列，跳过该文件\n")
        next
      }
      
      # 提取当前基因的数据（第2列匹配基因名）
      gene_subset <- bed_data[bed_data$V2 == gene, ]
      # 若有数据：计算第7列的均值（甲基化水平），添加到gene_data
      if (nrow(gene_subset) > 0) {
        mean_meth <- mean(gene_subset$V7, na.rm = TRUE)  # na.rm排除NA值
        gene_data <- rbind(gene_data, data.frame(
          Context = context,
          Sample = sample,
          Mean_Methylation = mean_meth
        ))
        has_data <- TRUE  # 标记当前基因有有效数据
      }
    }
  }
  
  # --------------------------------------------------------------------------
  # 子步骤1：处理文件缺失错误（若有bed文件缺失，终止程序并提示）
  # --------------------------------------------------------------------------
  if (length(error_files) > 0) {
    stop(paste0(
      "错误：以下bed文件缺失，请补充后重新运行：\n",
      paste(error_files, collapse = "\n")
    ))
  }
  
  # --------------------------------------------------------------------------
  # 子步骤2：当前基因无数据→跳过；有数据→生成图表
  # --------------------------------------------------------------------------
  if (!has_data) {
    skipped_genes <- c(skipped_genes, gene)
    cat("→ 基因", gene, "：无数据，已跳过\n")
  } else {
    valid_genes <- c(valid_genes, gene)
    
    # 数据预处理：因子化Context和Sample，确保绘图顺序正确（按设定的contexts/samples顺序）
    gene_data$Context <- factor(gene_data$Context, levels = contexts)
    gene_data$Sample <- factor(gene_data$Sample, levels = samples)
    
    # ------------------------------------------------------------------------
    # 绘制单个基因的柱状图（ggplot2，科研级风格）
    # ------------------------------------------------------------------------
    p <- ggplot(gene_data, aes(x = Context, y = Mean_Methylation, fill = Sample)) +
      # 柱状图核心：黑色细边框（linewidth=0.2）、分组排列（dodge）、合适宽度
      geom_col(
        position = position_dodge(width = 0.8),
        width = 0.7,
        color = "black",
        linewidth = 0.2
      ) +
      # 配色：使用预设的sample_colors，图例标签用sample_labels
      scale_fill_manual(
        values = sample_colors,
        labels = sample_labels
      ) +
      # Y轴范围：从0开始（甲基化水平不能为负），顶部留15%空白（避免柱子顶到图边）
      scale_y_continuous(
        expand = c(0, 0),
        limits = c(0, max(gene_data$Mean_Methylation, na.rm = TRUE) * 1.15)
      ) +
      # 主题：简洁风格，去除多余网格线，调整字体大小（适合期刊）
      theme_bw() +
      theme(
        # 面板设置：去除网格线，添加浅灰色边框
        panel.grid = element_blank(),
        panel.border = element_rect(color = "gray70", linewidth = 0.2),
        # 坐标轴：字体大小11，X轴标签水平（不旋转）
        axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title = element_text(
          size = 13,
          color = "black",
          face = "bold",
          margin = margin(t = 0, r = 10, b = 0, l = 0)  # 标题与轴间距
        ),
        # 图例：放在右侧，字体大小11，标题加粗
        legend.position = "right",
        legend.title = element_text(size = 12, color = "black", face = "bold"),
        legend.text = element_text(size = 11, color = "black"),
        legend.key.width = unit(0.8, "cm"),  # 统一图例宽度
        legend.background = element_rect(fill = NA),  # 透明背景
        # 图标题：突出基因名，居中，加粗
        plot.title = element_text(
          size = 14,
          color = "black",
          face = "bold",
          hjust = 0.5,
          margin = margin(b = 15)  # 标题与图间距
        )
      ) +
      # 图表标签：X轴、Y轴、图例标题、图标题
      labs(
        x = axis_x,
        y = axis_y,
        fill = "Sample Group",  # 图例标题
        title = paste0(plot_title_prefix, gene)  # 图标题（含当前基因名）
      )
    
    # ------------------------------------------------------------------------
    # 保存图表：PNG（300DPI，适合展示）和PDF（矢量图，适合论文发表）
    # ------------------------------------------------------------------------
    # 构建输出文件名（基因名+后缀，如：ANR1_methylation.png）
    png_path <- paste0(gene, "_methylation.png")
    pdf_path <- paste0(gene, "_methylation.pdf")
    
    # 保存PNG（高分辨率，背景白色）
    ggsave(
      png_path,
      plot = p,
      width = 6,    # 固定宽度6英寸
      height = 5,   # 固定高度5英寸
      dpi = 300,    # 300DPI（印刷级分辨率）
      bg = "white", # 白色背景（避免透明）
      limitsize = FALSE
    )
    
    # 保存PDF（矢量图，无分辨率限制，可无限放大）
    ggsave(
      pdf_path,
      plot = p,
      width = 6,
      height = 5,
      device = "pdf",
      bg = "white",
      limitsize = FALSE
    )
    
    cat("→ 基因", gene, "：已生成 ", png_path, " 和 ", pdf_path, "\n", sep = "")
  }
}


# ============================================================================
# 第五步：输出分析总结（统计有效/跳过基因数，便于核对）
# ============================================================================
cat("\n=== 分析完成 ===")
cat("\n1. 处理结果统计：")
cat("\n   - 总基因数：", length(gene_list))
cat("\n   - 有效基因数（已绘图）：", length(valid_genes))
cat("\n   - 跳过基因数（无数据）：", length(skipped_genes))

if (length(skipped_genes) > 0) {
  cat("\n2. 跳过的基因列表：")
  cat("\n   ", paste(skipped_genes, collapse = ", "))
}

cat("\n3. 输出文件说明：")
cat("\n   - PNG文件：300DPI，适合PPT/报告展示")
cat("\n   - PDF文件：矢量图，适合期刊论文发表")
cat("\n   - 所有文件均保存在当前目录：", getwd(), "\n")##############################################################################
# 脚本名称：single_gene_methylation_plots.R
# 功能描述：读取9个甲基化bed文件，为每个有效基因单独生成甲基化水平柱状图
# 输出文件：每个有效基因对应1个PNG（300DPI）和1个PDF（矢量图），存于当前目录
# 依赖包：tidyverse（包含ggplot2、dplyr等绘图与数据处理包）
# 使用方式：1. 确保gene_list和9个bed文件在工作目录；2. Rscript single_gene_methylation_plots.R
##############################################################################

# ============================================================================
# 第一步：加载依赖包并检查
# ============================================================================
# 检查tidyverse是否已安装，未安装则报错提示
if (!requireNamespace("tidyverse", quietly = TRUE)) {
	  stop("错误：需要安装'tidyverse'包！请先运行：install.packages('tidyverse')")
}
# 加载包（suppressMessages隐藏加载时的冲突提示，避免干扰）
suppressMessages(library(tidyverse))


# ============================================================================
# 第二步：配置基础参数（样本、甲基化类型、配色等，便于后续修改）
# ============================================================================
# 1. 样本列表（R/L/S，与bed文件名前缀对应）
samples <- c("R", "L", "S")
# 2. 甲基化类型（CpG/CHG/CHH，与bed文件名后缀对应）
contexts <- c("CpG", "CHG", "CHH")
# 3. 科研级配色（红色/R、蓝色/L、绿色/S，区分度高，适合期刊图）
sample_colors <- c(
		     "R" = "#E74C3C",
		       "L" = "#3498DB",
		       "S" = "#2ECC71"
		       )
# 4. 图例标签（比样本名更完整，提升图的可读性）
sample_labels <- c(
		     "R" = "Sample R",
		       "L" = "Sample L",
		       "S" = "Sample S"
		       )
# 5. 图表固定标签（坐标轴、标题前缀）
axis_x <- "Methylation Context"       # X轴标签（甲基化类型）
axis_y <- "Mean Methylation Level"    # Y轴标签（平均甲基化水平）
plot_title_prefix <- "Methylation Level of Gene: "  # 图标题前缀


# ============================================================================
# 第三步：读取基因列表并初始化记录变量
# ============================================================================
# 读取gene_list文件（无表头，第一列即为基因名）
gene_list_path <- "gene_list"
if (!file.exists(gene_list_path)) {
	  stop(paste0("错误：基因列表文件'gene_list'不存在！请检查路径：", getwd()))
}
gene_list <- read.table(
			  gene_list_path,
			    header = FALSE,        # 无表头
			    sep = "\t",            # 制表符分隔（兼容纯文本列）
			      stringsAsFactors = FALSE  # 基因名按字符处理，不转因子
			    )$V1  # 取第一列作为基因列表

# 初始化记录变量（跟踪处理结果）
valid_genes <- c()       # 成功绘图的有效基因
skipped_genes <- c()     # 无数据的跳过基因
error_files <- c()       # 缺失的bed文件（用于报错）


# ============================================================================
# 第四步：遍历每个基因，提取数据并生成独立图表
# ============================================================================
cat("=== 开始分析：共", length(gene_list), "个基因 ===\n")

for (gene in gene_list) {
	  # 临时变量：存储当前基因的甲基化数据（Context/Sample/Mean_Methylation）
	  gene_data <- data.frame(
				      Context = character(),
				          Sample = character(),
				          Mean_Methylation = numeric(),
					      stringsAsFactors = FALSE
					    )
  # 标记：当前基因是否有有效数据
  has_data <- FALSE
    
    # 嵌套循环：遍历所有9个bed文件（3样本×3甲基化类型）
    for (sample in samples) {
	        for (context in contexts) {
			      # 构建当前bed文件路径（如：R.CpG.bed、L.CHG.bed）
			      bed_path <- paste0(sample, ".", context, ".bed")
        
        # 检查bed文件是否存在：不存在则记录，后续统一报错
        if (!file.exists(bed_path)) {
		        error_files <- unique(c(error_files, bed_path))
	        next  # 跳过当前文件，继续检查其他文件
		      }
	      
	      # 读取bed文件（按制表符分隔，确保列识别正确）
	      bed_data <- read.table(
				             bed_path,
					             header = FALSE,
					             sep = "\t",
						             quote = "",          # 忽略文件中的引号，避免列分裂错误
						             na.strings = "",     # 空值按NA处理
							             stringsAsFactors = FALSE
							           )
	      
	      # 检查bed文件列数：至少7列（需提取第2列基因名、第7列甲基化水平）
	      if (ncol(bed_data) < 7) {
		              cat("警告：文件", bed_path, "列数不足7列，跳过该文件\n")
	              next
		            }
	            
	            # 提取当前基因的数据（第2列匹配基因名）
	            gene_subset <- bed_data[bed_data$V2 == gene, ]
	            # 若有数据：计算第7列的均值（甲基化水平），添加到gene_data
	            if (nrow(gene_subset) > 0) {
			            mean_meth <- mean(gene_subset$V7, na.rm = TRUE)  # na.rm排除NA值
		            gene_data <- rbind(gene_data, data.frame(
								               Context = context,
									                 Sample = sample,
									                 Mean_Methylation = mean_meth
											         ))
			            has_data <- TRUE  # 标记当前基因有有效数据
			          }
		        }
    }
    
    # --------------------------------------------------------------------------
    # 子步骤1：处理文件缺失错误（若有bed文件缺失，终止程序并提示）
    # --------------------------------------------------------------------------
    if (length(error_files) > 0) {
	        stop(paste0(
			          "错误：以下bed文件缺失，请补充后重新运行：\n",
				        paste(error_files, collapse = "\n")
				      ))
      }
      
      # --------------------------------------------------------------------------
      # 子步骤2：当前基因无数据→跳过；有数据→生成图表
      # --------------------------------------------------------------------------
      if (!has_data) {
	          skipped_genes <- c(skipped_genes, gene)
          cat("→ 基因", gene, "：无数据，已跳过\n")
	    } else {
		        valid_genes <- c(valid_genes, gene)
	      
	      # 数据预处理：因子化Context和Sample，确保绘图顺序正确（按设定的contexts/samples顺序）
	      gene_data$Context <- factor(gene_data$Context, levels = contexts)
	          gene_data$Sample <- factor(gene_data$Sample, levels = samples)
	          
	          # ------------------------------------------------------------------------
	          # 绘制单个基因的柱状图（ggplot2，科研级风格）
	          # ------------------------------------------------------------------------
	          p <- ggplot(gene_data, aes(x = Context, y = Mean_Methylation, fill = Sample)) +
			        # 柱状图核心：黑色细边框（linewidth=0.2）、分组排列（dodge）、合适宽度
			        geom_col(
					         position = position_dodge(width = 0.8),
						         width = 0.7,
						         color = "black",
							         linewidth = 0.2
							       ) +
				      # 配色：使用预设的sample_colors，图例标签用sample_labels
				      scale_fill_manual(
							        values = sample_colors,
								        labels = sample_labels
								      ) +
				            # Y轴范围：从0开始（甲基化水平不能为负），顶部留15%空白（避免柱子顶到图边）
				            scale_y_continuous(
							               expand = c(0, 0),
								               limits = c(0, max(gene_data$Mean_Methylation, na.rm = TRUE) * 1.15)
								             ) +
					          # 主题：简洁风格，去除多余网格线，调整字体大小（适合期刊）
					          theme_bw() +
						        theme(
							              # 面板设置：去除网格线，添加浅灰色边框
							              panel.grid = element_blank(),
								              panel.border = element_rect(color = "gray70", linewidth = 0.2),
								              # 坐标轴：字体大小11，X轴标签水平（不旋转）
								              axis.text = element_text(size = 11, color = "black"),
									              axis.text.x = element_text(angle = 0, hjust = 0.5),
									              axis.title = element_text(
														          size = 13,
															            color = "black",
															            face = "bold",
																              margin = margin(t = 0, r = 10, b = 0, l = 0)  # 标题与轴间距
																            ),
							              # 图例：放在右侧，字体大小11，标题加粗
							              legend.position = "right",
								              legend.title = element_text(size = 12, color = "black", face = "bold"),
								              legend.text = element_text(size = 11, color = "black"),
									              legend.key.width = unit(0.8, "cm"),  # 统一图例宽度
									              legend.background = element_rect(fill = NA),  # 透明背景
										              # 图标题：突出基因名，居中，加粗
										              plot.title = element_text(
															          size = 14,
																            color = "black",
																            face = "bold",
																	              hjust = 0.5,
																	              margin = margin(b = 15)  # 标题与图间距
																		              )
										            ) +
							      # 图表标签：X轴、Y轴、图例标题、图标题
							      labs(
								           x = axis_x,
									           y = axis_y,
									           fill = "Sample Group",  # 图例标题
										           title = paste0(plot_title_prefix, gene)  # 图标题（含当前基因名）
										         )
							          
							          # ------------------------------------------------------------------------
							          # 保存图表：PNG（300DPI，适合展示）和PDF（矢量图，适合论文发表）
							          # ------------------------------------------------------------------------
							          # 构建输出文件名（基因名+后缀，如：ANR1_methylation.png）
							          png_path <- paste0(gene, "_methylation.png")
							          pdf_path <- paste0(gene, "_methylation.pdf")
								      
								      # 保存PNG（高分辨率，背景白色）
								      ggsave(
									           png_path,
										         plot = p,
										         width = 6,    # 固定宽度6英寸
											       height = 5,   # 固定高度5英寸
											       dpi = 300,    # 300DPI（印刷级分辨率）
											             bg = "white", # 白色背景（避免透明）
											             limitsize = FALSE
												         )
								      
								      # 保存PDF（矢量图，无分辨率限制，可无限放大）
								      ggsave(
									           pdf_path,
										         plot = p,
										         width = 6,
											       height = 5,
											       device = "pdf",
											             bg = "white",
											             limitsize = FALSE
												         )
								          
								          cat("→ 基因", gene, "：已生成 ", png_path, " 和 ", pdf_path, "\n", sep = "")
								        }
}


# ============================================================================
# 第五步：输出分析总结（统计有效/跳过基因数，便于核对）
# ============================================================================
cat("\n=== 分析完成 ===")
cat("\n1. 处理结果统计：")
cat("\n   - 总基因数：", length(gene_list))
cat("\n   - 有效基因数（已绘图）：", length(valid_genes))
cat("\n   - 跳过基因数（无数据）：", length(skipped_genes))

if (length(skipped_genes) > 0) {
	  cat("\n2. 跳过的基因列表：")
  cat("\n   ", paste(skipped_genes, collapse = ", "))
}

cat("\n3. 输出文件说明：")
cat("\n   - PNG文件：300DPI，适合PPT/报告展示")
cat("\n   - PDF文件：矢量图，适合期刊论文发表")
cat("\n   - 所有文件均保存在当前目录：", getwd(), "\n")##############################################################################
# 脚本名称：single_gene_methylation_plots.R
# 功能描述：读取9个甲基化bed文件，为每个有效基因单独生成甲基化水平柱状图
# 输出文件：每个有效基因对应1个PNG（300DPI）和1个PDF（矢量图），存于当前目录
# 依赖包：tidyverse（包含ggplot2、dplyr等绘图与数据处理包）
# 使用方式：1. 确保gene_list和9个bed文件在工作目录；2. Rscript single_gene_methylation_plots.R
##############################################################################

# ============================================================================
# 第一步：加载依赖包并检查
# ============================================================================
# 检查tidyverse是否已安装，未安装则报错提示
if (!requireNamespace("tidyverse", quietly = TRUE)) {
	  stop("错误：需要安装'tidyverse'包！请先运行：install.packages('tidyverse')")
}
# 加载包（suppressMessages隐藏加载时的冲突提示，避免干扰）
suppressMessages(library(tidyverse))


# ============================================================================
# 第二步：配置基础参数（样本、甲基化类型、配色等，便于后续修改）
# ============================================================================
# 1. 样本列表（R/L/S，与bed文件名前缀对应）
samples <- c("R", "L", "S")
# 2. 甲基化类型（CpG/CHG/CHH，与bed文件名后缀对应）
contexts <- c("CpG", "CHG", "CHH")
# 3. 科研级配色（红色/R、蓝色/L、绿色/S，区分度高，适合期刊图）
sample_colors <- c(
		     "R" = "#E74C3C",
		       "L" = "#3498DB",
		       "S" = "#2ECC71"
		       )
# 4. 图例标签（比样本名更完整，提升图的可读性）
sample_labels <- c(
		     "R" = "Sample R",
		       "L" = "Sample L",
		       "S" = "Sample S"
		       )
# 5. 图表固定标签（坐标轴、标题前缀）
axis_x <- "Methylation Context"       # X轴标签（甲基化类型）
axis_y <- "Mean Methylation Level"    # Y轴标签（平均甲基化水平）
plot_title_prefix <- "Methylation Level of Gene: "  # 图标题前缀


# ============================================================================
# 第三步：读取基因列表并初始化记录变量
# ============================================================================
# 读取gene_list文件（无表头，第一列即为基因名）
gene_list_path <- "gene_list"
if (!file.exists(gene_list_path)) {
	  stop(paste0("错误：基因列表文件'gene_list'不存在！请检查路径：", getwd()))
}
gene_list <- read.table(
			  gene_list_path,
			    header = FALSE,        # 无表头
			    sep = "\t",            # 制表符分隔（兼容纯文本列）
			      stringsAsFactors = FALSE  # 基因名按字符处理，不转因子
			    )$V1  # 取第一列作为基因列表

# 初始化记录变量（跟踪处理结果）
valid_genes <- c()       # 成功绘图的有效基因
skipped_genes <- c()     # 无数据的跳过基因
error_files <- c()       # 缺失的bed文件（用于报错）


# ============================================================================
# 第四步：遍历每个基因，提取数据并生成独立图表
# ============================================================================
cat("=== 开始分析：共", length(gene_list), "个基因 ===\n")

for (gene in gene_list) {
	  # 临时变量：存储当前基因的甲基化数据（Context/Sample/Mean_Methylation）
	  gene_data <- data.frame(
				      Context = character(),
				          Sample = character(),
				          Mean_Methylation = numeric(),
					      stringsAsFactors = FALSE
					    )
  # 标记：当前基因是否有有效数据
  has_data <- FALSE
    
    # 嵌套循环：遍历所有9个bed文件（3样本×3甲基化类型）
    for (sample in samples) {
	        for (context in contexts) {
			      # 构建当前bed文件路径（如：R.CpG.bed、L.CHG.bed）
			      bed_path <- paste0(sample, ".", context, ".bed")
        
        # 检查bed文件是否存在：不存在则记录，后续统一报错
        if (!file.exists(bed_path)) {
		        error_files <- unique(c(error_files, bed_path))
	        next  # 跳过当前文件，继续检查其他文件
		      }
	      
	      # 读取bed文件（按制表符分隔，确保列识别正确）
	      bed_data <- read.table(
				             bed_path,
					             header = FALSE,
					             sep = "\t",
						             quote = "",          # 忽略文件中的引号，避免列分裂错误
						             na.strings = "",     # 空值按NA处理
							             stringsAsFactors = FALSE
							           )
	      
	      # 检查bed文件列数：至少7列（需提取第2列基因名、第7列甲基化水平）
	      if (ncol(bed_data) < 7) {
		              cat("警告：文件", bed_path, "列数不足7列，跳过该文件\n")
	              next
		            }
	            
	            # 提取当前基因的数据（第2列匹配基因名）
	            gene_subset <- bed_data[bed_data$V2 == gene, ]
	            # 若有数据：计算第7列的均值（甲基化水平），添加到gene_data
	            if (nrow(gene_subset) > 0) {
			            mean_meth <- mean(gene_subset$V7, na.rm = TRUE)  # na.rm排除NA值
		            gene_data <- rbind(gene_data, data.frame(
								               Context = context,
									                 Sample = sample,
									                 Mean_Methylation = mean_meth
											         ))
			            has_data <- TRUE  # 标记当前基因有有效数据
			          }
		        }
    }
    
    # --------------------------------------------------------------------------
    # 子步骤1：处理文件缺失错误（若有bed文件缺失，终止程序并提示）
    # --------------------------------------------------------------------------
    if (length(error_files) > 0) {
	        stop(paste0(
			          "错误：以下bed文件缺失，请补充后重新运行：\n",
				        paste(error_files, collapse = "\n")
				      ))
      }
      
      # --------------------------------------------------------------------------
      # 子步骤2：当前基因无数据→跳过；有数据→生成图表
      # --------------------------------------------------------------------------
      if (!has_data) {
	          skipped_genes <- c(skipped_genes, gene)
          cat("→ 基因", gene, "：无数据，已跳过\n")
	    } else {
		        valid_genes <- c(valid_genes, gene)
	      
	      # 数据预处理：因子化Context和Sample，确保绘图顺序正确（按设定的contexts/samples顺序）
	      gene_data$Context <- factor(gene_data$Context, levels = contexts)
	          gene_data$Sample <- factor(gene_data$Sample, levels = samples)
	          
	          # ------------------------------------------------------------------------
	          # 绘制单个基因的柱状图（ggplot2，科研级风格）
	          # ------------------------------------------------------------------------
	          p <- ggplot(gene_data, aes(x = Context, y = Mean_Methylation, fill = Sample)) +
			        # 柱状图核心：黑色细边框（linewidth=0.2）、分组排列（dodge）、合适宽度
			        geom_col(
					         position = position_dodge(width = 0.8),
						         width = 0.7,
						         color = "black",
							         linewidth = 0.2
							       ) +
				      # 配色：使用预设的sample_colors，图例标签用sample_labels
				      scale_fill_manual(
							        values = sample_colors,
								        labels = sample_labels
								      ) +
				            # Y轴范围：从0开始（甲基化水平不能为负），顶部留15%空白（避免柱子顶到图边）
				            scale_y_continuous(
							               expand = c(0, 0),
								               limits = c(0, max(gene_data$Mean_Methylation, na.rm = TRUE) * 1.15)
								             ) +
					          # 主题：简洁风格，去除多余网格线，调整字体大小（适合期刊）
					          theme_bw() +
						        theme(
							              # 面板设置：去除网格线，添加浅灰色边框
							              panel.grid = element_blank(),
								              panel.border = element_rect(color = "gray70", linewidth = 0.2),
								              # 坐标轴：字体大小11，X轴标签水平（不旋转）
								              axis.text = element_text(size = 11, color = "black"),
									              axis.text.x = element_text(angle = 0, hjust = 0.5),
									              axis.title = element_text(
														          size = 13,
															            color = "black",
															            face = "bold",
																              margin = margin(t = 0, r = 10, b = 0, l = 0)  # 标题与轴间距
																            ),
							              # 图例：放在右侧，字体大小11，标题加粗
							              legend.position = "right",
								              legend.title = element_text(size = 12, color = "black", face = "bold"),
								              legend.text = element_text(size = 11, color = "black"),
									              legend.key.width = unit(0.8, "cm"),  # 统一图例宽度
									              legend.background = element_rect(fill = NA),  # 透明背景
										              # 图标题：突出基因名，居中，加粗
										              plot.title = element_text(
															          size = 14,
																            color = "black",
																            face = "bold",
																	              hjust = 0.5,
																	              margin = margin(b = 15)  # 标题与图间距
																		              )
										            ) +
							      # 图表标签：X轴、Y轴、图例标题、图标题
							      labs(
								           x = axis_x,
									           y = axis_y,
									           fill = "Sample Group",  # 图例标题
										           title = paste0(plot_title_prefix, gene)  # 图标题（含当前基因名）
										         )
							          
							          # ------------------------------------------------------------------------
							          # 保存图表：PNG（300DPI，适合展示）和PDF（矢量图，适合论文发表）
							          # ------------------------------------------------------------------------
							          # 构建输出文件名（基因名+后缀，如：ANR1_methylation.png）
							          png_path <- paste0(gene, "_methylation.png")
							          pdf_path <- paste0(gene, "_methylation.pdf")
								      
								      # 保存PNG（高分辨率，背景白色）
								      ggsave(
									           png_path,
										         plot = p,
										         width = 6,    # 固定宽度6英寸
											       height = 5,   # 固定高度5英寸
											       dpi = 300,    # 300DPI（印刷级分辨率）
											             bg = "white", # 白色背景（避免透明）
											             limitsize = FALSE
												         )
								      
								      # 保存PDF（矢量图，无分辨率限制，可无限放大）
								      ggsave(
									           pdf_path,
										         plot = p,
										         width = 6,
											       height = 5,
											       device = "pdf",
											             bg = "white",
											             limitsize = FALSE
												         )
								          
								          cat("→ 基因", gene, "：已生成 ", png_path, " 和 ", pdf_path, "\n", sep = "")
								        }
}


# ============================================================================
# 第五步：输出分析总结（统计有效/跳过基因数，便于核对）
# ============================================================================
cat("\n=== 分析完成 ===")
cat("\n1. 处理结果统计：")
cat("\n   - 总基因数：", length(gene_list))
cat("\n   - 有效基因数（已绘图）：", length(valid_genes))
cat("\n   - 跳过基因数（无数据）：", length(skipped_genes))

if (length(skipped_genes) > 0) {
	  cat("\n2. 跳过的基因列表：")
  cat("\n   ", paste(skipped_genes, collapse = ", "))
}

cat("\n3. 输出文件说明：")
cat("\n   - PNG文件：300DPI，适合PPT/报告展示")
cat("\n   - PDF文件：矢量图，适合期刊论文发表")
cat("\n   - 所有文件均保存在当前目录：", getwd(), "\n")##############################################################################
# 脚本名称：single_gene_methylation_plots.R
# 功能描述：读取9个甲基化bed文件，为每个有效基因单独生成甲基化水平柱状图
# 输出文件：每个有效基因对应1个PNG（300DPI）和1个PDF（矢量图），存于当前目录
# 依赖包：tidyverse（包含ggplot2、dplyr等绘图与数据处理包）
# 使用方式：1. 确保gene_list和9个bed文件在工作目录；2. Rscript single_gene_methylation_plots.R
##############################################################################

# ============================================================================
# 第一步：加载依赖包并检查
# ============================================================================
# 检查tidyverse是否已安装，未安装则报错提示
if (!requireNamespace("tidyverse", quietly = TRUE)) {
	  stop("错误：需要安装'tidyverse'包！请先运行：install.packages('tidyverse')")
}
# 加载包（suppressMessages隐藏加载时的冲突提示，避免干扰）
suppressMessages(library(tidyverse))


# ============================================================================
# 第二步：配置基础参数（样本、甲基化类型、配色等，便于后续修改）
# ============================================================================
# 1. 样本列表（R/L/S，与bed文件名前缀对应）
samples <- c("R", "L", "S")
# 2. 甲基化类型（CpG/CHG/CHH，与bed文件名后缀对应）
contexts <- c("CpG", "CHG", "CHH")
# 3. 科研级配色（红色/R、蓝色/L、绿色/S，区分度高，适合期刊图）
sample_colors <- c(
		     "R" = "#E74C3C",
		       "L" = "#3498DB",
		       "S" = "#2ECC71"
		       )
# 4. 图例标签（比样本名更完整，提升图的可读性）
sample_labels <- c(
		     "R" = "Sample R",
		       "L" = "Sample L",
		       "S" = "Sample S"
		       )
# 5. 图表固定标签（坐标轴、标题前缀）
axis_x <- "Methylation Context"       # X轴标签（甲基化类型）
axis_y <- "Mean Methylation Level"    # Y轴标签（平均甲基化水平）
plot_title_prefix <- "Methylation Level of Gene: "  # 图标题前缀


# ============================================================================
# 第三步：读取基因列表并初始化记录变量
# ============================================================================
# 读取gene_list文件（无表头，第一列即为基因名）
gene_list_path <- "gene_list"
if (!file.exists(gene_list_path)) {
	  stop(paste0("错误：基因列表文件'gene_list'不存在！请检查路径：", getwd()))
}
gene_list <- read.table(
			  gene_list_path,
			    header = FALSE,        # 无表头
			    sep = "\t",            # 制表符分隔（兼容纯文本列）
			      stringsAsFactors = FALSE  # 基因名按字符处理，不转因子
			    )$V1  # 取第一列作为基因列表

# 初始化记录变量（跟踪处理结果）
valid_genes <- c()       # 成功绘图的有效基因
skipped_genes <- c()     # 无数据的跳过基因
error_files <- c()       # 缺失的bed文件（用于报错）


# ============================================================================
# 第四步：遍历每个基因，提取数据并生成独立图表
# ============================================================================
cat("=== 开始分析：共", length(gene_list), "个基因 ===\n")

for (gene in gene_list) {
	  # 临时变量：存储当前基因的甲基化数据（Context/Sample/Mean_Methylation）
	  gene_data <- data.frame(
				      Context = character(),
				          Sample = character(),
				          Mean_Methylation = numeric(),
					      stringsAsFactors = FALSE
					    )
  # 标记：当前基因是否有有效数据
  has_data <- FALSE
    
    # 嵌套循环：遍历所有9个bed文件（3样本×3甲基化类型）
    for (sample in samples) {
	        for (context in contexts) {
			      # 构建当前bed文件路径（如：R.CpG.bed、L.CHG.bed）
			      bed_path <- paste0(sample, ".", context, ".bed")
        
        # 检查bed文件是否存在：不存在则记录，后续统一报错
        if (!file.exists(bed_path)) {
		        error_files <- unique(c(error_files, bed_path))
	        next  # 跳过当前文件，继续检查其他文件
		      }
	      
	      # 读取bed文件（按制表符分隔，确保列识别正确）
	      bed_data <- read.table(
				             bed_path,
					             header = FALSE,
					             sep = "\t",
						             quote = "",          # 忽略文件中的引号，避免列分裂错误
						             na.strings = "",     # 空值按NA处理
							             stringsAsFactors = FALSE
							           )
	      
	      # 检查bed文件列数：至少7列（需提取第2列基因名、第7列甲基化水平）
	      if (ncol(bed_data) < 7) {
		              cat("警告：文件", bed_path, "列数不足7列，跳过该文件\n")
	              next
		            }
	            
	            # 提取当前基因的数据（第2列匹配基因名）
	            gene_subset <- bed_data[bed_data$V2 == gene, ]
	            # 若有数据：计算第7列的均值（甲基化水平），添加到gene_data
	            if (nrow(gene_subset) > 0) {
			            mean_meth <- mean(gene_subset$V7, na.rm = TRUE)  # na.rm排除NA值
		            gene_data <- rbind(gene_data, data.frame(
								               Context = context,
									                 Sample = sample,
									                 Mean_Methylation = mean_meth
											         ))
			            has_data <- TRUE  # 标记当前基因有有效数据
			          }
		        }
    }
    
    # --------------------------------------------------------------------------
    # 子步骤1：处理文件缺失错误（若有bed文件缺失，终止程序并提示）
    # --------------------------------------------------------------------------
    if (length(error_files) > 0) {
	        stop(paste0(
			          "错误：以下bed文件缺失，请补充后重新运行：\n",
				        paste(error_files, collapse = "\n")
				      ))
      }
      
      # --------------------------------------------------------------------------
      # 子步骤2：当前基因无数据→跳过；有数据→生成图表
      # --------------------------------------------------------------------------
      if (!has_data) {
	          skipped_genes <- c(skipped_genes, gene)
          cat("→ 基因", gene, "：无数据，已跳过\n")
	    } else {
		        valid_genes <- c(valid_genes, gene)
	      
	      # 数据预处理：因子化Context和Sample，确保绘图顺序正确（按设定的contexts/samples顺序）
	      gene_data$Context <- factor(gene_data$Context, levels = contexts)
	          gene_data$Sample <- factor(gene_data$Sample, levels = samples)
	          
	          # ------------------------------------------------------------------------
	          # 绘制单个基因的柱状图（ggplot2，科研级风格）
	          # ------------------------------------------------------------------------
	          p <- ggplot(gene_data, aes(x = Context, y = Mean_Methylation, fill = Sample)) +
			        # 柱状图核心：黑色细边框（linewidth=0.2）、分组排列（dodge）、合适宽度
			        geom_col(
					         position = position_dodge(width = 0.8),
						         width = 0.7,
						         color = "black",
							         linewidth = 0.2
							       ) +
				      # 配色：使用预设的sample_colors，图例标签用sample_labels
				      scale_fill_manual(
							        values = sample_colors,
								        labels = sample_labels
								      ) +
				            # Y轴范围：从0开始（甲基化水平不能为负），顶部留15%空白（避免柱子顶到图边）
				            scale_y_continuous(
							               expand = c(0, 0),
								               limits = c(0, max(gene_data$Mean_Methylation, na.rm = TRUE) * 1.15)
								             ) +
					          # 主题：简洁风格，去除多余网格线，调整字体大小（适合期刊）
					          theme_bw() +
						        theme(
							              # 面板设置：去除网格线，添加浅灰色边框
							              panel.grid = element_blank(),
								              panel.border = element_rect(color = "gray70", linewidth = 0.2),
								              # 坐标轴：字体大小11，X轴标签水平（不旋转）
								              axis.text = element_text(size = 11, color = "black"),
									              axis.text.x = element_text(angle = 0, hjust = 0.5),
									              axis.title = element_text(
														          size = 13,
															            color = "black",
															            face = "bold",
																              margin = margin(t = 0, r = 10, b = 0, l = 0)  # 标题与轴间距
																            ),
							              # 图例：放在右侧，字体大小11，标题加粗
							              legend.position = "right",
								              legend.title = element_text(size = 12, color = "black", face = "bold"),
								              legend.text = element_text(size = 11, color = "black"),
									              legend.key.width = unit(0.8, "cm"),  # 统一图例宽度
									              legend.background = element_rect(fill = NA),  # 透明背景
										              # 图标题：突出基因名，居中，加粗
										              plot.title = element_text(
															          size = 14,
																            color = "black",
																            face = "bold",
																	              hjust = 0.5,
																	              margin = margin(b = 15)  # 标题与图间距
																		              )
										            ) +
							      # 图表标签：X轴、Y轴、图例标题、图标题
							      labs(
								           x = axis_x,
									           y = axis_y,
									           fill = "Sample Group",  # 图例标题
										           title = paste0(plot_title_prefix, gene)  # 图标题（含当前基因名）
										         )
							          
							          # ------------------------------------------------------------------------
							          # 保存图表：PNG（300DPI，适合展示）和PDF（矢量图，适合论文发表）
							          # ------------------------------------------------------------------------
							          # 构建输出文件名（基因名+后缀，如：ANR1_methylation.png）
							          png_path <- paste0(gene, "_methylation.png")
							          pdf_path <- paste0(gene, "_methylation.pdf")
								      
								      # 保存PNG（高分辨率，背景白色）
								      ggsave(
									           png_path,
										         plot = p,
										         width = 6,    # 固定宽度6英寸
											       height = 5,   # 固定高度5英寸
											       dpi = 300,    # 300DPI（印刷级分辨率）
											             bg = "white", # 白色背景（避免透明）
											             limitsize = FALSE
												         )
								      
								      # 保存PDF（矢量图，无分辨率限制，可无限放大）
								      ggsave(
									           pdf_path,
										         plot = p,
										         width = 6,
											       height = 5,
											       device = "pdf",
											             bg = "white",
											             limitsize = FALSE
												         )
								          
								          cat("→ 基因", gene, "：已生成 ", png_path, " 和 ", pdf_path, "\n", sep = "")
								        }
}


# ============================================================================
# 第五步：输出分析总结（统计有效/跳过基因数，便于核对）
# ============================================================================
cat("\n=== 分析完成 ===")
cat("\n1. 处理结果统计：")
cat("\n   - 总基因数：", length(gene_list))
cat("\n   - 有效基因数（已绘图）：", length(valid_genes))
cat("\n   - 跳过基因数（无数据）：", length(skipped_genes))

if (length(skipped_genes) > 0) {
	  cat("\n2. 跳过的基因列表：")
  cat("\n   ", paste(skipped_genes, collapse = ", "))
}

cat("\n3. 输出文件说明：")
cat("\n   - PNG文件：300DPI，适合PPT/报告展示")
cat("\n   - PDF文件：矢量图，适合期刊论文发表")
cat("\n   - 所有文件均保存在当前目录：", getwd(), "\n")##############################################################################
# 脚本名称：single_gene_methylation_plots.R
# 功能描述：读取9个甲基化bed文件，为每个有效基因单独生成甲基化水平柱状图
# 输出文件：每个有效基因对应1个PNG（300DPI）和1个PDF（矢量图），存于当前目录
# 依赖包：tidyverse（包含ggplot2、dplyr等绘图与数据处理包）
# 使用方式：1. 确保gene_list和9个bed文件在工作目录；2. Rscript single_gene_methylation_plots.R
##############################################################################

# ============================================================================
# 第一步：加载依赖包并检查
# ============================================================================
# 检查tidyverse是否已安装，未安装则报错提示
if (!requireNamespace("tidyverse", quietly = TRUE)) {
	  stop("错误：需要安装'tidyverse'包！请先运行：install.packages('tidyverse')")
}
# 加载包（suppressMessages隐藏加载时的冲突提示，避免干扰）
suppressMessages(library(tidyverse))


# ============================================================================
# 第二步：配置基础参数（样本、甲基化类型、配色等，便于后续修改）
# ============================================================================
# 1. 样本列表（R/L/S，与bed文件名前缀对应）
samples <- c("R", "L", "S")
# 2. 甲基化类型（CpG/CHG/CHH，与bed文件名后缀对应）
contexts <- c("CpG", "CHG", "CHH")
# 3. 科研级配色（红色/R、蓝色/L、绿色/S，区分度高，适合期刊图）
sample_colors <- c(
		     "R" = "#E74C3C",
		       "L" = "#3498DB",
		       "S" = "#2ECC71"
		       )
# 4. 图例标签（比样本名更完整，提升图的可读性）
sample_labels <- c(
		     "R" = "Sample R",
		       "L" = "Sample L",
		       "S" = "Sample S"
		       )
# 5. 图表固定标签（坐标轴、标题前缀）
axis_x <- "Methylation Context"       # X轴标签（甲基化类型）
axis_y <- "Mean Methylation Level"    # Y轴标签（平均甲基化水平）
plot_title_prefix <- "Methylation Level of Gene: "  # 图标题前缀


# ============================================================================
# 第三步：读取基因列表并初始化记录变量
# ============================================================================
# 读取gene_list文件（无表头，第一列即为基因名）
gene_list_path <- "gene_list"
if (!file.exists(gene_list_path)) {
	  stop(paste0("错误：基因列表文件'gene_list'不存在！请检查路径：", getwd()))
}
gene_list <- read.table(
			  gene_list_path,
			    header = FALSE,        # 无表头
			    sep = "\t",            # 制表符分隔（兼容纯文本列）
			      stringsAsFactors = FALSE  # 基因名按字符处理，不转因子
			    )$V1  # 取第一列作为基因列表

# 初始化记录变量（跟踪处理结果）
valid_genes <- c()       # 成功绘图的有效基因
skipped_genes <- c()     # 无数据的跳过基因
error_files <- c()       # 缺失的bed文件（用于报错）


# ============================================================================
# 第四步：遍历每个基因，提取数据并生成独立图表
# ============================================================================
cat("=== 开始分析：共", length(gene_list), "个基因 ===\n")

for (gene in gene_list) {
	  # 临时变量：存储当前基因的甲基化数据（Context/Sample/Mean_Methylation）
	  gene_data <- data.frame(
				      Context = character(),
				          Sample = character(),
				          Mean_Methylation = numeric(),
					      stringsAsFactors = FALSE
					    )
  # 标记：当前基因是否有有效数据
  has_data <- FALSE
    
    # 嵌套循环：遍历所有9个bed文件（3样本×3甲基化类型）
    for (sample in samples) {
	        for (context in contexts) {
			      # 构建当前bed文件路径（如：R.CpG.bed、L.CHG.bed）
			      bed_path <- paste0(sample, ".", context, ".bed")
        
        # 检查bed文件是否存在：不存在则记录，后续统一报错
        if (!file.exists(bed_path)) {
		        error_files <- unique(c(error_files, bed_path))
	        next  # 跳过当前文件，继续检查其他文件
		      }
	      
	      # 读取bed文件（按制表符分隔，确保列识别正确）
	      bed_data <- read.table(
				             bed_path,
					             header = FALSE,
					             sep = "\t",
						             quote = "",          # 忽略文件中的引号，避免列分裂错误
						             na.strings = "",     # 空值按NA处理
							             stringsAsFactors = FALSE
							           )
	      
	      # 检查bed文件列数：至少7列（需提取第2列基因名、第7列甲基化水平）
	      if (ncol(bed_data) < 7) {
		              cat("警告：文件", bed_path, "列数不足7列，跳过该文件\n")
	              next
		            }
	            
	            # 提取当前基因的数据（第2列匹配基因名）
	            gene_subset <- bed_data[bed_data$V2 == gene, ]
	            # 若有数据：计算第7列的均值（甲基化水平），添加到gene_data
	            if (nrow(gene_subset) > 0) {
			            mean_meth <- mean(gene_subset$V7, na.rm = TRUE)  # na.rm排除NA值
		            gene_data <- rbind(gene_data, data.frame(
								               Context = context,
									                 Sample = sample,
									                 Mean_Methylation = mean_meth
											         ))
			            has_data <- TRUE  # 标记当前基因有有效数据
			          }
		        }
    }
    
    # --------------------------------------------------------------------------
    # 子步骤1：处理文件缺失错误（若有bed文件缺失，终止程序并提示）
    # --------------------------------------------------------------------------
    if (length(error_files) > 0) {
	        stop(paste0(
			          "错误：以下bed文件缺失，请补充后重新运行：\n",
				        paste(error_files, collapse = "\n")
				      ))
      }
      
      # --------------------------------------------------------------------------
      # 子步骤2：当前基因无数据→跳过；有数据→生成图表
      # --------------------------------------------------------------------------
      if (!has_data) {
	          skipped_genes <- c(skipped_genes, gene)
          cat("→ 基因", gene, "：无数据，已跳过\n")
	    } else {
		        valid_genes <- c(valid_genes, gene)
	      
	      # 数据预处理：因子化Context和Sample，确保绘图顺序正确（按设定的contexts/samples顺序）
	      gene_data$Context <- factor(gene_data$Context, levels = contexts)
	          gene_data$Sample <- factor(gene_data$Sample, levels = samples)
	          
	          # ------------------------------------------------------------------------
	          # 绘制单个基因的柱状图（ggplot2，科研级风格）
	          # ------------------------------------------------------------------------
	          p <- ggplot(gene_data, aes(x = Context, y = Mean_Methylation, fill = Sample)) +
			        # 柱状图核心：黑色细边框（linewidth=0.2）、分组排列（dodge）、合适宽度
			        geom_col(
					         position = position_dodge(width = 0.8),
						         width = 0.7,
						         color = "black",
							         linewidth = 0.2
							       ) +
				      # 配色：使用预设的sample_colors，图例标签用sample_labels
				      scale_fill_manual(
							        values = sample_colors,
								        labels = sample_labels
								      ) +
				            # Y轴范围：从0开始（甲基化水平不能为负），顶部留15%空白（避免柱子顶到图边）
				            scale_y_continuous(
							               expand = c(0, 0),
								               limits = c(0, max(gene_data$Mean_Methylation, na.rm = TRUE) * 1.15)
								             ) +
					          # 主题：简洁风格，去除多余网格线，调整字体大小（适合期刊）
					          theme_bw() +
						        theme(
							              # 面板设置：去除网格线，添加浅灰色边框
							              panel.grid = element_blank(),
								              panel.border = element_rect(color = "gray70", linewidth = 0.2),
								              # 坐标轴：字体大小11，X轴标签水平（不旋转）
								              axis.text = element_text(size = 11, color = "black"),
									              axis.text.x = element_text(angle = 0, hjust = 0.5),
									              axis.title = element_text(
														          size = 13,
															            color = "black",
															            face = "bold",
																              margin = margin(t = 0, r = 10, b = 0, l = 0)  # 标题与轴间距
																            ),
							              # 图例：放在右侧，字体大小11，标题加粗
							              legend.position = "right",
								              legend.title = element_text(size = 12, color = "black", face = "bold"),
								              legend.text = element_text(size = 11, color = "black"),
									              legend.key.width = unit(0.8, "cm"),  # 统一图例宽度
									              legend.background = element_rect(fill = NA),  # 透明背景
										              # 图标题：突出基因名，居中，加粗
										              plot.title = element_text(
															          size = 14,
																            color = "black",
																            face = "bold",
																	              hjust = 0.5,
																	              margin = margin(b = 15)  # 标题与图间距
																		              )
										            ) +
							      # 图表标签：X轴、Y轴、图例标题、图标题

