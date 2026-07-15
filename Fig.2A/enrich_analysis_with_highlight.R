args=commandArgs(TRUE)
if (length(args) < 1){
    print ("usage: all <DE_results> <GO2id> <GO2name> <KO2id> <KO2name> <DiffGroup> <Pvalue> <Qvalue>")
    print ("usage: kegg <DE_results> <KO2id> <KO2name> <DiffGroup> <Pvalue> <Qvalue>")
    print ("usage: go <DE_results> <GO2id> <GO2name> <DiffGroup> <Pvalue> <Qvalue>")
    q()
}

if (args[1] %in% c('-h', '-help', '--h', '--help')){
    print ("usage: all <DE_results> <GO2id> <GO2name> <KO2id> <KO2name> <DiffGroup> <Pvalue> <Qvalue>")
    print ("usage: kegg <DE_results> <KO2id> <KO2name> <DiffGroup> <Pvalue> <Qvalue>")
    print ("usage: go <DE_results> <GO2id> <GO2name> <DiffGroup> <Pvalue> <Qvalue>")
    q()
}


library(ggplot2)
library(stringr)
library(clusterProfiler)  # 移到前面统一加载

# 辅助函数：确保特定通路包含在要展示的结果中
ensure_pathways_included <- function(enrich_result, specific_ids, top_n = 20) {
    if (is.null(enrich_result) || nrow(enrich_result) == 0) {
        return(enrich_result)
    }
    
    # 获取top N结果
    top_results <- head(enrich_result, top_n)
    
    # 检查特定通路是否已在top结果中
    missing_ids <- setdiff(specific_ids, top_results$ID)
    
    # 如果有缺失的特定通路，且它们存在于完整结果中，则添加进来
    if (length(missing_ids) > 0) {
        # 从完整结果中找到这些特定通路
        specific_entries <- enrich_result[enrich_result$ID %in% missing_ids, ]
        
        if (nrow(specific_entries) > 0) {
            # 合并top结果和特定通路
            combined <- rbind(top_results, specific_entries)
            # 去重（以防特定通路恰好也在top结果中）
            combined <- combined[!duplicated(combined$ID), ]
            return(combined)
        }
    }
    
    return(top_results)
}

go_enrich <- function(subset, TERM2GENE_file, TERM2NAME_file, file_prefix, pvalue, qvalue){
    go_anno <- read.delim(TERM2GENE_file, colClasses = 'character', header = T, stringsAsFactors = FALSE)[,1:2]
    names(go_anno) <- c('gene_id', 'go_id')
    go_class <- read.delim(TERM2NAME_file, colClasses = 'character', header = FALSE, stringsAsFactors = FALSE)
    names(go_class) <- c('go_id', 'go_description', 'go_class')
    diff_gene <- read.delim(subset, header = FALSE, stringsAsFactors = FALSE)
    gene_select <- diff_gene$V1
    subset_go <- merge(diff_gene, go_anno, by.x = "V1", by.y = "gene_id")
    
    if ((nrow(subset_go)) > 0){
        go_rich_1 <- enricher(gene = gene_select,
            TERM2GENE = go_anno[,c('go_id', 'gene_id')],
            TERM2NAME = go_class[,c('go_id', 'go_description')],
            pvalueCutoff = 1,
            qvalueCutoff = 1)
        go_rich <- enricher(gene = gene_select,
            TERM2GENE = go_anno[,c('go_id', 'gene_id')],
            TERM2NAME = go_class[,c('go_id', 'go_description')],
            pvalueCutoff = pvalue,
            qvalueCutoff = qvalue)
        
        go_rich_1_rs <- as.data.frame(go_rich_1)
        write.table(go_rich_1_rs, paste(file_prefix, '.go_enrich.xls', sep=""), sep = '\t', row.names = FALSE, quote = FALSE)
        
        # 确保GO:0009813通路被包含在要展示的结果中
        go_rich_df <- as.data.frame(go_rich)
        if (nrow(go_rich_df) > 0) {
            specific_go <- "GO:0009813"
            display_go <- ensure_pathways_included(go_rich_df, specific_go, 20)
            
            # 转换回enrichResult对象以便使用clusterProfiler的绘图函数
            go_rich_display <- go_rich
            go_rich_display@result <- display_go
            
            # barplot
            bar_plot <- barplot(go_rich_display, drop=TRUE, color = "p.adjust", showCategory=nrow(display_go)) + 
                        scale_y_discrete(labels = function(x) str_wrap(x, width = 35))
            name <- paste(file_prefix, ".go_barplot.pdf", sep="")
            pdf(name, width=12, height=12)
            print(bar_plot)
            dev.off()  # 确保关闭设备
            
            name <- paste(file_prefix, ".go_barplot.png", sep="")
            png(name,width=12*600,height=12*600, res=600)
            print(bar_plot)
            dev.off()
            
            # dotplot
            dot_plot <- dotplot(go_rich_display, color = "p.adjust", showCategory=nrow(display_go)) + 
                        scale_y_discrete(labels = function(x) str_wrap(x, width = 35))
            name <- paste(file_prefix, ".go_dotplot.pdf", sep="")
            pdf(name, width=12, height=12)
            print(dot_plot)
            dev.off()
            
            name <- paste(file_prefix, ".go_dotplot.png", sep="")
            png(name,width=12*600,height=12*600, res=600)
            print(dot_plot)
            dev.off()
        } else {
            print(paste(file_prefix, " no significant GO enrichment results", sep = ''))
        }
    } else {
        warn <- paste(file_prefix, ' no go annotation', sep = '')
        print(warn)
    }
}

kegg_enrich <- function(subset, TERM2GENE_file, TERM2NAME_file, file_prefix, pvalue, qvalue){
    kegg_anno <- read.delim(TERM2GENE_file, colClasses = 'character', header = T, stringsAsFactors = FALSE)[,1:2]
    names(kegg_anno) <- c('gene_id', 'pathway')
    kegg_class <- read.delim(TERM2NAME_file, colClasses = 'character', header = FALSE, stringsAsFactors = FALSE)
    names(kegg_class) <- c('class_A', 'class_B', 'pathway', 'pathway_description')
    diff_gene <- read.delim(subset, header = FALSE, stringsAsFactors = FALSE)
    gene_select <- diff_gene$V1
    subset_ko <- merge(diff_gene, kegg_anno, by.x = "V1", by.y = "gene_id")
    
    if ((nrow(subset_ko)) > 0){
        kegg_rich_1 <- enricher(gene = gene_select,
            TERM2GENE = kegg_anno[,c('pathway', 'gene_id')],
            TERM2NAME = kegg_class[,c('pathway', 'pathway_description')],
            pvalueCutoff = 1,
            qvalueCutoff = 1)
        kegg_rich <- enricher(gene = gene_select,
            TERM2GENE = kegg_anno[,c('pathway', 'gene_id')],
            TERM2NAME = kegg_class[,c('pathway', 'pathway_description')],
            pvalueCutoff = pvalue,
            qvalueCutoff = qvalue)
        
        kegg_rich_1_rs <- as.data.frame(kegg_rich_1)
        write.table(kegg_rich_1_rs, file = paste(file_prefix, '.ko_enrich.xls', sep = ''), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        
        # 确保ko00940通路被包含在要展示的结果中
        kegg_rich_df <- as.data.frame(kegg_rich)
        if (nrow(kegg_rich_df) > 0) {
            specific_kegg <- "ko00940"
            display_kegg <- ensure_pathways_included(kegg_rich_df, specific_kegg, 20)
            
            # 转换回enrichResult对象以便使用clusterProfiler的绘图函数
            kegg_rich_display <- kegg_rich
            kegg_rich_display@result <- display_kegg
            
            # barplot
            bar_plot <- barplot(kegg_rich_display, drop=TRUE, color = "p.adjust", showCategory=nrow(display_kegg)) + 
                        scale_y_discrete(labels = function(x) str_wrap(x, width = 35))
            name <- paste(file_prefix, ".ko_barplot.pdf", sep="")
            pdf(name, width=12, height=12)
            print(bar_plot)
            dev.off()
            
            name <- paste(file_prefix, ".ko_barplot.png", sep="")
            png(name,width=12*600,height=12*600, res=600)
            print(bar_plot)
            dev.off()
            
            # dotplot
            dot_plot <- dotplot(kegg_rich_display, color = "p.adjust", showCategory=nrow(display_kegg)) + 
                        scale_y_discrete(labels = function(x) str_wrap(x, width = 35))
            name <- paste(file_prefix, ".ko_dotplot.pdf", sep="")
            pdf(name,width=12,height=12)
            print(dot_plot)
            dev.off()
            
            name <- paste(file_prefix, ".ko_dotplot.png", sep="")
            png(name,width=12*600,height=12*600, res=600)
            print(dot_plot)
            dev.off()
        } else {
            print(paste(file_prefix, " no significant KEGG enrichment results", sep = ''))
        }
    } else {
        warn <- paste(file_prefix, ' no ko annotation', sep = '')
        print(warn)
    }
}

if (args[1] == 'all'){
    if (length(args) != 9) {
        print ("usage: all <DE_results> <GO2id> <GO2name> <KO2id> <KO2name> <DiffGroup> <Pvalue> <Qvalue>")
        q()
    }else{
        run_type <- as.character(args[1])
        DE_result <- args[2]
        GO2id <- args[3]
        GO2name <- args[4]
        KO2id <- args[5]
        KO2name <- args[6]
        prefix <- as.character(args[7])
        pvalue <- as.double(args[8])
        qvalue <- as.double(args[9])
        
        go_enrich(DE_result,GO2id,GO2name,prefix,pvalue,qvalue)
        kegg_enrich(DE_result,KO2id,KO2name,prefix,pvalue,qvalue)
    }
}else if (args[1] == 'kegg'){
    if (length(args) != 7) {
        print ("usage: kegg <DE_results> <KO2id> <KO2name> <DiffGroup> <Pvalue> <Qvalue>")
        q()
    }else{
        run_type <- as.character(args[1])
        DE_result <- args[2]
        KO2id <- args[3]
        KO2name <- args[4]
        prefix <- as.character(args[5])
        pvalue <- as.double(args[6])
        qvalue <- as.double(args[7])
        
        kegg_enrich(DE_result,KO2id,KO2name,prefix,pvalue,qvalue)
    }
}else if (args[1] == 'go'){
    if (length(args) != 7) {
        print ("usage: go <DE_results> <GO2id> <GO2name> <DiffGroup> <Pvalue> <Qvalue>")
        q()
    }else{
        run_type <- as.character(args[1])
        DE_result <- args[2]
        GO2id <- args[3]
        GO2name <- args[4]
        prefix <- as.character(args[5])
        pvalue <- as.double(args[6])
        qvalue <- as.double(args[7])
        
        go_enrich(DE_result,GO2id,GO2name,prefix,pvalue,qvalue)
    }
}else{
    print ("usage: all <DE_results> <GO2id> <GO2name> <KO2id> <KO2name> <DiffGroup> <Pvalue> <Qvalue>")
    print ("usage: kegg <DE_results> <KO2id> <KO2name> <DiffGroup> <Pvalue> <Qvalue>")
    print ("usage: go <DE_results> <GO2id> <GO2name> <DiffGroup> <Pvalue> <Qvalue>")
    q()
}
