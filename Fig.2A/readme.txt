2C：barplot.py  图的颜色和坐标不一致，提供的代码仅供参考
2A：02.methy_3group_circ_plot.R
2D：enrich_analysis_with_highlight.R	图的标记和显示通路不一致，提供的代码仅供参考
2B：test_plot.R			x轴和箱体表示的不一致，提供的代码仅供参考
2G、H：run.R			y轴刻度和图的颜色等不一致，提供的代码仅供参考
2I：TF-pie.py		        颜色,图例的位置不一致，提供的代码仅供参考
2J：代码如下			颜色和模块名称不一样，提供的代码仅供参考
filename = paste(prefix, "_module_trait.pdf", sep = "")
pdf(file = filename)
#sizeGrWindow(10,5)  #调整图片长宽
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

write.table(moduleTraitCor, file="moduleTraitCor.xls", quote=FALSE, sep="\t",col.names =NA)
write.table(moduleTraitPvalue, file="moduleTraitPvalue.xls", quote=FALSE, sep="\t",col.names =NA)

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.3,
               zlim = c(-1,1),
               cex.lab = 0.5,
               main = paste("Module-trait relationships"))
dev.off()			
2K：Correlation.pearson.R			颜色不一致，提供的代码仅供参考	
