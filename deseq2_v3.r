#http://folk.uio.no/jonbra/MBV-INF4410_2017/exercises/2017-12-07_R_DESeq2_exercises_without_results.html - got script outline/figure ideas from this
#https://github.com/hbctraining/DGE_workshop/tree/master/lessons - much clearer understanding from this.
#############################################
##############   PACKAGES   #################
#############################################
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#install.packages("pheatmap")
#install.packages("reshape2")
#install.packages("gplots")
#install.packages("RColorBrewer")
#install.packages("ggplot2")
#install.packages('cowplot')
#install.packages('tidyverse')
#install.packages('ggrepel')
#install.packages('DEGreport')

library('tidyverse')
library("DESeq2")
library("pheatmap")
library("reshape2")
library("gplots")
library("RColorBrewer")
library("ggplot2")
library('gridExtra')
library('cowplot')
library('ggrepel')

#library('dplyr')
#library('tibble')
#######################################################
################## command line arguments #############
#######################################################
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Arguments required (order dependent): 1) count file, 2) condition file", call.=FALSE)
} else if (length(args)==1) {
# default output file
  args[3] = "out.txt"
}
#######################################################
setwd('/Users/ie3/Desktop/LMB_ONE_MONTH/rna_seq/deseq2_analysis/aerobic_no_contaminated_with_int_regions_analysis')
#######################################################

count_data = read.table(file = "BIR1049_aerobic_WT_and_KO_counts_with_int.txt", header = T, sep = "\t")
count_data = data.frame(count_data[,-1], row.names=count_data[,1])
col_data = read.table(file = "aerobic_conditions.txt", header = T, sep = "\t")

#######################################################
################### read in data ######################
#######################################################
#read in count data
count_data = read.table(file = args[1], header = T, sep = "\t")
#make gene names row names
count_data = data.frame(count_data[,-1], row.names=count_data[,1])
#read in metadata - the condition names must be same as column names
col_data = read.table(file = args[2], header = T, sep = "\t")


############################################
############ Quality Assessment ############ 
############################################
# requires data to be normalised - for the de analysis this is done within the deseq command
# convert count_data table to dataframe.
count_data = as.data.frame(count_data)
# ESSENTIALLY TO IMPROVE VISUALISATION we log-transform to make numbers on scale, improves the distances/clustering of the normalised counts. DESeq2 uses a regularized log transform (rlog) of the normalised counts for sample level QC as it moderates the variance across the mean, improving clustering.

############################################
############ NORMALIZATION #################
############################################
###### to compare normalised and unnormalised counts ######
not_norm_pseudoCount = log2(count_data + 1)
not_norm_pseudoCount <- cbind(Row.Names = rownames(not_norm_pseudoCount), not_norm_pseudoCount)
#making quality plots with ggplot: 
not_norm_df = melt(not_norm_pseudoCount, id=c('Row.Names')) # reshape the matrix 
# adding condition column
not_norm_df = data.frame(not_norm_df, Condition= substr(not_norm_df$variable, 9,10))
# distribution of raw read counts 
not_norm_count <- ggplot(not_norm_df, aes(x = not_norm_df$variable, y = value,  fill = Condition)) + geom_boxplot() + xlab("") + ylab(expression(log[2](count + 1))) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
not_norm_count
#Normalization is done by DESeq2. It will estimate a size factor (scaling factor) which all the genes in a sample will be multiplied with. Ideally the size factor should be 1, which means that no normalization will take place.
# creating DESeqDataSet
dds = DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ condition) # we're testing for the different condidtions
dds = estimateSizeFactors(dds)
sizeFactors(dds)
# plot normalized counts
norm_counts = counts(dds, normalized = TRUE) # Extract the normalized counts
pseudoCount = log2(norm_counts + 1) # convert to log-scale for visualization
df = melt(pseudoCount) # transpose the matrix
df = data.frame(df, Condition = substr(df$Var2, 9, 10))
norm_count <- ggplot(df, aes(x = df$Var2, y = value, fill = Condition)) + geom_boxplot() + xlab("") +ylab(expression(log[2](count + 1))) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
norm_count
# Create plot file 
pdf("aerobic_int_not_nomr_vs_norm_counts.pdf")
grid.arrange(not_norm_count, norm_count, ncol=2)
dev.off()
# store for normalised counts for future use
write.table(norm_counts, file="aerobic_int_normalized_counts.txt", sep="\t", quote=F, col.names=NA)
####################################################
#########     Quality Control with     #############
####### Unsupervised Clustering Methods ############
####################################################

#rlog transform deseq normalised read counts for visualisation
rld = rlogTransformation(dds)

# pca plot - do the samples with the same conditions cluster.
pca_plot <- plotPCA(rld, intgroup=c("condition"))
pca_plot
# Calculate distances using transformed (and normalized) counts
distsRL <- dist(t(assay(rld))) 
mat <- as.matrix(distsRL) # convert to matrix
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition)) # set rownames in the matrix
colnames(mat) = NULL # remove column names
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) # set colors
dist_mat <- pheatmap(mat,clustering_distance_rows=distsRL,clustering_distance_cols=distsRL,col=colors)
dist_mat
# create quality plots #
pdf("aerobic_int_quality_plots.pdf")
grid.arrange(grobs=list(pca_plot, dist_mat[[4]]))
dev.off()

########################################################
############## Differential Expression #################
#############           Analysis       #################
########################################################
# running differential expression analysis
#Everything from normalization to linear modeling was carried out by the use of a single function
dds <- DESeq(dds)
# plot dispersion estimates useful to determine how well deseq2 model fits. #
plotDispEsts(dds)

#building results table:
contrast <- c('condition','WT', 'KO')
# alpha is the significance thrsehold - default is 0.1.
#results table with unshrunken values
res_table_unshrunken <- results(dds, contrast=contrast, alpha=0.05,lfcThreshold=2)
# lcf shrinkage is to correct for high dispersion of some genes. coef is to tell lfc to shrink the values associated with this contrast.
res_table_shrunken <- lfcShrink(dds, conntrast=contrast, res=res_table_unshrunken, coef='condition_WT_vs_KO', lfcThreshold=2)

#MA plot shows the mean of the normalised counts versus the log2fold changes for all genes tested. plot also allows us to evaluate the magnitude of fold changes and how they are distributed to mean expression.
#unshrunken
MA_unshrunk <- plotMA(res_table_unshrunken, ylim=c(-3,3))
#shrunken
MA_shrunk <- plotMA(res_table_shrunken, ylim=c(-3,3)) 

#save ma plots to file 
#pdf('aerobic_MA_plots.pdf')
#grid.arrange(grobs=list(MA_shrunk, MA_unshrunk), newpage = T)
#dev.off()

#summarise results
# unshrunken
summary(res_table_unshrunken)
# shrunken
summary(res_table_shrunken)

#extract significantly differentially expressed genes #
# variables with cirterial
padj.cutoff <- 0.05

#get unshrunken results table:
res_table_us_tb <- res_table_unshrunken %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble()
# filter based on p=value and log2fold change
sig_res_table_us <- res_table_us_tb %>% filter(padj < padj.cutoff & abs(log2FoldChange) > 2)

#get shrunken results table
res_table_s_tb <- res_table_shrunken %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble()
# filter based on p=value and log2fold change
sig_res_table_s <- res_table_s_tb %>% filter(padj < padj.cutoff & abs(log2FoldChange) > 2)

#write tables to file
write.table(res_table_us_tb, file='aerobic_int_us_all_de_stats.txt', sep='\t', row.names=F, quote=F)
write.table(sig_res_table_us, file='aerobic_int_us_signif_de_stats.txt', sep='\t', row.names=F, quote=F)
write.table(res_table_s_tb, file='aerobic_int_s_all_de_stats.txt', sep='\t', row.names=F, quote=F)
write.table(sig_res_table_s, file='aerobic_int_s_signif_de_stats.txt', sep='\t', row.names=F, quote=F)

#plotCounts(dds, gene="Mycobacterium_abscessus_massiliense_BIR1049_v0.1_01187", intgroup="condition")

#################################################
#############  VISUALISATION   ##################
##################################################
# tibbles for metadata and normalized counts, apparently 
metadata <- col_data %>% rownames_to_column(var="samplename") %>% as_tibble()
normalized_counts <- norm_counts %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble()

########## creating heatmap ###############
# create table for significant genes for heatmap
norm_sig_s <- normalized_counts[,c(1,2:6)] %>% filter(gene %in% sig_res_table_s$gene) %>% data.frame() %>% column_to_rownames(var = "gene")


scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}



write.table(norm_sig_s,file = 'normalised_counts_with_int.txt')
#annotate the heatmap
annotation <- metadata %>% select(sample, condition) %>% data.frame(row.names = "sample")
#set a color palette
heat_colors <- brewer.pal(9, "Blues")
breaksList = seq(-2, 3, by = 0.25)
#heatmap
pheatmap(norm_sig_s, color = colorRampPalette(brewer.pal(n=9, name='Blues'))(length(breaksList)), breaks= breaksList, cluster_rows = T, show_rownames = T, annotation = annotation, border_color = NA,fontsize = 10, scale = 'row',fontsize_row = 10,height = 20)


#volcano plot
## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 2 in either direction
res_table_s_tb <- res_table_s_tb %>% mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 2)

## Volcano plot
ggplot(res_table_s_tb) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold, alpha=0.1)) +
  ggtitle("aerobic dpnM presence vs absence") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,400)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))


