
###Load this layout function before running the script for the first time
theme_Publication_blank <- function(base_size=10, base_family="") { #12 For ALDR paper
  require(grid)
  require(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA),
            panel.border = element_rect(colour = NA, fill = "transparent"),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,margin=margin(0,10,0,0)),
            axis.title.x = element_text(margin=margin(10,0,0,0)),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(size = 0.3),
            axis.line.x = element_line(size = 0.3, linetype = "solid", colour = "black"),
            axis.line.y = element_line(size = 0.3, linetype = "solid", colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA, fill="transparent"),
            legend.position = "bottom",
            #legend.direction = "horizontal",
            #legend.box = "horizontal",
            #legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(10, "pt"),
            #legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#d8d8d8",fill="#d8d8d8"),
            strip.text = element_text(face="bold")
    ))
  
} 

setwd('~/Desktop/PF4 Nature Rev Data/OMICS/bulk RNAseq/')

library('ggplot2')
library('factoextra')
library("genefilter")
library(DESeq2)
library(scales)
library(ggrepel)
library(dplyr)
library(pheatmap)
library(data.table)

#Load old dataset
load('dds_Hip_PF4.bin')

#### I) PCA analysis in volcano plots of the DEG's

##Draw PCA plot of all samples on the genes chaging with aging 
dds_temp <- dds_PF4_list$Hip

select <- results_list_PF4$Male$Old_Veh_vs_Young_Veh$ressig$gene_symbol
vsd <- varianceStabilizingTransformation(dds_temp)       #Variance stabilizing transformation

groups <- as.factor(dds_temp$age_treatment)
res.pca <-prcomp(t(assay(vsd)[select,]), scale = F)
myplot <- fviz_pca_ind(res.pca,
                       col.ind = groups, # color by groups
                       addEllipses = T, # Concentration ellipses
                       ellipse.type = "confidence",
                       legend.title = "Groups",
                       repel = T, label="none", axes = c(1,2)
)
myplot

#### II) Aging DEGs, PF4 trt DEGs, and Cognitive Signatures

##Aging DEGs
#Volcanoplot for DEG results - Aging vs Young 
plot_df <- as.data.frame(results_list_PF4$Male$Old_Veh_vs_Young_Veh$resall)
plot_df <- plot_df[order(plot_df$padj, decreasing = F),]
myplot <- ggplot(plot_df, aes(x = log2FoldChange, y=-log10(padj), color=ifelse(padj<0.05, ifelse(log2FoldChange>0, 'Upregulated', 'Downregulated'), 'n.s.'), label=gene_symbol)) + 
  geom_point() + theme_Publication_blank() +
  scale_color_manual(values = c('Upregulated'='#ea5430', 'Downregulated'='#6181d1',  'n.s.'='lightgrey')) + scale_x_continuous(limits = c(-2,2), oob = squish, expand = c(0,0))  +
  geom_text_repel(data= head(plot_df,10), color='black', point.size=NA, max.overlaps = Inf) # highlight the most significant genes
  
myplot

plot_df

#Save aging DEGs Aging vs Young to csv file
library(data.table)
plot_df <- as.data.frame(results_list_PF4$Male$Old_Veh_vs_Young_Veh$resall)
plot_df <- plot_df[order(plot_df$padj, decreasing = F),]
data_to_write_out <-plot_df 
fwrite(x= data_to_write_out, file ="670_aging_deg.csv")

#Heatmap for DEGs Aging vs Young 
plot_df <- as.data.frame(results_list_PF4$Male$Old_Veh_vs_Young_Veh$resall)
plot_df <- plot_df[order(plot_df$padj, decreasing = F),]
data <-as.data.frame(results_list_PF4$Male$Old_Veh_vs_Young_Veh$ressig %>% filter(padj < 0.05))$gene_symbol
dds_test <- dds_PF4_list$Hip[, dds_PF4_list$Hip$treatment == "Veh"]
vsd <- varianceStabilizingTransformation(dds_test)  

mat  <- assay(vsd)[data, ]
mat <- as.data.frame(t(scale(t(mat))))

anno <- as.data.frame(colData(vsd)[, 'age', drop=F])
mat2<-mat[,c(5,14:25,1:4,6:13)]
colnames(anno) <-'age'

annot_colors=list(age=c(Young="#3CC1C8", Old="#37658B"))

pheatmap(mat2, annotation_col = anno, show_colnames = F, cluster_rows = T, cluster_cols =F, annotation_colors=annot_colors, border_color="NA")


##PF4 trt DEGs
#Volcanoplot for DEG restuls - Veh vs PF4 in Aging 
plot_df <- as.data.frame(results_list_PF4$Male$Old_PF4_vs_Old_Veh$resall)
plot_df <- plot_df[order(plot_df$padj, decreasing = F),]

myplot <- ggplot(plot_df, aes(x = log2FoldChange, y=-log10(padj), color=ifelse(padj<0.05, ifelse(log2FoldChange>0, 'Upregulated', 'Downregulated'), 'n.s.'), label=gene_symbol)) + 
  geom_point() + theme_Publication_blank() +
  scale_color_manual(values = c('Upregulated'='#ea5430', 'Downregulated'='#6181d1',  'n.s.'='lightgrey')) + scale_x_continuous(limits = c(-2,2), oob = squish, expand = c(0,0))  +
  geom_text_repel(data = subset(plot_df, padj<0.05), color='black', point.size=NA, max.overlaps = Inf) # highlight the significant genes

myplot

#Save 24 pf4 edg to csv
library(data.table)
data_to_write_out <-plot_df 
fwrite(x= data_to_write_out, file ="24_pf4_deg.csv")

#Heatmap of PF4 trt DEGs
plot_df <- as.data.frame(results_list_PF4$Male$Old_PF4_vs_Old_Veh$resall)
plot_df <- plot_df[order(plot_df$padj, decreasing = F),]
data <-as.data.frame(results_list_PF4$Male$Old_PF4_vs_Old_Veh$ressig %>% filter(padj < 0.05))$gene_symbol
dds_test <- dds_PF4_list$Hip[, dds_PF4_list$Hip$age == "Old"]
vsd <- varianceStabilizingTransformation(dds_test)  

mat  <- assay(vsd)[data, ]
mat <- as.data.frame(t(scale(t(mat))))

anno <- as.data.frame(colData(vsd)[, 'treatment', drop=F])
mat <-mat[, order(anno$treatment)]
mat2<- mat[,c(11:22, 1:10)]
colnames(anno) <-'treatment'

annot_colors=list(treatment=c(Veh="grey", PF4="#9e1f63"))

pheatmap(mat2, annotation_col = anno, show_colnames = F, cluster_rows = T, cluster_cols =F, annotation_colors=annot_colors, border_color="NA")

##Cognitive signatures
#Analysis for WorkingMemory
load('YMaze_WorkingMemory_CorrelationResults.bin')
load('dds_Hip_PF4.bin')

dds_tissue <- dds_PF4_list$Hip

output_table <- counts(dds_tissue, normalized=T)
output_table <- as.data.frame(output_table)

#Get the top 500 genes
genes_to_check <- head(na.omit(resOrdered_Ymaze_HighLowDEGs), 500)$gene_symbol
Score_to_compare <- as.numeric(as.character(as.data.frame(colData(dds_tissue))[,'YoungAging_WorkingMemory']))

#Perform gene-vs-cogniton score correlataion
correlation_frame <- data_frame(gene='A', cor_spear=1, p_spear=1, cor_pears=1, p_pears=1)
i <- 1
for (gene in genes_to_check) {
  cor_frame <- data.frame(expr=t(output_table[gene,]), concentration= Score_to_compare)
  colnames(cor_frame) <- c('expr', 'concentration')
  cor_frame <- cor_frame %>% filter(expr > 0)
  cor_vector_spearman <- cor.test(cor_frame$expr, cor_frame$concentration, method = 'spearman')
  cor_vector_pearson <- cor.test(cor_frame$expr, cor_frame$concentration, method = 'pearson')
  correlation_frame[i, 'gene'] <- gene
  correlation_frame[i, 'cellsExpr'] <- nrow(cor_frame)
  correlation_frame[i, 'cor_spear'] <- cor_vector_spearman$estimate
  correlation_frame[i, 'p_spear'] <- cor_vector_spearman$p.value
  correlation_frame[i, 'cor_pears'] <- cor_vector_pearson$estimate
  correlation_frame[i, 'p_pears'] <- cor_vector_pearson$p.value
  i <- i + 1
}
#Cleanup results and run the multiple testing correction
correlation_frame <- na.omit(correlation_frame)
correlation_frame$padj_spear <- p.adjust(correlation_frame$p_spear, method = 'BH')
correlation_frame$padj_pear <- p.adjust(correlation_frame$p_pears, method = 'BH')
#Relabel
WorkingMemory_CorrelationGenes_fullTable <- correlation_frame
rm(correlation_frame)
#Inspect the results and select the top correlating genes
dim(WorkingMemory_CorrelationGenes_fullTable %>% filter(padj_pear < 0.1))
#Filter for high correlates
WorkingMemory_CorrelationGenes_topCorrelates <- (WorkingMemory_CorrelationGenes_fullTable %>% filter(padj_pear < 0.1)  %>% filter(abs(cor_spear) > .55))
#Create a transcription score
library(VISION)
color_code=c("Veh"="grey", "PF4"="#ea5430")
#Set the type of correlation you want to filter on 
correaltion_type <- 'pearson' #pearson; spearman
correaltion_column <- 'cor_pears' #cor_pears; cor_spear

dds_temp <- dds_PF4_list$Hip

signature_vector <- WorkingMemory_CorrelationGenes_topCorrelates %>% pull(correaltion_column)
names(signature_vector) <- toupper(WorkingMemory_CorrelationGenes_topCorrelates %>% pull('gene'))
#Creat signature vector
signature_vector[which(signature_vector<0)] <- -1
signature_vector[which(signature_vector>0)] <- 1
sigData <- signature_vector
sig <- createGeneSignature(name = "Young_LearningMemory", sigData = sigData)
#Run varianceStabilizedTransformation to obtain normalised counts
vsd <- varianceStabilizingTransformation(dds_temp)
scale_counttable <- assay(vsd)
row.names(scale_counttable) <- toupper(row.names(scale_counttable))

pos_sig_genes <- signature_vector[which(signature_vector>0)]
neg_sig_genes <- signature_vector[which(signature_vector<0)]
#Run for each sample the score calculation and store in a table
sample_table <- data.frame(sample='A', score=1, stringsAsFactors = F)
i <- 1
for (sample in colnames(scale_counttable)) {
  temp_df <- scale_counttable[,sample,drop=F]
  sample_score <- (sum(temp_df[names(signature_vector),] * signature_vector)) / length(signature_vector)
  sample_table[i,'sample'] <- as.character(sample)
  sample_table[i, "score"] <- sample_score
  i <- i + 1
}

plot_df <- sample_table
#Extract some metadata from the sample names
plot_df$age <- unlist(lapply(strsplit(as.character(plot_df$sample), split = "_"), function(x) x[[1]]))
plot_df$treatment <- unlist(lapply(strsplit(as.character(plot_df$sample), split = "_"), function(x) x[[2]]))
plot_df$SeqID <- plot_df$sample
plot_df$Learning_metric <- plyr::mapvalues(plot_df$SeqID, 
                                           from = as.character(as.data.frame(colData(dds_temp))$seqID),
                                           to = as.numeric(as.character(as.data.frame(colData(dds_temp))[,'YoungAging_WorkingMemory']))
                                           
)
plot_df$Learning_metric <- as.numeric(plot_df$Learning_metric)

#Plot the scatteer plot of score vs cognition 
(myplot <- ggplot(plot_df, aes(x=score, y=Learning_metric)) + geom_point(aes(shape=treatment, color=age)) + geom_smooth(method = 'lm') +
    theme_Publication_blank() + scale_shape_manual(values = c('Veh'=16, 'PF4'=3)) #+ scale_color_manual(values=c("Veh"="grey", "PF4"="#ea5430"))
)

#Test for correlation
plot_df_Transcription_score <- plot_df
cor.test(plot_df_Transcription_score$Learning_metric, plot_df_Transcription_score$score)

#Plot heatmap of correlate genes
topCorGenes <- WorkingMemory_CorrelationGenes_topCorrelates$gene
dds_test <- dds_PF4_list$Hip

#Plot heatmap of top 20 variable genes acrosss whole dataset
vsd <- varianceStabilizingTransformation(dds_test)       #Variance stabilizing transformation

mat  <- assay(vsd)[ topCorGenes, ]
mat <- as.data.frame(t(scale(t(mat))))
#mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("YoungAging_WorkingMemory", 'treatment', 'age', 'AgeMonths'), drop=F])
anno[,1] <- as.numeric(as.character(anno[,1]))
anno[,4] <- as.numeric(as.character(anno[,4]))

colnames(anno) <- c('score', 'treatment', 'age', 'Months')
anno$RNAscore <- plyr::mapvalues(row.names(anno), from = plot_df_Transcription_score$sample, plot_df_Transcription_score$score)
anno$RNAscore <- as.numeric(anno$RNAscore)
anno$MontsAge <- anno$Months
anno$Months <- NULL
mat <- mat[,order(anno$score)]
anno <- anno[order(anno$score),, drop=F]

#Plot heaatmap ordered by increasing cognition score
pheatmap(mat, annotation_col = anno, show_colnames = T, cluster_cols =F)

#Save 44 cognitive signatures to csv
library(data.table)
data_to_write_out <- mat
fwrite(x= data_to_write_out, file ="44_cognitive sinature.csv")


### III)venn diagram and box plots of normalized counts
#Plot Venn diagram
load('YMaze_WorkingMemory_CorrelationResults.bin')
library(Vennerable)
Vstem <- Venn( list('GenesCorrelatingWithCognition'=WorkingMemory_CorrelationGenes_topCorrelates$gene,
                    'AgingDEGs'=(results_list_PF4$Male$Old_Veh_vs_Young_Veh$ressig %>% filter(padj < 0.05))$gene_symbol,
                    'PF4DEGs'=(results_list_PF4$Male$Old_PF4_vs_Old_Veh$ressig %>% filter(padj < 0.05))$gene_symbol
)

)
plot(Vstem, doWeights=T)


### IV)Box plots of normalized counts
## Normalized counts from below are all Log2 transformed, and then plotted and analyzed.
#Akap11
gene <- 'Akap11'
data <- plotCounts(dds_PF4_list$Hip, gene=gene, intgroup=c("age_treatment"), returnData=T,normalized = T)

data$age_treatment <- factor(data$age_treatment,
                             levels=c('Young_Veh', 'Young_PF4',
                                      'Old_Veh', 'Old_PF4'), ordered = T)

data$age <- unlist(lapply(strsplit(as.character(row.names(data)), split = "_"), function(x) x[[1]]))
data$treatment <- unlist(lapply(strsplit(as.character(row.names(data)), split = "_"), function(x) x[[2]]))
data$age <- factor(data$age, levels = c('Young', 'Old'))
data$treatment <- factor(data$treatment, levels = c('Veh', 'PF4'))

(myplot <- ggplot(data, aes(x=age, y=count, color=treatment)) +
    geom_boxplot(outlier.colour = NA, notch = F) + xlab("") + ylab("Normalised counts")  +
    # geom_errorbar(aes(ymin=count-se, ymax=count+se),
    #               width=.2,                    # Width of the error bars
    #               position=position_dodge(.9)) +  
    theme_Publication_blank() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank()) + 
    geom_point(data=data, position = position_jitterdodge()) +
    scale_color_manual(values = c("Veh"="grey", "PF4"="#ea5430")) +   scale_y_continuous(expand = expansion(mult = c(.01, .01)), oob=squish) 
) 

myplot


#MeCP2
gene <- 'Mecp2'
data <- plotCounts(dds_PF4_list$Hip, gene=gene, intgroup=c("age_treatment"), returnData=T,normalized = T)

data$age_treatment <- factor(data$age_treatment,
                             levels=c('Young_Veh', 'Young_PF4',
                                      'Old_Veh', 'Old_PF4'), ordered = T)

data$age <- unlist(lapply(strsplit(as.character(row.names(data)), split = "_"), function(x) x[[1]]))
data$treatment <- unlist(lapply(strsplit(as.character(row.names(data)), split = "_"), function(x) x[[2]]))
data$age <- factor(data$age, levels = c('Young', 'Old'))
data$treatment <- factor(data$treatment, levels = c('Veh', 'PF4'))

(myplot <- ggplot(data, aes(x=age, y=count, color=treatment)) +
    geom_boxplot(outlier.colour = NA, notch = F) + xlab("") + ylab("Normalised counts")  +
    # geom_errorbar(aes(ymin=count-se, ymax=count+se),
    #               width=.2,                    # Width of the error bars
    #               position=position_dodge(.9)) +  
    theme_Publication_blank() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank()) + 
    geom_point(data=data, position = position_jitterdodge()) +
    scale_color_manual(values = c("Veh"="grey", "PF4"="#ea5430")) +   scale_y_continuous(expand = expansion(mult = c(.01, .01)), oob=squish) 
) 

myplot

#Arhgef9
gene <- 'Arhgef9'
data <- plotCounts(dds_PF4_list$Hip, gene=gene, intgroup=c("age_treatment"), returnData=T,normalized = T)

data$age_treatment <- factor(data$age_treatment,
                             levels=c('Young_Veh', 'Young_PF4',
                                      'Old_Veh', 'Old_PF4'), ordered = T)

data$age <- unlist(lapply(strsplit(as.character(row.names(data)), split = "_"), function(x) x[[1]]))
data$treatment <- unlist(lapply(strsplit(as.character(row.names(data)), split = "_"), function(x) x[[2]]))
data$age <- factor(data$age, levels = c('Young', 'Old'))
data$treatment <- factor(data$treatment, levels = c('Veh', 'PF4'))

(myplot <- ggplot(data, aes(x=age, y=count, color=treatment)) +
    geom_boxplot(outlier.colour = NA, notch = F) + xlab("") + ylab("Normalised counts")  +
    # geom_errorbar(aes(ymin=count-se, ymax=count+se),
    #               width=.2,                    # Width of the error bars
    #               position=position_dodge(.9)) +  
    theme_Publication_blank() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank()) + 
    geom_point(data=data, position = position_jitterdodge()) +
    scale_color_manual(values = c("Veh"="grey", "PF4"="#ea5430")) +   scale_y_continuous(expand = expansion(mult = c(.01, .01)), oob=squish) 
) 

myplot

