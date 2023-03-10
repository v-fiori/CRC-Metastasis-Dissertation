###################### METASTASIS POTENTIAL BIOMARKERS ##########################

library(rlang)
library(gplots)
library(ggplot2)
library(gtsummary)
library(edgeR)
library(DESeq2)
library(fgsea)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(pamr)
library(caret)

#import the data
counts <- read.delim("Kallisto_Gene_Counts_Early_Primaries_21_02_2022.tab")
meta <- read.delim("Kallisto_Meta_Early_Primaries_21_02_2022.tab")
genes <- rownames(counts)

#remove duplicated samples
colnames(counts) <- substr(colnames(counts),1,5)
counts <- counts[,!duplicated(colnames(counts))]
meta <- meta[!duplicated(substr(rownames(meta),1,5)),]
rownames(meta) <- meta$sample

#summarize the metadata
for (i in 1:nrow(meta)) {
  if (meta$PrimLocation[i] == "rectum" || meta$PrimLocation[i] == "colostomy" ||  meta$PrimLocation[i] == "transverse" || meta$PrimLocation[i] == "right/left" || is.na(meta$PrimLocation[i])) {
    meta$PrimLocation[i] = "Other"
  }}


sum_meta <- subset(meta, select = c("organ","stage","MetPot","CMS","PrimLocation","isdead","gender"))
sum_meta$MetPot <- ifelse(meta$MetPot == "Metastatic",1,0)
theme_gtsummary_journal(journal = "jama")
table_sum <- tbl_summary(sum_meta,by = MetPot) %>% add_overall() %>% bold_labels() %>% add_p()
table_sum2 <- tbl_summary(sum_meta,by = CMS) %>% add_overall() %>% bold_labels() %>% add_p()
library(flextable)
table_sum %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = "sum.docx")

#filtering and normalization with EdgeR
groups <- factor(meta$group,levels = c("PNM","PM"))
y <- DGEList(counts = counts, genes = genes, group = groups)
dim(y)

keep <- filterByExpr(y)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
genes_keep <- row.names(y)
counts_y <- counts[genes_keep,]

########## Principal Component Analysis
library(PCAtools)
y <- calcNormFactors(y)
p <- PCAtools::pca(cpm(y, normalized.lib.sizes = T, log=T), removeVar = 0.10, metadata = meta)
screeplot(p,axisLabSize = 18, titleLabSize = 22,components = getComponents(p,1:20))
biplot(p,
       #x = 'PC2',
       #y = 'PC4',
       showLoadings = F,
       colby = 'group',
       shape = 'cohort',
       #selectLab = rownames(meta)[c(60,61,68,69,75,76,86,87,91,92)],
       selectLab = c("CR379.1","CR438.1","CR512.1","CR566.1","CR585.1","CR379.2","CR438.2","CR512.2","CR566.2","CR585.2"),
       #colkey = c('CMS1' = 'forestgreen', 'CMS2' = 'purple', 'CMS3' = "red", 'CMS4' = 'darkgoldenrod3','Unknown' = "black"),
       hline = 0, vline = 0,
       encircle = FALSE,
       encircleFill = FALSE,
       encircleAlpha = 1, encircleLineSize = 3,
       legendPosition = 'right')



##### DESeq2 Differential Expression Analysis
dds <- DESeqDataSetFromMatrix(countData = round(counts_y) ,colData = meta, design = ~group)
dds$group <- relevel(dds$group, "PNM")
dds = estimateSizeFactors(dds)
ddsHTSeq <- DESeq(dds)

##### DEGs Results 
resultsNames(ddsHTSeq)
res <- results(ddsHTSeq)
res <- res[ order(res$padj),]
table(res$padj<0.05)

##### import genes' names
gene_info <- read.table("e100_biomart_gene_positions.tab", sep="\t", header = T)

###### filter DEGs for log2FC > 1 or < -1 and padj < 0.05
res_filter <- res[abs(res$log2FoldChange) > 1 & res$padj < 0.05, ]
##### Conver the gene codes into their respective names
res$gene_id <- row.names((res))
res <- merge(as.data.frame(res), gene_info, by.x="gene_id",by.y="Gene.stable.ID")
res <- res[ order(res$padj),]

res_filter$gene_id <- row.names((res_filter))
res_filter <- merge(as.data.frame(res_filter), gene_info, by.x="gene_id",by.y="Gene.stable.ID")
res_filter <- res_filter[ order(res_filter$padj),]

# check if there are names missing
sum(is.na(res_filter$Gene.name))
sum(is.na(res$Gene.name))


######## Heatmap dos top DEGs
#normalize the counts from DESeq2
normCounts <- counts(ddsHTSeq, normalized = TRUE)
#res_stage_filt <- res_stage[abs(res_stage$log2FoldChange) > 1 & res_stage$padj < 0.05, ]
#selecionar os 20 genes
select <- res_filter$gene_id[1:21]
log10_normCounts <- log10(normCounts + 1)
values <- log10_normCounts[select,]
values_2 <- as.data.frame(values)
values_2$id <- row.names((values_2))
values_2 <- merge(values_2, gene_info, by.x="id",by.y="Gene.stable.ID")
row.names(values_2) <- values_2$Gene.name
values_2$id <- NULL
values <- values_2[,1:ncol(values)]

#reorder the samples to appear first the ones with metastasis
meta2 <- meta[order(meta$MetPot),]

values <- values[,match(row.names(meta2),colnames(values))]
values <- values[-21,]


metpot_colors <- unlist(lapply(meta2$MetPot,function(x){
  if(grepl("NonMetastatic",x)) '#53FA0B'
  else if(grepl("Metastatic",x)) '#000000'
}))

heatmap.2(as.matrix(values),
          labRow = row.names(values),
          col = bluered(20),
          margins = c(8,9),
          trace = "none",
          density.info = "none",
          labCol = FALSE,
          scale = "row",
          ColSideColors = metpot_colors,
          dendrogram='none',
          #main = "Heatmap of the 20 top DEGs expression levels",
          Colv = FALSE)

#distfun = function(x) as.dist(1 - cor(t(x))))
legend("left", title = "Metastatic Potential",legend=c("Metastatic","NonMetastatic"), 
       fill=c("#000000","#53FA0B"), cex=0.8, box.lty=1)

####### box plot of the 20 top DEGs ##########

library(ggplot2)

met_box <- rep(c("Met"),31) 
nmet_box <- rep(c("NonMet"),118)
metpot_box <- c(met_box,nmet_box)
t_values <- as.data.frame(t(values))
data_box <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(data_box) <- c("Gene","MetPot","Gene Expression")
for (i in colnames(t_values)) {
  genes_box <- rep(c(i),149)
  express_box <- c(t_values[,i])
  data_box2 <- data.frame(genes_box,metpot_box,express_box)
  colnames(data_box2) <- c("Gene","MetPot","Gene Expression")
  data_box <- rbind(data_box,data_box2)
  
}

ggplot(data_box, aes(x = Gene, y = `Gene Expression`,fill = MetPot)) + 
  geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  #scale_fill_manual(breaks = data_box$MetPot, values=c("red", "green")) +
  scale_fill_manual(values=c("#69b3a2", "grey")) +
  labs(fill = "Metastatic Potential")


####### Volcano Plot

keyvals <- ifelse(
  res$log2FoldChange < -1 & res$padj < 0.05, 'blue',
  ifelse(res$log2FoldChange > 1 & res$padj < 0.05, 'red',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red'] <- 'Up Regulated'
names(keyvals)[keyvals == 'black'] <- 'Not DEG'
names(keyvals)[keyvals == 'blue'] <- 'Down Regulated'

EnhancedVolcano(res,
                lab = res$Gene.name,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Metastatic vs NonMetastatic',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 2.0,
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5)



#################Gene Set Enrichment Analysis
##REACTOME
res$gene_id <- row.names((res))
gene_info <- read.table("e100_biomart_gene_positions.tab", sep="\t", header = T)
res <- merge(as.data.frame(res), gene_info, by.x="gene_id",by.y="Gene.stable.ID")
res <- res[ order(res$padj),]

ranks <- res$log2FoldChange
names(ranks) <- res$Gene.name

library(fgsea)
reactome <-  gmtPathways("c2.cp.reactome.v7.1.symbols.gmt")
set.seed(101)
fgseaRes <- fgsea(pathways = reactome, 
                  stats    = ranks,
                  minSize  = 15,
                  maxSize  = 500)
collapsedPathways <- collapsePathways(fgseaRes[order(padj)][padj < 0.05], reactome, ranks)
mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]

down_pathways <- head(fgseaRes[(fgseaRes$NES < 0) & (fgseaRes$pathway %in% mainPathways),], n=10)
down_pathways <- down_pathways[order(down_pathways$NES, decreasing = T),]

up_pathways <- head(fgseaRes[(fgseaRes$NES > 0) & (fgseaRes$pathway %in% mainPathways),], n=10)
up_pathways <- up_pathways[order(up_pathways$NES, decreasing = T),]

select_pathways <- rbind(up_pathways,down_pathways)
select_pathways$pathway <- sub("REACTOME_","",select_pathways$pathway)
select_pathways$pathway <- gsub("_"," ",select_pathways$pathway)


nup <- nrow(up_pathways)
ndown <- nrow(down_pathways)

bp <- barplot(select_pathways$NES, names.arg = "", xlim=c(-4,4),
              horiz = T, col = c(rep("darkgreen",nup),rep("darkred",ndown)),main ="Gene Set Enrichment Analysis Met vs NonMet", xlab="Normalized Enrichment Score (NES)")
text(x=-0.1, y=bp[1:nup],select_pathways$pathway[1:nup], cex=0.55, adj=1)
text(x=0.1, y=bp[(nup+1):(nup+1+ndown)],select_pathways$pathway[(nup+1):(nup+1+ndown)], cex=0.55, adj=0)

#Get the DEGs that are on the top enriched pathways
enriched_pathways <- cbind(up_pathways$pathway,down_pathways$pathway)
enr_pathways_info <- reactome[enriched_pathways]
degs <- res_filter$Gene.name
gsea_degs <- c()
gsea_pathways <- c()
for (i in degs){
  for (j in names(enr_pathways_info)){
    if (i %in% unlist(enr_pathways_info[j])){
      gsea_degs <- cbind(gsea_degs,i)
      gsea_pathways <- cbind(gsea_pathways,j)
    }
  }
}
names(gsea_degs) <- gsea_pathways

########### CANCER HALLMARKS
paths <-  gmtPathways("h.all.v7.4.symbols.gmt")
set.seed(101)
fgseaRes <- fgsea(pathways = paths, 
                  stats    = ranks,
                  minSize  = 15,
                  maxSize  = 500)
collapsedPathways <- collapsePathways(fgseaRes[order(padj)][padj < 0.05], paths, ranks)
mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]

down_pathways <- head(fgseaRes[(fgseaRes$NES < 0) & (fgseaRes$pathway %in% mainPathways),], n=10)
down_pathways <- down_pathways[order(down_pathways$NES, decreasing = T),]

up_pathways <- head(fgseaRes[(fgseaRes$NES > 0) & (fgseaRes$pathway %in% mainPathways),], n=10)
up_pathways <- up_pathways[order(up_pathways$NES, decreasing = T),]

select_pathways <- rbind(up_pathways,down_pathways)
select_pathways$pathway <- sub("REACTOME_","",select_pathways$pathway)
select_pathways$pathway <- gsub("_"," ",select_pathways$pathway)


nup <- nrow(up_pathways)
ndown <- nrow(down_pathways)

bp <- barplot(select_pathways$NES, names.arg = "", xlim=c(-4,4),
              horiz = T, col = c(rep("darkgreen",nup),rep("darkred",ndown)),main = "Gene Set Enrichment Analysis Met vs NonMet", xlab="Normalized Enrichment Score (NES)")
text(x=-0.1, y=bp[1:nup],select_pathways$pathway[1:nup], cex=0.55, adj=1)
text(x=0.1, y=bp[(nup+1):(nup+1+ndown)],select_pathways$pathway[(nup+1):(nup+1+ndown)], cex=0.55, adj=0)


#Get the DEGs that are on the top enriched pathways
enriched_pathways_2 <- cbind(up_pathways$pathway,down_pathways$pathway)
enr_pathways_info_2 <- paths[enriched_pathways_2]
degs_2 <- res_filter$Gene.name
gsea_degs_2 <- c()
gsea_pathways_2 <- c()
for (i in degs_2){
  for (j in names(enr_pathways_info_2)){
    if (i %in% unlist(enr_pathways_info_2[j])){
      gsea_degs_2 <- cbind(gsea_degs_2,i)
      gsea_pathways_2 <- cbind(gsea_pathways_2,j)
    }
  }
}
names(gsea_degs_2) <- gsea_pathways_2



#################     Predictive Model with PAMR
logCPMs <- cpm(y,normalized.lib.sizes = T, log = T)
train_data <- logCPMs
train_data_y <- factor(meta$group)
set.seed(101)
#train with 70% of data
trainIndex <- createDataPartition(train_data_y, p = 0.7, list = FALSE, times = 1)
pamr.data.train <- list(x=train_data[,trainIndex],y=train_data_y[trainIndex],
                        geneid=rownames(train_data),
                        genenames=rownames(train_data))
pamr.train <- pamr.train(pamr.data.train)
pamr.results<- pamr.cv(fit = pamr.train, data = pamr.data.train)
pamr.plotcv(pamr.results)
pamr.results
pamr.confusion(pamr.results, threshold = 1.8)

#test on the other 30%
pamr.data.test <- list(x=train_data[,-trainIndex],y=train_data_y[-trainIndex],
                       geneid=rownames(train_data),
                       genenames=rownames(train_data))
pamr.results.test <- pamr.predict(pamr.train, pamr.data.test$x,
                                  threshold = 1.8,
                                  prior = pamr.train$prior,
                                  threshold.scale = pamr.train$threshold.scale)
table("prediction" = pamr.results.test, "reality" = meta$group[-trainIndex])

# select the genes used that are also DEGs 
pamr.results.train.genes <- pamr.listgenes(pamr.train, pamr.data.train, threshold=1.8)
pamr.results.train.genes <- as.data.frame(pamr.results.train.genes)
head(pamr.results.train.genes$id)
CPMs_pam <- counts[rownames(counts) %in% pamr.results.train.genes$id,]

CPMs_pam$id <- rownames(CPMs_pam)
gene_info <- read.table("e100_biomart_gene_positions.tab", sep="\t", header = T)
CPMs_pam <- merge(as.data.frame(CPMs_pam), gene_info, by.x="id",by.y="Gene.stable.ID")
#Heatmap for the intersection with the DEGs
#reorder the samples to appear first the ones with metastasis
legend <- meta[order(meta$MetPot),]
values <- CPMs_pam[,match(row.names(legend),colnames(CPMs_pam))]
values$id <- CPMs_pam$id
values <- merge(values, gene_info, by.x="id",by.y="Gene.stable.ID")

metpot_colors <- unlist(lapply(legend$MetPot,function(x){
  if(grepl("NonMetastatic",x)) '#53FA0B'
  else if(grepl("Metastatic",x)) '#000000'
}))
heatmap.2(as.matrix(values[2:150]),
          labRow = values$Gene.name,
          col = bluered(20),
          margins = c(8,9),
          trace = "none",
          density.info = "none",
          scale = "row",
          labCol = FALSE,
          ColSideColors = metpot_colors,
          dendrogram='none',
          main = "Heatmap of the PAM classificator genes",
          Colv = FALSE)


legend("top", title = "Metastatic Potential",legend=c("Metastatic","NonMetastatic"), 
       fill=c("#000000","#53FA0B"), cex=0.8, box.lty=1)

######SNP impact
genes_m <- table(mutations_filt$Gene.refGene[mutations_filt$sample %in% metast])
genes_nm <- table(mutations_filt$Gene.refGene[mutations_filt$sample %in% n_metast])
genes_unfilt_m <- table(mutations_unfilt$Gene.refGene[mutations_unfilt$sample %in% metast])
genes_unfilt_nm <- table(mutations_unfilt$Gene.refGene[mutations_unfilt$sample %in% n_metast])
genes_m <- as.data.frame(genes_m)
genes_nm <- as.data.frame(genes_nm)
genes_unfilt_m <- as.data.frame(genes_unfilt_m)
genes_unfilt_nm <- as.data.frame(genes_unfilt_nm)
total_genes <- merge.data.frame(genes_m,genes_unfilt_m, by.x = "Var1", by.y = "Var1", all.x = TRUE)
total_genes$Freq.y <- total_genes$Freq.y - total_genes$Freq.x
total_genes <- merge.data.frame(total_genes,genes_nm, by.x = "Var1", by.y = "Var1", all.x = TRUE)
colnames(total_genes) <- c("Genes", "M_impact","M_no_impact","NM_impact")
total_genes <- merge.data.frame(total_genes,genes_unfilt_nm, by.x = "Genes", by.y = "Var1", all.x = TRUE)
colnames(total_genes) <- c("Genes", "M_impact","M_no_impact","NM_impact","NM_no_impact")
total_genes$NM_no_impact <- total_genes$NM_no_impact - total_genes$NM_impact
total_genes <- total_genes[order(total_genes$M_impact, decreasing = TRUE),]
total_genes <- na.omit(total_genes)


######### fisher's test

library(stats)
total_genes$p.value <- NA
total_genes$odds_ratio <- NA
for (i in 1:nrow(total_genes)) {
  df <-  data.frame("Met" = c(total_genes$M_impact[i],total_genes$M_no_impact[i]),
                    "Non_Met" = c(total_genes$NM_impact[i],total_genes$NM_no_impact[i]),
                    row.names = c("Impact","No_Impact"))
  fisher <- fisher.test(df)
  total_genes$p.value[i] <- fisher$p.value
  total_genes$odds_ratio[i] <- fisher$estimate
}
filt_total <- total_genes#subset(total_genes, rowSums(total_genes[,2:5])>5)
filt_total <- filt_total[order(filt_total$p.value),]
filt_total$adj.p.value <- p.adjust(filt_total$p.value, method = "fdr", n = length(filt_total$p.value))

#################################CNVS
#############genes
sam1 <- unique(cnvs_cohort1$sample)
sam2 <- unique(cnvs_cohort2$sample)
same <- intersect(sam1,sam2)

unique_cnvs1 <- unique(cnvs_cohort1[c("sample","geneName","CNV")])
unique_cnvs2 <- unique(cnvs_cohort2[c("sample","geneName","CNV")])

genes_cnv1_met <- as.data.frame(table(unique_cnvs1$geneName[unique_cnvs1$sample %in% metast]))
genes_cnv2_met <- as.data.frame(table(unique_cnvs2$geneName[unique_cnvs2$sample %in% metast]))
genes_cnv1_nonmet <- as.data.frame(table(unique_cnvs1$geneName[unique_cnvs1$sample %in% n_metast]))
genes_cnv2_nonmet <- as.data.frame(table(unique_cnvs2$geneName[unique_cnvs2$sample %in% n_metast]))
genes_cnv1_met <- genes_cnv1_met[order(genes_cnv1_met$Freq, decreasing = TRUE),]
genes_cnv2_met <- genes_cnv2_met[order(genes_cnv2_met$Freq,decreasing = TRUE),]
genes_cnv1_nonmet <- genes_cnv1_nonmet[order(genes_cnv1_nonmet$Freq,decreasing = TRUE),]
genes_cnv2_nonmet <- genes_cnv2_nonmet[order(genes_cnv2_nonmet$Freq,decreasing = TRUE),]
samples <- unique(unique_cnvs1$sample)
genes_cnv1_met$No_cnv <- sum(samples %in% metast, na.rm = TRUE)
genes_cnv1_met$No_cnv <- genes_cnv1_met$No_cnv - genes_cnv1_met$Freq
genes_cnv1_nonmet$No_cnv <- sum(samples %in% n_metast, na.rm = TRUE)
genes_cnv1_nonmet$No_cnv <- genes_cnv1_nonmet$No_cnv - genes_cnv1_nonmet$Freq
samples_2 <- unique(unique_cnvs2$sample)
genes_cnv2_met$No_cnv <- sum(samples_2 %in% metast, na.rm = TRUE)
genes_cnv2_met$No_cnv <- genes_cnv2_met$No_cnv - genes_cnv2_met$Freq
genes_cnv2_nonmet$No_cnv <- sum(samples_2 %in% n_metast, na.rm = TRUE)
genes_cnv2_nonmet$No_cnv <- genes_cnv2_nonmet$No_cnv - genes_cnv2_nonmet$Freq

#fisher for cohort1
total_cnvs1 <- merge.data.frame(genes_cnv1_met,genes_cnv1_nonmet, by.x = "Var1", by.y = "Var1", all.x = TRUE)
colnames(total_cnvs1) <- c("Genes", "M_cnv","M_no_cnv","NM_cnv","NM_no_cnv")
total_cnvs1$NM_cnv[is.na(total_cnvs1$NM_cnv)] =  0
total_cnvs1$NM_no_cnv[is.na(total_cnvs1$NM_no_cnv)] = sum(samples %in% n_metast, na.rm = TRUE)


library(stats)
total_cnvs1$p.value <- NA
total_cnvs1$odds_ratio <- NA
for (i in 1:nrow(total_cnvs1)) {
  df <-  data.frame("Met" = c(total_cnvs1$M_cnv[i],total_cnvs1$M_no_cnv[i]),
                    "Non_Met" = c(total_cnvs1$NM_cnv[i],total_cnvs1$NM_no_cnv[i]),
                    row.names = c("Cnv","No_cnv"))
  fisher <- fisher.test(df)
  total_cnvs1$p.value[i] <- fisher$p.value
  total_cnvs1$odds_ratio[i] <- fisher$estimate
}

total_cnvs1$adj.p.value <- p.adjust(total_cnvs1$p.value, method = "fdr", n = length(total_cnvs1$p.value))
total_cnvs1 <- total_cnvs1[order(total_cnvs1$p.value),]



#cohort2
total_cnvs2 <- merge.data.frame(genes_cnv2_met,genes_cnv2_nonmet, by.x = "Var1", by.y = "Var1", all.x = TRUE)
colnames(total_cnvs2) <- c("Genes", "M_cnv","M_no_cnv","NM_cnv","NM_no_cnv")
total_cnvs2$NM_cnv[is.na(total_cnvs2$NM_cnv)] =  0
total_cnvs2$NM_no_cnv[is.na(total_cnvs2$NM_no_cnv)] = sum(samples %in% n_metast, na.rm = TRUE)


library(stats)
total_cnvs2$p.value <- NA
total_cnvs2$odds_ratio <- NA
for (i in 1:nrow(total_cnvs2)) {
  df <-  data.frame("Met" = c(total_cnvs2$M_cnv[i],total_cnvs2$M_no_cnv[i]),
                    "Non_Met" = c(total_cnvs2$NM_cnv[i],total_cnvs2$NM_no_cnv[i]),
                    row.names = c("Cnv","No_cnv"))
  fisher <- fisher.test(df)
  total_cnvs2$p.value[i] <- fisher$p.value
  total_cnvs2$odds_ratio[i] <- fisher$estimate
}
total_cnvs2 <- total_cnvs2[order(total_cnvs2$p.value),]
total_cnvs2$adj.p.value <- p.adjust(total_cnvs2$p.value, method = "fdr", n = length(total_cnvs2$p.value))
total_cnvs2 <- total_cnvs2[order(total_cnvs2$adj.p.value),]
