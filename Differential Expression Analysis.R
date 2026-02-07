library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(org.EcK12.eg.db)

C <- read.csv(":/Cfeaturecount_R2-counts.csv")
D <- read.csv(":/Dfeaturecount_R2-counts.csv")
E <- read.csv(":/Efeaturecount_R2-counts.csv")
F <- read.csv(":/Ffeaturecount_R2-counts.csv")


CD_EF = data.frame(C = C[,7], D = D[,7],E = E[,7], F = F[,7])
rownames(CD_EF) <- C[,1]
CD_EF_2 = CD_EF[rowSums(CD_EF) > 2,]


setcolor <- colorRampPalette(brewer.pal(9,"GnBu"))(100)
cordata <- cor(as.matrix(CD_EF_2 ), method="pearson")
clusdata <- hcluster(t(CD_EF_2 ), method="pearson")
heatmap.2(cordata, Rowv = as.dendrogram(clusdata), symm=T, trace="none",col=setcolor,
          margins=c(11,11), main="Pearson correlation")


conditions <- c('control','control','exp','exp')
batch <- factor(c(1,1,1,1))
sample <- data.frame(conditions,batch)

ddsFullCountTable<- DESeqDataSetFromMatrix(countData = CD_EF_2,
                                           colData = sample,  design = ~ + conditions)
dds <- DESeq(ddsFullCountTable)
normalized_counts <- counts(dds, normalized=TRUE)

sampleB = 'control'
sampleA = 'exp'
contrastV <- c("conditions",sampleA,sampleB)#sampleA/sampleB
res <- results(dds,contrast = contrastV)
res <- cbind(CD_EF_2 , normalized_counts, ID = rownames(res), as.data.frame(res))
res <- cbind.data.frame(counts(dds, normalized = TRUE),res)

res$ID<-CD_EF_2[,1]
res_de <- subset(res, res$pva<0.05)
res_de_up <- subset(res_de, res_de$log2FoldChange >= 1)
res_de_down <- subset(res_de, res_de$log2FoldChange<=(-1)*1)

res_de_up_id = data.frame(ID=res_de_up$ID)
res_de_down_id = data.frame(ID=res_de_down$ID)
de_id = rbind(res_de_up_id, res_de_down_id)
FDR <- res$padj


logCounts <- log2( res$baseMean+1)
logFC <- res$log2FoldChange
pvalue <- res$pvalue
cut_off_qvalue = 0.05
cut_off_logFC = 2
res$Sig <- ifelse(res$pvalue < cut_off_qvalue & 
                    abs(res$log2FoldChange) >= cut_off_logFC,
                  ifelse(res$log2FoldChange > cut_off_logFC ,'UP','DOWN'),'NS')
res <- data.frame(res)


ggplot(res, aes(x = logFC, y = -1*log10(pvalue), colour=Sig)) +
  geom_point(alpha=0.4, size=2) +
  scale_color_manual(values=c("#546de5","gray", "#ff4757"))+
  geom_vline(xintercept=c(-2,2),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_qvalue),
             lty=4,col="black",lwd=0.8) +
  labs(x="log2(Fold Change)",
       y="-log10 (pvalue)")+
  theme_bw()+
  ggtitle("Volcano Plot")+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank(),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank()   
  )



up_name <- row.names(res_de_up)
keggup <- enrichKEGG(gene = up_name , organism ='eco' , 
                     pAdjustMethod = 'fdr',
                     pvalueCutoff = 0.1 , 
                     qvalueCutoff = 0.1,) 
barplot(keggup,)

down_name <- row.names(res_de_down)
keggdown <- enrichKEGG(gene = down_name , organism ='eco' , 
                       pAdjustMethod = 'fdr',
                       pvalueCutoff = 0.5 , 
                       qvalueCutoff = 0.5,) 
barplot(keggdown,)

