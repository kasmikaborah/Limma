# Limma
**Statistical significance analysis for RNA-seq data**

**load required libaries**
library(edgeR)
library(RColorBrewer)
library(gplots)
library(limma)
library(ggplot2)
library(ggrepel)

**load input data**
count1 <- read.csv("filte_low_geo.csv", header = TRUE, sep = ",", row.names = 1)
count1
dim(count1)#(38442   388)

**transpose your data**

dd1_t<-t(count1)
dd1_t
dim(dd1_t)

**Normalize the data to fit into ML model for downstream analysis**
min_max_normalization <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

normalized_data <- apply(dd1_t, 2, min_max_normalization)
head(normalized_data)
dim(normalized_data)
counts<-t(normalized_data)
dim(counts)
write.table(counts,file="tcga_normal.txt",quote=F,sep="\t",row.names=TRUE,col.name=T)

dim(normalized_data)
#View (normalized_data)
#abline(normalized_data=0.1)


#counts<-t(normalized_data)
#dim(counts)
#head(counts)
#class(dge)
 

#write.table(counts,file="normalized_data_geodata.txt",quote=F,sep="\t",row.names=TRUE,col.name=T)


d0 <- DGEList(count1)
d0

metadata <- read.csv("metadata.csv")
head(metadata)
dim(metadata)

identical(metadata$Run, colnames(count1))
counts <- df_clean[,metadata$Run]#if rows are not same than run this line)

**preprocessing**
d0 <- calcNormFactors(d0)
d0$samples

#plotMDS(d0, col = as.numeric(metadata))

#write.table(d0,"normalizationd0.txt",sep = "\t")
mm <- model.matrix(~0 + Stage, data = metadata)
mm

#logcpm transformation
logcpm<-cpm(d0, log=TRUE)
write.table(logcpm, file = 'logcpm_transformation_input200.txt', sep = "\t", quote = F)
#Voom transformation
####increasing the memory limit#############
memory.limit()###16180
memory.limit(size=800000)###41000
y <- voom(d0, mm, plot = T)
write.table(d0, file = 'voom_transformation_input200.txt', sep = "\t", quote = F)

#Fitting linear models in limma
fit <- lmFit(y, mm)
head(coef(fit))

#Specify which groups to compare using contrasts:
contr <- makeContrasts(StageNormal- StageTumor, levels = colnames(coef(fit)))
contr

#Estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)
tmp
tmp1 <- eBayes(tmp)
tmp1

#toptable with adjustment like BH...
top.table <- topTable(tmp1, adjust.method = "BH", sort.by = "P", n = Inf)
head(top.table,20)
length(which(top.table$adj.P.Val < 0.05))#(34961)
DEGlist<-top.table[which(top.table$adj.P.Val<0.05),]
head(DEGlist)

**save output data**
write.table(DEGlist,file="geo_sig.txt",quote=F,sep="\t",row.names=TRUE,col.name=T)

**Visualization of data**

# Define subsets with consistent thresholds
upregulated <- subset(DEGlist, logFC > 1.5)      # up-regulated
downregulated <- subset(DEGlist, logFC < -1.5)   # down-regulated (note: negative)

# Write to files
write.table(upregulated, file="upregulated.txt", quote=F, sep="\t", row.names=TRUE, col.names=T)
write.table(downregulated, file="downregulated.txt", quote=F, sep="\t", row.names=TRUE, col.names=T)

# Create volcano plot data frame
df <- data.frame(
  gene = rownames(DEGlist),
  logFC = DEGlist$logFC,
  adjPval = DEGlist$adj.P.Val,
  status = ifelse(DEGlist$logFC > 1.5, "Up", 
           ifelse(DEGlist$logFC < -1.5, "Down", "Unchanged"))
)

# Plot
ggplot(df, aes(x = logFC, y = -log10(adjPval))) +
  geom_point(aes(color = status)) +
  scale_color_manual(values = c("Up" = "red", "Down" = "green", "Unchanged" = "gray"))

  
**OR**

volcanoplot(DEGlist, highlight = 0L, hl.col = 'BLUE', coef = 1L, names = tmp1$genes$ID, xlab = 'Log2 Fold Change' , ylab = NULL, pch = 16, cex =0.35,values=c("U","D"),hl.col=c("green","red" ))
volcanoplot(DEGlist, highlight=8, names = rownames(tmp1), main="Tumor vs Normal")
upregulated<-subset(DEGlist, logFC > 1.5) # up-regulated
dim(upregulated)
write.table(downregulated,file="downregulated.txt",quote=F,sep="\t",row.names=TRUE,col.name=T)
downregulated<-subset(DEGlist, logFC < 1) # down-regulated
dim(downregulated)
write.table(upregulated,file="downregulated.txt",quote=F,sep="\t",row.names=TRUE,col.name=T)
df<-data.frame(upregulated,downpregulated)
ggplot(data = df, aes(y = y, x = x)) +
  geom_point(size=2, aes(colour=z)) +
  scale_color_manual("Status", colour= color)

plotWithHighlights(downregulated,upregulated,status=df,values=c("U","D"),hl.col=c("green","red"))
