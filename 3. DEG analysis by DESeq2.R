library(DESeq2)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)

# read.table
X_data <- read.table("gene-count-matrix.txt", header = T, row.names = 1)
names(X_data)<- c("T01", "T02", "T03", "T04", "T05", "T06", "T07", "T08", "T09", "T10")

# 
input_data <- X_data[c("T01", "T02", "T03", "T04", "T05", "T06")]
condition <- factor(c(rep("WT", 3), rep("KO", 3)))

# 
input_data <- input_data[which(rowSums(input_data) > 0), ]
input_data <- round(input_data, digits = 0)
input_data <- as.matrix(input_data)

# 
coldata <- data.frame(row.names = colnames(input_data), condition)
head(coldata)

# build deseq input matrix
dds <- DESeqDataSetFromMatrix(countData = input_data, colData = coldata, design = ~condition)
dds$condition <- relevel(dds$condition, ref = "WT")

# DESeq2 
dds <- DESeq(dds)

# 
res <- results(dds, alpha = 0.05)
summary(res)
res <- res[order(res$padj), ]
res
table(res$padj < 0.05)

# 
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized = T)),
                 by = "row.names", sort = F)
names(resdata)[1] <- "ID_Gene"
resdata <- separate(resdata, ID_Gene, into = c("ID", "Gene"), sep = "_")
head(resdata)

######
# Volcano plot

# padj
resdata$significant <- ifelse(resdata$padj < 0.05 & resdata$log2FoldChange > 1, "UP",
                              ifelse(resdata$padj < 0.05 & resdata$log2FoldChange < -1,
                                     "DOWN", "FALSE"))

ggplot(data = resdata, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point() + xlim(-10, 10) + #ylim(0, 6) +
  scale_color_manual(values = c("blue", "grey40", "red")) +
  geom_hline(yintercept = -log10(0.05), lty = 4, lwd = 0.6, alpha = 0.8, colour = "grey") +
  geom_vline(xintercept = c(1,-1), lty = 4, lwd = 0.6, alpha = 0.8, colour = "grey") +
  theme_bw() + guides(color = F) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(title = "KO vs WT", x = "log2(fold change)", y = "-log10(adjusted P)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text_repel(data = subset(resdata, -log10(padj) > 5 & abs(log2FoldChange) > 1), 
                  aes(label = Gene), col = "black", alpha = 0.8)

# pvalue
resdata$significant <- ifelse(resdata$pvalue < 0.05 & resdata$log2FoldChange > 1, "UP",
                              ifelse(resdata$pvalue < 0.05 & resdata$log2FoldChange < -1,
                                     "DOWN", "FALSE"))

ggplot(data = resdata, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point() + #xlim(-10, 10) + #ylim(0, 6) +
  scale_color_manual(values = c("blue", "grey40", "red")) +
  geom_hline(yintercept = -log10(0.05), lty = 4, lwd = 0.6, alpha = 0.8, colour = "grey") +
  geom_vline(xintercept = c(1,-1), lty = 4, lwd = 0.6, alpha = 0.8, colour = "grey") +
  theme_bw() + guides(color = F) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(title = "KO vs WT", x = "log2(fold change)", y = "-log10(P-value)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text_repel(data = subset(resdata, -log10(pvalue) > 10 & abs(log2FoldChange) > 1), 
                  aes(label = Gene), col = "black", alpha = 0.8)

# output results
write.table(resdata, file = "X.KO vs WT.diffexpr-results.txt", 
            sep = "\t", quote = F, row.names = F)

