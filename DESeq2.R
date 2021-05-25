library(DESeq2)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)

# read.table
X_data <- read.table("gene-count-matrix.txt", header = T, row.names = 1)
names(X_data)<- c("T01", "T02", "T03", "T04", "T05", "T06", "T07", "T08", "T09", "T10")


# 控制条件：因子
input_data <- X_data[c("T01", "T02", "T03", "T04", "T05", "T06")]
condition <- factor(c(rep("WT", 3), rep("KO", 3)))

# 取整,函数round
input_data <- input_data[which(rowSums(input_data) > 0), ]
input_data <- round(input_data, digits = 0)

# 准备文件
# as.matrix 将输入文件转换为表达矩阵；
input_data <- as.matrix(input_data)

# input_data根据控制条件构建data.frame
coldata <- data.frame(row.names = colnames(input_data), condition)
head(coldata)

# build deseq input matrix 构建输入矩阵
# countData作为矩阵的input_data；colData Data.frame格式；控制条件design;
dds <- DESeqDataSetFromMatrix(countData = input_data, colData = coldata, design = ~condition)
dds$condition <- relevel(dds$condition, ref = "WT")

# DESeq2 进行差异分析
dds <- DESeq(dds)
# 实际包含3步

# 提取结果
# dds dataset格式转换为DESeq2 中result数据格式，矫正值默认0.1
res <- results(dds, alpha = 0.05)
# 查看res(DESeqresults格式),可以看到上下调基因
summary(res)

# res(resultset)按照P值排序
res <- res[order(res$padj), ]
res

table(res$padj < 0.05)

# 将进过矫正后的表达量结果加进去；
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized = T)),
                 by = "row.names", sort = F)
names(resdata)[1] <- "ID_Gene"

resdata <- separate(resdata, ID_Gene, into = c("ID", "Gene"), sep = "_")
# 查看(resdata)
head(resdata)

######
# 画火山图

# 按照padj分类
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

# 按照pvalue分类
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

# output result 输出结果
write.table(resdata, file = "X.KO vs WT.diffexpr-results.txt", 
            sep = "\t", quote = F, row.names = F)
