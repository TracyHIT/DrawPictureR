library("DESeq2")
cts.count <- read.delim("~/Documents/DrawPictureR/data/dfr_0.01.count.txt")
condition <- factor(c("A.a","A.a","A.a","A.b","A.b","A.b",
                      "B.a","B.a","B.a","B.b","B.b","B.b",
                      "C.a","C.a","C.a","C.b","C.b","C.b",
                      "CK","CK","CK"))
colData <- data.frame(row.names=colnames(cts.count)[-1],condition)
database <- as.matrix(cts.count[,-1])
dds <- DESeqDataSetFromMatrix(countData = database, colData = colData, design= ~ condition)
dds <- DESeq(dds)
res <- results(dds)
table(res$padj < 0.05)
res <- res[order(res$padj),]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata,file = "test.csv")

#
test_brench