zxcvbnk
code for edgeR
------------------

library(edgeR) # For an alternative analysis approach
library(AnnotationDbi)


countData <- as.matrix(read.csv("count_GSE135_2.csv", row.names = 1))
# Create a DataFrame with sample information
sampleInfo <- read.csv("GSE58135.csv")


# Define the "group" vector
group <- factor(c("1", "1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2"))



# Create a DGEList object
dgeGlm <- DGEList(counts = countData, group = group)
str(dgeGlm)
str(group)
group




# Filter lowly expressed genes
keep <- rowSums(cpm(dgeGlm)>2) >= 10 #count per million (CPM)
dgeGlm <- dgeGlm[keep,] 

#Perform TMM normalization
dgeGlm <- calcNormFactors(dgeGlm)

names(dgeGlm)
dgeGlm[["samples"]]


design <- model.matrix(~group)
design


# fitting a generalized linear model (GLM)
dgeGlmComDisp <- estimateGLMCommonDisp(dgeGlm,design, verbose = TRUE)

dgeGlmTrendDisp <- estimateGLMTrendedDisp(dgeGlmComDisp, design)

dgeGlmTagDisp <- estimateGLMTagwiseDisp(dgeGlmTrendDisp,design)


plotBCV(dgeGlmTagDisp)

fit<- glmFit(dgeGlmTagDisp,design)

colnames(coef(fit))

lrt <- glmLRT(fit, coef = 2) #performing likelihood ratio tests (LRT)

ttGlm <- topTags(lrt, n=Inf)
head(ttGlm)
class(ttGlm)
summary(decideTests(lrt))

summary(deGlm <- decideTestsDGE(lrt, p=0.01, adjust= "fdr"))
tagsGlm <- rownames(dgeGlmTagDisp)[as.logical(deGlm)]


plotSmear(lrt, de.tags = tagsGlm)
hits2 <- ttGlm$table[ttGlm$table$FDR < 0.01, ]
head(hits2)
####### Sorting
df <- as.data.frame(hits2, logFC  = A1:A11410, FDR = E1:E11410)

subdif <- subset(df, logFC >1 & FDR < 0.01 )
subdif1 <- subset(df, logFC < -1 & FDR < 0.01 )
write.csv(df,"diff_edgeR_GSE58135.csv")
write.csv(subdif, "edgeR_up.csv")
write.csv(subdif1, "edgeR_down.csv")
