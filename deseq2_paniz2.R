#this is the workflow for tximport and DESeq2 to detect the effect of IGF2BP3 and IL-13 in BEAS-2B cells
R.version
#4.3.2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("refGenome")
#install.packages('refGenome')
install.packages('tximeta')
BiocManager::install("tximport")
# a and yes
BiocManager::install('GenomicFeatures')
# a and yes
BiocManager::install('DESeq2')
BiocManager::install('DelayedArray')
BiocManager::install("txdbmaker")
BiocManager::install("pheatmap")

library(tximport)
library(tidyverse)
library(GenomicFeatures)
library(DelayedArray) 
library(DESeq2)
library(magrittr)
library(txdbmaker)
library("pheatmap")

setwd(dir = '~/Documents/seq_data_paniz/input_salmon/')
getwd()

#need to make sure the files are in the same folder as where you're working for this next bit to work
sample_names <- list.files()
#want to convert data type to a matrix from 2D...use the function as. for it to happen
sample_names <- as.data.frame(sample_names)
sample_names <- sample_names[1:60,]
#put a comma to let it know its a dataframe(table) (not a list)
class(sample_names)
dim(sample_names)

sample_names

#lOAD METADTA CSV (that we made) 
samples <- read.csv('~/Documents/seq_data_paniz/metadata.csv')
dim(samples) 

#change donor_for_df to a factor rather than an interger
samples$donor <- as.factor(samples$donor)

# make three df, one per fraction

df_total <- filter(samples, fraction == "input")
df_free <- filter(samples, fraction == "free_mRNP")
df_polysome <- filter(samples, fraction == "polysome")

# Make txdb object from GFF file, make sure to use the same version as the one for alignment
txdb <- makeTxDbFromGFF("~/Documents/sp_perucha/gencode.v43.annotation.gtf.gz", 
                        format = "gtf")

saveDb(txdb, file="gencode.v43.sqlite")
# next time you can just load with this line (no need to makeTxDb...)
#txdb <- loadDb("~/Documents/seq_data_paniz/gencode.v43.sqlite")

# select the columns that we need now
columns(txdb)
k <- keys(txdb, "GENEID")
tx2gene <- select(txdb, keys = k, keytype = 'GENEID', columns = 'TXNAME')

#Important, TXNAME has to be the first column
res <- AnnotationDbi::select(txdb, k, "TXNAME", "GENEID")
head(res) # for transcript-level analysis
tx2gene <- res[,2:1] # for gene-level analysis
head(tx2gene)
write.csv(tx2gene, file='tx2gene.csv', row.names = FALSE, quote = FALSE )


#use tximport to read our count files. Since we are using the same gtf file, the versions of the transcripts will be the same, hence the argument ignoreTxVersion = FALSE. If not, we set it to TRUE
txi <- tximport(files = sample_names$sample_names, type = 'salmon', tx2gene = tx2gene, ignoreTxVersion = FALSE, ignoreAfterBar = TRUE)
summary(txi)
head(txi$counts)

# we set the column names to the sample names using our samples object
colnames(txi$counts) <- df_total$Sample_name
head(txi$counts)

#save the count file (abundance, counts, length)
write.csv(txi, "input_analysis/countmatrix.csv") # together

count_raw <- txi$counts # separately
write.csv(count_raw, "input_analysis/raw_counts_matrix.csv")

# we now import the txi object into DESeq (~ tells DESEQ what to compare)
#design considers the effect of condition, treatment, and then the interaction b/w condition with treatment
dds <- DESeqDataSetFromTximport(txi = txi, 
                                colData = df_total, 
                                design = ~ donor)
#The first thing we notice is that both our counts and average transcript length were used to construct the DESeq object. We also see a warning message, where our condition was converted to a factor. Both of these messages are ok to see!

vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, c("fraction", 'Sample_name'))
#saved image in folder plots, PCA_total_donors

### We make 3 analyses, one per fraction

total_files <- list.files()
#want to convert data type to a matrix from 2D...use the function as. for it to happen
total_files <- as.data.frame(total_files)
#total_files <- total_files[1:20,]
#put a comma to let it know its a dataframe(table) (not a list)
class(total_files)
dim(total_files)
head(total_files)

#use tximport to read our count files. Since we are using the same gtf file, the versions of the transcripts will be the same, hence the argument ignoreTxVersion = FALSE. If not, we set it to TRUE
txi <- tximport(files = total_files$total_files, type = 'salmon', tx2gene = tx2gene, ignoreTxVersion = FALSE, ignoreAfterBar = TRUE)
summary(txi)
head(txi$counts)

# we set the column names to the sample names using our df_total object
colnames(txi$counts) <- df_total$Sample_name
head(txi$counts)

# save the total counts raw
count_raw <- txi$counts # separately
write.csv(count_raw, "raw_counts_input_matrix.csv")


# we now import the txi object into DESeq (~ tells DESEQ what to compare)
#design considers the effect of condition, treatment, and then the interaction b/w condition with treatment
dds <- DESeqDataSetFromTximport(txi = txi, 
                                colData = df_total, 
                                design = ~ donor + genotype + condition + genotype:condition)

#to fix the ensembl annotations (incl. only the first 15 characters ignoring the versions)
head(rownames(dds))
table(duplicated(substr(rownames(dds),1,15)))
rownames(dds) <- make.unique(substr(rownames(dds),1,15))

dim(dds)

#Now we can run our differential expression pipeline. First, it is sometimes convenient to remove genes where all the samples have very small counts. It's less of an issue for the statistical methods, and mostly just wasted computation, as it is not possible for these genes to exhibit statistical significance for differential expression. Here we count how many genes (out of those with at least a single count) have 3 samples with a count of 10 or more:
dds <- dds[rowSums(counts(dds)) > 0,]
keep <- rowSums(counts(dds) >= 10) >= 5
table(keep)
# filter them out of dds matrix
dds[keep, ]
dds <- dds[keep, ]

m <- counts(dds)
write.csv(m, 'input_analysis/gene_input_norm_counts.csv')

# now we do differential gene expression

# as we have two different things to compare (IGF2BP3 effects and IL-13 effects), we make design for interactions
# we specify our references

dds$condition <- relevel(dds$condition, ref = "IL13")
dds$genotype <- relevel(dds$genotype, ref = "siC")

dds <- DESeq(dds)
dim(dds)
head(dds)
resultsNames(dds)

# the effect of IL13 for siC (the main effect) In 'contrast' the order is factor, numerator, denominator
res_controlIL13 <- results(dds, contrast=c("condition", "unstim", "IL13"))
res.sort <- res_controlIL13[order(res_controlIL13$padj),]
res.sort 
res_controlIL13$padj < 0.1

write.csv(res_controlIL13, '~/Documents/seq_data_paniz/input_salmon/input_analysis/siC_IL13_reordered.csv')

# the effect of IGF2BP3

res_IGF2BP3 <- results(dds, contrast = c('genotype', 'siIGF2BP3', 'siC'))
res.sort3 <- res_IGF2BP3[order(res_IGF2BP3$padj), ]
res.sort3$padj < 0.1

write.csv(res_IGF2BP3, '~/Documents/seq_data_paniz/input_salmon/input_analysis/siIGF2BP3_vs_siC.csv')

# the effect of IGF2BP3 KD on IL-13
res_IGF2BP3IL13 <- results(dds, contrast = list( c("condition_unstim_vs_IL13","genotypesiIGF2BP3.conditionunstim")))
res.sort2 <- res_IGF2BP3IL13[order(res_IGF2BP3IL13$padj),]
res.sort2

write.csv(res_IGF2BP3IL13, '~/Documents/seq_data_paniz/input_salmon/input_analysis/siIGF2BP3_IL-13.csv')

# the interaction term, answering: is the IL-13 effect *different* across genotypes?
res.IG.IL13 <- results(dds, name="genotypesiIGF2BP3.conditionunstim")
res.IG.IL13 <- res_IGF2BP3IL13[order(res_IGF2BP3IL13$padj),]
res.IG.IL13

write.csv(res.IG.IL13, '~/Documents/seq_data_paniz/input_salmon/input_analysis/IL-13_effects_IGF2BP3dependent.csv')

# heatmap of count matrix

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","genotype")])
ntd <- normTransform(dds)
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
