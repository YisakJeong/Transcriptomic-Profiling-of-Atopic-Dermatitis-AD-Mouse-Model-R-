BiocManager::install("tximport")

library(tximport)
library(biomaRt)
library(DESeq2)
library(dplyr)


##### data load
base_dir <- "C:/Users/home/Documents/choilab 인턴/7월/프로젝트/20260724_salmon_JE36_output_ss/salmon_JE36"
files <- list.files(path=base_dir,pattern="quant.sf",recursive=TRUE,full.names = TRUE)
names(files) <- basename(dirname(files)) # 각 파일의 이름으로 샘플명 자동 부여
sample_names <- names(files) # 샘플 이름 목록


# group 이름 지정
group <- case_when(
  grepl("Cont", sample_names)     ~ "Control",
  grepl("Nor", sample_names)      ~ "Normal",
  grepl("JE36-300", sample_names) ~ "JE36-300", # 긴 문자열을 위로!
  grepl("JE36-100", sample_names) ~ "JE36-100",
  grepl("JE36-30", sample_names)  ~ "JE36-30",
  grepl("Posit", sample_names)    ~ "Positive",
  TRUE               ~ "Postl"
)

table(sample_names, group)



samples <- data.frame(
  
  group=factor(group),
  row.names=sample_names
)


#BiocManager::install("EnsDb.Mmusculus.v79")
library(EnsDb.Mmusculus.v79)
library(ensembldb)



tx2gene <- transcripts(
  EnsDb.Mmusculus.v79,
  columns = c("tx_id", "gene_name")
) |>
  as.data.frame() |>
  dplyr::select(tx_id, gene_name)


head(tx2gene)
txi <- tximport(files,type="salmon",tx2gene=tx2gene,ignoreTxVersion = TRUE) # txi에 abundance(TPM), counts(estimated), Length(유효길이)
dds <- DESeqDataSetFromTximport(txi,colData=samples,design=~group)
keep <- rowSums(counts(dds))>=10 # 모든 샘플의 readcount 합이 10 이상인 애들만 남기기
dds <- dds[keep,]


# 3. 두 그룹으로만 새로 vst 정규화 진행
vsd <- vst(dds, blind = TRUE)
colnames(vsd) <- gsub("^quant_", "", colnames(vsd))


# 히트맵 위해서 
# 1. Control과 Normal 그룹에 해당하는 샘플만 필터링  
keep_samples <- dds$group %in% c("Normal", "Control")
dds_sub <- dds[, keep_samples]

# 2. 안 쓰는 그룹 라벨(JE36, Positive 등) 제거 (Factor Refactor)
dds_sub$group <- factor(dds_sub$group, levels = c("Normal", "Control"))

vsd_sub <-vst(dds_sub,blind=TRUE)
colnames(vsd_sub) <- gsub("^quant_", "", colnames(vsd_sub))


library(ggplot2)
library(ggrepel)

p_raw <- plotPCA(vsd,intgroup="group")+
  geom_text_repel(aes(label=name),size=3,max.overlaps=20)+
  ggtitle("Raw Data PCA (Before Batch Corrrection)")+
  theme_minimal()
print(p_raw)


# PC plot for only AD vs Normal 
p_ad_norm <- plotPCA(vsd_sub,intgroup="group")+
  geom_text_repel(aes(label=name),size=3,max.overlaps=20)+
  ggtitle("PCA plot for AD vs Normal")+
  theme_minimal()
print(p_ad_norm)


# Sample Clustering
corMatrix <- cor(assay(vsd_sub))
library(RColorBrewer)
# 2. 색상 테마 설정
my_colors <- colorRampPalette(c("navy", "white", "red"))(255)

annot_info <- data.frame(
  group = vsd_sub$group,
  row.names = colnames(vsd_sub) # vsd_blind에 적용된 진짜 샘플명
)
library(pheatmap)
pheatmap(corMatrix,
         color = my_colors,
         annotation_col = annot_info, # 상단에 group 표시
         )


## batch correction
library(sva)
library(limma)

library(tibble)
library(tidyr)
target_genes <- c("Flg")
dds_sub <- estimateSizeFactors(dds_sub)
norm_counts <- counts(dds_sub,normalized=TRUE)
write.csv(as.data.frame(norm_counts),file="Normalized_counts.csv",row.names = TRUE)


library(EnsDb.Mmusculus.v79)


# DESeq2

dds1 <- DESeq(dds_sub)
resultsNames(dds1)
res <- results(dds1)

# LFC shrinkage 
reslfc <- lfcShrink(dds1, coef=resultsNames(dds1)[2],type="apeglm")
res_df<- reslfc %>% as.data.frame
res_df$stat <- res$stat



# Volcano Plot
keyvals <- ifelse(
  !is.na(res_df$padj) & res_df$padj <0.05 & res_df$log2FoldChange>1.0 ,"#e41a1c",
  ifelse(!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange< -1.0, "#377eb8",
         "#999999"))

names(keyvals)[keyvals=="#e41a1c"] <- "Up-regulated"
names(keyvals)[keyvals=="#377eb8"] <- "Down-regulated"
names(keyvals)[keyvals=="#999999"] <- "NS"

library(EnhancedVolcano)
EnhancedVolcano(
  res_df,
  lab = rownames(res_df),
  x = 'log2FoldChange',
  y = 'padj',
  pCutoff = 0.05,
  FCcutoff = 1.0,
  colCustom = keyvals,
  ylab = bquote(~-Log[10] ~ italic(P)[adj]),
  title = 'Normal vs AD',
  caption = 'FC cutoff: 1.0, Adjusted p-value cutoff: 0.05',
  legendPosition = 'right'
)



# Up-regulated genes ranked by STAT
up_genes <- res_df[res_df$log2FoldChange > 1 & res_df$padj < 0.05&!is.na(res_df$padj),]
nrow(up_genes)
up_gene_ranked <- up_genes[order(up_genes$stat,decreasing=TRUE),]


# Down-regulated genes ranked by STAT
down_genes <- res_df[res_df$log2FoldChange < -1 & res_df$padj < 0.05&!is.na(res_df$padj),]
nrow(down_genes)
down_genes_ranked <- down_genes[order(down_genes$stat),]
down_genes_ranked <- down_genes %>%
  mutate(
    rank_l2fc = rank(-log2FoldChange), # L2FC가 클수록 높은 순위 (1위, 2위...)
    rank_padj = rank(padj)             # padj가 작을수록(유의미할수록) 높은 순위
  )


# DEG selection
sig_genes <- rownames(res_df)[which(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1& !is.na(res_df$padj))]


# DEG 발현 패턴 heatmap으로 시각화
vsd <- vst(dds1,blind=FALSE)
mat <- assay(vsd)[sig_genes,]
                         
colnames(vsd) <- gsub("^quant_", "", colnames(vsd))
colnames(mat) <- gsub("^quant_", "", colnames(mat))
annot_info <- data.frame(
  group = vsd$group,
  row.names = colnames(vsd)
)

my_colors <- colorRampPalette(c("navy", "white", "red"))(255)
pheatmap(mat,
         scale = "row",                       # 유전자별 Z-score 계산 (패턴 비교의 핵심!)
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         annotation_col = annot_info,
         show_rownames = FALSE,              # 유전자가 많으면 FALSE, 상위 30~50개면 TRUE
         show_colnames = TRUE,
         color = my_colors,
         main = "Gene Expression Heatmap of DEGs (Normal vs AD)")



# GO Analysis
library(clusterProfiler)
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

# up gene GO
up_gene_symbol <- rownames(up_genes)
up_entrez <- bitr(up_gene_symbol,
                  fromType="SYMBOL",
                  toType="ENTREZID",
                  OrgDb = org.Mm.eg.db
                  )
ego_up <- enrichGO(gene=up_entrez$ENTREZID,
                   OrgDb = org.Mm.eg.db,
                   ont="BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)
head(ego_up)
dotplot(ego_up,showCategory=15, title="GO BP Enrichment (Up-regulated)")




# down-regulated genes GO
down_gene_symbol <- rownames(down_genes)
down_entrez <- bitr(down_gene_symbol,
                    fromType = "SYMBOL",
                    toType="ENTREZID",
                    OrgDb = org.Mm.eg.db)
ego_down <- enrichGO(gene=down_entrez$ENTREZID,
                     OrgDb=org.Mm.eg.db,
                     ont="BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)

dotplot(ego_down,showCategory=15,title="GO BP Enrichment (Down-regulated)")



# log file 
writeLines(capture.output(sessionInfo()), "logs/sessionInfo.txt")

