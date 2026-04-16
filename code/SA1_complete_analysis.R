
# =============================================================
# P6 SA1 — COMPLETE ANALYSIS SCRIPT
# SARS-CoV-2 Brain Senescence and Age-Dependent Cognitive Impairment
# de Erausquin Lab | UTHSCSA Glenn Biggs Institute
# Author: Kelechi Elechi
# Date: 2025
# =============================================================

# ── 0. Setup ─────────────────────────────────────────────
library(here)
setwd("/Users/kelechielechi/Documents/COVID PASC")
here::i_am(".here")

library(GEOquery); library(tidyverse); library(tximport)
library(DESeq2); library(edgeR); library(GSVA)
library(lme4); library(lmerTest); library(vegan)
library(org.Hs.eg.db); library(AnnotationDbi)
library(ComplexHeatmap); library(circlize)
library(ggplot2); library(here); library(janitor)
set.seed(42)

# ── 1. Load metadata ──────────────────────────────────────
gse188847 <- getGEO("GSE188847", GSEMatrix=TRUE,
                     destdir=here("data/raw/GEO/GSE163005"))
meta <- pData(gse188847[[1]])

meta_clean <- meta %>%
  filter(source_name_ch1 == "Frontal cortex of brain") %>%
  dplyr::select(geo_accession,
    age=`age:ch1`, group=`group:ch1`,
    sex=`Sex:ch1`, source=source_name_ch1) %>%
  mutate(
    age          = as.numeric(age),
    sex_male     = as.integer(sex == "M"),
    age_group_60 = factor(ifelse(age>=60,">=60","<60"),
                          levels=c("<60",">=60")),
    diagnosis    = case_when(
      group %in% c("COVID-19","ICU/VENT") ~ "COVID19",
      group == "Control"                  ~ "Control",
      group == "Old"                      ~ "Old_Control",
      group == "Young"                    ~ "Young_Control",
      TRUE ~ NA_character_),
    batch = ifelse(grepl("^GSM5690",geo_accession),
                   "Batch1","Batch2")
  )

write.csv(meta_clean,
          here("data/processed/GSE188847_metadata_clean.csv"),
          row.names=FALSE)

# ── 2. Load expression data (tximport) ────────────────────
tx2gene_raw <- read.csv(
  "/Users/kelechielechi/Documents/covid-brain-aging-thesis/data/processed/GSE188847/tx2gene.csv")
tx2gene_clean <- tx2gene_raw[,1:2]
colnames(tx2gene_clean) <- c("tx","gene")
tx2gene_clean$tx <- sub("\\.\\d+$","",tx2gene_clean$tx)

salmon_dir <- "/Users/kelechielechi/Documents/GSE188847/"
all_gz <- list.files(salmon_dir,pattern="\\.gz$",full.names=TRUE)
sample_names <- basename(all_gz) %>%
  sub("^GSM\\d+_","",.) %>%
  sub("\\.salmon.*\\.gz$","",.)
names(all_gz) <- sample_names

txi <- tximport(all_gz, type="salmon",
                tx2gene=tx2gene_clean,
                ignoreTxVersion=TRUE,
                countsFromAbundance="lengthScaledTPM")

# ── 3. Build DESeq2 object ────────────────────────────────
geo_titles <- data.frame(
  geo_accession = meta$geo_accession,
  sample_name   = meta$title,
  stringsAsFactors=FALSE)

sample_table <- geo_titles %>%
  left_join(meta_clean, by="geo_accession") %>%
  filter(!is.na(diagnosis))

keep_samples <- sample_table$sample_name
txi_fc <- list(
  counts    = txi$counts[,keep_samples],
  abundance = txi$abundance[,keep_samples],
  length    = txi$length[,keep_samples],
  countsFromAbundance = txi$countsFromAbundance)

sample_table <- sample_table[
  match(colnames(txi_fc$counts), sample_table$sample_name),]
rownames(sample_table) <- sample_table$sample_name

dds <- DESeqDataSetFromTximport(txi_fc, sample_table,
                                 design=~diagnosis)

saveRDS(dds,          here("data/processed/dds_raw.rds"))
saveRDS(sample_table, here("data/processed/sample_table.rds"))
saveRDS(txi_fc,       here("data/processed/txi_fc.rds"))

# ── 4. QC and gene filtering ──────────────────────────────
keep <- filterByExpr(counts(dds), group=dds$diagnosis,
                     min.count=10, min.prop=0.2)
dds_filtered <- dds[keep,]

vst_data <- vst(dds_filtered, blind=TRUE)

# Batch PERMANOVA
dist_obj   <- dist(t(assay(vst_data)))
batch_vec  <- sample_table$batch[match(colnames(vst_data),
                                        sample_table$sample_name)]
diag_vec   <- sample_table$diagnosis[match(colnames(vst_data),
                                            sample_table$sample_name)]
set.seed(42)
permanova  <- adonis2(dist_obj ~ batch_vec + diag_vec,
                      permutations=999, by="terms")
# Batch R2=0.020, p=0.072 — no batch correction needed

# ── 5. Gene symbol conversion ─────────────────────────────
id_map <- AnnotationDbi::select(org.Hs.eg.db,
  keys=rownames(dds_filtered),
  columns=c("ENSEMBL","SYMBOL"), keytype="ENSEMBL")

id_map_clean <- id_map %>%
  filter(!is.na(SYMBOL), SYMBOL!="") %>%
  dplyr::group_by(ENSEMBL) %>%
  dplyr::slice(1) %>% ungroup()

matched_idx <- match(rownames(dds_filtered),id_map_clean$ENSEMBL)
has_symbol  <- !is.na(matched_idx)
new_names   <- rownames(dds_filtered)
new_names[has_symbol] <- id_map_clean$SYMBOL[matched_idx[has_symbol]]
new_names_unique <- make.unique(new_names, sep=".")

dds_symbol <- dds_filtered
rownames(dds_symbol) <- new_names_unique

saveRDS(dds_symbol, here("data/processed/dds_symbol.rds"))
saveRDS(id_map_clean, here("data/processed/ensembl_to_symbol.rds"))

# ── 6. Load gene sets ─────────────────────────────────────
senmayo_genes <- read.csv(here("data/raw/genesets/SenMayo_genes.csv"))$gene
cellage_genes <- read.csv(here("data/raw/genesets/CellAge_genes.csv"))$gene
genage_genes  <- read.csv(here("data/raw/genesets/GenAge_genes.csv"))$gene
expressed_symbols <- rownames(dds_symbol)

gene_sets <- list(
  SenMayo = senmayo_genes[senmayo_genes %in% expressed_symbols],
  CellAge = cellage_genes[cellage_genes %in% expressed_symbols],
  GenAge  = genage_genes[genage_genes   %in% expressed_symbols])

# ── 7. ssGSEA senescence scoring ─────────────────────────
vst_data <- vst(dds_symbol, blind=TRUE)
expr_mat  <- assay(vst_data)

gsva_param <- gsvaParam(exprData=expr_mat,
                         geneSets=gene_sets,
                         kcdf="Gaussian",
                         minSize=10, maxSize=1000)
set.seed(42)
sen_scores <- gsva(gsva_param, verbose=TRUE)

scores_df <- as.data.frame(t(sen_scores)) %>%
  rownames_to_column("sample_name") %>%
  left_join(sample_table, by="sample_name") %>%
  mutate(SenMayo_z=scale(SenMayo)[,1],
         CellAge_z=scale(CellAge)[,1],
         GenAge_z =scale(GenAge)[,1])

write.csv(scores_df,
          here("results/tables/senescence_scores.csv"),
          row.names=FALSE)

# ── 8. Statistical models ─────────────────────────────────
# Primary: COVID19 vs Control, additive model (all 73 samples)
scores_df$diagnosis <- factor(scores_df$diagnosis,
  levels=c("Control","COVID19","Old_Control","Young_Control"))

results_secondary <- lapply(c("SenMayo_z","CellAge_z","GenAge_z"),
  function(sig) {
    fit <- lm(as.formula(paste0(sig,"~ age + diagnosis + sex_male")),
              data=scores_df)
    as.data.frame(coef(summary(fit))) %>%
      rownames_to_column("term") %>%
      mutate(signature=sig,
             p_adj=p.adjust(`Pr(>|t|)`,method="BH"))
  }) %>% bind_rows()
# KEY RESULT: SenMayo COVID19 beta=+0.954, p_adj=0.0022

# Severity gradient ANOVA
scores_severity <- scores_df %>%
  filter(diagnosis %in% c("Control","COVID19") |
         grepl("ICUVENT",sample_name)) %>%
  mutate(severity=case_when(
    diagnosis=="Control" ~ "Control",
    grepl("ICUVENT",sample_name) ~ "ICU/VENT",
    diagnosis=="COVID19" ~ "COVID-19"))
aov_fit <- aov(SenMayo_z ~ severity, data=scores_severity)
# KEY RESULT: F=9.223, p=0.0004; Control < COVID-19 < ICU/VENT

# ── 9. DESeq2 differential expression ────────────────────
dds_cv <- dds_symbol[,dds_symbol$diagnosis %in%
                       c("Control","COVID19")]
dds_cv$diagnosis <- relevel(droplevels(dds_cv$diagnosis),
                              ref="Control")
design(dds_cv) <- ~ sex_male + age + diagnosis
dds_cv  <- DESeq(dds_cv)
res_cv  <- results(dds_cv,
                    name="diagnosis_COVID19_vs_Control",
                    alpha=0.05)
# KEY RESULT: 6694 DEGs; 62/93 SenMayo genes upregulated

write.csv(as.data.frame(res_cv) %>% rownames_to_column("gene"),
          here("results/tables/DESeq2_COVID_vs_Control_full.csv"),
          row.names=FALSE)

