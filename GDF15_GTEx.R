library(tidyverse)
library(edgeR)
library(corrr)

# read data
raw_gene_counts <- read.delim(gzfile("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"), skip = 2)
gene_symbols <- raw_gene_counts[,1:2]

gene_counts <- raw_gene_counts %>%
  column_to_rownames("Name")
gene_counts <- gene_counts[,-1]


# read metadata and prepare for join with expression data
annotations_subject_phenotypes <- read_tsv("GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.tsv")
annotations_sample_attributes <- read_tsv("GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.tsv") %>%
  filter(SMAFRZE == "RNASEQ") %>%
  mutate(SUBJID = substring(SAMPID,1,10))

annotations_sample_attributes <- annotations_sample_attributes %>%
  mutate(SUBJID_tmp = case_when(
    substr(SUBJID, nchar(SUBJID), nchar(SUBJID)) == "-" ~ substr(SUBJID, 1, nchar(SUBJID)-1),
    .default = SUBJID
  )) %>%
  select(-SUBJID) %>%
  rename(SUBJID = SUBJID_tmp)

annotations_sample_attributes$SAMPID <- gsub("\\-", ".", annotations_sample_attributes$SAMPID)


# TMM normalization
dge <- DGEList(gene_counts[,annotations_sample_attributes$SAMPID], group = factor(annotations_sample_attributes$SMTSD))
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
exprs <- cpm(dge, log = FALSE)

# join with gene symbols
exprs <- exprs %>%
  as.data.frame() %>%
  rownames_to_column("Name") %>%
  full_join(gene_symbols, by = "Name") %>%
  na.omit() %>%
  select(-Name) %>%
  rename("Gene" = Description)

# filter for GDF15
exprs_gdf15 <- exprs %>%
  as.data.frame() %>%
  filter(Gene %in% c("GDF15")) %>%
  column_to_rownames("Gene") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SAMPID") %>%
  unique()


# join with sample metadata
exprs_gdf15 <- exprs_gdf15 %>%
  full_join(annotations_sample_attributes, by = "SAMPID") %>%
  select(GDF15, SMTSD, SUBJID) %>%
  unique() 

# join with subject phenotypes
exprs_gdf15 <- exprs_gdf15 %>%
  full_join(annotations_subject_phenotypes, by = "SUBJID") %>%
  select(GDF15, SMTSD, SUBJID, AGE) %>%
  unique() %>% 
  na.omit() %>%
  pivot_wider(names_from = SMTSD, values_from = GDF15) %>%
  select(- `Cells - Cultured fibroblasts`, - `Cells - EBV-transformed lymphocytes`, - `Cervix - Ectocervix`,
         - `Cervix - Endocervix`, - `Kidney - Medulla`, - `Fallopian Tube`)


# calculate tissue means
tissue_means_gdf15 <- exprs_gdf15 %>%
  select(-AGE) %>%
  column_to_rownames("SUBJID") %>%
  colMeans(na.rm = T) %>%
  as.data.frame() %>%
  rownames_to_column("Tissue")  %>%
  arrange(desc(`.`)) %>%
  rename("Mean" = ".")


# compute correlation with age in each tissue
gdf15_age_cor <- exprs_gdf15 %>%
  correlate(.,method = "spearman", use = "pairwise.complete.obs") %>%
  focus(AGE) %>%
  arrange(desc(AGE)) %>%
  rename("Tissue" = "term", "Spearman r - Age" = "AGE") 

