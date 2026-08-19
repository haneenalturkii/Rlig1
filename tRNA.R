library("readxl")
library("QFeatures")
library("tidyr")
library("dplyr")
library("ggplot2")
library("patchwork")
library("tidyverse")
library("gridExtra")
library("grid")
library("ComplexHeatmap")
library("magick")
library("circlize")
library("imputeLCMD")
library("FactoMineR")
library("factoextra")
library("limma")
library("ggrepel")
library("clusterProfiler")
library("org.Hs.eg.db")
library("RColorBrewer")
library("matrixStats")
library("STRINGdb")
library("UpSetR")

setwd("~/Desktop")

###### Reading data #####

prot_df <- read.csv("20240710_145530_240710_FS_WithoutRNA_Paper_Report (1).csv", 
                    check.names = FALSE)
names(prot_df)

sample_metadata <- read.delim("metadata.tsv")
print(sample_metadata)

### Shortening ids
# Keep only the Quantity columns.

id_cols <- c("PG.ProteinGroups", "PG.Genes")
quant_cols <- grep("\\.raw\\.PG\\.Quantity$", names(prot_df), value = TRUE)
quant_cols
prot_long <- prot_df %>%
  dplyr::select(all_of(id_cols), all_of(quant_cols)) %>%
  pivot_longer(
    cols = all_of(quant_cols),
    names_to = "sample_col",
    values_to = "quantity"
  ) %>%
  
  mutate(
    sample_id = sample_col %>%
      str_remove("^\\[\\d+\\]\\s*") %>%
      str_remove("\\.raw\\.PG\\.Quantity$")
  )

# Join metadata on sample_id
prot_annotated <- prot_long %>%
  left_join(sample_metadata, by = "sample_id")



##### Reading to a QFeatures object #####

prot_qf <- readQFeatures(
  assayData = prot_df,
  quantCols = 46:63,
  fnames = "PG.ProteinGroups",
  name = "proteins"
)

sid <- rownames(colData(prot_qf)) %>% str_remove("\\.raw\\.PG\\.Quantity$")

### find matching row in the metadata for each id
idx <- match(sid, sample_metadata$sample_id)

### for safety, stop if any of the 18 didn't match
stopifnot(!anyNA(idx))

### copy the labels into colData in their original order
colData(prot_qf)$sample_id <- sid
colData(prot_qf)$group     <- sample_metadata$group[idx]
colData(prot_qf)$stress    <- sample_metadata$stress[idx]

# check
colData(prot_qf)

unique(colData(prot_qf)$group)
unique(colData(prot_qf)$stress)

### Shorten the stress label
colData(prot_qf)$stress_short <- dplyr::recode(colData(prot_qf)$stress,
                                               "Stressed (40 µM menadione, 1 h)" = "str",
                                               "No stress" = "unst")
# Build the combined condition label from group + short stress.
colData(prot_qf)$condition <- paste(colData(prot_qf)$group,
                                    colData(prot_qf)$stress_short,
                                    sep = "_")

# Sanity check: 6 conditions, 3 samples each
table(colData(prot_qf)$condition)

##### Remove contaminants #####
# 57 contaminant proteins (trypsin, keratin, BSA, casein, etc.) are flagged 
prot_ids <- rownames(prot_qf[["proteins"]])
is_contaminant <- grepl("^Cont_", prot_ids)   # matches the "Cont_" prefix exactly
sum(is_contaminant)                            # expect 57

prot_qf <- prot_qf[!is_contaminant, , ]        # drop them from the whole object
dim(prot_qf[["proteins"]])                    

##### Dealing with missing data #####

prot_qf <- addAssay(prot_qf, prot_qf[["proteins"]], name = "prot_with_missing")
prot_qf <- addAssayLink(prot_qf, from = "proteins", to = "prot_with_missing")
plot(prot_qf)
na_report <- nNA(prot_qf, i = "prot_with_missing")
print(na_report)

na_report$nNAcols %>%
  as_tibble() %>%
  ggplot(aes(x = name, y = pNA)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Proportion missing values") +
  coord_flip() +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic()


##### Log transform #####

prot_qf <- logTransform(object = prot_qf,
                        base   = 2,
                        i = "prot_with_missing",
                        name = "log_proteins")


##### Imputing #####
set.seed(118)

prot_qf <- QFeatures::impute(object = prot_qf,
                             method = "knn",
                             MARGIN = 1,
                             i = "log_proteins",
                             name = "log2_imp_prot")

## check with histohram 
obs <- assay(prot_qf[["log_proteins"]])   
imp <- assay(prot_qf[["log2_imp_prot"]])     

was_na <- is.na(obs)  

hist(imp[!was_na], breaks = 50, col = rgb(0,0,1,0.4),
     main = "Observed (blue) vs imputed (red)", xlab = "log2 quantity")
hist(imp[was_na],  breaks = 50, col = rgb(1,0,0,0.4), add = TRUE)



##### Normalization #####

prot_qf <- normalize(prot_qf,
                     i      = "log2_imp_prot",
                     name   = "log_norm_prot",
                     method = "diff.median")



before_norm <- prot_qf[["log2_imp_prot"]]   %>% 
  assay() %>% 
  longForm() %>% 
  mutate(Condition = strsplit(as.character(colname), split = "_") %>% 
           sapply("[[", 1)) %>% 
  ggplot(aes(x = colname, y = value, fill = Condition)) + 
  geom_boxplot() + guides(x = guide_axis(angle = 90)) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text.y=element_blank(), legend.position = "none")

after_norm <- prot_qf[["log_norm_prot"]] %>% 
  assay() %>% 
  longForm() %>% 
  mutate(Condition = strsplit(as.character(colname), split = "_") %>% 
           sapply("[[", 1)) %>% 
  ggplot(aes(x = colname, y = value, fill = Condition)) + 
  geom_boxplot() + guides(x = guide_axis(angle = 90)) + 
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_blank(), 
        axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        legend.position = "none")


norm_plots <- grid.arrange(before_norm, after_norm, ncol= 2)

## Other way to see it 

pre_norm <- longForm(prot_qf[,, "log2_imp_prot"], colvars = 'group') %>%
  ggplot(aes(x = value, colour = group, group = colname)) +
  geom_density() +
  theme_classic() +
  xlab("log2 (Abundance)") +
  ggtitle('Pre-normalisation')

post_norm <- longForm(prot_qf[,, 'log_norm_prot'], colvars = 'group') %>%
  ggplot(aes(x = value, colour = group, group = colname)) +
  geom_density() +
  theme_classic() +
  xlab("log2 (Abundance)") +
  ggtitle('Post-normalisation')

pre_norm + post_norm

##### PCA #####

# Build PCA object
PCA_obj <- PCA(
  X = prot_qf[["log_norm_prot"]] %>%     
    assay() %>%
    as.data.frame() %>%
    t(),                                
  ncp        = 10,       
  scale.unit = TRUE,     
  graph      = FALSE
)

# scree plot to see how much variance each PC explains
PCA_scree <- fviz_screeplot(
  PCA_obj,
  choice    = "variance",
  geom      = c("bar", "line"),
  barfill   = "dodgerblue4",
  addlabels = TRUE,
  ncp       = 8,
  main      = "Scree plot"
)
PCA_scree


# PCA plot 

scores <- as.data.frame(PCA_obj$ind$coord)
scores$group  <- colData(prot_qf)$group         
scores$stress <- colData(prot_qf)$stress_short

# Naming the stress legend understandable names
scores$stress <- dplyr::recode(scores$stress,
                               "str"  = "Stressed",
                               "unst" = "Unstressed"
)

# plot
ggplot(scores, aes(x = Dim.1, y = Dim.2, colour = group, shape = stress)) +
  geom_point(size = 3, stroke = 1.5) +
  scale_colour_brewer(palette = "Dark2") +
  scale_shape_manual(values = c("Stressed" = 4, "Unstressed" = 16)) + 
  labs(title = "PC1 vs PC2",
       colour = "Bait",         
       shape  = "Stress state",  
       x = "Dim1 (35.6%)", y = "Dim2 (19.7%)") +
  theme_bw()


##### Statistical Analysis #####

## limma
# Set the level order 
condition <- factor(colData(prot_qf)$condition)
levels(condition)  

# Build the design
limma_design <- model.matrix(~ 0 + condition)

# Rename columns to the bare condition names 
colnames(limma_design) <- levels(condition)

# Verify: 18 rows (samples), 6 columns (conditions), a single 1 per row.
limma_design

# Contrasts 
colnames(limma_design) <- make.names(colnames(limma_design))
colnames(limma_design) 

contrasts_matrix <- makeContrasts(
  AMP_vs_Beads_unst = AMP.Rlig1_unst - Bead.control_unst,
  AMP_vs_Beads_str  = AMP.Rlig1_str  - Bead.control_str,
  K57A_vs_Beads_unst = K57A.Rlig1_unst - Bead.control_unst,
  K57A_vs_Beads_str  = K57A.Rlig1_str  - Bead.control_str,
  AMP_vs_K57A_unst = AMP.Rlig1_unst - K57A.Rlig1_unst,
  AMP_vs_K57A_str  = AMP.Rlig1_str  - K57A.Rlig1_str,
  levels = limma_design
)

contrasts_matrix

# Fit model
expr <- assay(prot_qf[["log_norm_prot"]])
dim(expr)

fit  <- lmFit(expr, design = limma_design)      # fit the model
fit2 <- contrasts.fit(fit, contrasts_matrix)           # apply my 6 contrasts
fit2 <- eBayes(fit2, robust = TRUE, trend = TRUE)

# Results as a summary across all six contrasts 
contrast_names <- colnames(contrasts_matrix)

summary_table <- sapply(contrast_names, function(cn) {
  res <- topTable(fit2, coef = cn, number = Inf, adjust.method = "BH")
  c(
    up   = sum(res$adj.P.Val < 0.05 & res$logFC >  1),   # enriched in first group
    down = sum(res$adj.P.Val < 0.05 & res$logFC < -1),   # enriched in second group
    total_sig = sum(res$adj.P.Val < 0.05)                 # any significant
  )
})

t(summary_table)

##### Full result tables per contrast #####

# Protein annotation
uniprot_annotations <- rowData(prot_qf[["log_norm_prot"]]) %>%
  as.data.frame() %>%
  dplyr::select(PG.Genes, PG.ProteinDescriptions)

# loop using my object names

limma_results_all <- vector("list", length(contrast_names))
names(limma_results_all) <- contrast_names

for (cn in contrast_names) {
  limma_results_all[[cn]] <- topTable(
    fit2,                 
    coef    = cn,
    number  = Inf,
    adjust.method = "BH",
    sort.by = "p"
  ) %>%
    data.frame() %>%
    rownames_to_column("UniprotID") %>%
    merge(uniprot_annotations, by.x = "UniprotID", by.y = "row.names",
          all.x = TRUE) %>%     # all.x = TRUE keeps all result rows
    mutate(contrast = cn) %>%
    arrange(P.Value)
}

limma_results_all_df <- dplyr::bind_rows(limma_results_all)

# Check

all_results_report <- table(
  status = ifelse(limma_results_all_df$adj.P.Val < 0.05 &
                    abs(limma_results_all_df$logFC) > 1,
                  ifelse(limma_results_all_df$logFC > 0, "Enriched_up",
                         "Enriched_down"),
                  "Not_sig"),
  contrast = limma_results_all_df$contrast
)
all_results_report


##### Checking step: Histogram for P-value distribution across all contrast

limma_results_all_df %>%
  ggplot(aes(P.Value)) +
  geom_histogram(bins = 50) +
  theme_classic() +
  labs(x = "P-value", y = "Frequency") +
  facet_wrap(~contrast)


##### GO enrichment #####

# Four contrasts 
contrasts_for_go <- c("AMP_vs_Beads_unst", "AMP_vs_Beads_str",
                      "K57A_vs_Beads_unst", "K57A_vs_Beads_str")

go_results <- list()

for (cn in contrasts_for_go) {
  
  # enriched (up) proteins for this contrast
  sig_up <- limma_results_all_df %>%
    filter(contrast == cn, adj.P.Val < 0.05, logFC > 1) %>%
    pull(UniprotID)
  
  # background = all proteins quantified in this contrast
  universe <- limma_results_all_df %>%
    filter(contrast == cn) %>%
    pull(UniprotID)
  
  go_results[[cn]] <- enrichGO(
    gene     = sig_up,
    universe = universe,
    OrgDb    = org.Hs.eg.db,
    keyType  = "UNIPROT",
    ont      = "ALL",           
    pvalueCutoff = 0.05,
    readable = TRUE
  )
  message(cn, ": ", length(sig_up), " enriched proteins analysed")
}


# take the AMP vs Beads result and keep only MF terms
go_mf <- go_results[["AMP_vs_Beads_unst"]]
go_mf@result <- go_mf@result %>% filter(ONTOLOGY == "MF")

dotplot(go_mf, x = "Count",
        showCategory = 12,
        color = "p.adjust",
        font.size = 13,
        label_format = 40) +
  ggtitle("GO Molecular Function: AMP vs Beads (unstressed)") +
  theme(axis.text.y  = element_text(size = 13),
        axis.text.x  = element_text(size = 12),
        axis.title.x = element_text(size = 13),
        plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text  = element_text(size = 11),
        legend.title = element_text(size = 11))

ggsave("~/Desktop/go_dotplot_MF.png", width = 9, height = 7, dpi = 300, bg = "white")




########## tRNA-PROCESSING MACHINERY ANALYSIS ##########

## Define the tRNA machinery, grouped by pathway step #####
trna_machinery <- tibble::tribble(
  ~step,                    ~gene,
  "1. Endonuclease (cut)",  "TSEN2",
  "1. Endonuclease (cut)",  "TSEN15",
  "1. Endonuclease (cut)",  "TSEN34",
  "1. Endonuclease (cut)",  "TSEN54",
  "2. 5' processing (RNase P)", "POP1",
  "2. 5' processing (RNase P)", "POP4",
  "2. 5' processing (RNase P)", "POP5",
  "2. 5' processing (RNase P)", "POP7",
  "2. 5' processing (RNase P)", "RPP14",
  "2. 5' processing (RNase P)", "RPP21",
  "2. 5' processing (RNase P)", "RPP25",
  "2. 5' processing (RNase P)", "RPP25L",
  "2. 5' processing (RNase P)", "RPP30",
  "2. 5' processing (RNase P)", "RPP38",
  "2. 5' processing (RNase P)", "RPP40",
  "3. 3' processing",       "ELAC2",
  "3. 3' processing",       "ELAC1",
  "4. End-healing",         "CLP1",
  "4. End-healing",         "RTCA",
  "4. End-healing",         "ANGEL2",
  "5. Ligation (RTCB cplx)","RTCB",
  "5. Ligation (RTCB cplx)","DDX1",
  "5. Ligation (RTCB cplx)","FAM98A",
  "5. Ligation (RTCB cplx)","FAM98B",
  "5. Ligation (RTCB cplx)","RTRAF",
  "5. Ligation (RTCB cplx)","C2orf49",
  "6. Maturation (CCA)",    "TRNT1"
)

### Pull each enzyme's result from limma output for AMP vs beads
trna_results <- limma_results_all_df %>%
  filter(contrast == "AMP_vs_Beads_unst",
         PG.Genes %in% trna_machinery$gene) %>%
  dplyr::select(PG.Genes, logFC, adj.P.Val) %>%
  right_join(trna_machinery, by = c("PG.Genes" = "gene")) %>%
  arrange(step, adj.P.Val)
head(trna_results, 30)

cat("\nSignificantly enriched (FDR<0.05, logFC>1): ",
    sum(trna_results$adj.P.Val < 0.05 & trna_results$logFC > 1, na.rm = TRUE),
    " of ", sum(!is.na(trna_results$logFC)), " detected\n")

### tRNA machinery under K57A vs Beads
trna_k57a <- limma_results_all_df %>%
  filter(contrast == "K57A_vs_Beads_unst",
         PG.Genes %in% trna_machinery$gene) %>%
  dplyr::select(PG.Genes, logFC, adj.P.Val) %>%
  arrange(adj.P.Val)
head(trna_k57a, 30)

cat("K57A enriched:", sum(trna_k57a$adj.P.Val < 0.05 & trna_k57a$logFC > 1, na.rm=TRUE),
    "of", sum(!is.na(trna_k57a$logFC)), "\n")

### check the same machinery under stress 
# Does the tRNA module hold under stress too? (supports "stress-robust")
trna_stress <- limma_results_all_df %>%
  filter(contrast == "AMP_vs_Beads_str",
         PG.Genes %in% trna_machinery$gene) %>%
  dplyr::select(PG.Genes, logFC, adj.P.Val) %>%
  arrange(adj.P.Val)
cat("\n--- Same machinery, STRESSED ---\n")
head(trna_stress, 30)


##### FIGURE 1: AMP vs beads pathway #####

plot_trna <- trna_results %>%
  filter(!is.na(logFC)) %>%
  mutate(
    sig = ifelse(adj.P.Val < 0.05 & logFC > 1, "Enriched (FDR<0.05)", "Not significant"),
    step = factor(step, levels = c(
      "1. Endonuclease (cut)", "2. 5' processing (RNase P)",
      "3. 3' processing", "4. End-healing",
      "5. Ligation (RTCB cplx)", "6. Maturation (CCA)"))
  )

fig_bars <- ggplot(plot_trna,
                   aes(x = reorder(PG.Genes, logFC), y = logFC, fill = sig)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  coord_flip() +
  facet_grid(step ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(
    values = c("Enriched (FDR<0.05)" = "#B2182B",
               "Not significant"     = "#B0B0B0"),
    name = NULL) +
  labs(
    title = "Rlig1 co-enriches with the tRNA-repair machinery",
    subtitle = "The tRNA-repair pathway co-enriches with active Rlig1",
    x = NULL,
    y = "Enrichment over bead control (log2 fold change)") +
  theme_bw(base_size = 14) +
  theme(
    strip.text.y.left = element_text(angle = 0, hjust = 0.5, face = "bold", size = 13),
    strip.placement   = "outside",
    plot.title        = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle     = element_text(size = 13, colour = "grey30", hjust = 0.5),
    legend.position   = "top",
    legend.text       = element_text(size = 13),          
    axis.text.y       = element_text(size = 12),         
    axis.text.x       = element_text(size = 12),          
    axis.title.x      = element_text(size = 14),
    panel.grid.major.y = element_blank())

fig_bars
ggsave("~/Desktop/trna_machinery.png", fig_bars,
       width = 10, height = 9, dpi = 300, bg = "white")


##### FIGURE 2: K57A vs Beads pathway #####

# Prep K57A data with pathway steps and significance flag
plot_trna_k57 <- trna_k57a %>%
  left_join(trna_machinery, by = c("PG.Genes" = "gene")) %>%
  filter(!is.na(logFC)) %>%
  mutate(
    sig = ifelse(adj.P.Val < 0.05 & logFC > 1, "Enriched (FDR<0.05)", "Not significant"),
    step = factor(step, levels = c(
      "1. Endonuclease (cut)", "2. 5' processing (RNase P)",
      "3. 3' processing", "4. End-healing",
      "5. Ligation (RTCB cplx)", "6. Maturation (CCA)"))
  )

fig_bars_k57 <- ggplot(plot_trna_k57,
                       aes(x = reorder(PG.Genes, logFC), y = logFC, fill = sig)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  coord_flip() +
  facet_grid(step ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(
    values = c("Enriched (FDR<0.05)" = "#B2182B",
               "Not significant"     = "#B0B0B0"),
    name = NULL) +
  labs(
    title = "The tRNA machinery is also enriched by inactive K57A-Rlig1",
    subtitle = "The catalytically dead mutant recovers the same pathway, at lower fold change",
    x = NULL,
    y = "Enrichment over bead control (log2 fold change)") +
  theme_bw(base_size = 14) +
  theme(
    strip.text.y.left = element_text(angle = 0, hjust = 0.5, face = "bold", size = 13),
    strip.placement   = "outside",
    plot.title        = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle     = element_text(size = 13, colour = "grey30", hjust = 0.5),
    legend.position   = "top",
    legend.text       = element_text(size = 13),          
    axis.text.y       = element_text(size = 12),         
    axis.text.x       = element_text(size = 12),          
    axis.title.x      = element_text(size = 14),
    panel.grid.major.y = element_blank())

fig_bars_k57
ggsave("~/Desktop/trna_machinery_k57a.png", fig_bars_k57,
       width = 10, height = 9, dpi = 300, bg = "white")

##### Targeted enrichment test: is the tRNA machinery over-represented? #####

# Define the full tRNA-processing gene set
trna_set <- c("TSEN2","TSEN15","TSEN34","TSEN54","POP1","POP4","POP5","POP7",
              "RPP14","RPP21","RPP25","RPP25L","RPP30","RPP38","RPP40",
              "ELAC2","ELAC1","CLP1","RTCA","ANGEL2","RTCB","DDX1",
              "FAM98A","FAM98B","RTRAF","C2orf49","TRNT1")

# Take limma results for AMP vs Beads (unstressed), testable proteins only
res <- limma_results_all_df %>%
  filter(contrast == "AMP_vs_Beads_unst", !is.na(adj.P.Val))

# Flag each protein: is it tRNA machinery? is it enriched?
res <- res %>%
  mutate(
    is_trna     = PG.Genes %in% trna_set,
    is_enriched = adj.P.Val < 0.05 & logFC > 1
  )

# Build the 2x2 contingency table
tab <- table(
  tRNA     = factor(res$is_trna,     levels = c(TRUE, FALSE)),
  enriched = factor(res$is_enriched, levels = c(TRUE, FALSE))
)
print(tab)

# Fisher's exact test (one-sided: is tRNA OVER-represented?)
ft <- fisher.test(tab, alternative = "greater")
cat("\ntRNA enzymes enriched:",
    tab["TRUE","TRUE"], "of", sum(tab["TRUE",]),
    "(", round(100*tab["TRUE","TRUE"]/sum(tab["TRUE",])), "%)\n")
cat("Background enrichment rate:",
    round(100*tab["FALSE","TRUE"]/sum(tab["FALSE",])), "%\n")
cat("Odds ratio:", round(ft$estimate, 1), "\n")
cat("Fisher exact p:", format(ft$p.value, scientific = TRUE, digits = 3), "\n")



##### heatmap of tRNA machinery genes accross all 18 samples #####

## Define the tRNA machinery genes 
trna_genes <- c("TSEN34","TSEN54","POP1","POP4","POP7","RPP14","RPP30",
                "RPP38","RPP40","RPP25L","ELAC2","CLP1","RTCA","ANGEL2",
                "RTCB","DDX1","FAM98A","FAM98B","RTRAF","C2orf49","TRNT1")

## Pull their abundance matrix #####
gene_col  <- rowData(prot_qf[["log_norm_prot"]])$PG.Genes
trna_rows <- gene_col %in% trna_genes
trna_ids  <- rownames(prot_qf[["log_norm_prot"]])[trna_rows]

quant_trna <- assay(prot_qf[["log_norm_prot"]])[trna_ids, ]
rownames(quant_trna) <- gene_col[trna_rows]

## Column annotation (bait + stress) 
col_annot <- data.frame(
  Bait   = colData(prot_qf)$group,
  Stress = ifelse(colData(prot_qf)$stress_short == "str", "Stressed", "Unstressed"),
  row.names = colnames(prot_qf[["log_norm_prot"]])
)

## Fixed column order: Unstressed -> Stressed; Beads -> AMP -> K57A 
col_order <- order(
  factor(col_annot$Stress, levels = c("Unstressed", "Stressed")),
  factor(col_annot$Bait,   levels = c("Bead control", "AMP-Rlig1", "K57A-Rlig1"))
)
quant_trna <- quant_trna[, col_order]          # reorder matrix columns
col_annot  <- col_annot[col_order, , drop = FALSE]  # reorder annotation

## Fix legend order via factor levels 
col_annot$Bait   <- factor(col_annot$Bait,
                           levels = c("Bead control", "AMP-Rlig1", "K57A-Rlig1"))
col_annot$Stress <- factor(col_annot$Stress,
                           levels = c("Unstressed", "Stressed"))


## Annotation bar colours 
annot_colors <- list(
  Bait   = c("Bead control" = "#1B9E77",
             "AMP-Rlig1"    = "#D14984",
             "K57A-Rlig1"   = "#7570B3"),
  Stress = c("Unstressed"   = "#4F797E",
             "Stressed"     = "#D19649")
)

## Row clustering (columns stay in fixed order) 
dist_rows_trna <- as.dist(1 - cor(t(quant_trna), use = "pairwise.complete.obs"))

## Draw + save 
png("~/Desktop/trna_heatmap.png", width = 11, height = 9, units = "in", res = 300)
pheatmap(
  quant_trna,
  scale = "row",                       # Z-score per enzyme
  cluster_cols = FALSE,                # keep bait x stress order
  cluster_rows = TRUE,                 # group similar enzymes
  clustering_distance_rows = dist_rows_trna,
  clustering_method = "ward.D2",
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = col_annot,
  annotation_colors = annot_colors,
  gaps_col = 9,                        # gap between Unstressed | Stressed
  border_color = NA,
  color = colorRampPalette(rev(brewer.pal(7, "RdBu")))(100),
  fontsize = 15,
  fontsize_row = 14,                   # readable — only ~21 rows
  fontsize_col = 13,
  annotation_legend = TRUE,
  main = "tRNA-repair machinery across baits and stress states"
)
dev.off()


##### Doors for Future work using Upset graphs#####

### See which proteins are in which contrast (if a protein is present in 6 it's baiscally noise)
prot_graph <- read.delim("transposed.tsv",sep="\t")
upset(prot_graph, sets = c("AMP_vs_Beads_str","AMP_vs_Beads_unst",
                           "AMP_vs_K57A_str", "AMP_vs_K57A_unst",
                           "K57A_vs_Beads_str","K57A_vs_Beads_unst"))

### See the same thing but without the noise proteins 

png("~/Desktop/upset_plot.png", width = 14, height = 9, units = "in", res = 300)

prot_graph_clean <- read.delim("clean.tsv",sep="\t")
upset(prot_graph_clean,
      sets = c("AMP_vs_Beads_str", "AMP_vs_Beads_unst",
               "AMP_vs_K57A_str",  "AMP_vs_K57A_unst",
               "K57A_vs_Beads_str","K57A_vs_Beads_unst"),
      text.scale = 2.2,               # bump ALL text up 
      point.size = 3,                 # bigger dots
      line.size = 1.2,                # thicker connecting lines
      mainbar.y.label = "Intersection Size",
      sets.x.label = "Set Size")

# Add a title on top 
grid.text("Shared enriched proteins across contrasts",
          x = 0.65, y = 0.97,
          gp = gpar(fontsize = 20, fontface = "bold"))

dev.off()   







