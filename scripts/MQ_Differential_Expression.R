# set to the current working directory
setwd("~/Desktop/swopnil.pradhan/MQ_RNA_seq/counts")

# install the required libraries 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("BiocManager") 
# Install only once if not already installed
BiocManager::install(c("EnhancedVolcano", "DESeq2"), force = T)
library(EnhancedVolcano)
library(DESeq2)
library(dplyr)
library(pheatmap)
library(ggplot2)

# import the counts matrix
feature_counts <- read.table(file="counts.txt", header=TRUE)

# remove the couluns 2:6 as it contains info which we do not need currently (Chr, Start, End, Strand and Length)
# We are only intrested in the Geneid and the Counts for each of the samples
feature_counts <- feature_counts[, -c(2:6)]
head(feature_counts)

# prepare metadata file in excel and import 
metadata <- read.csv("metadata.csv", header = T, sep = ",")

# set each of the colnames as factor to include in the design
metadata$Strain <- factor(metadata$Strain)
metadata$Group <- factor(metadata$Group)
metadata$Muscle <- factor(metadata$Muscle)

# set the Geneid as rownames to ensure the matrix of feature_counts and metadata matches 12*12
rownames(feature_counts) <- feature_counts$Geneid
feature_counts <- feature_counts[ , -1]

# check result
head(feature_counts)

# now we start the differential expression
# We have also set group as one of the factors

dds <- DESeqDataSetFromMatrix(countData = feature_counts,
                              colData = metadata,
                              design = ~Group)


dds <- DESeq(dds)
res <- results(dds)

resultsNames(dds)

# we can do pairwise comparisons to compare each of the groups
groups <- levels(dds$Group)
pairs <- combn(groups, 2, simplify = F)

# results_list is an empty list to store results of each pairwise comparisons. 
results_list <- list()

# Loops through each pair of groups and stores results in results_list
# We have 6 possible comparisons for our analysis using the 4 groups.

for (p in pairs) {
  res <- results(dds, contrast = c("Group", p[1], p[2]))
  comparison_name <- paste(p[1], "_vs_", p[2], sep = "")
  results_list[[comparison_name]] <- res
}

# Now we can identify the diffrentially expressed genes which are also significant for each of the study. 
# Setting up the thereshold levels. 
padj_cutoff <- 0.05
log2fc_cutoff <- 1

# Again setting DEGs_list as an empty list to store the data
DEGs_list <- list()

for (name in names(results_list)) {
  res <- results_list[[name]]
  
  # Filter significant genes
  DEGs <- subset(res, padj < padj_cutoff & abs(log2FoldChange) >= log2fc_cutoff)
  
  # Store in list
  DEGs_list[[name]] <- DEGs
}


DEG_genes <- lapply(DEGs_list, function(x) rownames(x))


# Saving the DEGs for each comparisons as .csv file
for (name in names(DEGs_list)) {
  write.csv(as.data.frame(DEGs_list[[name]]), paste0(name, "_DEGs.csv"))
}

# Combine all DEGs into a single vector
all_DEGs <- unique(unlist(lapply(DEGs_list, rownames)))

# Get normalized counts
norm_counts <- counts(dds, normalized = TRUE)
DEG_counts <- norm_counts[all_DEGs, ]

# Optionally, scale by gene for better heatmap visualization
scaled_counts <- t(scale(t(DEG_counts)))

# Create heatmap
library(grid)

heat_map <- grid.grabExpr(pheatmap(scaled_counts,
                            cluster_rows = TRUE,
                            cluster_cols = TRUE,
                            show_rownames = FALSE,
                            show_colnames = TRUE,
                            main = "Heatmap of DEGs across all comparisons"))


ggsave("DEG_heatmap.pdf", heat_map, width = 8, height = 10)

install.packages("UpSetR")
library(UpSetR)

# Create a list of gene sets per comparison
DEG_gene_sets <- lapply(DEGs_list, rownames)

# Convert to UpSetR format
DEG_upset <- fromList(DEG_gene_sets)

# Plot
upset_plot <- upset(DEG_upset, 
      order.by = "freq",
      main.bar.color = "steelblue",
      sets.bar.color = "tomato",
      text.scale = c(2, 2, 1.5, 1.5, 1.5, 1.2))

upset_plot

tiff("Upset_plot.tiff", width = 3000, height = 2400, res = 300)
upset(DEG_upset, 
      order.by = "freq",
      main.bar.color = "steelblue",
      sets.bar.color = "tomato",
      text.scale = c(2, 2, 1.5, 1.5, 1.5, 1.2))
dev.off()


# Volcano Plots
# Create a folder if it doesn't exist
if(!dir.exists("Volcano_Plots")){
  dir.create("Volcano_Plots")
}

plot_volcano <- function(res, title = "Volcano Plot", log2fc_cutoff = 1, padj_cutoff = 0.05) {
  
  res_df <- as.data.frame(res)
  res_df$Gene <- rownames(res_df)
  
  # Assign significance
  res_df$Significant <- "Not Significant"
  res_df$Significant[res_df$padj < padj_cutoff & res_df$log2FoldChange >= log2fc_cutoff] <- "Up"
  res_df$Significant[res_df$padj < padj_cutoff & res_df$log2FoldChange <= -log2fc_cutoff] <- "Down"
  
  # Count Up/Down genes
  n_up <- sum(res_df$Significant == "Up")
  n_down <- sum(res_df$Significant == "Down")
  
  # Create volcano plot
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Down" = "blue", "Up" = "red", "Not Significant" = "gray")) +
    geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "black") +
    annotate("text", x = max(res_df$log2FoldChange, na.rm = TRUE)*0.7, 
             y = max(-log10(res_df$padj), na.rm = TRUE)*0.9, 
             label = paste0("Up: ", n_up), color = "red", size = 5) +
    annotate("text", x = min(res_df$log2FoldChange, na.rm = TRUE)*0.7, 
             y = max(-log10(res_df$padj), na.rm = TRUE)*0.9, 
             label = paste0("Down: ", n_down), color = "blue", size = 5) +
    theme_minimal() +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10 Adjusted P-value")
  
  return(p)
}

for (name in names(results_list)) {
  res <- results_list[[name]]
  
  # Generate annotated volcano plot
  volcano_plot <- plot_volcano(res, title = paste0("Volcano: ", name))
  
  # Define file name
  file_name <- paste0("Volcano_Plots/", name, "_volcano.pdf")
  
  # Save PDF
  ggsave(file_name, plot = volcano_plot, width = 7, height = 5)
}

library(clusterProfiler)
library(org.Gg.eg.db)  # change to your organism (e.g., org.Gg.eg.db for chicken)
library(enrichplot)
library(dplyr)


# Create folder for GO results
if(!dir.exists("GO_Results")){
  dir.create("GO_Results")
}

# Function for GO analysis with grouped BP/MF/CC barplot
run_GO_combined_barplot_grouped <- function(res, name, padj_cutoff = 0.05, log2fc_cutoff = 1) {
  
  res_df <- as.data.frame(res)
  res_df$Gene <- rownames(res_df)
  
  # Filter significant genes
  sig_genes <- res_df$Gene[res_df$padj < padj_cutoff & abs(res_df$log2FoldChange) >= log2fc_cutoff]
  if (length(sig_genes) == 0) {
    message(paste("⚠️ No significant genes for:", name))
    return(NULL)
  }
  
  # Map gene symbols to ENTREZ IDs (chicken)
  entrez_ids <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Gg.eg.db)
  if (is.null(entrez_ids) || nrow(entrez_ids) == 0) {
    message(paste("⚠️ No gene mappings for:", name))
    return(NULL)
  }
  
  # Run GO enrichment for BP, MF, CC
  go_terms <- list()
  for (ont in c("BP", "MF", "CC")) {
    ego <- enrichGO(
      gene = entrez_ids$ENTREZID,
      OrgDb = org.Gg.eg.db,
      keyType = "ENTREZID",
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      readable = TRUE
    )
    
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      df <- as.data.frame(ego)
      df$Ontology <- ont
      go_terms[[ont]] <- df
    }
  }
  
  if (length(go_terms) == 0) {
    message(paste("⚠️ No enriched GO terms for:", name))
    return(NULL)
  }
  
  # Combine BP, MF, CC results
  go_combined <- dplyr::bind_rows(go_terms)
  
  # Save combined table
  write.csv(go_combined, paste0("GO_Results/", name, "_GO_combined.csv"), row.names = FALSE)
  
  # Pick top 10 terms from each ontology
  go_plot <- go_combined %>%
    group_by(Ontology) %>%
    top_n(-10, p.adjust) %>%
    arrange(Ontology, p.adjust) %>%   # group by Ontology, sort by p.adjust
    ungroup()
  
  # Create factor ordering: BP terms first, then MF, then CC
  go_plot$Description <- factor(go_plot$Description, levels = rev(unique(go_plot$Description)))
  
  # --- GROUPED BARPLOT ---
  p <- ggplot(go_plot, aes(x = Description, y = Count, fill = Ontology)) +
    geom_bar(stat = "identity") +
    facet_wrap(~Ontology, scales = "free_y", ncol = 1, strip.position = "top") +  # separate panels stacked vertically
    coord_flip() +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("GO Enrichment by Ontology (BP, MF, CC):", name),
      x = "GO Term",
      y = "Gene Count",
      fill = "Ontology"
    ) +
    theme(
      axis.text.y = element_text(size = 9),
      strip.text = element_text(face = "bold", size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # Save plot
  ggsave(paste0("GO_Results/", name, "_GO_combined_grouped_barplot.pdf"), plot = p, width = 9, height = 10)
  
  return(go_combined)
}

# 🔁 Run for all DE result sets
GO_combined_results <- list()

for (name in names(results_list)) {
  res <- results_list[[name]]
  GO_combined_results[[name]] <- run_GO_combined_barplot_grouped(res, name)
}


# KEGG pathway enrichment analysis

# Create output folder
if(!dir.exists("KEGG_Results")){
  dir.create("KEGG_Results", recursive = TRUE)
}

# Function for KEGG enrichment + Rich Factor bubble plot
run_KEGG_richfactor <- function(res, name, padj_cutoff = 0.05, log2fc_cutoff = 1) {
  
  res_df <- as.data.frame(res)
  res_df$Gene <- rownames(res_df)
  
  # Select significant genes
  sig_genes <- res_df$Gene[res_df$padj < padj_cutoff & abs(res_df$log2FoldChange) >= log2fc_cutoff]
  if (length(sig_genes) == 0) {
    message(paste("⚠️ No significant genes for:", name))
    return(NULL)
  }
  
  # Convert gene symbols to ENTREZ IDs
  entrez_ids <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Gg.eg.db)
  if (is.null(entrez_ids) || nrow(entrez_ids) == 0) {
    message(paste("⚠️ No mapped genes for:", name))
    return(NULL)
  }
  
  # Run KEGG enrichment for chicken (gga)
  ekegg <- enrichKEGG(
    gene = entrez_ids$ENTREZID,
    organism = "gga",
    keyType = "kegg",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  if (is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0) {
    message(paste("⚠️ No enriched KEGG pathways for:", name))
    return(NULL)
  }
  
  # Make results readable
  ekegg <- setReadable(ekegg, OrgDb = org.Gg.eg.db, keyType = "ENTREZID")
  kegg_df <- as.data.frame(ekegg)
  
  # Calculate Rich Factor
  # Example GeneRatio format: "5/89" → 5 DE genes out of 89 annotated
  kegg_df <- kegg_df %>%
    mutate(
      GeneRatioNum = as.numeric(sub("/.*", "", GeneRatio)),
      BgRatioNum = as.numeric(sub(".*/", "", BgRatio)),
      RichFactor = GeneRatioNum / BgRatioNum
    )
  
  # Save full KEGG results table with RichFactor
  out_csv <- file.path("KEGG_Results", paste0(name, "_KEGG_richfactor.csv"))
  write.csv(kegg_df, out_csv, row.names = FALSE)
  
  # Select top pathways (by adjusted p-value)
  kegg_plot <- kegg_df %>%
    arrange(p.adjust) %>%
    head(20)
  
  # Bubble plot using Rich Factor
  p <- ggplot(kegg_plot, aes(x = RichFactor,
                             y = reorder(Description, RichFactor),
                             size = Count,
                             color = -log10(p.adjust))) +
    geom_point(alpha = 0.8) +
    scale_color_gradient(low = "blue", high = "red") +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("KEGG Pathway Enrichment (Rich Factor):", name),
      x = "Rich Factor",
      y = "Pathway",
      color = "-log10(adj. p-value)",
      size = "Gene Count"
    ) +
    theme(
      axis.text.y = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # Save bubble plot
  out_pdf <- file.path("KEGG_Results", paste0(name, "_KEGG_richfactor_bubbleplot.pdf"))
  ggsave(out_pdf, plot = p, width = 9, height = 7)
  return(kegg_df)
}

# 🔁 Run KEGG analysis for all datasets in results_list
KEGG_results_list <- list()

for (name in names(results_list)) {
  res <- results_list[[name]]
  KEGG_results_list[[name]] <- run_KEGG_richfactor(res, name)
}
