library(readxl)
library(tidyverse)
library(ggrepel)
library(pheatmap)
library(data.table)
library(dplyr)

network <- read.delim(
  "string_interactions_short_15_nodes.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)
head(network)

nodes <- unique(c(network$'#node1', network$node2))
length(nodes)
nodes

mapping <- read_excel(
  "protein_STRING_Mapping.xlsx")
head(mapping)
colnames(mapping)
mapping
setdiff(nodes, mapping$STRING_ID)

edges <- network %>%
  select(
    node1 = `#node1`,
    node2
  )

edges <- edges %>%
  rowwise() %>%
  mutate(
    protein_A = min(node1, node2),
    protein_B = max(node1, node2)
  ) %>%
  ungroup() %>%
  distinct(protein_A, protein_B)

degree_table <- bind_rows(
  edges %>% count(protein_A, name = "degree") %>%
    rename(STRING_ID = protein_A),
  
  edges %>% count(protein_B, name = "degree") %>%
    rename(STRING_ID = protein_B)
) %>%
  group_by(STRING_ID) %>%
  summarise(
    degree = sum(degree),
    .groups = "drop"
  )

gene_network <- mapping %>%
  left_join(
    degree_table,
    by = "STRING_ID"
  )
gene_network

enrichment <- read.delim(
  "enrichment.all.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

colnames(enrichment) #check
enrichment <- enrichment %>%
  rename(
    category = `#category`,
    term_id = `term ID`,
    description = `term description`,
    gene_count = `observed gene count`,
    background_count = `background gene count`,
    strength = strength,
    FDR = `false discovery rate`,
    matching_ids = `matching proteins in your network (IDs)`,
    matching_genes = `matching proteins in your network (labels)`
  )

total_genes <- 14

GO_BP <- enrichment %>%
  filter(
    category == "GO Process",
    FDR < 0.05
  ) %>%
  mutate(
    GeneRatio = gene_count / total_genes,
    neg_log10_FDR = -log10(FDR)
  )
nrow(GO_BP)

GO_MF <- enrichment %>%
  filter(
    category == "GO Function",
    FDR < 0.05
  )
nrow(GO_MF)
GO_CC <- enrichment %>%
  filter(
    category == "GO Component",
    FDR < 0.05
  )
nrow(GO_CC)

GO_BP_plot <- GO_BP %>%
  arrange(FDR) %>%
  mutate(
    description = factor(
      description,
      levels = description
    )
  )

ggplot(
  GO_BP_plot,
  aes(
    x = GeneRatio,
    y = description,
    size = gene_count,
    color = neg_log10_FDR
  )
) +
  geom_point(alpha = 0.9) +
  scale_size_continuous(
    name = "Gene Count",
    range = c(4, 10)
  ) +
  scale_color_gradient(
    name = "-log10(FDR)",
    low = "steelblue4",
    high = "skyblue"
  ) +
  labs(
    title = "GO Biological Process Enrichment",
    x = "Gene Ratio",
    y = "Biological Process"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    axis.text.y = element_text(
      size = 10
    ),
    axis.text.x = element_text(
      size = 10
    ),
    axis.title = element_text(
      face = "bold"
    ),
    legend.title = element_text(
      face = "bold"
    )
  )

GO_MF_plot <- GO_MF %>%
  arrange(FDR) %>%
  slice_head(n = 15) %>%
  mutate(
    term = reorder(description, strength),
    neg_log10_FDR = -log10(FDR)
  )
ggplot(
  GO_MF_plot,
  aes(
    x = strength,
    y = term,
    size = gene_count,
    color = neg_log10_FDR
  )
) +
  geom_point() +
  labs(
    title = "GO Molecular Function Enrichment",
    x = "Enrichment strength",
    y = "Molecular Function",
    size = "Gene count",
    color = "-log10(FDR)"
  ) +
  theme_bw()

GO_CC_plot <- GO_CC %>%
  arrange(FDR) %>%
  slice_head(n = 15) %>%
  mutate(
    term = reorder(description, strength),
    neg_log10_FDR = -log10(FDR)
  )

ggplot(
  GO_CC_plot,
  aes(
    x = strength,
    y = term,
    size = gene_count,
    color = neg_log10_FDR
  )
) +
  geom_point() +
  labs(
    title = "GO Cellular Component Enrichment",
    x = "Enrichment strength",
    y = "Cellular Component",
    size = "Gene count",
    color = "-log10(FDR)"
  ) +
  theme_bw()
#heatmap generation
annotation <- read.delim(
  "string_functional_annotations_15_nodes.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
annotation <- annotation %>%
  rename(
    term_id = `term ID`,
    description = `term description`
  )
colnames(annotation)

significant_BP <- enrichment %>%
  filter(
    category == "GO Process",
    FDR < 0.05
  ) %>%
  select(
    `term_id`,
    `description`
  )
colnames(significant_BP)

BP_annotations <- annotation %>%
  filter(category == "GO Process") %>%
  semi_join(
    significant_BP,
    by = c("term_id", "description")
  )
BP_annotations <- BP_annotations %>%
  left_join(
    mapping %>% select(STRING_ID, 'Protein name'),
    by = c("#node" = "STRING_ID")
  )
#heatmap generation
heatmap_data <- BP_annotations %>%
  select(
    'Protein name',
    description
  ) %>%
  distinct() %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = description,
    values_from = value,
    values_fill = 0
  )
heatmap_matrix <- as.data.frame(heatmap_data)

rownames(heatmap_matrix) <- heatmap_matrix$'Protein name'

heatmap_matrix$'Protein name' <- NULL

heatmap_matrix <- as.matrix(heatmap_matrix)
#heatmap
pheatmap(
  heatmap_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "GO Biological Process – Network Protein Associations",
  border_color = NA
)
