library(phyloseq)
library(ggplot2)
library(vegan)
library(dplyr)
library(DESeq2)
#Load the BIOM file
pseqbr <- import_biom("salmonella_sample.biom")
#Rename the columns
colnames(tax_table(pseqbr)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
sample_names(pseqbr) <- gsub("_bracken_species", "", sample_names(pseqbr))
#Load the Metadata
meta <- read.csv("metadata.csv", row.names = 1) # Sets the SRR column as the ID
sample_data(pseqbr) <- sample_data(meta)
sample_sums(pseqbr)
#Rarefaction
otu_table <- as.data.frame(t(otu_table(pseqbr)))
rare_curve <- rarecurve(otu_table, step = 100000)
#relative abundance
physeq_rel <- transform_sample_counts(pseqbr, function(x) x / sum(x))

#A plot at phylum level for cleaner visualization
physeq_phy <- tax_glom(physeq_rel, taxrank = "Phylum")
dfphy <- psmelt(physeq_phy)
ggplot(dfphy, aes(x = Abbreviation, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(y = "Relative Abundance", x = "Sample")
#Getting a better visualization
#Taking top 20 phyla
top20 <- names(sort(taxa_sums(physeq_rel), decreasing = TRUE)[1:20])
pseqtop20  <- prune_taxa(top20, physeq_rel)
df_top20 <- psmelt(pseqtop20)
#Plot
ggplot(df_top20, aes(x = Abbreviation, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(y = "Relative Abundance", x = "Sample")

#Checking out at the family level
physeq_fam <- tax_glom(physeq_rel, taxrank = "Family")
dffam <- psmelt(physeq_fam)
ggplot(dffam, aes(x = Abbreviation, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(y = "Relative Abundance", x = "Sample")
#taking top 20 family
top20fam <- names(sort(taxa_sums(physeq_rel), decreasing = TRUE)[1:20])
pseqtop20fam  <- prune_taxa(top20fam, physeq_rel)
df_top20fam <- psmelt(pseqtop20fam)
#Plot
ggplot(df_top20fam, aes(x = Abbreviation, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(y = "Relative Abundance", x = "Sample")
#Checking out at the genus level
physeq_gen <- tax_glom(physeq_rel, taxrank = "Genus")
dfgen <- psmelt(physeq_gen)
ggplot(dfgen, aes(x = Abbreviation, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(title = "Genus-Level Relative Abundance by Group", y = "Relative Abundance", x = "Sample")
#taking top 20 genus
top10gen <- names(sort(taxa_sums(physeq_rel), decreasing = TRUE)[1:10])
pseqtop10gen  <- prune_taxa(top10gen, physeq_rel)
df_top10gen <- psmelt(pseqtop10gen)
#Plot
ggplot(df_top10gen, aes(x = Abbreviation, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(y = "Relative Abundance", x = "Sample")
#Checking out at the species level
physeq_spe <- tax_glom(physeq_rel, taxrank = "Species")
dfspe <- psmelt(physeq_spe)
ggplot(dfspe, aes(x = Abbreviation, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(y = "Relative Abundance", x = "Sample")
#taking top 10 species
top10spe <- names(sort(taxa_sums(physeq_rel), decreasing = TRUE)[1:10])
pseqtop10spe  <- prune_taxa(top10spe, physeq_rel)
df_top10spe <- psmelt(pseqtop10spe)
#Plot
ggplot(df_top10spe, aes(x = Abbreviation, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(y = "Relative Abundance", x = "Sample")
#Alpha Diversity
plot_richness(pseqbr, x="Diet", color="Abbreviation", measures=c("Observed", "Shannon", "Chao1", "Ace", "Simpson", "InvSimpson", "Fisher")) +
  geom_boxplot(aes(group=samples), alpha=0.1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  labs(title = "Alpha Diversity",
       x = "Samples",
       y = "Alpha Diversity Measure")

#Beta Diversity
# Here's a PCoA with bray-curtis
  ord.pcoa.bray <- ordinate(pseqbr, method="PCoA", distance="bray")
  plot_ordination(pseqbr, ord.pcoa.bray, color="Abbreviation", title="Bray PCoA: Individual Samples") + 
    geom_point(size = 3) +
    theme_minimal()
  
# Here's an NMDS with the same distance measure. 
  ord.nmds.bray <- ordinate(pseqbr, method="NMDS", distance="bray")
  plot_ordination(pseqbr, ord.nmds.bray, color="Abbreviation", title="Bray NMDS") + geom_point(size = 4)
  
# What about Jaccard distance?
  ord.pcoa.jaccard <- ordinate(pseqbr, method="PCoA", distance="jaccard")
  plot_ordination(pseqbr, ord.pcoa.jaccard, color="Abbreviation", title="Bray PCoA") + geom_point(size = 4)

#Permanova
  meta <- as(sample_data(pseqbr), "data.frame")
  permanova <- adonis2(phyloseq::distance(pseqbr, method = "bray") ~ Diet, data = meta)
 print(permanova)

#Differential Abundance Analysis
dds <- phyloseq_to_deseq2(pseqbr, ~ Diet) #Converting phyloseq to deseq format
dds_res <- DESeq(dds, test="Wald", fitType="parametric")
res <- results(dds_res, cooksCutoff = FALSE)
res_df <- as.data.frame(res)
sig_results <- res_df[which(res_df$padj < 0.05), ] #Filter for padj<0.05
tax_table_df <- as.data.frame(tax_table(pseqbr)) #Merging the tax table to the results
sig_plus_tax <- merge(sig_results, tax_table_df, by="row.names")
sig_plus_tax <- sig_plus_tax[order(sig_plus_tax$log2FoldChange, decreasing = TRUE), ] # Rearranging by Log2FoldChange (highest difference first)
head(sig_plus_tax) 
 
#Volcano plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
   geom_point(aes(color = padj < 0.05), alpha = 0.5) +
   scale_color_manual(values = c("black", "red")) +
   theme_bw() +
   labs(title = "Volcano Plot: Vegan vs. Omnivore",
        x = "Log2 Fold Change",
        y = "-Log10 Adjusted P-value")
 