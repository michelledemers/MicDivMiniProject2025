if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")
library(dada2)
packageVersion("dada2")

setwd("fastqFiles/")

fnFs <- sort(list.files(getwd(), pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(getwd(), pattern="_R2_001.fastq", full.names = TRUE))
dada2::plotQualityProfile(fnFs[1:2])
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#filter step
filtFs <- file.path(getwd(), "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(getwd(), "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
#rerun this, but remove the file "25-MAD-B1S-NO3-T1_S58_L001_R1_001.fastq"
sample.names2 <- sample.names[-grep("25-MAD-B1S-NO3-T1", sample.names)]
filtFs <- file.path(getwd(), "filtered", paste0(sample.names2, "_F_filt.fastq.gz"))
filtRs <- file.path(getwd(), "filtered", paste0(sample.names2, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names2
names(filtRs) <- sample.names2
out_filt <- dada2::filterAndTrim(fnFs[-grep("25-MAD-B1S-NO3-T1", sample.names)], filtFs, fnRs[-grep("25-MAD-B1S-NO3-T1", sample.names)], filtRs, truncLen=c(150,150),
                            maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                            compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out_filt)

errF_filt <- dada2::learnErrors(filtFs, multithread=TRUE)
errR_filter <- dada2::learnErrors(filtRs, multithread=TRUE)
dadaFs_filt <- dada2::dada(filtFs, err=errF_filt, multithread=TRUE)
dadaRs_filt <- dada2::dada(filtRs, err=errR_filter, multithread=TRUE)
mergers_filt <- dada2::mergePairs(dadaFs_filt, filtFs, dadaRs_filt, filtRs, justConcatenate = TRUE)
head(mergers_filt[2])
seqtab_filt <- dada2::makeSequenceTable(mergers_filt)
dim(seqtab_filt)

errF <- dada2::learnErrors(fnFs, multithread=TRUE)
errR <- dada2::learnErrors(fnRs, multithread=TRUE)
dadaFs <- dada2::dada(fnFs, err=errF, multithread=TRUE)
dadaRs <- dada2::dada(fnRs, err=errR, multithread=TRUE)
mergers <- dada2::mergePairs(dadaFs, fnFs, dadaRs, fnRs, justConcatenate = TRUE)
head(mergers[2])
seqtab <- dada2::makeSequenceTable(mergers)
dim(seqtab)

## remove chimeras
seqtab_filt.nochim <- dada2::removeBimeraDenovo(seqtab_filt, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab_filt.nochim)
sum(seqtab_filt.nochim)/sum(seqtab_filt)
seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

##track reads
getN <- function(x) sum(dada2::getUniques(x))
track_filt <- cbind(out_filt, sapply(dadaFs_filt, getN), sapply(dadaRs_filt, getN), sapply(mergers_filt, getN), rowSums(seqtab_filt.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track_filt) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_filt) <- sample.names2
head(track_filt)
out[,2] <- out[,1]
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

## taxonomy
taxa_filt <- dada2::assignTaxonomy(seqtab_filt.nochim, "../silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa_filt.print <- taxa_filt # Removing sequence rownames for display only
rownames(taxa_filt.print) <- NULL
head(taxa_filt.print)

taxa <- dada2::assignTaxonomy(seqtab.nochim, "../silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

library(phyloseq)
library(ggplot2)
library(Biostrings)
library(dplyr)

samdf <- read.table("16smetadata_phyloseq.txt", header = 1)
samdf_filt <- samdf[-grep("25-MAD-B1S-NO3-T1", samdf$sample.id),]
row.names(samdf) <- NULL
samdf <- samdf %>%
  tibble::column_to_rownames("sample.id") 
row.names(samdf_filt) <- NULL
samdf_filt <- samdf_filt %>%
  tibble::column_to_rownames("sample.id") 
ps_filt <- phyloseq(otu_table(seqtab_filt.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf_filt), 
               tax_table(taxa_filt))
row.names(seqtab.nochim) <- sapply(strsplit(basename(row.names(seqtab.nochim)), "_"), `[`, 1)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

## import from qiime2 pipeline

library(decontam)
library(qiime2R)
library(tidyr)
asvs <- read_qza("03-dada2-feature-table_forward.qza")
taxonomy <- read_qza("05-taxonomy-blast-90-1.qza")
tree <- read_qza("07-rooted-tree.qza")
asv_df <- asvs$data # we can view and save the actual data by specifying '$'
taxonomy_df <- taxonomy$data # save taxonomy info as a dataframe
phylo_tree <- tree$data

taxonomy_fixed_df <- taxonomy_df %>% separate_wider_delim(Taxon, delim = ";", names_sep = "", too_few = "align_start")
taxonomy_fixed_df[is.na(taxonomy_fixed_df)] <- "Unassigned" # rename NAs into unassigned
head(taxonomy_fixed_df)

taxonomy_fixed_df <- as.data.frame(taxonomy_fixed_df) # force into a dataframe
row.names(taxonomy_fixed_df) <- taxonomy_fixed_df$Feature.ID # make first column into row names
taxonomy_fixed_df$Feature.ID <- NULL # remove the first column
taxonomy_matrix <- as.matrix(taxonomy_fixed_df) # convert to a matrix 

physeq_asv <- otu_table(asv_df, taxa_are_rows = T) # convert into phyloseq object
physeq_tax <- tax_table(taxonomy_matrix) # convert into phyloseq object
physeq_meta <- sample_data(samdf) # convert into phyloseq object

ps_forward <- phyloseq(physeq_asv, physeq_tax, physeq_meta) # merge into phyloseq object
#ps_forward_tree <- merge_phyloseq(ps_forward, phylo_tree) # add tree into phyloseq object
# the tree command isn't really working so we'll move on
ps_forward # work with this one going forward
ps_filt
ps

# look at size of libraries in controls versus true samples
df <- as.data.frame(sample_data(ps_forward)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps_forward)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

sample_data(ps_forward)$is.neg <- sample_data(ps_forward)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(ps_forward, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)


ps_forward_sans_controls <- subset_samples(ps_forward, Sample_or_Control != "Control_sample")

## alpha and beta diversity
ps_forward_ra <- transform_sample_counts(ps_forward_sans_controls, function(x) x / sum(x)) # convert to relative abundance
p <- plot_bar(ps_forward_ra, "Sample", fill="Taxon2") + # Sample is x-axis
  geom_bar(aes(color=Taxon2), stat="identity") + # specifies the type of graph and color
  facet_grid(~Timepoint, scales = "free", space = "free_x") # lets separate our samples by Site
p
ggsave(plot = p, filename = "phyla_taxonomy_plot.pdf", width = 20, height = 10, units = "in")

p <- plot_bar(ps_forward_ra, "Sample", fill="Taxon3") + # Sample is x-axis
  geom_bar(aes(color=Taxon3), stat="identity") + # specifies the type of graph and color
  facet_grid(~Timepoint, scales = "free", space = "free_x") # lets separate our samples by Site
p
ggsave(plot = p, filename = "class_taxonomy_plot.pdf", width = 20, height = 10, units = "in")

#filter to top 15 class and phyla
ps_phyla <- tax_glom(ps_forward_ra, taxrank = "Taxon2", NArm = FALSE)
rel_abund_df <- psmelt(ps_phyla)
names(rel_abund_df)[13] <- "Domain"
names(rel_abund_df)[14] <- "Phylum"
taxaSums <- data.frame(tax_table(ps_phyla)[,"Taxon2"],
                      taxa_sums = taxa_sums(ps_phyla)) %>%
  arrange(desc(taxa_sums))
top15 = head(rownames(taxaSums), 15)
y = prune_taxa(top15, ps_phyla)
df1 = data.frame(ID = c(taxa_names(y), "Other"), Phylum = c(tax_table(y)[,"Taxon2"], "Other"))
df2 = t(cbind(otu_table(t(y)), data.frame(Other = 1 - sample_sums(y))))
df = cbind(df1, df2) %>%
  pivot_longer(-c(ID, Phylum), names_to = "SampleType", values_to = "Abundance") %>%
  as.data.frame
df$Phylum <- gsub("p__", "", df$Phylum)
df$Phylum = as.factor(df$Phylum)
df <- merge(df, samdf, by.x = "SampleType", by.y = "row.names")
df$Phylum = factor(df$Phylum, levels = c(gsub("p__", "",taxaSums[top15,"Taxon2"]), "Other"))

p <- ggplot(df[-grep("SVS",df$SampleType),], aes(SampleType, y = Abundance, fill = Phylum)) +
  geom_bar(aes(color=Phylum), stat="identity") +
  facet_grid(~Timepoint, scales = "free", space = "free_x") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(plot = p, filename = "phyla_taxonomy_plot_top15.pdf", width = 20, height = 10, units = "in")


ps_class <- tax_glom(ps_forward_ra, taxrank = "Taxon3", NArm = FALSE)
rel_abund_df <- psmelt(ps_class)
names(rel_abund_df)[13] <- "Domain"
names(rel_abund_df)[14] <- "Phylum"
names(rel_abund_df)[15] <- "Class"
taxaSums <- data.frame(tax_table(ps_class)[,"Taxon3"],
                       taxa_sums = taxa_sums(ps_class)) %>%
  arrange(desc(taxa_sums))
top15 = head(rownames(taxaSums), 15)
y = prune_taxa(top15, ps_class)
df1 = data.frame(ID = c(taxa_names(y), "Other"), Class = c(tax_table(y)[,"Taxon3"], "Other"))
df2 = t(cbind(otu_table(t(y)), data.frame(Other = 1 - sample_sums(y))))
df = cbind(df1, df2) %>%
  pivot_longer(-c(ID, Class), names_to = "SampleType", values_to = "Abundance") %>%
  as.data.frame
df$Class <- gsub("c__", "", df$Class)
df$Class = as.factor(df$Class)
df <- merge(df, samdf, by.x = "SampleType", by.y = "row.names")
df$Class = factor(df$Class, levels = c(gsub("c__", "",taxaSums[top15,"Taxon3"]), "Other"))

p <- ggplot(df[-grep("SVS",df$SampleType),], aes(SampleType, y = Abundance, fill = Class)) +
  geom_bar(aes(color=Class), stat="identity") +
  facet_grid(~Timepoint, scales = "free", space = "free_x") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(plot = p, filename = "class_taxonomy_plot_top15.pdf", width = 20, height = 10, units = "in")



alpha_div_observed <- phyloseq::estimate_richness(ps_forward_sans_controls, measures = "Observed") # calculate alpha diversity
head(alpha_div_observed) # view the first few lines of our R object
metadata_df <- phyloseq::sample_data(ps_forward_sans_controls) # access metadata from the phyloseq object
head(metadata_df) 
alpha_div_observed_metadata <- data.frame(
  phyloseq::sample_data(ps_forward_sans_controls), # get metadata
  "Reads" = phyloseq::sample_sums(ps_forward_sans_controls), # number of reads
  "Observed" = phyloseq::estimate_richness(ps_forward_sans_controls, measures = "Observed"), # count observed ASVs
  "Shannon" = phyloseq::estimate_richness(ps_forward_sans_controls, measures = "Shannon"), # Calculate Shannon Diversity
  "InvSimpson" = phyloseq::estimate_richness(ps_forward_sans_controls, measures = "InvSimpson")) # calculate InvSimpson
alpha_div_observed_plot <- ggplot(alpha_div_observed_metadata[-which(is.na(alpha_div_observed_metadata$Treatment)),], aes(x=Treatment, y=Observed)) + # what are we plotting
  geom_boxplot(aes(color = Treatment)) + # make it into a boxplot plot
  theme_minimal() +
  facet_grid(~Timepoint, scales = "free", space = "free_x")# and lets make sure its "pretty"  
alpha_div_observed_plot
alpha_div_shannon_plot <- ggplot(alpha_div_observed_metadata[-which(is.na(alpha_div_observed_metadata$Treatment)),], aes(x=Treatment, y=Shannon)) + # what are we plotting
  geom_boxplot(aes(color = Treatment)) + # make it into a boxplot plot
  theme_minimal() +
  facet_grid(~Timepoint, scales = "free", space = "free_x")# and lets make sure its "pretty"  
alpha_div_shannon_plot
ggsave(plot = alpha_div_shannon_plot, filename = "alpha_div_shannon_plot.pdf", width = 10, height = 5, units = "in")

alpha_div_observed_metadata$Evenness <- alpha_div_observed_metadata$Shannon/log(alpha_div_observed_metadata$Observed)
head(alpha_div_observed_metadata)
alpha_div_Evenness_plot <- ggplot(alpha_div_observed_metadata[-which(is.na(alpha_div_observed_metadata$Treatment)),], aes(x=Treatment, y=Evenness)) + # what are we plotting
  geom_boxplot(aes(color = Treatment)) + # make it into a boxplot plot
  theme_minimal() +
  facet_grid(~Timepoint, scales = "free", space = "free_x")# and lets make sure its "pretty"  
alpha_div_Evenness_plot
ggsave(plot = alpha_div_Evenness_plot, filename = "alpha_div_Evenness_plot.pdf", width = 10, height = 5, units = "in")

alpha_div_shannon_plot <- ggplot(alpha_div_observed_metadata[-which(is.na(alpha_div_observed_metadata$Treatment)),], aes(x=Treatment, y=Shannon)) + # what are we plotting
  geom_boxplot(aes(color = Treatment)) + # make it into a boxplot plot
  theme_minimal() +
  facet_grid(Site~Timepoint, scales = "free", space = "free_x")# and lets make sure its "pretty"  
alpha_div_shannon_plot
ggsave(plot = alpha_div_shannon_plot, filename = "alpha_div_shannon_plot_site.pdf", width = 10, height = 10, units = "in")

alpha_div_observed_metadata$Evenness <- alpha_div_observed_metadata$Shannon/log(alpha_div_observed_metadata$Observed)
head(alpha_div_observed_metadata)
alpha_div_Evenness_plot <- ggplot(alpha_div_observed_metadata[-which(is.na(alpha_div_observed_metadata$Treatment)),], aes(x=Treatment, y=Evenness)) + # what are we plotting
  geom_boxplot(aes(color = Treatment)) + # make it into a boxplot plot
  theme_minimal() +
  facet_grid(Site~Timepoint, scales = "free", space = "free_x")# and lets make sure its "pretty"  
alpha_div_Evenness_plot
ggsave(plot = alpha_div_Evenness_plot, filename = "alpha_div_Evenness_plot_site.pdf", width = 10, height = 10, units = "in")

alpha_div_shannon_plot <- ggplot(alpha_div_observed_metadata[-which(alpha_div_observed_metadata$Timepoint == "TE"),], aes(x=Treatment, y=Shannon)) + # what are we plotting
  geom_boxplot(aes(color = Treatment)) + # make it into a boxplot plot
  theme_minimal() +
  facet_grid(Site~Timepoint, scales = "free", space = "free_x")# and lets make sure its "pretty"  
alpha_div_shannon_plot
ggsave(plot = alpha_div_shannon_plot, filename = "alpha_div_shannon_plot_faceted.pdf", width = 10, height = 10, units = "in")

alpha_div_observed_metadata$Evenness <- alpha_div_observed_metadata$Shannon/log(alpha_div_observed_metadata$Observed)
head(alpha_div_observed_metadata)
alpha_div_Evenness_plot <- ggplot(alpha_div_observed_metadata[-which(alpha_div_observed_metadata$Timepoint == "TE"),], aes(x=Treatment, y=Evenness)) + # what are we plotting
  geom_boxplot(aes(color = Treatment)) + # make it into a boxplot plot
  theme_minimal() +
  facet_grid(Site~Timepoint, scales = "free", space = "free_x")# and lets make sure its "pretty"  
alpha_div_Evenness_plot
ggsave(plot = alpha_div_Evenness_plot, filename = "alpha_div_Evenness_plot_faceted.pdf", width = 10, height = 10, units = "in")


alpha_div_observed_metadata <- data.frame(
  phyloseq::sample_data(ps_forward_sans_controls), # get metadata
  "Reads" = phyloseq::sample_sums(ps_forward_sans_controls), # number of reads
  "Observed" = phyloseq::estimate_richness(ps_forward_sans_controls, measures = "Observed"), # count observed ASVs
  "Shannon" = phyloseq::estimate_richness(ps_forward_sans_controls, measures = "Shannon"), # Calculate Shannon Diversity
  "InvSimpson" = phyloseq::estimate_richness(ps_forward_sans_controls, measures = "InvSimpson")) # calculate InvSimpson
ttest <- t(sapply(alpha_div_observed_metadata, function(x) unlist(t.test(x~sample_data(ps_forward_sans_controls)$Treatment)[c("estimate","p.value","statistic","conf.int")])))
ttest
alpha_div_observed_plot <- ggplot(alpha_div_observed_metadata[-which(is.na(alpha_div_observed_metadata$Treatment)),], aes(x=Treatment, y=Observed)) + # what are we plotting
  geom_boxplot(aes(color = Treatment)) + # make it into a boxplot plot
  theme_minimal() +
  facet_grid(~Timepoint, scales = "free", space = "free_x")# and lets make sure its "pretty"  
alpha_div_observed_plot
alpha_div_shannon_plot <- ggplot(alpha_div_observed_metadata[-which(is.na(alpha_div_observed_metadata$Treatment)),], aes(x=Treatment, y=Shannon)) + # what are we plotting
  geom_boxplot(aes(color = Treatment)) + # make it into a boxplot plot
  theme_minimal() +
  facet_grid(Site~Timepoint, scales = "free", space = "free_x")# and lets make sure its "pretty"  
alpha_div_shannon_plot
ggsave(plot = alpha_div_shannon_plot, filename = "alpha_div_shannon_plot2.pdf", width = 10, height = 5, units = "in")

alpha_div_observed_metadata$Evenness <- alpha_div_observed_metadata$Shannon/log(alpha_div_observed_metadata$Observed)
head(alpha_div_observed_metadata)
alpha_div_Evenness_plot <- ggplot(alpha_div_observed_metadata[-which(is.na(alpha_div_observed_metadata$Treatment)),], aes(x=Treatment, y=Evenness)) + # what are we plotting
  geom_boxplot(aes(color = Treatment)) + # make it into a boxplot plot
  theme_minimal() +
  facet_grid(Site~Timepoint, scales = "free", space = "free_x")# and lets make sure its "pretty"  
alpha_div_Evenness_plot
ggsave(plot = alpha_div_Evenness_plot, filename = "alpha_div_Evenness_plot2.pdf", width = 10, height = 5, units = "in")

## beta diversity
# look at all my samples, over time, space and treatment
ps_forward_sans_controls_sans_SVS <- subset_samples(ps_forward_sans_controls, Site == "A" | Site == "B")

ps_forward_sans_controls_relab <- transform_sample_counts(ps_forward_sans_controls_sans_SVS, function(x) x / sum(x) )
bray_dist <- phyloseq::distance(ps_forward_sans_controls_relab, method="bray") # calculate bray-curtis metric
ordination <- ordinate(ps_forward_sans_controls_relab, method="PCoA", distance=bray_dist) # Perform ordination using bray-curtis metric
p <-plot_ordination(ps_forward_sans_controls_relab, ordination, color="Treatment", shape = "Timepoint") + theme(aspect.ratio=1) +
  facet_wrap(~Site, nrow = 2)# and lets make sure its "pretty"  
# plot the ordination using ggplot2
ggsave(plot = p, filename = "PCoA_2.pdf", width = 4, height = 5, units = "in")


library(vegan)
adonis2(bray_dist ~ sample_data(ps_forward_sans_controls_sans_SVS)$Site)
adonis2(bray_dist ~ sample_data(ps_forward_sans_controls_sans_SVS)$Timepoint)
adonis2(bray_dist ~ sample_data(ps_forward_sans_controls_sans_SVS)$Treatment)
adonis2(bray_dist ~ sample_data(ps_forward_sans_controls_sans_SVS)$SampleEnvironment)

# make this prettier
library(ggordiplots)
gg_ordiplot(ordination, groups = ps_forward_sans_controls_relab$Management, pt.size = 3)

#Look at just water column samples over space
ps_forward_sans_controls_water <- subset_samples(ps_forward_sans_controls, Timepoint == "T0" | Timepoint == "TE")
ps_forward_sans_controls_water <- subset_samples(ps_forward_sans_controls_water, SampleEnvironment == "Water")
ps_forward_sans_controls_water
ps_forward_sans_controls_relab <- transform_sample_counts(ps_forward_sans_controls_water, function(x) x / sum(x) )
bray_dist <- phyloseq::distance(ps_forward_sans_controls_relab, method="bray") # calculate bray-curtis metric
ordination <- ordinate(ps_forward_sans_controls_relab, method="PCoA", distance=bray_dist) # Perform ordination using bray-curtis metric
plot_ordination(ps_forward_sans_controls_relab, ordination, color="Lon") + theme(aspect.ratio=1) # and lets make sure its "pretty"  
# nothing of real significance here

#11
## differential abundance
library(ALDEx2)
library(tidyr)
library(tibble)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
# compare how the communities vary at time t2, across sample location if we can
barrel_one_0_2 <- subset_samples(ps_forward_sans_controls_sans_SVS) # subset our data 
barrel_one_0_2 <- subset_samples(barrel_one_0_2, SampleEnvironment != "Control") # subset our data 
barrel_one_0_2 <- prune_taxa(taxa_sums(barrel_one_0_2) > 0, barrel_one_0_2) # remove ASVs not present
barrel_one_0_2
barrel_one_0_2 <- tax_glom(barrel_one_0_2, taxrank="Taxon3")
aldex2_sample_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(barrel_one_0_2)),
                                    phyloseq::sample_data(barrel_one_0_2)$Site,
                                    test="kw")
aldex2_sample_site$OTU <- row.names(aldex2_sample_site)
row.names(aldex2_sample_site) <- NULL
aldex_taxa_info <- data.frame(tax_table(barrel_one_0_2)) # convert taxonomy table to dataframe
aldex_taxa_info <- aldex_taxa_info %>% rownames_to_column(var = "OTU") # change rownames to column named OTU
aldex2_sample_site_sig <- aldex2_sample_site %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
dim(aldex2_sample_site_sig)
sig_aldex2_gen_result_location <- left_join(aldex2_sample_site_sig, aldex_taxa_info) # append taxonomy info
gen_otu_table_location <- data.frame(phyloseq::otu_table(barrel_one_0_2)) # extract otu table
gen_otu_table_location <- rownames_to_column(gen_otu_table_location, var = "OTU") # covert rownames to column

sig_aldex2_gen_result_location_tax <- left_join(sig_aldex2_gen_result_location, gen_otu_table_location) # combine the tables
rownames(sig_aldex2_gen_result_location_tax) <- sig_aldex2_gen_result_location_tax$Taxon3 # change Taxon22 (Genus) to rownames
sig_aldex2_gen_result_location_tax <- sig_aldex2_gen_result_location_tax[, -c(1:11)] # remove tax info from otu table
shsk_gen_czm_location <- (apply(sig_aldex2_gen_result_location_tax, 1, function(x){log(x+1) - mean(log(x+1))})) # log transformation
Z.Score.gen_location <- scale(t(shsk_gen_czm_location)) # scale the dataset

df <- sample_data(barrel_one_0_2) # extract metadata from phyloseq object
list=df[order(df$Site, decreasing=T),] # order based on Site so we can group them together
sample=rownames(list) # extract the rownames from the list (SampleID)
heatmap_annotation_top = HeatmapAnnotation(
  Location = df$Site) # We will annotate by site 
col_matrix <- brewer.pal(6, "BrBG")
hm_gen_location <- Heatmap(Z.Score.gen_location)
hm_gen_location
hm_gen_location <- Heatmap(Z.Score.gen_location, name = "Z-score", col = col_matrix, 
                           top_annotation = heatmap_annotation_top,
                           cluster_columns = F, 
                           column_names_gp = gpar(fontsize = 6))
ggsave(hm_gen_location, "diff_abundance_based_on_site.pdf")

#12
# looking for danger bugs
write.table(taxonomy_fixed_df, file = "taxonomy_fixed.txt", quote = F, sep = ",")
tax_coliform <- read.delim("taxonomy_fixed_coliforms.txt")
taxa_table_new <- tax_coliform %>% 
  tibble::column_to_rownames("ASV")
taxa_mat_new <- as.matrix(taxa_table_new)
ps_forward_new <- phyloseq(physeq_asv, tax_table(taxa_mat_new), physeq_meta)
ps_forward_new
ps_forward_new <- subset_samples(ps_forward_new, Timepoint == "T0" | Timepoint == "TE")
# need to determine where the bugs occur
ps_forward_ra <- transform_sample_counts(ps_forward_new, function(x) x / sum(x)) # convert to relative abundance
p <- plot_bar(ps_forward_ra, "Sample", fill="Putative.Fecal.Coliform") + # Sample is x-axis
  geom_bar(aes(color=Putative.Fecal.Coliform), stat="identity") + # specifies the type of graph and color
  facet_grid(~Timepoint, scales = "free", space = "free_x") # lets separate our samples by Site
p
ggsave(plot = p, filename = "phyla_taxonomy_Putative.Fecal.Coliform.pdf", width = 5, height = 5, units = "in")

#filter to top 15 class and phyla
rel_abund_df <- psmelt(ps_forward_ra)
rel_abund_df <- rel_abund_df %>% 
  filter(SampleEnvironment != "Control") %>%
  filter(Putative.Fecal.Coliform == "Yes")
df2 <- rel_abund_df %>% 
  group_by(Sample) %>%
  summarise(sum = sum(Abundance))
df2$Site <- c("A","A","A","A","B","B","B","B","1","2","3","4","5")
df2$Site <- factor(df2$Site, levels = rev(c("A", "1", "2", "3", "4", "B", "5")))
p <- ggplot(df2, aes(x = Site, y = sum)) + geom_point() +
  stat_summary(fun.y=mean, aes(group=1), geom="line", colour="red") +
  stat_summary(fun.y=mean, aes(group=1), geom="point", colour="red", size=3, shape=4)
ggsave(plot = p, "PutativeFecalColiforms.pdf", width = 5, height = 5, units = "in")

## from galaxy
library(phyloseq)
asv_table <- read.delim("Galaxy1339-[dada2_ removeBimeraDenovo on data 1336_ sequencetable].dada2_sequencetable", sep =  "\t", header = 1, check.names=F)
taxa_table <- read.delim("Galaxy1341-[dada2_ assignTaxonomy and addSpecies on data 1339].tabular", sep =  "\t", header = 1, check.names=F)
metadata <- read.delim("Galaxy823-[16smetadata.tsv].tabular", sep =  "\t", header = 1, check.names=F)

names(asv_table)[1] <- "X"
names(taxa_table)[1] <- "X"

asv_table <- asv_table %>%
  tibble::column_to_rownames("X") 
taxa_table <- taxa_table %>% 
  tibble::column_to_rownames("X")
metadata <- metadata %>% 
  tibble::column_to_rownames("SampleID") 

asv_mat <- as.matrix(asv_table)
taxa_mat <- as.matrix(taxa_table)

OTU = otu_table(asv_mat, taxa_are_rows = TRUE)
TAX = tax_table(taxa_mat)
samples = sample_data(metadata)
test <- phyloseq(OTU, TAX, samples)
test


# making maps
library(leaflet)
library(htmlwidgets)
library(webshot)
m = leaflet() %>% addTiles()
m
m = m %>% setView(-70.635, 41.5755, zoom = 16)
m
m <- m %>% addMarkers(-70.63537, 41.576007, popup = "Site A") %>% addMarkers(-70.640487, 41.577018, popup = "Site B")
m
saveWidget(m, "SiteMap.html", selfcontained = FALSE)
webshot("SiteMap.html", file = "SiteMap.png",
        cliprect = "viewport", vwidth = 750, vheight = 500)

m = leaflet() %>% addTiles()
m
m = m %>% setView(-70.635, 41.5755, zoom = 16)
m
m <- m %>% addMarkers(-70.63537, 41.576007, popup = "Site A") %>% addMarkers(-70.640487, 41.577018, popup = "Site B") %>% addMarkers(-70.636115, 41.575615, popup = "Site 1") %>% addMarkers(-70.636297,41.576419, popup = "Site 2") %>% addMarkers(-70.636759,41.5757759, popup = "Site 3") %>% addMarkers(-70.6381255, 41.5766886, popup = "Site 4") %>% addMarkers(-70.6417338, 41.5759944, popup = "Site 5")
m
saveWidget(m, "SiteMap_2.html", selfcontained = FALSE)
webshot("SiteMap_2.html", file = "SiteMap_2.png",
        cliprect = "viewport", vwidth = 750, vheight = 500)


library(decontam)
library(qiime2R)
library(tidyr)
getwd()
list.files() 
asvs <- read_qza("../BT/results/03-dada2-feature-table_forward.qza")
taxonomy <- read_qza("05-taxonomy-blast-90-1.qza")
tree <- read_qza("07-rooted-tree.qza")
asv_df <- asvs$data # we can view and save the actual data by specifying '$'
taxonomy_df <- taxonomy$data # save taxonomy info as a dataframe
phylo_tree <- tree$data

taxonomy_fixed_df <- taxonomy_df %>% separate_wider_delim(Taxon, delim = ";", names_sep = "", too_few = "align_start")
taxonomy_fixed_df[is.na(taxonomy_fixed_df)] <- "Unassigned" # rename NAs into unassigned
head(taxonomy_fixed_df)

taxonomy_fixed_df <- as.data.frame(taxonomy_fixed_df) # force into a dataframe
row.names(taxonomy_fixed_df) <- taxonomy_fixed_df$Feature.ID # make first column into row names
taxonomy_fixed_df$Feature.ID <- NULL # remove the first column
taxonomy_matrix <- as.matrix(taxonomy_fixed_df) # convert to a matrix 

physeq_asv <- otu_table(asv_df, taxa_are_rows = T) # convert into phyloseq object
physeq_tax <- tax_table(taxonomy_matrix) # convert into phyloseq object
physeq_meta <- sample_data(
