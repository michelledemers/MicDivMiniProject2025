# Bioinformatics Protocol for 16S Amplicon Sequencing
## Michelle DeMers
### August 2, 2025 - August 9, 2025

This workflow is working with a demultiplexed dataset of 16S amplicon sequences, having amplified the 16S V4 region. This is from a small pertubation pilot experiment, where water column and sediment samples were taken from a salt marsh environment and then varying nitrogen sources were added over a short period of time to assess how 'untreated' wastewater pertubation might affect microbial communities in these tow spaces, relative to their initial composition and, possibly, predicted function.

## Setting up the conda environment and all my data/directories

From [here](https://library.qiime2.org/quickstart/amplicon).

		srun --pty  --cpus-per-task=1 --job-name=interact --ntasks=1 --nodes=1 --partition=batch --time=02:00:00 --mem=12GB /bin/bash -l
		module load Miniforge3
		mkdir /home/md89517/conda-env
		conda env create \
	  --name qiime2-amplicon-2025.7 \
	  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.7/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml
		
Then we'll have the folder we are working in/with.
		
		cd /scratch/md89517/
		mkdir MD2025-amplicon-project
		cd MD2025-amplicon-project/
		mkdir data databases metadata results scripts
		
## Scripts for running the pipeline

### 01-Assessing-Data-Quality

Ahead of all of the below code, we need to activate our conda environment that we created above.

		module load Miniforge3
		source activate qiime2-amplicon-2025.7

Here we want to convert all of out fastqc files into qiime2 artifacts, and then we'll generate a report of the quality of the sequencing run within that.

		qiime tools import \
		--type 'SampleData[PairedEndSequencesWithQuality]' \
		--input-path fq-manifest.tsv \
		--output-path demux.qza \
		--input-format PairedEndFastqManifestPhred33V2
		 
The above command requires a manifest. You can use this OR, you can use the below code:

		qiime tools import \
	  --type 'SampleData[PairedEndSequencesWithQuality]' \
	  --input-path ${INPUT1} \
	  --output-path ${OUTPUT1} \
	  --input-format CasavaOneEightSingleLanePerSampleDirFmt
	  
Then, we need to visualize the quality of the data with the below code.

		qiime demux summarize \
	  --i-data ${INPUT} \
	  --o-visualization ${OUTPUT}

### 02-Assessing-Data-Quality with FastQC and MultiQC

We'll also use fastQC and multiQC to assess data quality, just to verify that the data quality is similar to what we get from qiime before trimming in the next step.

	module load FastQC
	INPUT=/scratch/md89517/MD2025-amplicon-project/data
	OUTPUT=/scratch/md89517/MD2025-amplicon-project/results/02-fastqc
	mkdir -p ${OUTPUT}
	fastqc ${INPUT}/* -o ${OUTPUT}

	module load MultiQC
	OUTPUT2=/scratch/md89517/MD2025-amplicon-project/results/02-multiqc
	mkdir -p ${OUTPUT2}
	multiqc --outdir ${OUTPUT2} ${OUTPUT}

Confirmed, data is okay. Wish it was better, but we move on.
	
### 03-DADA2-Denoising and 04-Denoising-Visualization

We want to trim adapters, primers, and low quality sequences, and then we want to merge our sequences. Being that the sequence lengths we have are 150 nts....we may not want to trim a lot, if anything. We also know the forward and reverse reads do not overlap, so that will be interesting in this next step. More information [here](https://docs.qiime2.org/jupyterbooks/cancer-microbiome-intervention-tutorial/020-tutorial-upstream/040-denoising.html) and [here](https://benjjneb.github.io/dada2/tutorial.html), for CLI and R, respectively.

For now, we will try to stay in command line and use:

		module load Miniforge3
		source activate /home/md89517/.conda/envs/qiime2-amplicon-2025.7
	
	  qiime dada2 denoise-paired \
	  --i-demultiplexed-seqs /scratch/md89517/MD2025-amplicon-project/results/01-data-quality-reports/qiime-imports.qza \
	  --p-trunc-len-f 150 \
	  --p-trunc-len-r 150 \
	  --o-representative-sequences /scratch/md89517/MD2025-amplicon-project/results/03-dada2-rep-seq.qza \
	  --o-table /scratch/md89517/MD2025-amplicon-project/results/03-dada2-feature-table.qza \
	  --o-denoising-stats /scratch/md89517/MD2025-amplicon-project/results/03-dada2-stats.qza
	  --p-min-overlap 4
	
		qiime metadata tabulate \
	  --m-input-file /scratch/md89517/MD2025-amplicon-project/results/03-dada2-stats.qza \
	  --o-visualization /scratch/md89517/MD2025-amplicon-project/results/04-dada2-stats-summ.qzv
	  
Based on `04-dada2-stats-summ.qzv`, nothing really merged well. It would seem that, on average, 25% passed the filter and then ~10% merged, and then really nothing remained after removing chimeras. So, we are 1) going to move into dada2 pipeline in R, which allows us to concatenate the reads as opposed to relying on an overlap, like the qiime2 pipeline requires in command line and 2) we will use only the forward reads. So, we have two pipelines going forward, which will hopefully merge at somepoint in R. (You can also do this step in Galaxy. To work with the phyloseq output from Galaxy, use [this](https://vaulot.github.io/tutorials/Phyloseq_tutorial.html).)

#### Pipeline 1 (R)

From above, we are working with the tutorial from [here](https://benjjneb.github.io/dada2/tutorial.html). This allows the reads to merge via concatenation, so you can use both forward and reverse reads. In R:

	setwd("fastqFiles/")
	fnFs <- sort(list.files(getwd(), pattern="_R1_001.fastq", full.names = TRUE))
	fnRs <- sort(list.files(getwd(), pattern="_R2_001.fastq", full.names = TRUE))
	# view qaultiy of first two samples, forward files
	dada2::plotQualityProfile(fnFs[1:2])
	sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
	#filter step
	filtFs <- file.path(getwd(), "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
	filtRs <- file.path(getwd(), "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
	names(filtFs) <- sample.names
	names(filtRs) <- sample.names
	filtFs <- file.path(getwd(), "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
	filtRs <- file.path(getwd(), "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
	names(filtFs) <- sample.names
	names(filtRs) <- sample.names
	out <- dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150),
	                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
	                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
	head(out)
	
	# First, rerun the filtering , but remove the file "25-MAD-B1S-NO3-T1_S58_L001_R1_001.fastq"
	sample.names2 <- sample.names[-grep("25-MAD-B1S-NO3-T1", sample.names)]
	filtFs <- file.path(getwd(), "filtered", paste0(sample.names2, "_F_filt.fastq.gz"))
	filtRs <- file.path(getwd(), "filtered", paste0(sample.names2, "_R_filt.fastq.gz"))
	names(filtFs) <- sample.names2
	names(filtRs) <- sample.names2
	out_filt <- dada2::filterAndTrim(fnFs[-grep("25-MAD-B1S-NO3-T1", sample.names)], filtFs, fnRs[-grep("25-MAD-B1S-NO3-T1", sample.names)], filtRs, truncLen=c(150,150),
	                            maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
	                            compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
	head(out_filt)
	# Then we make our ASVs using the filtered and the non-filtered reads, just to see what we get:
	
	# Make ASVs and ASV table from filtered results:
	errF_filt <- dada2::learnErrors(filtFs, multithread=TRUE)
	errR_filter <- dada2::learnErrors(filtRs, multithread=TRUE)
	dadaFs_filt <- dada2::dada(filtFs, err=errF_filt, multithread=TRUE)
	dadaRs_filt <- dada2::dada(filtRs, err=errR_filter, multithread=TRUE)
	mergers_filt <- dada2::mergePairs(dadaFs_filt, filtFs, dadaRs_filt, filtRs, justConcatenate = TRUE)
	head(mergers_filt[2])
	seqtab_filt <- dada2::makeSequenceTable(mergers_filt)
	dim(seqtab_filt)
	
	# Make ASVs and ASV table from NON-filtered results:
	errF <- dada2::learnErrors(fnFs, multithread=TRUE)
	errR <- dada2::learnErrors(fnRs, multithread=TRUE)
	dadaFs <- dada2::dada(fnFs, err=errF, multithread=TRUE)
	dadaRs <- dada2::dada(fnRs, err=errR, multithread=TRUE)
	mergers <- dada2::mergePairs(dadaFs, fnFs, dadaRs, fnRs, justConcatenate = TRUE)
	head(mergers[2])
	seqtab <- dada2::makeSequenceTable(mergers)
	dim(seqtab)
	
This does create some ASVs, we'll see what we get moving forward.

#### Pipeline 2 (CLI)

Uses forward reads only, and trims all forward reads at nt 140.

		module load Miniforge3
		source activate /home/md89517/.conda/envs/qiime2-amplicon-2025.7

		qiime dada2 denoise-single \
		--i-demultiplexed-seqs /scratch/md89517/MD2025-amplicon-project/results/01-data-quality-reports/qiime-imports.qza \
		--p-trunc-len 140 \
		--p-trim-left 0 \
		--o-table /scratch/md89517/MD2025-amplicon-project/results/03-dada2-feature-table_forward.qza \
		--o-representative-sequences /scratch/md89517/MD2025-amplicon-project/results/03-dada2-rep-seq_forward.qza \
		--o-denoising-stats /scratch/md89517/MD2025-amplicon-project/results/03-dada2-stats_forward.qza \
		--verbose
	
		qiime metadata tabulate \
	  --m-input-file /scratch/md89517/MD2025-amplicon-project/results/03-dada2-stats_forward.qza \
	  --o-visualization /scratch/md89517/MD2025-amplicon-project/results/04-dada2-stats-summ_forward.qzv
	  
This maintained about 50% of all my reads. We will see how this turns out going forward.

### 05-Taxonomy-Classification

We'll be using Silva in both pipelines.

#### Pipeline 1 (R)

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
	
	## taxa
	taxa_filt <- assignTaxonomy(seqtab_filt.nochim, "~/tax/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
	taxa_filt.print <- taxa_filt # Removing sequence rownames for display only
	rownames(taxa_filt.print) <- NULL
	head(taxa_filt.print)
	
	taxa <- assignTaxonomy(seqtab.nochim, "~/tax/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
	taxa.print <- taxa # Removing sequence rownames for display only
	rownames(taxa.print) <- NULL
	head(taxa.print)

Well, there's something.

#### Pipeline 2 (CLI)

In an ideal world, our server would have a qiime classifier (from any given source, like GTDB or Silva) built for our region of interest, and we wouldn't have to install it on our own user. However, that is not the case, so we will go get a pre-built classifier. There are caveats and not on a time crunch, I might build my own, but for now, using the qiime pipeline, we would find one that is publicly available. We would use a Silva classifier, from [here](https://forum.qiime2.org/t/silva-v4-classifier/32410/6), which links to [here](https://bwsyncandshare.kit.edu/s/zK5zjAbsTFQRpbo?path=%2Fq2_24-10%2FSilva138_2). From this link, we are using the *V4.qza files. This is technically an older version of a qiime object, so let's see if it will even run with our newer version of qiime2. The three classifier files are currently installed in `/scratch/md89517/MD2025-amplicon-project/databases/`. 

		module load Miniforge3
		source activate qiime2-amplicon-2025.7
	
		INPUT=/scratch/md89517/MD2025-amplicon-project/results/03-dada2-rep-seq_forward.qza
		REFSEQ=/scratch/md89517/MD2025-amplicon-project/databases/seqs-V4.qza
		REFTAX=/scratch/md89517/MD2025-amplicon-project/databases/taxa-V4.qza
		CLASS=/scratch/md89517/MD2025-amplicon-project/results/05-taxonomy-blast-90-1.qza
		SEARCH=/scratch/md89517/MD2025-amplicon-project/results/05-taxonomy-seach-results
	
		qiime feature-classifier classify-consensus-blast \
		  --i-query ${INPUT} \
		  --i-reference-taxonomy ${REFTAX} \
		  --i-reference-reads ${REFSEQ} \
		  --p-maxaccepts 1 \
		  --p-perc-identity 0.90 \
		  --o-classification ${CLASS} \
		  --o-search-results ${SEARCH} \
		  --p-num-threads 12
		  
The output here suggests we have...diversity? At the very least it classified things. We continue on. :)

### 06-Importing-Metadata 

#### Pipeline 1 (R)

Metadata is saved, currently, in the MiniProject folder, as a txt. We'll properly import when the two pipelines converge in a couple steps.

#### Pipeline 2 (CLI)

We'll add the metadata to the table summary we generated before, and visualize.

	module load Miniforge3
	source activate qiime2-amplicon-2025.7
	
	INPUT=/scratch/md89517/MD2025-amplicon-project/results/03-dada2-feature-table_forward.qza
	OUTPUT=/scratch/md89517/MD2025-amplicon-project/results
	METADATA=/scratch/md89517/MD2025-amplicon-project/metadata/16smetadata.txt
	
	qiime feature-table summarize \
  --i-table ${INPUT} \
  --o-visualization ${OUTPUT}/06-feature-table.qzv \
  --m-sample-metadata-file ${METADATA}

### 07-Phylogeny-Building

#### Pipeline 1

There is no tree building step in this pipeline. TBH not even sure why it exists because you'd think you'd want different trees per sample/treatment/site, anyway? Unclear.

#### Pipeline 2 (CLI)

This is the last step before we merge and decide what reads/filtering we will be using! This is very exciting.

	module load Miniforge3
	source activate qiime2-amplicon-2025.7

	INPUT=/scratch/md89517/MD2025-amplicon-project/results/03-dada2-rep-seq.qza
	OUTPUT=/scratch/md89517/MD2025-amplicon-project/results

	qiime alignment mafft \
	  --i-sequences ${INPUT} \
	  --o-alignment ${OUTPUT}/07-alignment.qza
	
	qiime phylogeny fasttree \
	  --i-alignment ${OUTPUT}/07-alignment.qza \
	  --o-tree  ${OUTPUT}/07-unrooted-tree.qza
	
	qiime phylogeny midpoint-root \
	  --i-tree ${OUTPUT}/07-unrooted-tree.qza \
	  --o-rooted-tree ${OUTPUT}/07-rooted-tree.qza
	  
### 08-Moving-to-Phyloseq

#### Pipeline 1

	library(phyloseq)
	library(ggplot2)
	library(Biostrings)
	library(dplyr)
	
	samdf <- read.table("16smetadata_qiime.txt", header = 1)
	samdf_filt <- samdf[-grep("25-MAD-B1S-NO3-T1", samdf$sample.id),]
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

Here we have our filtered one and our non-filtered one.

#### Pipeline 2 (CLI)

	library(decontam)
	library(qiime2R)
	library(tidyr)
	samdf <- read.table("16smetadata_qiime.txt", header = 1)
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
	

The rest of this will be in R, and as specified above, we'll be using our forward reads only.

### 09-Decontamination

	# look at size of libraries in controls versus true samples
	df <- as.data.frame(sample_data(ps_forward)) # Put sample_data into a ggplot-friendly data.frame
	df$LibrarySize <- sample_sums(ps_forward)
	df <- df[order(df$LibrarySize),]
	df$Index <- seq(nrow(df))
	ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()
	
	sample_data(ps_forward)$is.neg <- sample_data(ps_forward)$Sample_or_Control == "Control Sample"
	contamdf.prev <- isContaminant(ps_forward, method="prevalence", neg="is.neg")
	table(contamdf.prev$contaminant)
	
	contamdf.prev05 <- isContaminant(ps_forward, method="prevalence", neg="is.neg", threshold=0.5)
	table(contamdf.prev05$contaminant)

No contaminants get removed. Lol. We keep moving.

### 10-Alpha-and-Beta-Diversity

	# start with relative abundance taxonomy plots
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
	
	# alpha diversity plots
		
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
	

	#Look at just water column samples over space
	ps_forward_sans_controls_water <- subset_samples(ps_forward_sans_controls, Timepoint == "T0" | Timepoint == "TE")
	ps_forward_sans_controls_water <- subset_samples(ps_forward_sans_controls_water, SampleEnvironment == "Water")
	ps_forward_sans_controls_water
	ps_forward_sans_controls_relab <- transform_sample_counts(ps_forward_sans_controls_water, function(x) x / sum(x) )
	bray_dist <- phyloseq::distance(ps_forward_sans_controls_relab, method="bray") # calculate bray-curtis metric
	ordination <- ordinate(ps_forward_sans_controls_relab, method="PCoA", distance=bray_dist) # Perform ordination using bray-curtis metric
	plot_ordination(ps_forward_sans_controls_relab, ordination, color="Lon") + theme(aspect.ratio=1) # and lets make sure its "pretty"  
	# nothing of real significance here

### 11-Differential-Abundance

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
	                                    phyloseq::sample_data(barrel_one_0_2)$Site, # across sample location!
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
