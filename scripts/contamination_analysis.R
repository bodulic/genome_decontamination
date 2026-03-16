#!/usr/bin/env Rscript
#Contamination analysis script
#Loading the required packages
library(data.table)
library(taxonomizr)
library(Biostrings)
library(compositions)
library(mclust)
library(umap)
library(ggplot2)
library(GenomicRanges)

#Importing external arguments
ext_args <- commandArgs(trailingOnly = T)
THREAD_N <- as.numeric(ext_args[1])
TAXONOMIZR_DB <- ext_args[2]
GENOME_FASTA_PATH <- ext_args[3]
SCAFF_LEN_MIN <- as.numeric(ext_args[4])
SEED <- as.numeric(ext_args[5])
MCLUST_G_MAX <- as.numeric(ext_args[6])
CHISQ_ALPHA <- as.numeric(ext_args[7])
GC_Z_THR <- as.numeric(ext_args[8])
DIAMOND_RES_PATH <- ext_args[9]
TARGET_TAXON <- ext_args[10]
DIA_WEIGHT_RATIO_MULT <- as.numeric(ext_args[11])
DIA_COV_FRACTION_THR <- as.numeric(ext_args[12])
MEGAN_LCA_RES_PATH <- ext_args[13]
KRAKEN_RES_PATH <- ext_args[14]
KRAKEN_KMER_FRACTION_RATIO_THR <- as.numeric(ext_args[15])
KRAKEN_KMER_FRACTION_THR <- as.numeric(ext_args[16])
KRAKEN_KMER_FRACTION_TAX_ASSIGN_THR <- as.numeric(ext_args[17])
KRAKEN_CLASSIF_KMER_FRACTION_TAX_ASSIGN_THR <- as.numeric(ext_args[18])
DIA_COV_FRACTION_TAX_ASSIGN_THR <- as.numeric(ext_args[19])
DIA_REL_COV_FRACTION_TAX_ASSIGN_THR <- as.numeric(ext_args[20])

setDTthreads(THREAD_N)

#Defining functions
map_tax_categories <- function(domain, phylum) {
 if (domain == "Viruses" | (is.na(domain) & is.na(phylum) == F)) return("Viruses")
 if (domain == "Bacteria") return("Bacteria")
 if (domain == "Archaea") return("Archaea")
 if (domain != "Eukaryota") return("Other")  
  
 metazoa <- c("Porifera", "Cnidaria", "Ctenophora", "Platyhelminthes" ,"Annelida", "Mollusca", "Nematoda", "Arthropoda", "Echinodermata", "Chordata", "Brachiopoda", "Bryozoa", "Hemichordata", "Tardigrada", "Onychophora", "Rotifera", "Xenacoelomorpha", "Nemertea", "Priapulida", "Placozoa", "Phoronida", "Acanthocephala", "Chaetognatha", "Nematomorpha", "Entoprocta", "Dicyemida", "Orthonectida", "Loricifera", "Kinorhyncha", "Gnathostomulida", "Micrognathozoa", "Gastrotricha", "Cycliophora")
 fungi <- c("Ascomycota", "Basidiomycota", "Glomeromycota", "Mucoromycota", "Zoopagomycota", "Neocallimastigomycota", "Blastocladiomycota", "Chytridiomycota", "Microsporidia", "Cryptomycota", "Olpidiomycota")
 viridiplantae <- c("Streptophyta", "Chlorophyta", "Prasinodermophyta")
 haptophyta <- c("Haptophyta")
 rhodophyta <- c("Rhodophyta")
 stramenopiles <- c("Ochrophyta", "Bigyra", "Hyphochytriomycota", "Oomycota", "Bacillariophyta")
 alveolata <- c("Dinoflagellata", "Apicomplexa", "Ciliophora", "Perkinsozoa")
 excavata <- c("Euglenozoa", "Metamonada", "Discoba", "Malawimonada", "Preaxostyla", "Heterolobosea", "Fornicata")
 rhizaria <- c("Retaria", "Cercozoa", "Endomyxa", "Foraminifera")
 choanoflagellates <- c("Choanoflagellata", "Choanoflagellatea")
 amoebozoa <- c("Amoebozoa", "Evosea", "Discosea") 
  
 if (phylum %chin% metazoa)           return("Metazoa")
 if (phylum %chin% fungi)             return("Fungi")
 if (phylum %chin% viridiplantae)     return("Viridiplantae")
 if (phylum %chin% haptophyta)        return("Haptophyta")
 if (phylum %chin% rhodophyta)        return("Rhodophyta")
 if (phylum %chin% stramenopiles)     return("Stramenopiles")
 if (phylum %chin% alveolata)         return("Alveolata")
 if (phylum %chin% excavata)          return("Excavata") 
 if (phylum %chin% rhizaria)          return("Rhizaria")
 if (phylum %chin% choanoflagellates) return("Choanoflagellata")
 if (phylum %chin% amoebozoa)         return("Amoebozoa")
 return("Eukaryote_other")
}

assign_taxon <- function(input_dt) {
 taxons <- getTaxonomy(input_dt[, taxid], sqlFile = TAXONOMIZR_DB)
 input_dt[, ':=' ("domain" = taxons[, 1], "phylum" = taxons[, 2], "class" = taxons[, 3], "order" = taxons[, 4], "family" = taxons[, 5], "genus" = taxons[, 6])]
 input_dt <- input_dt[is.na(domain) == F & is.na(phylum) == F]
 input_dt[, "taxon" := mapply(map_tax_categories, domain, phylum)]
 input_dt[, "taxid" := NULL]
 return(input_dt)
}

#K-mer profile analysis
genomic_scaffolds <- readDNAStringSet(GENOME_FASTA_PATH, format = "fasta")
genomic_scaffolds <- genomic_scaffolds[width(genomic_scaffolds) >= SCAFF_LEN_MIN]
names(genomic_scaffolds) <- sub(" .*", "", names(genomic_scaffolds))
scaffold_len_dt <- data.table("scaffold_name" = names(genomic_scaffolds), "scaffold_len" = width(genomic_scaffolds))

kmer_count <- as.data.table(oligonucleotideFrequency(genomic_scaffolds, width = 5, step = 1, as.prob = F))
kmer_count[, "scaffold" := scaffold_len_dt[, scaffold_name]]
kmer_count <- melt(kmer_count, id.vars = "scaffold",  variable.name = "kmer", value.name = "count")
kmer_count[, "represent_kmer" := pmin(as.character(kmer), reverseComplement(DNAStringSet(kmer)))]

kmer_count <- kmer_count[, .("count" = sum(count)), by = c("scaffold", "represent_kmer")]
kmer_count <- dcast(kmer_count, scaffold ~ represent_kmer, value.var = "count", fill = 0)
kmer_count_scaff <- kmer_count[, scaffold]
kmer_count <- as.matrix(kmer_count[, !"scaffold"])
rownames(kmer_count) <- kmer_count_scaff
rm(kmer_count_scaff)

kmer_count <- kmer_count + 0.5        
kmer_prop <- kmer_count / rowSums(kmer_count)
kmer_prop_clr <- clr(acomp(kmer_prop))
rm(kmer_count, kmer_prop)

kmer_prop_clr_pcs <- prcomp(kmer_prop_clr, center = T, scale. = F)
var_exp <- kmer_prop_clr_pcs$sdev^2 / sum(kmer_prop_clr_pcs$sdev^2)
pc_n <- which(cumsum(var_exp) >= 0.9)[1]
pc_n <- min(pc_n, 20)
kmer_prop_clr_pcs <- kmer_prop_clr_pcs$x[, 1 : pc_n]
rm(kmer_prop_clr, var_exp, pc_n)

set.seed(SEED)
mc_fits <- mclustICL(kmer_prop_clr_pcs, G = 1 : MCLUST_G_MAX, modelNames = c("EII", "EEI", "EEE", "VII", "VEI", "VVI"))
best_model_indices <- which(mc_fits == max(mc_fits, na.rm = T), arr.ind = T)
best_model_indices <- best_model_indices[which.min(best_model_indices[, "row"]), ] 
best_g <- as.integer(best_model_indices[1])

if(best_g > 1) {
 best_model_name  <- colnames(mc_fits)[best_model_indices[2]]
 set.seed(SEED)
 mc_fit_final <- Mclust(kmer_prop_clr_pcs, G = best_g, modelNames = best_model_name)
 if(is.null(mc_fit_final)) {
  stop("Mclust failed - retry with a lower max G value (1-30). Exiting.")
 }
  
 mc_clusters <- mc_fit_final$classification
 mc_umap <- umap(kmer_prop_clr_pcs, n_neighbors = 30, min_dist = 0.1)
 mc_umap <- as.data.table(cbind("scaffold_name" = rownames(mc_umap$data), "UMAP1" = mc_umap$layout[, 1], "UMAP2" = mc_umap$layout[, 2], "Cluster" = as.factor(mc_clusters)))
 mc_umap[, "Cluster" := as.character(Cluster)]
  
 mc_umap_plot <- ggplot(data = mc_umap, aes(UMAP1, UMAP2, color = Cluster)) +
  geom_point(size = 0.25, alpha = 0.7) +
  theme_minimal() +
  theme(panel.grid = element_blank())  +
  theme(panel.background = element_rect(fill = NA, color = "grey70")) +
  theme(legend.title = element_text(size = 7)) +
  theme(legend.text = element_text(size = 6)) +
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(axis.title = element_text(size = 7)) +
  theme(axis.text = element_text(size = 6))
  
 ggsave(mc_umap_plot, file = "mc_umap_plot.tiff", height = 2.5, width = 3.5, dpi = 300, bg = "white")
 write.table(mc_umap, file = "mc_umap_dt.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
 rm(best_model_name, mc_umap, mc_umap_plot)
  
 mc_clusters <- data.table("scaffold_name" = names(mc_clusters), "cluster" = mc_clusters)
 scaffold_info_dt <- merge(scaffold_len_dt, mc_clusters)
 mc_cluster_tot_len <- scaffold_info_dt[, .("scaff_len_sum" = sum(scaffold_len)), by = "cluster"]
 setorder(mc_cluster_tot_len, -scaff_len_sum)
 rm(mc_clusters)
  
 centers <- mc_fit_final$parameters$mean
 ref_cluster <- mc_cluster_tot_len[1, cluster]
 ref_center <- centers[, ref_cluster]
 cov_matrix <- mc_fit_final$parameters$variance$sigma[,, ref_cluster, drop = F][,, 1] 
  
 distance_to_ref_center <- sapply(1 : ncol(centers), function(i) {
  sqrt(mahalanobis(x = centers[, i], center = ref_center, cov = cov_matrix))
 })
 rm(mc_fit_final, mc_cluster_tot_len, ref_cluster, ref_center, cov_matrix)
  
 distance_thr <- sqrt(qchisq(1 - CHISQ_ALPHA, df = nrow(centers)))
 outlier_clusters <- which(distance_to_ref_center > distance_thr)
 scaffold_info_dt[, "outlier_kmer" := fifelse(cluster %in% outlier_clusters, T, F)]
 rm(centers, distance_to_ref_center, distance_thr, outlier_clusters)
} else {
 scaffold_info_dt <- scaffold_len_dt
 scaffold_info_dt[, ':=' ("cluster" = 1, "outlier_kmer" = F)]
}
rm(kmer_prop_clr_pcs, mc_fits, best_model_indices, best_g)

#GC content analysis
nucl_freq <- letterFrequency(genomic_scaffolds, letters = c("A", "T", "G", "C"), as.prob = T)
gc_fraction <- nucl_freq[, "G"] + nucl_freq[, "C"]
scaffold_info_dt[, "gc_content" := gc_fraction]
scaffold_info_dt[, "gc_z" := (gc_content - mean(gc_content)) / sd(gc_content)]
scaffold_info_dt[, "outlier_gc" := fifelse(abs(gc_z) > GC_Z_THR, T, F)]
rm(nucl_freq, gc_fraction)

#DIAMOND consensus analysis
diamond_hits <- fread(DIAMOND_RES_PATH, header = F, select = c(1, 2, 4, 7, 12))
setnames(diamond_hits, c("scaffold_name", "subject_name", "scaffold_start", "scaffold_end", "bitscore"))
diamond_hits <- diamond_hits[scaffold_name %in% scaffold_info_dt[, scaffold_name]]
diamond_hits[scaffold_start > scaffold_end, ':=' ("scaffold_start" = scaffold_end, "scaffold_end" = scaffold_start)]
diamond_hits[, "taxid" := accessionToTaxa(subject_name, sqlFile = TAXONOMIZR_DB)]
diamond_hits <- assign_taxon(diamond_hits)

if(diamond_hits[, .N] == 0) {
 stop("None of the supplied DIAMOND hits are associated with a known taxonomy. Exiting.")
}

diamond_hits_gr <- makeGRangesFromDataFrame(diamond_hits, seqnames.field = "scaffold_name", start.field = "scaffold_start", end.field = "scaffold_end", ignore.strand = T, keep.extra.columns = T)
diamond_hit_sets <- unlist(reduce(split(diamond_hits_gr, diamond_hits_gr$taxon), ignore.strand = T, min.gapwidth = 0), use.names = T)
diamond_hit_sets$taxon <- names(diamond_hit_sets)
diamond_hit_sets$dia_hit_set_id <- 1 : length(diamond_hit_sets)

diamond_hit_sets_dt <- as.data.table(diamond_hit_sets)[, c("seqnames", "width", "taxon")]
setnames(diamond_hit_sets_dt, c("scaffold_name", "dia_hit_set_len", "taxon"))
diamond_hit_sets_dt <- diamond_hit_sets_dt[, .("dia_hit_set_cov_sum" = sum(dia_hit_set_len)), by = c("scaffold_name", "taxon")]
diamond_hit_sets_dt <- merge(diamond_hit_sets_dt, scaffold_len_dt, by = "scaffold_name")
diamond_hit_sets_dt[, "dia_hit_set_cov_fraction" := dia_hit_set_cov_sum / scaffold_len]

diamond_hits_gr_reduced_dt <- as.data.table(reduce(diamond_hits_gr, ignore.strand = T, min.gapwidth = 0))[, c("seqnames", "width")]
setnames(diamond_hits_gr_reduced_dt, c("scaffold_name", "dia_red_hits_len"))
diamond_hits_gr_reduced_dt <- diamond_hits_gr_reduced_dt[, .("dia_red_hits_len_sum" = sum(dia_red_hits_len)), by = "scaffold_name"]
diamond_hit_sets_dt <- merge(diamond_hit_sets_dt, diamond_hits_gr_reduced_dt, by = "scaffold_name")
diamond_hit_sets_dt[, "rel_dia_hit_set_cov_fraction" := dia_hit_set_cov_sum / dia_red_hits_len_sum]

diamond_hits_sets_ovl <- findOverlaps(diamond_hits_gr, diamond_hit_sets, ignore.strand = T)
diamond_hits_gr_ovl_dt <- as.data.table(diamond_hits_gr[queryHits(diamond_hits_sets_ovl)])[, c("seqnames", "bitscore", "taxon")]
diamond_hit_sets_ovl_dt <- as.data.table(diamond_hit_sets[subjectHits(diamond_hits_sets_ovl)])[, c("dia_hit_set_id", "width", "taxon")]
diamond_hits_sets_ovl_dt <- cbind(diamond_hits_gr_ovl_dt, diamond_hit_sets_ovl_dt)
setnames(diamond_hits_sets_ovl_dt, c("scaffold_name", "bitscore", "taxon", "dia_hit_set_id", "dia_hit_set_width", "taxon_red"))
diamond_hits_sets_ovl_dt <- diamond_hits_sets_ovl_dt[taxon == taxon_red]
diamond_hits_sets_ovl_dt[, "taxon_red" := NULL]
rm(diamond_hits_gr, diamond_hit_sets, diamond_hits_sets_ovl, diamond_hits_gr_ovl_dt, diamond_hit_sets_ovl_dt)

diamond_hit_sets_weighted <- diamond_hits_sets_ovl_dt[, .("dia_hit_set_weight" = mean(bitscore)), by = c("dia_hit_set_id", "dia_hit_set_width", "scaffold_name", "taxon")]
diamond_hit_sets_weighted[, "dia_hit_set_weight" := dia_hit_set_weight * sqrt(dia_hit_set_width)]
diamond_hit_sets_weighted <- diamond_hit_sets_weighted[, .("dia_hit_set_weight_sum" = sum(dia_hit_set_weight)), by = c("scaffold_name", "taxon")]
diamond_hit_sets_weighted <- merge(diamond_hit_sets_weighted, diamond_hit_sets_dt, by = c("scaffold_name", "taxon"))
diamond_hit_sets_weighted[, "taxon_cat" := fifelse(taxon == TARGET_TAXON, TARGET_TAXON, paste("Non", TARGET_TAXON, sep = "_"))]
setorder(diamond_hit_sets_weighted, -dia_hit_set_weight_sum)
diamond_hit_sets_weighted <- unique(diamond_hit_sets_weighted, by = c("scaffold_name", "taxon_cat"))
diamond_hit_sets_weighted[, "dia_hit_set_weight_ratio" := dia_hit_set_weight_sum / sum(dia_hit_set_weight_sum), by = "scaffold_name"]
diamond_hit_sets_weighted_target <- setnames(diamond_hit_sets_weighted[taxon_cat == TARGET_TAXON, c("scaffold_name", "dia_hit_set_weight_ratio", "dia_hit_set_cov_fraction", "rel_dia_hit_set_cov_fraction")], c("scaffold_name", "dia_hit_set_weight_ratio_target", "dia_hit_set_cov_fraction_target", "rel_dia_hit_set_cov_fraction_target"))
scaffold_info_dt <- merge(scaffold_info_dt, diamond_hit_sets_weighted_target, by = "scaffold_name", all.x = T)
diamond_hit_sets_weighted_non_target <- setnames(diamond_hit_sets_weighted[taxon_cat == paste("Non", TARGET_TAXON, sep = "_"), c("scaffold_name", "dia_hit_set_weight_ratio", "dia_hit_set_cov_fraction", "rel_dia_hit_set_cov_fraction")], c("scaffold_name", "dia_hit_set_weight_ratio_non_target", "dia_hit_set_cov_fraction_non_target", "rel_dia_hit_set_cov_fraction_non_target"))
scaffold_info_dt <- merge(scaffold_info_dt, diamond_hit_sets_weighted_non_target, by = "scaffold_name", all.x = T)
scaffold_info_dt[is.na(dia_hit_set_weight_ratio_target), ':=' ("dia_hit_set_weight_ratio_target" = 0, "dia_hit_set_cov_fraction_target" = 0, "rel_dia_hit_set_cov_fraction_target" = 0)]
scaffold_info_dt[is.na(dia_hit_set_weight_ratio_non_target), ':=' ("dia_hit_set_weight_ratio_non_target" = 0, "dia_hit_set_cov_fraction_non_target" = 0, "rel_dia_hit_set_cov_fraction_non_target" = 0)]
scaffold_info_dt[, "outlier_diamond" := fifelse(dia_hit_set_weight_ratio_non_target > DIA_WEIGHT_RATIO_MULT * dia_hit_set_weight_ratio_target & dia_hit_set_cov_fraction_non_target > DIA_COV_FRACTION_THR & dia_hit_set_cov_fraction_non_target > dia_hit_set_cov_fraction_target, T, F)]
rm(diamond_hit_sets_dt, diamond_hits_sets_ovl_dt, diamond_hit_sets_weighted, diamond_hit_sets_weighted_target, diamond_hit_sets_weighted_non_target)

#MEGAN LCA analysis
megan_lca_results <- fread(MEGAN_LCA_RES_PATH, header = F, select = c(1, 3))
setnames(megan_lca_results, c("scaffold_name", "taxid"))
megan_lca_results <- assign_taxon(megan_lca_results)

if(megan_lca_results[, .N] == 0) {
 stop("None of the analysed scaffolds were assigned a known taxonomy by MEGAN. Exiting.")
}

megan_lca_results <- setnames(megan_lca_results[, c("scaffold_name", "taxon")], c("scaffold_name", "taxon_lca"))
scaffold_info_dt <- merge(scaffold_info_dt, megan_lca_results, by = "scaffold_name", all.x = T)
scaffold_info_dt[, "outlier_lca" := fifelse(taxon_lca != TARGET_TAXON & is.na(taxon_lca) == F, T, F)]
rm(megan_lca_results)

#Kraken2 k-mer composition analysis
kraken_results <- fread(KRAKEN_RES_PATH, header = F, select = c(1 : 3, 5))
setnames(kraken_results, c("class_result", "scaffold_name", "scaffold_len", "kraken_taxid_kmer_counts"))
kraken_results <- kraken_results[class_result == "C"]
kraken_results[, "class_result" := NULL]

if(kraken_results[, .N] == 0) {
 stop("None of the analysed scaffolds were assigned a known taxonomy by Kraken2. Exiting.")
}

kraken_results_chunks <- split(kraken_results, cut(seq_len(kraken_results[, .N]), 10, labels = F))
kraken_results <- list()
for(i in 1 : length(kraken_results_chunks)) {
 kraken_results_chunk <- kraken_results_chunks[[i]]
 kraken_results_chunk <- kraken_results_chunk[, .("kraken_taxid_kmer_count" = unlist(strsplit(kraken_taxid_kmer_counts, " ", fixed = T))), by = "scaffold_name"]
 kraken_results_chunk[, c("taxid", "kraken_kmer_count") := tstrsplit(kraken_taxid_kmer_count, ":", fixed = T)]
 kraken_results_chunk[, "kraken_kmer_count" := as.numeric(kraken_kmer_count)]
 kraken_results_chunk[, "kraken_kmer_fraction" := kraken_kmer_count / sum(kraken_kmer_count), by = "scaffold_name"]
 kraken_results_chunk <- assign_taxon(kraken_results_chunk)
 kraken_results_chunk[, "kraken_classif_kmer_fraction" := kraken_kmer_count / sum(kraken_kmer_count), by = "scaffold_name"]
 kraken_results[[i]] <- kraken_results_chunk
 rm(kraken_results_chunk)
}
rm(kraken_results_chunks)
kraken_results <- rbindlist(kraken_results)

if(kraken_results[, .N] == 0) {
 stop("None of the supplied Kraken2 hits are associated with a known taxonomy. Exiting.")
}

kraken_kmer_fraction_scaffold <- kraken_results[, .("kraken_kmer_fraction_tax" = sum(kraken_kmer_fraction), "kraken_classif_kmer_fraction_tax" = sum(kraken_classif_kmer_fraction)), by = c("scaffold_name", "taxon")]
kraken_kmer_fraction_scaffold[, "taxon_cat" := fifelse(taxon == TARGET_TAXON, TARGET_TAXON, paste("Non", TARGET_TAXON, sep = "_"))]
setorder(kraken_kmer_fraction_scaffold, -kraken_kmer_fraction_tax)
kraken_kmer_fraction_scaffold <- unique(kraken_kmer_fraction_scaffold, by = c("scaffold_name", "taxon_cat"))
kraken_kmer_fraction_scaffold_target <- setnames(kraken_kmer_fraction_scaffold[taxon_cat == TARGET_TAXON, c("scaffold_name", "kraken_kmer_fraction_tax", "kraken_classif_kmer_fraction_tax")], c("scaffold_name", "kraken_kmer_fraction_target", "kraken_classif_kmer_fraction_target"))
scaffold_info_dt <- merge(scaffold_info_dt, kraken_kmer_fraction_scaffold_target, by = "scaffold_name", all.x = T)
kraken_kmer_fraction_scaffold_non_target <- setnames(kraken_kmer_fraction_scaffold[taxon_cat == paste("Non", TARGET_TAXON, sep = "_"), c("scaffold_name", "kraken_kmer_fraction_tax", "kraken_classif_kmer_fraction_tax")], c("scaffold_name", "kraken_kmer_fraction_non_target", "kraken_classif_kmer_fraction_non_target"))
scaffold_info_dt <- merge(scaffold_info_dt, kraken_kmer_fraction_scaffold_non_target, by = "scaffold_name", all.x = T)
scaffold_info_dt[is.na(kraken_kmer_fraction_target), ':=' ("kraken_kmer_fraction_target" = 0, "kraken_classif_kmer_fraction_target" = 0)]
scaffold_info_dt[is.na(kraken_kmer_fraction_non_target), ':=' ("kraken_kmer_fraction_non_target" = 0, "kraken_classif_kmer_fraction_non_target" = 0)]
scaffold_info_dt[, "outlier_kraken" := fifelse(kraken_kmer_fraction_non_target > KRAKEN_KMER_FRACTION_RATIO_THR * kraken_kmer_fraction_target & kraken_kmer_fraction_non_target > KRAKEN_KMER_FRACTION_THR, T, F)]
rm(kraken_kmer_fraction_scaffold, kraken_kmer_fraction_scaffold_target, kraken_kmer_fraction_scaffold_non_target)

#Integrating scores
scaffold_info_dt[, "compositional_score" := fifelse(outlier_kmer == T | outlier_gc == T, 1, 0)]
scaffold_info_dt[, "prot_taxonomy_score" := fifelse(outlier_diamond == T | (outlier_lca == T & dia_hit_set_weight_ratio_non_target > dia_hit_set_weight_ratio_target & dia_hit_set_cov_fraction_non_target > dia_hit_set_cov_fraction_target), 1, 0)]
scaffold_info_dt[, "nucl_taxonomy_score" := fifelse(outlier_kraken == T, 1, 0)]
scaffold_info_dt[, "contamination_score" := compositional_score + prot_taxonomy_score + nucl_taxonomy_score]
scaffold_info_dt[, "contamination" := fifelse(contamination_score >= 2, T, F)]
contam_scaffolds <- scaffold_info_dt[contamination == T, scaffold_name]

#Analysing contaminants
diamond_hits_contam <- diamond_hits[scaffold_name %in% contam_scaffolds]
kraken_results_contam <- kraken_results[scaffold_name %in% contam_scaffolds]
rm(diamond_hits, kraken_results)

if(diamond_hits_contam[, .N] > 0 & kraken_results_contam[, .N] > 0) {
 diamond_hits_contam_gr <- makeGRangesFromDataFrame(diamond_hits_contam, seqnames.field = "scaffold_name", start.field = "scaffold_start", end.field = "scaffold_end", ignore.strand = T, keep.extra.columns = T)
  
 contam_scaff_taxonomy <- scaffold_info_dt[contamination == T, c("scaffold_name", "scaffold_len")]
 for (tax_rank in c("genus", "family", "order", "class", "phylum", "domain")) {
  tax_rank_prot <- paste(tax_rank, "prot", sep = "_")
  tax_rank_nucl <- paste(tax_rank, "nucl", sep = "_")
    
  diamond_hit_contam_sets <- unlist(reduce(split(diamond_hits_contam_gr, mcols(diamond_hits_contam_gr)[[tax_rank]]), ignore.strand = T, min.gapwidth = 0), use.names = T)
  mcols(diamond_hit_contam_sets)[[tax_rank_prot]] <- names(diamond_hit_contam_sets)
  diamond_hit_contam_sets_dt <- as.data.table(diamond_hit_contam_sets)[, c("seqnames", "width", ..tax_rank_prot)]
  setnames(diamond_hit_contam_sets_dt, c("scaffold_name", "dia_hit_set_len", tax_rank_prot))
  diamond_hit_contam_sets_dt <- diamond_hit_contam_sets_dt[, .("dia_hit_set_cov_sum" = sum(dia_hit_set_len)), by = c("scaffold_name", tax_rank_prot)]
  diamond_hit_contam_sets_dt <- merge(diamond_hit_contam_sets_dt, scaffold_len_dt, by = "scaffold_name")
  diamond_hit_contam_sets_dt[, paste("dia_hit_set_cov_fraction", tax_rank, sep = "_") := dia_hit_set_cov_sum / scaffold_len]
  diamond_hit_contam_sets_dt <- merge(diamond_hit_contam_sets_dt, diamond_hits_gr_reduced_dt, by = "scaffold_name")
  diamond_hit_contam_sets_dt[, paste("rel_dia_hit_set_cov_fraction", tax_rank, sep = "_") := dia_hit_set_cov_sum / dia_red_hits_len_sum]
  diamond_hit_contam_sets_dt <- diamond_hit_contam_sets_dt[, c("scaffold_name", paste("dia_hit_set_cov_fraction", ..tax_rank, sep = "_"), paste("rel_dia_hit_set_cov_fraction", ..tax_rank, sep = "_"), ..tax_rank_prot)]
  diamond_hit_contam_sets_dt <- diamond_hit_contam_sets_dt[order(-diamond_hit_contam_sets_dt[[2]])]
  diamond_hit_contam_sets_dt <- unique(diamond_hit_contam_sets_dt, by = "scaffold_name")
  contam_scaff_taxonomy <- merge(contam_scaff_taxonomy, diamond_hit_contam_sets_dt, by = "scaffold_name", all.x = T)
    
  kraken_contam_kmer_fraction_scaffold <- kraken_results_contam[, setNames(list(sum(kraken_kmer_fraction), sum(kraken_classif_kmer_fraction)), c(paste("kraken_kmer_fraction", tax_rank, sep = "_"), paste("kraken_classif_kmer_fraction", tax_rank, sep = "_"))), by = c("scaffold_name", tax_rank)]
  kraken_contam_kmer_fraction_scaffold <- kraken_contam_kmer_fraction_scaffold[is.na(kraken_contam_kmer_fraction_scaffold[[2]]) == F]
  colnames(kraken_contam_kmer_fraction_scaffold)[2] <- tax_rank_nucl
  kraken_contam_kmer_fraction_scaffold <- kraken_contam_kmer_fraction_scaffold[order(-kraken_contam_kmer_fraction_scaffold[[3]])]
  kraken_contam_kmer_fraction_scaffold <- unique(kraken_contam_kmer_fraction_scaffold, by = "scaffold_name")
  contam_scaff_taxonomy <- merge(contam_scaff_taxonomy, kraken_contam_kmer_fraction_scaffold[, c(1, 3, 4, 2)], by = "scaffold_name", all.x = T)
 }
 rm(scaffold_len_dt, diamond_hits_gr_reduced_dt, diamond_hits_contam, kraken_results_contam, diamond_hits_contam_gr,  diamond_hit_contam_sets, diamond_hit_contam_sets_dt, kraken_contam_kmer_fraction_scaffold)
  
 numeric_cols <- names(contam_scaff_taxonomy)[sapply(contam_scaff_taxonomy, is.numeric)]
 contam_scaff_taxonomy[, (numeric_cols) := lapply(.SD, function(x) fifelse(is.na(x), 0, x)), .SDcols = numeric_cols]
 contam_scaff_taxonomy <- merge(contam_scaff_taxonomy, scaffold_info_dt[, .(scaffold_name, prot_taxonomy_score, nucl_taxonomy_score)], by = "scaffold_name")
 
 for (tax_rank in c("domain", "phylum", "class", "order", "family", "genus")) {
  dia_hit_set_cov_fraction_col <- paste("dia_hit_set_cov_fraction", tax_rank, sep = "_")
  rel_dia_hit_set_cov_fraction_col <- paste("rel_dia_hit_set_cov_fraction", tax_rank, sep = "_")
  tax_rank_prot_col <- paste(tax_rank, "prot", sep = "_")
  kraken_kmer_fraction_col <- paste("kraken_kmer_fraction", tax_rank, sep = "_")
  kraken_classif_kmer_fraction_col <- paste("kraken_classif_kmer_fraction", tax_rank, sep = "_")
  tax_rank_nucl_col <- paste(tax_rank, "nucl", sep = "_")
  
  contam_scaff_taxonomy[nucl_taxonomy_score == 1 & get(kraken_kmer_fraction_col) > KRAKEN_KMER_FRACTION_TAX_ASSIGN_THR & get(kraken_classif_kmer_fraction_col) > KRAKEN_CLASSIF_KMER_FRACTION_TAX_ASSIGN_THR, "assigned_tax" := get(tax_rank_nucl_col)]
  contam_scaff_taxonomy[nucl_taxonomy_score == 0 & get(dia_hit_set_cov_fraction_col) > DIA_COV_FRACTION_TAX_ASSIGN_THR & get(rel_dia_hit_set_cov_fraction_col) > DIA_REL_COV_FRACTION_TAX_ASSIGN_THR, "assigned_tax" := get(tax_rank_prot_col)]
  contam_scaff_taxonomy[nucl_taxonomy_score == 1 & prot_taxonomy_score == 1 & get(kraken_kmer_fraction_col) <= KRAKEN_KMER_FRACTION_TAX_ASSIGN_THR & get(kraken_classif_kmer_fraction_col) <= KRAKEN_CLASSIF_KMER_FRACTION_TAX_ASSIGN_THR & get(dia_hit_set_cov_fraction_col) > DIA_COV_FRACTION_TAX_ASSIGN_THR & get(rel_dia_hit_set_cov_fraction_col) > DIA_REL_COV_FRACTION_TAX_ASSIGN_THR, "assigned_tax" := get(tax_rank_prot_col)]
 }
  
 contam_scaff_taxonomy[, ':=' ("prot_taxonomy_score" = NULL, "nucl_taxonomy_score" = NULL)]
 write.table(contam_scaff_taxonomy, file = "contaminant_scaffolds_taxonomy.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
 contamination_sequences <- genomic_scaffolds[names(genomic_scaffolds) %in% contam_scaffolds]
 writeXStringSet(contamination_sequences, file = "contaminant_sequences.fa", format = "fasta")
}

#Writing main results to files
write.table(scaffold_info_dt, file = "scaffold_contamination_info.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
decontaminated_genome <- genomic_scaffolds[names(genomic_scaffolds) %in% contam_scaffolds == F]
writeXStringSet(decontaminated_genome, file = "decontaminated_genome.fa", format = "fasta")
