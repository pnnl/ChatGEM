rm(list = ls())

library(magrittr)
library(ggplot2)
library(pmartR)
library(here)
library(stringr)
library(dplyr)
library(tidyr)
library(patchwork)
library(tibble)

here::i_am("Code/OPAL_Pp64_iBAQ_processing_2026_05_11_script.R")

# the first part of this script is built from Kelly Stratton's
# code for processing this data
# however, unlike that script we do not run reference rollup and statistical
# analyses at the protein level and instead focuse on
# running analyses for iBAQ processing

## Import/Processing -------------------------------------------------------

pep_edata <- readxl::read_xlsx(path = here("Data", "EMSL_62080_Opal_Pp_B2.xlsx"),
                              sheet = "Precursor")
pep_edata <- pep_edata %>%
  dplyr::rename(ProteinID = `Protein.Ids`,
                Modified_Sequence = `Modified.Sequence`,
                PrecursorID = `Precursor.Id`)
any(duplicated(pep_edata$PrecursorID)) # FALSE

# create e_meta and finalize e_data
pep_emeta <- pep_edata[, 1:3]
pep_edata <- pep_edata[, -c(1:2)] # PrecursorID is edata_cname

geneIDs <- read.csv(here("Data", "Pseudomonas_putida_KT2440_UniProt_2017-04-07_wAnnotation.csv"), header = TRUE)
geneIDs <- geneIDs %>%
  dplyr::rename(ProteinID = `Uniprot._ID`,
                Entry_Name = `Entry.Name`,
                Protein_Names = `Protein.names`,
                Gene_Names = `Gene.Names`,
                EC_number = `EC.number`,
                Gene_Ontology_biological_process = `Gene.Ontology..biological.process.`,
                Gene_Ontology_cellular_component = `Gene.Ontology..cellular.component.`,
                Gene_Ontology_GO = `Gene.Ontology..GO.`,
                Gene_Ontology_molecular_function = `Gene.Ontology..molecular.function.`,
                Gene_Ontology_IDs = `Gene.Ontology.IDs`)
# pep_emeta <- merge(x = pep_emeta, y = geneIDs, by = "ProteinID", all.x = TRUE)

# create f_data
pep_fdata <- readxl::read_xlsx(path = here("Data", "EMSL_62080_Opal_Pp_B2.xlsx"),
                               sheet = "Info")

pep_fdata <- pep_fdata[, c(1:5)]
pep_fdata <- pep_fdata %>%
  dplyr::rename(SampleID = Sample)

pep_fdata$Strain_short <- pep_fdata$Strain
pep_fdata$Strain_short <- str_replace_all(pep_fdata$Strain_short, "PNNL ", "")
pep_fdata$Strain_short <- str_replace_all(pep_fdata$Strain_short, "Sage landing pads", "Sage LPs")
pep_fdata$Strain_short <- str_replace_all(pep_fdata$Strain_short, "inducible, crystal violet", "ind, crys viol")
pep_fdata$Strain_short <- str_replace_all(pep_fdata$Strain_short, "constitutive", "const")

pep_fdata$Strain_Carbon_Time <- paste(pep_fdata$Strain_short, pep_fdata$Carbon_Source, pep_fdata$Time_Point, sep = "_")
pep_fdata$Carbon_Time <- paste(pep_fdata$Carbon_Source, pep_fdata$Time_Point, sep = "_")

table(pep_fdata$Strain_Carbon_Time)

all(pep_fdata$SampleID %in% names(pep_edata)[-1])
all(names(pep_edata)[-1] %in% pep_fdata$SampleID)
any(duplicated(pep_fdata$SampleID)) # FALSE

# create pmartR-compatible object
my_pepdata <- as.pepData(e_data = pep_edata,
                         f_data = pep_fdata,
                         e_meta = pep_emeta,
                         edata_cname = "PrecursorID",
                         fdata_cname = "SampleID",
                         emeta_cname = "ProteinID",
                         data_scale = "abundance")

# table of strain by inducer
table(my_pepdata$f_data$Strain_short, my_pepdata$f_data$Carbon_Source, useNA = "ifany")

# Transform data to the log2 scale
my_pepdata <- edata_transform(my_pepdata, data_scale = "log2")

my_pepdata <- group_designation(my_pepdata, main_effects = c("Strain_Carbon_Time"))

pep_gdf <- attr(my_pepdata, "group_DF") %>%
  dplyr::select(-SampleID) %>% # Remove replicate identifier; only retain design information
  dplyr::distinct()

pep_info <- get_data_info(my_pepdata) # to be used in the report

# Comparisons of interest are:
# 1. Compare Strain for each condition (time point + carbon source) --> 24 comparisons
#       a. KT 2440_Glucose_12 hr vs KT 3 x Sage LPs_Glucose_12 hr
#       b. KT 2440_Glucose_12 hr vs KT 3 x Sage LPs (ind, crys viol)_Glucose_12 hr
#       c. KT 2440_Glucose_12 hr vs KT 3 x Sage LPs, dCas12a (const)_Glucose_12 hr
#       d. KT 3 x Sage LPs_Glucose_12 hr vs KT 3 x Sage LPs (ind, crys viol)_Glucose_12 hr
#       e. KT 3 x Sage LPs_Glucose_12 hr vs KT 3 x Sage LPs, dCas12a (const)_Glucose_12 hr
#       f. KT 3 x Sage LPs (ind, crys viol)_Glucose_12 hr vs KT 3 x Sage LPs, dCas12a (const)_Glucose_12 hr
comps1_set1 <- data.frame(Control = c("KT 2440_Glucose_12 hr", "KT 2440_Glucose_12 hr",
                                      "KT 2440_Glucose_12 hr", "KT 3 x Sage LPs_Glucose_12 hr",
                                      "KT 3 x Sage LPs_Glucose_12 hr", "KT 3 x Sage LPs (ind, crys viol)_Glucose_12 hr"),
                          Test = c("KT 3 x Sage LPs_Glucose_12 hr", "KT 3 x Sage LPs (ind, crys viol)_Glucose_12 hr",
                                   "KT 3 x Sage LPs, dCas12a (const)_Glucose_12 hr", "KT 3 x Sage LPs (ind, crys viol)_Glucose_12 hr",
                                   "KT 3 x Sage LPs, dCas12a (const)_Glucose_12 hr", "KT 3 x Sage LPs, dCas12a (const)_Glucose_12 hr"))

#       g. KT 2440_Gluconate_12 hr vs KT 3 x Sage LPs_Gluconate_12 hr
#       h. KT 2440_Gluconate_12 hr vs KT 3 x Sage LPs (ind, crys viol)_Gluconate_12 hr
#       i. KT 2440_Gluconate_12 hr vs KT 3 x Sage LPs, dCas12a (const)_Gluconate_12 hr
#       j. KT 3 x Sage LPs_Gluconate_12 hr vs KT 3_Gluconate_12 x Sage LPs (ind, crys viol)_Gluconate_12 hr
#       k. KT 3 x Sage LPs_Gluconate_12 hr vs KT 3_Gluconate_12 x Sage LPs, dCas12a (const)_Gluconate_12 hr
#       l. KT 3 x Sage LPs (ind, crys viol)_Gluconate_12 hr vs KT 3 x Sage LPs, dCas12a (const)_Gluconate_12 hr
comps1_set2 <- data.frame(Control = c("KT 2440_Gluconate_12 hr", "KT 2440_Gluconate_12 hr",
                                      "KT 2440_Gluconate_12 hr", "KT 3 x Sage LPs_Gluconate_12 hr",
                                      "KT 3 x Sage LPs_Gluconate_12 hr", "KT 3 x Sage LPs (ind, crys viol)_Gluconate_12 hr"),
                          Test = c("KT 3 x Sage LPs_Gluconate_12 hr", "KT 3 x Sage LPs (ind, crys viol)_Gluconate_12 hr",
                                   "KT 3 x Sage LPs, dCas12a (const)_Gluconate_12 hr", "KT 3 x Sage LPs (ind, crys viol)_Gluconate_12 hr",
                                   "KT 3 x Sage LPs, dCas12a (const)_Gluconate_12 hr", "KT 3 x Sage LPs, dCas12a (const)_Gluconate_12 hr"))

#       m. KT 2440_Glucose_24 hr vs KT 3 x Sage LPs_Glucose_24 hr
#       n. KT 2440_Glucose_24 hr vs KT 3 x Sage LPs (ind, crys viol)_Glucose_24 hr
#       o. KT 2440_Glucose_24 hr vs KT 3 x Sage LPs, dCas12a (const)_Glucose_24 hr
#       p. KT 3 x Sage LPs hr vs KT 3 x Sage LPs (ind, crys viol)_Glucose_24 hr
#       q. KT 3 x Sage LPs hr vs KT 3 x Sage LPs, dCas12a (const)_Glucose_24 hr
#       r. KT 3 x Sage LPs (ind, crys viol)_Glucose_24 hr vs KT 3 x Sage LPs, dCas12a (const)_Glucose_24 hr
comps1_set3 <- data.frame(Control = c("KT 2440_Glucose_24 hr", "KT 2440_Glucose_24 hr",
                                      "KT 2440_Glucose_24 hr", "KT 3 x Sage LPs_Glucose_24 hr",
                                      "KT 3 x Sage LPs_Glucose_24 hr", "KT 3 x Sage LPs (ind, crys viol)_Glucose_24 hr"),
                          Test = c("KT 3 x Sage LPs_Glucose_24 hr", "KT 3 x Sage LPs (ind, crys viol)_Glucose_24 hr",
                                   "KT 3 x Sage LPs, dCas12a (const)_Glucose_24 hr", "KT 3 x Sage LPs (ind, crys viol)_Glucose_24 hr",
                                   "KT 3 x Sage LPs, dCas12a (const)_Glucose_24 hr", "KT 3 x Sage LPs, dCas12a (const)_Glucose_24 hr"))

#       s. KT 2440_Gluconate_24 hr vs KT 3 x Sage LPs_Gluconate_24 hr
#       t. KT 2440_Gluconate_24 hr vs KT 3 x Sage LPs (ind, crys viol)_Glucose_24 hr
#       u. KT 2440_Gluconate_24 hr vs KT 3 x Sage LPs, dCas12a (const)_Glucose_24 hr
#       v. KT 3 x Sage LPs_Gluconate_24 hr vs KT 3_Gluconate_12 x Sage LPs (ind, crys viol)_Gluconate_24 hr
#       w. KT 3 x Sage LPs_Gluconate_24 hr vs KT 3_Gluconate_12 x Sage LPs, dCas12a (const)_Gluconate_24 hr
#       x. KT 3 x Sage LPs (ind, crys viol)_Gluconate_24 hr vs KT 3 x Sage LPs, dCas12a (const)_Gluconate_24 hr
comps1_set4 <- data.frame(Control = c("KT 2440_Gluconate_24 hr", "KT 2440_Gluconate_24 hr",
                                      "KT 2440_Gluconate_24 hr", "KT 3 x Sage LPs_Gluconate_24 hr",
                                      "KT 3 x Sage LPs_Gluconate_24 hr", "KT 3 x Sage LPs (ind, crys viol)_Gluconate_24 hr"),
                          Test = c("KT 3 x Sage LPs_Gluconate_24 hr", "KT 3 x Sage LPs (ind, crys viol)_Gluconate_24 hr",
                                   "KT 3 x Sage LPs, dCas12a (const)_Gluconate_24 hr", "KT 3 x Sage LPs (ind, crys viol)_Gluconate_24 hr",
                                   "KT 3 x Sage LPs, dCas12a (const)_Gluconate_24 hr", "KT 3 x Sage LPs, dCas12a (const)_Gluconate_24 hr"))

# 2. Compare Time_Point for each condition (strain + carbon source) --> 8 comparisons
#       a. KT 2440_Glucose_24 hr vs KT 2440_Glucose_12 hr
#       b. KT 2440_Gluconate_24 hr vs KT 2440_Gluconate_12 hr
#       c. KT 3 x Sage LPs_Glucose_24 hr vs KT 3 x Sage LPs_Glucose_12 hr
#       d. KT 3 x Sage LPs_Gluconate_24 hr vs KT 3 x Sage LPs_Gluconate_12 hr
#       e. KT 3 x Sage LPs (ind, crys viol)_Glucose_24 hr vs KT 3 x Sage LPs (ind, crys viol)_Glucose_12 hr
#       f. KT 3 x Sage LPs (ind, crys viol)_Gluconate_24 hr vs KT 3 x Sage LPs (ind, crys viol)_Gluconate_12 hr
#       g. KT 3 x Sage LPs, dCas12a (const)_Glucose_24 hr vs KT 3 x Sage LPs, dCas12a (const)_Glucose_12 hr
#       h. KT 3 x Sage LPs, dCas12a (const)_Gluconate_24 hr vs KT 3 x Sage LPs, dCas12a (const)_Gluconate_12 hr
comps2 <- data.frame(Control = c("KT 2440_Glucose_12 hr", "KT 2440_Gluconate_12 hr",
                                 "KT 3 x Sage LPs_Glucose_12 hr", "KT 3 x Sage LPs_Gluconate_12 hr",
                                 "KT 3 x Sage LPs (ind, crys viol)_Glucose_12 hr", "KT 3 x Sage LPs (ind, crys viol)_Gluconate_12 hr",
                                 "KT 3 x Sage LPs, dCas12a (const)_Glucose_12 hr", "KT 3 x Sage LPs, dCas12a (const)_Gluconate_12 hr"),
                          Test = c("KT 2440_Glucose_24 hr", "KT 2440_Gluconate_24 hr",
                                   "KT 3 x Sage LPs_Glucose_24 hr", "KT 3 x Sage LPs_Gluconate_24 hr",
                                   "KT 3 x Sage LPs (ind, crys viol)_Glucose_24 hr", "KT 3 x Sage LPs (ind, crys viol)_Gluconate_24 hr",
                                   "KT 3 x Sage LPs, dCas12a (const)_Glucose_24 hr", "KT 3 x Sage LPs, dCas12a (const)_Gluconate_24 hr"))

# 3. Compare Carbon_Source for each condition (strain + time point) --> 8 comparisons
#       a. KT 2440_Gluconate_12 hr vs KT 2440_Glucose_12 hr
#       b. KT 2440_Gluconate_24 hr vs KT 2440_Glucose_24 hr
#       c. KT 3 x Sage LPs_Gluconate_12 hr vs KT 3 x Sage LPs_Glucose_12 hr
#       d. KT 3 x Sage LPs_Gluconate_24 hr vs KT 3 x Sage LPs_Glucose_24 hr
#       e. KT 3 x Sage LPs (ind, crys viol)_Gluconate_12 hr vs KT 3 x Sage LPs (ind, crys viol)_Glucose_12 hr
#       f. KT 3 x Sage LPs (ind, crys viol)_Gluconate_24 hr vs KT 3 x Sage LPs (ind, crys viol)_Glucose_24 hr
#       g. KT 3 x Sage LPs, dCas12a (const)_Gluconate_12 hr vs KT 3 x Sage LPs, dCas12a (const)_Glucose_12 hr
#       h. KT 3 x Sage LPs, dCas12a (const)_Gluconate_24 hr vs KT 3 x Sage LPs, dCas12a (const)_Glucose_24 hr
comps3 <- data.frame(Control = c("KT 2440_Glucose_12 hr", "KT 2440_Glucose_24 hr",
                                 "KT 3 x Sage LPs_Glucose_12 hr", "KT 3 x Sage LPs_Glucose_24 hr",
                                 "KT 3 x Sage LPs (ind, crys viol)_Glucose_12 hr", "KT 3 x Sage LPs (ind, crys viol)_Glucose_24 hr",
                                 "KT 3 x Sage LPs, dCas12a (const)_Glucose_12 hr", "KT 3 x Sage LPs, dCas12a (const)_Glucose_24 hr"),
                     Test = c("KT 2440_Gluconate_12 hr", "KT 2440_Gluconate_24 hr",
                              "KT 3 x Sage LPs_Gluconate_12 hr", "KT 3 x Sage LPs_Gluconate_24 hr",
                              "KT 3 x Sage LPs (ind, crys viol)_Gluconate_12 hr", "KT 3 x Sage LPs (ind, crys viol)_Gluconate_24 hr",
                              "KT 3 x Sage LPs, dCas12a (const)_Gluconate_12 hr", "KT 3 x Sage LPs, dCas12a (const)_Gluconate_24 hr"))


pep_comps <- data.frame(Control = c(comps1_set1$Control, comps1_set2$Control, comps1_set3$Control, comps1_set4$Control,
                                    comps2$Control, comps3$Control),
                        Test = c(comps1_set1$Test, comps1_set2$Test, comps1_set3$Test, comps1_set4$Test,
                                 comps2$Test, comps3$Test))
pep_comps <- pep_comps %>% dplyr::distinct() # Ensure that there are no redundant comparisons


# Define imdanova filter
gp_imdanova_filt <- imdanova_filter(my_pepdata)
summary(gp_imdanova_filt,
        omicsData = my_pepdata,
        comparisons = pep_comps,
        min_nonmiss_anova = 2,
        min_nonmiss_gtest = 3,
        remove_singleton_groups = TRUE)
# Apply imdanova filter based on specified comparisons
my_pepdata <- applyFilt(gp_imdanova_filt,
                        omicsData = my_pepdata,
                        comparisons = pep_comps,
                        min_nonmiss_anova = 2,
                        min_nonmiss_gtest = 3,
                        remove_singleton_groups = TRUE)

# -------------------------------------------------------------------------

## QC / Outlier Detection --------------------------------------------------

# molecule filter - not something Javi was doing, but he does do the IMD-ANOVA filter above which covers this
# molfilt <- molecule_filter(my_pepdata)
# summary(molfilt, min_num = 2)
# # Minimum Number:2
# # Filtered:641
# # Not Filtered:53306
# my_pepdata <- applyFilt(molfilt, my_pepdata, min_num = 2)


# proteomics filtering - require at least 2 peps per prot - not something Javi was doing
# protfilt <- proteomics_filter(my_pepdata)
# summary(protfilt, min_num_peps = 2)


# Apply rmd_filter
pep_rmd_filt <- rmd_filter(omicsData = my_pepdata,
                            ignore_singleton_groups = TRUE,
                            metrics = c("MAD", "Skewness", "Correlation", "Proportion_Missing"))
pep_cand_out <- pep_rmd_filt$SampleID[which(pep_rmd_filt$pvalue < 0.0001)]

pep_corr_res <- cor_result(my_pepdata)
pep_corr_res_mat <- reshape2::melt(pep_corr_res) %>%
  dplyr::rename(Correlation = value)
outcol <- ifelse(unique(pep_corr_res_mat$Var1) %in% pep_cand_out,"red", "black")

set.seed(1)
pep_pca_res <- dim_reduction(my_pepdata, k = 2)
temp_pcs_df <- data.frame(SampleID = pep_pca_res$SampleID,
                          PC1 = pep_pca_res$PC1,
                          PC2 = pep_pca_res$PC2)


# I'm curious to look at missingness for these data
## Missingness -------------------------------------------------
# this was added in this "with_gtest" version, to see what the patters of
# missing data are, how much, etc.
miss_filt <- missingval_result(omicsData = my_pepdata)
pep_missing_bar <- plot(miss_filt, omicsData = my_pepdata, plot_type = "bar", color_by = "Strain")
pep_missing_scatter <- plot(miss_filt, omicsData = my_pepdata, plot_type = "scatter")

# 2 outliers flagged/observed:
# "OPAL_02_Gluconate_24_h_B_147"
# "OPAL_02_Gluconate_24_h_E_157"
### Leaning toward removing both - on the basis of missing data, correlation, and standing out more from all the others in the PCA plot
pep_out <- pep_cand_out[c(grep("_147", pep_cand_out), grep("_157", pep_cand_out))]
myfilt <- custom_filter(my_pepdata, f_data_remove = pep_out)
my_pepdata <- applyFilt(myfilt, my_pepdata)

# -------------------------------------------------------------------------

## Normalization ----------------------------------------------------------

# Test median centering
my_pepdata_norm_obj <- normalize_global(omicsData = my_pepdata,
                                    subset_fn = "all",
                                    norm_fn = "median",
                                    apply_norm = FALSE)
normRes_tests(my_pepdata_norm_obj)
# $p_location
# 0.4890403

# Normalize data using median centering
my_pepdata_norm <- normalize_global(omicsData = my_pepdata,
                                    subset_fn = "all",
                                    norm_fn = "median",
                                    apply_norm = TRUE,
                                    backtransform = TRUE)

pep_corr_res2 <- cor_result(my_pepdata_norm)
pep_corr_res_mat2 <- reshape2::melt(pep_corr_res2) %>%
  dplyr::rename(Correlation = value)

set.seed(1)
pep_pca_res2 <- dim_reduction(my_pepdata_norm, k = 2)
temp_pcs_df2 <- data.frame(SampleID = pep_pca_res2$SampleID,
                           PC1 = pep_pca_res2$PC1,
                           PC2 = pep_pca_res2$PC2)


## iBAQ Processing ---------------------------------------------------------
digest_fasta2 <- function(fasta_df) {
   if (missing(fasta_df) || is.null(fasta_df)) stop("ERROR: Need fasta_df")
   if (!all(colnames(fasta_df) == c("Protein.Group", "Sequence"))) {
      stop("ERROR: FASTA format not recognised. Expect columns: Protein.Group, Sequence")
   }

   fasta_dt <- data.table::data.table(fasta_df)
   fasta_dt$Sequence <- gsub("((?<=[KR])(?=[^P]))", ",", fasta_dt$Sequence, perl = TRUE)
   digest_dt <- splitstackshape::cSplit(fasta_dt, "Sequence", sep = ",", direction = "long", drop = FALSE, fixed = FALSE)
   digest_dt <- digest_dt[nchar(Sequence) > 6 & nchar(Sequence) < 31, ]
   return(as.data.frame(digest_dt))
}



# instead of standard rollup like with other analyse
# we are going to run iBAQ calculations instead

# before rollup, calculate the number of unique peptides per protein and add that to emeta,
# then make sure to use the "emeta_cols" argument to protein_quant() so that column is retained
my_pepdata_norm$e_meta <- my_pepdata_norm$e_meta %>%
   group_by(ProteinID) %>%
   mutate(unique_PrecursorID_per_pro = n_distinct(PrecursorID),
          unique_ModifiedSequence_per_pro = n_distinct(Modified_Sequence)) %>%
   ungroup()

# pull out the edata information
precursor_df <- my_pepdata_norm$e_data %>%
   # want to add in protein information as well
   dplyr::left_join(my_pepdata_norm$e_meta %>% dplyr::select(ProteinID, Modified_Sequence, PrecursorID), by = "PrecursorID") %>%
   dplyr::relocate(ProteinID)

# create a list for ibaq setup with the following pieces of information:
# ProteinID, sample_list, quant, protein_list, id, and stripped.sequence
ibaq_rollup_list <- precursor_df %>%
   select(ProteinID, PrecursorID, Modified_Sequence,contains("OPAL")) %>%
   pivot_longer(-c(ProteinID, PrecursorID,Modified_Sequence),
                names_to = "sample_id", values_to = "intensity") %>%
   filter(!is.na(intensity)) %>%
   mutate(protein_group_seq = paste(ProteinID, PrecursorID, sep = "_")) %>%
   arrange(sample_id) %>%
   rename(protein_list = protein_group_seq,
          sample_list = sample_id,
          id = PrecursorID,
          Stripped.Sequence = Modified_Sequence,
          quant = intensity) %>%
   as.list()

# run fast_MaxLFQ and get matrix result
ibaq_rollup_mat <- iq::fast_MaxLFQ(ibaq_rollup_list)[[1]]

# convert to tibble with Protein.Group and Stripped.Sequence columns
ibaq_rollup <- as.data.frame(2^ibaq_rollup_mat) %>%    # back-transform to linear intensity
   # rename row names to protein group seq (which is just ProteinID for us in this use case)
   rownames_to_column(var = "protein_group_seq") %>%
   # this is not super relevant for us as stripped sequence was NA anyways for us
   mutate(Protein.Group = str_extract(protein_group_seq, "^[^_]+"),
          Stripped.Sequence = str_extract(protein_group_seq, "(?<=_).*")) %>%
   select(-protein_group_seq)

# include_ibaq2 <- my_pepdata_norm$e_meta %>%
#    # pull out number of peptides that map to a protein
#    dplyr::select(ProteinID,unique_PrecursorID_per_pro) %>%
#    dplyr::rename(Protein.Group = ProteinID,
#                  n = unique_PrecursorID_per_pro) %>%
#    dplyr::distinct()

# digest and count peptides per protein (n)
fasta_path = here("Data","ID_008893_17C184F0.fasta")
seq_fasta <- importFasta(fasta_path) %>%
   rename(Protein.Group = accession) %>%
   filter(!is.na(Protein.Group)) %>%
   select( -name) %>%
   rename(Sequence = sequence) %>%
   select(Protein.Group, Sequence)

# molecular weight data
pmw_df_init <- seq_fasta %>%
   dplyr::mutate(Sequence = toupper(Sequence),
                 MW = sapply(Sequence, function(seq){
                    seqinr::pmw(seqaa = unlist(stringr::str_split(seq, "")))
                 }))
pmw_df_mid <- pmw_df_init %>%
   dplyr::filter(!stringr::str_detect(Protein.Group,"\\|"))

pmw_df <- pmw_df_init %>%
   dplyr::mutate(Protein.Group = stringr::str_split_i(Protein.Group,"\\|",2)) %>%
   dplyr::filter(!is.na(Protein.Group)) %>%
   dplyr::bind_rows(pmw_df_mid)

include_ibaq_init <- digest_fasta2(seq_fasta) %>%
   # filter(!Protein.Group %in% Histones) %>%
   dplyr::group_by(Protein.Group) %>% tally(name = "n")

include_ibaq_mid <- digest_fasta2(seq_fasta) %>%
   # filter(!Protein.Group %in% Histones) %>%
   dplyr::group_by(Protein.Group) %>% tally(name = "n") %>%
   dplyr::filter(!stringr::str_detect(Protein.Group,"\\|"))

include_ibaq <- include_ibaq_init %>%
   dplyr::mutate(Protein.Group = stringr::str_split_i(Protein.Group,"\\|",2)) %>%
   dplyr::filter(!is.na(Protein.Group)) %>%
   dplyr::bind_rows(include_ibaq_mid)

split_prots <- stringr::str_split(ibaq_rollup$Protein.Group,"\\;")

max_mw <- purrr::map_df(split_prots, function(prots){
   seq_sub <- pmw_df[pmw_df$Protein.Group %in% prots,]
   max_mw = max(seq_sub$MW, na.rm = T)
   leading_protein = seq_sub$Protein.Group[which(seq_sub$MW == max_mw)]

   # if more than one have the same molecular weight max, pick the first one
   leading_protein <- leading_protein[1]

   seq_sub <- include_ibaq[include_ibaq$Protein.Group %in% prots,]
   seq_sub <- seq_sub[seq_sub$Protein.Group == leading_protein,]
   num_tryptic <- seq_sub$n[1]

   return(data.frame(LeadingProtein = leading_protein,
                     MW = max_mw,
                     n = num_tryptic))
})

# data frame with molecular weight for each protein
mw_df <- ibaq_rollup
mw_df$MW <- max_mw$MW
mw_df$LeadingProtein <- max_mw$LeadingProtein
mw_df$n <- max_mw$n
mw_df_all <- mw_df %>% dplyr::select(Protein.Group, LeadingProtein, MW,n) %>% dplyr::distinct()
mw_df <- mw_df_all %>% dplyr::select(-n)

# data frame with molecular weight for each protein
num_tryptic_df <- ibaq_rollup
num_tryptic_df$n <- max_mw$n
include_ibaq <- num_tryptic_df %>% dplyr::select(Protein.Group, n) %>% dplyr::distinct()

# Compute iBAQ per dataset
ibaq_int <- ibaq_rollup %>%
   select(Protein.Group, Stripped.Sequence, contains("OPAL")) %>%
   # pivot longer so we have the following columns:
   # Protein.Group, Stripped.Sequence, dataset (which is our sample ids), and intensity
   pivot_longer(-c(Protein.Group, Stripped.Sequence), names_to = "dataset", values_to = "intensity") %>%
   # filter out missing intensity values
   filter(!is.na(intensity)) %>%
   # add in information about peptides mapping to protein numbers
   inner_join(include_ibaq, by = "Protein.Group") %>%
   # group by dataset/sample ID, protein.group, and number of peptides mapping to protein
   # technically grouping by n is likely redundant
   group_by(dataset, Protein.Group, n) %>%
   # take the sum of intensity for each of those groups
   # here this feels unnecessary as the values have already been consolidated to one value
   summarise(sum_int = sum(intensity), .groups = "drop") %>%
   # divide sum by number of peptides used to map to protein
   mutate(iBAQ = sum_int/n) %>%
   # filter out missing ibaq values
   #filter(!is.na(iBAQ)) %>%
   # pull out only distinct values
   distinct(dataset, iBAQ, Protein.Group)

# put back into wide data.frame to match edata setup
ibaq_df <- ibaq_int %>%
   pivot_wider(names_from = dataset, values_from = iBAQ) %>%
   arrange(desc(rowMeans(select(., 2:ncol(.)), na.rm = TRUE)))


ibaq_df <- ibaq_df %>%
   dplyr::left_join(mw_df_all, by = "Protein.Group") %>%
   dplyr::relocate(Protein.Group,LeadingProtein,MW,n) %>%
   dplyr::rename(num_tryptic = n)


# save files
# without gCDW conversion
write.csv(ibaq_df, file = here("Processed_Data","05122026_OPAL_iBAQ_leadingprotein_no_gCDW_conversion_putida.csv"),
          row.names = F)
