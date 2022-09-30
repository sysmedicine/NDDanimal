########################################################################## FILE REQUIREMENTS ####
#
# [1] combineCounts/Data/metadata.xlsx
#     Curated metadata file.
#
##################################################################################### INPUTS ####
#
#
# [1] cts_ROSMAP_bulk_brain_RNA_seq/Results/
# [2] cts_ROSMAP_RNA_array/Results/
# [3] cts_ROSMAP_microglia_RNA_seq/Results/
# [4] cts_ROSMAP_microglia_single_cell_RNA_seq/Results/Expression/
# [5] cts_ROSMAP_single_nucleus_RNA_seq/Results/Expression/
# [6] cts_ROSMAP_blood_RNAseq/Results/
#      Expression counts in tsv format.
#
#################################################################################### OUTPUTS ####
#
# [1] combineCounts/Results/combineCounts.tsv
# [2] combineCounts/Results/combineCounts_entrez.tsv
#     Combined normalised counts files in tsv format.
#
# [3] vis/Data/combineCounts_impute.tsv
#     Combined normalised imputed counts in tsv format.
#
############################################################################# INITIALISATION ####
require(dplyr)
require(data.table)
require(edgeR)
require(foreach)
require(limma)
require(biomaRt)
require(tibble)
require(readxl)
require(Rmagic)
require(MetabolAnalyze)
require(Rtsne)
require(ggplot2)

## Define sources

userIn <- data.frame(row.names = c("ROSMAP_bulk_brain_RNA_seq",
                                   "ROSMAP_RNA_array",
                                   "ROSMAP_microglia_RNA_seq",
                                   "ROSMAP_microglia_single_cell_RNA_seq",
                                   "ROSMAP_single_nucleus_RNA_seq",
                                   "ROSMAP_blood_RNAseq"),
                     dirs = c("cts_ROSMAP_bulk_brain_RNA_seq/Results/",
                              "cts_ROSMAP_RNA_array/Results/",
                              "cts_ROSMAP_microglia_RNA_seq/Results/",
                              "cts_ROSMAP_microglia_single_cell_RNA_seq/Results/Expression/",
                              "cts_ROSMAP_single_nucleus_RNA_seq/Results/Expression/",
                              "cts_ROSMAP_blood_RNAseq/Results/"),
                     genesCol = c(1, 1, 1, 1, 1, 1),
                     countsCol = c(2, 2, 2, 2, 2, 2, 2))

userOut <- "combineCounts/Results/"
userOut_vis <- "vis/Data/combineCounts_impute.tsv"

####################################################################################### MAIN ####

inExcel <- "combineCounts/Data/metadata.xlsx"                                                # Prototype code for loading colData from the master counts excel book

colData <- read_excel(inExcel, sheet = 1)


## Counts matrix

## Read individual counts matrices
cts <- foreach (f = 1:nrow(userIn),                                                                         # For each source
                .final = function(f) setNames(f, rownames(userIn))
) %do% {
    foreach (g = list.files(path = as.character(userIn$dirs[f]),                                            # For each sample
                            full.names = T,
                            pattern = ".tsv"),
             .final = function(g) setNames(g,
                                           list.files(path = as.character(userIn$dirs[f]),
                                                      pattern = ".tsv")
             )
    ) %do% {
        c <- fread(g) %>%
            dplyr::select(gene = userIn$genesCol[f],                                                        # Read in expression and gene symbol data
                   counts = userIn$countsCol[f])
        c <- c[which(!(duplicated(c$gene))), ]                                                              # Remove duplicates
        names(c) <- c("gene",
                      paste0(rownames(userIn)[f],".",sub(x = g, pattern = "^.+/", replacement = "")))
        c
    }
} %>%
    unlist(recursive = F)

batch <- foreach (f = 1:nrow(userIn),                                                                       # store batch data
                    .final = function(f) setNames(f, rownames(userIn))
) %do% {
    foreach (g = list.files(path = as.character(userIn$dirs[f]),
                            full.names = T,
                            pattern = ".tsv"),
             .final = function(g) setNames(g,
                                           list.files(path = as.character(userIn$dirs[f]),
                                                      pattern = ".tsv")
             )
    ) %do% {
        f
    }
} %>% unlist()

## List all genes in all counts matrices and prepare to merge counts matrices
genes <- foreach (f = names(cts)) %do% {                                                                    # Get list of all genes
    cts[[f]]$gene
} %>% unlist() %>%
    unique()

cts <- foreach (f = names(cts),
                .final = function(f) setNames(f, names(cts))
) %do% {                                                                                                    # Add new rows for NA genes
    absentGenesIdx <- which(!genes %in% cts[[f]]$gene)
    newRows <- data.frame(gene = genes[absentGenesIdx],
                          counts = NA)
    rbind(unique(cts[[f]]),
          newRows,
          use.names = F)
}

## Merge counts matrices
for (f in 1:length(names(cts))) {                                                                           # Merge all samples into a single counts matrix
    if (f == 1){
        cts3 <- cts[[f]]                                                                                    # cts3 is counts matrix with ensembl ids
    }
    else {
        cts3 <- cts3 %>% full_join(cts[[f]],
                                   by = "gene")
    }
}

rownames(cts3) <- cts3$gene                                                                                 # Clean cts3 genes column to rownames
cts3$gene <- NULL

genes_idx_all_NA <- which(is.na(rowMeans(cts3, na.rm = T)))                                                 # Remove genes with all NA
if (length(genes_idx_all_NA) != 0) {
    cts3 <- cts3[-genes_idx_all_NA,]
    genes <- genes[-genes_idx_all_NA]
}

genes_idx_all_zero <- which(rowMeans(cts3, na.rm = T) == 0) %>% unlist() %>% as.numeric()                   # Remove genes with all zero
if (length(genes_idx_all_zero) != 0) {
    cts3 <- cts3[-genes_idx_all_zero,]
    genes <- genes[-genes_idx_all_zero]
}

## Normalisation

## Scale counts to 1 million per sample

scaleMillion <- function(x) {                                                                               # Function to scale counts so that each sample sums to 1 million counts
    scaleFactor <- 10^6 / sum(x, na.rm = T)
    x * scaleFactor
}

cts3 <- foreach (f = cts3,
                 .final = function(f) setNames(f, names(cts3))) %do% {                                      # Do
                    scaleMillion(f)
                } %>% as.data.frame()
rownames(cts3) <- genes

## Missing data imputation

cts3_MAGIC <- t(cts3)                                                                                       # Make a copy of cts3. Impute to remove NAs
cts3_MAGIC[is.na(cts3_MAGIC)] <- 0
cts3_MAGIC <- magic(cts3_MAGIC, npca = 20)
cts3_MAGIC <- t(as.matrix(cts3_MAGIC$result))
cts3_MAGIC <- as_tibble(cts3_MAGIC)
rownames(cts3_MAGIC) <- genes

## TMM normalise between samples within sources: normalisation factors

cts3_MAGIC_DGEList <- DGEList(cts3_MAGIC, group = colData$Source)                                           # On the imputed data, calculate TMM normalisation factors
TMM_norm_factors <- calcNormFactors(cts3_MAGIC_DGEList, method = "TMM")
TMM_norm_factors <- TMM_norm_factors$samples$norm.factors

normalise_TMM_Pareto_batch <- function(x) {

    ## TMM normalise

    cts3 <- foreach(i = 1:length(x),                                                                         # Apply TMM normalisation factors on the original (unimputed) data
                   .final = function(f) setNames(f, names(x))) %do% {
                       x[[i]] * TMM_norm_factors[i]
                   } %>% as.data.frame()
    rownames(x) <- genes

    ## Pareto scale each gene
    cts3_genes_stdevs <- apply(cts3, 1, sd, na.rm = T)
    cts3 <- sweep(cts3, 1, sqrt(cts3_genes_stdevs), "/")

    ## Remove batch effects between sources
    x <- removeBatchEffect(x = x,                                                                         # Remove batch effects (limma)
                           batch = batch) %>%
        as.data.frame()

    return(x)
}

cts3 <- normalise_TMM_Pareto_batch(cts3)
rownames(cts3) <- genes

## IO: Data exploration/vis to file

cts3_MAGIC_vis <- normalise_TMM_Pareto_batch(cts3_MAGIC)
cts3_MAGIC_vis <- t(cts3_MAGIC_vis)
fwrite(x = cts3_MAGIC_vis, file = userOut_vis, sep = "\t")


## Output

## Save ensembl counts matrix to disk
fwrite(x = cts3,                                                                                            # Write master counts files
       file = paste0(userOut, "combineCounts.tsv"),
       sep = "\t",
       row.names = T)
