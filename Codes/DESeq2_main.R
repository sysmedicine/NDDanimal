########################################################################## FILE REQUIREMENTS ####
#
# [1] combineCounts/Data/metadata.xlsx
#     Curated metadata file.
#
##################################################################################### INPUTS ####
#
# [1] combineCounts/Results/combineCounts.tsv
#     Combined normalised counts file.
#
#################################################################################### OUTPUTS ####
#
# [1] main/Results/DESeq2/
#     DESeq2 results.
#
############################################################################# INITIALISATION ####
require(data.table)
require(DESeq2)
require(dplyr)
require(foreach)
require(ggplot2)
require(readxl)
require(scales)

inExcel <- "combineCounts/Data/metadata.xlsx"                                                # colData
inCounts <- "combineCounts/Results/combineCounts.tsv"                      # Counts matrix
user_alpha <- 1e-10                                                                                         # Set alpha for determining significant DEGs
user_l2fc_thresh <- 0                                                                                       # Set l2fc threshold for determining significant DEGs
outdir <- "main/Results/DESeq2/"                                                                   # Output dir for saved files

dir.create(outdir,
           showWarnings = F)

################################################################################## FUNCTIONS ####

see <- function(x) {                                                                                        # Helper function to tabulate DESeq results
    tibble(gene = x@rownames,
           baseMean = x$baseMean,
           l2fc = x$log2FoldChange,
           l2fc_se = x$lfcSE,
           stat = x$stat,
           pval = x$pvalue,
           padj = x$padj)
}

reverselog_trans <- function(base = exp(1)) {                                                               # Inverse log scale for ggplot2
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv,
              log_breaks(base = base),
              domain = c(1e-100, Inf))
}

####################################################################################### MAIN ####

# File I/O

colData <- read_excel(inExcel, sheet = 1)                                                                   # Read colData

counts <- fread(inCounts)                                                                                   # Read counts matrix
genes <- counts[[1]]                                                                                        # Record gene names
counts <- counts[,-1]                                                                                       # Remove gene column from actual matrix
rownames(counts) <- genes                                                                                   # Add back gene names as row headers

# Data selection

colData <- colData %>% filter(!is.na(Phenotype))
samples_to_select <- colData$Colname

counts <- counts %>% dplyr::select(all_of(samples_to_select))

# Data formatting

counts <- counts %>%
    as.matrix() %>%
    exp() %>%
    apply(MARGIN = c(1,2), FUN = function(x) as.integer(x))

excess_NA_threshold <- dim(counts)[2] * 0.4
genes_idx_excess_NA <- which(rowCounts(x = counts, value = NA) >= excess_NA_threshold)               # Remove genes with too many NA
if (length(genes_idx_excess_NA) != 0) {
    counts <- counts[-genes_idx_excess_NA,]
    genes <- genes[-genes_idx_excess_NA]
}

#genes_idx_all_NA <- which(is.na(rowMeans(counts_select, na.rm = T)))                                       # Remove genes with all NA
#if (length(genes_idx_all_NA) != 0) {
#    counts_select <- counts_select[-genes_idx_all_NA,]
#    genes <- genes[-genes_idx_all_NA]
#}

genes_idx_all_zero <- which(rowMeans(counts, na.rm = T) == 0) %>% unlist() %>% as.numeric()          # Remove genes with all zero
if (length(genes_idx_all_zero) != 0) {
    counts <- counts[-genes_idx_all_zero,]
    genes <- genes[-genes_idx_all_zero]
}

#counts_select <- sapply(X = counts_select, FUN = as.integer)                                               # Convert counts matrix to integer
#counts_select <- exp(counts_select)                                                                        # Adjust counts so that no values are negative
counts[is.na(counts)] <- 0                                                                    # NA to zero (DESeq assumes not detected)

rownames(counts) <- genes                                                                            # Restore gene names

# DESeq2

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ 0 + Phenotype)

sizeFactors(dds) <- 1                                                                                       # Tell DESeq that the counts are already normalised

dds <- DESeq(dds)
# cluster b (differential), cluster a (reference)

phenotypes <- colData$Phenotype %>% unique()

for (a in phenotypes) {
    for (b in phenotypes) {
        if (a == b) next
        else {
            res <- results(dds,
                           contrast = c("Phenotype", b, a))


            resData <- see(res)

            # Volcano plot

            resData <- resData %>%
                mutate(colour = case_when(abs(l2fc) < user_l2fc_thresh | padj > user_alpha ~ "grey",
                                          l2fc >= user_l2fc_thresh & padj <= user_alpha ~ "red",
                                          l2fc <= -user_l2fc_thresh & padj <= user_alpha ~ "blue")) %>%
                arrange(padj)

            ggplot(data = resData,
                   mapping = aes(x = l2fc, y = padj, colour = colour, alpha = colour)) +
                geom_point() +
                scale_alpha_manual(values = c(0.3, 0.1, 0.3, 0)) +
                scale_colour_manual(values = c("#0000FF", "#AAAAAA", "#FF0000", "#000000")) +
                scale_y_continuous(trans = reverselog_trans(10)) +
                #geom_vline(xintercept = -user_l2fc_thresh,
                #           col = "black",
                #           linetype = "dotted",
                #           size = 0.7) +
                #geom_vline(xintercept = user_l2fc_thresh,
                #           col = "black",
                #           linetype = "dotted",
                #           size = 0.7) +
                geom_hline(yintercept = user_alpha,
                           col = "black",
                           linetype = "dotted",
                           size = 0.7) +
                theme(legend.position = "none") +
                labs(title = element_blank(),
                     x = element_blank(),
                     y = element_blank())
            ggsave(paste0(outdir, "DESeq2_", b, "_vs_", a, ".png"), width = 3, height = 3, dpi = 600)

            # Significant DEGs

            resData_sig <- resData %>%
                filter(colour == "red")

            # Save to file

            fwrite(x = resData,
                   file = paste0(outdir,
                                 "allDEG_",
                                 b,
                                 "_vs_",
                                 a,
                                 ".tsv"),
                   sep = "\t")

            fwrite(x = resData_sig,
                   file = paste0(outdir,
                                 "sigDEG_",
                                 b,
                                 "_vs_",
                                 a,
                                 ".tsv"),
                   sep = "\t")
        }
    }
}
