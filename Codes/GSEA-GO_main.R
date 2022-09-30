##################################################################################### INPUTS ####
#
# [1] main/Results/DESeq2/
#     DESeq2 results.
#
# [-] main/Data/GSEA/BioMartforGSC.tsv
#     KEGG dictionary for building  the gene set collections.  Optional file, will be dynamically
#     created if not supplied.
#
#################################################################################### OUTPUTS ####
#
# [1] main/Results/GSEA/
#     GSEA results.
#
############################################################################# INITIALISATION ####

require(biomaRt)
require(dplyr)
require(data.table)
require(piano)
require(foreach)
require(snowfall)

phenotypes <- c("Alzheimer's disease", "Control", "MCI", "Other dementia")
targetprefixes <- c("main/Results/DESeq2/allDEG_", "main/Results/GSEA-GO/GSEA_", "main/Results/GSEA-GO/GSEA_")
targetsuffixes <- c(".tsv", ".tsv", ".png")

i_o <-
    foreach (t = 1:length(targetprefixes), .final = function(t) setNames(t, c("in_DEG", "out_GSE_tsv", "out_GSE_png"))) %do% {
        foreach (a = phenotypes) %do% {
            foreach (b = phenotypes[which(phenotypes != a)]) %do% {
                paste0(targetprefixes[[t]], a, "_vs_", b, targetsuffixes[[t]])
            }
        } %>% unlist(F)
    }

infile_GSEASets <- "main/Data/GSEA/BioMartforGSC_GO.tsv"


ntop <- 0.05        # will take the top ntop*100 % of DEGs to determine GSEA


#### Functions #############################################################################################

### Get biomaRt results
fetchBM <- function(x,                                                                                      # Character or vector to submit to biomaRt
                    file = NULL                                                                             # File to check for reading or writing
) {
    # If file is specified, see if biomaRt results already exist, and if so, load them
    if (!is.null(file)) {                                                                                   # Check if file is specified
        try(
            suppressWarnings(res <- read.delim(file = file,                                                 # Try to load the file
                                               colClasses = "character")),
            silent = TRUE
        )
    }
    # If biomaRt results don't exist locally, fetch them from biomaRt and save them to file if specified
    while (!exists("res")) {                                                                                # Check if biomaRt results were loaded successfully
        try(                                                                                                # If no results were loaded successfuly, keep querying biomaRt until successful
            res <- getBM(attributes = c("ensembl_gene_id",
                                        "name_1006"),
                         filters = "ensembl_gene_id",
                         values = x,
                         mart = useMart("ensembl",
                                        dataset = "hsapiens_gene_ensembl")))
        if (exists("res") & (!is.null(file))) {                                                             # If results were obtained from bioMart (not locally), save to file
            write.table(x = res,
                        file = file,
                        sep = "\t")
        }
    }
    res
}

#### Main ##################################################################################################

doGSEA <- function(infile, outHeatmap, outTable) {

    data <- fread(infile)

    # Collect directions

    data <- data %>%
        filter(!is.na(pval)) %>%                                                                                # Remove na in pval
        mutate(dir = sign(l2fc),
               absstat = abs(stat)) %>%
        top_n(round(nrow(data) * ntop, 0), absstat)

    #dir <- data$dir
    #names(dir) <- data$gene

    # Collect deseq stats

    stat <- data$stat
    names(stat) <- data$gene

    # Query biomaRt to make gene set collections

    sets <- fetchBM(x = data$gene %>% unique(),
                    file = infile_GSEASets) %>%
        filter(name_1006 != "") %>%
        loadGSC()

    # Run GSEA

    gsaRes <- runGSA(geneLevelStats = stat,
                     directions = dir,
                     gsc = sets,
                     geneSetStat = "gsea",
                     ncpus = 5)

    # Save tables
    fwrite(x = GSAsummaryTable(gsaRes),
           file = outTable,
           sep = "\t")


}

for (i in 1:length(i_o$in_DEG)) {
    doGSEA(i_o$in_DEG[[i]], i_o$out_GSE_png[[i]], i_o$out_GSE_tsv[[i]])
}
