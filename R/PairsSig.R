#' Find significantly highly expressed LR pairs
#'
#' @description Find significantly highly expressed LR pairs for pairwise clusters
#' @param clu_pairs A list generated from \code{\link{FindPairs}}
#' @param per_num Number of repeat times for permutation test. Default is 1000
#' @param pvalue Include the significantly highly expressed LR pairs
#' with this cutoff of p value from permutation test. Default is \code{0.05}.
#' @return The same list containing the significant LR pairs stored
#' in \code{clu_pairs}. For p value, 0 means p < 1/per_num, namely p < 0.001 by default.
#' @import progress
#' @importFrom crayon red cyan green
#' @export PairsSig
#' @examples clu_pairs <- PairsSig(clu_pairs = FindPairs(object = mouse_kidney_203_Seurat,species = 'Mouse'))

PairsSig <- function(clu_pairs = NULL, per_num = 1000, pvalue = 0.05) {
    # check clu_pairs
    if (is.null(clu_pairs)) {
        stop("Please input the list from the function of FindPairs")
    }
    res_pair <- clu_pairs[["clu_pairs"]]
    if (nrow(res_pair) == 0) {
        stop("No available LR pairs in clu_pairs")
    }
    # check per_num
    if (per_num < 100) {
        stop("The minimum is 100!")
    }
    Sys.sleep(2)
    cat("Finding significant LR pairs between pairwise clusters", "\n")
    Sys.sleep(2)
    # Performing permutation test by random cell labels
    cat(cyan(paste0("Performing permutation test (repeating times: ", per_num, ")"), "\n"))
    Sys.sleep(2)
    cat(red("Note: this process might take a long time"))
    Sys.sleep(2)
    pb <- progress_bar$new(format = "[:bar] Finished::percent Remaining::eta", total = per_num, 
        clear = FALSE, width = 60, complete = "+", incomplete = "-")
    ndata <- clu_pairs[["ndata"]]
    clu_info <- clu_pairs[["clu_info"]]
    res_pairs <- clu_pairs[["clu_pairs"]]
    res_pairs_per <- as.data.frame(matrix(0, nrow = nrow(res_pairs), ncol = per_num))
    # clu_num
    clu_num <- unique(clu_info[, 2])
    clu_num1 <- clu_num[nchar(clu_num) == 1]
    clu_num2 <- clu_num[nchar(clu_num) == 2]
    clu_num1 <- clu_num1[order(clu_num1)]
    clu_num2 <- clu_num2[order(clu_num2)]
    clu_num <- c(clu_num1, clu_num2)
    # clu_pair
    clu_pair <- NULL
    for (i in 1:length(clu_num)) {
        d1 <- data.frame(cluster1 = rep(clu_num[i], length(clu_num)), cluster2 = clu_num, 
            stringsAsFactors = F)
        clu_pair <- rbind(clu_pair, d1)
    }
    clu_pair$clu1_clu2 <- paste(clu_pair$cluster1, clu_pair$cluster2, sep = "_")
    # add a function to sort pairwise clusters
    clu_sort <- function(clu_data) {
        clu_data <- clu_data[order(clu_data)]
        clu_data <- paste(clu_data[1], clu_data[2], sep = "_")
        return(clu_data)
    }
    for (k in 1:per_num) {
        set.seed(k)
        clu_info_per <- clu_info
        cell_id <- sample(c(1:ncol(ndata)), ncol(ndata))
        cell_id <- clu_info$cell[cell_id]
        clu_info_per$cell <- cell_id
        score_per <- NULL
        for (i in 1:nrow(clu_pair)) {
            clu_pair1 <- clu_pair[i, ]
            if (clu_pair1$clu1_clu2 %in% res_pairs$clu1_clu2) {
                res_pairs1 <- res_pairs[res_pairs$clu1_clu2 == clu_pair1$clu1_clu2, ]
                ndata1 <- ndata[, clu_info_per[clu_info_per$cluster == clu_pair1$cluster1, 
                  ]$cell]
                ndata2 <- ndata[, clu_info_per[clu_info_per$cluster == clu_pair1$cluster2, 
                  ]$cell]
                ndata1 <- ndata1[rownames(ndata1) %in% res_pairs1$ligand_gene_symbol, ]
                ndata2 <- ndata2[rownames(ndata2) %in% res_pairs1$receptor_gene_symbol, ]
                ndata1 <- as.data.frame(t(as.matrix(ndata1)))
                ndata2 <- as.data.frame(t(as.matrix(ndata2)))
                ndata1 <- ndata1[, res_pairs1$ligand_gene_symbol]
                ndata2 <- ndata2[, res_pairs1$receptor_gene_symbol]
                score_per1 <- as.numeric(apply(ndata1, 2, mean)) * as.numeric(apply(ndata2, 
                  2, mean))
                score_per <- c(score_per, score_per1)
            }
        }
        res_pairs_per[, k] <- score_per
        pb$tick()
    }
    cat(green("***Done***", "\n"))
    Sys.sleep(2)
    # calculating p value
    res_pairs_per$real <- res_pairs$lr_score
    # add a function to calculate p value
    clu_p_value <- function(scdata) {
        scdata_real <- scdata[length(scdata)]
        scdata <- scdata[-length(scdata)]
        scdata_num <- scdata[scdata < scdata_real]
        scdata_num <- length(scdata) - length(scdata_num)
        return(scdata_num/length(scdata))
    }
    res_pairs$p_value <- apply(res_pairs_per, 1, clu_p_value)
    res_pairs$type <- "Unsignificant"
    res_pairs[res_pairs$p_value < pvalue, ]$type <- "Significant"
    clu_pairs[["clu_pairs"]] <- res_pairs
    return(clu_pairs)
}
