#' Find significant crosstalk between pairwise clusters
#'
#' @description Find significant crosstalk between pairwise clusters with permutation test, which might take a long time.
#' @param clu_pairs A list generated from \code{\link{FindPairs}}
#' @param per_num Number of repeat times for permutation test. Default is 1000
#' @param pvalue Include the significantly cell-cell communications
#' with this cutoff of p value from permutation test. Default is 0.05
#' @return The same list containing the significant cell-cell communications stored
#' in \code{clu_clu_num} and in \code{clu_clu_score}. For p value, 0 means p < 1/per_num, namely p < 0.001 by default.
#' @import progress
#' @importFrom crayon red cyan green
#' @importFrom stats sd
#' @export CrossTalkSig
#' @examples
#' clu_pairs <- clu_pairs <- FindPairs(object = mouse_kidney_203_Seurat,species = 'Mouse')
#' clu_pairs <- CrossTalkSig(clu_pairs = clu_pairs, per_num = 100)

CrossTalkSig <- function(clu_pairs = NULL, per_num = 1000, pvalue = 0.05) {
    # check clu_pairs
    if (is.null(clu_pairs)) {
        stop("Please input the list from the function of FindPairs")
    }
    if (!is.data.frame(clu_pairs[["clu_crosstalk"]])) {
        stop("No available clu_crosstalk in clu_pairs")
    }
    if (!is.data.frame(clu_pairs[["clu_crosstalk_undir"]])) {
        stop("No available clu_crosstalk_undir in clu_pairs")
    }
    # check per_num
    if (per_num < 100) {
        stop("The minimum is 100!")
    }
    # para
    para <- clu_pairs[["para"]]
    species <- para[["species"]]
    use_LRdb <- para[["use_LRdb"]]
    cell_min_pct <- para[["cell_min_pct"]]
    p_value <- para[["p_value"]]
    # check LRdb
    if (!is.data.frame(use_LRdb)) {
        if (!use_LRdb %in% c("LRdb", "LRdb1", "LRdb2", "LRdb3", "LRdb4", "LRdb5", "LRdb6")) {
            stop("Please input the right LRdb! 'LRdb','LRdb1','LRdb2','LRdb3','LRdb4','LRdb5','LRdb6'")
        }
        if (use_LRdb == "LRdb") {
            LR_database <- LR_database[LR_database$species == species & LR_database$LRdb == "CellTalkDB", ]
        }
        if (use_LRdb %in% c("LRdb1", "LRdb2", "LRdb3")) {
            if (species != "Human") {
                stop("Please set species as 'Human'")
            }
            if (use_LRdb == "LRdb1") {
                LR_database <- LR_database[LR_database$LRdb == "CellTalkDB" & LR_database$species %in% c("Human", "Mouse2Human"),
                  ]
            }
            if (use_LRdb == "LRdb2") {
                LR_database <- LR_database[LR_database$LRdb %in% c("CellTalkDB", "SingleCellSignalR") & LR_database$species ==
                  "Human", ]
            }
            if (use_LRdb == "LRdb3") {
                LR_database <- LR_database[LR_database$LRdb %in% c("CellTalkDB", "SingleCellSignalR") & LR_database$species %in%
                  c("Human", "Mouse2Human"), ]
            }
        }
        if (use_LRdb %in% c("LRdb4", "LRdb5", "LRdb6")) {
            if (species != "Mouse") {
                stop("Please set species as 'Mouse'")
            }
            if (use_LRdb == "LRdb4") {
                LR_database <- LR_database[LR_database$LRdb == "CellTalkDB" & LR_database$species %in% c("Mouse", "Human2Mouse"),
                  ]
            }
            if (use_LRdb == "LRdb5") {
                LR_database <- LR_database[LR_database$LRdb %in% c("CellTalkDB", "SingleCellSignalR") & LR_database$species ==
                  "Human2Mouse", ]
            }
            if (use_LRdb == "LRdb6") {
                LR_database <- LR_database[LR_database$LRdb %in% c("CellTalkDB", "SingleCellSignalR") & LR_database$species %in%
                  c("Mouse", "Human2Mouse"), ]
            }
        }
        LR_database <- LR_database[, -c(4, 5)]
    }
    if (is.data.frame(use_LRdb)) {
        cat(red("Please ensure the custom LR database is a dataframe containing 2 columns: first is 'ligand', second is 'receptor'",
            "\n"))
        if (ncol(use_LRdb) != 2) {
            stop("Please ensure the custom LR database is a dataframe containing 2 columns: first is 'ligand', second is 'receptor'")
        }
        Sys.sleep(2)
        LR_database <- use_LRdb
        colnames(LR_database) <- c("ligand_gene_symbol", "receptor_gene_symbol")
        LR_database$lr_pair <- paste(LR_database$ligand_gene_symbol, LR_database$receptor_gene_symbol, sep = "_")
        LR_database <- LR_database[, c(3, 1, 2)]
        use_LRdb <- "custom LRdb"
    }
    Sys.sleep(2)
    cat("Finding significant cell-cell communications between pairwise clusters", "\n")
    Sys.sleep(2)
    # Performing permutation test by random cell labels
    cat(cyan(paste0("Performing permutation test (repeating times: ", per_num, ")"), "\n"))
    Sys.sleep(2)
    cat(red("Note: this process might take a long time"))
    pb <- progress_bar$new(format = "[:bar] Finished::percent Remaining::eta", total = per_num, clear = FALSE, width = 60,
        complete = "+", incomplete = "-")
    ndata <- clu_pairs[["ndata"]]
    clu_info <- clu_pairs[["clu_info"]]
    clu_crosstalk <- clu_pairs[["clu_crosstalk"]]
    clu_crosstalk_undir <- clu_pairs[["clu_crosstalk_undir"]]
    clu_clu_num_per <- as.data.frame(matrix(0, nrow = nrow(clu_crosstalk), ncol = per_num))
    clu_clu_score_per <- as.data.frame(matrix(0, nrow = nrow(clu_crosstalk), ncol = per_num))
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
        d1 <- data.frame(cluster1 = rep(clu_num[i], length(clu_num)), cluster2 = clu_num, stringsAsFactors = F)
        clu_pair <- rbind(clu_pair, d1)
    }
    # add a function to calculate gene percentage in scdata
    gene_pct <- function(scdata) {
        scdata1 <- scdata[scdata > 0]
        return(length(scdata1)/length(scdata))
    }
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
        for (i in 1:nrow(clu_pair)) {
            clu_pair1 <- clu_pair[i, ]
            ndata1 <- ndata[, clu_info_per[clu_info_per$cluster == clu_pair1$cluster1, ]$cell]
            ndata2 <- ndata[, clu_info_per[clu_info_per$cluster == clu_pair1$cluster2, ]$cell]
            res_pairs1 <- LR_database
            res_pairs1$ligand_clu <- clu_pair1$cluster1
            res_pairs1$receptor_clu <- clu_pair1$cluster2
            res_pairs1$ligand_pct <- 0
            res_pairs1$ligand_exp_avg <- 0
            res_pairs1$receptor_pct <- 0
            res_pairs1$receptor_exp_avg <- 0
            res_pairs1$lr_score <- 0
            # match ligand
            ligand_pct <- apply(ndata1, 1, gene_pct)
            ligand_mean <- apply(ndata1, 1, mean)
            ligand_info <- data.frame(gene = names(ligand_pct), pct = as.numeric(ligand_pct), exp_avg = as.numeric(ligand_mean),
                stringsAsFactors = F)
            ligand_info$z_score <- as.numeric(scale(ligand_info$exp_avg))
            ligand_info$p_value <- pnorm(q = ligand_info$exp_avg, mean = mean(ligand_info$exp_avg), sd = sd(ligand_info$exp_avg),
                lower.tail = F)
            ligand_info <- ligand_info[ligand_info$p_value < p_value, ]
            if (nrow(ligand_info) > 0) {
                colnames(ligand_info)[1] <- "ligand_gene_symbol"
                res_pairs1 <- base::merge(res_pairs1, ligand_info)
                res_pairs1$ligand_pct <- res_pairs1$pct
                res_pairs1$ligand_exp_avg <- res_pairs1$exp_avg
                res_pairs1 <- res_pairs1[, -c(11:14)]
                res_pairs1 <- res_pairs1[, c(2, 1, 3:10)]
            }
            # match receptor
            receptor_pct <- apply(ndata2, 1, gene_pct)
            receptor_mean <- apply(ndata2, 1, mean)
            receptor_info <- data.frame(gene = names(receptor_pct), pct = as.numeric(receptor_pct), exp_avg = as.numeric(receptor_mean),
                stringsAsFactors = F)
            receptor_info$z_score <- as.numeric(scale(receptor_info$exp_avg))
            receptor_info$p_value <- pnorm(q = receptor_info$exp_avg, mean = mean(receptor_info$exp_avg), sd = sd(receptor_info$exp_avg),
                lower.tail = F)
            receptor_info <- receptor_info[receptor_info$p_value < p_value, ]
            if (nrow(receptor_info) > 0) {
                colnames(receptor_info)[1] <- "receptor_gene_symbol"
                res_pairs1 <- base::merge(res_pairs1, receptor_info)
                res_pairs1$receptor_pct <- res_pairs1$pct
                res_pairs1$receptor_exp_avg <- res_pairs1$exp_avg
                res_pairs1 <- res_pairs1[, -c(11:14)]
                res_pairs1 <- res_pairs1[, c(2, 3, 1, 4:10)]
            }
            res_pairs1$lr_score <- res_pairs1$ligand_exp_avg * res_pairs1$receptor_exp_avg
            res_pairs1 <- res_pairs1[res_pairs1$ligand_pct > cell_min_pct & res_pairs1$receptor_pct > cell_min_pct, ]
            clu_clu_num_per[i, k] <- nrow(res_pairs1)
            clu_clu_score_per[i, k] <- sum(res_pairs1$lr_score)
        }
        pb$tick()
    }
    cat(green("***Done***", "\n"))
    Sys.sleep(2)
    # calculating p value
    clu_clu_num_per$real <- clu_crosstalk$LR_number
    clu_clu_score_per$real <- clu_crosstalk$LR_score
    # add a function to calculate p value
    clu_p_value <- function(scdata) {
        scdata_real <- scdata[length(scdata)]
        scdata <- scdata[-length(scdata)]
        scdata_num <- scdata[scdata < scdata_real]
        scdata_num <- length(scdata) - length(scdata_num)
        return(scdata_num/length(scdata))
    }
    clu_crosstalk$LR_number_pvalue <- apply(clu_clu_num_per, 1, clu_p_value)
    clu_crosstalk$LR_score_pvalue <- apply(clu_clu_score_per, 1, clu_p_value)
    clu_crosstalk$LR_number_type <- "Unsignificant"
    if (nrow(clu_crosstalk[clu_crosstalk$LR_number_pvalue < pvalue, ]) > 0) {
        clu_crosstalk[clu_crosstalk$LR_number_pvalue < pvalue, ]$LR_number_type <- "Significant"
    }
    clu_crosstalk$LR_score_type <- "Unsignificant"
    if (nrow(clu_crosstalk[clu_crosstalk$LR_score_pvalue < pvalue, ]) > 0) {
        clu_crosstalk[clu_crosstalk$LR_score_pvalue < pvalue, ]$LR_score_type <- "Significant"
    }
    clu_com <- clu_crosstalk_undir$clu_com
    clu_crosstalk_undir$LR_number_pvalue <- 1
    clu_crosstalk_undir$LR_score_pvalue <- 1
    for (i in 1:length(clu_com)) {
        clu_crosstalk1 <- which(clu_crosstalk$clu_com == clu_com[i])
        # LR_num_pvalue
        clu_clu_num_per1 <- clu_clu_num_per[clu_crosstalk1, ]
        clu_clu_num_per1 <- as.data.frame(matrix(colSums(clu_clu_num_per1), nrow = 1, byrow = T))
        clu_crosstalk_undir$LR_number_pvalue[i] <- apply(clu_clu_num_per1, 1, clu_p_value)
        # LR_score_pvalue
        clu_clu_score_per1 <- clu_clu_score_per[clu_crosstalk1, ]
        clu_clu_score_per1 <- as.data.frame(matrix(colSums(clu_clu_score_per1), nrow = 1, byrow = T))
        clu_crosstalk_undir$LR_score_pvalue[i] <- apply(clu_clu_score_per1, 1, clu_p_value)
    }
    clu_crosstalk_undir$LR_number_type <- "Unsignificant"
    clu_crosstalk_undir[clu_crosstalk_undir$LR_number_pvalue < pvalue, ]$LR_number_type <- "Significant"
    clu_crosstalk_undir$LR_score_type <- "Unsignificant"
    clu_crosstalk_undir[clu_crosstalk_undir$LR_score_pvalue < pvalue, ]$LR_score_type <- "Significant"
    clu_pairs[["clu_crosstalk"]] <- clu_crosstalk
    clu_pairs[["clu_crosstalk_undir"]] <- clu_crosstalk_undir
    return(clu_pairs)
}



