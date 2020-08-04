#' Plot crosstalk between pairwise clusters
#'
#' @description Plot crosstalk between pairwise clusters by sum of LR pairs number and scores for pairwise clusters with heatmap
#' @param clu_pairs A list generated from \code{\link{FindPairs}}
#' @param show_type Which to show, sum of LR pairs number or LR score
#' Input 'number' to show sum of LR pairs number; Input 'score' to show sum of LR pairs score
#' @param if_directed Considering the direction of LR pairs between two clusters. Default is TRUE
#' @param show_sig Whether to show significant crosstalk between pairwise clusters. Default is FALSE
#' To show significant crosstalk, please run \code{\link{CrossTalkSig}}
#' @param color_low Color of the lowest score. Default is 'white'
#' @param color_high Color of the highest score. Default is 'red'
#' @param border_color Color of cell borders on heatmap, use NA if no border should be drawn.Default is 'grey60'
#' @param cluster_rows Boolean values determining if rows should be clustered. Default is TRUE
#' @param cluster_cols Boolean values determining if columns should be clustered. Default is TRUE
#' @param symbol Significant symbol to plot. Default is '*'
#' @param symbol_col Color of the significant symbol. Default is 'black'
#' @param symbol_size Size of the significant symbol. Default is 12
#' @return A corrplot containing sum of LR pairs number or scores. When
#' if_directed is TRUE, rows mean the source(ligand),columns means the target(receptor)
#' @examples # Directed plot
#' plotpairsCor(clu_pairs = FindPairs(object = mouse_kidney_203_Seurat,species = 'Mouse'),
#'              if_directed = T)
#'
#' @examples # Undirected plot
#' plotpairsCor(clu_pairs = FindPairs(object = mouse_kidney_203_Seurat,species = 'Mouse'),
#'              if_directed = F)
#' @import pheatmap
#' @export PlotCrossTalkHeat

PlotCrossTalkHeat <- function(clu_pairs = NULL, show_type = NULL, if_directed = T, show_sig = F, color_low = "white", 
    color_high = "red", border_color = "grey60", cluster_rows = T, cluster_cols = T, symbol = "*", symbol_col = "black", 
    symbol_size = 12) {
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
    clu_crosstalk <- clu_pairs[["clu_crosstalk"]]
    clu_crosstalk_undir <- clu_pairs[["clu_crosstalk_undir"]]
    if (!"LR_number_pvalue" %in% colnames(clu_crosstalk) & show_sig == T) {
        stop("To show significant crosstalk, please run CrossTalkSig")
    }
    if (show_type != "number" & show_type != "score") {
        stop("Please input the right show_type!")
    }
    # add a function
    clu_sort <- function(clu_data) {
        clu_data <- clu_data[order(clu_data)]
        clu_data <- paste(clu_data[1], clu_data[2], sep = "_")
        return(clu_data)
    }
    # get clusters
    clu_info <- clu_pairs[["clu_info"]]
    clu_num <- unique(clu_info[, 2])
    clu_num1 <- clu_num[nchar(clu_num) == 1]
    clu_num2 <- clu_num[nchar(clu_num) == 2]
    clu_num1 <- clu_num1[order(clu_num1)]
    clu_num2 <- clu_num2[order(clu_num2)]
    clu_num <- c(clu_num1, clu_num2)
    # calculating sum of LR pairs and scores
    if (show_type == "number") {
        main_title <- "Sum of LR pairs number"
        cc_pairs <- as.data.frame(matrix(0, nrow = length(clu_num), ncol = length(clu_num)))
        colnames(cc_pairs) <- clu_num
        rownames(cc_pairs) <- clu_num
        cc_pairs_sig <- cc_pairs
        # row:source(ligand) col:target(receptor)
        for (i in 1:nrow(cc_pairs)) {
            for (j in 1:ncol(cc_pairs)) {
                clu_clu_name <- paste(rownames(cc_pairs)[i], colnames(cc_pairs)[j], sep = "_")
                if (if_directed == F) {
                  if (clu_clu_name %in% clu_crosstalk_undir$clu1_clu2) {
                    clu_crosstalk_undir1 <- clu_crosstalk_undir[clu_crosstalk_undir$clu1_clu2 == clu_clu_name, ]
                    cc_pairs[i, j] <- clu_crosstalk_undir1$LR_number
                    cc_pairs[j, i] <- clu_crosstalk_undir1$LR_number
                    if (show_sig == T) {
                      cc_pairs_sig[i, j] <- ""
                      cc_pairs_sig[j, i] <- ""
                      if (clu_crosstalk_undir1$LR_number_type == "Significant") {
                        cc_pairs_sig[i, j] <- symbol
                        cc_pairs_sig[j, i] <- symbol
                      }
                    }
                  }
                }
                if (if_directed == T) {
                  if (clu_clu_name %in% clu_crosstalk$clu1_clu2) {
                    clu_crosstalk1 <- clu_crosstalk[clu_crosstalk$clu1_clu2 == clu_clu_name, ]
                    cc_pairs[i, j] <- clu_crosstalk1$LR_number
                    if (show_sig == T) {
                      cc_pairs_sig[i, j] <- ""
                      if (clu_crosstalk1$LR_number_type == "Significant") {
                        cc_pairs_sig[i, j] <- symbol
                      }
                    }
                  }
                }
            }
        }
    }
    if (show_type == "score") {
        main_title <- "Sum of LR pairs score"
        cc_pairs <- as.data.frame(matrix(0, nrow = length(clu_num), ncol = length(clu_num)))
        colnames(cc_pairs) <- clu_num
        rownames(cc_pairs) <- clu_num
        cc_pairs_sig <- cc_pairs
        # row:source(ligand) col:target(receptor)
        for (i in 1:nrow(cc_pairs)) {
            for (j in 1:ncol(cc_pairs)) {
                clu_clu_name <- paste(rownames(cc_pairs)[i], colnames(cc_pairs)[j], sep = "_")
                if (if_directed == F) {
                  if (clu_clu_name %in% clu_crosstalk_undir$clu1_clu2) {
                    clu_crosstalk_undir1 <- clu_crosstalk_undir[clu_crosstalk_undir$clu1_clu2 == clu_clu_name, ]
                    cc_pairs[i, j] <- clu_crosstalk_undir1$LR_score
                    cc_pairs[j, i] <- clu_crosstalk_undir1$LR_score
                    if (show_sig == T) {
                      cc_pairs_sig[i, j] <- ""
                      cc_pairs_sig[j, i] <- ""
                      if (clu_crosstalk_undir1$LR_score_type == "Significant") {
                        cc_pairs_sig[i, j] <- symbol
                        cc_pairs_sig[j, i] <- symbol
                      }
                    }
                  }
                }
                if (if_directed == T) {
                  if (clu_clu_name %in% clu_crosstalk$clu1_clu2) {
                    clu_crosstalk1 <- clu_crosstalk[clu_crosstalk$clu1_clu2 == clu_clu_name, ]
                    cc_pairs[i, j] <- clu_crosstalk1$LR_score
                    if (show_sig == T) {
                      cc_pairs_sig[i, j] <- ""
                      if (clu_crosstalk1$LR_score_type == "Significant") {
                        cc_pairs_sig[i, j] <- symbol
                      }
                    }
                  }
                }
            }
        }
    }
    if (show_sig == T) {
        pheatmap(mat = cc_pairs, color = colorRampPalette(c(color_low, color_high))(100), border_color = border_color, 
            cluster_rows = cluster_rows, cluster_cols = cluster_cols, main = main_title, fontsize = 12, angle_col = 0, 
            display_numbers = as.matrix(cc_pairs_sig), number_color = symbol_col, fontsize_number = symbol_size)
    }
    if (show_sig == F) {
        pheatmap(mat = cc_pairs, color = colorRampPalette(c(color_low, color_high))(100), border_color = border_color, 
            cluster_rows = cluster_rows, cluster_cols = cluster_cols, main = main_title, fontsize = 12, angle_col = 0)
    }
}
