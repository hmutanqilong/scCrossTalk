#' Plot LR pairs between pairwise clusters
#'
#' @description Plot LR pairs between pairwise clusters with a heatmap of LR scores for each directed pairwise clusters
#' @param clu_pairs A list generated from \code{\link{FindPairs}}
#' @param LR_pairs Which LR pairs to show. Default is 'all'. Or custom LR pairs, e.g., App_Tspan12, App_Notch2
#' @param show_clusters Which clusters related LR pairs to show. Default is 'all'
#' @param show_sig Whether to show significant highly expressed LR pairs between pairwise clusters. Default is FALSE
#' To show significant LR pairs, please run \code{\link{PairsSig}}
#' @param color_low Color of the lowest score. Default is 'white'
#' @param color_high Color of the highest score. Default is 'red'
#' @param border_color Color of cell borders on heatmap, use NA if no border should be drawn.Default is 'grey60'
#' @param cluster_rows Boolean values determining if rows should be clustered. Default is TRUE
#' @param cluster_cols Boolean values determining if columns should be clustered. Default is TRUE
#' @param symbol Significant symbol to plot. Default is '*'
#' @param symbol_col Color of symbol. Default is 'black'
#' @param symbol_size Size of symbol. Default is 12
#' @return A heatmap of LR scores for each directed pairwise clusters
#' @examples
#' clu_pairs <- FindPairs(object = mouse_kidney_203_Seurat,species = 'Mouse')
#'
#' # Plot all LR pairs related with cluster 1
#' plotpairsHeat (clu_pairs = clu_pairs,
#'                LR_pairs = 'all',
#'                show_clusters = '1')
#'
#' # Plot App_Tspan12 pair related with all clusters
#' plotpairsHeat (clu_pairs = clu_pairs,
#'                LR_pairs = 'App_Tspan12',
#'                show_clusters = 'all')
#'
#' # Plot App_Tspan12 and App_Notch2 pairs related with cluster 1 and 2
#' plotpairsHeat (clu_pairs = clu_pairs,
#'                LR_pairs = c('App_Tspan12','App_Notch2'),
#'                show_clusters = c('1','2'))
#' @import pheatmap
#' @importFrom grDevices colorRampPalette
#' @export PlotPairsHeat

PlotPairsHeat <- function(clu_pairs = NULL, LR_pairs = "all", show_clusters = "all", show_sig = F, 
    color_low = "white", color_high = "red", border_color = "grey60", cluster_rows = T, cluster_cols = T, 
    symbol = "*", symbol_col = "black", symbol_size = 12) {
    # check
    if (is.null(clu_pairs)) {
        stop("Please input the list from the function of FindPairs")
    }
    res_pairs <- clu_pairs[["clu_pairs"]]
    if (nrow(res_pairs) == 0) {
        stop("There is no sigficant LR pairs in clu_pairs")
    }
    if (!"p_value" %in% colnames(res_pairs) & show_sig == T) {
        stop("To show significantly highly expressed LR pairs, please run PairsSig")
    }
    if (!is.character(LR_pairs)) {
        stop("Please input the right LR_pairs (character)")
    }
    if (!is.character(show_clusters)) {
        stop("Please input the right show_clusters (character)")
    }
    if (!is.logical(show_sig)) {
        stop("Please input the right show_sig (logical, TRUE or FALSE)")
    }
    if (!is.character(color_low)) {
        stop("Please input the right color_low (character)")
    }
    if (!is.character(color_high)) {
        stop("Please input the right color_high (character)")
    }
    if (!is.character(border_color)) {
        stop("Please input the right border_color (character)")
    }
    if (!is.logical(cluster_rows)) {
        stop("Please input the right cluster_rows (logical, TRUE or FALSE)")
    }
    if (!is.logical(cluster_cols)) {
        stop("Please input the right cluster_cols (logical, TRUE or FALSE)")
    }
    if (!is.character(symbol)) {
        stop("Please input the right symbol (character)")
    }
    if (!is.character(symbol_col)) {
        stop("Please input the right symbol_col (character)")
    }
    if (!is.numeric(symbol_size) | symbol_size < 0) {
        stop("Please input the right symbol_size (number, > 0)")
    }
    # check LR_pairs check show_clusters
    clu_num <- unique(c(res_pairs$ligand_clu, res_pairs$receptor_clu))
    if (!"all" %in% show_clusters) {
        if (!all(show_clusters %in% clu_num)) {
            show_clusters1 <- which(!show_clusters %in% clu_num)
            show_clusters <- show_clusters[show_clusters1]
            show_clusters1 <- show_clusters[1]
            if (length(show_clusters) > 1) {
                for (i in 2:length(show_clusters)) {
                  show_clusters1 <- paste(show_clusters1, show_clusters[i], sep = ",")
                }
            }
            stop(paste0(show_clusters1, ": not in clusters of clu_pairs list!"))
        }
        res_pairs <- res_pairs[res_pairs$ligand_clu %in% show_clusters | res_pairs$receptor_clu %in% 
            show_clusters, ]
        lr_pairs <- unique(res_pairs$lr_pair)
    }
    lr_pairs <- unique(res_pairs$lr_pair)
    if (!"all" %in% LR_pairs) {
        if (!all(LR_pairs %in% lr_pairs)) {
            LR_pairs1 <- which(!LR_pairs %in% lr_pairs)
            LR_pairs <- LR_pairs[LR_pairs1]
            LR_pairs1 <- LR_pairs[1]
            if (length(LR_pairs) > 1) {
                for (i in 2:length(LR_pairs)) {
                  LR_pairs1 <- paste(LR_pairs1, LR_pairs[i], sep = ",")
                }
            }
            stop(paste0(LR_pairs1, ": not in lr_pair of clu_pairs list!"))
        }
        lr_pairs <- lr_pairs[lr_pairs %in% LR_pairs]
        if (length(LR_pairs) == 1) {
            cluster_rows <- F
        }
    }
    clu_com <- unique(res_pairs$clu_com)
    res_pairs$clu_clu <- paste0(res_pairs$ligand_clu, "-", res_pairs$receptor_clu)
    clu_clu <- NULL
    for (i in 1:length(clu_com)) {
        clu_clu <- c(clu_clu, unique(res_pairs[res_pairs$clu_com == clu_com[i], ]$clu_clu))
    }
    clu_plot <- as.data.frame(matrix(0, nrow = length(lr_pairs), ncol = length(clu_clu)))
    colnames(clu_plot) <- clu_clu
    rownames(clu_plot) <- lr_pairs
    clu_plot_sig <- as.data.frame(matrix("", nrow = length(lr_pairs), ncol = length(clu_clu)), 
        stringsAsFactors = F)
    colnames(clu_plot_sig) <- clu_clu
    rownames(clu_plot_sig) <- lr_pairs
    for (i in 1:ncol(clu_plot)) {
        res_pairs1 <- res_pairs[res_pairs$clu_clu == clu_clu[i], ]
        res_pairs1 <- res_pairs1[res_pairs1$lr_pair %in% lr_pairs, ]
        if (nrow(res_pairs1) > 0) {
            for (j in 1:nrow(res_pairs1)) {
                clu_plot[res_pairs1$lr_pair[j], i] <- res_pairs1$lr_score[j]
                if (res_pairs1$type[j] == "Significant") {
                  clu_plot_sig[res_pairs1$lr_pair[j], i] <- symbol
                }
            }
        }
    }
    lr_pairs <- rownames(clu_plot)
    for (i in 1:length(lr_pairs)) {
        lr_pairs1 <- res_pairs[res_pairs$lr_pair == lr_pairs[i], ]
        lr_pairs1 <- lr_pairs1[1, c("ligand_gene_symbol", "receptor_gene_symbol")]
        lr_pairs1 <- paste0(lr_pairs1$ligand_gene_symbol, "-", lr_pairs1$receptor_gene_symbol)
        lr_pairs1 <- unique(lr_pairs1)
        lr_pairs[i] <- lr_pairs1
    }
    rownames(clu_plot) <- lr_pairs
    # plot
    if (show_sig == T) {
        pheatmap(mat = clu_plot, color = colorRampPalette(c(color_low, color_high))(50), border_color = border_color, 
            cluster_rows = cluster_rows, cluster_cols = cluster_cols, display_numbers = as.matrix(clu_plot_sig), 
            number_color = symbol_col, fontsize_number = symbol_size)
    }
    if (show_sig == F) {
        pheatmap(mat = clu_plot, color = colorRampPalette(c(color_low, color_high))(50), border_color = border_color, 
            cluster_rows = cluster_rows, cluster_cols = cluster_cols)
    }
}
