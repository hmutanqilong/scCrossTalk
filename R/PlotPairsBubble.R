#' Plot LR pairs between pairwise clusters
#'
#' @description Plot LR pairs between pairwise clusters with a bubble plot of LR pairs for pairwise clusters
#' @param clu_pairs A list generated from \code{\link{FindPairs}}
#' @param LR_pairs Which LR pairs to show. Default is 'all'. Or custom LR pairs, e.g., App_Tspan12, App_Notch2
#' @param show_clusters Which clusters related LR pairs to show. Default is 'all'
#' @param if_directed Considering the direction of LR pairs between two clusters. Default is TRUE
#' @param show_sig Whether to show significant highly expressed LR pairs between pairwise clusters. Default is FALSE
#' To show significant LR pairs, please run \code{\link{PairsSig}}
#' @param bubble_col Color of bubble. Default is 'black'
#' @param bubble_alpha Transparency of bubble. Default is 0.6
#' @param bubble_max_size Maximum size of bubble. Default is 20
#' @param show_text_cutoff Cutoff to show text label (LR score). Default is 0, namely to show all texts
#' @return A bubble plot of LR pairs for pairwise clusters
#' @examples
#' clu_pairs <- FindPairs(object = mouse_kidney_203_Seurat,species = 'Mouse')
#'
#' # Plot all LR pairs related with cluster 1
#' PlotPairsBubble (clu_pairs = clu_pairs,
#'                  LR_pairs = 'all',
#'                  show_clusters = '1')
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @export PlotPairsBubble

PlotPairsBubble <- function(clu_pairs = NULL, LR_pairs = "all", show_clusters = "all", if_directed = T, show_sig = F, 
    bubble_col = "black", bubble_alpha = 0.6, bubble_max_size = 20, show_text_cutoff = 1) {
    # check
    if (is.null(clu_pairs)) {
        stop("Please input the list from the function of FindPairs")
    }
    res_pairs <- clu_pairs[["clu_pairs"]]
    if (nrow(res_pairs) == 0) {
        stop("No available LR pairs in clu_pairs")
    }
    if (!"p_value" %in% colnames(res_pairs) & show_sig == T) {
        stop("To show significantly highly expressed LR pairs, please run PairsSig")
    }
    # check show_clusters
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
        res_pairs <- res_pairs[res_pairs$ligand_clu %in% show_clusters | res_pairs$receptor_clu %in% show_clusters, ]
        lr_pairs <- unique(res_pairs$lr_pair)
    }
    # check LR_pairs
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
    res_pairs$LR_pairs <- paste0(res_pairs$ligand_gene_symbol, "-", res_pairs$receptor_gene_symbol)
    if (show_sig == T) {
        res_pairs <- res_pairs[res_pairs$type == "Significant", ]
        if (nrow(res_pairs) == 0) {
            stop("No available significant LR pairs in clu_pairs")
        }
    }
    if (if_directed == F) {
        clu_sort <- function(clu_data) {
            clu_data <- clu_data[order(clu_data)]
            clu_data <- paste(clu_data[1], clu_data[2], sep = "-")
            return(clu_data)
        }
        res_pairs$clu_clu <- apply(res_pairs[, c("ligand_clu", "receptor_clu")], 1, clu_sort)
    }
    text_label <- ifelse(res_pairs$lr_score > show_text_cutoff, res_pairs$lr_pair, "")
    ggplot(data = res_pairs, aes(x = ligand_exp_avg, y = receptor_exp_avg)) + geom_point(aes(size = lr_score, fill = clu_clu), 
        shape = 21, colour = bubble_col, alpha = bubble_alpha) + geom_text_repel(label = text_label) + scale_size_area(max_size = bubble_max_size) + 
        guides(size = guide_legend((title = "LR score")), fill = guide_legend((title = "Crosstalk"))) + theme(legend.text = element_text(size = 10, 
        color = "black"), axis.title = element_text(size = 10, color = "black"), axis.text = element_text(size = 10, 
        color = "black"), legend.position = "right") + xlab("logCount (ligand average expression)") + ylab("logCount (receptor average expression)") + 
        theme_bw() + theme(panel.grid = element_blank())
}
