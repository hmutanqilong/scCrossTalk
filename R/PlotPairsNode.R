#' Plot LR pairs between pairwise clusters
#'
#' @description Plot LR pairs between pairwise clusters with a node plot of LR pairs between pairwise clusters
#' @param clu_pairs A list generated from \code{\link{FindPairs}}
#' @param show_sig Whether to show significant highly expressed LR pairs between pairwise clusters. Default is FALSE
#' To show significant LR pairs, please run \code{\link{PairsSig}}
#' @param edge_width Width of edge. Default is 1
#' @param edge_alpha Transparency of edge. Default is 0.5
#' @param node_size_min Minimum size of node. Default is 1
#' @param node_size_max Maximum size of node. Default is 10
#' @param text_size Size of text labels. Default is 3
#' @examples
#' clu_pairs <- FindPairs(object = mouse_kidney_203_Seurat,species = 'Mouse')
#'
#' PlotPairsNode(clu_pairs = clu_pairs, show_sig = F)
#' @return A node plot of LR pairs between pairwise clusters. The number behind gene is the cluster,
#' L, means ligand;R, means receptor. Size of gene node is the mean expression of this gene in this cluster
#' @import ggraph
#' @import ggplot2
#' @importFrom scales hue_pal
#' @importFrom igraph graph_from_data_frame
#' @export PlotPairsNode

PlotPairsNode <- function(clu_pairs = NULL, show_sig = F, edge_width = 1, edge_alpha = 0.5, node_size_min = 1, node_size_max = 10, 
    text_size = 3) {
    # check
    if (is.null(clu_pairs)) {
        stop("Please input the list from the function of FindPairs")
    }
    if (!is.list(clu_pairs)) {
        stop("Please input the right clu_pairs (generated from FindPairs)")
    }
    res_pairs <- clu_pairs[["clu_pairs"]]
    if (nrow(res_pairs) == 0) {
        stop("There is no sigficantly highly expressed LR in clu_pairs")
    }
    if (!is.logical(show_sig)) {
        stop("Please input the right show_sig (logical, TRUE or FALSE)")
    }
    if (!"p_value" %in% colnames(res_pairs) & show_sig == T) {
        stop("To show significantly highly expressed LR pairs, please run PairsSig")
    }
    if (edge_width < 0 | !is.numeric(edge_width)) {
        stop("Please input the right edge_width (number, > 0)")
    }
    if (edge_alpha < 0 | edge_alpha > 1 | !is.numeric(edge_alpha)) {
        stop("Please input the right edge_alpha (number, 0 < edge_alpha < 1)")
    }
    if (node_size_min < 0 | !is.numeric(node_size_min)) {
        stop("Please input the right node_size_min (number, > 0)")
    }
    if (node_size_max < 0 | !is.numeric(node_size_max)) {
        stop("Please input the right node_size_max (number, > 0)")
    }
    if (node_size_min >= node_size_max) {
        stop("Please input the right node_size_min and node_size_max (node_size_min < node_size_max)")
    }
    if (text_size < 0 | !is.numeric(text_size)) {
        stop("Please input the right text_size (number, > 0)")
    }
    res_pairs$ligand_gene_symbol <- paste0(res_pairs$ligand_gene_symbol, "(", res_pairs$ligand_clu, "L)")
    res_pairs$receptor_gene_symbol <- paste0(res_pairs$receptor_gene_symbol, "(", res_pairs$receptor_clu, "R)")
    pairs_sig <- 1:nrow(res_pairs)
    if (show_sig == T) {
        pairs_sig <- which(res_pairs$type == "Significant")
    }
    plot_res_con <- res_pairs[pairs_sig, c("ligand_gene_symbol", "receptor_gene_symbol", "ligand_clu")]
    colnames(plot_res_con) <- c("from", "to", "cluster")
    plot_res_con$value <- 1
    # raw clu_num
    clu_info <- clu_pairs[["clu_info"]]
    clu_num <- unique(clu_info[, 2])
    clu_num1 <- clu_num[nchar(clu_num) == 1]
    clu_num2 <- clu_num[nchar(clu_num) == 2]
    clu_num1 <- clu_num1[order(clu_num1)]
    clu_num2 <- clu_num2[order(clu_num2)]
    clu_num_raw <- c(clu_num1, clu_num2)
    plot_res <- data.frame(name = "0", size = 0, cluster = "0", stringsAsFactors = F)
    for (i in 1:length(clu_num_raw)) {
        clu_num_raw1 <- clu_num_raw[i]
        if (clu_num_raw1 %in% res_pairs$ligand_clu) {
            res_pairs1 <- res_pairs[res_pairs$ligand_clu == clu_num_raw1, ]
            res_pairs1 <- res_pairs1[, c("ligand_gene_symbol", "ligand_exp_avg", "ligand_clu")]
            res_pairs1 <- unique(res_pairs1)
            colnames(res_pairs1) <- c("name", "size", "cluster")
            res_pairs1 <- res_pairs1[order(-res_pairs1$size), ]
            plot_res <- rbind(plot_res, res_pairs1)
        }
        if (clu_num_raw1 %in% res_pairs$receptor_clu) {
            res_pairs1 <- res_pairs[res_pairs$receptor_clu == clu_num_raw1, ]
            res_pairs1 <- res_pairs1[, c("receptor_gene_symbol", "receptor_exp_avg", "receptor_clu")]
            res_pairs1 <- unique(res_pairs1)
            colnames(res_pairs1) <- c("name", "size", "cluster")
            res_pairs1 <- res_pairs1[order(-res_pairs1$size), ]
            plot_res <- rbind(plot_res, res_pairs1)
        }
    }
    plot_res <- plot_res[-1, ]
    # factor cluster
    clu_num <- unique(plot_res$cluster)
    clu_num <- clu_num[clu_num %in% clu_num_raw]
    clu_num1 <- clu_num[nchar(clu_num) == 1]
    clu_num2 <- clu_num[nchar(clu_num) == 2]
    clu_num1 <- clu_num1[order(clu_num1)]
    clu_num2 <- clu_num2[order(clu_num2)]
    clu_num <- c(clu_num1, clu_num2)
    plot_res$cluster <- factor(plot_res$cluster, levels = clu_num)
    plot_res$id <- 1:nrow(plot_res)
    # angle
    angle <- 360 * (plot_res$id - 0.5)/nrow(plot_res)
    plot_res$angle <- ifelse(angle > 180, 90 - angle + 180, 90 - angle)
    plot_res$hjust <- ifelse(angle > 180, 1, 0)
    # plot_clu_col
    clu_col <- hue_pal()(length(clu_num_raw))
    plot_clu_col <- NULL
    for (i in 1:length(clu_num)) {
        plot_clu_col1 <- which(clu_num_raw == clu_num[i])
        plot_clu_col1 <- clu_col[plot_clu_col1]
        plot_clu_col <- c(plot_clu_col, plot_clu_col1)
    }
    mygraph <- graph_from_data_frame(plot_res_con, vertices = plot_res, directed = FALSE)
    ggraph(mygraph, layout = "linear", circular = TRUE) + geom_edge_arc(aes(edge_colour = cluster), edge_alpha = edge_alpha, 
        edge_width = edge_width) + geom_node_point(aes(size = size, fill = cluster), shape = 21, color = "black", alpha = 0.9) + 
        scale_size_continuous(range = c(node_size_min, node_size_max)) + scale_fill_manual(values = plot_clu_col) + geom_node_text(aes(x = x * 
        1.06, y = y * 1.06, label = name, angle = angle, hjust = hjust, color = cluster), size = text_size) + scale_color_manual(values = plot_clu_col) + 
        scale_edge_color_manual(values = plot_clu_col) + expand_limits(x = c(-1.6, 1.6), y = c(-1.6, 1.6)) + coord_fixed() + 
        theme_minimal() + theme(legend.position = "none", panel.grid = element_blank(), axis.line = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), plot.margin = unit(c(0, 
            0, 0, 0), "null"), panel.spacing = unit(c(0, 0, 0, 0), "null"))
}
