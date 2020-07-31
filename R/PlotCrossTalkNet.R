#' Plot crosstalk between pairwise clusters
#'
#' @description Plot crosstalk between pairwise clusters by sum of LR pairs number or score between pairwise clusters with net plot
#' @param clu_pairs A list generated from \code{\link{FindPairs}}
#' @param show_type Which to show, sum of LR pairs number or LR score.
#' Input 'number' to show sum of LR pairs number; Input 'score' to show sum of LR pairs score
#' @param layout Layout of net plot, e.g., 'fr','mds',''randomly','dh',
#' 'gem','graphopt','grid','sphere','kk','lgl'. Default is 'nicely'
#' @param show_sig Whether to show significant crosstalk between pairwise clusters. Default is FALSE
#' To show significant crosstalk, please run \code{\link{CrossTalkSig}}
#' @param edge_col Color of edge. Default is 'black'
#' @param edge_alpha Transparency of edge. Default is 0.1
#' @param node_size_min Minimum size of node. Default is 5
#' @param node_size_max Maximum size of node. Default is 10
#' @param text_size Size of text labels. Default is 3
#' @param text_col Color of text labels Default is 'black'
#' @examples
#' clu_pairs <- FindPairs(object = mouse_kidney_203_Seurat, species = 'Mouse')
#'
#' PlotCrossTalkNet(clu_pairs = clu_pairs, show_type = 'number')
#'
#' PlotCrossTalkNet(clu_pairs = clu_pairs, show_type = 'score', show_sig = T)
#' @return Sankey plot of LR pairs between pairwise clusters
#' @import ggalluvial ggplot2
#' @importFrom scales hue_pal
#' @export PlotCrossTalkNet

PlotCrossTalkNet <- function(clu_pairs = NULL, show_type = NULL, layout = "nicely", show_sig = F, 
    edge_col = "black", edge_alpha = 0.1, node_size_min = 5, node_size_max = 10, text_size = 3, 
    text_col = "black") {
    # check clu_pairs
    if (is.null(clu_pairs) | !is.list(clu_pairs)) {
        stop("Please input the list from the function of FindPairs")
    }
    if (!is.data.frame(clu_pairs[["clu_crosstalk"]])) {
        stop("No available clu_crosstalk in clu_pairs")
    }
    clu_crosstalk_undir <- clu_pairs[["clu_crosstalk_undir"]]
    if (!"LR_number_pvalue" %in% colnames(clu_crosstalk_undir) & show_sig == T) {
        stop("To show significant crosstalk, please run CrossTalkSig")
    }
    if (is.null(show_type) | !is.character(show_type)) {
        stop("Please input the right show_type (character)")
    }
    if (show_type != "number" & show_type != "score") {
        stop("Please input the right show_type!")
    }
    if (!is.character(layout)) {
        stop("Please input the right layout (character)")
    }
    if (!is.logical(show_sig)) {
        stop("Please input the right show_sig (logical, TRUE or FALSE)")
    }
    if (!is.character(edge_col)) {
        stop("Please input the right edge_col (character)")
    }
    if (!is.numeric(edge_alpha) | edge_alpha < 0) {
        stop("Please input the right edge_alpha (number, > 0)")
    }
    if (!is.numeric(node_size_min) | node_size_min < 0) {
        stop("Please input the right node_size_min (number, > 0)")
    }
    if (!is.numeric(node_size_max) | node_size_max < 0) {
        stop("Please input the right node_size_max (number, > 0)")
    }
    if (node_size_min >= node_size_max) {
        stop("Please input the right node_size_min and node_size_max (node_size_min < node_size_max)")
    }
    if (!is.numeric(text_size) | text_size < 0) {
        stop("Please input the right text_size (number, > 0)")
    }
    if (!is.character(text_col)) {
        stop("Please input the right text_col (character)")
    }
    # check layout
    if (!layout %in% c("nicely", "mds", "fr", "randomly", "dh", "gem", "graphopt", "grid", 
        "sphere", "kk", "lgl")) {
        stop("Please input the right layout, 'nicely','mds','fr','randomly', etc.")
    }
    # raw clu_num
    clu_info <- clu_pairs[["clu_info"]]
    clu_num <- unique(clu_info[, 2])
    clu_num1 <- clu_num[nchar(clu_num) == 1]
    clu_num2 <- clu_num[nchar(clu_num) == 2]
    clu_num1 <- clu_num1[order(clu_num1)]
    clu_num2 <- clu_num2[order(clu_num2)]
    clu_num_raw <- c(clu_num1, clu_num2)
    plot_res <- as.data.frame(table(clu_info$cluster), stringsAsFactors = F)
    colnames(plot_res) <- c("name", "size")
    plot_res$cluster <- plot_res$name
    plot_res$cluster <- factor(plot_res$cluster, levels = clu_num_raw)
    # plot_clu_col
    clu_col <- hue_pal()(length(clu_num_raw))
    # match
    if (show_type == "number") {
        if (show_sig == F) {
            plot_res_con <- clu_crosstalk_undir[, c("cluster1", "cluster2", "LR_number", "cluster1")]
            colnames(plot_res_con) <- c("from", "to", "score", "cluster")
        }
        if (show_sig == T) {
            clu_crosstalk_undir <- clu_crosstalk_undir[clu_crosstalk_undir$LR_number_type == 
                "Significant", ]
            if (nrow(clu_crosstalk_undir) == 0) {
                stop("No significant crosstalk in clu_pairs")
            }
            plot_res_con <- clu_crosstalk_undir[, c("cluster1", "cluster2", "LR_number", "cluster1")]
            colnames(plot_res_con) <- c("from", "to", "score", "cluster")
        }
        
    }
    if (show_type == "score") {
        if (show_sig == F) {
            plot_res_con <- clu_crosstalk_undir[, c("cluster1", "cluster2", "LR_score", "cluster1")]
            colnames(plot_res_con) <- c("from", "to", "score", "cluster")
        }
        if (show_sig == T) {
            clu_crosstalk_undir <- clu_crosstalk_undir[clu_crosstalk_undir$LR_score_type == 
                "Significant", ]
            if (nrow(clu_crosstalk_undir) == 0) {
                stop("No significant crosstalk in clu_pairs")
            }
            plot_res_con <- clu_crosstalk_undir[, c("cluster1", "cluster2", "LR_score", "cluster1")]
            colnames(plot_res_con) <- c("from", "to", "score", "cluster")
        }
    }
    mygraph <- graph_from_data_frame(plot_res_con, vertices = plot_res, directed = TRUE)
    ggraph(mygraph, layout = layout) + geom_edge_link(aes(edge_width = score), edge_colour = edge_col, 
        edge_alpha = edge_alpha) + geom_node_point(aes(size = size, fill = cluster), shape = 21, 
        color = "black") + scale_size_continuous(range = c(node_size_min, node_size_max)) + 
        scale_fill_manual(values = clu_col) + geom_node_text(aes(label = name), size = text_size, 
        color = text_col) + expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2)) + theme_minimal() + 
        theme(legend.position = "none", panel.grid = element_blank(), axis.line = element_blank(), 
            axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
    
}
