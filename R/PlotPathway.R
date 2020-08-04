#' Plot a pathway for a cluster
#'
#' @description Plot a pathway for a cluster showing the expression of each gene with a net plot
#' @param clu_pairs A list generated from \code{\link{FindPairs}}
#' @param pathway_name Name of pathway, e.g., 'VEGFA-VEGFR2 Pathway'
#' @param cluster Pathway in which cluster to plot, e.g., '1'
#' @param layout Layout of net plot, e.g., 'nicely','fr','mds',''randomly','dh',
#' 'gem','graphopt','grid','sphere','kk','lgl'. Default is 'nicely'
#' @param show_text_cutoff Cutoff to show text label. Default is -1, namely to show all texts
#' @param node_size Minimum size of node. Default is 10
#' @param text_size Size of text labels. Default is 3
#' @param text_col Color of text labels. Default is 'black'
#' @param edge_width Width of edge. Default is 0.5
#' @param edge_col Color of edge.. Default is 'black'
#' @param edge_alpha Transparency of edge. Default is 0.2
#' @examples
#' clu_pairs <- FindPairs(object = mouse_kidney_203_Seurat,species = 'Mouse')
#'
#' PlotPathway(clu_pairs = clu_pairs,
#'             pathway_name = 'VEGFA-VEGFR2 Pathway',
#'             cluster = '1')
#' @return A net plot of a pathway containg the expression of genes for a cluster
#' @import ggraph ggplot2
#' @importFrom scales hue_pal
#' @importFrom igraph graph_from_data_frame
#' @export PlotPathway

PlotPathway <- function(clu_pairs = NULL, pathway_name = NULL, cluster = NULL, layout = "nicely", show_text_cutoff = -1, 
    node_size = 10, text_size = 3, text_col = "black", edge_width = 0.5, edge_col = "black", edge_alpha = 0.2) {
    # check
    if (is.null(clu_pairs)) {
        stop("Please input the list from the function of FindPairs")
    }
    if (!is.list(clu_pairs)) {
        stop("Please input the right clu_pairs (generated from FindPairs)")
    }
    ndata <- clu_pairs[["ndata"]]
    clu_info <- clu_pairs[["clu_info"]]
    pathway_name_test <- grep(pattern = pathway_name, x = kegg_pathway$pathway)
    if (is.null(pathway_name) | !is.character(pathway_name) | length(pathway_name_test) == 0) {
        stop("Please input the right pathway_name (character)")
    }
    if (is.null(cluster) | !is.character(cluster) | !cluster %in% clu_info$cluster) {
        stop("Please input the right cluster (character)")
    }
    if (!is.character(layout)) {
        stop("Please input the right layout (character)")
    }
    if (!is.logical(show_clu_node)) {
        stop("Please input the right show_clu_node (logical, TRUE or FALSE)")
    }
    if (!is.numeric(node_size) | node_size < 0) {
        stop("Please input the right node_size (number, >0)")
    }
    if (!is.numeric(text_size) | text_size < 0) {
        stop("Please input the right text_size (number, >0)")
    }
    if (!is.character(text_col)) {
        stop("Please input the right text_col (character)")
    }
    if (!is.numeric(edge_width) | edge_width < 0) {
        stop("Please input the right edge_width (number, >0)")
    }
    if (!is.character(edge_col)) {
        stop("Please input the right edge_col (character)")
    }
    if (!is.numeric(edge_alpha) | edge_alpha < 0) {
        stop("Please input the right edge_alpha (number, >0)")
    }
    # check layout
    if (!layout %in% c("nicely", "mds", "fr", "randomly", "dh", "gem", "graphopt", "grid", "sphere", "kk", "lgl")) {
        stop("Please input the right layout, 'nicely','mds','fr','randomly', etc.")
    }
    kegg_pathway <- kegg_pathway[pathway_name_test, ]
    para <- clu_pairs[["para"]]
    if (para[["species"]] == "Mouse") {
        kegg_pathway$a.gn <- kegg_pathway$a.gn_mouse
        kegg_pathway$b.gn <- kegg_pathway$b.gn_mouse
    }
    kegg_pathway <- kegg_pathway[, c(1, 2)]
    kegg_pathway <- unique(kegg_pathway)
    kegg_pathway <- kegg_pathway[kegg_pathway$a.gn != "NO", ]
    kegg_pathway <- kegg_pathway[kegg_pathway$b.gn != "NO", ]
    plot_res <- data.frame(name = unique(c(kegg_pathway$a.gn, kegg_pathway$b.gn)), color = 0, stringsAsFactors = F)
    genename <- plot_res$name
    genename <- genename[genename %in% rownames(ndata)]
    ndata <- ndata[rownames(ndata) %in% genename, ]
    clu_info <- clu_info[clu_info$cluster == cluster, ]
    ndata <- ndata[, colnames(ndata) %in% clu_info$cell]
    gene_avg_exp <- apply(ndata, 1, mean)
    gene_avg_exp <- data.frame(name = names(gene_avg_exp), avg_exp = as.numeric(gene_avg_exp), stringsAsFactors = F)
    plot_res1 <- merge(plot_res, gene_avg_exp)
    plot_res <- plot_res[!plot_res$name %in% plot_res1$name, ]
    if (nrow(plot_res) > 0) {
        plot_res$avg_exp <- 0
        plot_res1 <- rbind(plot_res, plot_res1)
    }
    plot_res <- plot_res1
    plot_res_con <- plot_res1[, c("color", "name")]
    colnames(plot_res_con) <- c("from", "to")
    plot_res_con$from <- cluster
    # plot_res
    plot_res1 <- plot_res1[1, ]
    plot_res1$name <- cluster
    plot_res1$avg_exp <- max(plot_res$avg_exp)
    plot_res <- rbind(plot_res, plot_res1)
    # raw clu_num
    clu_info <- clu_pairs[["clu_info"]]
    clu_num <- unique(clu_info[, 2])
    clu_num1 <- clu_num[nchar(clu_num) == 1]
    clu_num2 <- clu_num[nchar(clu_num) == 2]
    clu_num1 <- clu_num1[order(clu_num1)]
    clu_num2 <- clu_num2[order(clu_num2)]
    clu_num_raw <- c(clu_num1, clu_num2)
    # color
    clu_col <- hue_pal()(length(clu_num_raw))
    plot_clu_col <- clu_col[which(cluster %in% clu_num_raw)]
    mygraph <- graph_from_data_frame(plot_res_con, vertices = plot_res, directed = FALSE)
    ggraph(mygraph, layout = layout) + geom_edge_link(edge_colour = edge_col, edge_alpha = edge_alpha, edge_width = edge_width) + 
        geom_node_point(aes(fill = avg_exp), shape = 21, color = "black", size = node_size) + scale_fill_gradient(low = "white", 
        high = plot_clu_col) + geom_node_text(aes(label = ifelse(avg_exp > show_text_cutoff, as.character(name), "")), 
        size = text_size, color = text_col) + expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2)) + theme_minimal() + theme(panel.grid = element_blank(), 
        axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank()) + 
        labs(title = paste0(pathway_name, " for cluster ", cluster))
}
