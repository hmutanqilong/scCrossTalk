#' Plot pathways between pairwise clusters
#'
#' @description Plot pathways between pairwise clusters with a net plot of pathways between pairwise clusters
#' @param clu_pairs A list generated from \code{\link{FindPairs}}
#' @param show_sig Whether to show significant highly expressed LR pairs between pairwise clusters. Default is FALSE
#' To show significant LR pairs, please run \code{\link{PairsSig}}
#' @param show_clu_node Whether to show cluster nodes. Default is TRUE
#' @param layout Layout of net plot, e.g., 'fr','mds',''randomly','dh',
#' 'gem','graphopt','grid','sphere','kk','lgl'. Default is 'nicely'
#' @param show_text_cutoff Cutoff to show text label. Default is 0, namely to show all texts
#' @param node_size_min Minimum size of node. Default is 5
#' @param node_size_max Maximum size of node. Default is 10
#' @param text_size Size of text labels. Default is 3
#' @param text_col Color of text labels. Default is 'black'
#' @param edge_width Width of edge. Default is 0.5
#' @param edge_col Color of edge.. Default is 'black'
#' @param edge_alpha Transparency of edge. Default is 0.2
#' @examples
#' clu_pairs <- FindPairs(object = mouse_kidney_203_Seurat,species = 'Mouse')
#'
#' PlotPathwayNet(clu_pairs = clu_pairs, show_sig = F)
#' @return A net plot of pathways between pairwise clusters.
#' L, means ligand;R, means receptor. Size of node is the mean expression of this gene in this cluster
#' @import ggraph
#' @import ggplot2
#' @importFrom scales hue_pal
#' @importFrom igraph graph_from_data_frame
#' @export PlotPathwayNet

PlotPathwayNet <- function(clu_pairs = NULL, show_sig = F, show_clu_node = T, layout = "nicely",
                         show_text_cutoff = 0, node_size_min = 5, node_size_max = 10, text_size = 3, text_col = "black",
                         edge_width = 0.5, edge_col = "black", edge_alpha = 0.2) {
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
  if (!"p_value" %in% colnames(res_pairs) & show_sig == T) {
    stop("To show significantly highly expressed LR pairs, please run PairsSig")
  }
  if (!is.logical(show_sig)) {
    stop("Please input the right show_sig (logical, TRUE or FALSE)")
  }
  if (!is.logical(show_clu_node)) {
    stop("Please input the right show_clu_node (logical, TRUE or FALSE)")
  }
  if (!is.character(layout)) {
    stop("Please input the right layout (character)")
  }
  if (!is.numeric(show_text_cutoff) | show_text_cutoff < 0) {
    stop("Please input the right show_text_cutoff (number, >=0)")
  }
  if (!is.numeric(node_size_min) | node_size_min < 0) {
    stop("Please input the right node_size_min (number, >0)")
  }
  if (!is.numeric(node_size_max) | node_size_max < 0) {
    stop("Please input the right node_size_max (number, >0)")
  }
  if (node_size_min >= node_size_max) {
    stop("Please input the right node_size_min and node_size_max (node_size_min < node_size_max)")
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
  if (!layout %in% c("nicely", "mds", "fr", "randomly", "dh", "gem", "graphopt", "grid",
                     "sphere", "kk", "lgl")) {
    stop("Please input the right layout, 'nicely','mds','fr','randomly', etc.")
  }
  res_pairs$ligand_gene_symbol <- paste0(res_pairs$ligand_gene_symbol, "(", res_pairs$ligand_clu,
                                         "L)")
  res_pairs$receptor_gene_symbol <- paste0(res_pairs$receptor_gene_symbol, "(", res_pairs$receptor_clu,
                                           "R)")
  pairs_sig <- 1:nrow(res_pairs)
  if (show_sig == T) {
    pairs_sig <- which(res_pairs$type == "Significant")
  }
  plot_res_con <- res_pairs[pairs_sig, c("ligand_gene_symbol", "receptor_gene_symbol", "ligand_clu")]
  colnames(plot_res_con) <- c("from", "to", "cluster")
  if (show_clu_node == T) {
    plot_res_con1 <- res_pairs[pairs_sig, c("ligand_clu", "ligand_gene_symbol", "ligand_clu")]
    plot_res_con1 <- unique(plot_res_con1)
    colnames(plot_res_con1) <- c("from", "to", "cluster")
    plot_res_con2 <- res_pairs[pairs_sig, c("receptor_clu", "receptor_gene_symbol", "receptor_clu")]
    plot_res_con2 <- unique(plot_res_con2)
    colnames(plot_res_con2) <- c("from", "to", "cluster")
    plot_res_con <- rbind(plot_res_con, plot_res_con1, plot_res_con2)
  }
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
  if (show_clu_node == T) {
    plot_res1 <- data.frame(name = clu_num, size = max(plot_res$size + 1), cluster = clu_num,
                            stringsAsFactors = F)
    plot_res <- rbind(plot_res, plot_res1)
  }
  plot_res$cluster <- factor(plot_res$cluster, levels = clu_num)
  # plot_clu_col
  clu_col <- hue_pal()(length(clu_num_raw))
  plot_clu_col <- NULL
  for (i in 1:length(clu_num)) {
    plot_clu_col1 <- which(clu_num_raw == clu_num[i])
    plot_clu_col1 <- clu_col[plot_clu_col1]
    plot_clu_col <- c(plot_clu_col, plot_clu_col1)
  }
  mygraph <- graph_from_data_frame(plot_res_con, vertices = plot_res, directed = FALSE)
  ggraph(mygraph, layout = layout) + geom_edge_link(edge_colour = edge_col, edge_alpha = edge_alpha,
                                                    edge_width = edge_width) + geom_node_point(aes(size = size, fill = cluster), shape = 21,
                                                                                               color = "black", alpha = 0.9) + scale_size_continuous(range = c(node_size_min, node_size_max)) +
    scale_fill_manual(values = plot_clu_col) + geom_node_text(aes(label = ifelse(size >
                                                                                   show_text_cutoff, as.character(name), "")), size = text_size, color = text_col) +
    expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2)) + theme_minimal() + theme(legend.position = "none",
                                                                                panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
                                                                                axis.text = element_blank(), axis.title = element_blank())
}
