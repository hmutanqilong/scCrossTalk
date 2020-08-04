#' Plot pathways of target cluster from source cluster
#'
#' @description Plot pathways of target cluster from source cluster with a net plot
#' @param clu_pairs A list generated from \code{\link{FindPairs}}
#' @param show_sig Whether to show significant highly expressed LR pairs between pairwise clusters. Default is FALSE
#' To show significant LR pairs, please run \code{\link{PairsSig}}
#' @param layout Layout of net plot, e.g., 'fr','mds',''randomly','dh',
#' 'gem','graphopt','grid','sphere','kk','lgl'. Default is 'nicely'
#' @param ligand_clu Cluster of the ligand, e.g., '1'
#' @param ligand Ligand gene name, e.g., 'App'
#' @param receptor_clu Cluster of the receptor, e.g., '2'
#' @param receptor Receptor gene name, e.g., 'Notch2'.
#' @param top_path_number Number of top pathways for ligand or receptor. Default is 3
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
#' PlotPathwayNet(clu_pairs = clu_pairs,
#'                ligand_clu = '1',
#'                receptor_clu = '2')
#'
#' PlotPathwayNet(clu_pairs = clu_pairs,
#'                ligand_clu = '1',
#'                ligand = 'App',
#'                receptor_clu = '2',
#'                receptor = 'all')
#'
#' PlotPathwayNet(clu_pairs = clu_pairs,
#'                ligand_clu = '1',
#'                ligand = 'App',
#'                receptor_clu = '2',
#'                receptor = 'Notch2')
#'
#' PlotPathwayNet(clu_pairs = clu_pairs,
#'                ligand_clu = '1',
#'                ligand = 'App',
#'                receptor_clu = '2',
#'                receptor = c('Notch2','Tspan12'))
#'
#' @return A net plot of pathways between pairwise clusters.
#' L, means ligand;R, means receptor. Size of node is the mean expression of this gene in this cluster
#' @import ggraph ggplot2
#' @importFrom scales hue_pal
#' @importFrom igraph graph_from_data_frame
#' @export PlotPathwayNet

PlotPathwayNet <- function(clu_pairs = NULL, show_sig = F, layout = "kk", ligand_clu = NULL, ligand = "all", receptor_clu = NULL, 
    receptor = "all", top_path_number = 3, show_text_cutoff = 0, node_size_min = 5, node_size_max = 10, text_size = 3, 
    text_col = "black", edge_width = 0.5, edge_col = "black", edge_alpha = 0.2) {
    # check
    if (is.null(clu_pairs)) {
        stop("Please input the list from the function of FindPairs")
    }
    if (!is.list(clu_pairs)) {
        stop("Please input the right clu_pairs (generated from FindPairs)")
    }
    res_pairs <- clu_pairs[["clu_pairs"]]
    LR_pathway <- clu_pairs[["LR_pathway"]]
    if (nrow(res_pairs) == 0) {
        stop("There is no highly expressed LR in clu_pairs")
    }
    if (nrow(LR_pathway) == 0) {
        stop("There is related pathways for clusters")
    }
    if (!"p_value" %in% colnames(res_pairs) & show_sig == T) {
        stop("To show significantly highly expressed LR pairs, please run PairsSig")
    }
    if (!is.logical(show_sig)) {
        stop("Please input the right show_sig (logical, TRUE or FALSE)")
    }
    if (!is.character(layout)) {
        stop("Please input the right layout (character)")
    }
    if (!is.character(ligand_clu) | !ligand_clu %in% res_pairs$ligand_clu) {
        stop("Please input the right ligand_clu (character)")
    }
    if (ligand[1] != "all") {
        if (!is.character(ligand[1]) | !ligand[1] %in% res_pairs$ligand_gene_symbol) {
            stop("Please input the right ligand (character)")
        }
    }
    if (!is.character(receptor_clu) | !receptor_clu %in% res_pairs$receptor_clu) {
        stop("Please input the right receptor_clu (character)")
    }
    if (receptor[1] != "all") {
        if (!is.character(receptor[1]) | !receptor[1] %in% res_pairs$receptor_gene_symbol) {
            stop("Please input the right receptor (character)")
        }
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
    if (!layout %in% c("nicely", "mds", "fr", "randomly", "dh", "gem", "graphopt", "grid", "sphere", "kk", "lgl")) {
        stop("Please input the right layout, 'nicely','mds','fr','randomly', etc.")
    }
    res_pairs <- res_pairs[res_pairs$ligand_clu == ligand_clu & res_pairs$receptor_clu == receptor_clu, ]
    if (nrow(res_pairs) == 0) {
        stop(paste0("No available LR pairs from cluster ", ligand_clu, " to cluster ", receptor_clu))
    }
    pairs_sig <- 1:nrow(res_pairs)
    if (show_sig == T) {
        pairs_sig <- which(res_pairs$type == "Significant")
    }
    if (length(pairs_sig) == 0) {
        stop(paste0("No significant LR pairs from cluster ", ligand_clu, " to cluster ", receptor_clu))
    }
    res_pairs <- res_pairs[pairs_sig, ]
    if (ligand[1] != "all") {
        res_pairs <- res_pairs[res_pairs$ligand_gene_symbol %in% ligand, ]
        if (nrow(res_pairs) == 0) {
            stop(paste0("No available LR pairs from cluster ", ligand_clu, " to cluster ", receptor_clu, " for ligand ", 
                ligand))
        }
    }
    if (receptor[1] != "all") {
        res_pairs <- res_pairs[res_pairs$receptor_gene_symbol %in% receptor, ]
        if (nrow(res_pairs) == 0) {
            stop(paste0("No available LR pairs from cluster ", ligand_clu, " to cluster ", receptor_clu, " for receptor ", 
                receptor))
        }
    }
    # plot_res_con
    plot_res_con1 <- res_pairs[, c("ligand_clu", "ligand_gene_symbol")]
    plot_res_con1 <- unique(plot_res_con1)
    plot_res_con2 <- res_pairs[, c("ligand_gene_symbol", "receptor_gene_symbol")]
    plot_res_con2 <- unique(plot_res_con2)
    plot_res_con3 <- res_pairs[, c("receptor_clu", "receptor_gene_symbol")]
    plot_res_con3 <- unique(plot_res_con3)
    colnames(plot_res_con1) <- c("from", "to")
    colnames(plot_res_con2) <- c("from", "to")
    colnames(plot_res_con3) <- c("from", "to")
    plot_res_con <- rbind(plot_res_con1, plot_res_con2, plot_res_con3)
    LR_pathway1 <- LR_pathway[LR_pathway$type == "ligand", ]
    LR_pathway1 <- LR_pathway1[LR_pathway1$gene %in% res_pairs$ligand_gene_symbol, ]
    if (nrow(LR_pathway1) == 0) {
        stop(paste0("No available pathways for ligands in cluster ", ligand_clu))
    }
    genename <- unique(LR_pathway1$gene)
    LR_pathway_tmp <- LR_pathway1[1, ]
    for (i in 1:length(genename)) {
        genename1 <- genename[i]
        LR_pathway_tmp1 <- LR_pathway1[LR_pathway1$gene == genename1, ]
        if (nrow(LR_pathway_tmp1) > top_path_number) {
            LR_pathway_tmp1 <- LR_pathway_tmp1[c(1:top_path_number), ]
            LR_pathway_tmp1 <- LR_pathway_tmp1[order(-LR_pathway_tmp1$Freq), ]
        }
        LR_pathway_tmp <- rbind(LR_pathway_tmp, LR_pathway_tmp1)
    }
    LR_pathway1 <- LR_pathway_tmp[-1, ]
    LR_pathway1 <- LR_pathway1[, c("gene", "pathway")]
    colnames(LR_pathway1) <- c("from", "to")
    LR_pathway1$to <- paste0(LR_pathway1$to, "(", ligand_clu, ")")
    LR_pathway2 <- LR_pathway[LR_pathway$type == "receptor", ]
    LR_pathway2 <- LR_pathway2[LR_pathway2$gene %in% res_pairs$receptor_gene_symbol, ]
    if (nrow(LR_pathway2) == 0) {
        stop(paste0("No available pathways for receptors in cluster ", receptor_clu))
    }
    genename <- unique(LR_pathway2$gene)
    LR_pathway_tmp <- LR_pathway2[1, ]
    for (i in 1:length(genename)) {
        genename1 <- genename[i]
        LR_pathway_tmp1 <- LR_pathway2[LR_pathway2$gene == genename1, ]
        if (nrow(LR_pathway_tmp1) > top_path_number) {
            LR_pathway_tmp1 <- LR_pathway_tmp1[c(1:top_path_number), ]
            LR_pathway_tmp1 <- LR_pathway_tmp1[order(-LR_pathway_tmp1$Freq), ]
        }
        LR_pathway_tmp <- rbind(LR_pathway_tmp, LR_pathway_tmp1)
    }
    LR_pathway2 <- LR_pathway_tmp[-1, ]
    LR_pathway2 <- LR_pathway2[, c("gene", "pathway")]
    colnames(LR_pathway2) <- c("from", "to")
    LR_pathway2$to <- paste0(LR_pathway2$to, "(", receptor_clu, ")")
    plot_res_con <- rbind(plot_res_con, LR_pathway1, LR_pathway2)
    LR_pathway1$from <- ligand_clu
    LR_pathway2$from <- receptor_clu
    plot_res_con <- rbind(plot_res_con, LR_pathway1, LR_pathway2)
    # plot_res
    plot_res1 <- res_pairs[, c("ligand_gene_symbol", "ligand_exp_avg", "ligand_clu")]
    plot_res1 <- unique(plot_res1)
    plot_res2 <- res_pairs[, c("receptor_gene_symbol", "receptor_exp_avg", "receptor_clu")]
    plot_res2 <- unique(plot_res2)
    colnames(plot_res1) <- c("name", "size", "cluster")
    colnames(plot_res2) <- c("name", "size", "cluster")
    plot_res <- rbind(plot_res1, plot_res2)
    plot_res1$name <- plot_res1$cluster
    plot_res1$size <- max(plot_res$size) + 2
    plot_res1 <- unique(plot_res1)
    plot_res2$name <- plot_res2$cluster
    plot_res2$size <- max(plot_res$size) + 2
    plot_res2 <- unique(plot_res2)
    plot_res <- rbind(plot_res, plot_res1, plot_res2)
    LR_pathway1 <- LR_pathway[LR_pathway$type == "ligand", ]
    LR_pathway1 <- LR_pathway1[LR_pathway1$gene %in% res_pairs$ligand_gene_symbol, ]
    LR_pathway1 <- LR_pathway1[LR_pathway1$gene %in% res_pairs$ligand_gene_symbol, ]
    if (nrow(LR_pathway1) == 0) {
        stop(paste0("No available pathways for ligands in cluster ", ligand_clu))
    }
    genename <- unique(LR_pathway1$gene)
    LR_pathway_tmp <- LR_pathway1[1, ]
    for (i in 1:length(genename)) {
        genename1 <- genename[i]
        LR_pathway_tmp1 <- LR_pathway1[LR_pathway1$gene == genename1, ]
        if (nrow(LR_pathway_tmp1) > top_path_number) {
            LR_pathway_tmp1 <- LR_pathway_tmp1[c(1:top_path_number), ]
            LR_pathway_tmp1 <- LR_pathway_tmp1[order(-LR_pathway_tmp1$Freq), ]
        }
        LR_pathway_tmp <- rbind(LR_pathway_tmp, LR_pathway_tmp1)
    }
    LR_pathway1 <- LR_pathway_tmp[-1, ]
    LR_pathway1 <- LR_pathway1[, c(3, 4, 2)]
    colnames(LR_pathway1) <- c("name", "size", "cluster")
    LR_pathway1$size <- max(plot_res$size)
    LR_pathway1$cluster <- ligand_clu
    LR_pathway1$name <- paste0(LR_pathway1$name, "(", LR_pathway1$cluster, ")")
    LR_pathway2 <- LR_pathway[LR_pathway$type == "receptor", ]
    LR_pathway2 <- LR_pathway2[LR_pathway2$gene %in% res_pairs$receptor_gene_symbol, ]
    if (nrow(LR_pathway2) == 0) {
        stop(paste0("No available pathways for receptors in cluster ", receptor_clu))
    }
    genename <- unique(LR_pathway2$gene)
    LR_pathway_tmp <- LR_pathway2[1, ]
    for (i in 1:length(genename)) {
        genename1 <- genename[i]
        LR_pathway_tmp1 <- LR_pathway2[LR_pathway2$gene == genename1, ]
        if (nrow(LR_pathway_tmp1) > top_path_number) {
            LR_pathway_tmp1 <- LR_pathway_tmp1[c(1:top_path_number), ]
            LR_pathway_tmp1 <- LR_pathway_tmp1[order(-LR_pathway_tmp1$Freq), ]
        }
        LR_pathway_tmp <- rbind(LR_pathway_tmp, LR_pathway_tmp1)
    }
    LR_pathway2 <- LR_pathway_tmp[-1, ]
    LR_pathway2 <- LR_pathway2[, c(3, 4, 2)]
    colnames(LR_pathway2) <- c("name", "size", "cluster")
    LR_pathway2$size <- max(plot_res$size)
    LR_pathway2$cluster <- receptor_clu
    LR_pathway2$name <- paste0(LR_pathway2$name, "(", LR_pathway2$cluster, ")")
    plot_res <- rbind(plot_res, LR_pathway1, LR_pathway2)
    # raw clu_num
    clu_info <- clu_pairs[["clu_info"]]
    clu_num <- unique(clu_info[, 2])
    clu_num1 <- clu_num[nchar(clu_num) == 1]
    clu_num2 <- clu_num[nchar(clu_num) == 2]
    clu_num1 <- clu_num1[order(clu_num1)]
    clu_num2 <- clu_num2[order(clu_num2)]
    clu_num_raw <- c(clu_num1, clu_num2)
    # factor cluster
    clu_num <- unique(plot_res$cluster)
    clu_num <- clu_num[clu_num %in% clu_num_raw]
    clu_num1 <- clu_num[nchar(clu_num) == 1]
    clu_num2 <- clu_num[nchar(clu_num) == 2]
    clu_num1 <- clu_num1[order(clu_num1)]
    clu_num2 <- clu_num2[order(clu_num2)]
    clu_num <- c(clu_num1, clu_num2)
    plot_res$cluster <- factor(plot_res$cluster, levels = clu_num)
    # plot_clu_col
    clu_col <- hue_pal()(length(clu_num_raw))
    plot_clu_col <- NULL
    for (i in 1:length(clu_num)) {
        plot_clu_col1 <- which(clu_num_raw == clu_num[i])
        plot_clu_col1 <- clu_col[plot_clu_col1]
        plot_clu_col <- c(plot_clu_col, plot_clu_col1)
    }
    plot_res <- unique(plot_res)
    mygraph <- graph_from_data_frame(plot_res_con, vertices = plot_res, directed = FALSE)
    ggraph(mygraph, layout = layout) + geom_edge_link(edge_colour = edge_col, edge_alpha = edge_alpha, edge_width = edge_width) + 
        geom_node_point(aes(size = size, fill = cluster), shape = 21, color = "black", alpha = 0.9) + scale_size_continuous(range = c(node_size_min, 
        node_size_max)) + scale_fill_manual(values = plot_clu_col) + geom_node_text(aes(label = ifelse(size > show_text_cutoff, 
        as.character(name), "")), size = text_size, color = text_col) + expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2)) + 
        theme_minimal() + theme(legend.position = "none", panel.grid = element_blank(), axis.line = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
}
