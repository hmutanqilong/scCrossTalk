#' Plot a LR pair between pairwise clusters
#'
#' @description Plot a LR pair between pairwise clusters with a point plot of gene expression for pairwise ligands and receptors
#' @param clu_pairs A list generated from \code{\link{FindPairs}}
#' @param ligand Ligand gene name, e.g., 'App'.
#' @param ligand_clu Cluster of the ligand, e.g., '1'
#' @param receptor Receptor gene name, e.g., 'Tspan12'.
#' @param receptor_clu Cluster of the receptor, e.g., '2'
#' @param reduction Which dimensionality reduction to use, umap, tsne, pca.
#' @param size Size of the point. Default is 1
#' @param text_size Size of the axis x and y text.
#' @examples
#' clu_pairs <- FindPairs(object = mouse_kidney_203_Seurat,species = 'Mouse')
#'
#' PlotPairsDim(clu_pairs = clu_pairs,
#'              ligand = 'App',
#'              ligand_clu = '1',
#'              receptor = 'Tspan12',
#'              receptor_clu = '2')
#' @return A point plot of gene expression for pairwise ligands and receptors
#' @import ggplot2
#' @importFrom scales hue_pal
#' @importFrom ggpubr ggarrange
#' @importFrom dplyr filter
#' @export PlotPairsDim

PlotPairsDim <- function(clu_pairs = NULL, ligand = NULL, ligand_clu = NULL, receptor = NULL, 
    receptor_clu = NULL, reduction = "umap", size = 1, text_size = 12) {
    # check clu_pairs
    if (is.null(clu_pairs) | !is.list(clu_pairs)) {
        stop("Please input the list from the function of FindPairs")
    }
    # check ligand
    if (is.null(ligand) | !is.character(ligand)) {
        stop("Please input the right ligand (character)")
    }
    # check ligand_clu
    if (is.null(ligand_clu) | !is.character(ligand_clu)) {
        stop("Please input the right ligand_clu (character)")
    }
    # check receptor
    if (is.null(receptor) | !is.character(receptor)) {
        stop("Please input the right receptor (character)")
    }
    # check receptor_clu
    if (is.null(receptor_clu) | !is.character(receptor_clu)) {
        stop("Please input the right receptor_clu (character)")
    }
    if (!reduction %in% names(clu_pairs[["reductions"]]) | !is.character(reduction)) {
        stop("Please input the right method name of dimensionality reduction, 'umap' or 'tsne', ensuring Seurat object have done tsne or umap")
    }
    if (size < 0 | !is.numeric(size)) {
        stop("Please input the right size (number, >0)")
    }
    if (text_size < 0 | !is.numeric(text_size)) {
        stop("Please input the right text_size (number, >0)")
    }
    cell_dim <- as.data.frame(clu_pairs[["reductions"]][[reduction]]@cell.embeddings)
    # check ligand gene
    genename <- rownames(clu_pairs[["ndata"]])
    if (!ligand %in% genename) {
        stop(paste0(ligand, " is not in revised ndata. Please input the right gene name of ligand"))
    }
    if (!receptor %in% genename) {
        stop(paste0(receptor, " is not in revised ndata. Please input the right gene name of ligand"))
    }
    # obtain clusters
    clu_info <- clu_pairs[["clu_info"]]
    clu_num <- unique(clu_info[, 2])
    clu_num1 <- clu_num[nchar(clu_num) == 1]
    clu_num2 <- clu_num[nchar(clu_num) == 2]
    clu_num1 <- clu_num1[order(clu_num1)]
    clu_num2 <- clu_num2[order(clu_num2)]
    clu_num <- c(clu_num1, clu_num2)
    # check ligand_clu and receptor_clu
    if (!ligand_clu %in% clu_num) {
        stop(paste0(ligand_clu, " is not an effective cluster. Please input the right cluster of ligand"))
    }
    if (!receptor_clu %in% clu_num) {
        stop(paste0(receptor_clu, " is not an effective cluster. Please input the right cluster of receptor"))
    }
    # create data frame for plot
    cell_dim <- cbind(cell_dim, clu_info)
    ndata <- clu_pairs[["ndata"]]
    ndata1 <- ndata[rownames(ndata) == ligand, ]
    ndata2 <- ndata[rownames(ndata) == receptor, ]
    cell_dim$ligand <- as.numeric(ndata1)
    cell_dim$receptor <- as.numeric(ndata2)
    ligand_col <- hue_pal()(length(clu_num))[which(ligand_clu == clu_num)]
    receptor_col <- hue_pal()(length(clu_num))[which(receptor_clu == clu_num)]
    # plotting
    p_ligand <- ggplot(cell_dim, aes(x = UMAP_1, y = UMAP_2)) + geom_point(data = filter(cell_dim, 
        cluster == ligand_clu), aes(col = ligand), size = size) + scale_colour_gradient(low = "black", 
        high = ligand_col) + geom_point(data = filter(cell_dim, cluster != ligand_clu), col = "grey", 
        size = size) + labs(title = paste0("Cluster ", ligand_clu, ": ", ligand, "(ligand)")) + 
        theme_bw() + theme(panel.grid = element_blank()) + theme(plot.title = element_text(hjust = 0.5), 
        legend.justification = c(1, 1), axis.text.x = element_text(size = text_size, color = "black"), 
        axis.text.y = element_text(size = text_size, color = "black")) + labs(col = "logCount")
    p_receptor <- ggplot(cell_dim, aes(x = UMAP_1, y = UMAP_2)) + geom_point(data = filter(cell_dim, 
        cluster == receptor_clu), aes(col = receptor), size = size) + scale_colour_gradient(low = "black", 
        high = receptor_col) + geom_point(data = filter(cell_dim, cluster != receptor_clu), 
        col = "grey", size = size) + labs(title = paste0("Cluster ", receptor_clu, ": ", receptor, 
        "(receptor)")) + theme_bw() + theme(panel.grid = element_blank()) + theme(plot.title = element_text(hjust = 0.5), 
        legend.justification = c(1, 1), axis.text.x = element_text(size = text_size, color = "black"), 
        axis.text.y = element_text(size = text_size, color = "black")) + labs(col = "logCount")
    ggarrange(p_ligand, p_receptor)
}
