#' Plot a LR pair between pairwise clusters
#'
#' @description Plot a LR pair between pairwise clusters with a violin plot of
#' gene expression distribution for pairwise ligand and receptor
#' @param clu_pairs A list generated from \code{\link{FindPairs}}
#' @param ligand Ligand gene name, e.g., 'App'
#' @param ligand_clu Cluster of the ligand, e.g., '1'
#' @param receptor Receptor gene name, e.g., 'Tspan12'.
#' @param receptor_clu Cluster of the receptor, e.g., '2'
#' @param show_jitter to show jitter. Default is TRUE
#' @param jitter_size Size of jitter. Default is 2
#' @examples
#' clu_pairs <- FindPairs(object = mouse_kidney_203_Seurat,species = 'Mouse')
#'
#' PlotPairsViolin(clu_pairs = clu_pairs,
#'                 ligand = 'App',
#'                 ligand_clu = '1',
#'                 receptor = 'Tspan12',
#'                 receptor_clu = '2')
#' @return A violin plot of gene expression distribution for pairwise ligand and receptor
#' @import ggplot2
#' @importFrom scales hue_pal
#' @export PlotPairsViolin

PlotPairsViolin <- function(clu_pairs = NULL, ligand = NULL, ligand_clu = NULL, receptor = NULL, 
    receptor_clu = NULL, show_jitter = T, jitter_size = 2) {
    # check clu_pairs
    if (is.null(clu_pairs)) {
        stop("Please input the list from the function of FindPairs")
    }
    if (!is.list(clu_pairs)) {
        stop("Please input the right clu_pairs (generated from FindPairs)")
    }
    if (!is.character(ligand)) {
        stop("Please input the right ligand (character)")
    }
    if (!is.character(ligand_clu)) {
        stop("Please input the right ligand_clu (character)")
    }
    if (!is.character(receptor)) {
        stop("Please input the right receptor (character)")
    }
    if (!is.character(receptor_clu)) {
        stop("Please input the right receptor_clu (character)")
    }
    if (!is.logical(show_jitter)) {
        stop("Please input the right receptor_clu (logical, TRUE or FALSE)")
    }
    if (jitter_size < 0 | !is.numeric(jitter_size)) {
        stop("Please input the right jitter_size (number, > 0)")
    }
    ndata <- clu_pairs[["ndata"]]
    clu_info <- clu_pairs[["clu_info"]]
    res_pairs <- clu_pairs[["clu_pairs"]]
    # clu_num
    clu_info <- clu_pairs[["clu_info"]]
    clu_num <- unique(clu_info[, 2])
    clu_num1 <- clu_num[nchar(clu_num) == 1]
    clu_num2 <- clu_num[nchar(clu_num) == 2]
    clu_num1 <- clu_num1[order(clu_num1)]
    clu_num2 <- clu_num2[order(clu_num2)]
    clu_num <- c(clu_num1, clu_num2)
    if (nrow(res_pairs) == 0) {
        stop("No available LR pairs in the result of clu_pairs")
    }
    # check ligand and receptor
    res_pairs1 <- paste0(res_pairs$ligand_gene_symbol, res_pairs$ligand_clu, res_pairs$receptor_gene_symbol, 
        res_pairs$receptor_clu)
    input_pair <- paste0(ligand, ligand_clu, receptor, receptor_clu)
    if (!input_pair %in% res_pairs1) {
        stop(paste0("No available LR pair in the result of clu_pairs from ", ligand_clu, " to ", 
            receptor_clu))
    }
    ndata1 <- ndata[, clu_info[clu_info$cluster == ligand_clu, ]$cell]
    ndata2 <- ndata[, clu_info[clu_info$cluster == receptor_clu, ]$cell]
    ndata1 <- as.numeric(ndata1[rownames(ndata1) == ligand, ])
    ndata2 <- as.numeric(ndata2[rownames(ndata2) == receptor, ])
    ligand <- paste0(ligand, "(L)")
    receptor <- paste0(receptor, "(R)")
    res_plot1 <- data.frame(logCount = ndata1, LR = ligand, cluster = ligand_clu, stringsAsFactors = F)
    res_plot2 <- data.frame(logCount = ndata2, LR = receptor, cluster = receptor_clu, stringsAsFactors = F)
    res_plot <- rbind(res_plot1, res_plot2)
    res_plot$LR <- factor(res_plot$LR, levels = c(ligand, receptor))
    res_plot$cluster <- factor(res_plot$cluster, levels = c(ligand_clu, receptor_clu))
    clu_col <- hue_pal()(length(clu_num))
    clu_col1 <- clu_col[which(clu_num == ligand_clu)]
    clu_col2 <- clu_col[which(clu_num == receptor_clu)]
    if (show_jitter == T) {
        jitter_alpha <- 0.5
        violin_alpha <- 0.9
    }
    if (show_jitter == F) {
        jitter_alpha <- 0
        violin_alpha <- 1
    }
    ggplot(res_plot, aes(x = LR, y = logCount, fill = cluster)) + geom_jitter(alpha = jitter_alpha, 
        size = jitter_size) + geom_violin(alpha = violin_alpha) + scale_fill_manual(values = c(clu_col1, 
        clu_col2)) + theme_bw() + theme(panel.grid = element_blank(), axis.text.x = element_text(color = "black", 
        size = 12), axis.text.y = element_text(color = "black", size = 12), axis.title.x = element_blank())
}
