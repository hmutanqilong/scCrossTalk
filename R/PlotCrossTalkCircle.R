#' Plot crosstalk between pairwise clusters
#'
#' @description Plot crosstalk between pairwise clusters by sum of LR pairs number or score between pairwise clusters with Circle plot
#' @param clu_pairs A list generated from \code{\link{FindPairs}}
#' @param show_type Which to show, sum of LR pairs number or LR score.
#' Input 'number' to show sum of LR pairs number; Input 'score' to show sum of LR pairs score
#' @param show_sig Whether to show significant crosstalk between pairwise clusters. Default is FALSE
#' To show significant crosstalk, please run \code{\link{CrossTalkSig}}
#' @examples
#' clu_pairs <- FindPairs(object = mouse_kidney_203_Seurat, species = 'Mouse')
#'
#' PlotCrossTalkCircle(clu_pairs = clu_pairs, show_type = 'number')
#'
#' PlotCrossTalkCircle(clu_pairs = clu_pairs, show_type = 'score', show_sig = F)
#' @return Circle plot of LR pairs between pairwise clusters
#' @import circlize
#' @importFrom scales hue_pal
#' @importFrom graphics par
#' @export PlotCrossTalkCircle

PlotCrossTalkCircle <- function(clu_pairs = NULL, show_type = NULL, show_sig = F) {
    # check clu_pairs
    if (is.null(clu_pairs) | !is.list(clu_pairs)) {
        stop("Please input the list from the function of FindPairs")
    }
    if (!is.data.frame(clu_pairs[["clu_crosstalk"]])) {
        stop("No available clu_crosstalk in clu_pairs")
    }
    clu_crosstalk <- clu_pairs[["clu_crosstalk"]]
    if (!"LR_number_pvalue" %in% colnames(clu_crosstalk) & show_sig == T) {
        stop("To show significant crosstalk, please run CrossTalkSig")
    }
    if (!is.character(show_type)) {
        stop("Please input the right show_type (character)")
    }
    if (!is.logical(show_sig)) {
        stop("Please input the right show_sig (logical, TRUE or FALSE)")
    }
    if (show_type != "number" & show_type != "score") {
        stop("Please input the right show_type!")
    }
    if (show_type == "number") {
        if (show_sig == F) {
            plot_res <- clu_crosstalk[, c("cluster1_source", "cluster2_target", "LR_number")]
        }
        if (show_sig == T) {
            clu_crosstalk <- clu_crosstalk[clu_crosstalk$LR_number_type == "Significant", ]
            if (nrow(clu_crosstalk) == 0) {
                stop("No significant crosstalk in clu_pairs")
            }
            plot_res <- clu_crosstalk[, c("cluster1_source", "cluster2_target", "LR_number")]
        }

    }
    if (show_type == "score") {
        if (show_sig == F) {
            plot_res <- clu_crosstalk[, c("cluster1_source", "cluster2_target", "LR_score")]
        }
        if (show_sig == T) {
            clu_crosstalk <- clu_crosstalk[clu_crosstalk$LR_score_type == "Significant", ]
            if (nrow(clu_crosstalk) == 0) {
                stop("No significant crosstalk in clu_pairs")
            }
            plot_res <- clu_crosstalk[, c("cluster1_source", "cluster2_target", "LR_score")]
        }
    }
    colnames(plot_res) <- c("from", "to", "value")
    # raw clu_num
    clu_info <- clu_pairs[["clu_info"]]
    clu_num <- unique(clu_info[, 2])
    clu_num1 <- clu_num[nchar(clu_num) == 1]
    clu_num2 <- clu_num[nchar(clu_num) == 2]
    clu_num1 <- clu_num1[order(clu_num1)]
    clu_num2 <- clu_num2[order(clu_num2)]
    clu_num_raw <- c(clu_num1, clu_num2)
    # factor cluster
    clu_num <- unique(plot_res$from)
    clu_num <- clu_num[clu_num %in% clu_num_raw]
    clu_num1 <- clu_num[nchar(clu_num) == 1]
    clu_num2 <- clu_num[nchar(clu_num) == 2]
    clu_num1 <- clu_num1[order(clu_num1)]
    clu_num2 <- clu_num2[order(clu_num2)]
    clu_num <- c(clu_num1, clu_num2)
    plot_res$from <- factor(plot_res$from, levels = clu_num)
    # plot_clu_col
    clu_col <- hue_pal()(length(clu_num_raw))
    plot_clu_col <- NULL
    for (i in 1:length(clu_num)) {
        plot_clu_col1 <- which(clu_num_raw == clu_num[i])
        plot_clu_col1 <- clu_col[plot_clu_col1]
        plot_clu_col <- c(plot_clu_col, plot_clu_col1)
    }
    chordDiagram(x = plot_res, grid.col = plot_clu_col, preAllocateTracks = 1, transparency = 0.25, directional = 1,
        direction.type = c("arrows", "diffHeight"), diffHeight = -0.04, annotationTrack = "grid", link.arr.type = "big.arrow",
        link.sort = TRUE, link.largest.ontop = TRUE)
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        sector.name = get.cell.meta.data("sector.index")
        circos.text(mean(xlim), ylim[1] + 0.1, sector.name, facing = "inside", niceFacing = FALSE, adj = c(0.5, 0), cex = 1.5)
        circos.axis(h = "top", labels = FALSE, major.tick = FALSE, major.tick.percentage = 0.1, sector.index = sector.name,
            track.index = 2)
    }, bg.border = NA)
}
