#' Plot crosstalk between pairwise clusters
#'
#' @description Plot crosstalk between pairwise clusters by sum of LR pairs number or score between pairwise clusters with Sankey plot
#' @param clu_pairs A list generated from \code{\link{FindPairs}}
#' @param show_type Which to show, sum of LR pairs number or LR score.
#' Input 'number' to show sum of LR pairs number; Input 'score' to show sum of LR pairs score
#' @param show_sig Whether to show significant crosstalk between pairwise clusters. Default is FALSE
#' To show significant crosstalk, please run \code{\link{CrossTalkSig}}
#' @examples
#' clu_pairs <- FindPairs(object = mouse_kidney_203_Seurat, species = 'Mouse')
#'
#' PlotCrossTalkSan(clu_pairs = clu_pairs, show_type = 'number')
#'
#' PlotCrossTalkSan(clu_pairs = clu_pairs, show_type = 'score', show_sig = T)
#' @return Sankey plot of LR pairs between pairwise clusters
#' @import ggalluvial ggplot2
#' @importFrom scales hue_pal
#' @export PlotCrossTalkSan

PlotCrossTalkSan <- function(clu_pairs = NULL, show_type = NULL, show_sig = F) {
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
            plot_res$group <- 1:nrow(plot_res)
            plot_res1 <- plot_res[, c("cluster1_source", "LR_number", "group")]
            plot_res2 <- plot_res[, c("cluster2_target", "LR_number", "group")]
            plot_res1$Crosstalk <- "source(ligand)"
            plot_res2$Crosstalk <- "target(receptor)"
            plot_res1 <- plot_res1[, c(2, 3, 4, 1)]
            colnames(plot_res1)[4] <- "cluster"
            plot_res2 <- plot_res2[, c(2, 3, 4, 1)]
            colnames(plot_res2)[4] <- "cluster"
            plot_res <- rbind(plot_res1, plot_res2)
        }
        if (show_sig == T) {
            clu_crosstalk <- clu_crosstalk[clu_crosstalk$LR_number_type == "Significant", ]
            if (nrow(clu_crosstalk) == 0) {
                stop("No significant crosstalk in clu_pairs")
            }
            plot_res <- clu_crosstalk[, c("cluster1_source", "cluster2_target", "LR_number")]
            plot_res$group <- 1:nrow(plot_res)
            plot_res1 <- plot_res[, c("cluster1_source", "LR_number", "group")]
            plot_res2 <- plot_res[, c("cluster2_target", "LR_number", "group")]
            plot_res1$Crosstalk <- "source(ligand)"
            plot_res2$Crosstalk <- "target(receptor)"
            plot_res1 <- plot_res1[, c(2, 3, 4, 1)]
            colnames(plot_res1)[4] <- "cluster"
            plot_res2 <- plot_res2[, c(2, 3, 4, 1)]
            colnames(plot_res2)[4] <- "cluster"
            plot_res <- rbind(plot_res1, plot_res2)
        }
        
    }
    if (show_type == "score") {
        if (show_sig == F) {
            plot_res <- clu_crosstalk[, c("cluster1_source", "cluster2_target", "LR_score")]
            plot_res$group <- 1:nrow(plot_res)
            plot_res1 <- plot_res[, c("cluster1_source", "LR_score", "group")]
            plot_res2 <- plot_res[, c("cluster2_target", "LR_score", "group")]
            plot_res1$Crosstalk <- "source(ligand)"
            plot_res2$Crosstalk <- "target(receptor)"
            plot_res1 <- plot_res1[, c(2, 3, 4, 1)]
            colnames(plot_res1)[4] <- "cluster"
            plot_res2 <- plot_res2[, c(2, 3, 4, 1)]
            colnames(plot_res2)[4] <- "cluster"
            plot_res <- rbind(plot_res1, plot_res2)
        }
        if (show_sig == T) {
            clu_crosstalk <- clu_crosstalk[clu_crosstalk$LR_score_type == "Significant", ]
            if (nrow(clu_crosstalk) == 0) {
                stop("No significant crosstalk in clu_pairs")
            }
            plot_res <- clu_crosstalk[, c("cluster1_source", "cluster2_target", "LR_score")]
            plot_res$group <- 1:nrow(plot_res)
            plot_res1 <- plot_res[, c("cluster1_source", "LR_score", "group")]
            plot_res2 <- plot_res[, c("cluster2_target", "LR_score", "group")]
            plot_res1$Crosstalk <- "source(ligand)"
            plot_res2$Crosstalk <- "target(receptor)"
            plot_res1 <- plot_res1[, c(2, 3, 4, 1)]
            colnames(plot_res1)[4] <- "cluster"
            plot_res2 <- plot_res2[, c(2, 3, 4, 1)]
            colnames(plot_res2)[4] <- "cluster"
            plot_res <- rbind(plot_res1, plot_res2)
        }
    }
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
    # factor Crosstalk
    plot_res$Crosstalk <- factor(plot_res$Crosstalk, levels = c("source(ligand)", "target(receptor)"))
    # plot_clu_col
    clu_col <- hue_pal()(length(clu_num_raw))
    plot_clu_col <- NULL
    for (i in 1:length(clu_num)) {
        plot_clu_col1 <- which(clu_num_raw == clu_num[i])
        plot_clu_col1 <- clu_col[plot_clu_col1]
        plot_clu_col <- c(plot_clu_col, plot_clu_col1)
    }
    y_title <- paste0("Sum of LR pairs ", show_type)
    colnames(plot_res)[1] <- "LR"
    ggplot(plot_res, aes(x = Crosstalk, y = LR, stratum = cluster, alluvium = group, fill = cluster, label = cluster)) + 
        geom_flow(alpha = 0.7, width = 0.25, color = "darkgray") + geom_stratum(alpha = 1, width = 0.25) + geom_text(stat = "stratum", 
        size = 3.5, angle = 0) + scale_fill_manual(values = plot_clu_col) + ylab(y_title) + theme_test() + theme(legend.position = "none", 
        axis.text.y = element_blank(), axis.line = element_blank(), axis.ticks = element_blank())
}
