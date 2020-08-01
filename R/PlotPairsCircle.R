#' Plot LR pairs between pairwise clusters
#'
#' @description Plot LR pairs between pairwise clusters with a circle plot of LR pairs between pairwise clusters
#' @param clu_pairs A list generated from \code{\link{FindPairs}}
#' @param show_sig Whether to show significant highly expressed LR pairs between pairwise clusters. Default is FALSE
#' To show significant LR pairs, please run \code{\link{PairsSig}}
#' @param ligand_clu Cluster of the ligand, e.g., '1'
#' @param receptor_clu Cluster of the receptor, e.g., '2'
#' @examples
#' clu_pairs <- FindPairs(object = mouse_kidney_203_Seurat,species = 'Mouse')
#'
#' PlotPairsCircle(clu_pairs = clu_pairs,
#'                 show_sig = T,
#'                 ligand_clu = '1',
#'                 receptor_clu = '2')
#' @return A circle plot of LR pairs between pairwise clusters
#' @import ggplot2
#' @importFrom scales hue_pal
#' @importFrom grDevices colorRampPalette
#' @export PlotPairsViolin

PlotPairsCircle <- function(clu_pairs = NULL, show_sig = F, ligand_clu = NULL, receptor_clu = NULL){
  # check clu_pairs
  if (is.null(clu_pairs)) {
    stop("Please input the list from the function of FindPairs")
  }
  if (!is.list(clu_pairs)) {
    stop("Please input the right clu_pairs (generated from FindPairs)")
  }
  if (!is.logical(show_sig)) {
    stop("Please input the right show_sig (logical, TRUE or FALSE)")
  }
  if (!is.character(ligand_clu)) {
    stop("Please input the right ligand_clu (character)")
  }
  if (!is.character(receptor_clu)) {
    stop("Please input the right receptor_clu (character)")
  }
  res_pairs <- clu_pairs[["clu_pairs"]]
  if (nrow(res_pairs) == 0) {
    stop("There is no sigficantly highly expressed LR in clu_pairs")
  }
  if (!"p_value" %in% colnames(res_pairs) & show_sig == T) {
    stop("To show significantly highly expressed LR pairs, please run PairsSig")
  }
  # raw clu_num
  clu_info <- clu_pairs[["clu_info"]]
  clu_num <- unique(clu_info[, 2])
  clu_num1 <- clu_num[nchar(clu_num) == 1]
  clu_num2 <- clu_num[nchar(clu_num) == 2]
  clu_num1 <- clu_num1[order(clu_num1)]
  clu_num2 <- clu_num2[order(clu_num2)]
  clu_num_raw <- c(clu_num1, clu_num2)
  if (!ligand_clu %in% clu_num_raw) {
    stop(paste0(ligand_clu, " is not an effective cluster. Please input the right cluster of ligand"))
  }
  if (!receptor_clu %in% clu_num_raw) {
    stop(paste0(receptor_clu, " is not an effective cluster. Please input the right cluster of receptor"))
  }
  res_pairs <- res_pairs[res_pairs$ligand_clu == ligand_clu & res_pairs$receptor_clu == receptor_clu,]
  if (nrow(res_pairs) == 0) {
    stop(paste0("There is no LR pairs in clu_pairs from cluster from ", ligand_clu,
                ' to ', receptor_clu))
  }
  if (show_sig == T) {
    res_pairs <- res_pairs[res_pairs$type == 'Significant',]
    if (nrow(res_pairs) == 0) {
      stop(paste0("There is no significant LR pairs in clu_pairs from cluster from ", ligand_clu,
                  ' to ', receptor_clu))
    }
  }
  LR_score<- res_pairs$lr_score
  LR_score_min<- min(LR_score)
  LR_score_len<-  max(LR_score)-min(LR_score)
  for (i in 1:length(LR_score)) {
    LR_score[i] <- (LR_score[i]-LR_score_min)*2/LR_score_len+1
  }
  res_pairs<- res_pairs[,c("ligand_gene_symbol","receptor_gene_symbol")]
  clu_col <- hue_pal()(length(clu_num_raw))
  clu_col1 <- clu_col[which(clu_num == ligand_clu)]
  clu_col2 <- clu_col[which(clu_num == receptor_clu)]
  ligand_col<- rep(clu_col1,nrow(res_pairs))
  names(ligand_col) <- res_pairs$ligand_gene_symbol
  receptor_col<- rep(clu_col2,nrow(res_pairs))
  names(receptor_col) <- res_pairs$receptor_gene_symbol
  clu_col <- c(ligand_col,receptor_col)
  chordDiagramFromDataFrame(res_pairs, annotationTrack = "grid", preAllocateTracks = 1, directional = 1,
                            direction.type = "arrows",
                            link.arr.length = 0.12,
                            link.arr.width = 0.12,link.arr.type = 'ellipse',
                            link.arr.lty = par("lty"),
                            link.arr.lwd = LR_score, link.arr.col = 'royalblue',
                            grid.col = clu_col,col= 'lightblue', big.gap = 5, small.gap = 0.2)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=ifelse(nrow(res_pairs)>50, yes = 0.5,no=1))
    circos.axis(h = "top",labels = FALSE,minor.ticks = FALSE, major.at = c(xlim), major.tick.length = 2, sector.index = sector.name, track.index = 2)
  }, bg.border = NA)
}
