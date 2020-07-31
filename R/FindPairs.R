#' Find highly expressed ligand-receptor pairs and cell-cell communications for clusters
#'
#' @description Find highly expressed ligands and receptors between pairwise clusters using Z score for a \code{Seurat}
#' object (>= 3.0.0) after \code{log1p} normalization, \code{cluster analysis}
#' and \code{tSNE} or \code{Umap} dimensionality reduction, followed by inferring significant
#' LR pairs and cell-cell communication with permutation test.
#' @param object Seurat object (>= 3.0.0) after log1p normalization and cluster
#' analysis. Please ensure data is log1p normalized data and data has been clustered
#' before running scCrossTalk pipeline.
#' @param species Species of cells. The species must be defined. \code{'Human'} or \code{'Mouse'}.
#' @param use_LRdb Which ligand-receptor pairs database to ues. Default is
#' \code{'LRdb'} using CellTalkDB, namely 3,398 human LR pairs for \code{'Human'};
#' 2,033 mouse LR pairs for \code{'Mouse'}. See \code{Details} for more available
#' LR databases or a custom LR database.
#' @param revise_gene Whether to revise genes according to NCBI Gene symbols
#' (updated in April 28, 2020). Default is \code{TRUE}. When using a custom LR
#' database for \code{use_LRdb}, please be cautious to decide whether to revise genes
#' @param cell_min_pct Include the gene detected in at least this many cells
#' in each cluster. Default is \code{0.25}
#' @param p_value Include the significantly highly expressed ligands and receptors
#' with this cutoff of p value from Z score. Default is \code{0.05}
#' @return A list include a new data matrix, a data frame containing highly expressed
#' ligand-receptor pairs, and two data frames containing sum of LR pairs and
#' scores between pairwise clusters, respectively.
#' @examples clu_pairs <- FindPairs(object = mouse_kidney_203_Seurat,species = 'Mouse')
#' @import Seurat progress
#' @importFrom crayon red cyan green
#' @importFrom stats pnorm
#' @export FindPairs
#' @details [1] Set \code{use_LRdb} as \code{'LRdb1'}: For \code{'Human'}
#' \code{CellTalkDB}(3,398 human LR pairs)+
#' \code{CellTalkDB}(2,033 mouse LR pairs to human LR pairs)
#' # According to NCBI orthologs (updated in April 28, 2020))
#' @details [2] Set \code{use_LRdb} as \code{'LRdb2'}: For \code{'Human'}
#' \code{CellTalkDB}(3,398 human LR pairs)+
#' \code{SingleCellSignalR}(3,251 human LR pairs)
#' @details [3] Set \code{use_LRdb} as \code{'LRdb3'}: For \code{'Human'}
#' \code{CellTalkDB}(3,398 human LR pairs)+
#' \code{SingleCellSignalR}(3,251 human LR pairs)+
#' \code{CellTalkDB}(2,033 mouse LR pairs to human LR pairs)
#' @details [4] Set \code{use_LRdb} as \code{'LRdb4'}: For \code{'Mouse'}
#' \code{CellTalkDB}(2,033 mouse LR pairs)+
#' \code{CellTalkDB}(3,398 human LR pairs to mouse LR pairs)
#' @details [5] Set \code{use_LRdb} as \code{'LRdb5'}: For \code{'Mouse'}
#' \code{CellTalkDB}(3,398 human LR pairs to mouse LR pairs)+
#' \code{SingleCellSignalR}(3,251 human LR pairs to mouse LR pairs)
#' @details [6] Set \code{use_LRdb} as \code{'LRdb6'}: For \code{'Mouse'}
#' \code{CellTalkDB}(2,033 mouse LR pairs)+
#' \code{CellTalkDB}(3,398 human LR pairs to mouse LR pairs)+
#' \code{SingleCellSignalR}(3,251 human LR pairs to mouse LR pairs)

FindPairs <- function(object, species = NULL, use_LRdb = "LRdb", revise_gene = T, cell_min_pct = 0.25, 
    p_value = 0.05) {
    clu_pairs <- list()
    # extract normalized data from Seurat object
    ndata <- object[["RNA"]]@data
    clu_pairs[[1]] <- ndata
    names(clu_pairs)[1] <- "ndata"
    # check tsne and umap
    if (!"umap" %in% names(object) & !"tsne" %in% names(object)) {
        stop("Please perform dimensionality reduction using umap or tsne!")
    }
    clu_pairs[[2]] <- "NA"
    names(clu_pairs)[2] <- "reductions"
    clu_pairs[[2]] <- object@reductions
    # extract cluster information of all single cells
    clu_info <- Seurat::Idents(object = object)
    clu_info <- data.frame(cell = names(clu_info), cluster = as.character(clu_info), stringsAsFactors = F)
    clu_info[, 1] <- as.character(clu_info[, 1])
    clu_info[, 2] <- as.character(clu_info[, 2])
    clu_num <- unique(clu_info[, 2])
    clu_num1 <- clu_num[nchar(clu_num) == 1]
    clu_num2 <- clu_num[nchar(clu_num) == 2]
    clu_num1 <- clu_num1[order(clu_num1)]
    clu_num2 <- clu_num2[order(clu_num2)]
    clu_num <- c(clu_num1, clu_num2)
    rm(object)
    cat("Raw data matrix includes", ncol(ndata), "cells and", nrow(ndata), "genes", "\n")
    # check number of clusters
    if (length(clu_num) == 1) {
        stop("There is only one cluster, please do clustering analysis or select more clusters for Seurat object!")
    }
    clu_pairs[[3]] <- clu_info
    names(clu_pairs)[3] <- "clu_info"
    Sys.sleep(2)
    # check species
    if (is.null(species)) {
        stop("Please define the species! 'Human' or 'Mouse'.")
    }
    if (!species %in% c("Human", "Mouse")) {
        stop("Please input the right species! 'Human' or 'Mouse'.")
    }
    geneinfo <- geneinfo[geneinfo$species == species, ]
    # check LRdb
    if (!is.data.frame(use_LRdb)) {
        if (!use_LRdb %in% c("LRdb", "LRdb1", "LRdb2", "LRdb3", "LRdb4", "LRdb5", "LRdb6")) {
            stop("Please input the right LRdb! 'LRdb','LRdb1','LRdb2','LRdb3','LRdb4','LRdb5','LRdb6'")
        }
        if (use_LRdb == "LRdb") {
            LR_database <- LR_database[LR_database$species == species & LR_database$LRdb == 
                "CellTalkDB", ]
        }
        if (use_LRdb %in% c("LRdb1", "LRdb2", "LRdb3")) {
            if (species != "Human") {
                stop("Please set species as 'Human'")
            }
            if (use_LRdb == "LRdb1") {
                LR_database <- LR_database[LR_database$LRdb == "CellTalkDB" & LR_database$species %in% 
                  c("Human", "Mouse2Human"), ]
            }
            if (use_LRdb == "LRdb2") {
                LR_database <- LR_database[LR_database$LRdb %in% c("CellTalkDB", "SingleCellSignalR") & 
                  LR_database$species == "Human", ]
            }
            if (use_LRdb == "LRdb3") {
                LR_database <- LR_database[LR_database$LRdb %in% c("CellTalkDB", "SingleCellSignalR") & 
                  LR_database$species %in% c("Human", "Mouse2Human"), ]
            }
        }
        if (use_LRdb %in% c("LRdb4", "LRdb5", "LRdb6")) {
            if (species != "Mouse") {
                stop("Please set species as 'Mouse'")
            }
            if (use_LRdb == "LRdb4") {
                LR_database <- LR_database[LR_database$LRdb == "CellTalkDB" & LR_database$species %in% 
                  c("Mouse", "Human2Mouse"), ]
            }
            if (use_LRdb == "LRdb5") {
                LR_database <- LR_database[LR_database$LRdb %in% c("CellTalkDB", "SingleCellSignalR") & 
                  LR_database$species == "Human2Mouse", ]
            }
            if (use_LRdb == "LRdb6") {
                LR_database <- LR_database[LR_database$LRdb %in% c("CellTalkDB", "SingleCellSignalR") & 
                  LR_database$species %in% c("Mouse", "Human2Mouse"), ]
            }
        }
        LR_database <- LR_database[, -c(4, 5)]
    }
    if (is.data.frame(use_LRdb)) {
        cat(red("Please ensure the custom LR database is a dataframe containing 2 columns: first is 'ligand', second is 'receptor'", 
            "\n"))
        if (ncol(use_LRdb) != 2) {
            stop("Please ensure the custom LR database is a dataframe containing 2 columns: first is 'ligand', second is 'receptor'")
        }
        Sys.sleep(2)
        LR_database <- use_LRdb
        colnames(LR_database) <- c("ligand_gene_symbol", "receptor_gene_symbol")
        LR_database$lr_pair <- paste(LR_database$ligand_gene_symbol, LR_database$receptor_gene_symbol, 
            sep = "_")
        LR_database <- LR_database[, c(3, 1, 2)]
        use_LRdb <- "custom LRdb"
    }
    # revise gene symbol
    if (revise_gene == T) {
        cat(cyan("Revising gene symbols", "\n"))
        Sys.sleep(1)
        genename <- rownames(ndata)
        genename1 <- genename[genename %in% geneinfo$Symbol]
        genename2 <- genename[!genename %in% geneinfo$Symbol]
        genename3 <- genename2[genename2 %in% geneinfo$Synonyms]
        genename4 <- rep("NA", length(genename3))
        for (i in 1:length(genename3)) {
            d1 <- geneinfo[geneinfo$Synonyms == genename3[i], ]$Symbol
            if (length(d1) == 1) {
                genename4[i] <- d1
            }
        }
        genename3 <- c(genename1, genename3)
        genename4 <- c(genename1, genename4)
        genedata <- data.frame(raw_name = genename3, new_name = genename4, stringsAsFactors = F)
        genedata <- genedata[!genedata$new_name == "NA", ]
        genedata1 <- as.data.frame(table(genedata$new_name), stringsAsFactors = F)
        genedata1 <- genedata1[genedata1$Freq == 1, ]
        genedata <- genedata[genedata$new_name %in% genedata1$Var1, ]
        ndata <- ndata[genedata$raw_name, ]
        rownames(ndata) <- genedata$new_name
        cat(green("***Done***", "\n"))
        Sys.sleep(2)
        if (nrow(ndata) == 0) {
            stop("New data contains zero row. Please verify the species information and right symbols of input genes!")
        }
        cat("New data matrix includes", ncol(ndata), "cells,", nrow(ndata), "genes and", length(clu_num), 
            "clusters", "\n")
    }
    clu_pairs[[1]] <- ndata
    # add a function to calculate gene percentage in scdata
    gene_pct <- function(scdata) {
        scdata1 <- scdata[scdata > 0]
        return(length(scdata1)/length(scdata))
    }
    clu_sort <- function(clu_data) {
        clu_data <- clu_data[order(clu_data)]
        clu_data <- paste(clu_data[1], clu_data[2], sep = "_")
        return(clu_data)
    }
    # generating pair-wise clusters
    clu_pair <- NULL
    for (i in 1:length(clu_num)) {
        d1 <- data.frame(cluster1 = rep(clu_num[i], length(clu_num)), cluster2 = clu_num, 
            stringsAsFactors = F)
        clu_pair <- rbind(clu_pair, d1)
    }
    Sys.sleep(2)
    cat(cyan(paste0("Finding significant LR pairs using ", use_LRdb, ": ", nrow(LR_database), 
        " ", species, " LR pairs"), "\n"))
    pb <- progress_bar$new(format = "[:bar] Finished::percent Remaining::eta", total = nrow(clu_pair), 
        clear = FALSE, width = 60, complete = "+", incomplete = "-")
    # generating result file
    res_pairs <- NULL
    # Find ligand-receptor pairs for clusters
    for (i in 1:nrow(clu_pair)) {
        clu_pair1 <- clu_pair[i, ]
        ndata1 <- ndata[, clu_info[clu_info$cluster == clu_pair1$cluster1, ]$cell]
        ndata2 <- ndata[, clu_info[clu_info$cluster == clu_pair1$cluster2, ]$cell]
        res_pairs1 <- LR_database
        res_pairs1$ligand_clu <- clu_pair1$cluster1
        res_pairs1$receptor_clu <- clu_pair1$cluster2
        res_pairs1$ligand_pct <- 0
        res_pairs1$ligand_exp_avg <- 0
        res_pairs1$receptor_pct <- 0
        res_pairs1$receptor_exp_avg <- 0
        res_pairs1$lr_score <- 0
        # match ligand
        ligand_pct <- apply(ndata1, 1, gene_pct)
        ligand_mean <- apply(ndata1, 1, mean)
        ligand_info <- data.frame(gene = names(ligand_pct), pct = as.numeric(ligand_pct), 
            exp_avg = as.numeric(ligand_mean), stringsAsFactors = F)
        ligand_info$z_score <- as.numeric(scale(ligand_info$exp_avg))
        ligand_info$p_value <- pnorm(q = ligand_info$exp_avg, mean = mean(ligand_info$exp_avg), 
            sd = sd(ligand_info$exp_avg), lower.tail = F)
        ligand_info <- ligand_info[ligand_info$p_value < p_value, ]
        if (nrow(ligand_info) > 0) {
            colnames(ligand_info)[1] <- "ligand_gene_symbol"
            res_pairs1 <- base::merge(res_pairs1, ligand_info)
            res_pairs1$ligand_pct <- res_pairs1$pct
            res_pairs1$ligand_exp_avg <- res_pairs1$exp_avg
            res_pairs1 <- res_pairs1[, -c(11:14)]
            res_pairs1 <- res_pairs1[, c(2, 1, 3:10)]
        }
        # match receptor
        receptor_pct <- apply(ndata2, 1, gene_pct)
        receptor_mean <- apply(ndata2, 1, mean)
        receptor_info <- data.frame(gene = names(receptor_pct), pct = as.numeric(receptor_pct), 
            exp_avg = as.numeric(receptor_mean), stringsAsFactors = F)
        receptor_info$z_score <- as.numeric(scale(receptor_info$exp_avg))
        receptor_info$p_value <- pnorm(q = receptor_info$exp_avg, mean = mean(receptor_info$exp_avg), 
            sd = sd(receptor_info$exp_avg), lower.tail = F)
        receptor_info <- receptor_info[receptor_info$p_value < p_value, ]
        if (nrow(receptor_info) > 0) {
            colnames(receptor_info)[1] <- "receptor_gene_symbol"
            res_pairs1 <- base::merge(res_pairs1, receptor_info)
            res_pairs1$receptor_pct <- res_pairs1$pct
            res_pairs1$receptor_exp_avg <- res_pairs1$exp_avg
            res_pairs1 <- res_pairs1[, -c(11:14)]
            res_pairs1 <- res_pairs1[, c(2, 3, 1, 4:10)]
        }
        res_pairs1$lr_score <- res_pairs1$ligand_exp_avg * res_pairs1$receptor_exp_avg
        res_pairs1 <- res_pairs1[res_pairs1$ligand_pct > cell_min_pct & res_pairs1$receptor_pct > 
            cell_min_pct, ]
        res_pairs[[i]] <- res_pairs1
        names(res_pairs)[i] <- paste0("cluster", clu_pair1$cluster1, "-", "cluster", clu_pair1$cluster2)
        pb$tick()
    }
    Sys.sleep(2)
    cat(green("***Done***", "\n"))
    Sys.sleep(2)
    # formatting results
    res_pairs1 <- res_pairs[[1]]
    for (i in 2:length(res_pairs)) {
        res_pairs1 <- rbind(res_pairs1, res_pairs[[i]])
    }
    res_pairs1$clu_com <- apply(res_pairs1[, c("ligand_clu", "receptor_clu")], 1, clu_sort)
    res_pairs1$clu1_clu2 <- paste(res_pairs1$ligand_clu, res_pairs1$receptor_clu, sep = "_")
    clu_pairs[[4]] <- res_pairs1
    names(clu_pairs)[4] <- "clu_pairs"
    cat(paste0("There are ", nrow(res_pairs1), " highly expressed LR pairs for ", length(unique(res_pairs1$clu1_clu2)), 
        " pairwise clusters"))
    clu_crosstalk <- "NA"
    clu_crosstalk_undir <- "NA"
    if (nrow(res_pairs1) > 0) {
        # formatting directed cell-cell communication table
        clu_crosstalk <- as.data.frame(matrix(0, nrow = nrow(clu_pair), ncol = 6))
        colnames(clu_crosstalk) <- c("cluster1_source", "cluster2_target", "clu_com", "clu1_clu2", 
            "LR_number", "LR_score")
        clu_crosstalk$cluster1_source <- clu_pair$cluster1
        clu_crosstalk$cluster2_target <- clu_pair$cluster2
        clu_crosstalk$clu_com <- apply(clu_crosstalk[, c(1, 2)], 1, clu_sort)
        clu_crosstalk$clu1_clu2 <- paste(clu_crosstalk$cluster1_source, clu_crosstalk$cluster2_target, 
            sep = "_")
        for (i in 1:nrow(clu_crosstalk)) {
            clu_clu_name <- clu_crosstalk$clu1_clu2[i]
            if (clu_clu_name %in% res_pairs1$clu1_clu2) {
                clu_clu <- res_pairs1[res_pairs1$clu1_clu2 == clu_clu_name, ]
                clu_crosstalk$LR_number[i] <- nrow(clu_clu)
                clu_crosstalk$LR_score[i] <- sum(clu_clu$lr_score)
            }
        }
        clu_crosstalk_undir <- data.frame(clu_com = clu_crosstalk$clu_com, clu_com = clu_crosstalk$clu_com, 
            stringsAsFactors = F)
        clu_crosstalk_undir <- unique(clu_crosstalk_undir)
        clu_crosstalk_undir <- clu_crosstalk[rownames(clu_crosstalk_undir), ]
        rownames(clu_crosstalk_undir) <- 1:nrow(clu_crosstalk_undir)
        colnames(clu_crosstalk_undir)[1:2] <- c("cluster1", "cluster2")
        clu_com <- clu_crosstalk_undir$clu_com
        for (i in 1:length(clu_com)) {
            clu_crosstalk1 <- clu_crosstalk[clu_crosstalk$clu_com == clu_com[i], ]
            clu_crosstalk_undir$LR_number[i] <- sum(clu_crosstalk1$LR_number)
            clu_crosstalk_undir$LR_score[i] <- sum(clu_crosstalk1$LR_score)
        }
    }
    clu_pairs[[5]] <- clu_crosstalk
    names(clu_pairs)[5] <- "clu_crosstalk"
    clu_pairs[[6]] <- clu_crosstalk_undir
    names(clu_pairs)[6] <- "clu_crosstalk_undir"
    para <- list(species = species, use_LRdb = use_LRdb, cell_min_pct = cell_min_pct, p_value = p_value)
    clu_pairs[[7]] <- para
    names(clu_pairs)[7] <- "para"
    Sys.sleep(2)
    return(clu_pairs)
}
