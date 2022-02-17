#' extracts QTL-like SNPs based on GWASpoly package
#' @importFrom GWASpoly GWASpoly set.threshold manhattan.plot get.QTL set.K
#' @importFrom grDevices dev.off pdf
#' @importFrom utils write.csv
#' @param data_for_gwas phenotype and genotype data for GWASpoly from load_data
#' @param phenotype_name trait's name
#' @param qtl_threshold_model threshold method to select QTLs, "p_value" or "Bonferroni" or "FDR" or "permute"
#' @param qtl_threshold_value threshold value to select QTLs depends on qtl_threshold_model (p-value or q-value)
#' @param output prefix of output files
#' @param core_num the number of CPU cores to use
extract_peak_qtls <-
  function(data_for_gwas,
           phenotype_name,
           qtl_threshold_model,
           qtl_threshold_value,
           output,
           core_num) {

    data.loco <- set.K(data_for_gwas, LOCO=TRUE, n.core=2)
    
    params <- set.params()

    gwaspoly_result <- GWASpoly(
      data = data.loco,
      models = c("additive"),
      traits = c(phenotype_name),
      params = params,
      n.core = core_num
    )

    if (qtl_threshold_model == "p_value") {
      scores <- 10 ** (-gwaspoly_result@scores[[phenotype_name]])
      qtl_ind <- (!is.na(scores)) & (scores < qtl_threshold_value)
      if (sum(qtl_ind) > 0) {
        Trait <- rep(phenotype_name, sum(qtl_ind))
        Model <- rep(qtl_threshold_model, sum(qtl_ind))
        Threshold <- rep(qtl_threshold_value, sum(qtl_ind))
        Scores <- gwaspoly_result@scores[[phenotype_name]][qtl_ind, ]
        Map <- gwaspoly_result@map[qtl_ind, ]
        qtls <- cbind(Trait, Model, Threshold, Map, Scores)
      } else {
        qtls <- data.frame()
      }
      gwaspoly_threshold <- set.threshold(
        gwaspoly_result,
        method = "FDR",
        level = 0.05
      )
    } else {
      gwaspoly_threshold <- set.threshold(
        gwaspoly_result,
        method = qtl_threshold_model,
        level = qtl_threshold_value
      )
      qtls <- get.QTL(
        gwaspoly_threshold,
        traits = phenotype_name,
      )
    }
    pdf(paste(output, "_QTL_manhattan_plot.pdf", sep = ""))
    manhattan.plot(gwaspoly_threshold, trait = phenotype_name, model = "additive")
    dev.off()

    write.csv(qtls, file = paste(output, "_QTL_score.csv", sep = ""))

    ### Extract SNP positons which have peak scores
    peak_qtls <- c()
    if (dim(qtls)[1] > 0) {

      before_marker <- qtls$Marker[1]
      before_chr <- qtls$Chrom[1]
      base_marker <- qtls$Marker[1]
      base_score <- 0

      for (i in 1:dim(qtls)[1]) {

        marker <- qtls[i, 4]
        chr <- qtls[i, 5]
        score <- qtls[i, 9]

        other_peak <- (marker - before_marker) > 500
        other_chr <- chr != before_chr
        last_qtl <- i == dim(qtls)[1]
        if (other_peak || other_chr || last_qtl) {
          peak_qtls <- c(peak_qtls, base_marker)
          base_score <- score
          base_marker <- marker
        }else{
          if (score > base_score) {
            base_score <- score
            base_marker <- marker
          }
        }
        before_marker <- marker
        before_chr <- chr
      }
    }
    print("finish GWAS process")
    print(paste("QTL marker:", peak_qtls))
    return(peak_qtls)
  }

#' Convert position data (chrXX_ZZZZZ) to Marker number
#' @param position position data of SNP
#' @param marker_data indices of all SNP positions
pos2marker <- function(position, marker_data) {
  chrom <- strsplit(position, "_")[[1]][1]
  position <- as.numeric(strsplit(position, "_")[[1]][2])
  index <- (marker_data$Chrom == chrom & marker_data$Position == position)
  return(marker_data[index, "Marker"])
}

#' Convert Marker number to position data (chrXX_ZZZZZ)
#' @param marker marker index of SNP
#' @param marker_data indices of all SNP positions
marker2pos <- function(marker, marker_data) {
  chrom <- marker_data[marker_data$Marker == marker, "Chrom"]
  position <- marker_data[marker_data$Marker == marker, "Position"]
  return(paste(chrom, position, sep = "_"))
}

#' Calculate BayesFactor for Comparison between Model1 (without Epistasis) vs Model2 (with Epistasis)
#' @importFrom BayesFactor lmBF
#' @importFrom stats as.formula na.omit xtabs
#' @param pair_index indices of SNP pair to check epistasis
#' @param genotype_phenotype_dataset dataset of genotype and phenotype
#' @param phenotype_name trait's name
#' @param peak_qtls QTL-like SNP index list from extract_peak_qtls() or user specify
#' @param pairs SNP pairs
check_epistasis <-
  function(pair_index,
           genotype_phenotype_dataset,
           phenotype_name,
           peak_qtls,
           pairs) {

    snp1 <- pairs[1, pair_index]
    snp2 <- pairs[2, pair_index]
    qtl_num <- length(peak_qtls)
    ### SNPs index + 2 because of row.names & phenotype columns
    use_columns <- as.character(c(phenotype_name, peak_qtls, snp1, snp2))
    extract_data <- genotype_phenotype_dataset[, use_columns]
    extract_data <- na.omit(extract_data)
    remove_row <- rowSums(extract_data[, 2:(3 + qtl_num)] == 1) > 0
    extract_data <- extract_data[!remove_row, ]
    extract_data[, 2:(3 + qtl_num)] <- extract_data[, 2:(3 + qtl_num)] / 2

    h_x <- extract_data[, qtl_num + 2] == 0
    f_x <- extract_data[, qtl_num + 2] == 1
    x_h <- extract_data[, qtl_num + 3] == 0
    x_f <- extract_data[, qtl_num + 3] == 1

    base <- rep(0, dim(extract_data)[1])
    base[f_x & x_f] <- 1
    extract_data <- transform(extract_data, Epi1 = base)
    base <- rep(0, dim(extract_data)[1])
    base[f_x & x_h] <- 1
    extract_data <- transform(extract_data, Epi2 = base)
    base <- rep(0, dim(extract_data)[1])
    base[h_x & x_f] <- 1
    extract_data <- transform(extract_data, Epi3 = base)

    # only QTL model
    exp_val <- paste(colnames(extract_data)[2:(3 + qtl_num)], collapse = " + ")
    model1_fomula <- as.formula(paste(phenotype_name, exp_val, sep = " ~ "))
    # with Epistasis model
    exp_val <- paste(c(colnames(extract_data)[2:(3 + qtl_num)], c("Epi1", "Epi2", "Epi3")), collapse = " + ")
    model2_fomula <- as.formula(paste(phenotype_name, exp_val, sep = " ~ "))
    mod1 <- lmBF(model1_fomula, data = extract_data)
    mod2 <- lmBF(model2_fomula, data = extract_data)
    result <- mod2 / mod1
    return(list(snp1, snp2, exp(result@bayesFactor$bf)))
  }

#' Calculate BayesFactor for Comparison between Model1 (without Epistasis) vs Model2 (with Epistasis)
#' @importFrom BayesFactor lmBF
#' @importFrom stats as.formula na.omit xtabs
#' @param pair_index indices of SNP pair to check epistasis
#' @param genotype_phenotype_dataset dataset of genotype and phenotype
#' @param phenotype_name trait's name
#' @param peak_qtls QTL-like SNP index list from extract_peak_qtls() or user specify
#' @param pairs SNP pairs
check_epistasis2 <-
  function(pair_index,
           genotype_phenotype_dataset,
           phenotype_name,
           peak_qtls,
           pairs) {

    snp1 <- pairs[1, pair_index]
    snp2 <- pairs[2, pair_index]
    qtl_num <- length(peak_qtls)
    ### SNPs index + 2 because of row.names & phenotype columns
    use_columns <- as.character(c(phenotype_name, peak_qtls, snp1, snp2))
    extract_data <- genotype_phenotype_dataset[, use_columns]
    extract_data <- na.omit(extract_data)

    h_x <- extract_data[, qtl_num + 2] == 0
    m_x <- extract_data[, qtl_num + 2] == 1
    f_x <- extract_data[, qtl_num + 2] == 2
    x_h <- extract_data[, qtl_num + 3] == 0
    x_m <- extract_data[, qtl_num + 3] == 1
    x_f <- extract_data[, qtl_num + 3] == 2

    base <- rep(0, dim(extract_data)[1])
    base[f_x & x_f] <- 1
    extract_data <- transform(extract_data, Epi1 = base)
    base <- rep(0, dim(extract_data)[1])
    base[f_x & x_m] <- 1
    extract_data <- transform(extract_data, Epi2 = base)
    base <- rep(0, dim(extract_data)[1])
    base[f_x & x_h] <- 1
    extract_data <- transform(extract_data, Epi3 = base)
    base <- rep(0, dim(extract_data)[1])
    base[m_x & x_f] <- 1
    extract_data <- transform(extract_data, Epi4 = base)
    base <- rep(0, dim(extract_data)[1])
    base[m_x & x_m] <- 1
    extract_data <- transform(extract_data, Epi5 = base)
    base <- rep(0, dim(extract_data)[1])
    base[m_x & x_h] <- 1
    extract_data <- transform(extract_data, Epi6 = base)
    base <- rep(0, dim(extract_data)[1])
    base[h_x & x_f] <- 1
    extract_data <- transform(extract_data, Epi7 = base)
    base <- rep(0, dim(extract_data)[1])
    base[h_x & x_m] <- 1
    extract_data <- transform(extract_data, Epi8 = base)

    # only QTL model
    exp_val <- paste(colnames(extract_data)[2:(3 + qtl_num)], collapse = " + ")
    model1_fomula <- as.formula(paste(phenotype_name, exp_val, sep = " ~ "))
    # with Epistasis model
    exp_val <- paste(c(colnames(extract_data)[2:(3 + qtl_num)], c("Epi1", "Epi2", "Epi3", "Epi4", "Epi5", "Epi6", "Epi7", "Epi8")), collapse = " + ")
    model2_fomula <- as.formula(paste(phenotype_name, exp_val, sep = " ~ "))
    mod1 <- lmBF(model1_fomula, data = extract_data)
    mod2 <- lmBF(model2_fomula, data = extract_data)
    result <- mod2 / mod1
    return(list(snp1, snp2, exp(result@bayesFactor$bf)))
  }

#' draw heatmap of bayes factor
#' @importFrom heatmap3 heatmap3
#' @importFrom grDevices dev.off pdf
#' @param result bayesfactor of all combination of SNPs
#' @param filename name of save file
#' @param marker_data indices of all SNP positions
draw_heatmap <- function(result, filename, marker_data) {
  result$score <- log10(result$score)
  result_for_heatmap <- xtabs(score ~ first + second, result)
  rownames(result_for_heatmap) <-
    sapply(rownames(result_for_heatmap), marker2pos, marker_data = marker_data)
  colnames(result_for_heatmap) <-
    sapply(colnames(result_for_heatmap), marker2pos, marker_data = marker_data)

  if (dim(result_for_heatmap)[1] > 24) {
    row_number <- dim(result_for_heatmap)[1]
    lab_row <- rep("", row_number)
    lab_row[seq(1, row_number, row_number %/% 24)] <-
      rownames(result_for_heatmap)[seq(1, row_number, row_number %/% 24)]
  } else {
    lab_row <- rownames(result_for_heatmap)
  }
  if (dim(result_for_heatmap)[2] > 24) {
    col_number <- dim(result_for_heatmap)[2]
    lab_col <- rep("", col_number)
    lab_col[seq(1, col_number, col_number %/% 24)] <-
      colnames(result_for_heatmap)[seq(1, col_number, col_number %/% 24)]
  } else {
    lab_col <- colnames(result_for_heatmap)
  }
  pdf(filename)
  heatmap3(result_for_heatmap,
           Rowv = NA,
           Colv = NA,
           scale = "none",
           labRow = lab_row,
           labCol = lab_col,
           margins = c(10, 10)
  )
  dev.off()
}

#' load dataset to use GWAS and RILStEp
#' @export
#' @importFrom utils read.csv
#' @importFrom GWASpoly read.GWASpoly
#' @param phenotype_path path to csv file of phenotype data
#' @param genotype_path path to csv file of genotype data
#' @param phenotype_name trait's name
load_data <- function (phenotype_path, genotype_path, phenotype_name) {
  ### loading phenotype & genotype data
  phenotype_data <-
    read.csv(phenotype_path, sep = ",", header = TRUE, row.names = 1)
  genotype_data <-
    read.csv(genotype_path, sep = ",", header = TRUE)
  gwas_data <- read.GWASpoly(
    ploidy = 2,
    pheno.file = phenotype_path,
    geno.file = genotype_path,
    format = "numeric",
    n.traits = 1,
    delim = ","
  )
  loaded_data <- list(phenotype_name, phenotype_data, genotype_data, gwas_data)
  names(loaded_data) <- c("phenotype_name", "phenotype", "genotype", "for_gwas")
  return(loaded_data)
}

#' detect epistasis candidates based on bayes factor
#' @export
#' @importFrom future availableCores plan
#' @importFrom data.table transpose
#' @importFrom furrr future_map
#' @importFrom future multiprocess
#' @importFrom utils combn write.csv
#' @param loaded_data phenotype and genotype dataset from load_data
#' @param output prefix of output files
#' @param qtl_threshold_model threshold method to select QTLs, "p_value" or "Bonferroni" or "FDR" or "permute"
#' @param qtl_threshold_value threshold value to select QTLs depends on qtl_threshold_model (p-value or q-value)
#' @param interval select 1 SNP in this interval number
#' @param region specify SNP regions for detecting epistasis
#' @param qtls specify QTL-like SNPs
#' @param core_num specify the number of CPU core to use
#' @param heterozygous consider heterozygous(TRUE) or not(FALSE)
rilstep <-
  function(loaded_data,
           output,
           qtl_threshold_model = "FDR",
           qtl_threshold_value = 0.05,
           interval = 1,
           region = NA,
           qtls = NA,
           core_num = NA,
           heterozygous = FALSE) {

    print("start RILStEp !!")

    ### set core_num
    if (is.na(core_num)) {
      core_num <- availableCores() - 1
    }

    ### loading phenotype & genotype data
    phenotype_data <- loaded_data$phenotype
    genotype_data <- loaded_data$genotype
    marker_data <- genotype_data[, 1:3]
    phenotype_name <- loaded_data$phenotype_name

    ### extract QTL candidate SNPs
    if (is.na(qtls[1])) {
      peak_qtls <- extract_peak_qtls(loaded_data$for_gwas, phenotype_name, qtl_threshold_model, qtl_threshold_value, output, core_num)
    } else if (qtls == FALSE){
      peak_qtls <- c()
    } else {
      peak_qtls <- sapply(qtls, pos2marker, marker_data = marker_data)
    }

    ### transpose genotype data
    tmp_genotype_data <-
      data.table::transpose(genotype_data[, 4:dim(genotype_data)[2]])
    rownames(tmp_genotype_data) <-
      colnames(genotype_data[, 4:dim(genotype_data)[2]])
    colnames(tmp_genotype_data) <- 1:dim(genotype_data)[1]

    ### extract QTL & check SNPs data & make check pairs
    no_region <- is.na(region[1])

    if (no_region) {
      check_snps <- seq(1, dim(tmp_genotype_data)[2], by = interval)
      pairs <- combn(x = check_snps, m = 2)
    }else{
      region1 <- unlist(strsplit(region[1], ":"))
      region1_markers <- sapply(region1, pos2marker, marker_data = marker_data)
      if (length(region1_markers) == 1) {
        region1_snps <- region1_markers
      } else if (length(region1_markers) == 2) {
        region1_snps <- seq(region1_markers[1], region1_markers[2], by = interval)
      }

      if (length(region) == 1) {
        region2_markers <- c(1, dim(marker_data)[1])
      }else{
        region2 <- unlist(strsplit(region[2], ":"))
        region2_markers <- sapply(region2, pos2marker, marker_data = marker_data)
      }
      if (length(region2_markers) == 1) {
        region2_snps <- region2_markers
      } else if (length(region2_markers) == 2) {
        region2_snps <- seq(region2_markers[1], region2_markers[2], by = interval)
      }
      check_snps <- c(region1_snps, region2_snps)
      pairs <- t(expand.grid(region1_snps, region2_snps))
      pairs <- pairs[, pairs[1, ] != pairs[2, ]]
    }
    tmp_genotype_data <- tmp_genotype_data[, c(peak_qtls, check_snps)]
    genotype_phenotype_dataset <-
      merge(phenotype_data, tmp_genotype_data, by = "row.names")

    ### calculate bayes factor
    if (heterozygous) {
      plan(multiprocess, workers = core_num)
      t <- proc.time()
      result <- future_map(
        1:dim(pairs)[2],
        check_epistasis2,
        genotype_phenotype_dataset = genotype_phenotype_dataset,
        phenotype_name = phenotype_name,
        peak_qtls = peak_qtls,
        pairs = pairs
      )
    } else {
      plan(multiprocess, workers = core_num)
      t <- proc.time()
      result <- future_map(
        1:dim(pairs)[2],
        check_epistasis,
        genotype_phenotype_dataset = genotype_phenotype_dataset,
        phenotype_name = phenotype_name,
        peak_qtls = peak_qtls,
        pairs = pairs
      )
    }
    result <- data.frame(matrix(unlist(result), nrow = length(result), byrow = T))
    colnames(result) <- c("first", "second", "score")

    # draw heatmap
    heatmap_name <- paste(output, "_heatmap.pdf", sep = "")
    if (no_region) {
      draw_heatmap(result, heatmap_name, marker_data)
    } else {
      if ((length(region1_markers) > 1) && (length(region2_markers) > 1)) {
        draw_heatmap(result, heatmap_name, marker_data)
      }
    }

    # add marker position
    write.csv(marker_data, file = paste(output, "_marker_data.csv", sep = ""))
    write.csv(result, file = paste(output, "_epistasis_BFscore.csv", sep = ""))

    print(proc.time() - t)
    return(result[order(-result$score), ])
  }
