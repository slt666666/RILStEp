% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RILStEp.R
\name{extract_peak_qtls}
\alias{extract_peak_qtls}
\title{extracts QTL-like SNPs based on GWASpoly package}
\usage{
extract_peak_qtls(
  data_for_gwas,
  phenotype_name,
  qtl_threshold_model,
  qtl_threshold_value,
  output,
  core_num
)
}
\arguments{
\item{data_for_gwas}{phenotype and genotype data for GWASpoly from load_data}

\item{phenotype_name}{trait's name}

\item{qtl_threshold_model}{threshold method to select QTLs, "p_value" or "Bonferroni" or "FDR" or "permute"}

\item{qtl_threshold_value}{threshold value to select QTLs depends on qtl_threshold_model (p-value or q-value)}

\item{output}{prefix of output files}

\item{core_num}{the number of CPU cores to use}
}
\description{
extracts QTL-like SNPs based on GWASpoly package
}
