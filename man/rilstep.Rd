% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RILStEp.R
\name{rilstep}
\alias{rilstep}
\title{detect epistasis candidates based on bayes factor}
\usage{
rilstep(
  loaded_data,
  output,
  qtl_threshold_model = "FDR",
  qtl_threshold_value = 0.05,
  interval = 1,
  region = NA,
  qtls = NA,
  core_num = NA,
  heterozygous = FALSE
)
}
\arguments{
\item{loaded_data}{phenotype and genotype dataset from load_data}

\item{output}{prefix of output files}

\item{qtl_threshold_model}{threshold method to select QTLs, "p_value" or "Bonferroni" or "FDR" or "permute"}

\item{qtl_threshold_value}{threshold value to select QTLs depends on qtl_threshold_model (p-value or q-value)}

\item{interval}{select 1 SNP in this interval number}

\item{region}{specify SNP regions for detecting epistasis}

\item{qtls}{specify QTL-like SNPs}

\item{core_num}{specify the number of CPU core to use}

\item{heterozygous}{consider heterozygous(TRUE) or not(FALSE)}
}
\description{
detect epistasis candidates based on bayes factor
}
