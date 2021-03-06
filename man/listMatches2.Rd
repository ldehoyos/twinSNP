% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/listMatches2.R
\name{listMatches2}
\alias{listMatches2}
\title{List the matches of each SNP.}
\usage{
listMatches2(matched_snps)
}
\arguments{
\item{matched_snps}{Data frame with columns "input_snp" and "matched_snp".}
}
\value{
List with all the SNPs and a vector of its matches.
}
\description{
This function reads the output from "SNPmatch" and outputs the same
information as a list, where each object is a new input SNP which has
inside a vector with all its matches. It is the same as "listMatches" but
without ordering the variants (this saves computation time).
}
\examples{
SNP <- listMatches2(matched_snps)

}
