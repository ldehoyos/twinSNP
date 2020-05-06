#' List the matches of each SNP.
#'
#' This function reads the output from "SNPmatch" and outputs the same
#' information as a list, where each object is a new input SNP which has
#' inside a vector with all its matches. It is the same as "listMatches" but
#' without ordering the variants (this saves computation time).
#'
#' @param matched_snps Data frame with columns "input_snp" and "matched_snp".
#'
#' @return List with all the SNPs and a vector of its matches.
#'
#' @examples
#' SNP <- listMatches2(matched_snps)
#'
#' @export
#'
listMatches2 <- function(matched_snps){
  ### Create a data frame with only the columns of the input and matched SNPs.
  matched_snps <- data.frame(input_snp = matched_snps$input_snp,
                             matched_snp = matched_snps$snpID,
                             stringsAsFactors = F)

  ### See how many matches are there per SNP
  snp_count <- table(matched_snps$input_snp)
  snp_count <- data.frame(snp_count); colnames(snp_count) <- c("input_snp", "matches")
  snp_count$input_snp <- as.character(snp_count$input_snp)

  SNPinfo <- list()
  for (i in 1:nrow(snp_count)){SNPinfo[[i]] <- findSNPs2(matched_snps, snp_count[i, 1])} #
  names(SNPinfo) <- snp_count$input_snp # Assign the input_snps to their vectors of matched snps.

  return(SNPinfo)
}
