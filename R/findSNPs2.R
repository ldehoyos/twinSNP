#' Find the SNPs that match with the input SNP.
#'
#' The same as "findSNPs" but without ordering the SNPs.
#'
#' @param x Data frame with two columns (input_snp and matched_snps)
#' @param y Character with the input snp (chr:pos)
#'
#' @return Vector with all the matched.
#'
#' @examples
#' SNP <- findSNPs2(matched_snps, "1:1361641")
#'
#' @export
#'
findSNPs2 <- function(x, y) {
  SNPmatch <- grep(pattern = y, x[, 1])
  matchedSNPs <- x[SNPmatch, 2]
  return(matchedSNPs)
}
