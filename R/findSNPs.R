#' Find the SNPs that match with the input SNP.
#'
#'
#' @param x Data frame with two columns (input_snp and matched_snps)
#' @param y Character with the input snp (chr:pos)
#'
#' @return Vector with all the matched SNPs ordered by chromosome and then by position.
#'
#' @examples
#' SNP <- findSNPs(matched_snps, "1:1361641")
#'
#' @export
findSNPs <- function(x, y) {
  SNPmatch <- grep(pattern = y, x[, 1])
  matchedSNPs <- x[SNPmatch, 2]
  matchedSNPs_ordered <- orderSNPs(matchedSNPs)
  return(matchedSNPs_ordered)
}
