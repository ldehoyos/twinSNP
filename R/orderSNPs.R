#' Order a vector of SNPs
#'
#'
#' @param SNPvector A character vector with the SNPs in format chr:pos.
#'
#' @return Vector with all SNPs ordered by chromosome and then by position.
#'
#' @examples
#' myVector_ordered <- orderSNPs(myVector)
#'
#' @export
#'
orderSNPs <- function(SNPvector) {
  div  <- data.table(chr = as.numeric(gsub(":.*", "", SNPvector)),
                     pos = as.numeric(gsub(".:", "", SNPvector)))
  div_ordered <- div[order(chr, pos), ]
  SNPvector_out <- paste(div_ordered$chr, div_ordered$pos, sep = ":")
  return(SNPvector_out)
}
