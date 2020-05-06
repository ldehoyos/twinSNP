#' Randomly choose SNPs with similar characteristics to your input SNPs.
#'
#' Based on the information provided of the matches, the function chooses randomly
#' the matches with similar characteristics to the input SNPs. In case the number of original
#' SNPs (numSNPs) is not the same as the length of SNPlist, the difference will be filled with
#' random variants of the already sampled available SNPs.
#'
#' @param SNPlist A list where each object is a input SNP (specified in the name of the object) with a vector of its matches.
#' @param numSNPs Numeric value with the information of the original number of SNPs.
#'
#' @return Vector with a randomly chosen number of variants with similar characteristics to the input SNPs.
#'
#' @examples
#' permVector <- permSNPs(SNP_list, 8000)
#'
#' @export
permSNPs <-  function(SNPlist, numSNPs){
  # numSNPs is the number of original input SNPs.
  # Set a vector to save the input variants and the matched ouput variants
  permSNPs_in  <-   data.frame(matrix(ncol= 2, nrow= numSNPs))
  permSNPs_out <- as.character(matrix(ncol= 1, nrow= numSNPs))

  # Save the input variants with a match and order them by number of matches.
  permSNPs_in[1:length(SNPlist),1] <- names(SNPlist)

  nMatches <- lengths(SNPlist)
  permSNPs_in[,2] <- as.numeric(nMatches[match(permSNPs_in[,1], names(nMatches))])
  permSNPs_in <- permSNPs_in[order(permSNPs_in[,2]),]

  # In case the number of variants is not the same because they have been variants
  #   that do not have a match, fill those places with random matches of the rest of the variants.
  # To do so, only those variants that have more matches than the threshold (thr) will be used
  #   for this purpose.
  #

  thr <- 10

  if(numSNPs-length(SNPlist) != 0){
    cat("Attention: Out of ", numSNPs, " SNPs, ", numSNPs-length(SNPlist) , " do not have matches.\n", sep = "")
    if (length(permSNPs_in[which(permSNPs_in[,2]<thr),1]) != 0){
      permSNPs_in[(length(SNPlist)+1):numSNPs,1] <-
        sample(x= names(SNPlist)[-match(permSNPs_in[which(permSNPs_in[,2]<thr),1], names(SNPlist))],
               size = numSNPs-length(SNPlist))
    }else{
      permSNPs_in[(length(SNPlist)+1):numSNPs,1] <- sample(x= names(SNPlist), size = numSNPs-length(SNPlist), replace = F)
    }
  }

  permSNPs_in[,2] <- as.numeric(nMatches[match(permSNPs_in[,1], names(nMatches))])
  permSNPs_in <- permSNPs_in[order(permSNPs_in[,2]),1]

  matxSNP <-c()
  for (i in 1:length(permSNPs_in)){
    matxSNP <- SNPlist[[grep(pattern = permSNPs_in[i], x = names(SNPlist))]]
    if(length(which(!is.na(match(matxSNP, permSNPs_out)))) != 0){
      matxSNP <- matxSNP[-which(!is.na(match(matxSNP, permSNPs_out)))]
    }
    permSNPs_out[i] <- sample(x= matxSNP, size = 1)
  }

  return(permSNPs_out)
}
