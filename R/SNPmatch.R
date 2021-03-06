#' Find SNPs with similar characteristics to the input SNPs.
#'
#' Based on the information provided in the database, the function matches
#' the input SNPs with those SNPs in the database based on the parameters of
#' minor allele frequency, gene density, distance to nearest gene and
#' LD buddies.
#'
#' It also writes the information on the included ("input_snps_annotated.txt")
#' and excluded ("input_snps_excluded.txt") SNPs depending on if they are available
#' on the database.
#'
#' Also, at the end, if there is any SNP that does not have any match,
#' it will also write it in a txt ("input_snps_annotated_unmatched.txt")
#' and exclude it from the final output ("matched_snps_annotation.txt").
#'
#' @param database Data frame with information of SNPs
#' @param inputSNPs Character vector with the input SNPs in format chr:pos.
#' @param MAF Numeric value with interval limits of minor allele frequency.
#' @param GD Numeric value with interval limits of gene density.
#' @param DNG Numeric value with interval limits of distance to nearest gene.
#' @param LDB Numeric value with interval limits of LD buddies.
#'
#' @return Data table with the input SNPs, their matches and information about the matches.
#'
#' @examples
#' myDatabase <- data.frame(snpID= c("10:10001753", "10:10001794"), snp_maf= c(0.074, 0.410), gene_count= c(0,0), dist_nearest_gene_snpsnap= c(98932,98891), friends_ld05= c(35,168))
#' myVars <- data.frame(input_snp= c("18:53093251", "12:57682956", "10:10001794"))
#' matcx <- SNPmatch(database = myDatabase, inputSNPs = myVars, MAF = 0.05, GD = 0.5, DNG = 0.5, LDB = 0.5)
#'
#' @export
SNPmatch <- function(database, inputSNPs, MAF, GD, DNG, LDB){
  # 1. INCLUDED AND EXCLUDED SNPs.
  matxIn <- database[match(inputSNPs$input_snp, database$snpID),]
  excluded_snps <- inputSNPs[which(is.na(matxIn$snpID)),]
  included_snps <-   matxIn[-which(is.na(matxIn$snpID)),]

  # Save the files
  write.table(excluded_snps, "input_snps_excluded.txt", row.names = F, quote = F, col.names = F)
  write.table(included_snps, "input_snps_annotated.txt", row.names = F, quote = F, col.names = T)


  # 2. DO MATCHING
  # Take out the SNPs included to match
  datB <- data.table(database[-as.numeric(rownames(included_snps)),])
  rownames(datB) <- c()

  ### SET THE INTERVALS

  # Set the MAF interval at 5%
  included_snps$maf_min <- included_snps$snp_maf - MAF
  included_snps$maf_max <- included_snps$snp_maf + MAF

  # Set the gene density interval
  included_snps$gdens_min <- included_snps$gene_count - GD*included_snps$gene_count
  included_snps$gdens_max <- included_snps$gene_count + GD*included_snps$gene_count

  # Set the distance to nearest gene interval
  included_snps$dist_interval_min <- included_snps$dist_nearest_gene_snpsnap - DNG*included_snps$dist_nearest_gene_snpsnap
  included_snps$dist_interval_max <- included_snps$dist_nearest_gene_snpsnap + DNG*included_snps$dist_nearest_gene_snpsnap

  # Set the LD interval
  included_snps$ld_interval_min <- included_snps$friends_ld05 - LDB*included_snps$friends_ld05
  included_snps$ld_interval_max <- included_snps$friends_ld05 + LDB*included_snps$friends_ld05

  # Do the analysis and matching of SNPs
  matched_snps <-list()

  for (i in 1:nrow(included_snps)){

    # Do the filtering to select the variants that meet the requirements.
    A1 <- datB[data.table::between(snp_maf, lower= included_snps$maf_min[i], upper= included_snps$maf_max[i])]
    A2 <-   A1[data.table::between(gene_count, lower= included_snps$gdens_min[i], upper=included_snps$gdens_max[i])]
    A3 <-   A2[data.table::between(dist_nearest_gene_snpsnap, lower= included_snps$dist_interval_min[i], upper= included_snps$dist_interval_max[i])]
    subSnp <- A3[data.table::between(friends_ld05, lower= included_snps$ld_interval_min[i], upper=included_snps$ld_interval_max[i])]
    rownames(subSnp) <- c()

    matched_snps[[i]] <- data.table(input_snp= included_snps$snpID[i], subSnp, stringsAsFactors = F)
  }

  matchedSNPs <-  do.call(rbind, matched_snps)

  # 3. NON-MATCHED SNPs
  ### Ensure that all the SNPs have  a match, and if not, exclude those SNPs.

  if(length(which(is.na(matchedSNPs$snpID))) != 0){
    unmatched_snps <- matchedSNPs[which(is.na(matchedSNPs$snpID)),]
    write.table(unmatched_snps, "input_snps_annotated_unmatched.txt", quote = F,
                col.names = T, row.names = F)

    matched_snps2  <- matchedSNPs[-which(is.na(matchedSNPs$snpID)),]
    write.table(matched_snps2, "matched_snps_annotation.txt", quote = F, col.names = T, row.names = F)
    return(matched_snps2)
  }else{
    write.table(matchedSNPs, "matched_snps_annotation.txt", quote = F, col.names = T, row.names = F)
    return(matchedSNPs)
  }
}
