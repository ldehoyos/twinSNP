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
#' It is different to the function "SNPmatch" as this function outputs a
#' list with the information of the SNPs and its matches, which is also done
#' with the functions "listMatches" and "listMatches2". Thus, this function
#' is better in terms of computation time.
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
#' SNPlist <- SNPmatch_list(database = myDatabase, inputSNPs = myVars, MAF = 0.05, GD = 0.5, DNG = 0.5, LDB = 0.5)
#'
#' @import data.table
#'
#' @export
SNPmatchXList <- function(database, inputSNPs, MAF, GD, DNG, LDB){
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
  maf_min <- included_snps$snp_maf - MAF
  maf_max <- included_snps$snp_maf + MAF

  # Set the gene density interval
  gdens_min <- included_snps$gene_count - GD*included_snps$gene_count
  gdens_max <- included_snps$gene_count + GD*included_snps$gene_count

  # Set the distance to nearest gene interval
  dist_interval_min <- included_snps$dist_nearest_gene_snpsnap - DNG*included_snps$dist_nearest_gene_snpsnap
  dist_interval_max <- included_snps$dist_nearest_gene_snpsnap + DNG*included_snps$dist_nearest_gene_snpsnap

  # Set the LD interval
  ld_interval_min <- included_snps$friends_ld05 - LDB*included_snps$friends_ld05
  ld_interval_max <- included_snps$friends_ld05 + LDB*included_snps$friends_ld05

  # Do the analysis and matching of SNPs
  matched_snps <- SNPinfo <- list()

  for (i in 1:nrow(included_snps)){

    # Do the filtering to select the variants that meet the requirements.
    A1 <- datB[data.table::between(snp_maf, lower= maf_min[i], upper= maf_max[i])]
    A2 <-   A1[data.table::between(gene_count, lower= gdens_min[i], upper= gdens_max[i])]
    A3 <-   A2[data.table::between(dist_nearest_gene_snpsnap, lower= dist_interval_min[i], upper=dist_interval_max[i])]
    subSnp <- A3[data.table::between(friends_ld05, lower= ld_interval_min[i], upper=ld_interval_max[i])]
    rownames(subSnp) <- c()

    if(file.exists("matched_snps_annotation.txt")== T){file.remove("matched_snps_annotation.txt")}
    if(file.exists("input_snps_annotated_unmatched.txt")== T){file.remove("input_snps_annotated_unmatched.txt")}

    matched_snps <- data.frame(input_snp= included_snps$snpID[i], subSnp, stringsAsFactors = F)
    if (is.na(matched_snps$snpID)){
      write.table(matched_snps, "input_snps_annotated_unmatched.txt", append=T, quote=F, col.names=T, row.names=F)
    }else{
      write.table(matched_snps, "matched_snps_annotation.txt", append=T, quote=F, col.names=T, row.names=F)
    }
    SNPinfo[[i]] <- subSnp$snpID
  }
  names(SNPinfo) <- included_snps$snpID

  return(SNPinfo)
}
