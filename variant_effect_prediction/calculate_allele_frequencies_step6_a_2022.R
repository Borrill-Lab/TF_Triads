#'Calculate allele frequencies from He et al data
#'05/08/20
#'Exome_capture_preparation_050820.R
#'Prepare summary tables which record the number of variants of each type
#'by both continent and cultivar/landrace. REMEMBER to keep the record of
#'number of 1/1, 1/0 and 0/0 lines so that MAF can be recalculated from
#'different groupings.
#'
#'20/08/20
#'EDIT: add 'cultivar_col' with combined Web and Supp groupings, and remove an NA
#'
#'17/08/21
#'Pared down massively
#'Calculate overall allele frequency only (no splitting by continent/cultivar), for large-scale analysis
#'Keep ID column, as this matches the VEP data
#'
#'18/08/21
#'Tried to streamline script using pipe so we don't cache multiple copies of the table.
#'
#'10/01/22
#'Filter for varieties in the middle 95% for mutation load, so that odd individuals with high load
#'don't skew the results and add a disproportionate number of SNPs to the graph
#'
#'21/01/22
#'Recalculate with extremeload2.txt
#'Calculate MINOR allele frequency, not REFERENCE based allele frequency
#'
#'Calculate allele frequency after this filter.


#RUN ONCE
require(readr)
require(reshape2)
require(lubridate)
require(dplyr)

#FUNCTION
#Calculate allele frequencies
#NOTE FORWARD SLASH rather than pipe, / not |
#NOTE I have got rid of `1/0` as it doesn't appear in this data format
return.MAF <- function(GenFreq){
  ###Returns the allele frequency of the less common allele in the population ###
  p <- (GenFreq$`1/1`+0.5*GenFreq$`0/1`)/(GenFreq$`0/0`+GenFreq$`0/1`+GenFreq$`1/1`)
  q <- 1-p
  allele.freq <- pmin(p, q) #pmin calculates parallel minima for two vectors
  return(allele.freq)
}

#RUN SECTION



#Variants with >25% missing values have NOT been filtered out.
#Choose a table containing ONLY variants of interest.
setwd("/rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs/4_analysis")
VarTable <- read_delim("He19F_coding_variant_in_triad_genes.vcf", delim = "\t")
colnames(VarTable)[1] <- "CHROM"

#FILTER OUT EXTREME MUTATION LOAD
setwd("/rds/projects/b/borrillp-phd-students/Catherine/expression_and_load_data")
extreme_load <- read_lines("extremeload2.txt")
extreme_load <- gsub("\\.", "-", extreme_load) #reformat to match
VarTable_v1.5 <- VarTable[which(!(colnames(VarTable) %in% extreme_load))]

print("before:")
print(ncol(VarTable))
print("after:")
print(ncol(VarTable_v1.5))

#fun.aggregate = length automatically counts the number of records in each category
#let's hope these commands work with the pipe
VarTable_v2 <- VarTable_v1.5 %>%
  melt(id.vars = c("CHROM", "POS", "ID", "REF", "ALT","QUAL", "FILTER", "INFO", "FORMAT")) %>% #All non-ID vars go in variable column
  dcast(CHROM + POS + ID + REF + ALT ~ value, fun.aggregate = length)

VarTable_v3 <- VarTable_v2 %>%
  mutate(MAF=return.MAF(VarTable_v2))

#Save everything we've done
setwd("/rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs/4_analysis")
write_csv(VarTable_v3, file = paste("allele_frequency_coding_variant_extreme_load_",today(),".csv", sep=""))
