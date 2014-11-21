# set wd
setwd("C:/Users/dave/Desktop/Manoj/desktop/Lab_Documents/2014_Nov")

# input file with all variants to be filtered for CH
var = read.csv("CVID_variants_ch_script.csv", header = TRUE)

# input file of the combinations of sample IDs in trios. Pedigree file
ped = read.csv("pedigree.txt", header = TRUE)


## PART 1: SAMPLE INDEPENDENT

# Sort variants by gene name
var[order("Gene.refgene")]

# Filter only nonsynonymous or stop variants
step1 <- subset(var, var$ExonicFunc.refgene == "nonsynonymous SNV" | var$ExonicFunc.refgene == "stoploss" | var$ExonicFunc.refgene == "stopgain")

# Filter only variants with contorls of <= 10
step2 <- subset(step1, step1$control_summaries.3..sum_controls_775max <= 10)

head(step2)

## PART 2: FAMILY SPECIFIC

step3 <- step2[,c("chr","Start","Ref","Alt","CVID_sum_217","Normal_sum_56","Gene.refgene","ExonicFunc.refgene","control_summaries.3..sum_controls_775max","CVID1156","X_1738_CVID","X_1739_CVID")]


step4 <-

for (i in 1:nrow(ped)) {
  df[,c(ped$pro,ped$p1,ped$p2)]  
  # do all subsequent steps here once ready
  

# Extract columns for trio (change for pro, p1 and p2)
  step3 <- step2[,c("chr","Start","Ref","Alt","CVID_sum_217","Normal_sum_56","Gene.refgene","ExonicFunc.refgene","control_summaries.3..sum_controls_775max","CVID1156","X_1738_CVID","X_1739_CVID")]


# output file
# write.csv(TADA_gene_list, file="Varscan_All-Trios_2_prepped-vars.csv")

