### Author: Manoj Kanagaraj
### Notes: usage compatible with ANNOVAR refgene variant annotation headers

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

# create new column that IDs proband
step3$proband <- colnames(step4)[ncol(step4)-3]

# prep for sum
  # convert blanks and dots to 0
step3$CVID1156[step3$CVID1156 == ""] <- 0 # fix step specification
step3$CVID1156[step3$CVID1156 == "."] <- 0

step3$CVID1156 <- as.numeric(levels(step3$CVID1156))[step3$CVID1156] #needed in order to sum afterwards

## same for other two samples .. change to step3
step4$X_1738_CVID[step4$X_1738_CVID == ""] <- 0
step4$X_1738_CVID[step4$X_1738_CVID == "."] <- 0
step4$X_1739_CVID[step4$X_1739_CVID == ""] <- 0
step4$X_1739_CVID[step4$X_1739_CVID == "."] <- 0
step4$X_1738_CVID <- as.numeric(levels(step4$X_1738_CVID))[step4$X_1738_CVID]
step4$X_1739_CVID <- as.numeric(levels(step4$X_1739_CVID))[step4$X_1739_CVID]
# do the above 3 steps for all three samples, not one at a time..

step4$triosum<-step4[,10]+step4[,11]+step4[,12]

# create new column that sums across pro, p1, p2
step5 <- apply(step4,1,function(row) sum(vec[ row[ncol(step4-2)] : row[ncol(step4)] ] ))

for (i in 1:nrow(ped)) {
  df[,c(ped$pro,ped$p1,ped$p2)]  
  # do all subsequent steps here once ready
  

# Extract columns for trio (change for pro, p1 and p2)
  step3 <- step2[,c("chr","Start","Ref","Alt","CVID_sum_217","Normal_sum_56","Gene.refgene","ExonicFunc.refgene","control_summaries.3..sum_controls_775max","CVID1156","X_1738_CVID","X_1739_CVID")]


# output file
# write.csv(TADA_gene_list, file="Varscan_All-Trios_2_prepped-vars.csv")

