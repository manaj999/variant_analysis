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
step4$CVID1156[step4$CVID1156 == ""] <- 0 # fix step specification
step4$CVID1156[step4$CVID1156 == "."] <- 0

step4$CVID1156 <- as.numeric(levels(step4$CVID1156))[step4$CVID1156] #needed in order to sum afterwards

## same for other two samples .. change to step3
step4$X_1738_CVID[step4$X_1738_CVID == ""] <- 0
step4$X_1738_CVID[step4$X_1738_CVID == "."] <- 0
step4$X_1739_CVID[step4$X_1739_CVID == ""] <- 0
step4$X_1739_CVID[step4$X_1739_CVID == "."] <- 0
step4$X_1738_CVID <- as.numeric(levels(step4$X_1738_CVID))[step4$X_1738_CVID]
step4$X_1739_CVID <- as.numeric(levels(step4$X_1739_CVID))[step4$X_1739_CVID]
# do the above 3 steps for all three samples, not one at a time..


# filter only variants with proband >=1, change step numbers
step4<-subset(step3,step3[,10]>0)

# sum over 3 samples. could potentially change this to sum the samples specifically instead of depending upon index
step4$triosum<-step4[,10]+step4[,11]+step4[,12]
step4$partial_sum<-step4[,10]+step4[,11]

# filter the three samples only for those with 2
step5<-subset(step4,step4$triosum ==2)

# filter only recurrent genes
step6 <- subset(step5, step5$Gene.refgene %in% step5$Gene.refgene[duplicated(step5$Gene.refgene)])

# filter only recurrent genes with a hit in both parents. if this works, then can get rid of partial_sum
aggregate(step7$X_1738_CVID ~ step7$Gene.refgene, step7, FUN=sum)
aggregate(step7$X_1739_CVID ~ step7$Gene.refgene, step7, FUN=sum)

aggp1 <- subset(aggregate(step7$X_1739_CVID ~ step7$Gene.refgene, step7, FUN=sum),aggregate(step7$X_1739_CVID ~ step7$Gene.refgene, step7, FUN=sum)[,2]>0)[,1]
aggp2 <- subset(aggregate(step7$X_1738_CVID ~ step7$Gene.refgene, step7, FUN=sum),aggregate(step7$X_1738_CVID ~ step7$Gene.refgene, step7, FUN=sum)[,2]>0)

step8 <- subset(step7, step7$Gene.refgene %in% aggp1[,1] & step7$Gene.refgene %in% aggp2[,1])

# create new column that sums across pro, p1, p2
step5 <- apply(step4,1,function(row) sum(vec[ row[ncol(step4-2)] : row[ncol(step4)] ] ))

for (i in 1:nrow(ped)) {
  df[,c(ped$pro,ped$p1,ped$p2)]  
  # do all subsequent steps here once ready
  

# Extract columns for trio (change for pro, p1 and p2)
  step3 <- step2[,c("chr","Start","Ref","Alt","CVID_sum_217","Normal_sum_56","Gene.refgene","ExonicFunc.refgene","control_summaries.3..sum_controls_775max","CVID1156","X_1738_CVID","X_1739_CVID")]


# output file
# write.csv(TADA_gene_list, file="Varscan_All-Trios_2_prepped-vars.csv")

