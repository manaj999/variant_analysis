### Author: Manoj Kanagaraj
### Notes: usage compatible with ANNOVAR refgene variant annotation headers

# set wd
setwd("C:/Users/dave/Desktop/Manoj/desktop/Lab_Documents/2014_Nov")

# input file with all variants to be filtered for CH
var = read.csv("CVID_variants_ch_script.csv", header = TRUE)

# input file of the combinations of sample IDs in trios. Pedigree file
ped = read.csv("pedigree.txt", header = TRUE, stringsAsFactors=FALSE)


## PART 1: SAMPLE INDEPENDENT

# Sort variants by gene name
var[order("Gene.refgene")]

# Filter only nonsynonymous or stop variants
step1 <- subset(var, var$ExonicFunc.refgene == "nonsynonymous SNV" | var$ExonicFunc.refgene == "stoploss" | var$ExonicFunc.refgene == "stopgain")

# Filter only variants with contorls of <= 10
step2 <- subset(step1, step1$control_summaries.3..sum_controls_775max <= 8 & step1$Normal_sum_56 <= 1)

#head(step2)

## PART 2: FAMILY SPECIFIC

#step3 <- step2[,c("chr","Start","Ref","Alt","CVID_sum_217","Normal_sum_56","Gene.refgene","ExonicFunc.refgene","control_summaries.3..sum_controls_775max","CVID1156","X_1738_CVID","X_1739_CVID")]


#step3$proband <- colnames(step4)[ncol(step4)-3]

# prep for sum
  # convert blanks and dots to 0
# step4$CVID1156[step4$CVID1156 == ""] <- 0 # fix step specification
# step4$CVID1156[step4$CVID1156 == "."] <- 0
# 
# step4$CVID1156 <- as.numeric(levels(step4$CVID1156))[step4$CVID1156] #needed in order to sum afterwards

## same for other two samples .. change to step3
# step4$X_1738_CVID[step4$X_1738_CVID == ""] <- 0
# step4$X_1738_CVID[step4$X_1738_CVID == "."] <- 0
# step4$X_1739_CVID[step4$X_1739_CVID == ""] <- 0
# step4$X_1739_CVID[step4$X_1739_CVID == "."] <- 0
# step4$X_1738_CVID <- as.numeric(levels(step4$X_1738_CVID))[step4$X_1738_CVID]
# step4$X_1739_CVID <- as.numeric(levels(step4$X_1739_CVID))[step4$X_1739_CVID]
# do the above 3 steps for all three samples, not one at a time..


# filter only variants with proband >=1, change step numbers
#step4<-subset(step3,step3[,10]>0)

# sum over 3 samples. could potentially change this to sum the samples specifically instead of depending upon index
#step4$triosum<-step4[,10]+step4[,11]+step4[,12]
#step4$partial_sum<-step4[,10]+step4[,11]

# filter the three samples only for those with 2
#step5<-subset(step4,step4$triosum ==2)

# filter only recurrent genes
#step6 <- subset(step5, step5$Gene.refgene %in% step5$Gene.refgene[duplicated(step5$Gene.refgene)])

# filter only recurrent genes with a hit in both parents. if this works, then can get rid of partial_sum
#aggregate(step7$X_1738_CVID ~ step7$Gene.refgene, step7, FUN=sum)
#aggregate(step7$X_1739_CVID ~ step7$Gene.refgene, step7, FUN=sum)

#aggp1 <- subset(aggregate(step7$X_1739_CVID ~ step7$Gene.refgene, step7, FUN=sum),aggregate(step7$X_1739_CVID ~ step7$Gene.refgene, step7, FUN=sum)[,2]>0)
#aggp2 <- subset(aggregate(step7$X_1738_CVID ~ step7$Gene.refgene, step7, FUN=sum),aggregate(step7$X_1738_CVID ~ step7$Gene.refgene, step7, FUN=sum)[,2]>0)

#step8 <- subset(step7, step7$Gene.refgene %in% aggp1[,1] & step7$Gene.refgene %in% aggp2[,1])




## PART 2: FAMILY SPECIFIC
# do all subsequent steps here once ready
for (i in 1:nrow(ped)) {
   
  # to refer to a pro, p1, p2, use ped[i,1:3]
 
  

  # Extract columns for trio (change for pro, p1 and p2)
  step3 <- step2[,c("chr","Start","Ref","Alt","CVID_sum_217","Normal_sum_56","Gene.refgene","ExonicFunc.refgene","control_summaries.3..sum_controls_775max",ped[1,1],ped[1,2],ped[1,3])]

  # create new column that IDs proband
  step3$proband <- colnames(step3)[ncol(step3)-2]

  # prep for sum
  for(j in 1:ncol(ped)){
    cur = step3[ped[i,j]]
    step3[ped[i,j]][step3[ped[i,j]] == "." | step3[ped[i,j]] == ""] <- 0
    step3[ped[i,j]] <- lapply(step3[ped[i,j]], function(x) as.numeric(as.character(x)))
    #cur <- as.numeric(levels(cur))[cur]
    #step3[ped[i,j]][step3[ped[i,j]] == "." | step3[ped[i,j]] == ""] <- 0
    #step3[ped[i,j]] <- as.numeric(levels(step3[ped[i,j]]))[step3[ped[i,j]]]
  }
 
  # filter only variants with proband >=1
  step4<-subset(step3,step3[ped[i,1]]>0)
  # sum over 3 samples. could potentially change this to sum the samples specifically instead of depending upon index
  step4$triosum<-step4[ped[i,1]]+step4[i,2]+step4[i,3]
  
  # filter the three samples only for those with 2
  step5<-subset(step4,step4$triosum ==2)
  
  # filter only recurrent genes
  step6 <- subset(step5, step5$Gene.refgene %in% step5$Gene.refgene[duplicated(step5$Gene.refgene)])
  
  aggp1 <- subset(aggregate(step6[ped[i,2]] ~ step6$Gene.refgene, step6, FUN=sum),aggregate(step6[ped[i,2]] ~ step6$Gene.refgene, step6, FUN=sum)[,2]>0)
  aggp2 <- subset(aggregate(step6[ped[i,3]] ~ step6$Gene.refgene, step6, FUN=sum),aggregate(step6[ped[i,3]] ~ step6$Gene.refgene, step6, FUN=sum)[,2]>0)
  
  # filter only compound heterozygous events
  step7 <- subset(step6, step6$Gene.refgene %in% aggp1[,1] & step6$Gene.refgene %in% aggp2[,1])
}

# figure out output 

# output file
# write.csv(TADA_gene_list, file="Varscan_All-Trios_2_prepped-vars.csv")

