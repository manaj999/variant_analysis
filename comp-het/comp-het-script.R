### Author: Manoj Kanagaraj
### Notes: usage compatible with ANNOVAR refgene variant annotation headers

# set wd
setwd("C:/Users/dave/Desktop/Manoj/desktop/Lab_Documents/2014_Nov")

# input file with all variants to be filtered for CH
var = read.csv("CVID_variants_ch_script.csv", header = TRUE)

# input file of the combinations of sample IDs in trios. Pedigree file
ped = read.csv("pedigree.csv", header = TRUE, stringsAsFactors=FALSE)

# temporary flag
fam_1 = 0

## PART 1: SAMPLE INDEPENDENT

# Sort variants by gene name
var[order("Gene.refgene")]

# Filter only nonsynonymous or stop variants
step1 <- subset(var, var$ExonicFunc.refgene == "nonsynonymous SNV" | var$ExonicFunc.refgene == "stoploss" | var$ExonicFunc.refgene == "stopgain")

# Filter only variants with contorls of <= 10
step2 <- subset(step1, step1$control_summaries.3..sum_controls_775max <= 8 & step1$Normal_sum_56 <= 1)



## PART 2: FAMILY SPECIFIC

# do all subsequent steps here once ready
for (i in 1:nrow(ped)) {
  # to refer to a pro, p1, p2, use ped[i,1:3]
 
  

  # Extract columns for trio (change for pro, p1 and p2)
  step3 <- step2[,c("chr","Start","Ref","Alt","CVID_sum_217","Normal_sum_56","Gene.refgene","ExonicFunc.refgene","control_summaries.3..sum_controls_775max",ped[i,1],ped[i,2],ped[i,3])]

  # create new column that IDs proband
  step3$proband <- colnames(step3)[ncol(step3)-2]

  # prep for sum
  for(j in 1:ncol(ped)){
    cur = step3[ped[i,j]]
    step3[ped[i,j]][step3[ped[i,j]] == "." | step3[ped[i,j]] == ""] <- 0
    step3[ped[i,j]] <- lapply(step3[ped[i,j]], function(x) as.numeric(as.character(x)))
  }
 
  # filter only variants with proband >=1
  step4<-subset(step3,step3[ped[i,1]]>0)
  # sum over 3 samples. could potentially change this to sum the samples specifically instead of depending upon index
  step4$triosum<-step4[ped[i,1]]+step4[ped[i,2]]+step4[ped[i,3]]
  
  # filter the three samples only for those with 2
  step5<-subset(step4,step4$triosum ==2)
  
  # filter only recurrent genes
  step6 <- subset(step5, step5$Gene.refgene %in% step5$Gene.refgene[duplicated(step5$Gene.refgene)])
  if(nrow(step6) != 0){ 
    aggp1 <- subset(aggregate(as.matrix(step6[ped[i,2]]) ~ step6$Gene.refgene, step6, FUN=sum),aggregate(as.matrix(step6[ped[i,2]]) ~ step6$Gene.refgene, step6, FUN=sum)[,2]>0)
    aggp2 <- subset(aggregate(as.matrix(step6[ped[i,3]]) ~ step6$Gene.refgene, step6, FUN=sum),aggregate(as.matrix(step6[ped[i,3]]) ~ step6$Gene.refgene, step6, FUN=sum)[,2]>0)
    
    # filter only compound heterozygous events
    step7 <- subset(step6, step6$Gene.refgene %in% aggp1[,1] & step6$Gene.refgene %in% aggp2[,1])
    
    
    # clean up output columns
    # Extract columns for trio (change for pro, p1 and p2)
    step8 <- step7[,c("chr","Start","Ref","Alt","CVID_sum_217","Normal_sum_56","Gene.refgene","ExonicFunc.refgene","control_summaries.3..sum_controls_775max",ped[i,1],ped[i,2],ped[i,3],"proband")]
    
    colnames(step8) <- c("chr","Start","Ref","Alt","CVID_sum_217","Normal_sum_56","Gene.refgene","ExonicFunc.refgene","control_summaries.3..sum_controls_775max","pro","p1","p2","proband")
    
    if(fam_1 == 0){
      output <- step8
      fam_1 = 1
    }
    else{
      output <- rbind(output,step8)
    }
  }
  
  #file_name <- paste("comp-het_out_",ped[i,1],".csv",sep="")
  #write.csv(step8, file=file_name, row.names=FALSE)
#   #write.csv(step8, file="comp-het_out.csv", row.names=FALSE)
#   if(fam_1 == 0){
#     output <- step8
#     fam_1 = 1
#   }
#   else{
#     output <- rbind(output,step8)
#   }
  # add col.names=FALSE for successive families
}

# figure out output 

# output file
write.table(output, file = "comp-het_out.csv",row.names=FALSE, sep=",")
# write.csv(TADA_gene_list, file="Varscan_All-Trios_2_prepped-vars.csv")

