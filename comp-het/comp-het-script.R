### Author: Manoj Kanagaraj
### Notes: Usage compatible with ANNOVAR refgene variant annotation headers [ chr Start Ref Alt Case_sum  Normal_sum  Gene.refgene  ExonicFunc.refgene  Control_sum ... Sample_IDs ]
### Notes: Additional filters may be used to pre-select variants for this script. Otherwise, this script only filters variant type and frequency based on normals and controls

### EXAMPLE CMD: R comp-het-script <ARGUMENT 1> <ARGUMENT 2> <ARGUMENT 3> <ARGUMENT 4>

### COMMAND LINE ARGUMENTS:
## ARGUMENT 1: Working Directory
## ARGUMENT 2: Variant File (headers in same order as example file)
## ARGUMENT 3: Pedigree File
## ARGUMENT 4: Output File

## SET FILTER PARAMETERS HERE
# Control variants from 6500 exomes, etc.
control_max = 10

# Variants from all patients' parents
parent_max = 10 # should be at least one since it is a heterozygous variant

# SET WORKING DIRECTORY
# setwd("PATH to WD")

# SET INPUT FILE CONTAINING ANNOVAR OUTPUT
var = read.csv("rare-var_for_comp-het.csv", header = TRUE, check.names = FALSE)

# SET PEDIGREE FILE
ped = read.csv("pedigree_01-22.csv", header = TRUE, stringsAsFactors=FALSE, check.names = FALSE)

# SET OTHER VARIABLES
output <- data.frame()

## PART 1: SAMPLE INDEPENDENT

# Filter only nonsynonymous or stopgain/stoploss variants
# **Edit this step if you wish to include additional variant types**
step1 <- var[var$ExonicFunc.refgene == "nonsynonymous SNV" | var$ExonicFunc.refgene == "stoploss" | var$ExonicFunc.refgene == "stopgain",,drop=F]
stopifnot(nrow(step1)>0)


# Filter variants based on rareness using parameters set above
step2 <- subset(step1, step1[,'Control_sum'] <= control_max & step1[,'Normal_sum'] <= parent_max)
stopifnot(nrow(step2)>0)


## PART 2: FAMILY SPECIFIC

# Repeat for every family to aggregate compound heterozygous events
for (i in 1:nrow(ped)) {

  # Extract columns for current trio
  header = colnames(step2)
  step3 <- step2[,c(header[1:9],ped[i,1],ped[i,2],ped[i,3])]
  
  # Create new column to identify trio
  step3$proband <- colnames(step3[,ped[i,1],drop=F])

  # Clean up data before summing in case there are any blanks (. or "")
  #step3[step3[,ped[i,]]=="."|is.na(step3[,ped[i,]]), c(ped[i,])] <- 0

  for(j in 1:ncol(ped)){
    cur = step3[ped[i,j]]
    step3[ped[i,j]][step3[ped[i,j]] == "." | step3[ped[i,j]] == ""] <- 0
    step3[ped[i,j]] <- lapply(step3[ped[i,j]], function(x) as.numeric(as.character(x)))
  }
 
  
  # Sort variants by gene name
  step3 <- step3[order(step3$Gene.refgene),]
  
  # Filter only variants that occur in proband
  step4<-subset(step3,step3[ped[i,1]]>0)
  if(nrow(step4)==0) {next}
  
  # Sum over the three samples
  step4$triosum<-step4[ped[i,1]]+step4[ped[i,2]]+step4[ped[i,3]] #apply
  
  # Filter the three samples only for those with 2. This suggests that it is a heterozygous variant.
  step5<-subset(step4,step4$triosum ==2)
  if(nrow(step5)==0) {next}
  
  # Filter only recurrent genes
  step6 <- subset(step5, step5$Gene.refgene %in% step5$Gene.refgene[duplicated(step5$Gene.refgene)])
  if(nrow(step6)==0) {next}
  
    
  # Extract genes with heterozygous variants in each parent
  aggp1 <- subset(aggregate(as.matrix(step6[ped[i,2]]) ~ step6$Gene.refgene, step6, FUN=sum), aggregate(as.matrix(step6[ped[i,2]]) ~ step6$Gene.refgene, step6, FUN=sum)[,2]>0)
  aggp2 <- subset(aggregate(as.matrix(step6[ped[i,3]]) ~ step6$Gene.refgene, step6, FUN=sum), aggregate(as.matrix(step6[ped[i,3]]) ~ step6$Gene.refgene, step6, FUN=sum)[,2]>0)
  
  # Filter for genes that have recurrent heterozygous events, with at least one from each parent. This suggests that it is a heterozygous variant.
  step7 <- subset(step6, step6$Gene.refgene %in% aggp1[,1] & step6$Gene.refgene %in% aggp2[,1])
  
  
  # Extract relevant columns (exclude "triosum" column)
  step7 <- step7[ - c(14)]
  
  colnames(step7)[10:12] <- c("pro","p1","p2")
  
  # Use output df to aggregate compound heterozygous events from all families
  output <- rbind(output,step7)
}  

# Output file
write.table(output, file = "comp-het_out.csv", row.names=FALSE, sep=",")
