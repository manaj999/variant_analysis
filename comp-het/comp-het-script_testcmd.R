### Author: Manoj Kanagaraj
### Notes: usage compatible with ANNOVAR refgene variant annotation headers

### EXAMPLE CMD: R comp-het-script <ARGUMENT 1> <ARGUMENT 2> <ARGUMENT 3> <ARGUMENT 4>

### COMMAND LINE ARGUMENTS:
## ARGUMENT 1: Working Directory
## ARGUMENT 2: Variant File (headers in same order as example file)
## ARGUMENT 3: Pedigree File
## ARGUMENT 4: Output File

## SET FILTER PARAMETERS HERE
# Control variants from 6500 exomes, etc.
control_max = 0

# Variants from all patients' parents
parent_max = 1 # should be at least one since it is a heterozygous variant

# SET WORKING DIRECTORY
# setwd("C:/Users/dave/Desktop/Manoj/desktop/Lab_Documents/2014_Nov")
print(commandArgs()[1])
print(commandArgs()[2])
print(commandArgs()[3])
print(commandArgs()[4])

# SET INPUT FILE CONTAINING ANNOVAR OUTPUT
var = read.csv("CVID_variants_ch_script.csv", header = TRUE)

# SET PEDIGREE FILE
ped = read.csv("pedigree.csv", header = TRUE, stringsAsFactors=FALSE)
