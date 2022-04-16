#!/usr/bin/env Rscript

# Author: Dat T Nguyen <n.dat@outlook.com>
# Date: 15 March 2021

# This script is use to compute imputation accuracy between true genotypes and imputed genotype dosages. The measurement is square Person's correlation.
#input parameters are:
    # imputation=path/to/imputation_result_directory
    # out=out_name

# ./compute_imputation_accuracy.R imputation=/media/datn/data/DatProjects/vn_array/wraped_codes/imputation_naive_chr18_20000 out=test_chr18.Rdata
options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly=TRUE)
syntax='\nUsage:\t./compute_imputation_accuracy.R imputation=path/to/imputation_result_directory out=out_name\n\nInput parameters are:\n\timputation=path/to/imputation_result_directory\n\tout=out_file_name.Rdata : output file name, defaut value is: "imputation_results.Rdata"\n\n'

imputation_dir = out_name = NA

if(length(args) == 0 ){
  cat("\nNo argument, Program stop! \n")
  cat(syntax)
  quit()
}

for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))
  if (res[1] == "imputation") imputation_dir = res[2]
  #if (res[1] == "size") array_size = as.numeric(res[2])
  if (res[1] == "out") out_name = res[2]   
}

if(is.na(imputation_dir)){
  cat(syntax)
  cat("imputation_dir is not found! Program stop!\n\n")
  quit()
}

if ( is.na(out_name)) out_name = "imputation_results.Rdata"

# if(is.na(array_size)){
#   cat(syntax)
#   cat("size is not found! Program stop!\n\n")
#   quit()
# }
current_dir = getwd()
setwd(imputation_dir)

#######################################################################

ref = read.delim("ref.encoded.txt")

#x = head(ref)

info_slipt = unlist(strsplit(ref$INFO[1],";", fixed = T))
pick_pos = grep("AF=", info_slipt, fixed = T)[1]


ref$INFO = sapply(ref$INFO, function(x) unlist(strsplit(x,";"))[pick_pos])
ref$INFO = sapply(ref$INFO, function(x) unlist(strsplit(x,"="))[2])
ref$ID = paste(ref$X.CHROM, ref$POS, ref$REF, ref$ALT, sep = ":")

sample_names = colnames(ref)              
sample_names = sample_names[c(10: length(sample_names))]

## filter non-imputed SNPs
file = paste0("imputed_", sample_names[1], ".txt")
tem = read.delim(file)
pick = ref$ID %in% tem$ID
ref = ref[pick,]
## load imputed files
imputed = data.frame(ID = ref$ID)
for(sample in sample_names){
  file = paste0("imputed_", sample, ".txt")
  tem = read.delim(file)
  imputed[,sample] = tem[match(imputed$ID, tem$ID) ,sample]
}

data = ref[, c(1:8)]
x = as.matrix(ref[, sample_names])
y =  as.matrix(imputed[, sample_names])
data$r_2 = 0
for( i in 1 : nrow(data)){
  data$r_2[i] = cor(x[i,], y[i,], method = "pearson")
}
data$r_2 = data$r_2^2

data$INFO = as.numeric(data$INFO)
data$MAF = data$INFO
data$MAF[data$MAF > 0.5] = 1 - data$INFO[data$MAF > 0.5]


cutt_off = list(
  c(0.00, 0.01),
  c(0.01, 0.02),
  c(0.02, 0.03),
  c(0.03 ,0.04),
  c(0.04, 0.05),
  c(0.05, 0.075),
  c(0.075, 0.1),
  c(0.1, 0.125),
  c(0.125, 0.15),
  c(0.15, 0.175),
  c(0.175, 0.2),
  c(0.2, 0.225),
  c(0.225, 0.25),
  c(0.25, 0.275),
  c(0.275, 0.3),
  c(0.3, 0.325),
  c(0.325, 0.35),
  c(0.35, 0.375),
  c(0.375, 0.4),
  c(0.4, 0.425),
  c(0.425, 0.45),
  c(0.45, 0.475),
  c(0.475, 0.5)
)


# cutt_off = list(
#   c(0.01, 0.02),
#   c(0.02, 0.03),
#   c(0.03 ,0.04),
#   c(0.04, 0.05),
#   c(0.05, 0.075),
#   c(0.075, 0.1),
#   c(0.1, 0.125),
#   c(0.125, 0.15),
#   c(0.15, 0.2),
#   c(0.2, 0.25),
#   c(0.25, 0.3),
#   c(0.3, 0.35),
#   c(0.35, 0.4),
#   c(0.4, 0.45),
#   c(0.45, 0.5)
# )

res = c()
for (i in 1:length(cutt_off) ){
  pick = data$MAF > cutt_off[[i]][1] &  data$MAF <= cutt_off[[i]][2]
  tem = data[pick,]
  tem_res = mean(tem$r_2, na.rm = T)
  names(tem_res) = paste( cutt_off[[i]][1],  cutt_off[[i]][2], sep = ":")
  res = c(res, tem_res)
}
# back to the defaut directory
setwd(current_dir)
save(data, res, file = out_name)


cat(paste0("\nOutput: ", out_name, "\n"))
cat("Done!\n\n")
