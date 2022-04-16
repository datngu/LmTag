#!/usr/bin/env Rscript

# Author: Dat T Nguyen <n.dat@outlook.com>
# Date: 15 March 2021

# This script is create a naive SNP slection set with pre-set number of target SNP. Procedure: 1. Ranking all SNP positions 2. Pick SNPs by odering with pre-set total number. Reading the implemented codes for more details.
#input parameters are:
    # vcf=path/to/vcf.gz
    # size=number_of_tagSNP : for example: size=20000
    # out=out_file_name.txt : output file name, defaut value is: "naive_array_size.txt"

# ./build_naive_array.R vcf=/media/datn/data/DatProjects/vn_array/data/vn_504_maf_0.01/MAF_0.01.chr22_vn504_Biallelic.SNP_PHASED_SHAPEIT4_no_chr_version_dbSNP_151.vcf.gz size=20000 out=chr18_naive_20000.txt

options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly=TRUE)
syntax='\nUsage:\t./build_naive_array.R vcf=path/to/vcf.gz size=size out=output_name.txt\n\nInput parameters are:\n\tvcf=path/to/vcf.gz\n\tsize=number_of_tagSNP : for example: size=20000\n\tout=out_file_name.txt : output file name, defaut value is: "naive_array_size.txt"\n\n'

vcf_path = array_size = out_name = NA

if(length(args) == 0 ){
  cat("\nNo argument, Program stop! \n")
  cat(syntax)
  quit()
}

for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))
  if (res[1] == "vcf") vcf_path = res[2]
  if (res[1] == "size") array_size = as.numeric(res[2])
  if (res[1] == "out") out_name = res[2]   
}

if(is.na(vcf_path)){
  cat(syntax)
  cat("vcf is not found! Program stop!\n\n")
  quit()
}

if(is.na(array_size)){
  cat(syntax)
  cat("size is not found! Program stop!\n\n")
  quit()
}

# extract snp positions
all_site = paste0(vcf_path, "_all_site.txt")
command = paste0("zcat ", vcf_path, " | grep -v ^# | awk '{print $1, $2}' > ", all_site)
system(command)

site = read.table(all_site, header = F)

all_snp = site$V2
n_snp = nrow(site)
all_snp = sort(all_snp)
array_size = array_size + 1
step = n_snp / array_size
pick = round(seq(1, n_snp, step))
pick = pick[-1]
# pick = sample(c(1:nrow(site)), array_size) ## random choosing - not use
array = site[pick,]

if ( is.na(out_name)) out_name = paste0( "naive_array_", (array_size-1), ".txt")

write.table(array , file = out_name, col.names = F, row.names = F, sep ="\t", quote = F)

cat(paste0("\nOutput: ", out_name, "\n"))
cat("Done!\n")