#!/bin/bash

# Author: Dat T Nguyen <n.dat@outlook.com>
# Date: 09 April 2022

# This script is to use to pre-processing vcf, used in docker container

syntax="./process_vcf.sh in.vcf.gz out_vcf.gz"
in_vcf=$1
out_vcf=$2

if [[ -z "$in_vcf" ]]; then
   echo ""
   echo "Usage:"
   echo $syntax
   echo ""
   echo "ERROR: no input vcf file !"
   echo ""
   exit
fi

if [[ -z "$out_vcf" ]]; then
   echo ""
   echo "Usage:"
   echo $syntax
   echo ""
   echo "ERROR: no output vcf file !"
   echo ""
   exit
fi

echo PROCESSING your input vcf.gz file : $in_vcf!

bcftools view $1 -m2 -M2 -v snps -Q 0.9999999999:major -q 0.01:minor -e 'ALT="."' | bcftools +fill-tags | sed 's/chr//g' | bgzip > $2

echo DONE! The processed vcf.gz file is: $out_vcf!

