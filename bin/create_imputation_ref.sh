#!/bin/bash

# Author: Dat T Nguyen <n.dat@outlook.com>
# Date: 15 March 2021

# This script is to use to build m3vcf reference files for leave one out cross validation.
# We assume that bcftools, vcftools, Minimac3, and Minimac4 are properly installed and can be called directly from bash shell
#
#input parameters are:
    #	-v|--vcf -- input vcf file to buil reference (should be in bgzip formatted file)
    #	-o|--out -- output directory that contains all m3vcf files, pre-created dir is not needed.
    #	-p|--thread -- number of thread.

# bash create_imputation_ref.sh -v /media/datn/data/LD_boost_paper/cross_pop_g3/Biallelic.SNP.chr10.AFR_fill_tag.vcf.gz -o AFR_ref2
###################################
syntax=" ./create_imputation_ref.sh -v [path to your vcf file] -o [path to your output dir] -p [number of thread]"

while [[ $# -gt 1 ]]
do
key="$1"

case $key in

     -v|--vcf) 
     in_ref=$(readlink -f $2)
     shift
     ;; 

     -p|--thread) 
     CPUNUM=$2
     shift
     ;; 

     -o|--out)
     ref_dir=$2
     shift
     ;;
     *)

esac
shift
done

if [[ -z "$in_ref" ]]; then
   echo ""
   echo "Usage:"
   echo $syntax
   echo ""
   echo "ERROR: no input vcf reference file !"
   echo ""
   exit
fi

if [[ -z "$ref_dir" ]]; then
	echo "No specific setting for output, so we create [mkdir ref_dir] in your current directory.."
	ref_dir="ref_dir"
fi


if [[ -z "$CPUNUM" ]]; then
   CPUNUM=8
   echo "No specific setting for thread number, so use the default setting (-p 8)..."
fi


mkdir $ref_dir

# remove chr in vcf file
zcat $in_ref | sed 's/chr//g' | bgzip > ${ref_dir}/ref.vcf.gz
# extract all sample names
bcftools view -h $in_ref | grep "^#CHROM" | cut -f10- > $ref_dir/all_sample.txt

# function to create references for minimac4
function create_impute_ref {
   in_ref=$1
   #in_SNP_array=$2
	sample=$2
   ref_dir=$3
	mkdir $ref_dir/$sample
	# sample name for reference, exclde the imputed sample...
	bcftools view -h $in_ref | grep "^#CHROM" | cut -f10- | sed 's/\t/\n/g' | grep -v "^${sample}" > $ref_dir/${sample}/sampleID_by_line.txt
	# extract reference panel
	bcftools view $in_ref  -S $ref_dir/${sample}/sampleID_by_line.txt -a -o $ref_dir/${sample}/${sample}_ref_panel.vcf.gz -Oz	
	# create m3vcf
	Minimac3 --refHaps $ref_dir/${sample}/${sample}_ref_panel.vcf.gz\
         --processReference \
         --prefix $ref_dir/${sample}/${sample}_ref_panel
   sleep 2s

   # delete to release storage
   rm $ref_dir/${sample}/${sample}_ref_panel.vcf.gz
   rm $ref_dir/${sample}/${sample}_ref_panel.erate
   rm $ref_dir/${sample}/${sample}_ref_panel.rec
}

# parallel through all samples

N=$CPUNUM
(
for i in $(cat $ref_dir/all_sample.txt); 
do 
   ((j=j%N)); ((j++==0)) && wait
   create_impute_ref ${ref_dir}/ref.vcf.gz $i $ref_dir & 
done
)
sleep 15m # make sure all process done!
echo "Done processing!"