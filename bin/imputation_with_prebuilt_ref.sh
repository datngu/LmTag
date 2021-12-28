#!/bin/bash

# Author: Dat T Nguyen <n.dat@outlook.com>
# Date: 15 March 2021

# This script is to use to run leave one out imputation with pre-built m3vcf reference files, input is path to output directory that created by [create_imputation_ref.sh].
# We assume that bcftools, vcftools, Minimac3, and Minimac4 are properly installed and can be called directly from bash shell
#
#input parameters are:
    # -t|--tag -- input is a tab separated file, with two columns contain chr number and position of SNP and the file should have no header.  
    # -r|--ref -- path to output directory that is output of [create_imputation_ref.sh].
    # -o|--out -- output directory that contains imputation results, pre-created dir is not needed.
    # -p|--thread -- number of thread.

# /media/datn/data/DatProjects/vn_array/wraped_codes/imputation_with_prebuilt_ref.sh -t /media/datn/data/DatProjects/vn_array/build_training_chr18/train_7.txt -r /media/datn/data/DatProjects/vn_array/reference_imputation/chr_18_ref_Biallelic_SNP_only_maf_0.01 -o testing_imputation

###################################
syntax=" ./imputation_with_prebuilt_ref.sh -t [path to your tagSNP file] -r [path to your reference directory (output of create_imputation_ref.sh) ] -o [path to your output dir] -p [number of thread]"

while [[ $# -gt 1 ]]
do
key="$1"

case $key in

     -t|--tag) 
     in_tag=$(readlink -f $2)
     shift
     ;;

     -r|--ref) 
     ref_dir=$(readlink -f $2)
     shift
     ;;  

     -p|--thread) 
     CPUNUM=$2
     shift
     ;; 

     -o|--out)
     res_dir=$2
     shift
     ;;
     *)

esac
shift
done

if [[ -z "$in_tag" ]]; then
   echo ""
   echo "Usage:"
   echo $syntax
   echo ""
   echo "ERROR: no tagSNP file !"
   echo ""
   exit
fi

if [[ -z "$ref_dir" ]]; then
   echo ""
   echo "Usage:"
   echo $syntax
   echo ""
   echo "ERROR: no reference directory !"
   echo ""
   exit
fi

#####################

if [[ -z "$res_dir" ]]; then
  echo "No specific setting for output, so we create [mkdir res_dir] in your current directory.."
  res_dir="res_dir"
fi


if [[ -z "$CPUNUM" ]]; then
   CPUNUM=8
   echo "No specific setting for thread number, so use the default setting (-p 8)..."
fi


########################

mkdir $res_dir
vcftools --gzvcf $ref_dir/ref.vcf.gz --positions $in_tag --recode --stdout | bgzip > $res_dir/array.vcf.gz



function leave_one_out_imputation {
	in_SNP_array=$1
	sample=$2
   ref_dir=$3
   res_dir=$4
   mkdir $res_dir/$sample
	bcftools view $in_SNP_array -s $sample -a -o $res_dir/$sample/${sample}_SNP_array.vcf.gz -Oz
	minimac4 --refHaps $ref_dir/${sample}/${sample}_ref_panel.m3vcf.gz \
         --ChunkLengthMb 200 \
         --ChunkOverlapMb 20 \
         --haps $res_dir/$sample/${sample}_SNP_array.vcf.gz \
         --prefix $res_dir/imputed_${sample} \
         --ignoreDuplicates \
         --cpus 1
}

# parallel through all samples

N=$CPUNUM
(
for i in $(cat $ref_dir/all_sample.txt); 
do 
   ((j=j%N)); ((j++==0)) && wait
   leave_one_out_imputation $res_dir/array.vcf.gz $i $ref_dir $res_dir & 
done
)

sleep 40s
# process ref to read in R
cd $res_dir
zcat $ref_dir/ref.vcf.gz | grep -v "##" | sed 's/0|0/0/g' | sed 's/0|1/1/g' | sed 's/1|0/1/g' | sed 's/1|1/2/g' > ref.encoded.txt
rm *.info
rm -r `ls -d */`


for file in *dose.vcf.gz
do
	out=`echo $file | cut -d "." -f1`
	zcat $file | grep -v "^##" | awk '//{printf "%s\t%s\n", $3, $10}' | sed 's/0|0://g' | sed 's/0|1://g' | sed 's/1|0://g' | sed 's/1|1://g' > ${out}.txt
done

rm *dose.vcf.gz

echo "Done processing!"


