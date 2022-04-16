#!/bin/bash

# Author: Dat T Nguyen <n.dat@outlook.com>
# Date: 06 April 2022

# This script is a user friendly pipeline to use to build model for LmTag


syntax=" ./model_pipeline.sh 
    \n\n\t
    -v|--vcf : input vcf file to build  (should be in bgzip formatted file, and preprocessed by our recommended protocol: bcftools view YOUR.vcf.gz -m2 -M2 -v snps -Q 0.9999999999:major -q 0.01:minor -e 'ALT="."' | bcftools +fill-tags | sed 's/chr//g' | bgzip > YOUR_processed.vcf.gz)
    \n\n\t
    -r|--ref : path to output directory that is output of [create_imputation_ref.sh].
    \n\n\t
    -n|--size : number of tagSNP in simulated array. Theoratically, the number is need to be closed with your targeted number of tag SNP in the chromosome. We recommend this number is the maximum number of tag SNP selected by a standard greedy algorithm such as TagIt.
    \n\n\t
    -o|--out : [optional] output model file name, defaut: LmTag_model.Rdata.
    \n\n\t
    -p|--thread : [optional] number of thread - defaut: 8.
    \n
"

while [[ $# -gt 1 ]]
do
key="$1"

case $key in

     -v|--vcf) 
     in_vcf=$(readlink -f $2)
     shift
     ;;

     -r|--ref) 
     ref_dir=$(readlink -f $2)
     shift
     ;;  

     -n|--size) 
     size=$2
     shift
     ;; 

     -o|--out)
     out_fn=$2
     shift
     ;;

     -p|--thread) 
     CPUNUM=$2
     shift
     ;;
     *)

esac
shift
done

if [[ -z "$in_vcf" ]]; then
   echo ""
   echo "Usage:"
   echo -e $syntax
   echo ""
   echo "ERROR: no input vcf file !"
   echo ""
   exit
fi

if [[ -z "$ref_dir" ]]; then
   echo ""
   echo "Usage:"
   echo -e $syntax
   echo ""
   echo "ERROR: no reference m3vcf directory!"
   echo ""
   exit
fi

if [[ -z "$size" ]]; then
   echo ""
   echo "Usage:"
   echo -e $syntax
   echo ""
   echo "ERROR: no simulated array size input!"
   echo ""
   exit
fi

if [[ -z "$out_fn" ]]; then
   out_fn="LmTag_model.Rdata"
   echo "No specific setting for output name, so use the default setting (--out LmTag_model.Rdata)..."
fi

if [[ -z "$CPUNUM" ]]; then
   CPUNUM=8
   echo "No specific setting for thread number, so use the default setting (-p 8)..."
fi


start_time=$SECONDS

tem_dir=${out_fn}_model_tem
mkdir ${tem_dir}
# compute LD
#mkdir plink
plink --vcf ${in_vcf} \
      --vcf-half-call 'haploid' \
      --make-bed  --const-fid --out ${tem_dir}/tem_plink \
      --threads 1 \
      --memory 2000


plink --bfile ${tem_dir}/tem_plink \
      --r --ld-window-r2 0.2 \
      --ld-window 10000 \
      --ld-window-kb 1000 \
      --out ${tem_dir}/LD_0.2 \
      --threads $CPUNUM \
      --memory 2000

rm ${tem_dir}/tem_plink*

# LD file should be: ${tem_dir}/LD_0.2.ld
echo "Building simulated array!"
build_naive_array.R vcf=${in_vcf} size=${size} out=${tem_dir}/simulated_array.txt

echo "Imputing simulated array!"
imputation_with_prebuilt_ref.sh -t ${tem_dir}/simulated_array.txt -r ${ref_dir} -o ${tem_dir}/simulated_array -p ${CPUNUM}

echo "Estimating imputation accuracy of simulated array!"
compute_imputation_accuracy.R imputation=${tem_dir}/simulated_array out=${tem_dir}/simulated_array.Rdata

echo "Finding best tagged SNPs of simulated array!"
LmTag find --tag ${tem_dir}/simulated_array.txt --ld ${tem_dir}/LD_0.2.ld -o ${tem_dir}/find_tagged_snp_output.txt

echo "Estimating model parameters!"
buid_imputation_model.R imputation_Rdata=${tem_dir}/simulated_array.Rdata find_snp=${tem_dir}/find_tagged_snp_output.txt out_Rdata=${out_fn}

echo "Done model pipline!"

elapsed=$(( SECONDS - start_time ))
eval "echo Running time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')" > ${out_fn}.time.log
echo Running time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')

rm -r ${tem_dir}