# Documents for LmTag

###########################################################
## Update news
#### 28 December 2021: v0.1.0
#### 31 August 2021: testing version v0.0.2
- Add tagged marker positions column in output.
#### 30 July 2021: testing version v0.0.1
- Add --vip and --exclude parameters: use to add VIP SNPS list and excluded SNPS list for the program.
- VIP SNPs are prioritized to be selected by the algorithm while excluded SNPs are elimiated from the selection.
__NOTE__ not sure that all VIP SNPs are selected and all excluded SNPs are excluded because the program weight imputation is the primary factor of choosing SNPs.

#### 23 June 2021: testing version 0.0.0
- First submission
## 1. Introduction

`LmTag` is a model based method to find tagSNP in SNP array desgin that maximize imputation coverage and functional score of tag SNPs. Full details of the method is described in the method manuscript.
 
### Software requirements
`LmTag` is implemented in `R` and `C++`. `LmTag` requires `vcftools, bcftools, minimac3, minimac4, plink` for model construction step.
- `VCFtools 0.1.17` : https://vcftools.github.io
- `bcftools 1.10.2` : https://github.com/samtools/bcftools
- `plink1.9` : https://www.cog-genomics.org/plink/
- `minimac3` : https://genome.sph.umich.edu/wiki/Minimac3
- `minimac4` : https://genome.sph.umich.edu/wiki/Minimac4
- A `C++-11` compliant compiler version of `GCC (g++ >= 4.8.2)`
- `R` packages version `3.6.0` or latter

## 2. Download and installation
### 2.1 Download
### 2.2 Installation

```sh
# unzip the software package
tar -xzvf LmTag.tar.gz
## Move to the *LmTag* directory and do configuration for LmTag
cd LmTag
# build the C++ program
make
# export to PATH
export PATH=$PWD/bin:$PATH
cd ..
```
The installation assumes that `vcftools, bcftools, minimac3, minimac4, plink` are available and can be call directly from your terminal promt.

If `vcftools, bcftools, minimac3, minimac4, plink` are not available, please install and add these tools to `PATH` variable.
__OTHERWISE, THERE WILL BE ERRORS__

## 3. Input requirement
### LmTag require a phased vcf file to perform tag SNP selection with following criteria:
- Only biallelic SNPS are considered.
- `MAF` are carefully check before runing (recommend to use `filltag` command by `bcftools` before input to LmTag pipeline) as LmTag extracts MAF directly from `INFO/AF` information in vcf file.
- Minimum `MAF` threshold are pre-determine in vcf filer, so you need to do filtering before providing the file to `LmTag` pipeline.
- `#CHROM` should be encoded without 'chr' character.
### Other possible input files for `LmTag`:
- Functional scores of SNPs candidates, typically extracted from the `CADD` databases [optional]
- List of VIP SNPs - they will be prioritized in tag SNP selection at highest levels; typically SNPs in `GWAS catalog` or `ClinVar` databases that you really want to select as tag SNPs [optional]
- List of bad SNPs - they will be tried to exluded in tag SNP selection; typically SNPs in repeated regions or low quality that you don't want selected as tag SNPs [optional]

## 4. Step by step instruction to run LmTag
### 4.0 Data reprocessing and compute needed information
#### Obtain tutorial data:
We provide in this tutorial based on chromosome 10, __East Asian population__ dataset:
- vcf file: https://zenodo.org/api/files/efb7a8bc-efce-4391-8f65-59b974328c2e/chr10_EAS.vcf.gz
- CADD score file: https://zenodo.org/api/files/efb7a8bc-efce-4391-8f65-59b974328c2e/chr10_EAS_CADD.txt - 3 first columns are fixed format.
- VIP SNP list: https://zenodo.org/api/files/efb7a8bc-efce-4391-8f65-59b974328c2e/VIP_GWAS_CLINVAR_ALL.txt - 3 first columns are fixed format.
- 
#### Recommended protocol for vcf preprocessing:
Assumed that you have a raw vcf file: `chr10_EAS.vcf.gz` in your current directory, recommended protocol to prepare vcf file is:
```sh
bcftools view chr10_EAS.vcf.gz -m2 -M2 -v snps -Q 0.9999999999:major -q 0.01:minor -e 'ALT="."' | bcftools +fill-tags | sed 's/chr//g' | bgzip > chr10_EAS_processed.vcf.gz
```


#### Compute needed information
Now we compute `LD` with `plink v1.9` and extract `MAF` with `bcftools`:
```sh
# compute LD
mkdir plink
plink --vcf chr10_EAS_processed.vcf.gz \
      --vcf-half-call 'haploid' \
      --make-bed  --const-fid --out ./plink/chr10_EAS_tem_file \
      --threads 1 \
      --memory 2000


plink --bfile ./plink/chr10_EAS_tem_file \
      --r --ld-window-r2 0.2 \
      --ld-window 10000 \
      --ld-window-kb 1000 \
      --out ./plink/chr10_EAS_ld_0.2 \
      --threads 8 \
      --memory 2000

mv ./plink/chr10_EAS_ld_0.2.ld ./
rm -r plink

# extract MAF
echo $'CHR\tPOS\tAF' > chr10_EAS_extracted_AF.txt
bcftools query -f '%CHROM\t%POS\t%AF\n' chr10_EAS_processed.vcf.gz >> chr10_EAS_extracted_AF.txt
```

### 4.1 Model construction

#### 4.1.1 build m3vcf reference for imputation
We need to build a reference directory for leave one out cross validation imputation with `create_imputation_ref.sh`. This step may take very long time because it will generate n reference m3vcf files with n-1 samples. n is number of sample in your `vcf.gz` file. You may download pre-built reference instead of generate it yourself.
```sh
create_imputation_ref.sh -v chr10_EAS_processed.vcf.gz -o chr10_EAS_hg38_high_cov -p 16
```
__NOTE:__  Pre-built imputation reference panel of populations are available for downloading at:
EAS: https://zenodo.org/api/files/efb7a8bc-efce-4391-8f65-59b974328c2e/chr10_EAS.tar.gz
EUR: https://zenodo.org/api/files/efb7a8bc-efce-4391-8f65-59b974328c2e/chr10_EUR.tar.gz
SAS: https://zenodo.org/api/files/efb7a8bc-efce-4391-8f65-59b974328c2e/chr10_SAS.tar.gz

Assumed that you have downloaded `chr10_EAS.tar.gz`, the upzip command is:
```sh
tar -xvzf chr10_EAS.tar.gz
```


#### 4.1.1 Create naive array and compute imputation accuracy for naive array.
We need to build a naive array to entablish relation between LD, MAF, and distance of SNPs with `build_naive_array.R`. The idea is to sampling n SNPs uniformlly based on their index after sorting by genomic position. Next, we impute with pre-built imputation reference m3vcf with `imputation_with_prebuilt_ref.sh` and compute imputaion accuracy with `compute_imputation_accuracy.R`.
```sh
# size=32970 is the size (number of tag SNP selected by obtain by TagIt with the same input vcf file - read main paper for further information).
build_naive_array.R vcf=chr10_EAS_processed.vcf.gz size=32970 out=chr10_EAS_naive.txt

imputation_with_prebuilt_ref.sh -t chr10_EAS_naive.txt -r chr10_EAS_hg38_high_cov -o naive_chr10_EAS -p 16

compute_imputation_accuracy.R imputation=naive_chr10_EAS out=naive_chr10_EAS.Rdata
```
### 4.1.2 Find best tagSNP
Finding best tagSNP (belong to chr10_EAS_naive.txt) by `LmTag find` command for all SNPs.
```sh
LmTag find --tag chr10_EAS_naive.txt --ld chr10_EAS_ld_0.2.ld -o chr10_EAS_find_snp_output.txt
```
### 4.1.3 Building model with computed inputs
Now we can build model with `buid_imputation_model.R`
```sh
buid_imputation_model.R imputation_Rdata=naive_chr10_EAS.Rdata find_snp=chr10_EAS_find_snp_output.txt out_Rdata=chr10_EAS_model.Rdata
```
OR you may download the model at: https://zenodo.org/api/files/efb7a8bc-efce-4391-8f65-59b974328c2e/chr10_EAS_model.Rdata

### 4.2 Tag SNP selection with LmTag
#### 4.2.1 Fitting model to generate input for LmTag
This step need ld file generated by `plink v1.9`, AF file generated by `bcftools`, and model build by `./buid_imputation_model.R`.
```sh
fit_imputation_model.R model_Rdata=chr10_EAS_model.Rdata af=chr10_EAS_extracted_AF.txt ld=chr10_EAS_ld_0.2.ld ld_cutoff=0.8 out_ld=chr10_EAS_ld_fitted_model.txt
```
#### 4.2.3 Run LmTag to slect tagSNPs
```sh
## testing with k = 200
LmTag tag --ld_model chr10_EAS_ld_fitted_model.txt \
      --eff chr10_EAS_CADD.txt \
      --vip VIP_GWAS_CLINVAR_ALL.txt \
      -k 200 \
      -o chr10_EAS_tagSNP.txt
```
Input argument:
| Name      | Description |
| ----------- | ----------- |
| --ld_model | fitted model ld file, generated by it_imputation_model.R |
| --eff | file provide infomation of effect score of input SNP |
| --exclude | list of SNPs that will be try to avoid to select by the algorithm |
| --vip | list of SNPs that will be try to prioritize to select by the algorithm |
| -k | k value of the beam search algorithm |
| -o | output file name |

## 5. Output

LmTag output file is __chr10_EAS_tagSNP.txt__ with following infomation

| Name      | Description |
| ----------- | ----------- |
| chr | chromosome of tag SNP|
| pos | position of tag SNP|
| id |  typically isrsID of tag SNP - but depends on input ld file|
| sum_score | sum of impuation score - used for weighting in tag SNP selection|
| degree | degree of tag SNP in the graph|
| effect_score | effect score of tag SNP|
| flag | source of tag SNP - it can be normal, vip (from vip list), excluded (from excluded list)|
| tagged_pos | positions of tagged SNPs |


## 6. Evaluation imputation accuracy performance

```sh
cat chr10_EAS_tagSNP.txt | grep -v "chr" | awk '//{printf "%s\t%s\n", $1, $2}'  > chr10_EAS_tagSNP_cleaned.txt


imputation_with_prebuilt_ref.sh -t chr10_EAS_tagSNP_cleaned.txt -r chr10_EAS_hg38_high_cov -o LmTag_chr10_EAS -p 16

compute_imputation_accuracy.R imputation=LmTag_chr10_EAS out=LmTag_chr10_EAS.Rdata
```
