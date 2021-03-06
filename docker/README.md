# Documents for running LmTag with docker

`LmTag` is a model based method to find tagSNP in SNP array desgin that maximizes imputation coverage and functional score of tag SNPs. Full details of the method is described in the method manuscript.
 
## 1. Pulling LmTag image and testing the installation

```sh
# pull lmtag docker image:
docker pull ndatth/lmtag:v0.2.0

# test your docker
docker run --rm ndatth/lmtag:v0.2.0 LmTag_pipeline.sh
```

Expected output is:

```text
Usage:
./LmTag_pipeline.sh 

       -v|--vcf : input vcf file to build (should be in bgzip formatted file, and preprocessed by our recommended protocol: bcftools view YOUR.vcf.gz -m2 -M2 -v snps -Q 0.9999999999:major -q 0.01:minor -e 'ALT=.' | bcftools +fill-tags | sed 's/chr//g' | bgzip > YOUR_processed.vcf.gz) 

       -m|--model : path to buit model generated by [model_pipeline.sh]. 

       -M|--model_setting : model setting can be [linear] that model the relationship of LD, MAF and distance as a simple linear model or [interaction] that consider interaction term of LD and MAF_taggedSNP. In our experiment, interaction model provides slightly better performance than simple linear model. However, we still prefer the simpler model as the performance difference is not significant.

       -s|--score : [optional] file provide infomation of effect score of input SNPs in the targeted chromosome. In practice, we use CADD scores for this option. 

       -V|--vip : [optional] list of SNPs that will be try to prioritize to select by the algorithm - genome-wide SNP list accepted. In practice, we used a combination of ClinVar and GWAS Catalog for this option. 

       -o|--out : [optional] output tagSNPs selected by LmTag algorithm, defaut: LmTag_tagSNP.txt. 

       -k|--beam : [optional] beamwith - defaut: 1. 


ERROR: no input vcf file !

```


## 2. Input requirement

### LmTag require a phased vcf file to perform tag SNP selection with following criteria:

- Only biallelic SNPS are considered.
- `MAF` are carefully check before runing (recommend to use `filltag` command by `bcftools` before input to LmTag pipeline) as LmTag extracts MAF directly from `INFO/AF` information in vcf file.
- Minimum `MAF` threshold are pre-determine in vcf filer, so you need to do filtering before providing the file to `LmTag` pipeline.
- `#CHROM` should be encoded without 'chr' character.

### Other possible input files for `LmTag`:

- Functional scores of SNPs candidates, typically extracted from the `CADD` databases [optional]
- List of VIP SNPs - they will be prioritized in tag SNP selection at highest levels; typically SNPs in `GWAS catalog` or `ClinVar` databases that you really want to select as tag SNPs [optional]

### For better understanding input format of LmTag, we prepare a real dataset and step by step copy-paste tutorial in the next sections.

## 3. Downloading the tutorial dataset

We provide in this tutorial based on chromosome 10, __East Asian population__ dataset:
- vcf file: https://zenodo.org/record/5807198/files/chr10_EAS.vcf.gz?download=1
- CADD score file: https://zenodo.org/record/5807198/files/chr10_EAS_CADD.txt?download=1
- VIP SNP list: https://zenodo.org/record/5807198/files/VIP_GWAS_CLINVAR_ALL.txt?download=1

We will create a folder name `lmtag_test` for the tutorial and using `wget` to download the data:

```sh
mkdir lmtag_test

wget https://zenodo.org/record/5807198/files/chr10_EAS.vcf.gz?download=1 -O lmtag_test/chr10_EAS.vcf.gz

wget https://zenodo.org/record/5807198/files/chr10_EAS_CADD.txt?download=1 -O lmtag_test/chr10_EAS_CADD.txt

wget https://zenodo.org/record/5807198/files/VIP_GWAS_CLINVAR_ALL.txt?download=1 -O lmtag_test/VIP_GWAS_CLINVAR_ALL.txt

# to save time in m3vcf reference building, the pre-buit m3vcf reference data can be downloaded and uncompressed by: 

wget https://zenodo.org/record/5807198/files/chr10_EAS.tar.gz?download=1 -O lmtag_test/chr10_EAS.tar.gz

# uncompress the file

cd lmtag_test
tar -xvzf chr10_EAS.tar.gz
cd ..
```


## 4. Step by step instruction to run LmTag

#### 4.1 Recommended protocol for vcf preprocessing [optional]:

Assumed that you have a raw vcf file: `chr10_EAS.vcf.gz` in `lmtag_test` directory, recommended protocol to prepare vcf file is, users can customize as their own needs.

```sh
# assumming you are running unix and $PWD/lmtag_test: is the path to lmtag_test directory that is mounted to docker containter as /lmtag_test:
docker run --rm -v $PWD/lmtag_test:/lmtag_test ndatth/lmtag:v0.2.0 process_vcf.sh \
      lmtag_test/chr10_EAS.vcf.gz \
      lmtag_test/chr10_EAS_processed.vcf.gz
```


#### 4.2. Building m3vcf reference for imputation

***Can be ignored in this tutorial as we already downloaded m3vcf reference***


We need to build a reference directory for leave one out cross validation imputation with `create_imputation_ref.sh`. This step may take long time because it will generate n reference m3vcf files with n-1 samples. n is number of sample in your `vcf.gz` file. If you already downloaded it, you can ignore this step.

```sh
# assumming you are running unix and $PWD/lmtag_test: is the path to lmtag_test directory that is mounted to docker containter as /lmtag_test:

docker run --rm -v $PWD/lmtag_test:/lmtag_test ndatth/lmtag:v0.2.0 create_imputation_ref.sh \
      -v lmtag_test/chr10_EAS_processed.vcf.gz \
      -o lmtag_test/chr10_EAS_hg38_high_cov \
      -p 32
```


#### 4.2. Estimating LmTag model parameter


We prepared a wrapper that comprises all steps needed to build model for LmTag. In practice, we recommend to build model for a typical chromosome (in human, we used chromosome 10, as it sizes is approximately average of all chromosome) and apply the same model to the remaining chromosomes.

```text
Usage:
./model_pipeline.sh 

       -v|--vcf : input vcf file to build (should be in bgzip formatted file, and preprocessed by our recommended protocol: bcftools view YOUR.vcf.gz -m2 -M2 -v snps -Q 0.9999999999:major -q 0.01:minor -e 'ALT=.' | bcftools +fill-tags | sed 's/chr//g' | bgzip > YOUR_processed.vcf.gz) 

       -r|--ref : path to output directory that is output of [create_imputation_ref.sh]. 

       -n|--size : number of tagSNP in simulated array. Theoratically, the number is need to be closed with your targeted number of tag SNP in the chromosome. We recommend this number is the maximum number of tag SNP selected by a standard greedy algorithm such as TagIt. 

       -o|--out : [optional] output model file name, defaut: LmTag_model.Rdata. 

       -p|--thread : [optional] number of thread - defaut: 8.
```

Running tutorial:

```sh
# assumming you are running unix and $PWD/lmtag_test: is the path to lmtag_test directory that is mounted to docker containter as /lmtag_test:

docker run --rm -v $PWD/lmtag_test:/lmtag_test ndatth/lmtag:v0.2.0 model_pipeline.sh \
      -v lmtag_test/chr10_EAS_processed.vcf.gz \
      -r lmtag_test/chr10_EAS_hg38_high_cov \
      -n 32970 \
      -p 32 \
      -o lmtag_test/LmTag_model_test.Rdata
```

#### 4.2. Tag SNP selection with LmTag


We prepared a wrapper that comprises all steps needed to select tag SNP with LmTag.

```text
Usage:
./LmTag_pipeline.sh 

       -v|--vcf : input vcf file to build (should be in bgzip formatted file, and preprocessed by our recommended protocol: bcftools view YOUR.vcf.gz -m2 -M2 -v snps -Q 0.9999999999:major -q 0.01:minor -e 'ALT=.' | bcftools +fill-tags | sed 's/chr//g' | bgzip > YOUR_processed.vcf.gz) 

       -m|--model : path to buit model generated by [model_pipeline.sh]. 

       -M|--model_setting : model setting can be [linear] that model the relationship of LD, MAF and distance as a simple linear model or [interaction] that consider interaction term of LD and MAF_taggedSNP. In our experiment, interaction model provides slightly better performance than simple linear model. However, we still prefer the simpler model as the performance difference is not significant.

       -s|--score : [optional] file provide infomation of effect score of input SNPs in the targeted chromosome. In practice, we use CADD scores for this option. 

       -V|--vip : [optional] list of SNPs that will be try to prioritize to select by the algorithm - genome-wide SNP list accepted. In practice, we used a combination of ClinVar and GWAS Catalog for this option. 

       -o|--out : [optional] output tagSNPs selected by LmTag algorithm, defaut: LmTag_tagSNP.txt. 

       -k|--beam : [optional] beamwith - defaut: 1. 

```

Running tutorial:

There are two models can be choose:

- The linear option [-M linear] models imputation accuracy of a linear funcion of: `Imputation_accuracy ~ LD + MAF_tagSNP + MAF_taggedSNP + distance`.


```sh
# assumming you are running unix and $PWD/lmtag_test: is the path to lmtag_test directory that is mounted to docker containter as /lmtag_test:

docker run --rm -v $PWD/lmtag_test:/lmtag_test ndatth/lmtag:v0.2.0 LmTag_pipeline.sh \
      -v lmtag_test/chr10_EAS_processed.vcf.gz \
      -m lmtag_test/LmTag_model_test.Rdata \
      -s lmtag_test/chr10_EAS_CADD.txt \
      -V lmtag_test/VIP_GWAS_CLINVAR_ALL.tx \
      -k 5 \
      -M linear \
      -o lmtag_test/test_LmTag_model_linear.txt
```

- The interaction option [-M interaction] models imputation accuracy of a linear funcion with interaction term of LD and tagged_MAF: `Imputation_accuracy ~ LD + MAF_tagSNP + MAF_taggedSNP + distance + LD:MAF_taggedSNP`.



```sh
# assumming you are running unix and $PWD/lmtag_test: is the path to lmtag_test directory that is mounted to docker containter as /lmtag_test:

docker run --rm -v $PWD/lmtag_test:/lmtag_test ndatth/lmtag:v0.2.0 LmTag_pipeline.sh \
      -v lmtag_test/chr10_EAS_processed.vcf.gz \
      -m lmtag_test/LmTag_model_test.Rdata \
      -s lmtag_test/chr10_EAS_CADD.txt \
      -V lmtag_test/VIP_GWAS_CLINVAR_ALL.txt \
      -k 5 \
      -M interaction \
      -o lmtag_test/test_LmTag_model_interaction.txt
```

In our experiments, the interaction model provides slightly better performance, but it is not significant. We still recommend users to use the linear model for better intepretbility.



## 5. Output

LmTag output files have following infomation

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


## 6. License

The Software is restricted to non-commercial research purposes.

## 7. Reference

LmTag: functional-enrichment and imputation-aware tag SNP selection for population-specific genotyping arrays
Dat Thanh Nguyen, Quan Nguyen, Duong Thuy Nguyen, Nam Sy Vo
bioRxiv 2022.01.28.478108; doi: https://doi.org/10.1101/2022.01.28.478108


