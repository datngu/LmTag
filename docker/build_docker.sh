# this script build docker image of LmTag manually

# run container
docker run -it --rm -h datn --name LmTag ubuntu:18.04

apt update -y
apt install -yq vcftools
apt install -yq bcftools
apt install -yq r-base
apt install -yq tabix

R -e 'install.packages("data.table")'
R -e 'install.packages("R.utils")'



# crtl+p, crtl+q
docker cp bin_biotools/* LmTag:/bin

docker attach LmTag
cp bin_biotools/* bin/
exit

docker commit ubuntu_18.04 lmtag:v.0.1.1


## BUILD DOCKER AUTOMATICLY

docker build -t lmtag:v0.2.0 -f Dockerfile .


## testing
#docker run -it --rm --name LmTag_test -h datn lmtag:v0.2.0

# v0.2.0
docker tag lmtag:v0.2.0 ndatth/lmtag:v0.2.0
docker push ndatth/lmtag:v0.2.0