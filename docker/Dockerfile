FROM ubuntu:18.04

ENV TZ=America/Los_Angeles
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt update -y
RUN apt install -yq r-base
RUN apt install -yq tabix
RUN apt install -yq vcftools
RUN apt install -yq bcftools

RUN R -e 'install.packages("data.table")'
RUN R -e 'install.packages("R.utils")'

ADD ./bin_biotools/ /bin_biotools
COPY /bin_biotools/* /bin
RUN rm -rf /bin_biotools/



