FROM  ubuntu:20.04 as builder
USER root

MAINTAINER  cgphelp@sanger.ac.uk

ENV ANNOTATEVCF_VER '1.0.0'

# install system tools
RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends \
locales \
g++ \
make \
gcc \
pkg-config \
python3 python3-dev python3-pip python3-setuptools python3-wheel \
zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev \
bcftools tabix git

ENV CGP_OPT /opt/wtsi-cgp
RUN mkdir $CGP_OPT
ENV PYTHONPATH $CGP_OPT/python-lib/lib/python3.6/site-packages

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

COPY requirements.txt $CGP_OPT/requirements.txt
RUN pip3 --no-cache-dir install -r $CGP_OPT/requirements.txt
# install annotatevcf
RUN pip3 --no-cache-dir install https://github.com/cancerit/annotateVCF/releases/download/${ANNOTATEVCF_VER}/annotateVcf-${ANNOTATEVCF_VER}-py3-none-any.whl

COPY ..

FROM ubuntu:20.04

LABEL uk.ac.sanger.cgp="Cancer Genome Project, Wellcome Sanger Institute" \
      version="1.0.0" \
      description="Tool to perform vcf file annotation"

### security upgrades and cleanup
RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends \
apt-transport-https \
locales \
ca-certificates \
time \
unattended-upgrades \
python3 \
zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev && \
unattended-upgrade -d -v && \
apt-get remove -yq unattended-upgrades && \
apt-get autoremove -yq

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

ENV CGP_OPT /opt/wtsi-cgp
ENV PATH $CGP_OPT/bin:$CGP_OPT/python-lib/bin:$PATH
ENV PYTHONPATH $CGP_OPT/python-lib/lib/python3.6/site-packages
ENV LD_LIBRARY_PATH $OPT/lib
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

RUN mkdir -p $CGP_OPT
COPY --from=builder $CGP_OPT $CGP_OPT

## USER CONFIGURATION 
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

USER ubuntu

WORKDIR /home/ubuntu

CMD ["/bin/bash"]
