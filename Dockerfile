FROM  ubuntu:20.04 as builder
USER root

MAINTAINER  cgphelp@sanger.ac.uk

ENV ANNOTATEVCF_VER '1.2.4'

# install system tools
RUN apt-get -yq update \
&& apt-get -yq install --no-install-recommends software-properties-common \
&& add-apt-repository ppa:deadsnakes/ppa \
&& locales \
&& g++ \
&& make \
&& gcc \
&& pkg-config \
&& python3.7 python3.7-dev python3-pip python3-setuptools \
&& zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev \
&& curl -sSL --retry 10 https://bootstrap.pypa.io/get-pip.py > get-pip.py \
&& python3.7 get-pip.py \
&& git

ENV CGP_OPT /opt/wtsi-cgp
RUN mkdir $CGP_OPT
ENV PATH $CGP_OPT/bin:$CGP_OPT/python-lib/bin:$PATH
ENV PYTHONPATH $CGP_OPT/python-lib/lib/python3.7/site-packages

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

# install annotatevcf
RUN pip3.7 install --upgrade setuptools

RUN pip3.7 install --install-option="--prefix=$CGP_OPT/python-lib" https://github.com/cancerit/annotateVCF/archive/${ANNOTATEVCF_VER}.tar.gz

COPY . .

RUN python3.7 setup.py sdist
RUN pip3.7 install --install-option="--prefix=$CGP_OPT/python-lib" dist/$(ls -1 dist/)

FROM ubuntu:20.04

LABEL uk.ac.sanger.cgp="Cancer Genome Project, Wellcome Sanger Institute" \
      version="1.2.4" \
      description="Tool to perform vcf file annotation"

### security upgrades and cleanup
RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends \
apt-transport-https \
locales \
ca-certificates \
time \
tabix \
bcftools \
unattended-upgrades \
python3.7 \
python3-setuptools \
python3-pip \
zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev && \
unattended-upgrade -d -v
RUN apt-get autoremove -yq

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8
ENV CGP_OPT /opt/wtsi-cgp
ENV PATH $CGP_OPT/bin:$CGP_OPT/python-lib/bin:$PATH
ENV PYTHONPATH $CGP_OPT/python-lib/lib/python3.7/site-packages
RUN pip3 install --upgrade setuptools
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
