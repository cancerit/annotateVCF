FROM  ubuntu:16.04
MAINTAINER  cgphelp@sanger.ac.uk

LABEL uk.ac.sanger.cgp="Cancer Genome Project, Wellcome Sanger Institute" \
      version="1.0.0" \
      description="Tool to perform vcf file annotation"



USER root

ENV VERSION 1.0.0
ENV OPT /opt/wtsi-cgp
ENV PATH $OPT/bin:$PATH
ENV LD_LIBRARY_PATH $OPT/lib

# install system tools
RUN apt-get update && \
  apt-get install -yq --no-install-recommends lsb-release && \
  apt-get update && \
  apt-get install bcftools \
  apt-get install tabix \
  apt-get install -qy --no-install-recommends \
    apt-transport-https \
    locales \
    libcairo2-dev \
    python3  \
    python3-dev \
    python3-setuptools \
    python3-pip \
    python3-wheel \
    unattended-upgrades && \
    unattended-upgrade -d -v && \
    apt-get remove -yq unattended-upgrades && \
    apt-get autoremove -yq

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

COPY requirements.txt $OPT/requirements.txt
RUN pip3 --no-cache-dir install -r $OPT/requirements.txt
# install annotatevcf
RUN pip3 --no-cache-dir install https://github.com/cancerit/pyCRISPRcleanR/releases/download/${VERSION}/annotateVcf-${VERSION}-py3-none-any.whl

### security upgrades and cleanup
RUN apt-get -yq update && \
    apt-get -yq install unattended-upgrades && \
    unattended-upgrades

RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

USER ubuntu

WORKDIR /home/ubuntu

CMD ["/bin/bash"]
