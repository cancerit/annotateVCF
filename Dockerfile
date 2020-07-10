FROM python:3.7-slim as builder

MAINTAINER  cgphelp@sanger.ac.uk

ENV ANNOTATEVCF_VER '1.0.0'

# install system tools

ENV CGP_OPT /opt/wtsi-cgp
RUN mkdir $CGP_OPT

RUN python -m venv $CGP_OPT/venv
# Make sure we use the virtualenv:
ENV PATH="$CGP_OPT/venv/bin:$PATH"

RUN pip3 install wheel
RUN pip3 install --upgrade setuptools
RUN pip3 --no-cache-dir install https://github.com/cancerit/annotateVCF/releases/download/${ANNOTATEVCF_VER}/annotateVcf-${ANNOTATEVCF_VER}-py3-none-any.whl
COPY . .

FROM python:3.7-slim

LABEL uk.ac.sanger.cgp="Cancer Genome Project, Wellcome Sanger Institute" \
      version="$ANNOTATEVCF_VER" \
      description="Tool to perform vcf file annotation"

# Make sure we use the virtualenv:
ENV CGP_OPT /opt/wtsi-cgp
RUN mkdir -p $CGP_OPT
COPY --from=builder $CGP_OPT $CGP_OPT
RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends \
    bcftools tabix

ENV PATH="$CGP_OPT/venv/bin:$PATH"

## USER CONFIGURATION 

RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

USER ubuntu

WORKDIR /home/ubuntu

CMD ["/bin/bash"]
