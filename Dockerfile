FROM python:3.7-slim as builder

MAINTAINER  cgphelp@sanger.ac.uk

ENV ANNOTATEVCF_VER '1.0.2'

# install system tools

ENV CGP_OPT /opt/wtsi-cgp/venv

RUN python -m venv $CGP_OPT
# Make sure we use the virtualenv:
ENV PATH $CGP_OPT/bin:$CGP_OPT/python-lib/bin:$PATH
ENV PYTHONPATH $CGP_OPT/python-lib/lib/python3.7/site-packages
RUN pip3 install --upgrade setuptools

RUN pip3 install --install-option="--prefix=$CGP_OPT/python-lib" https://github.com/cancerit/annotateVCF/archive/${ANNOTATEVCF_VER}.tar.gz

COPY . .

RUN python3 setup.py sdist
RUN pip3 install --install-option="--prefix=$CGP_OPT/python-lib" dist/$(ls -1 dist/)

FROM python:3.7-slim

LABEL uk.ac.sanger.cgp="Cancer Genome Project, Wellcome Sanger Institute" \
      version="$ANNOTATEVCF_VER" \
      description="Tool to perform vcf file annotation"

# Make sure we use the virtualenv:
ENV CGP_OPT /opt/wtsi-cgp/venv
RUN mkdir -p $CGP_OPT
COPY --from=builder $CGP_OPT $CGP_OPT
RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends \
    bcftools tabix
#set PATHS 
ENV PATH $CGP_OPT/bin:$CGP_OPT/python-lib/bin:$PATH
ENV PYTHONPATH $CGP_OPT/python-lib/lib/python3.7/site-packages
## USER CONFIGURATION 
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

USER ubuntu

WORKDIR /home/ubuntu

CMD ["/bin/bash"]
