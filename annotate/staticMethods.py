import sys
import os
import tarfile
import json
import io
import tempfile
from subprocess import Popen, PIPE, STDOUT
import logging
import pkg_resources

log = logging.getLogger(__name__)

BCFTOOLS_ALL_PASS_CMD='bcftools view -f PASS {} | bgzip -c >{} && tabix -p vcf {}'

BCFTOOLS_GENE_CMD='bcftools view -f PASS {} | bcftools annotate -a {} -i \'INFO/VC="stop_lost" || INFO/VC="start_lost" ' \
             '|| INFO/VC="ess_splice" || INFO/VC="nonsense" || INFO/VC="frameshift" \' -h {} -c CHROM,FROM,TO,INFO/DRV >{}'

BCFTOOLS_MUT_CMD='bcftools view -f PASS {} | ' \
                 'bcftools annotate -a {} -h {} -c CHROM,FROM,TO,INFO/DRV |' \
                 'bcftools annotate  -i \'INFO/DRV==INFO/VC  && INFO/DRV!="."\' | bgzip -c >{} && tabix -p vcf {}'

BGZIP_TABIX_CMD='cat {} | bgzip -c > {} && tabix -p vcf {}'

BCFTOOLS_CONCAT_VCF='bcftools concat -a --rm-dups all {} {} {} | bgzip -c >{} && tabix -p vcf {}'

class StaticMthods(object):
    """ Static methosds for common tasks """

    def __init__(self):
        super().__init__()

    def input_checker(infile):
        """
          checks user input file and returns it's type
        """
        try:
            if tarfile.is_tarfile(infile):
                log.info(("input is an archive:", infile))
                return 'tar'
            else:
                log.info(("input is a file:", infile))
                return 'file'
        except IsADirectoryError:
            return 'dir'
        except IOError as ioe:
            sys.exit('Error in reading input file:{}'.format(ioe.args[0]))

    # ------------------------------Analysis methods---------------------------------
		@staticmethod
    def combine_count_n_library(countfile, libfile, outdir='./'):
				print("hello")
		return None
