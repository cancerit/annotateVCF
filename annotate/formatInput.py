import logging
import json
import os
import sys
from annotate.abstractAnnotate import AbstractAnnotate
from annotate.staticMethods import StaticMthods as SM
log = logging.getLogger(__name__)

'''
  This code run's bcftools to annoate driver gene and varinat sites
'''


class AnnotateVcf(AbstractAnnotate):
    """
        Main class , loads user defined parameters and files
    """

    def check_input(self):
        """
           check input type and presence of user supplied
           input files
        """
        super().check_input()
        input_type = []
        for infile in (self.vcf_file):
            input_type.append(SM.input_checker(infile))
        return input_type

    def run_analysis(self):
        """
          method to run the analysis
        """
        controls = self.ncontrols
        min_read_count = self.minreads
        min_target_genes = self.mingenes
        ignored_genes = self.ignored_genes
        cpus = self.num_processors
        outdir = self.outdir
        gene_sig_dir = self.gene_sig_dir
        result_cfg = self.results_cfg
        iter = self.numiter

        global MAGECK_CMD
        if outdir:
            os.makedirs(outdir + '/mageckOut', exist_ok=True)
            os.makedirs(outdir + '/bagelOut', exist_ok=True)

        # check input files
        (input1, input2) = self.check_input()

        if input1 and input2:
            log.info("Running analysis, input file checks DONE.....")
            cldf = SM.combine_count_n_library(self.countfile, self.libfile, outdir=outdir)
        else:
            sys.exit('Input data is not in required format, see inputFormat in README file')

