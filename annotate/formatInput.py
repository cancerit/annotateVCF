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

        for infile in ([self.vcf_file]):
            input_type.append(SM.input_checker(infile))
        return input_type

    def run_analysis(self):
        """
          method to run the analysis
        """
        outdir = self.outdir
        keepTmp = self.keepTmp
        vcf_file = self.vcf_file
        drv_data = self.drv_data
        drv_json = self.drv_json
        outputPath = outdir + '/out_anotatevcf/'
        os.makedirs(outputPath, exist_ok=True)
        # check input files
        (vcf_status) = self.check_input()

        if vcf_status[0] == 'y':
            log.info("input file checks PASSED, running analysis, .....")
            # seprate extension and file name
            (file_name, ext) = SM.get_file_metadata(vcf_file)
            # prepare reference input data
            (drv_genes, drv_genes_prev, drv_muts, header_file,
             genome_loc, lof_consequences) = SM.prepare_ref_data(drv_json, drv_data)
            with SM.tempdir(outputPath) as base_dir:
                logging.info(base_dir)
                outfile_name = base_dir + '/' + file_name
                final_output = outputPath + file_name
                if keepTmp:
                    outfile_name = outputPath + file_name
                # get lof consequences
                info_vc_prm = SM.format_consequences(lof_consequences)
                # create list of lof and ambiguous genes...
                lof_gene_list = SM.get_drv_gene_list(drv_genes)
                # create dict of previous gene symbol
                prev_gene_dict = SM.get_drv_prev_gene_dict(drv_genes_prev)
                # get passed  records from VCF file
                filtered_vcf = SM.get_filtered_vcf(vcf_file, outfile_name + '.pass' + ext)
                # get records with LoF effect ....
                drv_gene_vcf = SM.map_drv_genes(filtered_vcf, info_vc_prm, genome_loc, header_file,
                                                outfile_name + '.genes.vcf')
                # filter LoF records with Driver genes [ LoF and ambiguious ]
                lof_drv_gene_vcf = SM.filter_lof_genes(drv_gene_vcf, lof_gene_list, prev_gene_dict,
                                                       outfile_name + '.genes.lof.vcf')
                # map driver mutatations...
                drv_mut_vcf = SM.map_drv_mutations(filtered_vcf, drv_muts,
                                                   header_file, outfile_name + '.muts' + ext)
                # combine results...
                SM.concat_results(drv_mut_vcf, lof_drv_gene_vcf,
                                  filtered_vcf, final_output + '.drv' + ext)
                logging.info("analysis completed successfully")
                logging.info("temporary files will be removed if flag -tmp is not set")

        else:
            sys.exit('Input data is not in required format OR input file does  \
                not exists, see inputFormat in README file')
