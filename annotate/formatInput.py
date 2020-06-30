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
        
        for infile in (self.vcf_file, self.drv_muts, self.genome_loc, self.lof_consequences , self.header_file ) :
          input_type.append(SM.input_checker(infile))
        return input_type

    def run_analysis(self):
        """
          method to run the analysis
        """
        vcf = self.vcf_file
        drv_muts = self.drv_muts
        drv_genes = self.drv_genes
        genome_loc = self.genome_loc
        lof_consequences = self.lof_consequences
        header_file = self.header_file			
        outdir = self.outdir

        global MAGECK_CMD
        if outdir:
            outputPath= outdir + '/anotatevcfOut/'
            os.makedirs(outputPath, exist_ok=True)
        # check input files
        (vcf_status) = self.check_input()
     
        if vcf_status[0] ==  'y':
            log.info("input file checks PASSED, running analysis, .....")
            (file_name, ext)=SM.get_file_metadata(vcf)
            with SM.tempdir(outputPath) as base_dir:
                logging.info(base_dir)
                outfile_name=base_dir+'/'+file_name
                final_output=outputPath+file_name
                # create list of lof genes...
                with open(drv_genes) as f:
                    lof_gene_list = f.read().splitlines()
                info_vcf_prm = SM.format_consequences(lof_consequences)
                filtered_vcf = SM.get_filtered_vars(vcf, outfile_name+'.pass'+ext)
                drv_gene_vcf = SM.map_drv_genes(filtered_vcf, info_vcf_prm, genome_loc, header_file, outfile_name+'.genes.vcf')
                lof_drv_gene_vcf = SM.filter_lof_genes(drv_gene_vcf, lof_gene_list, outfile_name+'.genes.lof.vcf' )
                drv_mut_vcf = SM.map_drv_mutations(filtered_vcf, drv_muts, header_file, outfile_name+'.muts'+ext)
                SM.concat_results(drv_mut_vcf, lof_drv_gene_vcf, filtered_vcf, final_output+'.drv'+ext)
                logging.info("analysis completed successfully")
                logging.info(base_dir)
                
        else:
            sys.exit('Input data is not in required format OR input file does not exists, see inputFormat in README file')

