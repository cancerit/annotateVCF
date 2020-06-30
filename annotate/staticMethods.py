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

FILTER_VARS='bcftools view -f PASS {} | bgzip -c >{} && tabix -f -p vcf {}'

MAP_GENES='bcftools annotate -a {} -i \'INFO/VC="stop_lost" || INFO/VC="start_lost" ' \
             '|| INFO/VC="ess_splice" || INFO/VC="nonsense" || INFO/VC="frameshift" \' -h {} -c CHROM,FROM,TO,INFO/DRV {} >{}'

MAP_MUTATIONS='bcftools annotate -a {} -h {} -c CHROM,FROM,TO,INFO/DRV  {}|' \
                 'bcftools annotate  -i \'INFO/DRV==INFO/VC  && INFO/DRV!="."\' | bgzip -c >{} && tabix -f -p vcf {}'

BGZIP_TABIX='cat {} | bgzip -c >{} && tabix -f -p vcf {}'

CONCAT_VCF='bcftools concat -a --rm-dups all {} {} {} | bgzip -c >{} && tabix -f -p vcf {}'

class StaticMthods(object):
    """ Static methosds for common tasks """

    def __init__(self):
        super().__init__()

    def input_checker(infile):
        """
          checks user input file and returns it's type
        """
        try:
            if os.path.exists(infile):
                log.info(("input file exists", infile))
                return 'y'
            else:
                return None
        except IsADirectoryError:
            return None
        except IOError as ioe:
            sys.exit('Error in reading input file{}:{}'.format(ioe.args[0], infile))

    # ------------------------------Analysis methods---------------------------------
    @staticmethod
    def get_file_metadata(full_file_name):
        """
          takes file path as input and gives its path and processed extension
          #  If there are two extensions adds second extensions as prefix
        """
        (_, name) = os.path.split(full_file_name)
        (name_no_ext, first_ext) = os.path.splitext(name)
        (name_no_ext2, second_ext) = os.path.splitext(name_no_ext)
        ext = second_ext + first_ext
        return name_no_ext2, ext

    @staticmethod
    def get_filtered_vars(vcf, outfile_name):
        """
        :param vcf: input vcf file
        :param filename: filename no extension
        :param ext: file extension
        :param outdir:
        :return: filtered var vcf outfile
        """
        global FILTER_VARS
        cmd=FILTER_VARS.format(vcf, outfile_name, outfile_name)
        StaticMthods.run_command(cmd)
        return outfile_name
    
    @staticmethod
    def map_drv_genes(filtered_vcf, genome_loc, header_file, outfile_name):
        """
        :param filtered_vcf:  read filtered vcf file
        :param filename: oufile name wiyhout extension
        :param ext: file extension
        :param genome_loc: genome location
        :param header_file: vcf custom header file
        :param outdir: output directory
        :return:
        """
        global MAP_GENES
        cmd=MAP_GENES.format(genome_loc, header_file, filtered_vcf, outfile_name)
        StaticMthods.run_command(cmd)
        return outfile_name

    @staticmethod
    def filter_lof_genes(drv_gene_vcf, lof_gene_list, out_filename):
        """
        :param drv_gene_vcf:
        :param drv_genes:
        :param out_filename:
        :return:
        """
        fh = open(out_filename, "w")
        with open(drv_gene_vcf) as gene_f:
            for line in gene_f:
                # write header lines
                if line.startswith('#'):
                    fh.write(line)
                else:
                    gene = (line.split(';VD=')[1]).split('|')[0]
                    # write matching LoF genes....
                    if gene in lof_gene_list:
                        fh.write(line)
        fh.close()
        compressed_output= out_filename + '.gz'
        cmd = BGZIP_TABIX.format(out_filename, compressed_output, compressed_output )
        StaticMthods.run_command(cmd)
        return compressed_output

    @staticmethod
    def map_drv_mutations(filtered_vcf, drv_muts, header_file, outfile_name):
        """

        :param filtered_vcf:
        :param drv_muts:
        :param header_file:
        :param outfile_name:
        :return:
        """
        global MAP_MUTATIONS
        cmd = MAP_MUTATIONS.format(drv_muts, header_file, filtered_vcf, outfile_name, outfile_name)
        StaticMthods.run_command(cmd)
        return outfile_name

    @staticmethod
    def concat_results(drv_muts, lof_vcf, filtered_vcf, outfile_name):
        """
        
        :param drv_muts: 
        :param lof_vcf: 
        :param filtered_vcf: 
        :param outfile_name: 
        :return: 
        """
        global CONCAT_VCF
        cmd = CONCAT_VCF.format(drv_muts, lof_vcf, filtered_vcf, outfile_name, outfile_name)
        StaticMthods.run_command(cmd)
        return

    @staticmethod
    def run_command(cmd):
        """
        : runs external command
        :param cmd:
        :return: command output
        """
        """ runs command in a shell, returns stdout and exit code"""
        if not len(cmd):
            raise ValueError("Must supply at least one argument")
        try:
            # To capture standard error in the result, use stderr=subprocess.STDOUT:
            cmd_obj = Popen(cmd, stdin=None, stdout=PIPE, stderr=PIPE,
                            shell=True, universal_newlines=True, bufsize=-1,
                            close_fds=True, executable='/bin/bash')
            logging.info("running command:{}".format(cmd))
            (out, error) = cmd_obj.communicate()
            exit_code = cmd_obj.returncode
            if (exit_code == 0):
                logging.info("bcftools run successfully")
            else:
                logging.debug("Error: bcftools exited with non zero exit status, please check log file more details")
                logging.error("OUT:{}:Error:{}:Exit:{}".format(out, error, exit_code))
            return 0
        except OSError as oe:
            logging.error("Unable to run command:{} Error:{}".format(cmd, oe.args[0]))
            sys.exit("Unable to run command:{} Error:{}".format(cmd, oe.args[0]))


