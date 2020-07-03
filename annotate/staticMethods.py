import sys
import os
import io
import json
import tempfile
from subprocess import Popen, PIPE, STDOUT
import logging
import pkg_resources
import shutil

from contextlib import contextmanager

log = logging.getLogger(__name__)

FILTER_VARS = 'bcftools view -f PASS {} | bgzip -c >{} && tabix -f -p vcf {}'

MAP_GENES = 'bcftools annotate -a {} -i \' {} \' -h {} -c CHROM,FROM,TO,INFO/DRV {} >{}'

MAP_MUTATIONS = 'bcftools annotate -a {} -h {} -c CHROM,FROM,TO,INFO/DRV  {}|' \
                'bcftools annotate  -i \'INFO/DRV==INFO/VC  && INFO/DRV!="."\' | ' \
                'bgzip -c >{} && tabix -f -p vcf {}'

BGZIP_TABIX = 'cat {} | bgzip -c >{} && tabix -f -p vcf {}'

CONCAT_VCF = 'bcftools concat --allow-overlaps --rm-dups all {} {} {} | bcftools sort | ' \
             'bgzip -c >{} && tabix -f -p vcf {}'

UNHEADER_VCF = 'bcftools view -H {} >{}'


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
                return 'y'
            else:
                return None
        except IsADirectoryError:
            return None
        except IOError as ioe:
            sys.exit('Error in reading input file{}:{}'.format(ioe.args[0], infile))

    # ------------------------------Analysis methods---------------------------------
    @staticmethod
    def prepare_ref_data(json_file, ref_dir):
        try:
            if json_file is None:
                sys.exit('Driver reference json configuration file must be provided')
            with open(json_file, 'r') as cfgfile:
                cfg = json.load(cfgfile)
                path_dict = StaticMthods._format_dir_input(ref_dir)
                if path_dict[cfg['drv_genes']]:
                    drv_genes = path_dict[cfg['drv_genes']]
                else:
                    sys.exit('paramater not found {}'.format(cfg['drv_genes']))
                if path_dict[cfg['drv_genes_prev']]:
                    prev_genes = path_dict[cfg['drv_genes_prev']]
                else:
                    sys.exit('paramater not found {}'.format(cfg['drv_genes_prev']))
                if path_dict[cfg['drv_mut']]:
                    drv_mut = path_dict[cfg['drv_mut']]
                else:
                    sys.exit('paramater not found {}'.format(cfg['drv_mut']))
                if path_dict[cfg['header_info']]:
                    header_info = path_dict[cfg['header_info']]
                else:
                    sys.exit('paramater not found {}'.format(cfg['header_info']))
                if path_dict[cfg['genome_loc']]:
                    genome_loc = path_dict[cfg['genome_loc']]
                else:
                    sys.exit('paramater not found {}'.format(cfg['genome_loc']))
                if cfg['lof_consequences']:
                    lof_consequences = cfg['lof_consequences']
                else:
                    sys.exit('paramater not found {}'.format(cfg['lof_consequences']))

        except json.JSONDecodeError as jde:
            sys.exit('json error:{}'.format(jde.args[0]))
        except FileNotFoundError as fne:
            sys.exit('Can not find json file:{}'.format(fne.args[0]))
        return (drv_genes, prev_genes, drv_mut, header_info, genome_loc, lof_consequences)

    @staticmethod
    def get_drv_gene_list(drv_genes):
        with open(drv_genes) as f_drv:
            lof_gene_list = f_drv.read().splitlines()
        return lof_gene_list

    @staticmethod
    def get_drv_prev_gene_dict(drv_genes_prev):
        prev_gene_dict = {}
        with open(drv_genes_prev) as f_prev:
            for line in f_prev:
                (key, val) = line.split("\t")
                prev_gene_dict[key] = val.strip()
        return prev_gene_dict

    @staticmethod
    def format_consequences(lof_con):
        effect_type = []
        for effect in lof_con:
            effect_type.append("INFO/VC=\"" + effect + "\"")
        return ' || '.join(effect_type)

    @staticmethod
    def _format_dir_input(file_path):
        """
          creates a diretory object of key = file name and
          values = [file paths, name, extension, size]
        """
        path_dict = {}
        for dirpath, _, files in os.walk(file_path):
            for filename in files:
                fullpath = os.path.join(dirpath, filename)
                (_, name) = os.path.split(fullpath)
                path_dict[name] = fullpath
        return path_dict

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
    def get_filtered_vcf(vcf, outfile_name):
        """
        :param vcf: input vcf file
        :param filename: filename no extension
        :param ext: file extension
        :param outdir:
        :return: filtered var vcf outfile
        """
        global FILTER_VARS
        cmd = FILTER_VARS.format(vcf, outfile_name, outfile_name)
        StaticMthods.run_command(cmd)
        return outfile_name

    @staticmethod
    def map_drv_genes(filtered_vcf, info_vcf_prm, genome_loc, header_file, outfile_name):
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
        cmd = MAP_GENES.format(genome_loc, info_vcf_prm, header_file, filtered_vcf, outfile_name)
        StaticMthods.run_command(cmd)
        return outfile_name

    @staticmethod
    def filter_lof_genes(drv_gene_vcf, lof_gene_list, prev_gene_dict, out_filename):
        """
        :param drv_gene_vcf:
        :param drv_genes:
        :param prev_gene_dict:
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
                    elif prev_gene_dict.get(gene, None):
                        fh.write(line)
        fh.close()
        compressed_output = out_filename + '.gz'
        cmd = BGZIP_TABIX.format(out_filename, compressed_output, compressed_output)
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
    def unheader_vcf(header_vcf, unheader_vcf):
        # only used for testing ...
        global UNHEADER_VCF
        cmd = UNHEADER_VCF.format(header_vcf, unheader_vcf)
        StaticMthods.run_command(cmd)
        return unheader_vcf

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
            # logging.info("running command:{}".format(cmd))
            (out, error) = cmd_obj.communicate()
            exit_code = cmd_obj.returncode
            if (exit_code == 0):
                logging.info("Command run successfully:\n{}\n".format(cmd))
            else:
                logging.debug("Error: bcftools exited with non zero exit \
                              status, please check log file more details")
                logging.error("OUT:{}:Error:{}:Exit:{}".format(out, error, exit_code))
            return 0
        except OSError as oe:
            logging.error("Unable to run command:{} Error:{}".format(cmd, oe.args[0]))
            sys.exit("Unable to run command:{} Error:{}".format(cmd, oe.args[0]))

    @contextmanager
    def tempdir(mypath):
        path = tempfile.mkdtemp(dir=mypath)
        try:
            yield path
        finally:
            try:
                shutil.rmtree(path)
            except IOError:
                sys.stderr.writie('Failed to clean up temp dir {}'.format(path))
