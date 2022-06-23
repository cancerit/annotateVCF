#
# Copyright (c) 2021
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
#
# This file is part of annotatevcf.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
# 2009, 2011-2012’ should be interpreted as being identical to a statement that
# reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
# statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
# identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012’.


import logging
import sys
import os
import re
import io
import copy
import vcf

from vcf import utils
from collections import defaultdict
from subprocess import Popen, PIPE, STDOUT, TimeoutExpired
import pkg_resources

'''
  This code runs bcftools and pyvcf to annotates user provided driver gene and variant sites
'''


class VcfAnnotator:
    """
        Annotator class, does all the processing and output generation
    """

    def __init__(self, f, basedir):
        """
        initialises required parameters..
        :param f: formatted input object
        :param basedir: tmpdir
        """
        self.outdir = basedir
        self.status_dict = f.format(['input_status'])
        self.input_data = f.format(['vcf_file'])
        # set input vcf parameters ...
        self._set_input_vcf(self.input_data)
        self.drv_header = f.header_line
        self.merge_vcf_dict = {}
        self.annotated_vcf_list = []
        self.keepTmp = f.keepTmp
        # run analysis
        self._runAnalysis(f)

    def _runAnalysis(self, f):
        """
        runs selected  the analysis steps based on user input
        :param f: formatted input object
        :return:
        """
        status = self.status_dict['input_status']
        vcf_filters = f.format(['vcf_filters'])
        vcf_filter_params = vcf_filters.get('vcf_filters', None)
        self.info_filters = vcf_filter_params['INFO']
        self.format_filters = vcf_filter_params['FORMAT']
        self.filter_filters = vcf_filter_params['FILTER']
        self.flag_germline = vcf_filter_params['INFO_FLAG_GERMLINE']
        drv_type = f.format(['drv_type'])
        self.drv_type_dict = drv_type.get('drv_type', None)

        run_status = False
        a_type = ['normal_panel', 'mutations', 'lof_genes', 'cancer_predisposition']

        for analysis in a_type:
            if status[analysis] and analysis == 'normal_panel':
                logging.info(f"Tagging germline variants with INFO tag:{self.flag_germline}")
                self.tag_germline_vars(f.np_vcf)
                run_status = True
            if status[analysis] and analysis == 'mutations':
                logging.info("Annotating driver mutations with INFO field:DRV=<consequence(s)>")
                self.annot_drv_muts(f.muts_file)
                run_status = True
            if status[analysis] and analysis == 'lof_genes':
                logging.info("Annotating LoF genes with INFO field:DRV=LoF")
                self.annotate_lof_genes(f.genes_file)
                run_status = True
            if status[analysis] and analysis == 'cancer_predisposition':
                logging.info("Annotating germline variants with INFO field:CPV=<consequence(s)")
                self.annotate_cpv(f.cpv_file)
                run_status = True
        if run_status:
            logging.info("concatenating results")
            self.concat_results()
            logging.info("Analysis completed successfully")
        else:
            logging.info("Input files not accessible analysis aborted")

    def _set_input_vcf(self, input_data):
        """
        sets input vcf data
        :param input_data:
        :return:
        """
        self.vcf_path = input_data['vcf_file']['path']
        self.vcf_ext = input_data['vcf_file']['ext']
        self.vcf_name = input_data['vcf_file']['name']
        self.outfile_name = self.outdir + '/' + self.vcf_name + '{}'

    def tag_germline_vars(self, np_vcf):
        """
        use normal panel to tag germline variants and create filtered vcf file to use in
        subsequent driver annotation steps ...
        sets filtered vcf as new user input parameter for downstream analysis
        add tagged and filtered vcf files to concat in final step
        :param np_vcf:
        :return:
        """
        tagged_vcf = self.outfile_name.format('_np.vcf.gz')
        filtered_vcf = self.outfile_name.format('_np_filtered.vcf.gz')

        cmd = f"bcftools annotate -a {np_vcf} -i '{self.filter_filters}'  -m '{self.flag_germline}'" \
              f" {self.vcf_path} | bgzip -c >{tagged_vcf} && tabix -p vcf {tagged_vcf}"
        _run_command(cmd)

        cmd = f"bcftools  view -i '{self.flag_germline}=0' {tagged_vcf} | " \
              f"bgzip -c >{filtered_vcf} && tabix -p vcf {filtered_vcf}"

        _run_command(cmd)

        # update input file for driver annotation step ...
        input_data = copy.deepcopy(self.input_data)
        input_data['vcf_file']['path'] = filtered_vcf
        self._set_input_vcf(input_data)

        # add vcf path to concat in the  order of priority...
        self.merge_vcf_dict['d'] = filtered_vcf
        # contains all record taggedwith NPGL
        self.merge_vcf_dict['e'] = tagged_vcf

    def annot_drv_muts(self, muts_file):
        """
        annotate vcf using driver mutation set
        add annotated vcf file to concat in final step
        :param muts_file:
        :return:
        """
        muts_outfile = self.outfile_name.format('_muts.vcf.gz')
        combined_filter = _combine_filters([self.filter_filters, self.format_filters])
        cmd = f"bcftools annotate -i '{combined_filter}'" \
              f" --merge-logic DRV:unique" \
              f" -a {muts_file} -h {self.drv_header} " \
              f"-c CHROM,FROM,TO,INFO/DRV {self.vcf_path} " \
              f" | bcftools annotate  -i 'INFO/DRV!=\".\" && INFO/DRV[*]==INFO/VC' " \
              f" | bgzip -c >{muts_outfile} && tabix -f -p vcf {muts_outfile}"
        _run_command(cmd)
        self.annotated_vcf_list.append(muts_outfile)

    def annotate_lof_genes(self, genes_file):
        """
        annotate vcf using known LoF genes
        :param genes_file: lof gene names file
        :return:
        """
        # create dummy genome locationo file to annoate LoF genes...
        get_gene = re.compile(r'\bVD=([-\w]+)')
        genome_loc_file = self.outdir + '/genome.tab.gz'
        create_dummy_genome(self.vcf_path, genome_loc_file)
        genes_outfile = self.outfile_name.format('_genes.vcf')
        lof_outfile = self.outfile_name.format('_genes_lof.vcf')
        lof_gene_list = get_drv_gene_list(genes_file)
        combined_filter = _combine_filters([self.filter_filters, self.format_filters, self.info_filters])
        cmd = f"bcftools annotate -a {genome_loc_file} -i '{combined_filter}' " \
              f"-h {self.drv_header} -c CHROM,FROM,TO,INFO/DRV {self.vcf_path} >{genes_outfile}"

        _run_command(cmd)
        with open(lof_outfile, "w") as lof_fh, open(genes_outfile, 'r') as gene_f:
            for line in gene_f:
                # write header lines
                if line.startswith('#'):
                    lof_fh.write(line)
                else:
                    gene = get_gene.search(line)[1]
                    # gene = _get_gene('VD', 0))
                    # write matching LoF genes....
                    if gene.upper() in lof_gene_list:
                        lof_fh.write(line)
        #compress and store vcf
        lof_outfile_gz = compress_vcf(lof_outfile)
        self.annotated_vcf_list.append(lof_outfile_gz)

    def annotate_cpv(self, cpv_file):
        """
        annotate vcf using cancer predisposition variant set
        add annotated vcf file to concat in final step
        :param cpv_file:
        :return:
        """
        cpv_outfile = self.outfile_name.format('_cpv.vcf.gz')
        combined_filter = _combine_filters([self.filter_filters, self.format_filters])
        cmd = f"bcftools annotate -i '{combined_filter}'" \
              f" --merge-logic CPV:unique" \
              f" -a {cpv_file} -h {self.drv_header} " \
              f"-c CHROM,FROM,TO,INFO/CPV {self.vcf_path} " \
              f" | bcftools annotate  -i 'INFO/CPV!=\".\" && INFO/CPV[*]==INFO/VC'  " \
              f" | bgzip -c >{cpv_outfile} && tabix -f -p vcf {cpv_outfile}"
        _run_command(cmd)
        self.annotated_vcf_list.append(cpv_outfile)

    def _get_vcf_readers(self):
        """
         create vcf file readers for each annotated vcf file
        :return: list of vcf file readers
        """
        vcf_readers = []
        for my_vcf in self.annotated_vcf_list:
            vcf_readers.append(vcf.Reader(filename=my_vcf))
        return  vcf_readers

    def walk_annotated_vcf(self, vcf_readers):
        """
        walk through each line of vcf files, matching records will be walked together
        However to circumvent occasional anomalous behaviour of pyVCF I have created separate 'vcf_key' to store matched records and their driver type
        records written in a collection for each 'vcf_key' as ke val pair of driver_type and vcf record
        e.g., 'chr17_7675139_C_[A]': {'somatic': <vcf.model._Record object at 0x7f82ca93ee50>, 'germline': <vcf.model._Record object at 0x7f82ca93ee80>}
        :return: collection as above
        """
        write_record = defaultdict(dict)

        for record in utils.walk_together(*vcf_readers):
            filtered_record = list(filter(None, record))
            for vcf_line in filtered_record:
                vcf_key = f"{vcf_line.CHROM}_{vcf_line.POS}_{vcf_line.REF}_{vcf_line.ALT}"
                (drv_type, drv_val) = self.get_driver_type(vcf_line)
                write_record[vcf_key][drv_val] = vcf_line
        return write_record

    def get_driver_type(self, vcf_line):
        """
        get driver type for given vcf line
        :return:
        """
        for drv_type in self.drv_type_dict.keys():
            if vcf_line.INFO.get(drv_type, None):
                drv_val = self.drv_type_dict.get(drv_type, 'NA')
                return drv_type, drv_val

    def set_driver_type(self, vcf_line):
        """
        set driver type for  vcf line to False
        :return:
        """
        for drv_type in self.drv_type_dict.keys():
          if vcf_line.INFO.get(drv_type, None):
            vcf_line.INFO[drv_type] = False


    def print_record(self, write_record, vcf_writer):
        """
        print driver record to vcf header template file generated from one of the annotated vcf file
        :return:
        """
        for records in write_record.values():
            drv_type = set()
            edit_line = None
            for my_drv, edit_line in records.items():
                drv_type.add(my_drv)
            self.set_driver_type(edit_line)
            edit_line.INFO['DRV'] = sorted(drv_type)
            vcf_writer.write_record(edit_line)
        vcf_writer.close()

    def concat_results(self):
        """
        process and concatenate the annotated vcf files...
        :return:
        """
        annotated_vcfs = self.annotated_vcf_list
        concat_drv_out = self.outfile_name.format('_drv.vcf.gz')
        # combine annotated vcf files...
        vcf_readers=self._get_vcf_readers()
        write_record=self.walk_annotated_vcf(vcf_readers)
        tmp_header = f"{self.outdir}/tmp_header_file"
        (vcf_writer, combined_vcf) =_get_vcf_writer(annotated_vcfs[0], tmp_header)
        self.print_record(write_record, vcf_writer)
        combined_vcf_gz = compress_vcf(combined_vcf)

        self.merge_vcf_dict['a'] = combined_vcf_gz
        self.merge_vcf_dict['f'] = self.input_data['vcf_file']['path']

        vcf_files = ([self.merge_vcf_dict[myvcf] for myvcf in sorted(self.merge_vcf_dict.keys())])

        cmd = f"bcftools concat --allow-overlaps --rm-dups all {' '.join(vcf_files)} | bcftools sort | " \
              f"bgzip -c >{concat_drv_out} && tabix -f -p vcf {concat_drv_out}"
        _run_command(cmd)
        _run_command('cp -p ' + concat_drv_out + '  ' + concat_drv_out + '.tbi ' + self.outdir + '/..')
        if self.keepTmp:
            _run_command('cp -rp ' + self.outdir + ' ' + self.outdir + '_' + self.vcf_name)


# generic methods ....
def _combine_filters(filter_array):
    """
    :param filter_array: filtring paramters
    :return: return formatted filtering parameters if present otherwise () equivalent to no filter...
    """
    if any(filter_array):
        return f"({') && ('.join(filter(None, filter_array))})"
    else:
        return "()"


def _get_gene(line, gene_field, field_loc):
    # Not used ... kept for future implementation of different annotation fields....
    # ANN=T|missense_variant|MODERATE|AGAP005273|AGAP005273| [ e.g. 'ANN', 3]
    # VD=TP5-TEST1-TEST2|CCDS11118.1|r.276_277insa|c.86_87insA|p.N29fs*14| [ e.g. 'VD', 0 ]
    info_list = line.split("\t")[7].split(';')
    info_dict = dict(f.split('=') for f in info_list if '=' in f)
    gene = info_dict[gene_field].split('|')[field_loc]
    return gene.upper()


def get_drv_gene_list(drv_genes):
    with open(drv_genes) as f_drv:
        lof_gene_list = [gene.upper() for gene in f_drv.read().splitlines()]
    return lof_gene_list


def compress_vcf(vcf):
    """
    :param vcf:
    :return:
    """
    outfile = vcf + '.gz'
    cmd = f"bgzip -c <{vcf} >{outfile} && tabix -f -p vcf {outfile}"
    _run_command(cmd)
    return outfile


def create_dummy_genome(input_vcf, genome_loc):
    """
    :param input_vcf:
    :param genome_loc:
    :return:
    """
    cmd = f"tabix -l {input_vcf} | xargs -I vcf_chr printf 'vcf_chr\t1\t400000000\tLoF\n' |" \
          f"bgzip -c >{genome_loc} && tabix -s1 -b2 -e3 {genome_loc}"
    _run_command(cmd)


def unheader_vcf(header_vcf, unheader_vcf):
    """
    This is used for testing purpose only...
    :param header_vcf:
    :param unheader_vcf:
    :return:
    """
    cmd = f"bcftools view -H {header_vcf} >{unheader_vcf}"
    _run_command(cmd)
    return unheader_vcf

def _get_vcf_writer(my_vcf_file, tmp_header):
    """
    :param my_vcf_file:
    :param my_vcf_header_file:
    :return: my_vcf_header file
    """
    cmd = f"bcftools view -h {my_vcf_file} >{tmp_header}.vcf"
    _run_command(cmd)
    vcf_reader = vcf.Reader(filename=f"{tmp_header}.vcf")
    vcf_writer = vcf.Writer(open(f"{tmp_header}_combined.vcf", 'w'), vcf_reader)
    return vcf_writer, f"{tmp_header}_combined.vcf"

def _run_command(cmd):
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
        (out, error) = cmd_obj.communicate()
        exit_code = cmd_obj.returncode
        if exit_code == 0:
            logging.debug(f"Command run successfully:\n{cmd}"
                          f"OUT:{out} Error:{error} Exit:{exit_code}\n")
        else:
            logging.info("Error: command exited with non zero exit \
                          status, please check logging file for more details")
            logging.error(f"OUT:{out}:Error:{error}:Exit:{exit_code}")
            if exit_code != 0:
                sys.exit("Exiting...")
        return
    except TimeoutExpired:
        cmd_obj.kill()
        (out, error) = cmd_obj.communicate()
        logging.error(f"Unable to run command:{cmd}: Out:{out} : Error:{error}")
        sys.exit(f"Unable to run command:{cmd}: Out:{out} : Error:{error}")
