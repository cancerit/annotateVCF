import logging
import sys
import os
import re
import io
import tempfile
import pkg_resources
import json
import shutil

from contextlib import contextmanager

'''
  This code runs bcftools to annotate driver gene and variant sites
'''


class IO_Formatter:
    """
        Main class , loads user defined parameters and files
    """

    def __init__(self, **kwargs):
        self.vcf_file = kwargs['vcf_file']
        self.genes_file = kwargs.get('lof_genes', None)
        self.muts_file = kwargs.get('mutations', None)
        self.cpv_file = kwargs.get('cancer_predisposition', None)
        self.np_vcf = kwargs.get('normal_panel', None)
        self.header_line = kwargs['header_line']
        self.outdir = kwargs['outdir']
        self.keepTmp = kwargs['keepTmp']
        self.json_file = kwargs.get('vcf_filters', None)
        # check input data ...

    def format(self, input_array):
        """

        :param input_array: array containing parameters to invoke formatter functions
        :return: formatter data for given parameter
        """
        format_store = {}
        for input_type in input_array:
            formatter = self._get_formatter(input_type)
            format_store[input_type] = formatter()
        return format_store

    def _get_formatter(self, input_type):
        if input_type == 'vcf_file':
            return self._format_vcf
        elif input_type == 'outdir':
            return self._get_outdir_path
        elif input_type == 'input_status':
            return self._check_input
        elif input_type == 'vcf_filters':
            return self._get_filters
        else:
            raise ValueError(input_type)

    def _check_input(self):
        """
        check user input files accessibility status
        :return: status dict
        """
        input_status = check_inputs({'vcf_file': self.vcf_file, 'normal_panel': self.np_vcf,
                                     'mutations': self.muts_file, 'lof_genes': self.genes_file, 'cancer_predisposition': self.cpv_file})
        if input_status['vcf_file'] is None:
            sys.exit("Please provide input vcf file")
        return input_status

    def _format_vcf(self):
        return get_file_metadata(self.vcf_file)

    def _get_filters(self):
        """
          load parameters from json config file
        """
        inc_filters = ['FORMAT', 'FILTER', 'INFO', 'INFO_FLAG_GERMLINE']
        formatted_filters = parse_filters(self.json_file, 'include', inc_filters)
        return formatted_filters

    def _get_outdir_path(self):
        """
         formatter function: add default output directory
         :return:
         """
        outputPath = self.outdir
        os.makedirs(outputPath, exist_ok=True)
        return outputPath

    # generic functions ....


def check_inputs(file_dict):
    """
    :param file_dict: check status of input files..
    :return: status dict
    """
    for f_type in file_dict.keys():
        if file_dict[f_type] and os.path.isfile(file_dict[f_type]):
            logging.debug(f"input: {f_type}  :{file_dict[f_type]}")
            file_dict[f_type] = True
        else:
            logging.debug(f"File not provided, skipping analysis step  : {f_type}")
            file_dict[f_type] = False
    return file_dict


def get_file_metadata(full_file_name):
    """
      takes file path as input and gives its path and processed extension
      #  If there are two extensions adds second extensions as prefix to first ext...
    """
    (_, name) = os.path.split(full_file_name)
    (name_no_ext, first_ext) = os.path.splitext(name)
    (name_no_ext2, second_ext) = os.path.splitext(name_no_ext)
    ext = second_ext + first_ext
    file_metadata = {'name': name_no_ext2, 'ext': ext, 'path': full_file_name}
    return file_metadata


def parse_filters(json_file, filter_type, filters):
    """
      load filtering parameters from json config file
    """
    filter_param_dict = {}
    try:
        if json_file is None:
            sys.exit('Json configuration file must be provided')
        with open(json_file, 'r') as cfgfile:
            filter_cfg = json.load(cfgfile)
            for filter in filters:
                filter_param_dict[filter] = filter_cfg[filter_type][filter]
    except json.JSONDecodeError as jde:
        sys.exit('json error:{}'.format(jde.args[0]))
    except FileNotFoundError as fne:
        sys.exit('Can not find json file:{}'.format(fne.args[0]))
    return filter_param_dict


@contextmanager
def tempdir(mypath):
    """

    :param mypath: path to create tmpdir
    :return:
    """
    path = tempfile.mkdtemp(dir=mypath)
    try:
        yield path
    finally:
        try:
            shutil.rmtree(path)
        except IOError:
            sys.stderr.write(f'Failed to clean up temp dir {path}')
