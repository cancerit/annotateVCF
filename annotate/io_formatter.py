import logging
import sys
import os
import re
import io
import json
import tempfile
from subprocess import Popen, PIPE, STDOUT
import pkg_resources
import shutil

from contextlib import contextmanager

log = logging.getLogger(__name__)

'''
  This code runs bcftools to annoate driver gene and variant sites
'''

class IO_Formatter:
    """
        Main class , loads user defined parameters and files
    """
    def __init__(self, **kwargs):
        self.vcf_file = kwargs['vcf_file']
        self.gene_file = kwargs.get('lof_genes',None)
        self.muts_file = kwargs.get('mutations',None)
        self.germline_files = kwargs.get('germline',None)
        self.lof_type = kwargs['lof_type']
        self.outdir = kwargs['outdir']
        self.keepTmp = kwargs['keepTmp']
        self.filter = kwargs.get('vcf_filter',None) 
        
    def format(self,input_array):
        format_store={}
        for input_type in input_array:
          formatter = self._get_formatter(input_type)     
          format_store[input_type]=formatter()
        return format_store

    def _get_formatter(self, input_type):
        if input_type == 'lof_genes':
          return self._format_header_genes
        elif input_type == 'mutations':
          return self._format_header_muts
        elif input_type == 'germline':
          return self._format_header_germline
        elif input_type == 'vcf_file':
          return self._format_vcf
        elif input_type == 'lof_type':
          return self._format_lof
        elif input_type == 'outdir':
          return self._get_outdir_path
        else:
          raise ValueError(input_type)

    def _format_header_genes(self):
         header_lines={}
         file_name, header_line, ext = format_header_driver(self.gene_file)
         header_lines[file_name]= [header_line, self.gene_file, ext]
         return header_lines

    def _format_header_muts(self):
         header_lines={}
         file_name, header_line, ext = format_header_driver(self.muts_file)
         header_lines[file_name]= [header_line, self.muts_file, ext]
         return header_lines

    def _format_header_germline(self):
         header_lines={}
         count = 0
         for file_path in self.germline_files:
            count += 1
            file_name, header_line, ext = format_header_germline(file_path, count)
            header_lines[file_name] = [header_line, file_path, ext]
         return header_lines
    
    def _format_vcf(self):
         f_filter=format_filter(self.filter, 'FILTER=')
         (file_name, ext)=get_file_metadata(self.vcf_file)
         return({'name':file_name, 'path':self.vcf_file, 'ext':ext, 'filter':f_filter})

    def _format_lof(self):
        return format_filter(self.lof_type, 'INFO/VC=')   

    def _get_outdir_path(self):
        outputPath = self.outdir + '/out_anotatevcf/'
        os.makedirs(outputPath, exist_ok=True)
        return outputPath 


def get_drv_gene_list(drv_genes):
    with open(drv_genes) as f_drv:
      lof_gene_list = f_drv.read().splitlines()
    return lof_gene_list

def compress_vcf(vcf):
  vcf_out=vcf+'.gz'
  BGZIP_TABIX = 'bgzip -c <{} >{} && tabix -f -p vcf {}'
  cmd = BGZIP_TABIX.format(vcf, vcf_out, vcf_out)
  _run_command(cmd)
  return vcf_out
  
def create_dummy_genome(input_vcf, genome_loc):
    DUMMY_GENOME_TAB='tabix -l {} | xargs -I vcf_chr printf \'vcf_chr\t1\t400000000\tLoF\n\' | bgzip -c >{} && tabix -s1 -b2 -e3 {}'
    cmd = DUMMY_GENOME_TAB.format(input_vcf, genome_loc, genome_loc)
    _run_command(cmd)
    
def format_filter(filter_vals,filter_type):
    format_store=[]
    for val in filter_vals:
        format_store.append(filter_type + "\"" + val + "\"" )
    return ' || '.join(format_store)

def get_parameters(param_dict,outdir):
    input_param = {}
    for file_name, (header_string, input_file)  in param_dict.items():
        header_file = outdir+'/'+file_name+'.header'
        with open(header_file, 'w') as hf:
           hf.write(header_string)
        input_param[input_file] = header_file
    return input_param

def format_header_germline(myfile,counter):
    """
       generate vcf header and INFO field
    """
    file_name, ext = get_file_metadata(myfile)
    #acronym = filter(str.isupper, file_name.title())
    #ID=''.join([i for i in acronym])
    header_line = "##INFO=<ID=GVF{},Number=., Type=String, Description=\"Germline Variant Filter ({})\">".format(counter, file_name)
    return (file_name, header_line, ext)

def format_header_driver(myfile):
    """
       generate vcf heade and INFO filed
    """
    file_name, ext = get_file_metadata(myfile)
    header_line = "##INFO=<ID=DRV,Number=., Type=String, Description=\"Driver Variant Class\">"
    return (file_name, header_line, ext)

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
        # logging.info("running command:{}".format(cmd))
        (out, error) = cmd_obj.communicate()
        exit_code = cmd_obj.returncode
        if (exit_code == 0):
            logging.info("Command run successfully:\n{} OUT:{} Error:{} Exit:{}\n".format(cmd, out, error, exit_code))
        else:
            logging.debug("Error: command exited with non zero exit \
                          status, please check log file for more details")
            logging.error("OUT:{}:Error:{}:Exit:{}".format(out, error, exit_code))
        return
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
        sys.stderr.write('Failed to clean up temp dir {}'.format(path))
