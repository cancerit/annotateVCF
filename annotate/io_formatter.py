import logging
import sys
import os
import re
import io
import tempfile
import pkg_resources
import shutil

from contextlib import contextmanager

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
        self.np_vcf = kwargs.get('normal_panel',None)
        self.germline_tag = kwargs['germline_tag']
        self.lof_type = kwargs['lof_type']
        self.outdir = kwargs['outdir']
        self.keepTmp = kwargs['keepTmp']
        self.filter = kwargs.get('vcf_filter',None) 
        # check input data ...

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
        elif input_type == 'normal_panel':
          return self._format_vcf_np
        elif input_type == 'vcf_file':
          return self._format_vcf
        elif input_type == 'lof_type':
          return self._format_lof
        elif input_type == 'outdir':
          return self._get_outdir_path
        elif input_type == 'input_status':
          return self._check_input
        else:
          raise ValueError(input_type)

    def _check_input(self):
         input_status = check_inputs( {'vcf_file':self.vcf_file, 'normal_panel':self.np_vcf, 
                                     'mutations':self.muts_file, 'lof_genes':self.gene_file})
         if input_status['vcf_file'] is None:
           sys.exit("Please provide input vcf file")
         return input_status
           
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

    def _format_vcf_np(self):
         f_filter=format_filter(self.filter, 'FILTER=')
         (file_name, ext)=get_file_metadata(self.np_vcf)
         return({'name':file_name, 'path':self.np_vcf, 'ext':ext, 'filter':f_filter, 'tag':self.germline_tag})

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

# generic functions ....
def check_inputs(file_dict):
  for f_type in file_dict.keys():
    if file_dict[f_type] and os.path.isfile(file_dict[f_type]):
      logging.info("{} input file :{}".format( f_type, file_dict[f_type]))
      file_dict[f_type] = True
    else:
      logging.info("File : {} not found, skipping analysis step for : {}".format(file_dict[f_type], f_type)) 
      file_dict[f_type] = False
  return file_dict

def format_filter(filter_vals, filter_type):
  format_store=[]
  for val in filter_vals:
    format_store.append(filter_type + "\"" + val + "\"" )
  return ' || '.join(format_store)

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

@contextmanager
def tempdir(mypath):
    path = tempfile.mkdtemp(dir=mypath)
    try:
      yield path
    finally:
      try:
        print("shutil.rmtree(path)")
      except IOError:
        sys.stderr.write('Failed to clean up temp dir {}'.format(path))
