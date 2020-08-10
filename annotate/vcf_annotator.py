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

class VcfAnnotator:
    """
        Main class , loads user defined parameters and files
    """
    def __init__(self, f, basedir):
       input_data = f.format(['vcf_file'])
       self.drv_mut = f.format(['mutations'])
       self.drv_genes = f.format(['lof_type', 'lof_genes'])
       
       self.germline_filter = f.format(['germline'])
       #set input vcf  paramaters and tmp_oudir.... required for downnstream functions....
       input_data = f.format(['vcf_file'])
       self.vcf_path = input_data['vcf_file']['path']
       self.vcf_ext = input_data['vcf_file']['ext']
       self.vcf_name = input_data['vcf_file']['name']
       self.vcf_filter = input_data['vcf_file']['filter']
       self.outdir = basedir
       self.outfile_name = basedir+'/'+self.vcf_name+'{}'

       print("input_data:{}\n drv_mut:{}\n drv_genes:{}\n germline_filter:{}\nBasedir:{}\n".format(input_data,self.drv_mut,self.drv_genes,self.germline_filter, basedir))
        
    def tag_germline_vars(self):
       print(self.germline_filter)   
    def annot_drv_muts(self):
       muts_outfile = self.outfile_name.format('_muts.vcf.gz')
       ANNOTATE_MUTS = 'bcftools annotate -i \'{}\' --merge-logic DRV:unique -a {} -h {} -c CHROM,FROM,TO,INFO/DRV  {}|' \
                'bcftools annotate  -i \'DRV!="." && DRV[*]==VC\' | ' \
                'bgzip -c >{} && tabix -f -p vcf {}'
       #{'driver_mutations_test': ['##INFO=<ID=DRV,Number=., Type=String, Description="Driver Variant Class">', 'tests/test_drvData/driver_mutations_test.tsv.gz']}
       input_param = get_parameters(self.drv_mut['mutations'], self.outdir)
       
       for drv_muts_file, header_file in input_param.items():
         cmd = ANNOTATE_MUTS.format(self.vcf_filter, drv_muts_file, header_file, self.vcf_path, muts_outfile, muts_outfile)
         _run_command(cmd)         
         #_run_command('cp -p '+self.outfile+' '+self.outdir+'/..') 
       return muts_outfile 

    

    def annotate_lof_genes(self):
      get_gene=re.compile('.*;VD=(\w+)|.*')
      ANNOTATE_GENES = 'bcftools annotate -i \'{}\' -a {} -i \' {} \' -h {} -c CHROM,FROM,TO,INFO/DRV {} >{}' 
      # create dummy genome locationo file to annoate LoF genes...
      genome_loc_file = self.outdir+'/genome.tab.gz'
      create_dummy_genome(self.vcf_path, genome_loc_file) 
      input_param = get_parameters(self.drv_genes['lof_genes'], self.outdir)
      genes_outfile = self.outfile_name.format('_genes.vcf')
      lof_outfile = self.outfile_name.format('_genes_lof.vcf')

      for drv_gene_file, header_file in input_param.items():
        lof_gene_list = get_drv_gene_list(drv_gene_file)
        #map lof effect types to pass variants...
        cmd = ANNOTATE_GENES.format(self.vcf_filter, genome_loc_file, self.drv_genes['lof_type'], header_file, self.vcf_path, genes_outfile)
        _run_command(cmd)
        with open(lof_outfile, "w") as lof_fh, open(genes_outfile, 'r') as gene_f:
          for line in gene_f:
            # write header lines
            if line.startswith('#'):
              lof_fh.write(line)
            else:
              gene = get_gene.match(line)[1]
              # write matching LoF genes....
              if gene in lof_gene_list:
                lof_fh.write(line)

      return compress_vcf(lof_outfile)       

    def concat_results(self, muts_vcf, lof_vcf):
      concat_drv_out = self.outfile_name.format('_drv.vcf.gz')  
      CONCAT_VCF = 'bcftools concat --allow-overlaps --rm-dups all {} {} {} | bcftools sort | ' \
      'bgzip -c >{} && tabix -f -p vcf {}'
      cmd=CONCAT_VCF.format(muts_vcf, lof_vcf, self.vcf_path, concat_drv_out, concat_drv_out)
      _run_command(cmd)
      _run_command('cp -p '+concat_drv_out+'  '+self.outdir+'/..')
      return concat_drv_out

def get_parameters(param_dict,outdir):
    input_param = {}
    print(param_dict)
    for file_name, (header_string, input_file, ext)  in param_dict.items():
        header_file = outdir+'/'+file_name+'.header'
        with open(header_file, 'w') as hf:
           hf.write(header_string)
        input_param[input_file] = header_file
    return input_param

def get_drv_gene_list(drv_genes):
    with open(drv_genes) as f_drv:
      lof_gene_list = f_drv.read().splitlines()
    return lof_gene_list

def compress_vcf(vcf):
  outfile=vcf+'.gz'
  BGZIP_TABIX = 'bgzip -c <{} >{} && tabix -f -p vcf {}'
  cmd = BGZIP_TABIX.format(vcf, outfile, outfile)
  _run_command(cmd)
  return outfile
  
def create_dummy_genome(input_vcf, genome_loc):
    DUMMY_GENOME_TAB='tabix -l {} | xargs -I vcf_chr printf \'vcf_chr\t1\t400000000\tLoF\n\' | bgzip -c >{} && tabix -s1 -b2 -e3 {}'
    cmd = DUMMY_GENOME_TAB.format(input_vcf, genome_loc, genome_loc)
    _run_command(cmd)
    
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
