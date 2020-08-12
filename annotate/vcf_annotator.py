import logging
import sys
import os
import re
import io
import copy
from subprocess import Popen, PIPE, STDOUT
import pkg_resources

'''
  This code runs bcftools to annoate driver gene and variant sites
'''

class VcfAnnotator:
    """
        Main class , loads user defined parameters and files
    """
    def __init__(self, f, basedir):
       self.outdir = basedir
       self.status_dict=f.format(['input_status'])
       self.input_data = f.format(['vcf_file'])
       # set input vcf parameters ...
       self._set_input_vcf(self.input_data)
       self.merge_vcf_dict={}
       # run analysis
       self._runAnalysis(f)
    
    def _runAnalysis(self, f):
        status=self.status_dict['input_status']
        run_status = False
        a_type= ['normal_panel', 'mutations', 'lof_genes']
        for analysis in a_type:
           if status[analysis] and analysis == 'normal_panel':
             np = f.format([analysis])
             self.tag_germline_vars(np[analysis])
             run_status = True
             logging.info("Analysis completed for : {}".format(analysis)) 
           if status[analysis] and analysis == 'mutations':
             drv_mut = f.format([analysis])
             self.annot_drv_muts(drv_mut[analysis])
             run_status = True
             logging.info("Analysis completed for : {}".format(analysis)) 
           if status[analysis] and analysis == 'lof_genes':
             drv_genes = f.format(['lof_type', 'lof_genes'])
             self.annotate_lof_genes(drv_genes)
             run_status = True
             logging.info("Analysis completed for : {}".format(analysis)) 
        if run_status:       
          logging.info("concatenating results") 
          self.concat_results() 
        else:
          logging.info("Input files not accessible analysis aborted")
     
    def _set_input_vcf(self,input_data):
       self.vcf_path = input_data['vcf_file']['path']
       self.vcf_ext = input_data['vcf_file']['ext']
       self.vcf_name = input_data['vcf_file']['name']
       self.vcf_filter = input_data['vcf_file']['filter']
       self.outfile_name = self.outdir + '/'+self.vcf_name+'{}'
   
    def tag_germline_vars(self, np):
         self.tagged_vcf = self.outfile_name.format('_np.vcf.gz') 
         self.filtered_vcf = self.outfile_name.format('_np_filtered.vcf.gz')         

         TAG_GERMLINE='bcftools annotate -a {} -i \'{}\'  -m \'{}\'  {} | bgzip -c >{} && tabix -p vcf {}'
         cmd=TAG_GERMLINE.format(np['path'], np['filter'], np['tag'], self.vcf_path, self.tagged_vcf, self.tagged_vcf)
         _run_command(cmd)
         
         FILTER_GERMLINE ='bcftools  view -i \'{} && {}=0\' {} | bgzip -c >{} && tabix -p vcf {}'
         cmd = FILTER_GERMLINE.format(np['filter'], np['tag'], self.tagged_vcf, self.filtered_vcf, self.filtered_vcf)
         
         _run_command(cmd)
         
         # update input file for driver annotation step ... 
         input_data = copy.deepcopy(self.input_data)
         input_data['vcf_file']['path'] = self.filtered_vcf
         self._set_input_vcf(input_data)
          
         # add vcf path to concat 
         self.merge_vcf_dict['c']= self.filtered_vcf
         self.merge_vcf_dict['d']= self.tagged_vcf
   
    def annot_drv_muts(self, drv_mut):
       muts_outfile = self.outfile_name.format('_muts.vcf.gz')
       ANNOTATE_MUTS = 'bcftools annotate -i \'{}\' --merge-logic DRV:unique -a {} -h {} -c CHROM,FROM,TO,INFO/DRV  {}|' \
                'bcftools annotate  -i \'DRV!="." && DRV[*]==VC\' | ' \
                'bgzip -c >{} && tabix -f -p vcf {}'
       #{'driver_mutations_test': ['##INFO=<ID=DRV,Number=., Type=String, Description="Driver Variant Class">', 'tests/test_drvData/driver_mutations_test.tsv.gz']}
       input_param = get_parameters(drv_mut, self.outdir)
       
       for drv_muts_file, header_file in input_param.items():
         cmd = ANNOTATE_MUTS.format(self.vcf_filter, drv_muts_file, header_file, self.vcf_path, muts_outfile, muts_outfile)
         _run_command(cmd)         
         #_run_command('cp -p '+self.outfile+' '+self.outdir+'/..')
       self.merge_vcf_dict['a']=muts_outfile 

    def annotate_lof_genes(self, drv_genes):
      get_gene = re.compile('.*;VD=(\w+)|.*')
      ANNOTATE_GENES = 'bcftools annotate -a {} -i \'{} && ({})\' -h {} -c CHROM,FROM,TO,INFO/DRV {} >{}' 
      # create dummy genome locationo file to annoate LoF genes...
      genome_loc_file = self.outdir+'/genome.tab.gz'
      create_dummy_genome(self.vcf_path, genome_loc_file) 
      input_param = get_parameters(drv_genes['lof_genes'], self.outdir)
      genes_outfile = self.outfile_name.format('_genes.vcf')
      lof_outfile = self.outfile_name.format('_genes_lof.vcf')

      for drv_gene_file, header_file in input_param.items():
        lof_gene_list = get_drv_gene_list(drv_gene_file)
        #map lof effect types to pass variants...
        cmd = ANNOTATE_GENES.format(genome_loc_file, self.vcf_filter, drv_genes['lof_type'], header_file, self.vcf_path, genes_outfile)
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
      self.merge_vcf_dict['b']=compress_vcf(lof_outfile)

    def concat_results(self):
      concat_drv_out = self.outfile_name.format('_drv.vcf.gz')  
      self.merge_vcf_dict['e']=self.input_data['vcf_file']['path']
 
      vcf_files=([self.merge_vcf_dict[myvcf] for myvcf in sorted(self.merge_vcf_dict.keys())])
      
      CONCAT_VCF = 'bcftools concat --allow-overlaps --rm-dups all {} | bcftools sort | ' \
      'bgzip -c >{} && tabix -f -p vcf {}'
      cmd=CONCAT_VCF.format(' '.join(vcf_files), concat_drv_out, concat_drv_out)
      _run_command(cmd)
      _run_command('cp -p '+concat_drv_out+'  '+concat_drv_out+'.tbi '+self.outdir+'/..')

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
                          status, please check logging file for more details")
            logging.error("OUT:{}:Error:{}:Exit:{}".format(out, error, exit_code))
        return
    except OSError as oe:
        logging.error("Unable to run command:{} Error:{}".format(cmd, oe.args[0]))
        sys.exit("Unable to run command:{} Error:{}".format(cmd, oe.args[0]))