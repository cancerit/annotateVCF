import annotate.io_formatter as formatter
import annotate.vcf_annotator as annotator
import pytest
import os
import filecmp
import unittest.mock

'''
written test to check codebase integrity
of annotateVcf
'''

class TestClass():
  configdir = os.path.join(os.path.dirname(os.path.realpath(__file__)))
  test_dir = configdir + '/test_input/'
  test_out = configdir + '/test_output/'

  options_vcf_celline = {'vcf_file': test_dir + 'celline.vcf.gz',
                   'vcf_filter': ['PASS'],
                   'lof_genes': test_dir + 'lof_genes_v1.0.txt',
                   'mutations': test_dir + 'driver_mutations_sorted.tsv.gz',
                   'lof_type': ["stop_lost","start_lost","ess_splice","frameshift","nonsense"],
                   'header_line': test_dir + 'info.header',
                   'outdir': test_dir + 'tmpout',
                   'keepTmp': False,
                   'germline_tag': 'NPGL',
                   'normal_panel': test_dir+'np.vcf.gz'
                   }

  file_metadata = {'vcf_file': {'name': 'celline',
                                'ext': '.vcf.gz',
                                'path': test_dir + 'celline.vcf.gz'
                                }}
  file_dict = {'input_status': {'mutations': True, 'lof_genes': True,
                                'normal_panel': True, 'vcf_file': True}}
  lof_type = {
      'lof_type': 'INFO/VC="stop_lost" || INFO/VC="start_lost" || INFO/VC="ess_splice" || INFO/VC="frameshift" || INFO/VC="nonsense"'}
  my_filter = {'format_filter': 'FILTER="PASS"'}
  # celline output
  muts_vcf = test_out + 'celline_muts.vcf.gz'
  lof_vcf = test_out + 'celline_genes_lof.vcf.gz'
  drv_vcf = test_out + 'celline_drv.vcf.gz'


  my_formatter=formatter.IO_Formatter(**options_vcf_celline) 
  outdir_path=my_formatter.format(['outdir'])

  def test_celline_vcf_input(self):
      # check input type function
    f = self.my_formatter
    file_metadata = self.file_metadata
    assert file_metadata == f.format(['vcf_file']),'test_vcf_input test OK'

  def test_celline_file_input(self):
    file_dict=self.file_dict
    f = self.my_formatter
    assert file_dict == f.format(['input_status']),'test_file_input test OK'

  def test_celline_lof_format(self):
    lof_type=self.lof_type
    f = self.my_formatter
    assert lof_type == f.format(['lof_type']),'test_lof_format test OK'

  def test_celline_filter_format(self):
    my_filter=self.my_filter
    f = self.my_formatter
    assert my_filter == f.format(['format_filter']),'test_filter_format test OK'
  def chek_celline_outdir(slef):
       self.test_dir + 'tmpout' == self.outdir_path 

  def test_celline_vcf_formatter(self):
    f = self.my_formatter
    outdir_path = self.outdir_path
    with formatter.tempdir(outdir_path['outdir']) as base_dir:
      obs_muts_vcf = base_dir + '/celline_muts.vcf.gz'
      obs_lof_vcf = base_dir + '/celline_genes_lof.vcf.gz'
      obs_drv_vcf = base_dir + '/celline_drv.vcf.gz'

      annotator.VcfAnnotator(f, base_dir)

      exp_muts_vcf_sub = annotator.unheader_vcf(self.muts_vcf, base_dir+'/exp_muts.vcf')
      obs_muts_vcf_sub = annotator.unheader_vcf(obs_muts_vcf, base_dir+'/obs_muts.vcf')
      assert filecmp.cmp(exp_muts_vcf_sub, obs_muts_vcf_sub, shallow=True), 'muts records in vcf files are identical OK'

      exp_lof_vcf_sub = annotator.unheader_vcf(self.lof_vcf, base_dir+'/exp_genes.lof.vcf')
      obs_lof_vcf_sub = annotator.unheader_vcf(obs_lof_vcf, base_dir+'/obs_genes.lof.vcf')
      assert filecmp.cmp(exp_lof_vcf_sub, obs_lof_vcf_sub, shallow=True), 'lof records in vcf files are identical OK'

      exp_drv_vcf_sub = annotator.unheader_vcf(self.drv_vcf, base_dir+'/exp_drv.vcf')
      obs_drv_vcf_sub = annotator.unheader_vcf(obs_drv_vcf, base_dir+'/obs_drv.vcf')
      assert filecmp.cmp(exp_drv_vcf_sub, obs_drv_vcf_sub, shallow=True), 'final drv records in vcf files are identical OK'
