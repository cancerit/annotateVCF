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
                   'lof_genes': test_dir + 'lof_genes_v1.0.txt',
                   'mutations': test_dir + 'driver_mutations_sorted.tsv.gz',
                   'vcf_filters': test_dir + 'filters.json',
                   'header_line': test_dir + 'info.header',
                   'outdir': test_dir + 'tmpout',
                   'keepTmp': False,
                   'normal_panel': test_dir+'np.vcf.gz'
                   }

  file_metadata = {'vcf_file': {'name': 'celline',
                                'ext': '.vcf.gz',
                                'path': test_dir + 'celline.vcf.gz'
                                }}
  file_dict = {'input_status': {'mutations': True, 'lof_genes': True,'cancer_predisposition': False,
                                'normal_panel': True, 'vcf_file': True}}

  info_filter="INFO/VC=\"stop_lost,start_lost,ess_splice,frameshift,nonsense\""
  format_filter= "FORMAT/VAF[*] > 0.15"
  filters_filter= "FILTER=\"PASS\""
  info_flag_germline="NPGL"
  # celline output
  muts_vcf = f"{test_out}/celline_muts.vcf.gz"
  lof_vcf = f"{test_out}/celline_genes_lof.vcf.gz"
  drv_vcf = f"{test_out}/celline_drv.vcf.gz"
  my_formatter=formatter.IO_Formatter(**options_vcf_celline) 
  outdir_path=my_formatter.format(['outdir'])

  def test_celline_vcf_input(self):
      # check input type function
    f = self.my_formatter
    assert self.file_metadata == f.format(['vcf_file']),'test_vcf_input test OK'

  def test_celline_file_input(self):
    file_dict=self.file_dict
    f = self.my_formatter
    assert file_dict == f.format(['input_status']),'test_file_input test OK'

  def test_celline_info_filter(self):
    f = self.my_formatter
    vcf_filters = f.format(['vcf_filters'])
    vcf_filter_params = vcf_filters.get('vcf_filters', None)
    assert self.info_filter == vcf_filter_params['INFO'],'test_INFO test OK'

  def test_celline_format_filter(self):
    f = self.my_formatter
    vcf_filters = f.format(['vcf_filters'])
    vcf_filter_params = vcf_filters.get('vcf_filters', None)
    assert self.format_filter == vcf_filter_params['FORMAT'],'test_FORMAT test OK'

  def test_celline_filter_filters(self):
    f = self.my_formatter
    vcf_filters = f.format(['vcf_filters'])
    vcf_filter_params = vcf_filters.get('vcf_filters', None)
    assert self.filters_filter == vcf_filter_params['FILTER'],'test_FILTER test OK'

  def test_celline_info_flag_germline(self):
    f = self.my_formatter
    vcf_filters = f.format(['vcf_filters'])
    vcf_filter_params = vcf_filters.get('vcf_filters', None)
    assert self.info_flag_germline == vcf_filter_params['INFO_FLAG_GERMLINE'],'test_INFO_FLAG_GERMLINE test OK'

  def chek_celline_outdir(slef):
       self.test_dir + 'tmpout' == self.outdir_path 

  def test_celline_vcf_formatter(self):
    f = self.my_formatter
    outdir_path = self.outdir_path
    with formatter.tempdir(outdir_path['outdir']) as base_dir:
      obs_muts_vcf = f"{base_dir}/celline_muts.vcf.gz"
      obs_lof_vcf = f"{base_dir}/celline_genes_lof.vcf.gz"
      obs_drv_vcf = f"{base_dir}/celline_drv.vcf.gz"

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
