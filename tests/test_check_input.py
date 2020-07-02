import annotate.formatInput as fc
import annotate.staticMethods as sm
import pytest
import os

'''
written test to check codebase integrity
of annotateVcf
'''

class TestClass():
    configdir = os.path.join(os.path.dirname(os.path.realpath(__file__)))
    driver_data = configdir + 'test_drvData.json'
    vcf_file= configdir + 'test_input/input.vcf.gz'

    def test_file_input(self):
        # check input type function
        #for infile in ([self.vcf_file, self.drv_mut, self.genome_tab, self.lof_consq , self.header_info] ) :
        my_input=fc.AnnotateVcf(vcf_file=self.vcf_file,driver_data=self.driver_data, outdir="../out_test", keepTmp=True)
        assert ['y'] == my_input.check_input(),'input_check test OK'

if __name__ == '__main__':
  mytests=TestClass()
  mytests()
