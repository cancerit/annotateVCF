import annotate.formatInput as fc
import annotate.staticMethods as sm
import pytest
import os
import unittest.mock

'''
written test to check codebase integrity
of annotateVcf
'''

class TestClass():
    configdir = os.path.join(os.path.dirname(os.path.realpath(__file__)))
    driver_data = configdir + '/test_drvData'
    driver_json = configdir + '/test_drvData.json'
    vcf_file= configdir + '/test_input_vcf/input.vcf.gz'

    def test_file_input(self):
        # check input type function
        my_input=fc.AnnotateVcf(vcf_file=self.vcf_file,driver_data=self.driver_data, driver_json=self.driver_json, outdir="../out_test", keepTmp=True)
        print(self.vcf_file)
        assert ['y'] == my_input.check_input(),'input_check test OK'
