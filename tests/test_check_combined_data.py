import annotate.formatInput as fc
import annotate.staticMethods as sm
import pytest
import os, tempfile
import filecmp

'''
written test to check codebase integrity
of archCompare
'''


class TestClass():
    pass

    def test_static_methods(self):
        
        test_drv_json = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_drvData.json')
        configdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_drvData')
        configdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_drvData')
        header_info = configdir + 'vcf_info.txt'
        genome_tab = configdir + 'genome_grch38.tab.gz'
        drv_gene = configdir + 'lof_genes_v1.0.txt'
        drv_gene_prev = configdir + 'lof_genes_previous_symbol.v1.0.txt'
        drv_mut = configdir + 'driver_mutations_grch38_v1.0.tsv.gz'
        lof_consq = configdir + 'lof_consequences.txt'
        vcf_file=testdir + 'input.vcf.gz'
        lof_con=['stop_lost','splice_site']
        # check input type function
        mystatic_obj = sm.StaticMthods()
        exp_con = mystatic_obj.format_consequences(lof_con)
        assert ('INFO/VC="stop_lost" || INFO/VC="splice_site"') == exp_con, 'format_consequences OK'
        # assert True == result, 'process_segments: check results'


if __name__ == '__main__':
    mytests = TestClass()
    mytests()
