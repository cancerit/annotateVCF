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
        configdir = os.path.join(os.path.dirname(os.path.realpath(__file__)))
        driver_data = configdir + '/test_drvData/'
        output_data = configdir + '/test_output/'
        # input files...
        driver_json = configdir + '/test_drvData.json'
        vcf_file = configdir + '/test_input_vcf/input.vcf.gz'
        drv_genes = driver_data + 'lof_genes_v1.0.txt' 
        drv_genes_prev = driver_data + 'lof_genes_previous_symbol.v1.0.txt' 
        drv_muts = driver_data + 'driver_mutations_test.tsv.gz' 
        header_file = driver_data + 'vcf_info.txt' 
        genome_loc = driver_data + 'genome_grch38.tab.gz' 
        lof_consequences = ["stop_lost","start_lost","ess_splice","frameshift","nonsense"]
        # output files ...
        pass_vcf = output_data + 'input.pass.vcf.gz'
        genes_lof_vcf = output_data + 'input.genes.lof.vcf.gz'
        muts_vcf = output_data + 'input.muts.vcf.gz'
        drv_vcf = output_data + 'input.drv.vcf.gz'
        outputPath='out_test'
        os.makedirs(outputPath, exist_ok=True)
        file_name='input'
        ext='.vcf.gz'
        outfile_name = outputPath + file_name

        exp_drv_gene_list=['APC', 'B2M']                                           
        exp_drv_gene_prev_dict={'ATD': 'ATM', 'ATC': 'ATM', 'ATDC': 'ATM', 'ATA': 'ATM'}
        # create object for static class ....
        so = sm.StaticMthods()

        exp_metadata =  (file_name, ext)
        obs_metadata = so.get_file_metadata(vcf_file)
        assert exp_metadata == obs_metadata, 'get_file_metadata test OK'

        exp_ref_data=(drv_genes, drv_genes_prev, drv_muts, header_file,genome_loc, lof_consequences)
        obs_ref_data = so.prepare_ref_data(driver_json,driver_data)
        assert exp_ref_data == obs_ref_data, 'prepare_ref_data test OK'

        exp_consequences = 'INFO/VC="stop_lost" || INFO/VC="start_lost" || INFO/VC="ess_splice" || INFO/VC="frameshift" || INFO/VC="nonsense"'
        obs_consequences = so.format_consequences(lof_consequences)
        assert exp_consequences == obs_consequences, 'format_consequence OK'
        
        obs_lof_gene_list = so.get_drv_gene_list(drv_genes) 
        assert exp_drv_gene_list == obs_lof_gene_list, 'get_drv_gene_list OK'

        obs_drv_gene_prev_dict = so.get_drv_prev_gene_dict(drv_genes_prev) 
        assert exp_drv_gene_prev_dict == obs_drv_gene_prev_dict, 'get_drv_prev_gene_dict OK'
        obs_pass_vcf = so.get_filtered_vcf(vcf_file, outfile_name+'.pass'+ext)
        (obs_vcf,exp_vcf)=so.test_compare_vcf(obs_pass_vcf,pass_vcf,outputPath+'/test_pass')
        assert filecmp.cmp(exp_vcf, obs_vcf,
                           shallow=True), 'passe records in vcf files are identical'
        
       

if __name__ == '__main__':
    mytests = TestClass()
    mytests()
