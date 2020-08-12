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
        Annotator class, does all the processing and output generation
    """
    def __init__(self, f, basedir):
        """
        initialises required parameters..
        :param f: formatted input object
        :param basedir: tmpdir
        """
        self.outdir = basedir
        self.status_dict=f.format(['input_status'])
        self.input_data = f.format(['vcf_file'])

        # set input vcf parameters ...
        self._set_input_vcf(self.input_data)
        self.drv_header = self.get_drv_header_file()
        self.merge_vcf_dict={}
        # run analysis
        self._runAnalysis(f)

    def _runAnalysis(self, f):
        """
        runs selected  the analysis steps based on user input
        :param f: formatted input object
        :return:
        """
        status=self.status_dict['input_status']
        vcf_filter_data = f.format(['format_filter'])
        self.vcf_filter = vcf_filter_data['format_filter']
        run_status = False
        a_type= ['normal_panel', 'mutations', 'lof_genes']
        
        for analysis in a_type:
            if status[analysis] and analysis == 'normal_panel':
                self.np_tag = f.np_tag
                self.np_vcf = f.np_vcf
                self.tag_germline_vars()
                run_status = True
                logging.info("Analysis completed for : {}".format(analysis))
            if status[analysis] and analysis == 'mutations':
                self.muts_file = f.muts_file
                self.annot_drv_muts()
                run_status = True
                logging.info("Analysis completed for : {}".format(analysis))
            if status[analysis] and analysis == 'lof_genes':
                self.genes_file = f.genes_file
                lof_types_data = f.format(['lof_type'])
                self.lof_types = lof_types_data['lof_type']
                self.annotate_lof_genes()
                run_status = True
                logging.info("Analysis completed for : {}".format(analysis))
        if run_status:
            logging.info("concatenating results")
            self.concat_results()
        else:
            logging.info("Input files not accessible analysis aborted")

    def _set_input_vcf(self,input_data):
        """
        sets input vcf data
        :param input_data:
        :return:
        """
        self.vcf_path = input_data['vcf_file']['path']
        self.vcf_ext = input_data['vcf_file']['ext']
        self.vcf_name = input_data['vcf_file']['name']
        self.outfile_name = self.outdir + '/'+self.vcf_name+'{}'

    def tag_germline_vars(self):
        """
        use normal panel to tag germline variants and created filtered vcf file to use in subsequent  driver annotation steps ...
        sets filtered vcf as new user input parameter for downstream analysis
        add tagged and filtered vcf files to concat in final step
        :param np: normal panel parameters
        :return:
        """
        
        self.tagged_vcf = self.outfile_name.format('_np.vcf.gz')
        self.filtered_vcf = self.outfile_name.format('_np_filtered.vcf.gz')

        TAG_GERMLINE='bcftools annotate -a {} -i \'{}\'  -m \'{}\'  {} | bgzip -c >{} && tabix -p vcf {}'
        cmd=TAG_GERMLINE.format(self.np_vcf, self.vcf_filter, self.np_tag, self.vcf_path, self.tagged_vcf, self.tagged_vcf)
        _run_command(cmd)

        FILTER_GERMLINE ='bcftools  view -i \'{} && {}=0\' {} | bgzip -c >{} && tabix -p vcf {}'
        cmd = FILTER_GERMLINE.format(self.vcf_filter, self.np_tag, self.tagged_vcf, self.filtered_vcf, self.filtered_vcf)

        _run_command(cmd)

        # update input file for driver annotation step ...
        input_data = copy.deepcopy(self.input_data)
        input_data['vcf_file']['path'] = self.filtered_vcf
        self._set_input_vcf(input_data)

        # add vcf path to concat
        self.merge_vcf_dict['c']= self.filtered_vcf
        self.merge_vcf_dict['d']= self.tagged_vcf

    def annot_drv_muts(self):
        """
        annotate vcf using driver mutation set
        add annotated vcf file to concat in final step
        :param drv_mut: reference mutation parameters
        :return:
        """
        muts_outfile = self.outfile_name.format('_muts.vcf.gz')
        ANNOTATE_MUTS = 'bcftools annotate -i \'{}\' --merge-logic DRV:unique -a {} -h {} -c CHROM,FROM,TO,INFO/DRV  {}|' \
                        'bcftools annotate  -i \'DRV!="." && DRV[*]==VC\' | ' \
                        'bgzip -c >{} && tabix -f -p vcf {}'
        #{'driver_mutations_test': ['##INFO=<ID=DRV,Number=., Type=String, Description="Driver Variant Class">', 'tests/test_drvData/driver_mutations_test.tsv.gz']}
        cmd = ANNOTATE_MUTS.format(self.vcf_filter, self.muts_file, self.drv_header, self.vcf_path, muts_outfile, muts_outfile)
        _run_command(cmd)
        self.merge_vcf_dict['a']=muts_outfile

    def annotate_lof_genes(self):
        """
        annotate vcf using known LoF genes
        :param drv_genes:
        :return:
        """
        get_gene = re.compile('.*;VD=(\w+)|.*')
        ANNOTATE_GENES = 'bcftools annotate -a {} -i \'{} && ({})\' -h {} -c CHROM,FROM,TO,INFO/DRV {} >{}'
        # create dummy genome locationo file to annoate LoF genes...
        genome_loc_file = self.outdir+'/genome.tab.gz'
        create_dummy_genome(self.vcf_path, genome_loc_file)
        genes_outfile = self.outfile_name.format('_genes.vcf')
        lof_outfile = self.outfile_name.format('_genes_lof.vcf')
        lof_gene_list = get_drv_gene_list(self.genes_file)
        #map lof effect types to pass variants...
        cmd = ANNOTATE_GENES.format(genome_loc_file, self.vcf_filter, self.lof_types, self.drv_header, self.vcf_path, genes_outfile)
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

    def get_drv_header_file(self):
        """
        create header annotation tag for driver annotations
        :param header_id:
        :param header_desc:
        :return:
        """
        header_file = self.outfile_name.format('_drv.header')
        with open(header_file, 'w') as hf:
            hf.write('##INFO=<ID=DRV,Number=., Type=String, Description="Driver Variant Class">')
        return header_file

def get_drv_gene_list(drv_genes):
    with open(drv_genes) as f_drv:
        lof_gene_list = f_drv.read().splitlines()
    return lof_gene_list

def compress_vcf(vcf):
    """
    :param vcf:
    :return:
    """
    outfile=vcf+'.gz'
    BGZIP_TABIX = 'bgzip -c <{} >{} && tabix -f -p vcf {}'
    cmd = BGZIP_TABIX.format(vcf, outfile, outfile)
    _run_command(cmd)
    return outfile

def create_dummy_genome(input_vcf, genome_loc):
    """
    :param input_vcf:
    :param genome_loc:
    :return:
    """
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
