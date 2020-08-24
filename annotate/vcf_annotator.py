import logging
import sys
import os
import re
import io
import copy
from subprocess import Popen, PIPE, STDOUT
import pkg_resources

'''
  This code runs bcftools to annotates driver gene and variant sites
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
        self.status_dict = f.format(['input_status'])
        self.input_data = f.format(['vcf_file'])

        # set input vcf parameters ...
        self._set_input_vcf(self.input_data)
        self.drv_header = f.header_line
        self.merge_vcf_dict = {}
        self.keepTmp = f.keepTmp
        # run analysis
        self._runAnalysis(f)

    def _runAnalysis(self, f):
        """
        runs selected  the analysis steps based on user input
        :param f: formatted input object
        :return:
        """
        status = self.status_dict['input_status']
        vcf_filters = f.format(['format_filter'])
        self.vcf_filter = vcf_filters['format_filter']
        run_status = False
        a_type = ['normal_panel', 'mutations', 'lof_genes']

        for analysis in a_type:
            if status[analysis] and analysis == 'normal_panel':
                logging.info(f"Tagging germline variants with INFO tag:{f.np_tag}")
                self.tag_germline_vars(f.np_tag, f.np_vcf)
                run_status = True
            if status[analysis] and analysis == 'mutations':
                logging.info("Annotating driver mutations with INFO field:DRV=<consequence(s)>")
                self.annot_drv_muts(f.muts_file)
                run_status = True
            if status[analysis] and analysis == 'lof_genes':
                logging.info("Annotating LoF genes with INFO filed:DRV=LoF")
                lof_types = f.format(['lof_type'])
                self.annotate_lof_genes(f.genes_file, lof_types['lof_type'])
                run_status = True
        if run_status:
            logging.info("concatenating results")
            self.concat_results()
            logging.info("Analysis completed successfully")
        else:
            logging.info("Input files not accessible analysis aborted")

    def _set_input_vcf(self, input_data):
        """
        sets input vcf data
        :param input_data:
        :return:
        """
        self.vcf_path = input_data['vcf_file']['path']
        self.vcf_ext = input_data['vcf_file']['ext']
        self.vcf_name = input_data['vcf_file']['name']
        self.outfile_name = self.outdir + '/' + self.vcf_name + '{}'

    def tag_germline_vars(self, np_tag, np_vcf):
        """
        use normal panel to tag germline variants and created filtered vcf file to use in
        subsequent driver annotation steps ...
        sets filtered vcf as new user input parameter for downstream analysis
        add tagged and filtered vcf files to concat in final step
        :param np_tag:
        :param np_vcf:
        :return:
        """
        tagged_vcf = self.outfile_name.format('_np.vcf.gz')
        filtered_vcf = self.outfile_name.format('_np_filtered.vcf.gz')

        cmd = f"bcftools annotate -a {np_vcf} -i '{self.vcf_filter}'  -m '{np_tag}'" \
              f" {self.vcf_path} | bgzip -c >{tagged_vcf} && tabix -p vcf {tagged_vcf}"
        _run_command(cmd)

        cmd = f"bcftools  view -i '({self.vcf_filter}) && {np_tag}=0' {tagged_vcf} | " \
              f"bgzip -c >{filtered_vcf} && tabix -p vcf {filtered_vcf}"

        _run_command(cmd)

        # update input file for driver annotation step ...
        input_data = copy.deepcopy(self.input_data)
        input_data['vcf_file']['path'] = filtered_vcf
        self._set_input_vcf(input_data)

        # add vcf path to concat in the  order of priority...
        self.merge_vcf_dict['c'] = filtered_vcf
        self.merge_vcf_dict['d'] = tagged_vcf

    def annot_drv_muts(self, muts_file):
        """
        annotate vcf using driver mutation set
        add annotated vcf file to concat in final step
        :param muts_file:
        :return:
        """
        muts_outfile = self.outfile_name.format('_muts.vcf.gz')
        cmd = f"bcftools annotate -i '{self.vcf_filter}' --merge-logic DRV:unique" \
              f" -a {muts_file} -h {self.drv_header} " \
              f"-c CHROM,FROM,TO,INFO/DRV {self.vcf_path} |" \
              f"bcftools annotate  -i 'DRV!=\".\" && DRV[*]==VC' | " \
              f"bgzip -c >{muts_outfile} && tabix -f -p vcf {muts_outfile}"
        _run_command(cmd)
        self.merge_vcf_dict['a'] = muts_outfile

    def annotate_lof_genes(self, genes_file, lof_types):
        """
        annotate vcf using known LoF genes
        :param genes_file: lof gene names file
        :param lof_types: lof consequences type string
        :return:
        """
        get_gene = re.compile(r'.*;VD=(\w+)|.*')
        # create dummy genome locationo file to annoate LoF genes...
        genome_loc_file = self.outdir + '/genome.tab.gz'
        create_dummy_genome(self.vcf_path, genome_loc_file)
        genes_outfile = self.outfile_name.format('_genes.vcf')
        lof_outfile = self.outfile_name.format('_genes_lof.vcf')
        lof_gene_list = get_drv_gene_list(genes_file)
        # map lof effect types to pass variants...
        cmd = f"bcftools annotate -a {genome_loc_file} -i '({self.vcf_filter}) && ({lof_types})' " \
              f"-h {self.drv_header} -c CHROM,FROM,TO,INFO/DRV {self.vcf_path} >{genes_outfile}"
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
        self.merge_vcf_dict['b'] = compress_vcf(lof_outfile)

    def concat_results(self):
        concat_drv_out = self.outfile_name.format('_drv.vcf.gz')
        self.merge_vcf_dict['e'] = self.input_data['vcf_file']['path']

        vcf_files = ([self.merge_vcf_dict[myvcf] for myvcf in sorted(self.merge_vcf_dict.keys())])

        cmd = f"bcftools concat --allow-overlaps --rm-dups all {' '.join(vcf_files)} | bcftools sort | " \
              f"bgzip -c >{concat_drv_out} && tabix -f -p vcf {concat_drv_out}"
        _run_command(cmd)
        _run_command('cp -p ' + concat_drv_out + '  ' + concat_drv_out + '.tbi ' + self.outdir + '/..')
        if self.keepTmp:
            _run_command('cp -rp ' + self.outdir + ' ' + self.outdir + '_' + self.vcf_name)


# generic methods ....
def get_drv_gene_list(drv_genes):
    with open(drv_genes) as f_drv:
        lof_gene_list = f_drv.read().splitlines()
    return lof_gene_list


def compress_vcf(vcf):
    """
    :param vcf:
    :return:
    """
    outfile = vcf + '.gz'
    cmd = f"bgzip -c <{vcf} >{outfile} && tabix -f -p vcf {outfile}"
    _run_command(cmd)
    return outfile


def create_dummy_genome(input_vcf, genome_loc):
    """
    :param input_vcf:
    :param genome_loc:
    :return:
    """
    cmd = f"tabix -l {input_vcf} | xargs -I vcf_chr printf 'vcf_chr\t1\t400000000\tLoF\n' |" \
          f"bgzip -c >{genome_loc} && tabix -s1 -b2 -e3 {genome_loc}"
    _run_command(cmd)


def unheader_vcf(header_vcf, unheader_vcf):
    """
    This is used for testing purpose only...
    :param header_vcf:
    :param unheader_vcf:
    :return:
    """
    cmd = f"bcftools view -H {header_vcf} >{unheader_vcf}"
    _run_command(cmd)
    return unheader_vcf


def _run_command(cmd):
    """
    : runs external command
    :param cmd:
    :return: command output
    """
    """ runs command in a shell, returns stdout and exit code"""
    if not len(cmd):
        raise ValueError("Must supply at least one argument")

    # To capture standard error in the result, use stderr=subprocess.STDOUT:
    cmd_obj = Popen(cmd, stdin=None, stdout=PIPE, stderr=PIPE,
                    shell=True, universal_newlines=True, bufsize=-1,
                    close_fds=True, executable='/bin/bash')
    try:
        # To capture standard error in the result, use stderr=subprocess.STDOUT:
        # logging.info("running command:{}".format(cmd))
        (out, error) = cmd_obj.communicate()
        exit_code = cmd_obj.returncode
        if (exit_code == 0):
            logging.debug(f"Command run successfully:\n{cmd}"
                          f"OUT:{out} Error:{error} Exit:{exit_code}\n")
        else:
            logging.info("Error: command exited with non zero exit \
                          status, please check logging file for more details")
            logging.error(f"OUT:{out}:Error:{error}:Exit:{exit_code}")
            if exit_code != 0:
                sys.exit("Exiting...")
        return
    except TimeoutExpired:
        cmd_obj.kill()
        (out, error) = cmd_obj.communicate()
        logging.error(f"Unable to run command:{cmd}: Out:{out} : Error:{error}")
        sys.exit(f"Unable to run command:{cmd}: Out:{out} : Error:{error}")
