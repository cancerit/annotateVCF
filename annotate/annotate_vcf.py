import os
import sys
import argparse
import itertools
import logging
from subprocess import Popen, PIPE, STDOUT

logging.basicConfig(filename='example.log',level=logging.DEBUG)

version = "1.0.0"

"""
    This script annotate vcf file using bcftools annotate functionality and custom python code...
    
"""

# print data frame

BCFTOOLS_ALL_PASS_CMD='bcftools view -f PASS {} | bgzip -c >{} && tabix -p vcf {}'

BCFTOOLS_GENE_CMD='bcftools view -f PASS {} | bcftools annotate -a {} -i \'INFO/VC="stop_lost" || INFO/VC="start_lost" ' \
             '|| INFO/VC="ess_splice" || INFO/VC="nonsense" || INFO/VC="frameshift" \' -h {} -c CHROM,FROM,TO,INFO/DRV >{}'

BCFTOOLS_MUT_CMD='bcftools view -f PASS {} | ' \
                 'bcftools annotate -a {} -h {} -c CHROM,FROM,TO,INFO/DRV |' \
                 'bcftools annotate  -i \'INFO/DRV==INFO/VC  && INFO/DRV!="."\' | bgzip -c >{} && tabix -p vcf {}'

BGZIP_TABIX_CMD='cat {} | bgzip -c > {} && tabix -p vcf {}'

BCFTOOLS_CONCAT_VCF='bcftools concat -a --rm-dups all {} {} {} | bgzip -c >{} && tabix -p vcf {}'

def merge_files(**opt):
    """
    :param opt: extension of file containing vaf values, vaf_cutoff
    :return:
    """
    vcf_file = opt['vcf_file'].strip()
    driver_mut_file = opt['driver_mutations'].strip()
    driver_gene_file = opt['driver_genes']
    genome_loc = opt['genomic_loc']
    info_header = opt['info_header']
    with open(driver_gene_file) as f:
        driver_gene_annot_list=f.read().splitlines()
        run_bcftools(vcf_file,driver_mut_file,driver_gene_annot_list,genome_loc,info_header)

    return

def run_bcftools(vcf_file, driver_mut_file, driver_gene_annot_list, genome_loc, info_header):
    """
    :param count_file:
    :return:
    """
    global BCFTOOLS_GENE_CMD
    global BCFTOOLS_MUT_CMD

    (out_file_name_no_ext,file_ext)=_get_file_metadata(vcf_file)
    vcf_out_mut=out_file_name_no_ext+'.mut'+file_ext
    cmd = BCFTOOLS_MUT_CMD.format(vcf_file, driver_mut_file, info_header,
                                  vcf_out_mut, vcf_out_mut)
    run_command(cmd)

    vcf_out_gene = out_file_name_no_ext + '.gene.vcf'
    cmd = BCFTOOLS_GENE_CMD.format(vcf_file, genome_loc, info_header,
                                   vcf_out_gene)

    run_command(cmd)

    vcf_out_gene_filtered = filter_lof_genes(vcf_out_gene, driver_gene_annot_list, out_file_name_no_ext)
    vcf_out_pass=out_file_name_no_ext+'.pass'+file_ext
    cmd=BCFTOOLS_ALL_PASS_CMD.format(vcf_file,vcf_out_pass,vcf_out_pass)

    run_command(cmd)
    # merge all files in order mut, gene and pass
    drv_annot_outfile = out_file_name_no_ext + '.drv' + file_ext
    cmd=BCFTOOLS_CONCAT_VCF.format(vcf_out_mut, vcf_out_gene_filtered, vcf_out_pass, drv_annot_outfile, drv_annot_outfile )

    run_command(cmd)

    return None

def filter_lof_genes(vcf_out_gene, driver_gene_annot_list, out_file_name_no_ext):
    filtered_outfile=out_file_name_no_ext+'.gene.filtered.vcf'
    fh = open(filtered_outfile,"w")
    with open(vcf_out_gene) as gene_f:
        for line in gene_f.readlines():
            #write headwr lines
            if line.startswith('#'):
                fh.write(line)
            else:
                gene=(line.split(';VD=')[1]).split('|')[0]
                #write matching LoF genes with LoF effect....
                if gene in driver_gene_annot_list:
                    fh.write(line)
    fh.close()
    filtered_outfile_gzip=filtered_outfile+'.gz'
    cmd=BGZIP_TABIX_CMD.format(filtered_outfile, filtered_outfile_gzip, filtered_outfile_gzip, )
    run_command(cmd)
    return filtered_outfile_gzip

        
def _get_file_metadata(full_file_name):
    """
      takes file path as input and gives its path and processed extension
      #  If there are two extensions adds second extensions as prefix
    """
    (_, name) = os.path.split(full_file_name)
    (name_no_ext, first_ext) = os.path.splitext(name)
    (name_no_ext2, second_ext) = os.path.splitext(name_no_ext)
    first_ext = second_ext + first_ext
    return name_no_ext2, first_ext

def run_command(cmd):
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
        logging.info("running command:{}".format(cmd))
        (out, error) = cmd_obj.communicate()
        exit_code = cmd_obj.returncode
        if (exit_code == 0):
            logging.info("bcftools run successfully")
        else:
            logging.debug("Error: bcftools exited with non zero exit status, please check log file more details")
            logging.error("OUT:{}:Error:{}:Exit:{}".format(out, error, exit_code))
        return 0
    except OSError as oe:
        logging.error("Unable to run command:{} Error:{}".format(cmd, oe.args[0]))
        sys.exit("Unable to run command:{} Error:{}".format(cmd, oe.args[0]))

def main():
    usage = "\n %prog [options] -vcf test.vcf -mut mut.tsv.gz -gene gene_names.txt -hl vcf_header.txt "

    optParser = argparse.ArgumentParser(prog='jaccard_index')
    optional = optParser._action_groups.pop()
    required = optParser.add_argument_group('required arguments')

    required.add_argument("-vcf", "--vcf_file", type=str, dest="vcf_file", required=True,
                          default=None, help="vcf_file to annotate")

    required.add_argument("-muts", "--driver_mutations", type=str, dest="driver_mutations", required=True,
                          default=None, help="tab separated driver mutation input file (Bgzip-compressed and tabix-indexed file with annotations), see http://samtools.github.io/bcftools/bcftools.html#annotate")

    required.add_argument("-genes", "--driver_genes", type=str, dest="driver_genes", required=True,
                          default=None, help=" file containing list of driver genes, one per line")

    required.add_argument("-gl", "--genomic_loc", type=str, dest="genomic_loc", required=True,
                          default=None, help=" tab separated file containing genomic locations, created from genome indedx file")

    required.add_argument("-ih", "--info_header", type=str, dest="info_header", required=True,
                          default=None,
                          help=" VCF INFO header field file, refer http://samtools.github.io/bcftools/bcftools.html#annotate")

    optional.add_argument("-o", "--outfile", type=str, dest="outfile",
                          default=None, help="path to outfile file, STOUT if not provided")

    optional.add_argument("-v", "--version", action='version', version='%(prog)s ' + version)
    optional.add_argument("-q", "--quiet", action="store_false", dest="verbose", required=False, default=True)

    optParser._action_groups.append(optional)
    if len(sys.argv) == 0:
        optParser.print_help()
        sys.exit(1)
    opts = optParser.parse_args()
    if not opts.vcf_file:
        sys.exit('\nERROR Arguments required\n\tPlease run: annotate_vcf.py --help\n')
    print("Annotating VCF files")
    # vars function returns __dict__ of Namespace instance
    merge_files(**vars(opts))

if __name__ == "__main__":
    main()
