from annotate.formatInput import AnnotateVcf
import sys
import os
import argparse
import pkg_resources
import logging.config

# loda config and reference files....
configdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config/')
header_info = configdir + 'vcf_info.txt'
log_config = configdir + 'logging.conf'
genome_tab = configdir + 'genome.tab.gz'
drv_gene = configdir + 'driver_genes.txt'
drv_mut = configdir + 'driver_mutations.tsv.gz'
lof_consq = configdir + 'lof_consequences.txt'
logging.config.fileConfig(log_config)
log = logging.getLogger(__name__)

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
version = pkg_resources.require("annotateVcf")[0].version

def main():
    usage = "\n %prog [options] -vcf test.vcf -mut mut.tsv.gz -gene gene_names.txt -hl vcf_header.txt "

    optParser = argparse.ArgumentParser(prog='annotateVcf',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional = optParser._action_groups.pop()
    required = optParser.add_argument_group('required arguments')

    required.add_argument("-vcf", "--vcf_file", type=str, dest="vcf_file", required=True,
                          default=None, help="vcf_file to annotate")

    optional.add_argument("-dm", "--driver_mutations", type=str, dest="driver_mutations", required=False,
                          default=drv_mut, help='tab separated driver mutation input file  \
                                                (Bgzip-compressed and tabix-indexed file with annotations),  \
                                                see http://samtools.github.io/bcftools/bcftools.html#annotate')

    optional.add_argument("-dg", "--driver_genes", type=str, dest="driver_genes", required=False,
                          default=drv_gene, help=" file containing list of driver genes, one per line")

    optional.add_argument("-gl", "--genomic_loc", type=str, dest="genomic_loc", required=False,
                          default=genome_tab, help=" tab separated file containing genomic locations,  \
                                                   created from genome indedx file")

    optional.add_argument("-hf", "--header_file", type=str, dest="header_file", required=False,
                          default=header_info,
                          help="VCF header INFO file, \
                               refer http://samtools.github.io/bcftools/bcftools.html#annotate")

    optional.add_argument("-lof", "--lof_consequences", type=str, dest="lof_consequences", required=False,
                          default=lof_consq, help="file containing list of consequences to be considered  under  \
                                                 loss of function category")

    optional.add_argument("-o", "--outdir", type=str, dest="outdir",
                          default=None, help="path to output directory")


    optional.add_argument("-v", "--version", action='version', version='%(prog)s ' + version)
    optional.add_argument("-q", "--quiet", action="store_false", dest="verbose", required=False, default=True)

    optParser._action_groups.append(optional)
    if len(sys.argv) == 0:
        optParser.print_help()
        sys.exit(1)
    opts = optParser.parse_args()
    if not opts.vcf_file:
        sys.exit('\nERROR Arguments required\n\tPlease run: annotateVcf --help\n')
    print("Annotating VCF files")
    # vars function returns __dict__ of Namespace instance
    annotater=AnnotateVcf(**vars(opts))
    annotater.run_analysis()
if __name__ == '__main__':
    main()

