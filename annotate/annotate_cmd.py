from annotate.formatInput import AnnotateVcf
import sys
import os
import argparse
import pkg_resources
import logging.config

configdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config/')
header_info = configdir + 'info.txt'
log_config = configdir + 'logging.conf'
genome_tab = configdir + 'genome.tab.gz'
drv_gene= configdir + 'drv_genes.txt'
drv_mut= configdir + 'drv_mutation.tab.gz'
logging.config.fileConfig(log_config)
log = logging.getLogger(__name__)

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
version = pkg_resources.require("annotateVcf")[0].version

def main():
    usage = "\n %prog [options] -vcf test.vcf -mut mut.tsv.gz -gene gene_names.txt -hl vcf_header.txt "

    optParser = argparse.ArgumentParser(prog='annotateVcf', formatter_class=argparse.ArgumentDefaultsHelpFormatter))
    optional = optParser._action_groups.pop()
    required = optParser.add_argument_group('required arguments')

    required.add_argument("-vcf", "--vcf_file", type=str, dest="vcf_file", required=True,
                          default=None, help="vcf_file to annotate")

    optional.add_argument("-muts", "--driver_mutations", type=str, dest="driver_mutations", required=True,
                          default=drv_mut, help="tab separated driver mutation input file (Bgzip-compressed and tabix-indexed file with annotations), see http://samtools.github.io/bcftools/bcftools.html#annotate")

    optional.add_argument("-genes", "--driver_genes", type=str, dest="driver_genes", required=True,
                          default=drv_gene, help=" file containing list of driver genes, one per line")

    optional.add_argument("-gl", "--genomic_loc", type=str, dest="genomic_loc", required=True,
                          default=genome_tab, help=" tab separated file containing genomic locations, created from genome indedx file")

    optional.add_argument("-ih", "--info_header", type=str, dest="info_header", required=True,
                          default=header_info,
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
    myannotate=AnnotateVcf(**vars(opts))
		myannotate.run_analysis()
if __name__ == '__main__':
    main()

