import annotate.io_formatter as formatter
import annotate.vcf_annotator as annotator
import sys
import os
import argparse
import pkg_resources
import logging.config

# load config and reference files....

configdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config/')
log_config = configdir + 'logging.conf'
info_header = configdir + 'info.header'
logging.config.fileConfig(log_config)
log = logging.getLogger(__name__)
version = pkg_resources.require("annotateVcf")[0].version

def main():

    usage = "\n %prog [options] -vcf input.vcf [-drv_json test.json -drv_data test_dir] "

    optParser = argparse.ArgumentParser(prog='annotateVcf',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional = optParser._action_groups.pop()
    required = optParser.add_argument_group('required arguments')

    required.add_argument("-vcf", "--vcf_file", type=str, dest="vcf_file", required=True,
                          default=None, help="vcf_file to annotate")

    required.add_argument("-filter", "--vcf_filter", type=str, dest="vcf_filter", nargs='+', required=False,
                          default=['PASS'], help="Include variant sites matching vcf FILTER flag(s)")

    optional.add_argument("-np", "--normal_panel", type=str, dest="normal_panel", required=False, 
                          default=None, help="normal panle file to flag germline variant sites")

    optional.add_argument("-gt", "--germline_tag", type=str, dest="germline_tag", required=False, 
                          default="NPGL", help="tag to dsiplay in normal panel filtered vcf header and INFO field")

    optional.add_argument("-g", "--lof_genes", type=str, dest="lof_genes", required=False,
                          default=None, help="LoF gene name file to use for Loss of function annotations")

    optional.add_argument("-m", "--mutations", type=str, dest="mutations", required=False,
                          default=None, help="driver mutations file to use for driver variant annotations")

    optional.add_argument("-lof", "--lof_type", type=str, dest="lof_type", nargs='+', metavar='N',
                          required=False, default=["stop_lost","start_lost","ess_splice",
                          "frameshift","nonsense"], help="Loss of function effect type")

    optional.add_argument("-hl", "--header_line", type=str, dest="header_line",
                          required=False, default=info_header, help="Loss of function effect type")

    optional.add_argument("-o", "--outdir", type=str, dest="outdir",
                          default="./", help="path to output directory")

    optional.add_argument("-tmp", "--keepTmp", action="store_true", dest="keepTmp",
                          default=False, help="Flag to keep tmporary files")

    optional.add_argument("-v", "--version", action='version', version='%(prog)s ' + version)
    optional.add_argument("-q", "--quiet", action="store_false", dest="verbose", required=False, default=True)

    optParser._action_groups.append(optional)
    if len(sys.argv) == 0:
        optParser.print_help()
        sys.exit(1)
    opts = optParser.parse_args()
    if not opts.vcf_file:
        sys.exit('\nERROR Arguments required\n\tPlease run: annotateVcf --help\n')
    print("Annotating VCF file")
    
    # vars function returns __dict__ of Namespace instance
    my_formatter = formatter.IO_Formatter(**vars(opts))
    outdir_path=my_formatter.format(['outdir'])
    with formatter.tempdir(outdir_path['outdir']) as base_dir:
      annotator.VcfAnnotator(my_formatter, base_dir)

if __name__ == '__main__':
    main()
