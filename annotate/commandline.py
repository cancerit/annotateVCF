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
filters_json = configdir + 'filters.json'
logging.config.fileConfig(log_config)
log = logging.getLogger(__name__)
version = pkg_resources.require("annotateVcf")[0].version


def main():
    usage = "\n %prog [options] -vcf input.vcf [-filter -np -gt -g -m -lof -hl -o ]"

    optParser = argparse.ArgumentParser(prog='annotateVcf',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional = optParser._action_groups.pop()
    required = optParser.add_argument_group('required arguments')

    required.add_argument("-vcf", "--vcf_file", type=str, dest="vcf_file", required=True,
                          default=None, help="vcf_file to annotate")

    optional.add_argument("-filters", "--vcf_filters", type=str, dest="vcf_filters", required=False,
                          default=filters_json, help="Include vcf filters file in json format \
                          [please refer bcftools documentation for more details \
                           : http://samtools.github.io/bcftools/bcftools.html#expressions]")

    optional.add_argument("-np", "--normal_panel", type=str, dest="normal_panel", required=False,
                          default=None, help="normal panel file to flag germline variant sites")

    optional.add_argument("-g", "--lof_genes", type=str, dest="lof_genes", required=False,
                          default=None, help="Known LoF genes file to annotate LoF variants")

    optional.add_argument("-m", "--mutations", type=str, dest="mutations", required=False,
                          default=None, help="driver mutations file to use for driver variant annotations")

    optional.add_argument("-cpv", "--cancer_predisposition", type=str, dest="cancer_predisposition",
                          required=False, default=None,
                          help="provide file to annotate germline predispostion variants")

    optional.add_argument("-hl", "--header_line", type=str, dest="header_line",
                          required=False, default=info_header, help="vcf info header line and info tag")

    optional.add_argument("-o", "--outdir", type=str, dest="outdir",
                          default="./out_annotatevcf", help="path to output directory")

    optional.add_argument("-tmp", "--keepTmp", action="store_true", dest="keepTmp",
                          default=False, help="Flag to keep temporary files")

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
    outdir_path = my_formatter.format(['outdir'])
    with formatter.tempdir(outdir_path['outdir']) as base_dir:
        annotator.VcfAnnotator(my_formatter, base_dir)


if __name__ == '__main__':
    main()
