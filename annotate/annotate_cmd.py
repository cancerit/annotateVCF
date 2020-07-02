from annotate.formatInput import AnnotateVcf
import sys
import os
import argparse
import pkg_resources
import logging.config

# loda config and reference files....
configdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config/')
log_config = configdir + 'logging.conf'
drv_json = configdir + 'drvData.json'
drv_data = configdir + 'drvData'
logging.config.fileConfig(log_config)
log = logging.getLogger(__name__)

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
version = pkg_resources.require("annotateVcf")[0].version

def main():
    usage = "\n %prog [options] -vcf input.vcf -drv drver.json "

    optParser = argparse.ArgumentParser(prog='annotateVcf',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional = optParser._action_groups.pop()
    required = optParser.add_argument_group('required arguments')

    required.add_argument("-vcf", "--vcf_file", type=str, dest="vcf_file", required=True,
                          default=None, help="vcf_file to annotate")

    optional.add_argument("-drv_json", "--driver_json", type=str, dest="driver_json", required=False,
                          default=drv_json, help="json file containing driver data file names")

    optional.add_argument("-drv_data", "--driver_data", type=str, dest="driver_data", required=False,
                          default=drv_data, help="directory path to driver annoatation files")

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
    print("Annotating VCF files")
    # vars function returns __dict__ of Namespace instance
    annotater=AnnotateVcf(**vars(opts))
    annotater.run_analysis()
if __name__ == '__main__':
    main()

