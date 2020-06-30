from abc import ABC, abstractmethod


class AbstractAnnotate(ABC):
    '''
      abstract class inilializes files to be comapred and implements two required methods
      to check user input and load config file
    '''

    def __init__(self, **kwargs):
        self.vcf_file = kwargs['vcf_file']
        self.drv_muts = kwargs['driver_mutations']
        self.drv_genes = kwargs.get('driver_genes')
        self.genome_loc = kwargs.get('genomic_loc')
        self.lof_consequences = kwargs.get('lof_consequences')
        self.header_file = kwargs.get('header_file')
        self.outdir = kwargs.get('outdir')
        self.keepTmp = kwargs.get('keepTmp')
        super().__init__()

    @abstractmethod
    def check_input(self):
        pass
