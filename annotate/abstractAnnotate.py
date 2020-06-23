from abc import ABC, abstractmethod


class AbstractAnnoate(ABC):
    '''
      abstract class inilializes files to be comapred and implements two required methods
      to check user input and load config file
    '''

    def __init__(self, **kwargs):
        self.vcffile = kwargs['vcf_file']
        self.drv_mut = kwargs['driver_mutations']
        self.drv_genes = kwargs.get('driver_genes')
        self.genome_tab = kwargs.get('genomic_loc')
        self.outdir = kwargs.get('outdir')
        self.info_header = kwargs.get('info_header')
        super().__init__()

    @abstractmethod
    def check_input(self):
        pass
