from abc import ABC, abstractmethod


class AbstractAnnotate(ABC):
    '''
      abstract class inilializes files to be comapred and implements two required methods
      to check user input and load config file
    '''

    def __init__(self, **kwargs):
        self.vcf_file = kwargs['vcf_file']
        self.drv_json = kwargs['driver_json']
        self.drv_data = kwargs['driver_data']
        self.outdir = kwargs['outdir']
        self.keepTmp = kwargs['keepTmp']
        super().__init__()

    @abstractmethod
    def check_input(self):
        pass
