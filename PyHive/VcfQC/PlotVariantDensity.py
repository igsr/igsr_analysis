import eHive
import os
import tempfile
from VcfQC import VcfQC

class PlotVariantDensity(eHive.BaseRunnable):
    """Plot Variant density using a VCF file"""
    
    def param_defaults(self):
        return {
        }

    def run(self):

        filepath=self.param_required('filepath')

        self.warning('Analysing file: %s'% filepath)

        vcfQC = VcfQC(vcf=filepath,bedtools_folder=self.param_required('bedtools_folder'),r_folder=self.param_required('r_folder'),r_scripts=self.param_required('r_scripts'))

        outprefix=""
        path=os.path.split(filepath)[0]
        file=os.path.split(filepath)[1]
        if self.param_is_defined('work_dir'):
            outprefix= tempfile.NamedTemporaryFile(dir=self.param('work_dir'),prefix=file)
        else:
            outprefix= tempfile.NamedTemporaryFile(dir=path,prefix=file)

        #getting the parameters for plot
        plot_params= {
            'height': self.param_required('height'),
            'chromName_cex': self.param_required('chromName_cex')
        }
        density_plot_f=vcfQC.plot_variant_density(length=self.param_required('length'), 
                                                  genome= self.param_required('genome'), 
                                                  outprefix=outprefix.name, 
                                                  plot_params= plot_params)

        self.param('density_plot_f', density_plot_f)

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 'density_plot_f' : self.param('density_plot_f') }, 1)
