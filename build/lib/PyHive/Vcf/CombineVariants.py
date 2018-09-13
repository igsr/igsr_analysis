import eHive
import subprocess
import os
import sys

class CombineVariants(eHive.BaseRunnable):
    """Combine a snp with an indel VCF"""
  
    def run(self):
        vcf_snps = self.param_required('vcf_snps')
        vcf_indels = self.param_required('vcf_indels')
        
        snps_filename=os.path.split(vcf_snps)[1]
        indels_filename=os.path.split(vcf_indels)[1]

        snps_workdir=os.path.split(vcf_snps)[0]
        indels_workdir=os.path.split(vcf_indels)[0]

        if snps_workdir!=indels_workdir:
             raise Exception("Folders containing the snp and indel VCFs are not the same")

        outfile=snps_workdir+"/"+snps_filename+"."+indels_filename

        command="{0}/bcftools concat {1} {2} -a -o {3} -O z".format(self.param_required('bcftools_folder'),vcf_snps, vcf_indels, outfile)

        if self.param('verbose')=="True":
            self.warning("Combining variants files")
            self.warning("Command used: %s" % command)
 
        try:
            subprocess.check_output(command,shell=True)
        except subprocess.CalledProcessError as e:
            self.warning("Something went wrong while combining the vcfs")
            print(e.output)
            sys.exit(0)

        self.param('out_vcf', outfile)

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 'out_vcf' : self.param('out_vcf') }, 1)




