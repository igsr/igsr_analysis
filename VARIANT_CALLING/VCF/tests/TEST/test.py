from VCFIntegration.Shapeit import Shapeit


shapeit_o=Shapeit(ligateHAPLOTYPES_folder = '~/bin/shapeit2_v2_12/bin/ligateHAPLOTYPES/bin/')
shapeit_o.ligate_shapeitchunks(vcf_f="/nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/VARCALL_ALLGENOME_13022017/COMBINING/DEVEL/EXAMPLEFROMWEB/DATA/example/GLs.vcf.gz",
                               scaffolded_samples="/nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/VARCALL_ALLGENOME_13022017/COMBINING/DEVEL/EXAMPLEFROMWEB/only_3.txt",
                               chunk_str="/nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/VARCALL_ALLGENOME_13022017/COMBINING/DEVEL/EXAMPLEFROMWEB/output.shapeit.22.20000000.20060000.haps.gz /nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/VARCALL_ALLGENOME_13022017/COMBINING/DEVEL/EXAMPLEFROMWEB/output.shapeit.22.20050000.20100000.haps.gz",
                               output_prefix="output.shapeit.22.ligated")


