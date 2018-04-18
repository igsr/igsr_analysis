package PyHive::PipeConfig::INTEGRATION::VCFIntegrationGATKINDEL;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');  # All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub default_options {
    my ($self) = @_;
    return {
        %{ $self->SUPER::default_options() },               # inherit other stuff from the base class

        'pipeline_name' => 'run_vcfintegration_indels',       # name used by the beekeeper to prefix job names on the farm

        # runnable-specific parameters' defaults:
        'hostname'   => 'mysql-g1kdcc-public',
        'username'     => 'g1krw',
        'port'        => 4197,
        'db' => undef,
        'pwd' => undef,
        'work_dir'    => undef,
        'final_dir' => undef,
	'faix' => undef,
	'newheader' => undef,
	'filelist' => undef, #  List of Bamfiles used for BAM Transposition. If more than one file then the transposition will be done in different runs
	'bedtools_folder' => '/homes/ernesto/bin/bedtools-2.25.0/bin/',
	'bcftools_folder' => '~/bin/bcftools-1.6/',
	'bgzip_folder' => '/nfs/production/reseq-info/work/ernesto/bin/anaconda3/bin/',
	'beagle_folder' => '~/bin/beagle/',
	'caller' => 'UG', # VariantRecalibrator
	'centromeres' => '/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/centromeres_locations.sort.bed', # BED file with centromeres in order to be considered by PyHive.Factories.CoordFactory
	'gatk_folder' => '~/bin/GATK/',
	'ginterval' => undef, # if defined, then do the integration for a certain genomic region
	'gmap_folder' => '/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/GENETIC_MAP/CHROS',
	'vcflib_folder' => '~/bin/vcflib/bin/', # folder containing the vcfallelicprimitives binary
	'makeBGLCHUNKS_folder' => '~/bin/shapeit2_v2_12/bin/makeBGLCHUNKS/bin/',
	'prepareGenFromBeagle4_folder' => '~/bin/shapeit2_v2_12/bin/prepareGenFromBeagle4/bin/',
	'ligateHAPLOTYPES_folder' => '~/bin/shapeit2_v2_12/bin/ligateHAPLOTYPES/bin/',
	'samtools_folder' => '/homes/ernesto/bin/samtools-1.6/bin/',
	'shapeit_folder' => '~/bin/shapeit2_v2_12/bin/',
	'tabix_folder' => '/nfs/production/reseq-info/work/ernesto/bin/anaconda3/bin/',
	'transposebam_folder' => '/homes/ernesto/lib/reseqtrack//c_code/transpose_bam/',
	'tranches' => '[100.0,99.9,99.0,98.0,97.0,96.0,95.0,92.0,90.0,85.0,80.0,75.0,70.0,65.0,60.0,55.0,50.0]', #VariantRecalibrator
        'resources_indels' => '/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/GATK_BUNDLE/resources_indels.json', # VariantRecalibrator
        'indels_annotations' => ['QD','DP','FS','SOR','ReadPosRankSum','MQRankSum','InbreedingCoeff'], #VQSR. annotations for recalibrating indels
	'input_scaffold_prefix' => ['/nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/VARCALL_ALLGENOME_13022017/COMBINING/PRODUCTION/HD_GENOTYPES/OMNI/PHASING/ALL.chip.omni_broad_sanger_combined.20140818.refcorr.biallelic.snps', 
				    '/nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/VARCALL_ALLGENOME_13022017/COMBINING/PRODUCTION/HD_GENOTYPES/AFFY/PHASING/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.ucsc.hg38.refcorr.biallelic.snps'],
# SHAPEIT. Specify here the prefix for the scaffolded microarray genotypes
	'inputthr' => 1.0, # SHAPEIT
	'window' => 0.1, # SHAPEIT
	'states' => 400, # SHAPEIT
	'statesrandom' => 200, # SHAPEIT
	'burn' => 0, # SHAPEIT
	'run' => 12, # SHAPEIT
	'prune' => 4, # SHAPEIT
	'main' => 20, # SHAPEIT
	'window_bglchnks' => 700, # makeBGLCHUNKS
	'overlap_bglchnks' => 200, # makeBGLCHUNKS
	'genome_file' => undef, #PyHive.Factories.CoordFactory. Used to generate the chunks
	'window_coordfactory_4transposebam' =>  undef, #PyHive.Factories.CoordFactory used for the transposebam analysis
	'window_coordfactory_4shapeit' =>  undef, #PyHive.Factories.CoordFactory used for the shapeit analysis
	'offset_coordfactory' => undef, #PyHive.Factories.CoordFactory
	'outprefix' => 'combined.all.chr20', # Prefix used for all output files
	'scaffolded_samples' => undef, #PyHive.VcfIntegration.run_ligateHAPLOTYPES
	'reference' => '/nfs/production/reseq-info/work/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa',
	'store_attributes' => 'False',
        'filelayout' => undef, #file layout for final phased file
	'newlayout' =>  undef, # new file layout for final phased file
	'lsf_queue'   => 'production-rh7'
    };
}

=head2 pipeline_create_commands

    Description : Implements pipeline_create_commands() interface method of Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that lists the commands that will create and set up the Hive database.
                  In addition to the standard creation of the database and populating it with Hive tables and procedures it also creates two pipeline-specific tables used by Runnables to communicate.

=cut

sub pipeline_create_commands {
    my ($self) = @_;
    return [
        @{$self->SUPER::pipeline_create_commands},  # inheriting database and hive tables' creation

        'mkdir -p '.$self->o('work_dir'),
        ];
}


sub resource_classes {
    my ($self) = @_;
    return {
        %{$self->SUPER::resource_classes},
	'500Mb' => { 'LSF' => '-C0 -M512 -q '.$self->o('lsf_queue').' -R"select[mem>512] rusage[mem=512]"' },
        '1Gb' => { 'LSF' => '-C0 -M1024 -q '.$self->o('lsf_queue').' -R"select[mem>1024] rusage[mem=1024]"' },
        '2Gb' => { 'LSF' => '-C0 -M2048 -q '.$self->o('lsf_queue').' -R"select[mem>2048] rusage[mem=2048]"' },
        '5Gb' => { 'LSF' => '-C0 -M5120 -q '.$self->o('lsf_queue').' -R"select[mem>5120] rusage[mem=5120]"' },
        '5Gb5cpus' => { 'LSF' => '-n 5 -C0 -M5120 -q '.$self->o('lsf_queue').' -R"select[mem>5120] rusage[mem=5120]"' },
        '8Gb' => { 'LSF' => '-C0 -M8192 -q '.$self->o('lsf_queue').' -R"select[mem>8192] rusage[mem=8192]"' },
        '10Gb5cpus' => { 'LSF' => '-n 5 -C0 -M10000 -q '.$self->o('lsf_queue').' -R"select[mem>10000] rusage[mem=10000]"' },
        '12Gb4cpus' => { 'LSF' => '-n 4 -C0 -M12288 -q '.$self->o('lsf_queue').' -R"select[mem>12288] rusage[mem=12288]"' },
        '12Gb' => { 'LSF' => '-C0 -M12288 -q '.$self->o('lsf_queue').' -R"select[mem>12288] rusage[mem=12288]"' },
        '15Gb' => { 'LSF' => '-n 20 -C0 -M15360 -q '.$self->o('lsf_queue').' -R"select[mem>15360] rusage[mem=15360]"' },
        '20GbUni' => { 'LSF' => '-C0 -M20000 -q '.$self->o('lsf_queue').' -R"select[mem>20000] rusage[mem=20000]"' },
        '20Gb5cpus' => { 'LSF' => '-n 5 -C0 -M20000 -q '.$self->o('lsf_queue').' -R"select[mem>20000] rusage[mem=20000]"' },
        '20Gb20cpus' => { 'LSF' => '-n 20 -C0 -M20000 -q '.$self->o('lsf_queue').' -R"select[mem>20000] rusage[mem=20000]"' },
        '10cpus' => { 'LSF' => '-n 10 -C0 -M1024 -q '.$self->o('lsf_queue').' -R"select[mem>1024] rusage[mem=1024]"' }
    };
}


sub hive_meta_table {
    my ($self) = @_;
    return {
        %{$self->SUPER::hive_meta_table},       # here we inherit anything from the base class

        'hive_use_param_stack'  => 1,           # switch on the new param_stack mechanism
    };
}

sub pipeline_analyses {
    my ($self) = @_;
    return [

	{   -logic_name => 'find_vcfs_to_combine',
            -module     => 'PyHive.Seed.SeedVCFIntegration',
	    -language   => 'python3',
            -parameters => {
                'filepath'     => '#filepath#'
            },
            -flow_into => {
                1 => ['combine_vcfs']
            },
        },

	{   -logic_name => 'combine_vcfs',
            -module     => 'PyHive.Vcf.VcfCombine',
            -language   => 'python3',
            -parameters => {
                'flist'     => '#flist#',
		'reference' => $self->o('reference'),
		'ginterval' => $self->o('ginterval'),
		'threads' => 4,
		'bcftools_folder' => $self->o('bcftools_folder'),
		'gatk_folder' => $self->o('gatk_folder'),
		'outprefix' => $self->o('outprefix'),
		'work_dir' => $self->o('work_dir')
            },
	    -rc_name => '12Gb4cpus',
	    -flow_into => {
		1 => {
		    'drop_INFO' => {'filepath' => '#out_vcf#'}
		},
	    },
        },

	{   -logic_name => 'drop_INFO',
            -module     => 'PyHive.Vcf.dropInfo',
            -language   => 'python3',
            -parameters => {
                'bcftools_folder' => $self->o('bcftools_folder'),
		'outprefix' => $self->o('outprefix'),
                'work_dir' => $self->o('work_dir')
            },
		    -flow_into => {
			1 => {'splitmultiallelic' => {'filepath' => '#out_vcf#'}}
		}
        },

	{   -logic_name => 'splitmultiallelic',
            -module     => 'PyHive.Vcf.BcftoolsVcfNorm',
            -language   => 'python3',
            -parameters => {
                'bcftools_folder' => $self->o('bcftools_folder'),
		'multiallelics' => 'split',
		'type' => 'both',
                'outprefix' => $self->o('outprefix'),
		'reference' => $self->o('reference'),
                'work_dir' => $self->o('work_dir')
            },
	    -flow_into => {
		1 => {'index_vcf1' => {'filepath' => '#out_vcf#'}}
	    }
        },

	{   -logic_name => 'index_vcf1',
            -module     => 'PyHive.Vcf.VcfIxByTabix',
            -language   => 'python3',
            -parameters => {
                'filepath' => '#filepath#',
                'tabix_folder' => $self->o('tabix_folder'),
                'work_dir' => $self->o('work_dir')
            },
	    -flow_into => {
		1 => {'run_VcfAllelicPrim' => {'filepath' => '#filepath#'}}
	    }
        },

	{   -logic_name => 'run_VcfAllelicPrim',
            -module     => 'PyHive.Vcf.VcfAllelicPrim',
            -language   => 'python3',
            -parameters => {
		'compress' =>1,
		'downstream_pipe' => '~/bin/vt/vt sort - | ~/bin/vt/vt uniq -',
                'filepath' => '#filepath#',
		'bgzip_folder' => $self->o('bgzip_folder'),
		'vcflib_folder' => $self->o('vcflib_folder'),
                'work_dir' => $self->o('work_dir')
            },
		    -flow_into => {
			1 => {'mergemultiallelic' => {'filepath' => '#out_vcf#'}}
		}
        },

	{   -logic_name => 'mergemultiallelic',
            -module     => 'PyHive.Vcf.BcftoolsVcfNorm',
            -language   => 'python3',
            -parameters => {
                'bcftools_folder' => $self->o('bcftools_folder'),
                'multiallelics' => 'merge',
                'type' => 'both',
                'outprefix' => $self->o('outprefix'),
                'reference' => $self->o('reference'),
                'work_dir' => $self->o('work_dir')
            },
	    -flow_into => {
		1 => {'select_indels' => {'filepath' => '#out_vcf#'}}
	    }
        },

	{   -logic_name => 'select_indels',
            -module     => 'PyHive.VcfFilter.SplitVariants',
            -language   => 'python3',
            -parameters => {
                'bcftools_folder' => $self->o('bcftools_folder'),
                'compress' => 'True',
                'type' => 'indels',
                'work_dir' => $self->o('work_dir')
            },
	    -flow_into => {
		1 => {'index_vcf2' => {
		    'filepath' => '#out_vcf#',
		      }
		}
	    }
        },

	{   -logic_name => 'index_vcf2',
            -module     => 'PyHive.Vcf.VcfIxByTabix',
            -language   => 'python3',
            -parameters => {
		'tabix_folder' => $self->o('tabix_folder'),
                'work_dir' => $self->o('work_dir')
            },
	    -flow_into => {
		1 => {'shorten_bamfiles' => {'out_vcf' => '#filepath#'}}
	    }
        },


	{   -logic_name => 'shorten_bamfiles',
            -module     => 'PyHive.File.ShortenFilePaths',
            -language   => 'python3',
            -parameters => {
                'filelist' => $self->o('filelist'),
                'work_dir' => $self->o('work_dir')
            },
	    -flow_into => {
		1 => {'coord_factory' => {
		    'short_flist' => '#short_flist#',
		    'out_vcf' => '#out_vcf#'
		      }}
	    }
        },

	{   -logic_name => 'coord_factory',
	    -module     => 'PyHive.Factories.CoordFactory',
            -language   => 'python3',
            -parameters => {
                'bedtools_folder' => $self->o('bedtools_folder'),
                'genome_file' => $self->o('genome_file'),
		'rextend' => '-1',
#		'ix' => 4,
                'window' => $self->o('window_coordfactory_4transposebam'),
                'verbose' => 1
            },
	    -flow_into => {
		'2->A' => { 'transpose_bam' => {
                    'out_vcf' => '#out_vcf#',
                    'region' => '#chunk#',
		    'ix' => '#ix#',
		    'filelist' => '#short_flist#'
			    }
		},
                'A->1' => [ 'merge_vcf'],
	    },
	},

        {   -logic_name => 'transpose_bam',
            -module     => 'PyHive.Factories.TransposeBam',
            -language   => 'python3',
            -parameters => {
                'filelist' => '#filelist#',
		'region' => '#region#',
                'outprefix' => 'test',
		'transposebam_folder' => $self->o('transposebam_folder'),
                'work_dir' => $self->o('work_dir')
            },
	    -flow_into => {
		1 => {'merge_transpose_bams' => {
		    'bamlist' => '#out_bamlist#',
		    'region' => '#region#',
		      }
		}
	    },
	    -rc_name => '12Gb'
        },

	{   -logic_name => 'merge_transpose_bams',
            -module     => 'PyHive.Bam.RunSamToolsMerge',
            -language   => 'python3',
            -parameters => {
                'bamlist' => '#bamlist#',
                'region' => '#region#',
                'outprefix' => 'test',
                'samtools_folder' => $self->o('samtools_folder'),
                'work_dir' => $self->o('work_dir')
            },
		    -flow_into => {
			1 => {'index_merge_transpose_bams' => {
                    'out_vcf' => '#out_vcf#',
                    'region' => '#region#',
                    'bamfile' => '#merged_bam#'
			      }
			}
		},
            -rc_name => '5Gb'
        },

	{   -logic_name => 'index_merge_transpose_bams',
            -module     => 'PyHive.Bam.RunSamToolsIndex',
            -language   => 'python3',
            -parameters => {
                'bamfile' => '#bamfile#',
                'samtools_folder' => $self->o('samtools_folder'),
            },
	    -flow_into => {
		1 => {'run_gatkug_indels' => {
                    'out_vcf' => '#out_vcf#',
                    'chunk' => '#region#',
                    'bamlist' => '#bamfile#'
		      }
		}
	    },
            -rc_name => '5Gb'
        },

	{   -logic_name => 'run_gatkug_indels',
            -module     => 'PyHive.VariantCalling.GATK_UG',
            -language   => 'python3',
            -parameters => {
		'genotyping_mode' => 'GENOTYPE_GIVEN_ALLELES',
		'glm' => 'INDEL',
		'alleles' => '#out_vcf#',
		'output_mode' => 'EMIT_ALL_SITES',
		'chunk' => '#chunk#',
		'dcov' => 1000,
                'gatk_folder' => $self->o('gatk_folder'),
		'bamlist' => '#bamlist#',
		'bgzip_folder' => $self->o('bgzip_folder'),
		'max_deletion_fraction' => 1.5, #set this parameter to >1 to disable it
		'work_dir' => $self->o('work_dir')."/gatk_ug",
		'reference' => $self->o('reference'),
		'outprefix' => '#out_vcf#',
		'threads' => 1,
                'verbose' => 1
            },
	    -rc_name => '8Gb',
	    -flow_into => {
                1 => [ '?accu_name=allchunks_files&accu_address=[]&accu_input_variable=out_vcf','?accu_name=allixs&accu_address=[]&accu_input_variable=ix']
	    },
        },

	{   -logic_name => 'merge_vcf',
            -module        => 'PyHive.Vcf.VcfConcat',
            -language   => 'python3',
            -parameters    => {
                'outprefix' => "#out_vcf#.merged.vcf.gz",
                'bcftools_folder' => $self->o('bcftools_folder'),
                'verbose' => 'True',
                'work_dir' => $self->o('work_dir')
            },
	    -flow_into => {
		1 => {
		    'index_vcf3' => {
                        'filepath' => '#merged_file#'
		    }
		},
	    },
            -analysis_capacity => 1,
            -rc_name => '500Mb'
        },

	{   -logic_name => 'index_vcf3',
            -module     => 'PyHive.Vcf.VcfIxByTabix',
            -language   => 'python3',
            -parameters => {
                'tabix_folder' => $self->o('tabix_folder'),
                'work_dir' => $self->o('work_dir')
            },
	    -flow_into => {
		1 => {
                    'run_variantrecalibrator_indels' => {
                        'filepath' => '#filepath#'
                    }
                },
	    }
        },

	{   -logic_name => 'run_variantrecalibrator_indels',
            -module        => 'PyHive.VcfFilter.VariantRecalibrator',
            -language   => 'python3',
            -parameters    => {
                'filepath' => '#filepath#',
                'work_dir' => $self->o('work_dir'),
                'caller' => $self->o('caller'),
		'annotations' => $self->o('indels_annotations'),
                'gatk_folder' => $self->o('gatk_folder'),
                'reference' => $self->o('reference'),
                'resources' => $self->o('resources_indels'),
                'tranches' => $self->o('tranches'),
		'intervals' => $self->o('ginterval'),
                'mode' => 'INDEL'
            },
	    -flow_into => {
		1 => {
		    'run_applyrecalibration_indels' => {
                        'filepath' => '#filepath#',
                        'recal_file' => '#recal_f#',
                        'tranches_file' => '#tranches_f#'
		    }},
	    },
            -analysis_capacity => 20,
            -rc_name => '12Gb',
        },

	{   -logic_name => 'run_applyrecalibration_indels',
            -module        => 'PyHive.VcfFilter.ApplyRecalibration',
            -language   => 'python3',
            -parameters    => {
                'filepath' => '#filepath#',
                'work_dir' => $self->o('work_dir'),
                'caller' => $self->o('caller'),
                'gatk_folder' => $self->o('gatk_folder'),
                'bgzip_folder' => $self->o('bgzip_folder'),
                'tabix_folder' => $self->o('tabix_folder'),
                'reference' => $self->o('reference'),
		'ts_filter_level' => 50.0,
                'recal_file' => '#recal_f#',
                'tranches_file' => '#tranches_f#',
                'mode' => 'SNP'
            },
	    -flow_into => {
		1 => {
		    'select_variants' => {
                        'filepath' => '#vcf_filt#'
		    }},
	     },
	     -rc_name => '5Gb',
	},

	{   -logic_name => 'select_variants',
            -module     => 'PyHive.VcfFilter.SelectVariants',
            -language   => 'python3',
            -parameters => {
                'filepath' => '#filepath#',
                'outprefix' => '#filepath#',
                'work_dir' => $self->o('work_dir'),
                'bcftools_folder' => $self->o('bcftools_folder')
            },
            -rc_name => '500Mb',
	    -flow_into => {
		1 => {
		    'convert_pl2gl' => {
                        'filepath' => '#out_vcf#'
		    }},
	   },
        },

	{   -logic_name => 'convert_pl2gl',
            -module     => 'PyHive.Vcf.convertPL2GL',
            -language   => 'python3',
            -parameters => {
		'filepath' => '#filepath#',
		'outprefix' => '#filepath#',
                'work_dir' => $self->o('work_dir'),
                'bcftools_folder' => $self->o('bcftools_folder')
            },
            -rc_name => '500Mb',
	    -flow_into => {
		1 => {
		    'index_vcf4' => {
                        'filepath' => '#out_vcf#'
		    }},
	    },
        },

	{   -logic_name => 'index_vcf4',
            -module     => 'PyHive.Vcf.VcfIxByTabix',
            -language   => 'python3',
            -parameters => {
		'filepath' => '#filepath#', 
                'tabix_folder' => $self->o('tabix_folder'),
                'work_dir' => $self->o('work_dir')
            },
	    -flow_into => {
		1 => {
		    'split_chr' => {
			'filepath' => '#filepath#'
		    }
		},
	    }
        },

	{   -logic_name => 'split_chr',
            -module        => 'PyHive.Factories.SplitVCFintoChros',
            -language   => 'python3',
            -parameters    => {
                'filepath' => '#filepath#',
                'bcftools_folder' => $self->o('bcftools_folder'),
                'faix' => $self->o('faix'),
                'threads' => 10,
		'filt_string' => 'PASS',
                'verbose' => 'True',
                'work_dir' => $self->o('work_dir')
            },
		    -flow_into => {
			2 => { 'rename_chros' => {
			        'filepath' => '#chr#',
				    'ix' => '#ix#',
				'chromname' => '#chromname#'}},
		},
            -analysis_capacity => 1,
            -rc_name => '10cpus'
        },


	{   -logic_name => 'rename_chros',
            -module     => 'PyHive.Vcf.VcfReplaceChrNames',
            -language   => 'python3',
            -parameters => {
                'chr_types' => 'ensembl',
                'work_dir' => $self->o('work_dir')."/#chromname#",
		'bgzip_folder' => $self->o('bgzip_folder')
            },
            -rc_name => '500Mb',
	    -flow_into => {
		1 => {'chunk_factory1' => {
                    'filepath'=> '#vcf_f#'
		      }
		}
	    },
        },


	{   -logic_name => 'chunk_factory1',
            -module     => 'PyHive.Factories.BeagleChunkFactory',
            -language   => 'python3',
            -parameters => {
		'makeBGLCHUNKS_folder' => $self->o('makeBGLCHUNKS_folder'),
		'work_dir' => $self->o('work_dir'),
		'correct' => 1,
		'chro' => '#chromname#',
		'window' => $self->o('window_bglchnks'),
		'overlap' => $self->o('overlap_bglchnks'),
		'verbose' => 1
            },
            -rc_name => '500Mb',
	    -flow_into => {
		'2->A' => {'run_beagle' => {
		    'vcf_file'=> '#filepath#',
		    'region_chunk' => '#chunk#'
			   }
		},
		'A->1' => { 'prepareGen_from_Beagle' => {'vcf_file' => '#filepath#'}}
	    },
	},

	{   -logic_name => 'run_beagle',
            -module     => 'PyHive.VcfIntegration.run_Beagle',
            -language   => 'python3',
            -parameters => {
                'beagle_folder' => $self->o('beagle_folder'),
                'work_dir' => $self->o('work_dir')."/#chromname#/beagle",
		'outprefix' => '#vcf_file#',
		'correct' => 1,
		'nthreads' => 5,
		'verbose' => 1 
            },
	    -flow_into => {
		1 => [ '?accu_name=allbeagle_files&accu_address=[]&accu_input_variable=vcf_f'],
	    },
	    -rc_name => '20Gb5cpus'
        },

	{   -logic_name => 'prepareGen_from_Beagle',
            -module     => 'PyHive.VcfIntegration.run_prepareGenFromBeagle4',
            -language   => 'python3',
            -parameters => {
                'prepareGenFromBeagle4_folder' => $self->o('prepareGenFromBeagle4_folder'),
                'work_dir' => $self->o('work_dir')."/#chromname#",
                'outprefix' => '#vcf_file#.shapeit_input',
		'prefix_in' => '#vcf_file#',
                'verbose' => 1
            },
	    -flow_into => {
		1 => {'chunk_factory2' => {
		    'vcf_file'=> '#vcf_file#',
		    'input_gen'=> '#input_gen#',
                    'input_init' => '#input_init#'		      
		      },
		}
 	    },
            -rc_name => '5Gb'
        },

	{   -logic_name => 'chunk_factory2',
            -module     => 'PyHive.Factories.CoordFactory',
            -language   => 'python3',
            -parameters => {
                'bedtools_folder' => $self->o('bedtools_folder'),
		'genome_file' => $self->o('genome_file'),
                'window' => $self->o('window_coordfactory_4shapeit'),
                'offset' => $self->o('offset_coordfactory'),
		'subtract' => $self->o('centromeres'),
                'verbose' => 1
            },
            -rc_name => '500Mb',
	    -flow_into => {
		'2->A' => {'run_shapeit' => {
                    'gen_gz'=> '#gen_gz#',
		    'gen_sample' => '#gen_sample#',
		    'hap_gz' => '#hap_gz#',
		    'hap_sample' => '#hap_sample#',
                    'chunk' => '#chunk#'
			   }
		},
		'A->1' => { 'run_ligate_haplotypes' => {'vcf_file' => '#vcf_file#'}}
	    },
        },

	{   -logic_name => 'run_shapeit',
            -module     => 'PyHive.VcfIntegration.run_Shapeit',
            -language   => 'python3',
            -parameters => {
                'filepath'     => '#filepath#',
		'gmap_folder' => $self->o('gmap_folder'),
                'shapeit_folder' => $self->o('shapeit_folder'),
		'inputthr' => $self->o('inputthr'),
		'window' => $self->o('window'),
		'states' => $self->o('states'),
		'statesrandom' => $self->o('statesrandom'),
		'burn' => $self->o('burn'),
		'run' => $self->o('run'),
		'prune' => $self->o('prune'),
		'main' => $self->o('main'),
		'outprefix' =>  '#vcf_file#',
		'input_scaffold_prefix' => $self->o('input_scaffold_prefix'),
                'newheader' => $self->o('newheader'),
                'work_dir' => $self->o('work_dir')."/#chr#/shapeit",
                'samplefile' => '#samplefile#'
            },
	    -rc_name => '15Gb',
	    -flow_into => {
                1 => [ '?accu_name=allchunks_files&accu_address=[]&accu_input_variable=hap_gz']
	    },
        },

	{   -logic_name => 'run_ligate_haplotypes',
            -module        => 'PyHive.VcfIntegration.run_ligateHAPLOTYPES',
            -language   => 'python3',
            -parameters    => {
		'hapgz_list' => '#allchunks_files#',
		'vcf_f' => '#vcf_file#',
		'outprefix' => '#vcf_file#.phased',
		'scaffolded_samples' => $self->o('scaffolded_samples'),
                'work_dir' => $self->o('work_dir')."/#chromname#",
                'ligateHAPLOTYPES_folder' => $self->o('ligateHAPLOTYPES_folder'),
                'verbose' => 'True'
            },
            -analysis_capacity => 1,
            -rc_name => '2Gb',
	    -flow_into => {
		1 => {'run_convert_vcf' => {
                    'hap_gz' => '#hap_gz#',
                    'hap_sample' => '#hap_sample#',
		    'vcf_file' => '#vcf_file#'
		      },
		}
	}
        },

	{   -logic_name => 'run_convert_vcf',
            -module        => 'PyHive.VcfIntegration.run_Shapeit_convert2vcf',
            -language   => 'python3',
            -parameters    => {
		'hap_gz' => '#hap_gz#',
		'hap_sample' => '#hap_sample#',
		'compress' => 1,
                'outprefix' => '#vcf_file#.phased',
                'work_dir' => $self->o('work_dir')."/#chromname#",
		'shapeit_folder' => $self->o('shapeit_folder'),
                'verbose' => 'True'
            },
            -analysis_capacity => 1,
            -rc_name => '500Mb',
	    -flow_into => {
		1 => {'store_phased_vcf' => {'filename' => '#out_vcf#'}}
            }
	},

	{   -logic_name => 'store_phased_vcf',
            -module        => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'type' => 'PHASED_VCF',
                'final_dir' => $self->o('final_dir'),
		'filelayout' => $self->o('filelayout'),
                'newlayout' => $self->o('newlayout'),
                'add_date' => 'True',
                'extension' => 'phased.vcf.gz',

            },
        }

	
	];
}

1;

