package PyHive::PipeConfig::INTEGRATION::VCFIntegration;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');  # All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub default_options {
    my ($self) = @_;
    return {
        %{ $self->SUPER::default_options() },               # inherit other stuff from the base class

        'pipeline_name' => 'run_snptools',       # name used by the beekeeper to prefix job names on the farm

        # runnable-specific parameters' defaults:
        'hostname'   => 'mysql-g1kdcc-public',
        'username'     => 'g1krw',
        'port'        => 4197,
        'db' => undef,
        'pwd' => undef,
        'work_dir'    => undef,
        'final_dir' => undef,
	'faix' => undef,
	'samplefile' => undef,
	'newheader' => undef,
	'bedtools_folder' => '/homes/ernesto/bin/bedtools-2.25.0/bin/',
	'bcftools_folder' => '~/bin/bcftools-1.6/',
	'beagle_folder' => '~/bin/beagle/',
	'gatk_folder' => '~/bin/GATK/',
	'ginterval' => 'chr20', # if defined, then do the integration for a certain genomic region
	'gmap_folder' => '/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/GENETIC_MAP/CHROS',
	'makeBGLCHUNKS_folder' => '~/bin/shapeit2_v2_12/bin/makeBGLCHUNKS/bin/',
	'prepareGenFromBeagle4_folder' => '~/bin/shapeit2_v2_12/bin/prepareGenFromBeagle4/bin/',
	'ligateHAPLOTYPES_folder' => '~/bin/shapeit2_v2_12/bin/ligateHAPLOTYPES/bin/',
	'shapeit_folder' => '~/bin/shapeit2_v2_12/bin/',
	'input_scaffold_prefix' => '/nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/VARCALL_ALLGENOME_13022017/COMBINING/PRODUCTION/HD_GENOTYPES/OMNI/PHASING/ALL.chip.omni_broad_sanger_combined.20140818.refcorr.biallelic.snps', # SHAPEIT. Specify here the prefix for the scaffolded microarray genotypes
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
	'genome_file' => '/nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/VARCALL_ALLGENOME_13022017/COMBINING/DEVEL/INTEGRATION_PIPELINE/chr20.genome', #
	'window_coordfactory' =>  '1400000', #PyHive.Factories.CoordFactory
	'offset_coordfactory' => '200000', #PyHive.Factories.CoordFactory
	'scaffolded_samples' => '', #PyHive.VcfIntegration.run_ligateHAPLOTYPES
	'reference' => '/nfs/production/reseq-info/work/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa',
	'store_attributes' => 'False',
        'filelayout' => [ 'prefix','coords','type1','type2','extension'],
	'newlayout' =>  [ 'prefix','coords','type1','type2'],
	'snptools_folder' => '~/bin/snptools/',
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
        '8Gb' => { 'LSF' => '-C0 -M8192 -q '.$self->o('lsf_queue').' -R"select[mem>8192] rusage[mem=8192]"' },
        '12Gb' => { 'LSF' => '-C0 -M12288 -q '.$self->o('lsf_queue').' -R"select[mem>12288] rusage[mem=12288]"' },
	'15Gb' => { 'LSF' => '-n 20 -C0 -M15360 -q '.$self->o('lsf_queue').' -R"select[mem>15360] rusage[mem=15360]"' },
	'20Gb' => { 'LSF' => '-n 20 -C0 -M20000 -q '.$self->o('lsf_queue').' -R"select[mem>20000] rusage[mem=20000]"' },
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
		'bcftools_folder' => $self->o('bcftools_folder'),
		'gatk_folder' => $self->o('gatk_folder'),
		'outprefix' => 'combined.all.chr20.10e620e6',
		'work_dir' => $self->o('work_dir')
            },
	    -rc_name => '12Gb',
	    -flow_into => {
		1 => {'select_biallelic_snps' => {'filepath' => '#out_vcf#'}}
	    }
	    
        },

	{   -logic_name => 'select_biallelic_snps',
            -module     => 'PyHive.VcfFilter.SplitVariants',
            -language   => 'python3',
            -parameters => {
                'bcftools_folder' => $self->o('bcftools_folder'),
		'biallelic' => 'True',
		'compress' => 'False',
		'type' => 'snps',
                'work_dir' => $self->o('work_dir')
            },
	    -flow_into => {
		1 => ['sample_factory']
	    }
        },

	{   -logic_name => 'sample_factory',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters => {
                'inputcmd'     => 'cat #samplefile#',
                'column_names' => [ 'samplename' ],
            },
            -flow_into => {
                '2->A' => [ 'run_SNPTools_bamodel' ],
		'A->1' => [ 'collect_raw_files' ]
            },
        },

	{   -logic_name => 'run_SNPTools_bamodel',
            -module     => 'PyHive.VcfIntegration.SNPTools_bamodel',
            -language   => 'python3',
            -parameters => {
                'vcf_file'     => '#out_vcf#',
                'snptools_folder' => $self->o('snptools_folder'),
		'sample' => '#samplename#',
		'bamlist' => '#bamlist#',
		'work_dir' => $self->o('work_dir'),
		'verbose' => 1
            },
	    -rc_name => '500Mb',
	    -analysis_capacity => 1000,
	    -flow_into => {
                1 => [ '?accu_name=allraws_files&accu_address=[]&accu_input_variable=raw_f']
	    },
        },

	{   -logic_name => 'collect_raw_files',
            -module     => 'PyHive.VcfIntegration.CollectRawF',
            -language   => 'python3',
            -parameters => {
		'work_dir' => $self->o('work_dir')
            },
            -rc_name => '500Mb',
	    -flow_into => {
                1 => {
		    'run_SNPTools_poprob' => {
			'rawlist'=> '#outfile#'
		    }
		}
	    },
        },

	{   -logic_name => 'run_SNPTools_poprob',
            -module     => 'PyHive.VcfIntegration.SNPTools_poprob',
            -language   => 'python3',
            -parameters => {
		'vcf_file' => '#out_vcf#',
                'snptools_folder' => $self->o('snptools_folder'),
                'rawlist' => '#rawlist#',
		'outprefix' => '#out_vcf#.pop_likelihood',
                'work_dir' => $self->o('work_dir'),
		'verbose' => 1
            },
	    -rc_name => '2Gb',
	    -flow_into => {
		1 => {
		    'chr_factory' => {
                        'probf'=> '#prob_f#'
		    }
		}
	    },
        },

	{   -logic_name => 'chr_factory',
            -module     => 'PyHive.Factories.ChrFactory',
	    -language   => 'python3',
            -parameters => {
		'faix' => $self->o('faix'),
                'column_names' => [ 'chro' ],
            },
	    -flow_into => {
		2 => [ 'run_SNPTools_prob2vcf' ]
            },
        },

	{   -logic_name => 'run_SNPTools_prob2vcf',
            -module     => 'PyHive.VcfIntegration.SNPTools_prob2vcf',
            -language   => 'python3',
            -parameters => {
		'vcf_file' => '#out_vcf#',
                'probf' => '#probf#',
                'snptools_folder' => $self->o('snptools_folder'),
                'outprefix' => '#out_vcf#',
		'chr' => '#chro#',
		'work_dir' => $self->o('work_dir'),
		'verbose' => 1
            },
            -rc_name => '500Mb',
	    -flow_into => {
		1 => {'rename_chros' => {
                    'filepath'=> '#vcf_f#'
		      }
		}
            },
	},

	{   -logic_name => 'rename_chros',
            -module     => 'PyHive.Vcf.VcfReplaceChrNames',
            -language   => 'python3',
            -parameters => {
                'chr_types' => 'ensembl',
                'work_dir' => $self->o('work_dir'),
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
		'chro' => '#chr#',
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
                'work_dir' => $self->o('work_dir'),
		'outprefix' => '#vcf_file#',
		'correct' => 1,
		'nthreads' => 1,
		'verbose' => 1 
            },
	    -flow_into => {
		1 => [ '?accu_name=allbeagle_files&accu_address=[]&accu_input_variable=vcf_f'],
	    },
	    -rc_name => '5Gb'
        },

	{   -logic_name => 'prepareGen_from_Beagle',
            -module     => 'PyHive.VcfIntegration.run_prepareGenFromBeagle4',
            -language   => 'python3',
            -parameters => {
                'prepareGenFromBeagle4_folder' => $self->o('prepareGenFromBeagle4_folder'),
                'work_dir' => $self->o('work_dir'),
                'outprefix' => '#vcf_file#shapeit_input',
		'prefix_in' => '#vcf_file#',
                'verbose' => 1
            },
	    -flow_into => {
		1 => {'chunk_factory2' => {
		    'input_gen'=> '#input_gen#',
                    'input_init' => '#input_init#'		      
		      },
		}
 	    },
            -rc_name => '500Mb'
        },

	{   -logic_name => 'chunk_factory2',
            -module     => 'PyHive.Factories.CoordFactory',
            -language   => 'python3',
            -parameters => {
                'bedtools_folder' => $self->o('bedtools_folder'),
		'genome_file' => $self->o('genome_file'),
                'window' => $self->o('window_coordfactory'),
                'offset' => $self->o('offset_coordfactory'),
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
		'A->1' => { 'run_ligate_haplotypes' => {'vcf_file' => '#filepath#'}}
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
                'work_dir' => $self->o('work_dir'),
                'samplefile' => '#samplefile#'
            },
	    -rc_name => '5Gb',
	    -flow_into => {
                1 => [ '?accu_name=allchunks_files&accu_address=[]&accu_input_variable=hap_gz']
	    },
        },

	{   -logic_name => 'run_ligate_haplotypes',
            -module        => 'PyHive.VcfIntegration.run_ligateHAPLOTYPES',
            -language   => 'python3',
            -parameters    => {
		'hapgz_list' => '#allchunks_files#',
		'vcf_f' => '',
		'scaffolded_samples' => $self->o('scaffolded_samples'),
                'work_dir' => $self->o('work_dir'),
                'ligateHAPLOTYPES_folder' => $self->o('ligateHAPLOTYPES_folder'),
                'verbose' => 'True'
            },
            -analysis_capacity => 1,
            -rc_name => '500Mb'
        }
	];
}

1;

