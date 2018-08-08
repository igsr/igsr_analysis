package PyHive::PipeConfig::INTEGRATION::PHASING.pm;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');  # All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub default_options {
    my ($self) = @_;
    return {
        %{ $self->SUPER::default_options() },               # inherit other stuff from the base class

        'pipeline_name' => 'run_vcfintegration_phasing',       # name used by the beekeeper to prefix job names on the farm

        # runnable-specific parameters' defaults:
        'hostname'   => 'mysql-g1kdcc-public',
        'username'     => 'g1krw',
        'port'        => 4197,
        'db' => undef,
        'pwd' => undef,
        'work_dir'    => undef,
        'final_dir' => undef,
	'log_dir' => undef,
	'faix' => undef,
	'newheader' => undef,
	'filelist' => undef, # List of Bamfiles used for BAM Transposition. If more than one file then the transposition will be done in different runs
	'bedtools_folder' => '/homes/ernesto/bin/bedtools-2.25.0/bin/',
	'bcftools_folder' => '~/bin/bcftools-1.6/',
	'bgzip_folder' => '/nfs/production/reseq-info/work/ernesto/bin/anaconda3/bin/',
	'beagle_folder' => '~/bin/beagle/',
	'beagle_jar' => 'beagle.08Jun17.d8b.jar',
	'caller' => 'UG', # UnifiedGenotyper
	'gmap_folder' => '/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/GENETIC_MAP/CHROS',
	'vcflib_folder' => '~/bin/vcflib/bin/', # folder containing the vcfallelicprimitives binary
	'makeBGLCHUNKS_folder' => '~/bin/shapeit2_v2_12/bin/makeBGLCHUNKS/bin/',
	'prepareGenFromBeagle4_folder' => '~/bin/shapeit2_v2_12/bin/prepareGenFromBeagle4/prepareGenFromBeagle4/bin/',
	'ligateHAPLOTYPES_folder' => '~/bin/shapeit2_v2_12/bin/ligateHAPLOTYPES/bin/',
	'samtools_folder' => '/homes/ernesto/bin/samtools-1.6/bin/',
	'shapeit_folder' => '~/bin/shapeit2_v2_12/bin/',
	'tabix_folder' => '/nfs/production/reseq-info/work/ernesto/bin/anaconda3/bin/',
	'transposebam_folder' => '/homes/ernesto/lib/reseqtrack//c_code/transpose_bam/',
	'tranches' => '[100.0,99.9,99.5,99.2,99.0,98.0,97.0,96.0,95.0,92.0,90.0,85.0,80.0,75.0,70.0,65.0,60.0,55.0,50.0]', #VariantRecalibrator
	'resources_snps' => '/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/GATK_BUNDLE/resources_snps.json', # VariantRecalibrator
	'snps_annotations' => undef, # VQSR. annotations for recalibrating snps
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
	'samplefile' => undef, # SHAPEIT
	'window_bglchnks' => undef, # makeBGLCHUNKS
        'overlap_bglchnks' => undef, # makeBGLCHUNKS
	'window_shapeitchnks' => undef, # makeBGLCHUNKS 4 Shapeit
        'overlap_shapeitchnks' => undef, # makeBGLCHUNKS 4 Shapeit
#	'window_bglchnks' => 12000, # recommended in supp of phase3 for makeBGLCHUNKS
#	'overlap_bglchnks' => 2000, # recommended in supp of phase 3 for makeBGLCHUNKS
	'genome_file' => undef, #PyHive.Factories.CoordFactory. Used to generate the chunks
	'window_coordfactory_4transposebam' =>  undef, #PyHive.Factories.CoordFactory used for the transposebam analysis
	'outprefix' => undef, # Prefix used for all output files
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
	'10Gb10cpus' => { 'LSF' => '-n 10 -C0 -M10000 -q '.$self->o('lsf_queue').' -R"select[mem>10000] rusage[mem=10000]"' },
        '12Gb4cpus' => { 'LSF' => '-n 4 -C0 -M12288 -q '.$self->o('lsf_queue').' -R"select[mem>12288] rusage[mem=12288]"' },
	'12Gb' => { 'LSF' => '-C0 -M12288 -q '.$self->o('lsf_queue').' -R"select[mem>12288] rusage[mem=12288]"' },
	'15Gb' => { 'LSF' => '-n 20 -C0 -M15360 -q '.$self->o('lsf_queue').' -R"select[mem>15360] rusage[mem=15360]"' },
	'20GbUni' => { 'LSF' => '-C0 -M20000 -q '.$self->o('lsf_queue').' -R"select[mem>20000] rusage[mem=20000]"' },
	'20Gb5cpus' => { 'LSF' => '-n 5 -C0 -M20000 -q '.$self->o('lsf_queue').' -R"select[mem>20000] rusage[mem=20000]"' },
	'20Gb10cpus' => { 'LSF' => '-n 10 -C0 -M20000 -q '.$self->o('lsf_queue').' -R"select[mem>20000] rusage[mem=20000]"' },
	'20Gb20cpus' => { 'LSF' => '-n 20 -C0 -M20000 -q '.$self->o('lsf_queue').' -R"select[mem>20000] rusage[mem=20000]"' },
	'75Gb20cpus' => { 'LSF' => '-n 20 -C0 -M75000 -q '.$self->o('lsf_queue').' -R"select[mem>75000] rusage[mem=75000]"' },
	'10cpus' => { 'LSF' => '-n 10 -C0 -M1024 -q '.$self->o('lsf_queue').' -R"select[mem>1024] rusage[mem=1024]"' },
	'20cpus' => { 'LSF' => '-n 20 -C0 -M1024 -q '.$self->o('lsf_queue').' -R"select[mem>1024] rusage[mem=1024]"' }
    };
}


sub hive_meta_table {
    my ($self) = @_;
    return {
        %{$self->SUPER::hive_meta_table},       # here we inherit anything from the base class

        'hive_use_param_stack'  => 0,           # switch on the new param_stack mechanism
    };
}

sub pipeline_analyses {
    my ($self) = @_;
    return [

	{   -logic_name => 'find_vcf_files',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters => {
                'inputcmd'     => 'cat #file#',
                'column_names' => [ 'initial_filename' ],
            },
            -flow_into => {
                2 => ['split_chr']
            },
        },

	{   -logic_name => 'split_chr',
            -module        => 'PyHive.Factories.SplitVCFintoChros',
            -language   => 'python3',
            -parameters    => {
                'filepath' => '#out_vcf#',
                'bcftools_folder' => $self->o('bcftools_folder'),
                'faix' => $self->o('faix'),
                'threads' => 10,
		'filt_string' => 'PASS',
                'verbose' => 'True',
                'work_dir' => $self->o('work_dir')
            },
	    -flow_into => {
		2 => [ 'rename_chros'] 
	    },
            -analysis_capacity => 1,
            -rc_name => '10cpus'
        },

	{   -logic_name => 'rename_chros',
            -module     => 'PyHive.Vcf.VcfReplaceChrNames',
            -language   => 'python3',
            -parameters => {
		'filepath' => '#chr#',
                'chr_types' => 'ensembl',
                'work_dir' => $self->o('work_dir')."/#chromname#",
		'bgzip_folder' => $self->o('bgzip_folder')
            },
            -rc_name => '500Mb',
	    -flow_into => {
		1 => {'chunk_factory1' => INPUT_PLUS() }
	    },
        },


	{   -logic_name => 'chunk_factory1',
            -module     => 'PyHive.Factories.BeagleChunkFactory',
            -language   => 'python3',
            -parameters => {
		'filepath' => '#vcf_f#',
		'makeBGLCHUNKS_folder' => $self->o('makeBGLCHUNKS_folder'),
		'work_dir' => $self->o('work_dir')."/#chromname#/beagle",
		'window' => $self->o('window_bglchnks'),
		'overlap' => $self->o('overlap_bglchnks'),
		'verbose' => 1
            },
            -rc_name => '500Mb',
	    -flow_into => {
		'2->A' => {'run_beagle_lowmem' => {
		    'vcf_file'=> '#filepath#',
		    'region_chunk' => '#chunk#'
			   }
		},
		'A->1' => { 'prepareGen_from_Beagle' => {'vcf_file' => '#filepath#'}}
	    },
	},

	{   -logic_name => 'run_beagle_lowmem',
            -module     => 'PyHive.VcfIntegration.run_Beagle',
            -language   => 'python3',
            -parameters => {
                'beagle_folder' => $self->o('beagle_folder'),
		'beagle_jar' => $self->o('beagle_jar'),
                'work_dir' => $self->o('work_dir')."/#chromname#/beagle",
		'outprefix' => '#vcf_file#',
		'niterations' => 15, #recommended in Supp P3
		'correct' => 1,
		'nthreads' => 5,
		'verbose' => 1 
            },
	    -flow_into => {
	       -1 => [ 'run_beagle_himem' ],
		1 => [ '?accu_name=allbeagle_files&accu_address=[]&accu_input_variable=vcf_f'],
	    },
	    -rc_name => '20Gb5cpus'
        },

	{   -logic_name => 'run_beagle_himem',
            -module     => 'PyHive.VcfIntegration.run_Beagle',
            -language   => 'python3',
            -parameters => {
                'beagle_folder' => $self->o('beagle_folder'),
		'beagle_jar' => $self->o('beagle_jar'),
                'work_dir' => $self->o('work_dir')."/#chromname#/beagle",
                'outprefix' => '#vcf_file#',
                'niterations' => 15, #recommended in Supp P3
                'correct' => 1,
                'nthreads' => 10,
                'verbose' => 1
            },
	    -flow_into => {
                1 => [ '?accu_name=allbeagle_files&accu_address=[]&accu_input_variable=vcf_f'],
	    },
            -rc_name => '20Gb10cpus'
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
		1 => {'chunk_factory2' => INPUT_PLUS() }
 	    },
            -rc_name => '5Gb'
        },

	{   -logic_name => 'chunk_factory2',
            -module     => 'PyHive.Factories.ShapeitChunkFactory',
            -language   => 'python3',
            -parameters => {
		'filepath' => '#vcf_file#',
                'makeBGLCHUNKS_folder' => $self->o('makeBGLCHUNKS_folder'),
                'work_dir' => $self->o('work_dir')."/#chromname#/shapeit/",
                'window' => $self->o('window_shapeitchnks'),
                'overlap' => $self->o('overlap_shapeitchnks'),
                'verbose' => 1
            },
            -rc_name => '500Mb',
	    -flow_into => {
		'2->A' => {'run_shapeit_lowmem' => {
                    'input_gen'=> '#input_gen#',
                    'input_init' => '#input_init#',
                    'chunk' => '#chunk#'
			   }
		},
		'A->1' => { 'run_ligate_haplotypes' => {'vcf_file' => '#vcf_file#'}}
	    },
        },


	{   -logic_name => 'run_shapeit_lowmem',
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
                'work_dir' => $self->o('work_dir')."/#chromname#/shapeit",
		'thread' => 10,
                'samplefile' => $self->o('samplefile')
            },
	    -rc_name => '10Gb10cpus',
	    -flow_into => {
	       -1 => [ 'run_shapeit_himem'],
                1 => [ '?accu_name=allshapeitoutput_files&accu_address=[]&accu_input_variable=hap_gz']
	    },
        },

	{   -logic_name => 'run_shapeit_himem',
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
                'work_dir' => $self->o('work_dir')."/#chromname#/shapeit",
                'thread' => 10,
                'samplefile' => $self->o('samplefile')
            },
            -rc_name => '20Gb10cpus',
	    -flow_into => {
                1 => [ '?accu_name=allshapeitoutput_files&accu_address=[]&accu_input_variable=hap_gz']
	    },
        },

	{   -logic_name => 'run_ligate_haplotypes',
            -module        => 'PyHive.VcfIntegration.run_ligateHAPLOTYPES',
            -language   => 'python3',
            -parameters    => {
		'hapgz_list' => '#allshapeitoutput_files#',
		'vcf_f' => '#vcf_file#',
		'outprefix' => '#vcf_file#.phased',
		'scaffolded_samples' => $self->o('scaffolded_samples'),
                'work_dir' => $self->o('work_dir')."/#chromname#",
                'ligateHAPLOTYPES_folder' => $self->o('ligateHAPLOTYPES_folder'),
                'verbose' => 'True'
            },
            -analysis_capacity => 1,
            -rc_name => '5Gb',
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
            -rc_name => '5Gb',
#	    -flow_into => {
#		1 => ['store_phased_vcf']
#            }
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

