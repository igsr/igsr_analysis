package PyHive::PipeConfig::INTEGRATION::VCFIntegrationGATKUG;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');  # All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub default_options {
    my ($self) = @_;
    return {
        %{ $self->SUPER::default_options() },               # inherit other stuff from the base class

        'pipeline_name' => 'run_vcfintegration_gatkug',       # name used by the beekeeper to prefix job names on the farm

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
	'variant_type' => undef, # What type of variant is going to be used with this pipeline. Possible values are 'snps' or 'indels'
	'filelist' => undef, # List of Bamfiles used for BAM Transposition. If more than one file then the transposition will be done in different runs
	'bedtools_folder' => '/homes/ernesto/bin/bedtools-2.25.0/bin/',
	'bcftools_folder' => '~/bin/bcftools-1.6/',
	'bgzip_folder' => '/nfs/production/reseq-info/work/ernesto/bin/anaconda3/bin/',
	'caller' => 'UG', # UnifiedGenotyper
	'glm' => undef, # UnifiedGenotyper options
	'dcov' => 250, # UnifiedGenotyper options
	'max_deletion_fraction' => 1.5, # UnifiedGenotyper options
	'centromeres' => '/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/centromeres_and_gaps.bed', # BED file with centromeres and gaps in order to be considered by PyHive.Factories.CoordFactory
	'gatk_folder' => '~/bin/GATK/',
	'java_tmpdir' => '/gpfs/nobackup/resequencing_informatics/ernesto/tmp', # necessary for GATK ApplyRecalibration not to crash 
	'ginterval' => undef, # if defined, then do the integration for a certain genomic region
	'gmap_folder' => '/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/GENETIC_MAP/CHROS',
	'vcflib_folder' => '~/bin/vcflib/bin/', # folder containing the vcfallelicprimitives binary
	'samtools_folder' => '/homes/ernesto/bin/samtools-1.6/bin/',
	'tabix_folder' => '/nfs/production/reseq-info/work/ernesto/bin/anaconda3/bin/',
	'transposebam_folder' => '/homes/ernesto/lib/reseqtrack//c_code/transpose_bam/',
	'tranches' => '[100.0,99.9,99.5,99.2,99.0,98.0,97.0,96.0,95.0,92.0,90.0,85.0,80.0,75.0,70.0,65.0,60.0,55.0,50.0]', #VariantRecalibrator
	'resources' => '/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/GATK_BUNDLE/resources_snps.json', # VariantRecalibrator
	'annotations' => undef, # VQSR. annotations for recalibrating variants
	'mode' => undef, # VQSR mode. Possible values are 'SNP' and 'INDEL'
	'ts_filter_level' => undef, # VQSR ApplyRecalibration 
	'genome_file' => undef, #PyHive.Factories.CoordFactory. Used to generate the chunks
	'outprefix' => undef, # Prefix used for all output files
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

	{   -logic_name => 'find_vcfs_to_combine',
            -module     => 'PyHive.Seed.SeedVCFIntegration',
	    -language   => 'python3',
            -parameters => {
                'filepath'     => '#filepath#'
            },
            -flow_into => {
		'2->A' => [ 'splitmultiallelic' ],
                'A->1' => [ 'combine_vcfs' ],
            },
        },

	{   -logic_name => 'splitmultiallelic',
            -module     => 'PyHive.Vcf.BcftoolsVcfNorm',
            -language   => 'python3',
            -parameters => {
                'filepath' => '#file#',
                'bcftools_folder' => $self->o('bcftools_folder'),
                'multiallelics' => 'split',
                'type' => 'both',
                'outprefix' => $self->o('outprefix'),
                'reference' => $self->o('reference'),
                'work_dir' => $self->o('work_dir')."/normalization"
            },
            -rc_name => '500Mb',
	    -flow_into => {
                1 => { 'index_vcf1' => INPUT_PLUS() }
	    }
        },

	{   -logic_name => 'index_vcf1',
            -module     => 'PyHive.Vcf.VcfIxByTabix',
            -language   => 'python3',
            -parameters => {
                'filepath' => '#out_vcf#',
                'tabix_folder' => $self->o('tabix_folder'),
                'work_dir' => $self->o('work_dir')."/normalization"
            },
            -rc_name => '500Mb',
	    -flow_into => {
		1 => {'run_VcfAllelicPrim'=> INPUT_PLUS() }
	    }
        },

        {   -logic_name => 'run_VcfAllelicPrim',
            -module     => 'PyHive.Vcf.VcfAllelicPrim',
            -language   => 'python3',
            -parameters => {
                'filepath' => '#out_vcf#',
                'compress' =>1,
                'downstream_pipe' => '~/bin/vt/vt sort - | ~/bin/vt/vt uniq -',
                'bgzip_folder' => $self->o('bgzip_folder'),
                'vcflib_folder' => $self->o('vcflib_folder'),
                'work_dir' => $self->o('work_dir')."/normalization"
            },
            -rc_name => '2Gb',
	    -flow_into => {
                1 => ['mergemultiallelic']
	    },
	    -analysis_capacity => 1
        },

        {   -logic_name => 'mergemultiallelic',
            -module     => 'PyHive.Vcf.BcftoolsVcfNorm',
            -language   => 'python3',
            -parameters => {
                'filepath' => '#out_vcf#',
                'bcftools_folder' => $self->o('bcftools_folder'),
                'multiallelics' => 'merge',
                'type' => 'both',
                'outprefix' => $self->o('outprefix'),
                'reference' => $self->o('reference'),
                'work_dir' => $self->o('work_dir')."/normalization"
            },
            -rc_name => '500Mb',
	    -flow_into => {
                1 => ['select_biallelic']
	    }
        },

	{   -logic_name => 'select_biallelic',
            -module     => 'PyHive.VcfFilter.SplitVariants',
            -language   => 'python3',
            -parameters => {
                'filepath' => '#out_vcf#',
                'bcftools_folder' => $self->o('bcftools_folder'),
                'compress' => 'True',
                'type' => $self->o('variant_type'),
                'biallelic' => 'True',
                'work_dir' => $self->o('work_dir')."/normalization"
            },
            -rc_name => '500Mb',
	    -flow_into => {
                1 => ['index_vcf2']
	     }
        },

	{   -logic_name => 'index_vcf2',
            -module     => 'PyHive.Vcf.VcfIxByTabix',
            -language   => 'python3',
            -parameters => {
                'filepath' => '#out_vcf#',
                'tabix_folder' => $self->o('tabix_folder'),
                'work_dir' => $self->o('work_dir')."/normalization"
            },
            -rc_name => '500Mb',
	    -flow_into => {
		1 => [ '?accu_name=allfiles2combine&accu_address=[]&accu_input_variable=out_vcf','?accu_name=alldatasets2combine&accu_address=[]&accu_input_variable=dataset']
	    }
        },

	{   -logic_name => 'combine_vcfs',
            -module     => 'PyHive.Vcf.VcfCombine',
            -language   => 'python3',
            -parameters => {
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
		1 => ['index_vcf3']
	    }
	    
        },

	{   -logic_name => 'index_vcf3',
            -module     => 'PyHive.Vcf.VcfIxByTabix',
            -language   => 'python3',
            -parameters => {
                'filepath' => '#out_vcf#',
                'tabix_folder' => $self->o('tabix_folder'),
                'work_dir' => $self->o('work_dir')
            },
            -rc_name => '500Mb',
	    -flow_into => {
                1 => { 'shorten_bamfiles' => INPUT_PLUS() },
	    }
        },

	{   -logic_name => 'shorten_bamfiles',
            -module     => 'PyHive.File.ShortenFilePaths',
            -language   => 'python3',
            -parameters => {
		'filepath' => '#out_vcf#',
                'filelist' => $self->o('filelist'),
                'work_dir' => $self->o('work_dir')
            },
	    -rc_name => '500Mb',
	    -flow_into => {
		1 => {'coord_factory' => INPUT_PLUS() }
	    }
        },

	{   -logic_name => 'coord_factory',
	    -module     => 'PyHive.Factories.CoordFactory',
            -language   => 'python3',
            -parameters => {
                'bedtools_folder' => $self->o('bedtools_folder'),
                'genome_file' => $self->o('genome_file'),
		'rextend' => '-1',
		'log_dir' => $self->o('log_dir'),
		'chunk_ixs' => '[223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260]',
                'window' => $self->o('window_coordfactory_4transposebam'),
                'verbose' => 1
            },
	    -rc_name => '500Mb',
	    -flow_into => {
		'2->A' => { 'transpose_bams_lowmem' => INPUT_PLUS() },
                'A->1' => [ 'merge_vcf'],
	    },
	},

        {   -logic_name => 'transpose_bams_lowmem',
            -module     => 'PyHive.Factories.TransposeBam',
            -language   => 'python3',
            -parameters => {
                'filelist' => '#filelist#',
		'region' => '#chunk#',
                'outprefix' => 'test',
		'transposebam_folder' => $self->o('transposebam_folder'),
                'work_dir' => $self->o('work_dir')
            },
	    -flow_into => {
	       -1 => {'transpose_bams_himem' => INPUT_PLUS() },
		1 => {'merge_transpose_bams' => INPUT_PLUS() }
	    },
	    -rc_name => '5Gb'
        },

	{   -logic_name => 'transpose_bams_himem',
            -module     => 'PyHive.Factories.TransposeBam',
            -language   => 'python3',
            -parameters => {
                'filelist' => '#filelist#',
                'region' => '#chunk#',
                'outprefix' => 'test',
                'transposebam_folder' => $self->o('transposebam_folder'),
                'work_dir' => $self->o('work_dir')
            },
	    -flow_into => {
		1 => {'merge_transpose_bams' => INPUT_PLUS() }
	    },
            -rc_name => '12Gb'
        },

	{   -logic_name => 'merge_transpose_bams',
            -module     => 'PyHive.Bam.RunSamToolsMerge',
            -language   => 'python3',
            -parameters => {
                'bamlist' => '#out_bamlist#',
                'region' => '#chunk#',
                'outprefix' => 'test',
                'samtools_folder' => $self->o('samtools_folder'),
                'work_dir' => $self->o('work_dir')
            },
	    -flow_into => {
		1 => {'index_merge_transpose_bams' => INPUT_PLUS() }
	    },
            -rc_name => '5Gb'
        },

	{   -logic_name => 'index_merge_transpose_bams',
            -module     => 'PyHive.Bam.RunSamToolsIndex',
            -language   => 'python3',
            -parameters => {
                'bamfile' => '#merged_bam#',
                'samtools_folder' => $self->o('samtools_folder'),
            },
	    -flow_into => {
		1 => {'run_gatkug_lowmem' => INPUT_PLUS() }
	    },
            -rc_name => '5Gb'
        },

	{   -logic_name => 'run_gatkug_lowmem',
            -module     => 'PyHive.VariantCalling.GATK_UG',
            -language   => 'python3',
            -parameters => {
		'genotyping_mode' => 'GENOTYPE_GIVEN_ALLELES',
		'glm' => $self->o('glm'),
		'alleles' => '#out_vcf#',
		'output_mode' => 'EMIT_ALL_SITES',
		'chunk' => '#chunk#',
		'dcov' => $self->o('dcov'),
                'gatk_folder' => $self->o('gatk_folder'),
		'bamlist' => '#merged_bam#',
		'bgzip_folder' => $self->o('bgzip_folder'),
		'work_dir' => $self->o('work_dir')."/gatk_ug",
		'log_file' => $self->o('work_dir')."/gatk_ug/gatk_ug",
		'max_deletion_fraction' => $self->o('max_deletion_fraction'), 
		'reference' => $self->o('reference'),
		'outprefix' => '#out_vcf#',
		'threads' => 1,
                'verbose' => 1
            },
	    -rc_name => '2Gb',
	    -flow_into => {
	       -1 => [ 'run_gatkug_himem1'], 
                1 => [ '?accu_name=allchunks_files&accu_address=[]&accu_input_variable=out_vcf','?accu_name=allixs&accu_address=[]&accu_input_variable=ix']
	    },
        },

	{   -logic_name => 'run_gatkug_himem1',
            -module     => 'PyHive.VariantCalling.GATK_UG',
            -language   => 'python3',
            -parameters => {
                'genotyping_mode' => 'GENOTYPE_GIVEN_ALLELES',
                'glm' => $self->o('glm'),
                'alleles' => '#out_vcf#',
                'output_mode' => 'EMIT_ALL_SITES',
                'chunk' => '#chunk#',
                'dcov' => $self->o('dcov'),
                'gatk_folder' => $self->o('gatk_folder'),
                'bamlist' => '#merged_bam#',
                'bgzip_folder' => $self->o('bgzip_folder'),
                'work_dir' => $self->o('work_dir')."/gatk_ug",
		'log_file' => $self->o('work_dir')."/gatk_ug/gatk_ug",
                'max_deletion_fraction' => $self->o('max_deletion_fraction'),
                'reference' => $self->o('reference'),
                'outprefix' => '#out_vcf#',
                'threads' => 1,
                'verbose' => 1
            },
            -rc_name => '5Gb',
	    -flow_into => {
	       -1 => [ 'run_gatkug_himem2' ],
                1 => [ '?accu_name=allchunks_files&accu_address=[]&accu_input_variable=out_vcf','?accu_name=allixs&accu_address=[]&accu_input_variable=ix']
	   },
        },
	
	{   -logic_name => 'run_gatkug_himem2',
            -module     => 'PyHive.VariantCalling.GATK_UG',
            -language   => 'python3',
            -parameters => {
                'genotyping_mode' => 'GENOTYPE_GIVEN_ALLELES',
                'glm' => $self->o('glm'),
                'alleles' => '#out_vcf#',
                'output_mode' => 'EMIT_ALL_SITES',
                'chunk' => '#chunk#',
                'dcov' => $self->o('dcov'),
                'gatk_folder' => $self->o('gatk_folder'),
                'bamlist' => '#merged_bam#',
                'bgzip_folder' => $self->o('bgzip_folder'),
                'work_dir' => $self->o('work_dir')."/gatk_ug",
		'log_file' => $self->o('work_dir')."/gatk_ug/gatk_ug",
                'max_deletion_fraction' => $self->o('max_deletion_fraction'),
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
		'log_dir' => $self->o('log_dir'),
		'threads' => 5,
                'work_dir' => $self->o('work_dir')
            },
	    -flow_into => {
		1 => ['index_vcf4']
	    },
            -analysis_capacity => 1,
            -rc_name => '500Mb'
        },

	{   -logic_name => 'index_vcf4',
            -module     => 'PyHive.Vcf.VcfIxByTabix',
            -language   => 'python3',
            -parameters => {
		'filepath' => '#merged_file#',
                'tabix_folder' => $self->o('tabix_folder'),
                'work_dir' => $self->o('work_dir')
            },
	    -rc_name => '500Mb',
	    -flow_into => {
		1 => {'run_variantrecalibrator' => INPUT_PLUS() }
	    }
        },

	{   -logic_name => 'run_variantrecalibrator',
            -module        => 'PyHive.VcfFilter.VariantRecalibrator',
            -language   => 'python3',
            -parameters    => {
                'filepath' => '#merged_file#',
                'work_dir' => $self->o('work_dir'),
		'log_file' => $self->o('work_dir')."/gatk_variantrecalibratior",
                'caller' => $self->o('caller'),
		'annotations' => $self->o('annotations'),
                'gatk_folder' => $self->o('gatk_folder'),
                'reference' => $self->o('reference'),
                'resources' => $self->o('resources'),
                'tranches' => $self->o('tranches'),
                'intervals' => $self->o('ginterval'),
		'verbose' => 'True',
                'mode' => $self->o('mode')
            },
	    -flow_into => {
		1 => {'run_applyrecalibration' => INPUT_PLUS() }
	    },
            -analysis_capacity => 1,
            -rc_name => '12Gb',
        },

	{   -logic_name => 'run_applyrecalibration',
            -module        => 'PyHive.VcfFilter.ApplyRecalibration',
            -language   => 'python3',
            -parameters    => {
                'filepath' => '#merged_file#',
                'work_dir' => $self->o('work_dir'),
		'log_file' => $self->o('work_dir')."/gatk_applyrecalibration",
		'tmp_dir' => $self->o('java_tmpdir'),
                'caller' => $self->o('caller'),
                'gatk_folder' => $self->o('gatk_folder'),
                'bgzip_folder' => $self->o('bgzip_folder'),
                'tabix_folder' => $self->o('tabix_folder'),
                'reference' => $self->o('reference'),
                'recal_file' => '#recal_f#',
		'threads' => 1,
		'ts_filter_level' => $self->o('ts_filter_level'),
                'tranches_file' => '#tranches_f#',
                'mode' => $self->o('mode')
            },
	    -flow_into => {
		1 => ['select_variants']
	    },
	    -analysis_capacity => 1,
	    -rc_name => '5Gb',
	},

        {   -logic_name => 'select_variants',
            -module     => 'PyHive.VcfFilter.SelectVariants',
            -language   => 'python3',
            -parameters => {
                'filepath' => '#vcf_filt#',
		'uncalled' => 'exclude', #it is necessary to exclude sites with missing genotypes in order for Beagle not to crash
                'outprefix' => '#vcf_filt#',
		'threads' => 20,
                'work_dir' => $self->o('work_dir'),
                'bcftools_folder' => $self->o('bcftools_folder')
            },
            -rc_name => '20cpus',
	    -flow_into => {
		1 => ['convert_pl2gl']
	    },
        },

	{   -logic_name => 'convert_pl2gl',
            -module     => 'PyHive.Vcf.convertPL2GL',
            -language   => 'python3',
            -parameters => {
		'filepath' => '#out_vcf#',
		'outprefix' => '#out_vcf#',
		'threads' => 1,
                'work_dir' => $self->o('work_dir'),
                'bcftools_folder' => $self->o('bcftools_folder')
            },
            -rc_name => '500Mb',
	    -flow_into => {
		1 => ['index_vcf5']
	    },
        },

	{   -logic_name => 'index_vcf5',
            -module     => 'PyHive.Vcf.VcfIxByTabix',
            -language   => 'python3',
            -parameters => {
		'filepath' => '#out_vcf#', 
                'tabix_folder' => $self->o('tabix_folder'),
                'work_dir' => $self->o('work_dir')
            },
	    -rc_name => '500Mb',
        }
	];
}

1;

