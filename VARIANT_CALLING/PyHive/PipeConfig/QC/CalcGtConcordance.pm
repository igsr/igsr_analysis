package PyHive::PipeConfig::CalcGtConcordance;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');  # All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;


sub default_options {
    my ($self) = @_;
    return {
        %{ $self->SUPER::default_options() },               # inherit other stuff from the base class

        'pipeline_name' => 'CalcGtConcordance',                   # name used by the beekeeper to prefix job names on the farm

        # runnable-specific parameters' defaults:
        'hostname'   => 'mysql-g1kdcc-public',
        'username'     => 'g1krw',
        'port'        => 4197,
        'db' => undef,
        'pwd' => undef,
        'work_dir'    => undef,
        'final_dir' => undef,
	'picard_folder' => '/homes/ernesto/bin/',
	'truth_vcf' => '/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/GENOTYPES/OMNI/REF_CORRECTED/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.hg38.vcf.ucsc.reheaded.refcorrected.final.bgz',
	'intervals' => '/nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/varcall_chr20_20122016/Low_Cov/CONCORDANCE/TEST/ifile.ival',
	'outprefix' => 'LCvsOMNI',
        'lsf_queue'   => 'production-rh7',
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
        {   -logic_name => 'find_files',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters => {
                'inputcmd'     => 'cat #file#',
                'column_names' => [ 'filename' ],
            },
            -flow_into => {
                2 => ['find_samples']
            },
        },

	{   -logic_name => 'find_samples',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters => {
                'inputcmd'     => 'cat #samples#',
                'column_names' => [ 'sample' ],
            },
            -flow_into => {
                2 => ['run_calc_concordance']
            },
        },
	{   -logic_name => 'run_calc_concordance',
            -module        => 'PyHive.VcfQC.GTConcordance',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
		'sample'=> '#sample#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'final_dir' => $self->o('final_dir'),
		'picard_folder' => $self->o('picard_folder'),
		'truth_vcf' => $self->o('truth_vcf'),
		'intervals' => $self->o('intervals'),
		'outprefix' => $self->o('outprefix'),		
            },
            -analysis_capacity => 20,
            -rc_name => '500Mb',
        },
	];
}

1;
