package PyHive::PipeConfig::QC::RunBcfToolsStats;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');  # All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub default_options {
    my ($self) = @_;
    return {
        %{ $self->SUPER::default_options() },               # inherit other stuff from the base class

        'pipeline_name' => 'run_filter_freebayes',       # name used by the beekeeper to prefix job names on the farm

        # runnable-specific parameters' defaults:
        'hostname'   => 'mysql-g1kdcc-public',
        'username'     => 'g1krw',
        'port'        => 4197,
        'db' => undef,
        'pwd' => undef,
        'work_dir'    => undef,
        'final_dir' => undef,
	'filelayout' => [ 'dataset','caller','date','sample','status','extension','compression'],
        'newlayout' =>  [ 'dataset','caller','sample','status'],
	'store_attributes' => 'True',
	'bcftools_folder' => '/nfs/software/ensembl/RHEL7-JUL2017-core2/linuxbrew/bin/',
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
        '8Gb' => { 'LSF' => '-C0 -M8192 -q '.$self->o('lsf_queue').' -R"select[mem>8192] rusage[mem=8192]"' }
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
                'column_names' => [ 'initial_filename' ],
            },
            -flow_into => {
                2 => {'chr_factory' => {
		    'filepath' => '#initial_filename#'
		      }
		}
            },
        },

	{   -logic_name => 'chr_factory',
            -module     => 'PyHive.Factories.ChrFactory',
	        -language   => 'python3',
            -parameters => {
		'faix' => $self->o('faix'),
                'column_names' => [ 'chr' ],
            },
	    -flow_into => {
		2 => ['run_bcftools_stats']
            },
        },

	
	{   -logic_name => 'run_bcftools_stats',
            -module        => 'PyHive.VcfQC.BcftoolsStats',
            -language   => 'python3',
            -parameters    => {
                'filepath' => '#filepath#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
		'region' => '#chr#',
		'verbose' => 'True',
		'filter_str' => 'PASS,.',
		'work_dir' => $self->o('work_dir'),
                'store_attributes' => $self->o('store_attributes'),
                'bcftools_folder' => $self->o('bcftools_folder')
            },
	    -analysis_capacity => 1, #important to prevent concurrent access to the same VCF file
            -rc_name => '500Mb'
	    
	}
	];
}

1;

