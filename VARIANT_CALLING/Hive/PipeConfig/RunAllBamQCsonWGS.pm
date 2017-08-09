=pod

=head1 NAME

    Hive::PipeConfig::RunAllQCsonWGS

=head1 SYNOPSIS

    init_pipeline.pl Hive::PipeConfig::RunAllQCsonWGS -inputfile file_list -work_dir dir_name -db reseqtrack_db_name -pwd db_pwd

    Run all QC tests on a WGS file. chkindel_rg, verifybamid, picard's CollectWgsMetrics

=cut

package Hive::PipeConfig::RunAllQCsonWGS;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');  # All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub default_options {
    my ($self) = @_;
    return {
        %{ $self->SUPER::default_options() },               # inherit other stuff from the base class

        'pipeline_name' => 'run_all_bamqcs_on_wgs',                   # name used by the beekeeper to prefix job names on the farm

        # runnable-specific parameters' defaults:
        'hostname'   => 'mysql-g1kdcc-public',
        'username'     => 'g1krw',
        'port'        => 4197,
        'db' => undef,
        'pwd' => undef,
        'work_dir'    => undef,
        'final_dir' => undef,
        'python_folder' => '/nfs/production/reseq-info/work/ernesto/bin/anaconda3/bin/',
        'script_folder' => '/homes/ernesto/lib/reseq-personal/ernesto/igsr/BamQC/src/',
        'chk_indel_rg_folder' => '/homes/ernesto/bin/',
	'verifybamid_folder' => '/nfs/production/reseq-info/work/bin/verifyBamID_1.1.3/verifyBamID/bin/',
	'java_folder'    => '/nfs/production/reseq-info/work/bin/java/jdk1.8.0_40/bin/',
        'picard_folder' => '/nfs/production/reseq-info/work/bin/picard-2.7.1/',
	'reference' => '/nfs/production/reseq-info/work/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa',
        'samtools_folder' => '/nfs/production/reseq-info/work/bin/samtools-1.3/',
        'lsf_queue'   => 'production-rh7',
    };
}


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
	%{$self->SUPER::hive_meta_table},
    'hive_use_param_stack' => 1,
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
		2 => WHEN(
                    '#run_chk_indel_rg#==1' => ['run_chk_indel_rg'],
                    '#run_verifybamid#==1' => ['run_verifybamid'],
                    '#run_picard_on_wgsfile#==1' => ['run_picard_on_wgsfile'],
                    ),
            },
        },

        {   -logic_name => 'run_chk_indel_rg',
            -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters    => {
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'python_folder' => $self->o('python_folder'),
                'script_folder' => $self->o('script_folder'),
                'chk_indel_rg_folder' => $self->o('chk_indel_rg_folder'),
                'final_dir' => $self->o('final_dir'),
                'cmd'       => '#python_folder#/python #script_folder#/run_chk_indel_rg.py --hostname #hostname# --username #username# --port #port# --pwd #pwd# --db #db# --exe #chk_indel_rg_folder# --filename #filename# --outdir #final_dir#',
            },
	    -analysis_capacity => 20,
            -rc_name => '500Mb',
        },

	{   -logic_name => 'run_verifybamid',
            -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters    => {
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'genotype_folder' => $self->o('genotype_folder'),
                'prefix' => $self->o('prefix'),
                'python_folder' => $self->o('python_folder'),
                'script_folder' => $self->o('script_folder'),
                'verifybamid_folder' => $self->o('verifybamid_folder'),
                'final_dir' => $self->o('final_dir'),
                'cmd'       => '#python_folder#/python #script_folder#/run_verifybamid.py --hostname #hostname# --username #username# --port #port# --pwd #pwd# --db #db# --exe #verifybamid_folder# --filename #filename# --outdir #final_dir# --genotypes #genotype_folder# --prefix #prefix#',
            },
	    -analysis_capacity => 20,
            -rc_name => '1Gb',
        },

	{   -logic_name => 'run_picard_on_wgsfile',
            -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters    => {
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'java_folder' => $self->o('java_folder'),
                'python_folder' => $self->o('python_folder'),
                'script_folder' => $self->o('script_folder'),
                'picard_folder' => $self->o('picard_folder'),
                'samtools_folder' => $self->o('samtools_folder'),
                'final_dir' => $self->o('final_dir'),
                'reference' => $self->o('reference'),
                'cmd'       => '#python_folder#/python #script_folder#/run_picard_on_WGSfile.py --hostname #hostname# --username #username# --port #port# --pwd #pwd# --db #db# --picard #picard_folder# --java #java_folder# --samtools #samtools_folder# --filename #filename# --outdir #final_dir# --reference #reference# ',
            },
            -analysis_capacity => 20,
            -rc_name => '12Gb'
        },

	];
}

1;
