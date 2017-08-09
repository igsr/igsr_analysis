=pod

=head1 NAME

    Hive::PipeConfig::RunPicardOnExome_conf

=head1 SYNOPSIS

    init_pipeline.pl Hive::PipeConfig::RunPicardOnExome_conf -inputfile file_list -work_dir dir_name -db reseqtrack_db_name -pwd db_pwd


=cut


package Hive::PipeConfig::RunPicardOnExome_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');  # All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly


sub default_options {
    my ($self) = @_;
    return {
        %{ $self->SUPER::default_options() },               # inherit other stuff from the base class

        'pipeline_name' => 'run_picard_on_exomebam',                   # name used by the beekeeper to prefix job names on the farm

        # runnable-specific parameters' defaults:
        'hostname'   => 'mysql-g1kdcc-public',
        'username'     => 'g1krw',
        'port'        => 4197,
	'db' => undef,
	'pwd' => undef,
        'work_dir'    => undef,
	'final_dir' => undef,
	'python_folder' => '/nfs/production/reseq-info/work/ernesto/bin/anaconda2/bin/',
	'script_folder' => '/homes/ernesto/lib/reseq-personal/ernesto/igsr/BamQC/src/',
	'java_folder'    => '/nfs/production/reseq-info/work/bin/java/jdk1.8.0_40/bin/',
        'picard_folder' => '/nfs/production/reseq-info/work/bin/picard-2.7.1/',
	'samtools_folder' => '/nfs/production/reseq-info/work/bin/samtools-1.3/',
	'targetfile' => '/nfs/production/reseq-info/work/zheng/grch38_bams_hsmetrics/g1k_exome_consensus.grch38.sam',
        'lsf_queue'   => 'production-rh6',
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
	  '2Gb' => { 'LSF' => '-C0 -M2048 -q '.$self->o('lsf_queue').' -R"select[mem>2048] rusage[mem=2048]"' },
          '5Gb' => { 'LSF' => '-C0 -M5120 -q '.$self->o('lsf_queue').' -R"select[mem>5120] rusage[mem=5120]"' },
          '8Gb' => { 'LSF' => '-C0 -M8192 -q '.$self->o('lsf_queue').' -R"select[mem>8192] rusage[mem=8192]"' },
          '12Gb' => { 'LSF' => '-C0 -M12288 -q '.$self->o('lsf_queue').' -R"select[mem>12288] rusage[mem=12288]"' }
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
                2  =>  [ 'analyze_bam' ]
            },
        },

        {   -logic_name => 'analyze_bam',
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
		'targetfile' => $self->o('targetfile'),
                'cmd'       => '#python_folder#/python #script_folder#/run_picard_on_exomefile.py --hostname #hostname# --username #username# --port #port# --pwd #pwd# --db #db# --picard #picard_folder# --java #java_folder# --samtools #samtools_folder# --filename #filename# --outdir #final_dir# --targetfile #targetfile# ',
            },
            -analysis_capacity => 30,
            -rc_name => '12Gb'
        },
    ];
}

1;

