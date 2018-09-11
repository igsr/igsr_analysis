=pod

=head1 NAME

    PyHive::PipeConfig::QC::RunBamQCsonWGS

=head1 SYNOPSIS

    Run all QC tests on a WGS file. chkindel_rg, verifybamid, picard's CollectWgsMetrics

Below in the 'default_options' function we can see the options that will control the behaviour
of the pipeline. The options that do not have a default value must be set when initializing the
pipeline using 'init_pipeline.pl'. Here are explanations for some of these options, modify them according
to your needs:

 -hostname, username, port, db, pwd control the connection details for the ReseqTrack database
 -work_dir: folder that will be used to put the intermediate files

 -store_attributes: Possible values are 'True'/'False'. If 'True' then the pipeline will store the QC stats
  on the BAM calculated by all programs used by this pipeline in the 'Atttribute' table of the ReseqTrack database
 -genotype_folder: Folder containing the chip high density microarray VCF that will be used by VerifyBAMID

=cut

package PyHive::PipeConfig::QC::RunBamQCsonWGS;

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
        'python_folder' => '/nfs/production/reseq-info/work/ernesto/bin/anaconda3/bin/',
        'script_folder' => '/homes/ernesto/lib/reseq-personal/ernesto/igsr/BamQC/src/',
	'store_attributes' => undef,
        'chk_indel_rg_folder' => '/homes/ernesto/bin/',
	'verifybamid_folder' => '/nfs/production/reseq-info/work/bin/verifyBamID_1.1.3/verifyBamID/bin/',
	'genotype_folder' => undef,
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
                'column_names' => [ 'filepath' ],
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
            -module        => 'PyHive.BamQC.RunChkIndelRg',
	    -language   => 'python3',
            -parameters    => {
		'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'python_folder' => $self->o('python_folder'),
                'chk_indel_rg_folder' => $self->o('chk_indel_rg_folder'),
                'work_dir' => $self->o('work_dir'),
            },
	    -analysis_capacity => 20,
            -rc_name => '500Mb',
	    -flow_into => {
		1 => ['store_chkindelrg_attribute']
	    },
        },

	{   -logic_name => 'store_chkindelrg_attribute',
            -module        => 'PyHive.Attribute.StoreAttribute',
            -language   => 'python3',
            -parameters    => {
		'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
		'store_attributes' => $self->o('store_attributes')
            },
        },

	{   -logic_name => 'run_verifybamid',
            -module        => 'PyHive.BamQC.RunVerifyBamId',
	    -language   => 'python3',
            -parameters    => {
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'genotype_folder' => $self->o('genotype_folder'),
                'python_folder' => $self->o('python_folder'),
                'verifybamid_folder' => $self->o('verifybamid_folder'),
                'work_dir' => $self->o('work_dir')
            },
	    -analysis_capacity => 1,
            -rc_name => '1Gb',
	    -flow_into => {
                1 => ['store_vfbamid_attribute']
	    },
        },

	{   -logic_name => 'store_vfbamid_attribute',
            -module        => 'PyHive.Attribute.StoreAttribute',
            -language   => 'python3',
            -parameters    => {
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'store_attributes' => $self->o('store_attributes')
            },
        },
	
	{   -logic_name => 'run_picard_on_wgsfile',
            -module        => 'PyHive.BamQC.RunPicardOnWGSfile',
	    -language   => 'python3',
            -parameters    => {
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'java_folder' => $self->o('java_folder'),
                'python_folder' => $self->o('python_folder'),
                'picard_folder' => $self->o('picard_folder'),
                'work_dir' => $self->o('work_dir'),
                'reference' => $self->o('reference'),
            },
            -analysis_capacity => 20,
            -rc_name => '12Gb',
	    -flow_into => {
                1 => ['store_picard_attribute']
            },
        },

	{   -logic_name => 'store_picard_attribute',
            -module        => 'PyHive.Attribute.StoreAttribute',
            -language   => 'python3',
            -parameters    => {
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'store_attributes' => $self->o('store_attributes')
            },
        }

	];
}

1;
=cut
