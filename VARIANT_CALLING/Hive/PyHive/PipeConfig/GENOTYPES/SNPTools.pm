package PyHive::PipeConfig::GENOTYPES::SNPTools;

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
	'bcftools_folder' => '~/bin/bcftools/',
	'beagle_folder' => '~/bin/beagle/',
	'store_attributes' => 'False',
        'filelayout' => [ 'name','coords','status','type','extension'],
	'newlayout' =>  [ 'name','coords','status','type'],
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
	'15Gb' => { 'LSF' => '-C0 -M15360 -q '.$self->o('lsf_queue').' -R"select[mem>15360] rusage[mem=15360]"' },
	'20Gb' => { 'LSF' => '-C0 -M20000 -q '.$self->o('lsf_queue').' -R"select[mem>20000] rusage[mem=20000]"' },
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
                2 => ['split_filename']
            },
        },

        {   -logic_name => 'split_filename',
            -module     => 'PyHive.File.SplitFile',
            -language   => 'python3',
            -parameters => {
                'filename'     => '#initial_filename#',
                'filelayout' => $self->o('filelayout'),
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
            },
            -flow_into => {
                1 => {'sample_factory' => {
                    'filepath'=> '#filepath#',
                    'layout_dict'=> '#layout_dict#'
                      }
                }
            },
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
            -module     => 'PyHive.VcfGenotype.SNPTools_bamodel',
            -language   => 'python3',
            -parameters => {
                'vcf_file'     => '#filepath#',
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
            -module     => 'PyHive.VcfGenotype.CollectRawF',
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
            -module     => 'PyHive.VcfGenotype.SNPTools_poprob',
            -language   => 'python3',
            -parameters => {
		'vcf_file' => '#filepath#',
                'snptools_folder' => $self->o('snptools_folder'),
                'rawlist' => '#rawlist#',
		'outprefix' => '#initial_filename#.pop_likelihood',
                'work_dir' => $self->o('work_dir')
            },
	    -rc_name => '500Mb',
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
		'2->A' => [ 'run_SNPTools_prob2vcf' ],
                'A->1' => [ 'merge_vcf' ]
            },
        },

	{   -logic_name => 'run_SNPTools_prob2vcf',
            -module     => 'PyHive.VcfGenotype.SNPTools_prob2vcf',
            -language   => 'python3',
            -parameters => {
		'vcf_file' => '#filepath#',
                'probf' => '#probf#',
                'snptools_folder' => $self->o('snptools_folder'),
                'outprefix' => '#initial_filename#',
		'chr' => '#chro#',
		'work_dir' => $self->o('work_dir')
            },
            -rc_name => '500Mb',
	    -flow_into => {
		1 => {'run_beagle' => {
                    'filepath'=> '#vcf_f#'
		      }
		}
            },
	},

	{   -logic_name => 'run_beagle',
            -module     => 'PyHive.VcfGenotype.run_Beagle',
            -language   => 'python3',
            -parameters => {
                'vcf_file'     => '#filepath#',
                'beagle_folder' => $self->o('beagle_folder'),
                'work_dir' => $self->o('work_dir'),
		'outprefix' => '#initial_filename#',
		'window' => 12000,
		'overlap' => 2000,
		'niterations' => 15,
		'verbose' => 'True' 
            },
	    -flow_into => {
		1 => {'vcf_reheader' => {
                    'filepath'=> '#vcf_f#'
                      }
                }
	    },
        },

	{   -logic_name => 'vcf_reheader',
            -module     => 'PyHive.Vcf.VcfReheader',
            -language   => 'python3',
            -parameters => {
                'filepath'     => '#filepath#',
                'bcftools_folder' => $self->o('bcftools_folder'),
                'newheader' => $self->o('newheader'),
                'work_dir' => $self->o('work_dir'),
                'samplefile' => '#samplefile#'
            },
	    -flow_into => {
                1 => [ '?accu_name=allchunks_files&accu_address=[]&accu_input_variable=vcf_f',
                       '?accu_name=allixs&accu_address=[]&accu_input_variable=ix']
	    },
        },

	{   -logic_name => 'merge_vcf',
            -module        => 'PyHive.Vcf.VcfConcat',
            -language   => 'python3',
            -parameters    => {
                'outprefix' => $self->o('work_dir')."/#initial_filename#.merged.vcf.gz",
                'bcftools_folder' => $self->o('bcftools_folder'),
                'verbose' => 'True',
                'work_dir' => $self->o('work_dir')
            },
	    -flow_into => {
		1 => {'store_vcf' => {'filename' => '#merged_file#'}}
	    },
            -analysis_capacity => 1,
            -rc_name => '500Mb'
        },

	{   -logic_name => 'store_vcf',
            -module        => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'type' => 'GT_VCF',
                'final_dir' => $self->o('final_dir'),
                'newlayout' => $self->o('newlayout'),
                'add_date' => 'True',
                'extension' => 'beagle.vcf.gz'
            },
        }

	];
}

1;

