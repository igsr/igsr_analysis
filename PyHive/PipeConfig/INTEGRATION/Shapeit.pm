package PyHive::PipeConfig::INTEGRATION::Shapeit;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf'); # All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub default_options {
    my ($self) = @_;
    return {
        %{$self->SUPER::default_options()}, # inherit other stuff from the base class

        'pipeline_name'      => 'run_shapeit', # name used by the beekeeper to prefix job names on the farm

        # runnable-specific parameters' defaults:
        'hostname'           => 'mysql-g1kdcc-public',
        'username'           => 'g1krw',
        'port'               => 4197,
        'db'                 => undef,
        'pwd'                => undef,
        'work_dir'           => undef,
        'final_dir'          => undef,
        'bcftools_folder'    => '/nfs/software/ensembl/RHEL7-JUL2017-core2/linuxbrew/bin/',
        'shapeit_folder'     => '~/bin/shapeit2_v2_12/bin/',
        'scaffolded_samples' => undef, #PyHive.VcfIntegration.run_ligateHAPLOTYPES
        'tabix_folder'       => '/nfs/software/ensembl/RHEL7-JUL2017-core2/linuxbrew/bin/',
        'newlayout'          => [ 'chip', 'chr' ],
        'lsf_queue'          => 'production-rh7'
    };
}

=head2 pipeline_create_commands

    Description : Implements pipeline_create_commands() interface method of Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that lists the commands that will create and set up the Hive database.
                  In addition to the standard creation of the database and populating it with Hive tables and procedures it also creates two pipeline-specific tables used by Runnables to communicate.

=cut

sub pipeline_create_commands {
    my ($self) = @_;
    return [
        @{$self->SUPER::pipeline_create_commands}, # inheriting database and hive tables' creation

        'mkdir -p ' . $self->o('work_dir'),
    ];
}


sub resource_classes {
    my ($self) = @_;
    return {
        %{$self->SUPER::resource_classes},
        '500Mb'  => { 'LSF' => '-C0 -M512 -q ' . $self->o('lsf_queue') . ' -R"select[mem>512] rusage[mem=512]"' },
        '1Gb'    => { 'LSF' => '-C0 -M1024 -q ' . $self->o('lsf_queue') . ' -R"select[mem>1024] rusage[mem=1024]"' },
        '2Gb'    => { 'LSF' => '-C0 -M2048 -q ' . $self->o('lsf_queue') . ' -R"select[mem>2048] rusage[mem=2048]"' },
        '5Gb'    => { 'LSF' => '-C0 -M5120 -q ' . $self->o('lsf_queue') . ' -R"select[mem>5120] rusage[mem=5120]"' },
        '8Gb'    => { 'LSF' => '-C0 -M8192 -q ' . $self->o('lsf_queue') . ' -R"select[mem>8192] rusage[mem=8192]"' },
        '12Gb'   => { 'LSF' => '-C0 -M12288 -q ' . $self->o('lsf_queue') . ' -R"select[mem>12288] rusage[mem=12288]"' },
        '15Gb'   => { 'LSF' => '-n 20 -C0 -M15360 -q ' . $self->o('lsf_queue') . ' -R"select[mem>15360] rusage[mem=15360]"' },
        '20Gb'   => { 'LSF' => '-n 20 -C0 -M20000 -q ' . $self->o('lsf_queue') . ' -R"select[mem>20000] rusage[mem=20000]"' },
        '10cpus' => { 'LSF' => '-n 10 -C0 -M1024 -q ' . $self->o('lsf_queue') . ' -R"select[mem>1024] rusage[mem=1024]"' },
        '20cpus' => { 'LSF' => '-n 20 -C0 -M5120 -q ' . $self->o('lsf_queue') . ' -R"select[mem>5120] rusage[mem=5120]"' }
    };
}


sub hive_meta_table {
    my ($self) = @_;
    return {
        %{$self->SUPER::hive_meta_table}, # here we inherit anything from the base class

        'hive_use_param_stack' => 1, # switch on the new param_stack mechanism
    };
}

sub pipeline_analyses {
    my ($self) = @_;
    return [

        { -logic_name   => 'seed_pipeline',
            -module     => 'PyHive.Seed.SeedShapeit',
            -language   => 'python3',
            -parameters => {
                'filepath' => '#file#',
                'prefix'   => 'OMNI'

            },
            -flow_into  => {
                1 => { 'run_shapeit' => {
                    'input_bed' => '#input_bed#',
                    'outprefix' => '#outprefix#'
                }
                }
            },
        },

        { -logic_name   => 'run_shapeit',
            -module     => 'PyHive.VcfIntegration.run_Shapeit',
            -language   => 'python3',
            -parameters => {
                'duohmm'         => 1,
                'window'         => 0.5,
                'states'         => 200,
                'burn'           => 10,
                'prune'          => 10,
                'main'           => 50,
                'thread'         => 20,
                'input_bed'      => '#input_bed#',
                'input_gen'      => undef,
                'outprefix'      => '#chro#',
                'shapeit_folder' => $self->o('shapeit_folder'),
                'work_dir'       => $self->o('work_dir'),
                'verbose'        => 1
            },
            -rc_name    => '20cpus',
            -flow_into  => {
                1 => { 'convert2vcf' => {
                    'hap_gz'     => '#hap_gz#',
                    'hap_sample' => '#hap_sample#',
                    'outprefix'  => '#outprefix#'
                }
                }
            },
        },

        { -logic_name   => 'convert2vcf',
            -module     => 'PyHive.VcfIntegration.run_Shapeit_convert2vcf',
            -language   => 'python3',
            -parameters => {
                'hap_gz'         => '#hap_gz#',
                'hap_sample'     => '#hap_sample#',
                'shapeit_folder' => $self->o('shapeit_folder'),
                'compress'       => 1,
                'verbose'        => 1,
                'work_dir'       => $self->o('work_dir'),
                'outprefix'      => '#outprefix#'

            },
            -rc_name    => '500Mb',
            -flow_into  => {
                1 => { 'store_hapgz_file' => {
                    'out_vcf'    => '#out_vcf#',
                    'filename'   => '#hap_gz#',
                    'hap_sample' => '#hap_sample#'
                }
                }
            },
        },

        { -logic_name   => 'store_hapgz_file',
            -module     => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters => {
                'filename'    => '#filename#',
                'hostname'    => $self->o('hostname'),
                'username'    => $self->o('username'),
                'port'        => $self->o('port'),
                'db'          => $self->o('db'),
                'pwd'         => $self->o('pwd'),
                'type'        => 'OMNI_SCAFFOLD_HAPGZ',
                'final_dir'   => $self->o('final_dir'),
                'oldlayout'   => [ 'chip', 'chr', 'extension', 'compression' ],
                'newlayout'   => $self->o('newlayout'),
                'extension'   => 'phased.haps',
                'add_date'    => 'True',
                'store'       => 1,
                'compression' => 'gz'
            },
            -flow_into  => {
                1 => { 'store_hapsample_file' => {
                    'filename' => '#hap_sample#',
                    'out_vcf'  => '#out_vcf#'
                }
                }
            },
        },

        { -logic_name   => 'store_hapsample_file',
            -module     => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters => {
                'filename'  => '#filename#',
                'hostname'  => $self->o('hostname'),
                'username'  => $self->o('username'),
                'port'      => $self->o('port'),
                'db'        => $self->o('db'),
                'pwd'       => $self->o('pwd'),
                'type'      => 'OMNI_SCAFFOLD_HAPSAMPLE',
                'final_dir' => $self->o('final_dir'),
                'oldlayout' => [ 'chip', 'chr', 'extension1', 'extension2' ],
                'newlayout' => $self->o('newlayout'),
                'add_date'  => 'True',
                'store'     => 1,
                'extension' => 'phased.haps.sample',
            },
            -flow_into  => {
                1 => { 'index_vcf' => { 'filepath' => '#out_vcf#' } }
            },
        },

        { -logic_name          => 'index_vcf',
            -module            => 'PyHive.Vcf.VcfIxByTabix',
            -language          => 'python3',
            -parameters        => {
                'filepath'     => '#filepath#',
                'tabix_folder' => $self->o('tabix_folder'),
            },
            -flow_into         => {
                1 => { 'store_vcf_file' => {
                    'filename' => '#filepath#',
                    'vcf_ix'   => '#vcf_ix#'
                }
                }
            },
            -analysis_capacity => 1,
            -rc_name           => '500Mb'
        },

        { -logic_name   => 'store_vcf_file',
            -module     => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters => {
                'filename'    => '#filename#',
                'hostname'    => $self->o('hostname'),
                'username'    => $self->o('username'),
                'port'        => $self->o('port'),
                'db'          => $self->o('db'),
                'pwd'         => $self->o('pwd'),
                'type'        => 'OMNI_SCAFFOLD_VCF',
                'final_dir'   => $self->o('final_dir'),
                'oldlayout'   => [ 'chip', 'chr', 'extension', 'compression' ],
                'newlayout'   => $self->o('newlayout'),
                'add_date'    => 'True',
                'store'       => 1,
                'extension'   => 'phased.vcf',
                'compression' => 'gz'
            },
            -flow_into  => {
                1 => { 'store_vcf_file_ix' => { 'filename' => '#vcf_ix#' } }
            },
        },

        { -logic_name   => 'store_vcf_file_ix',
            -module     => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters => {
                'filename'  => '#filename#',
                'hostname'  => $self->o('hostname'),
                'username'  => $self->o('username'),
                'port'      => $self->o('port'),
                'db'        => $self->o('db'),
                'pwd'       => $self->o('pwd'),
                'type'      => 'OMNI_SCAFFOLD_VCF_IX',
                'final_dir' => $self->o('final_dir'),
                'oldlayout' => [ 'chip', 'chr', 'extension', 'compression', 'extension' ],
                'newlayout' => $self->o('newlayout'),
                'add_date'  => 'True',
                'store'     => 1,
                'extension' => 'phased.vcf.gz.tbi',
            },
        },

    ];
}

1;

