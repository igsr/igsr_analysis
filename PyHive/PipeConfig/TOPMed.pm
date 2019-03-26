package PyHive::PipeConfig::TOPMed;
use strict;
use warnings;

# All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly
use base ('ReseqTrack::Hive::PipeConfig::ReseqTrackGeneric_conf');
#use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;



sub default_options {
    my ($self) = @_;
    return {

        # inherit other stuff from the base class
        %{$self->SUPER::default_options()},

        # name used by the beekeeper to prefix job names on the farm
        pipeline_name                 => 'run_topmed',

        seeding_module                => 'ReseqTrack::Hive::PipeSeed::BasePipeSeed',
        seeding_options               => {
            output_columns     => $self->o('sample_columns'),
            output_attributes  => $self->o('sample_attributes'),
            require_columns    => $self->o('require_sample_columns'),
            exclude_columns    => $self->o('exclude_sample_columns'),
            require_attributes => $self->o('require_sample_attributes'),
            exclude_attributes => $self->o('exclude_sample_attributes'),
        },
        require_run_attributes        => {},
        require_experiment_attributes => {},
        require_study_attributes      => {},
        require_sample_attributes     => {},
        exclude_run_attributes        => {},
        exclude_experiment_attributes => {},
        exclude_study_attributes      => {},
        exclude_sample_attributes     => {},
        require_experiment_columns    => { instrument_platform => [ 'ILLUMINA' ], },
        require_run_columns           => { status => [ 'public', 'private' ], },
        require_study_columns         => {},
        require_sample_columns        => {},
        exclude_sample_columns        => {},

        singularity_cache             => $self->o('ENV', 'SINGULARITY_CACHEDIR'),
        singularity_image             => 'broadinstitute/gtex_rnaseq:V8',

        # runnable-specific parameters' defaults:

        'lsf_queue'                   => 'production-rh74',

        'RGSM'                        => '#sample_source_id#',
        'RGPU'                        => '#run_source_id#',

        'sample_attributes'           => [],
        'sample_columns'              => [ 'sample_id', 'sample_source_id', 'sample_alias' ],
        'run_attributes'              => [],
        'run_columns'                 => [ 'run_source_id', 'center_name', 'run_alias' ],
        'study_attributes'            => [],
        'study_columns'               => [ 'study_source_id' ],
        'experiment_attributes'       => [],
        'experiment_columns'          => [ 'instrument_platform', 'paired_nominal_length' ],


        final_output_layout           => '#sample_source_id#/alignment',
        name_file_module              => 'ReseqTrack::Hive::NameFile::BaseNameFile',
        name_file_method              => 'basic',
        name_file_params              => {
            new_dir       => '#final_output_dir#/#final_output_layout#',
            new_basename  => '#sample_source_id#.bwa',
            add_datestamp => 1,
            suffix        => '.bam',
        },

    };
}

sub pipeline_create_commands {
    my ($self) = @_;

    return [
        @{$self->SUPER::pipeline_create_commands},
    ];
}

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        %{$self->SUPER::pipeline_wide_parameters},

        dir_label_params => [ 'study_source_id', 'sample_source_id', 'run_source_id' ],
    };
}


sub resource_classes {
    my ($self) = @_;
    return {
        %{$self->SUPER::resource_classes},
        '200Mb'  => { 'LSF' => '-C0 -M200 -q ' . $self->o('lsf_queue') . ' -R"select[mem>200] rusage[mem=200]"' },
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


sub pipeline_analyses {
    my ($self) = @_;

    my @analyses;

    push(@analyses, {
        -logic_name  => 'get_seeds',
        -module      => 'ReseqTrack::Hive::Process::SeedFactory',
        -meadow_type => 'LOCAL',
        -parameters  => {
            seeding_module  => $self->o('seeding_module'),
            seeding_options => $self->o('seeding_options'),
        },
        -flow_into   => {
            2 => [ 'libraries_factory' ],
        },
    });

    push(@analyses, {
        -logic_name  => 'libraries_factory',
        -module      => 'ReseqTrack::Hive::Process::RunMetaInfoFactory',
        -meadow_type => 'LOCAL',
        -parameters  => {
            factory_type                  => 'library',
            require_experiment_columns    => $self->o('require_experiment_columns'),
            require_study_columns         => $self->o('require_study_columns'),
            require_experiment_attributes => $self->o('require_experiment_attributes'),
            require_study_attributes      => $self->o('require_study_attributes'),
            exclude_experiment_attributes => $self->o('exclude_experiment_attributes'),
            exclude_study_attributes      => $self->o('exclude_study_attributes'),
        },
        -flow_into   => {
            2 => [ 'runs_factory' ],
        },
    });

    push(@analyses, {
        -logic_name  => 'runs_factory',
        -module      => 'ReseqTrack::Hive::Process::RunMetaInfoFactory',
        -meadow_type => 'LOCAL',
        -parameters  => {
            factory_type                  => 'run',
            require_experiment_columns    => $self->o('require_experiment_columns'),
            require_study_columns         => $self->o('require_study_columns'),
            require_run_columns           => $self->o('require_run_columns'),
            require_experiment_attributes => $self->o('require_experiment_attributes'),
            require_study_attributes      => $self->o('require_study_attributes'),
            require_run_attributes        => $self->o('require_run_attributes'),
            exclude_experiment_attributes => $self->o('exclude_experiment_attributes'),
            exclude_study_attributes      => $self->o('exclude_study_attributes'),
            exclude_run_attributes        => $self->o('exclude_run_attributes'),

            output_run_columns            => $self->o('run_columns'),
            output_study_columns          => $self->o('study_columns'),
            output_experiment_columns     => $self->o('experiment_columns'),
            output_run_attributes         => $self->o('run_attributes'),
            output_study_attributes       => $self->o('study_attributes'),
            output_experiment_attributes  => $self->o('experiment_attributes'),
        },
        -flow_into   => {
            2 => [ 'find_source_fastqs' ],
        },
    });

    push(@analyses, {
        -logic_name  => 'find_source_fastqs',
        -module      => 'ReseqTrack::Hive::Process::ImportCollection',
        -meadow_type => 'LOCAL',
        -parameters  => {
            collection_type    => $self->o('type_fastq'),
            collection_name    => '#run_source_id#',
            output_param       => 'fastq',
            reseqtrack_options => {
                flows_do_count_param => 'fastq',
                flows_do_count       => { 1 => '1+', },
            }
        },
        -flow_into   => {
            1 => [ 'star_align' ],
        },
    });

    push(@analyses, {
        -logic_name        => 'star_align',
        -language          => 'python3',
        -module            => 'PyHive.TOPMed.RunStar',
        -parameters        => {
            singularity_cache => $self->o('singularity_cache'),
            singularity_image => $self->o('singularity_image'),
            singularity_exe   => $self->o('singularity_exe'),
            star_index        => $self->o('star_index'),
            #output_directory  => $self->o('output_dir'),
            prefix            => 'star',
            num_threads       => 4
        },
        -rc_name           => '200Mb',
        -analysis_capacity => 4,
        -hive_capacity     => 200,
        #-flow_into         => { }
    });

    return \@analyses;
}

1;
