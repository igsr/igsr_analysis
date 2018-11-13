=head1 NAME

 PyHive::PipeConfig::FILTER::VCFilterFreeBayes

=head1 SYNOPSIS

This pipeline is used to filter a VCF generated using Freebayes caller.

Below in the 'default_options' function we can see the options that will control the behaviour
of the pipeline. The options that do not have a default value must be set when initializing the
pipeline using 'init_pipeline.pl'. Here are explanations for some of these options, modify them according
to your needs:

-hostname, username, port, db, pwd control the connection details for the ReseqTrack database
-work_dir: folder that will be used to put the intermediate files
-final_dir: folder that will be used to put the final pipeline files

-store_attributes: Possible values are 'True'/'False'. If 'True' then the pipeline will store the stats
 on the VCF calculated by some of the programs used by this pipeline (i.e. Picard CollectVariantCallingMetrics,
 BCFTools stats) in the 'Atttribute' table of the ReseqTrack database
-filelayout: this is used by the pipeline in order to know how to construct the final filename. It will label
 each bit in the initial filename with a certain name that will be used by the 'newlayout' parameter in order
 to construct the final filename. For example: 'dataset,caller,date,extension,compression'
 for a file named 'lc_ex_bam.freebayes.20170911.chr20_10e6_11e6.vcf.gz'
-newlayout: if newlayout is 'dataset,caller,extension,compression' and given the value used for 'filelayout'
 then the final filename will be 'lc_ex_bam.freebayes.vcf.gz'

-exclude_bed: Path to .BED file containing regions that will be filtered from the final VCF file (i.e. centromeres,
 chrX and chrY non PAR regions, etc.)
-faix: This file it is in the .faix format and controls the chromosomes in the in initial VCF that will be analyzed
-bcftools_stats_region: This option will define the regions used by BCFTools stats: i.e. chr20
-intervals_f: This points to the Picard-style interval_list file used by Picard CollectVariantCallingMetrics to set
 the chromosome used to calculate the stats
-truth_vcf: This points to the truth VCF used by Picard CollectVariantCallingMetrics in order to make calculations

Options used by 'bcftools filter -e' to do the filtering:
-f_expression: Filter expression that will be used with the bcftools filter -e parameter for filtering. Example:
 'QUAL<1'
-f_name: Used to specify the filter label used for SNPs and INDELs in the FILTER column of the VCF 

=cut

package PyHive::PipeConfig::FILTER::VCFilterFreeBayes;

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
	'store_attributes' => 'False',
	'newlayout' => undef, # new layout used for the generated files. i.e. [ 'dataset','caller']
        'filelayout' => undef, # layout of file that is analyzed by the pipeline. i.e. [ 'dataset','caller','date','extension','compression']
	'bcftools_folder' => '/nfs/production/reseq-info/work/bin/bcftools-1.3/',
	'exclude_bed' => '/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/exclude_nonvalid.bed',
        'bgzip_folder' => '/nfs/software/ensembl/RHEL7/linuxbrew/bin/',
	'tabix_folder' => '/nfs/software/ensembl/RHEL7/linuxbrew/bin/',
	'f_expression' => undef,
	'f_name' => undef,
	'bcftools_stats_region' => undef, # Define what chro will be analyzed by bcftools stats
	'picard_folder' => '/homes/ernesto/bin', # CollectVariantCallingMetrics
	'truth_vcf' => '/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/GATK_BUNDLE/dbsnp_146.hg38.vcf.gz', # CollectVariantCallingMetrics
	'intervals_f' => undef, # Define what chr to analyze by CollectVariantCallingMetrics
	'faix' => undef, # this controls what chros will be analyzed
	'vt_folder' => '/homes/ernesto/bin/vt/',
	'gatk_folder' => '/homes/ernesto/bin/GATK/',
	'reference' => '/nfs/production/reseq-info/work/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa',
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
	'15Gb' => { 'LSF' => '-C0 -M15360 -q '.$self->o('lsf_queue').' -R"select[mem>15360] rusage[mem=15360]"' },
	'20Gb' => { 'LSF' => '-C0 -M20000 -q '.$self->o('lsf_queue').' -R"select[mem>20000] rusage[mem=20000]"' },
	'10cpus' => { 'LSF' => '-n 10 -C0 -M1024 -q '.$self->o('lsf_queue').' -R"select[mem>1024] rusage[mem=1024]"' },
	'20cpus' => { 'LSF' => '-n 20 -C0 -M1024 -q '.$self->o('lsf_queue').' -R"select[mem>1024] rusage[mem=1024]"' },
	'30cpus' => { 'LSF' => '-n 30 -C0 -M1024 -q '.$self->o('lsf_queue').' -R"select[mem>1024] rusage[mem=1024]"' }
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
                1 => {'run_bcftools_stats1' => {
		    'filepath'=> '#filepath#', 
		    'layout_dict'=> '#layout_dict#'
		      }
		}
            },
        },

	{   -logic_name => 'run_bcftools_stats1',
            -module        => 'PyHive.VcfQC.BcftoolsStats',
            -language   => 'python3',
            -parameters    => {
                'filepath' => '#filepath#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
		'region' => $self->o('bcftools_stats_region'),
		'verbose' => 'True',
		'work_dir' => $self->o('work_dir'),
                'store_attributes' => $self->o('store_attributes'),
                'bcftools_folder' => $self->o('bcftools_folder')
            },
            -analysis_capacity => 20,
            -rc_name => '500Mb',
	    -flow_into => {
		1 => {'store_stats1' => {'filename' => '#stats_file#',
					 'initial_file' => '#filepath#',
					 'layout_dict'=> '#layout_dict#'
		      }
		}
	    },
	},

        {   -logic_name => 'store_stats1',
            -module        => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'type' => 'UNFILT_VCF_STATS',
                'final_dir' => $self->o('final_dir'),
                'newlayout' => $self->o('newlayout'),
                'add_date' => 'True',
                'extension' => 'unfilt.stats'
            },
            -analysis_capacity => 20,
            -rc_name => '500Mb',
            -flow_into => {
		1 => {'run_cvcmetrics1' => {'filepath' => '#initial_file#'}} 
	    },
        },

	{   -logic_name => 'run_cvcmetrics1',
            -module        => 'PyHive.VcfQC.CollectVariantCallingMetrics',
            -language   => 'python3',
            -parameters    => {
               'filepath' => '#filepath#',
               'hostname' => $self->o('hostname'),
               'username' => $self->o('username'),
               'port' => $self->o('port'),
               'db' => $self->o('db'),
               'pwd' => $self->o('pwd'),
               'store_attributes' => $self->o('store_attributes'),
               'picard_folder' => $self->o('picard_folder'),
               'truth_vcf' => $self->o('truth_vcf'),
               'intervals' => $self->o('intervals_f')
            },
	    -flow_into => {
               1 => {'store_cvcmetrics_detailmetrics1' => {
		   'filename' => '#detail_metrics#', 
		   'summary_metrics' => '#summary_metrics#', 
		   'initial_file' => '#filepath#' 
		     }
	       }
            },
            -analysis_capacity => 20,
            -rc_name => '12Gb',
        },

	{   -logic_name => 'store_cvcmetrics_detailmetrics1',
            -module        => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'type' => 'UNFILT_VCF_DETAILMETRICS',
                'final_dir' => $self->o('final_dir'),
                'newlayout' => $self->o('newlayout'),
                'add_date' => 'True',
                'extension' => 'unfilt.variant_calling_detail_metrics'
            },
	    -flow_into => {
		1 => {
		    'store_cvcmetrics_summarymetrics1' => {
			'filename' => '#summary_metrics#',
			'initial_file' => '#initial_file#'
		    }},
	    },
	    -analysis_capacity => 20,
            -rc_name => '500Mb',
        },

	{   -logic_name => 'store_cvcmetrics_summarymetrics1',
            -module        => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'type' => 'UNFILT_VCF_SUMMARYMETRICS',
                'final_dir' => $self->o('final_dir'),
                'newlayout' => $self->o('newlayout'),
                'add_date' => 'True',
                'extension' => 'unfilt.variant_calling_summary_metrics'
            },
	    -flow_into => {
		1 => {'split_chr' => {'filepath' => '#initial_file#'}}
            },
            -analysis_capacity => 20,
            -rc_name => '500Mb',
        },

	{   -logic_name => 'split_chr',
            -module        => 'PyHive.Factories.SplitVCFintoChros',
            -language   => 'python3',
            -parameters    => {
                'filepath' => '#filepath#',
                'bcftools_folder' => $self->o('bcftools_folder'),
		'faix' => $self->o('faix'),
                'threads' => 10,
                'verbose' => 'True',
		'work_dir' => $self->o('work_dir')
            },
	    -flow_into => {
                '2->A' => { 'norm_vcf' => { 
		    'filepath' => '#chr#',
		    'ix' => '#ix#'
			    }
		},
                'A->1' => [ 'merge_vcf'],
		},
            -analysis_capacity => 1,
            -rc_name => '10cpus'
        },
	
	{   -logic_name => 'norm_vcf',
            -module        => 'PyHive.Vcf.NormVcf',
            -language   => 'python3',
            -parameters    => {
                'filepath' => '#filepath#',
		'vt_folder' => $self->o('vt_folder'),
		'bgzip_folder' => $self->o('bgzip_folder'),
		'compress' => 1,
		'reference' => $self->o('reference'),
            },
	    -flow_into => {
		1 => {'index_vcf1' => {
		    'filepath' => '#vcf_file#',
		    'ix' => '#ix#'
		      }
		},
  	    },
            -analysis_capacity => 400,
            -rc_name => '500Mb'
        },

	{   -logic_name => 'index_vcf1',
            -module        => 'PyHive.Vcf.VcfIxByTabix',
            -language   => 'python3',
            -parameters    => {
                'filepath' => '#filepath#',
                'tabix_folder' => $self->o('tabix_folder'),
            },
            -flow_into => {
                1 => {'vcf_allelicprimitives' => {
                    'filepath' => '#filepath#',
                    'ix' => '#ix#'
                      }
                }
            },
            -analysis_capacity => 400,
            -rc_name => '500Mb'
        },

	{   -logic_name => 'vcf_allelicprimitives',
            -module        => 'PyHive.Vcf.VariantsToAllelicPrimitives',
            -language   => 'python3',
            -parameters    => {
                'filepath' => '#filepath#',
                'gatk_folder' => $self->o('gatk_folder'),
                'bgzip_folder' => $self->o('bgzip_folder'),
		'reference' => $self->o('reference'),
		'compress' => 'True'
            },
	    -flow_into => {
		1 => {'index_vcf2' => {
                    'filepath' => '#vcf_file#',
                    'ix' => '#ix#'
		      }
		},
	    },
            -analysis_capacity => 400,
            -rc_name => '20Gb'
        },

	{   -logic_name => 'index_vcf2',
            -module        => 'PyHive.Vcf.VcfIxByTabix',
            -language   => 'python3',
            -parameters    => {
                'filepath' => '#filepath#',
                'tabix_folder' => $self->o('tabix_folder'),
            },
            -flow_into => {
		1 => {'run_exclude_regions' => {
                    'filepath' => '#filepath#',
                    'ix' => '#ix#'
                      }
                }
            },
            -analysis_capacity => 400,
            -rc_name => '500Mb'
        },


	{   -logic_name => 'run_exclude_regions',
            -module        => 'PyHive.VcfFilter.SubsetVcfWithBed',
            -language   => 'python3',
            -parameters    => {
                'filepath' => '#filepath#',
                'bcftools_folder' => $self->o('bcftools_folder'),
                'threads' => 1,
		'verbose' => 'True',
                'action' => 'exclude',
                'bed' => $self->o('exclude_bed'),
            },
	    -flow_into => {
		1  =>  {'bcftools_filter' => {
		    'filepath' => '#subset_file#'
			}
		},
	    },
            -analysis_capacity => 400,
            -rc_name => '2Gb'
        },


        {   -logic_name => 'bcftools_filter',
            -module        => 'PyHive.VcfFilter.BcftoolsFilter',
            -language   => 'python3',
            -parameters    => {
                'filepath' => '#filepath#',
                'bcftools_folder' => $self->o('bcftools_folder'),
                'filter_expression' => $self->o('f_expression'),
                'filter_name' => $self->o('f_name')
            },
            -analysis_capacity => 200,
            -rc_name => '500Mb',
	    -flow_into => {
		1 => {'index_vcf3' => {
		    'filepath' => '#out_vcf#',
		    'ix' => '#ix#'
		      }
		},
	    },
        },

        {   -logic_name => 'index_vcf3',
            -module        => 'PyHive.Vcf.VcfIxByTabix',
            -language   => 'python3',
            -parameters    => {
                'filepath' => '#filepath#',
                'tabix_folder' => $self->o('tabix_folder'),
            },
            -flow_into => {
		1 => [ '?accu_name=allchunks_files&accu_address=[]&accu_input_variable=filepath', '?accu_name=allixs&accu_address=[]&accu_input_variable=ix']
            },
            -analysis_capacity => 200,
            -rc_name => '500Mb'
        },

	{   -logic_name => 'merge_vcf',
            -module        => 'PyHive.Vcf.VcfConcat',
            -language   => 'python3',
            -parameters    => {
		'outprefix' => $self->o('work_dir')."/#initial_filename#.merged.vcf.gz",
                'bcftools_folder' => $self->o('bcftools_folder'),
		'log_dir' => $self->o('log_dir'),
		'work_dir'=> $self->o('work_dir'),
		'verbose' => "True"
            },
	    -flow_into => {
                1 => {'index_vcf4' => {'filepath' => '#merged_file#'}}
            },
            -analysis_capacity => 1,
            -rc_name => '500Mb'
        },

	{   -logic_name => 'index_vcf4',
            -module        => 'PyHive.Vcf.VcfIxByTabix',
            -language   => 'python3',
            -parameters    => {
                'filepath' => '#filepath#',
                'tabix_folder' => $self->o('tabix_folder'),
            },
            -flow_into => {
                1 => {'run_bcftools_stats2' => {
		    'filepath' => '#filepath#',
		    'vcf_ix' => '#vcf_ix#'
		      }
		}
            },
            -analysis_capacity => 1,
            -rc_name => '500Mb'
        },

	{   -logic_name => 'run_bcftools_stats2',
            -module        => 'PyHive.VcfQC.BcftoolsStats',
            -language   => 'python3',
            -parameters    => {
                'filepath' => '#filepath#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
		'verbose' => 'True',
		'region' => $self->o('bcftools_stats_region'),
		'filter_str' => 'PASS,.',
                'work_dir' => $self->o('work_dir'),
                'store_attributes' => $self->o('store_attributes'),
                'bcftools_folder' => $self->o('bcftools_folder')
            },
            -analysis_capacity => 20,
            -rc_name => '500Mb',
	    -flow_into => {
		1 => {'store_stats2' => {
		    'filename' => '#stats_file#',
		    'merged_file' => '#filepath#',
		    'vcf_ix' => '#vcf_ix#'
		      }
		}
	    },
        },

	{   -logic_name => 'store_stats2',
            -module        => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'type' => 'FILT_VCF_STATS',
                'final_dir' => $self->o('final_dir'),
                'newlayout' => $self->o('newlayout'),
                'add_date' => 'True',
                'extension' => 'filt.stats'
            },
	    -flow_into => {
		1 => {'run_cvcmetrics2' => {
		    'merged_file' => '#merged_file#',
		    'vcf_ix' => '#vcf_ix#'
		      }
		}
	    },
            -analysis_capacity => 20,
            -rc_name => '500Mb'
        },

	{   -logic_name => 'run_cvcmetrics2',
            -module        => 'PyHive.VcfQC.CollectVariantCallingMetrics',
            -language   => 'python3',
            -parameters    => {
               'filepath' => '#filepath#',
               'hostname' => $self->o('hostname'),
               'username' => $self->o('username'),
               'port' => $self->o('port'),
               'db' => $self->o('db'),
               'pwd' => $self->o('pwd'),
               'store_attributes' => $self->o('store_attributes'),
               'work_dir' => $self->o('work_dir'),
               'picard_folder' => $self->o('picard_folder'),
               'truth_vcf' => $self->o('truth_vcf'),
               'intervals' => $self->o('intervals_f')
            },
	    -flow_into => {
	       1 => {'store_cvcmetrics_detailmetrics2' => {'filename' => '#detail_metrics#',
                                                           'summary_metrics' => '#summary_metrics#',
                                                           'merged_file' => '#filepath#',
                                                           'vcf_ix' => '#vcf_ix#'
		     }},
	    },
            -analysis_capacity => 20,
            -rc_name => '12Gb',
        },

	{   -logic_name => 'store_cvcmetrics_detailmetrics2',
            -module        => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'type' => 'FILT_VCF_DETAILMETRICS',
                'final_dir' => $self->o('final_dir'),
                'newlayout' => $self->o('newlayout'),
                'add_date' => 'True',
                'extension' => 'filt.variant_calling_detail_metrics'
            },
	    -flow_into => {
		1 => {
		    'store_cvcmetrics_summarymetrics2' => {
                        'filename' => '#summary_metrics#',
                        'merged_file' => '#merged_file#',
			'vcf_ix' => '#vcf_ix#'
		    }
		},
	    },
            -analysis_capacity => 20,
            -rc_name => '500Mb',
        },

        {   -logic_name => 'store_cvcmetrics_summarymetrics2',
            -module        => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'type' => 'FILT_VCF_SUMMARYMETRICS',
                'final_dir' => $self->o('final_dir'),
                'newlayout' => $self->o('newlayout'),
                'add_date' => 'True',
                'extension' => 'filt.variant_calling_summary_metrics'
            },
	    -flow_into => {
		1  =>  {
		    'store_filtered_file' => {
			'filename' => '#merged_file#',
			'vcf_ix' => '#vcf_ix#'
		    }
		},
	     },
            -analysis_capacity => 1,
            -rc_name => '500Mb'
        },

	
	{   -logic_name => 'store_filtered_file',
            -module        => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'type' => 'FILT_VCF',
                'final_dir' => $self->o('final_dir'),
                'newlayout' => $self->o('newlayout'),
                'add_date' => 'True',
                'extension' => 'filt.vcf.gz',

            },
            -flow_into => {
                1 => {
		    'store_filtered_file_ix' => {
			'filename' => '#vcf_ix#',
		    }}
            },
            -analysis_capacity => 1,
            -rc_name => '500Mb'
        },

	{   -logic_name => 'store_filtered_file_ix',
            -module        => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'type' => 'FILT_VCF_IX',
                'final_dir' => $self->o('final_dir'),
                'newlayout' => $self->o('newlayout'),
                'add_date' => 'True',
                'extension' => 'filt.vcf.gz.tbi',

            },
            -analysis_capacity => 1,
            -rc_name => '500Mb'
	}
	];
}

1;

