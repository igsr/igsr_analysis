=head1 NAME

 PyHive::PipeConfig::FILTER::VCFilterGATK

=head1 SYNOPSIS

This pipeline is used to filter a VCF generated using the Genome Analysis ToolKit (GATK).

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
-faix: This file is in the .faix format and controls the chromosomes in the initial VCF that will be analyzed
-bcftools_stats_region: This option will define the regions used by BCFTools stats: i.e. chr20
-intervals_f: This points to the Picard-style interval_list file used by Picard CollectVariantCallingMetrics to set 
 the chromosome used to calculate the stats
-truth_vcf: This points to the truth VCF used by Picard CollectVariantCallingMetrics in order to make calculations

GATK VQSR options:
-tranches: String representing a list used by GATK VQSR VariantRecalibrator in order to set the levels of truth 
 sensitivity at which to slice the data
-resources_snps: Path to json file with the location and prior probs of the sites used by GATK VariantRecalibrator
 to build the model for SNPs
-resources_indels: Path to json file with the location and prior probs of the sites used by GATK VariantRecalibrator
 to build the model for INDELs
-indels_annotations: List with annotations used by GATK VariantRecalibrators to build the model

=cut

package PyHive::PipeConfig::FILTER::VCFilterGATK;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');  # All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;


sub default_options {
    my ($self) = @_;
    return {
        %{ $self->SUPER::default_options() },               # inherit other stuff from the base class

        'pipeline_name' => 'run_filtering_on_gatk',                   # name used by the beekeeper to prefix job names on the farm

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
	'bcftools_folder' => '~/bin/bcftools-1.6/',
	'bgzip_folder' => '/nfs/production/reseq-info/work/ernesto/bin/anaconda3/bin/',
	'tabix_folder' => '/nfs/production/reseq-info/work/ernesto/bin/anaconda3/bin/',
	'exclude_bed' => '/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/exclude_nonvalid.bed',
	'faix' => undef, # faix file with chros that will be analyzed by the pipeline
	'bcftools_stats_region' => undef, # Define what chro will be analyzed by bcftools stats
	'picard_folder' => '/homes/ernesto/bin', # CollectVariantCallingMetrics
	'truth_vcf' => '/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/GATK_BUNDLE/dbsnp_146.hg38.vcf.gz', # CollectVariantCallingMetrics
	'intervals_f' => undef, # CollectVariantCallingMetrics
	'caller' => 'UG', # VariantRecalibrator
	'gatk_folder' => '/homes/ernesto/bin/GATK/', # VariantRecalibrator
	'java_tmpdir' => '/gpfs/nobackup/resequencing_informatics/ernesto/tmp', # necessary for GATK ApplyRecalibration not to crash 
	'tranches' => '[100.0,99.9,99.0,98.0,97.0,96.0,95.0,92.0,90.0,85.0,80.0,75.0,70.0,65.0,60.0,55.0,50.0]',
	'reference' => '/nfs/production/reseq-info/work/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa', # VariantRecalibrator
	'resources_snps' => '/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/GATK_BUNDLE/resources_snps.json', # VariantRecalibrator
	'resources_indels' => '/nfs/production/reseq-info/work/ernesto/isgr/SUPPORTING/REFERENCE/GATK_BUNDLE/resources_indels.json', # VariantRecalibrator
	'indels_annotations' => ['QD','DP','FS','SOR','ReadPosRankSum','MQRankSum','InbreedingCoeff'], #annotations for recalibrating indels
	'vt_folder' => '/homes/ernesto/bin/vt/',
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
	'50Gb' => { 'LSF' => '-C0 -M50000 -q '.$self->o('lsf_queue').' -R"select[mem>50000] rusage[mem=50000]"' },
	'10cpus' => { 'LSF' => '-n 10 -C0 -M1024 -q '.$self->o('lsf_queue').' -R"select[mem>1024] rusage[mem=1024]"' },
	'20cpus' => { 'LSF' => '-n 20 -q '.$self->o('lsf_queue') },
        '30cpus' => { 'LSF' => '-n 30 -q '.$self->o('lsf_queue') }
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
		1 => {'store_stats1' => {
		    'filename' => '#stats_file#',
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
               'work_dir' => $self->o('work_dir'),
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
		1 => {'split_chr' => {
		    'filepath' => '#initial_file#'
		      }
		}
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
                1 => {'run_exclude_regions' => {
		    'filepath' => '#vcf_file#',
		    'ix' => '#ix#'
		      }
		},
            },
            -analysis_capacity => 400,
            -rc_name => '2Gb'
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
                1 => [ '?accu_name=allchunks_files&accu_address=[]&accu_input_variable=subset_file', '?accu_name=allixs&accu_address=[]&accu_input_variable=ix']
            },
            -analysis_capacity => 400,
            -rc_name => '2Gb'
        },

	{   -logic_name => 'merge_vcf',
            -module        => 'PyHive.Vcf.VcfConcat',
            -language   => 'python3',
            -parameters    => {
                'outprefix' => $self->o('work_dir')."/#initial_filename#.merged.vcf.gz",
                'bcftools_folder' => $self->o('bcftools_folder'),
		'log_dir' => $self->o('log_dir'),
                'verbose' => 'True',
                'work_dir' => $self->o('work_dir')
            },

	    -flow_into => {
		1 => {'index_vcf1' => {'filepath' => '#merged_file#'}}
	    },
            -analysis_capacity => 1,
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
                1 => {'run_variantrecalibrator_snps' => {
		        'filepath' => '#filepath#'
		      }
		}
            },
            -analysis_capacity => 1,
            -rc_name => '500Mb'
        },

	{   -logic_name => 'run_variantrecalibrator_snps',
            -module        => 'PyHive.VcfFilter.VariantRecalibrator',
            -language   => 'python3',
            -parameters    => {
                'filepath' => '#filepath#',
                'work_dir' => $self->o('work_dir'),
                'caller' => $self->o('caller'),
		'log_file' => $self->o('log_dir')."/gatk_variantrecalibratior",
                'gatk_folder' => $self->o('gatk_folder'),
                'reference' => $self->o('reference'),
		'resources' => $self->o('resources_snps'),
		'tranches' => $self->o('tranches'),
#		'intervals' => 'chr20:10000000-20000000',
		'mode' => 'SNP'
            },
	    -flow_into => {
		1 => {
		    'run_applyrecalibration_snps' => {
			'filepath' => '#filepath#',
			'recal_file' => '#recal_f#',
			'tranches_file' => '#tranches_f#'
		    }},
	    },
            -analysis_capacity => 20,
            -rc_name => '12Gb',
        },

	{   -logic_name => 'run_applyrecalibration_snps',
            -module        => 'PyHive.VcfFilter.ApplyRecalibration',
            -language   => 'python3',
            -parameters    => {
                'filepath' => '#filepath#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'work_dir' => $self->o('work_dir'),
                'caller' => $self->o('caller'),
                'gatk_folder' => $self->o('gatk_folder'),
		'bgzip_folder' => $self->o('bgzip_folder'),
		'tabix_folder' => $self->o('tabix_folder'),
		'tmp_dir' => $self->o('java_tmpdir'),
                'reference' => $self->o('reference'),
                'recal_file' => '#recal_f#',
                'tranches_file' => '#tranches_f#',
		'ts_filter_level' => 99.9,
                'mode' => 'SNP'
            },
	    -flow_into => {
		1 => {
		    'run_variantrecalibrator_indels' => {
			'filepath' => '#vcf_filt#',
			'recal_file' => '#recal_f#',
                        'tranches_file' => '#tranches_f#'
		    }},
	    },
            -analysis_capacity => 20,
            -rc_name => '12Gb',
        },

	{   -logic_name => 'run_variantrecalibrator_indels',
            -module        => 'PyHive.VcfFilter.VariantRecalibrator',
            -language   => 'python3',
            -parameters    => {
                'filepath' => '#filepath#',
                'work_dir' => $self->o('work_dir'),
                'caller' => $self->o('caller'),
                'gatk_folder' => $self->o('gatk_folder'),
                'reference' => $self->o('reference'),
                'resources' => $self->o('resources_indels'),
		'tranches' => $self->o('tranches'),
		'log_file' => $self->o('log_dir')."/gatk_variantrecalibratior",
#		'intervals' => 'chr20:10000000-20000000',
                'mode' => 'INDEL',
                'annotations' => $self->o('indels_annotations'),
                'max_gaussians' => 4,
	    },
	    -flow_into => {
		1 => {
		    'run_applyrecalibration_indels' => {
			'filepath' => '#filepath#',
			'recal_file' => '#recal_f#',
			'tranches_file' => '#tranches_f#'
		    }},
	    },
            -analysis_capacity => 20,
            -rc_name => '12Gb',
        },

	{   -logic_name => 'run_applyrecalibration_indels',
            -module        => 'PyHive.VcfFilter.ApplyRecalibration',
            -language   => 'python3',
            -parameters    => {
                'filepath' => '#filepath#',
                'work_dir' => $self->o('work_dir'),
                'caller' => $self->o('caller'),
                'gatk_folder' => $self->o('gatk_folder'),
		'bgzip_folder' => $self->o('bgzip_folder'),
                'tabix_folder' => $self->o('tabix_folder'),
                'reference' => $self->o('reference'),
		'tmp_dir' => $self->o('java_tmpdir'),
                'ts_filter_level' => 99.5,
                'recal_file' => '#recal_f#',
                'tranches_file' => '#tranches_f#',
                'mode' => 'INDEL'
            },
	    -flow_into => {
		1 => {
		    'run_bcftools_stats2' => {
                        'filepath' => '#vcf_filt#',
			 'vcf_ix' => '#vcf_filt_ix#'
		    }},
	    },
            -analysis_capacity => 20,
            -rc_name => '12Gb',
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
		'region' => $self->o('bcftools_stats_region'),
		'verbose' => 'True',
		'filter_str' => 'PASS,.',
                'work_dir' => $self->o('work_dir'),
                'store_attributes' => $self->o('store_attributes'),
                'bcftools_folder' => $self->o('bcftools_folder')
            },
	    -flow_into => {
		1 => {'store_stats2' => {
		    'filename' => '#stats_file#',
		    'filt_file' => '#filepath#',
		    'vcf_ix' => '#vcf_ix#'
		      }
		}
	    },
            -analysis_capacity => 20,
            -rc_name => '500Mb',
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
		    'filepath' => '#filt_file#',
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
		1 => {'store_cvcmetrics_detailmetrics2' => {
		    'filename' => '#detail_metrics#', 
		    'summary_metrics' => '#summary_metrics#', 
		    'filt_file' => '#filepath#',
		    'vcf_ix' => '#vcf_ix#'
		      }
		}
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
                       'filt_file' => '#filt_file#',
		       'vcf_ix' => '#vcf_ix#'
		    }},
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
		1 => {'store_filt_file' => {
		    'filename' => '#filt_file#',
		    'vcf_ix' => '#vcf_ix#'
		      }
		}
	    },
            -analysis_capacity => 20,
            -rc_name => '500Mb',
        },

	{   -logic_name => 'store_filt_file',
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
                    'store_filt_file_ix' => {
                        'filename' => '#vcf_ix#',
                    }}
            },
            -analysis_capacity => 20,
            -rc_name => '500Mb'
        },

	{   -logic_name => 'store_filt_file_ix',
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
