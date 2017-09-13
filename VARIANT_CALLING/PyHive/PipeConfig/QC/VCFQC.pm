package PyHive::PipeConfig::QC::VCFQC;

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
	'r_folder' => '/nfs/software/ensembl/RHEL7/linuxbrew/bin/', # for running PyHive.VcfQC.PlotVariantDensity
        'r_scripts' => '/nfs/gns/homes/ernesto/lib/igsr_analysis_master/igsr_analysis/VARIANT_CALLING/VCF/scripts/R/', # for running PyHive.VcfQC.PlotVariantDensity
	'genome_file' => undef, # for running PyHive.VcfQC.PlotVariantDensity
	'windows' => undef, # for running PyHive.VcfQC.PlotVariantDensity
	'length' => undef, # for running PyHive.VcfQC.PlotVariantDensity
	'genome_canonical' => '/nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/VARCALL_ALLGENOME_13022017/VCF_QC/DEVEL_PIPELINE/SMALL_INPUT/WGS_LC_bcftools/FILES/chr_canonicals.genome', # for running PyHive.VcfQC.PlotVariantDensity
	'genome_unlocalized_scaffolds' => '/nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/VARCALL_ALLGENOME_13022017/VCF_QC/DEVEL_PIPELINE/SMALL_INPUT/WGS_LC_bcftools/FILES/unlocalized_scaffolds.genome', # for running PyHive.VcfQC.PlotVariantDensity
	'genome_unplaced_scaffolds_part1' => '/nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/VARCALL_ALLGENOME_13022017/VCF_QC/DEVEL_PIPELINE/SMALL_INPUT/WGS_LC_bcftools/FILES/unplaced_scaffolds.part1.genome', # for running PyHive.VcfQC.PlotVariantDensity
	'genome_unplaced_scaffolds_part2' => '/nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/VARCALL_ALLGENOME_13022017/VCF_QC/DEVEL_PIPELINE/SMALL_INPUT/WGS_LC_bcftools/FILES/unplaced_scaffolds.part2.genome', # for running PyHive.VcfQC.PlotVariantDensity
	'genome_alt_scaffolds_part1' => '/nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/VARCALL_ALLGENOME_13022017/VCF_QC/DEVEL_PIPELINE/SMALL_INPUT/WGS_LC_bcftools/FILES/alt_scaffolds.part1.genome', # for running PyHive.VcfQC.PlotVariantDensity
	'genome_alt_scaffolds_part2' => '/nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/VARCALL_ALLGENOME_13022017/VCF_QC/DEVEL_PIPELINE/SMALL_INPUT/WGS_LC_bcftools/FILES/alt_scaffolds.part2.genome', # for running PyHive.VcfQC.PlotVariantDensity
	'genome_alt_scaffolds_part3' => '/nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/VARCALL_ALLGENOME_13022017/VCF_QC/DEVEL_PIPELINE/SMALL_INPUT/WGS_LC_bcftools/FILES/alt_scaffolds.part3.genome', # for running PyHive.VcfQC.PlotVariantDensity
	'chr_file' => undef, # for running PyHive.VcfQC.ChrosInVcf
	'store_attributes' => 'False',
	'bcftools_folder' => '/nfs/production/reseq-info/work/bin/bcftools-1.3/',
	'bedtools_folder' => '/nfs/gns/homes/ernesto/bin/bedtools-2.25.0/bin/',
	'region_bed' => undef, #for running PyHive.VcfQC.VariantsInRegions
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
                1 => {'chros_in_vcf' => {
		    'filepath'=> '#filepath#',
		    'layout_dict'=> '#layout_dict#'
		      }
		}
            },
        },

	{   -logic_name => 'chros_in_vcf',
            -module     => 'PyHive.VcfQC.ChrosInVcf',
            -language   => 'python3',
            -parameters => {
                'filepath'     => '#filepath#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
		'bcftools_folder' => $self->o('bcftools_folder'),
		'chr_file' => $self->o('chr_file'),
		'store_attributes' => $self->o('store_attributes'),
		'filter_str' => 'PASS,.'
	    },
	    -analysis_capacity => 5,
            -rc_name => '500Mb',
	    -flow_into => {
		1 => {'variants_inexcluded_regions' => {
		    'filepath'=> '#filepath#' ,
		    'layout_dict'=> '#layout_dict#'
		      }
		}
	    }
	},

	{   -logic_name => 'variants_inexcluded_regions',
            -module     => 'PyHive.VcfQC.VariantsInRegions',
            -language   => 'python3',
            -parameters => {
                'filepath'     => '#filepath#',
                'bedtools_folder' => $self->o('bedtools_folder'),
                'region' => $self->o('region_bed'),
                'work_dir' => $self->o('work_dir')
            },
	    -analysis_capacity => 1,
            -rc_name => '5Gb',
	    -flow_into => {
		1 => {'store_inexcluded_regions_counts' => {
		    'filepath'=> '#filepath#', 
		    'filename' => '#outfile#',
		    'layout_dict'=> '#layout_dict#'
		      }
		}
	    }
        },

	{   -logic_name => 'store_inexcluded_regions_counts',
            -module        => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'type' => 'VCF_INEXCLUDED_COUNTS',
                'final_dir' => $self->o('final_dir'),
                'newlayout' => $self->o('newlayout'),
                'add_date' => 'True',
                'extension' => 'counts.txt'
            },
	    -flow_into => {
		1 => {
		    'select_snps' => {'filepath' => '#filepath#'},
		    'select_indels' => {'filepath' => '#filepath#'}
		}
	    },
        },
	
	{   -logic_name => 'select_snps',
            -module        => 'PyHive.VcfFilter.SplitVariants',
            -language   => 'python3',
            -parameters    => {
                'filepath' => '#filepath#',
                'bcftools_folder' => $self->o('bcftools_folder'),
                'type' => 'snps',
		'work_dir' => $self->o('work_dir')
            },
	    -analysis_capacity => 20,
            -rc_name => '1Gb',
	    -flow_into => {
		1 => {
		    'plot_density_canonical_snps' => {'filepath' => '#out_vcf#'},
		    'plot_density_alt_scaffolds1_snps' => {'filepath' => '#out_vcf#'},
		    'plot_density_alt_scaffolds2_snps' => {'filepath' => '#out_vcf#'},
		    'plot_density_alt_scaffolds3_snps' => {'filepath' => '#out_vcf#'},
		    'plot_density_unloc_scaffolds_snps' => {'filepath' => '#out_vcf#'},
		    'plot_density_unplc_scaffolds1_snps' => {'filepath' => '#out_vcf#'},
		    'plot_density_unplc_scaffolds2_snps' => {'filepath' => '#out_vcf#'}
		},
	    },
        },

	{   -logic_name => 'plot_density_canonical_snps',
            -module     => 'PyHive.VcfQC.PlotVariantDensity',
            -language   => 'python3',
            -parameters => {
                'filepath'     => '#filepath#',
                'bcftools_folder' => $self->o('bcftools_folder'),
                'bedtools_folder' => $self->o('bedtools_folder'),
                'r_folder' => $self->o('r_folder'),
                'r_scripts' => $self->o('r_scripts'),
                'length' => 1000000,
		'height' => 1800,
		'chromName_cex' => 0.8,
                'genome' => $self->o('genome_canonical'),
                'work_dir' => $self->o('work_dir')
            },
	    -analysis_capacity => 5,
            -rc_name => '5Gb',
	    -flow_into => {
		1 => {'store_density_canonical_snps' => {'filename'=>'#density_plot_f#'}}
	    },
        },

	{   -logic_name => 'store_density_canonical_snps',
            -module        => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'type' => 'VCF_DENSITYPLOT',
                'final_dir' => $self->o('final_dir'),
                'newlayout' => $self->o('newlayout'),
                'add_date' => 'True',
                'extension' => 'snps.canonical.png'
            },
        },

	{   -logic_name => 'plot_density_alt_scaffolds1_snps',
            -module     => 'PyHive.VcfQC.PlotVariantDensity',
            -language   => 'python3',
            -parameters => {
                'filepath'     => '#filepath#',
                'bcftools_folder' => $self->o('bcftools_folder'),
                'bedtools_folder' => $self->o('bedtools_folder'),
                'r_folder' => $self->o('r_folder'),
                'r_scripts' => $self->o('r_scripts'),
                'length' => 10000,
		'height' => 11000,
                'chromName_cex' => 0.8,
                'genome' => $self->o('genome_alt_scaffolds_part1'),
                'work_dir' => $self->o('work_dir')
            },
	    -analysis_capacity => 5,
            -rc_name => '5Gb',
	    -flow_into => {
		1 => {'store_density_alt1_snps' => {'filename'=>'#density_plot_f#'}},
	    },
        },

	{   -logic_name => 'store_density_alt1_snps',
            -module        => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'type' => 'VCF_DENSITYPLOT',
                'final_dir' => $self->o('final_dir'),
                'newlayout' => $self->o('newlayout'),
                'add_date' => 'True',
                'extension' => 'snps.alt1.png'
            },
        },

	{   -logic_name => 'plot_density_alt_scaffolds2_snps',
            -module     => 'PyHive.VcfQC.PlotVariantDensity',
            -language   => 'python3',
            -parameters => {
                'filepath'     => '#filepath#',
                'bcftools_folder' => $self->o('bcftools_folder'),
                'bedtools_folder' => $self->o('bedtools_folder'),
                'r_folder' => $self->o('r_folder'),
                'r_scripts' => $self->o('r_scripts'),
		'length' => 10000,
                'height' => 11000,
                'chromName_cex' => 0.8,
                'genome' => $self->o('genome_alt_scaffolds_part2'),
                'work_dir' => $self->o('work_dir')
            },
	    -analysis_capacity => 5,
            -rc_name => '5Gb',
	    -flow_into => {
		1 => {'store_density_alt2_snps' => {'filename'=>'#density_plot_f#'}},
	    },
        },

	{   -logic_name => 'store_density_alt2_snps',
            -module        => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
                'username' => $self->o('username'),
                'hostname' => $self->o('hostname'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'type' => 'VCF_DENSITYPLOT',
                'final_dir' => $self->o('final_dir'),
                'newlayout' => $self->o('newlayout'),
                'add_date' => 'True',
                'extension' => 'snps.alt2.png'
            },
        },

	{   -logic_name => 'plot_density_alt_scaffolds3_snps',
            -module     => 'PyHive.VcfQC.PlotVariantDensity',
            -language   => 'python3',
            -parameters => {
                'filepath'     => '#filepath#',
                'bcftools_folder' => $self->o('bcftools_folder'),
                'bedtools_folder' => $self->o('bedtools_folder'),
                'r_folder' => $self->o('r_folder'),
                'r_scripts' => $self->o('r_scripts'),
		'length' => 10000,
                'height' => 11000,
                'chromName_cex' => 0.8,
                'genome' => $self->o('genome_alt_scaffolds_part3'),
                'work_dir' => $self->o('work_dir')
            },
            -analysis_capacity => 5,
            -rc_name => '5Gb',
	    -flow_into => {
		1 => {'store_density_alt3_snps' => {'filename'=>'#density_plot_f#'}},
	    },
        },

	{   -logic_name => 'store_density_alt3_snps',
            -module        => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
                'username' => $self->o('username'),
                'hostname' => $self->o('hostname'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'type' => 'VCF_DENSITYPLOT',
                'final_dir' => $self->o('final_dir'),
                'newlayout' => $self->o('newlayout'),
                'add_date' => 'True',
                'extension' => 'snps.alt3.png'
            },
        },

	{   -logic_name => 'plot_density_unloc_scaffolds_snps',
            -module     => 'PyHive.VcfQC.PlotVariantDensity',
            -language   => 'python3',
            -parameters => {
                'filepath'     => '#filepath#',
                'bcftools_folder' => $self->o('bcftools_folder'),
                'bedtools_folder' => $self->o('bedtools_folder'),
                'r_folder' => $self->o('r_folder'),
                'r_scripts' => $self->o('r_scripts'),
                'length' => 10000,
		'height' => 7000,
                'chromName_cex' => 0.8,
                'genome' => $self->o('genome_unlocalized_scaffolds'),
                'work_dir' => $self->o('work_dir')
            },
	    -analysis_capacity => 5,
            -rc_name => '5Gb',
	    -flow_into => {
		1 => {'store_density_unloc_scaffolds_snps' => {'filename'=>'#density_plot_f#'}},
	    },
        },
	
	{   -logic_name => 'store_density_unloc_scaffolds_snps',
            -module        => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'type' => 'VCF_DENSITYPLOT',
                'final_dir' => $self->o('final_dir'),
                'newlayout' => $self->o('newlayout'),
                'add_date' => 'True',
                'extension' => 'snps.unloc.png'
            },
            -analysis_capacity => 20,
            -rc_name => '500Mb'
        },

	{   -logic_name => 'plot_density_unplc_scaffolds1_snps',
            -module     => 'PyHive.VcfQC.PlotVariantDensity',
            -language   => 'python3',
            -parameters => {
                'filepath'     => '#filepath#',
                'bcftools_folder' => $self->o('bcftools_folder'),
                'bedtools_folder' => $self->o('bedtools_folder'),
                'r_folder' => $self->o('r_folder'),
                'r_scripts' => $self->o('r_scripts'),
		'length' => 10000,
                'height' => 7000,
                'chromName_cex' => 0.8,
                'genome' => $self->o('genome_unplaced_scaffolds_part1'),
                'work_dir' => $self->o('work_dir')
            },
	    -analysis_capacity => 5,
            -rc_name => '5Gb',
	    -flow_into => {
		1 => {'store_density_unplc_scaffolds1_snps' => {'filename'=>'#density_plot_f#'}},
	    },
        },

	{   -logic_name => 'store_density_unplc_scaffolds1_snps',
            -module        => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'type' => 'VCF_DENSITYPLOT',
                'final_dir' => $self->o('final_dir'),
                'newlayout' => $self->o('newlayout'),
                'add_date' => 'True',
                'extension' => 'snps.unplc1.png'
            },
            -analysis_capacity => 20,
            -rc_name => '500Mb'
        },

	{   -logic_name => 'plot_density_unplc_scaffolds2_snps',
            -module     => 'PyHive.VcfQC.PlotVariantDensity',
            -language   => 'python3',
            -parameters => {
                'filepath'     => '#filepath#',
                'bcftools_folder' => $self->o('bcftools_folder'),
                'bedtools_folder' => $self->o('bedtools_folder'),
                'r_folder' => $self->o('r_folder'),
                'r_scripts' => $self->o('r_scripts'),
		'length' => 10000,
                'height' => 7000,
                'chromName_cex' => 0.8,
                'genome' => $self->o('genome_unplaced_scaffolds_part2'),
                'work_dir' => $self->o('work_dir')
            },
            -analysis_capacity => 5,
            -rc_name => '5Gb',
	    -flow_into => {
		1 => {'store_density_unplc_scaffolds2_snps' => {'filename'=>'#density_plot_f#'}},
	    },
        },

	{   -logic_name => 'store_density_unplc_scaffolds2_snps',
            -module        => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'type' => 'VCF_DENSITYPLOT',
                'final_dir' => $self->o('final_dir'),
                'newlayout' => $self->o('newlayout'),
                'add_date' => 'True',
                'extension' => 'snps.unplc2.png'
            },
            -analysis_capacity => 20,
            -rc_name => '500Mb'
        },

        {   -logic_name => 'select_indels',
            -module        => 'PyHive.VcfFilter.SplitVariants',
            -language   => 'python3',
            -parameters    => {
                'filepath' => '#filepath#',
                'bcftools_folder' => $self->o('bcftools_folder'),
                'type' => 'indels',
		'work_dir' => $self->o('work_dir')
            },
            -analysis_capacity => 5,
            -rc_name => '1Gb',
	    -flow_into => {
		1 => {
		    'plot_density_canonical_indels' => {'filepath' => '#out_vcf#'},
		    'plot_density_alt_scaffolds1_indels' => {'filepath' => '#out_vcf#'},
		    'plot_density_alt_scaffolds2_indels' => {'filepath' => '#out_vcf#'},
		    'plot_density_alt_scaffolds3_indels' => {'filepath' => '#out_vcf#'},
		    'plot_density_unloc_scaffolds_indels' => {'filepath' => '#out_vcf#'},
		    'plot_density_unplc_scaffolds1_indels' => {'filepath' => '#out_vcf#'},
		    'plot_density_unplc_scaffolds2_indels' => {'filepath' => '#out_vcf#'}
		},
	    },
	},

	{   -logic_name => 'plot_density_canonical_indels',
            -module     => 'PyHive.VcfQC.PlotVariantDensity',
            -language   => 'python3',
            -parameters => {
                'filepath'     => '#filepath#',
                'bcftools_folder' => $self->o('bcftools_folder'),
                'bedtools_folder' => $self->o('bedtools_folder'),
                'r_folder' => $self->o('r_folder'),
                'r_scripts' => $self->o('r_scripts'),
                'length' => 1000000,
		'height' => 1800,
		'chromName_cex' => 0.8,
                'genome' => $self->o('genome_canonical'),
                'work_dir' => $self->o('work_dir')
            },
	    -analysis_capacity => 5,
            -rc_name => '5Gb',
	    -flow_into => {
		1 => {'store_density_canonical_indels' => {'filename'=>'#density_plot_f#'}}
	    }
        },

	{   -logic_name => 'store_density_canonical_indels',
            -module        => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'type' => 'VCF_DENSITYPLOT',
                'final_dir' => $self->o('final_dir'),
                'newlayout' => $self->o('newlayout'),
                'add_date' => 'True',
                'extension' => 'indels.canonical.png'
            }
        },

	{   -logic_name => 'plot_density_alt_scaffolds1_indels',
            -module     => 'PyHive.VcfQC.PlotVariantDensity',
            -language   => 'python3',
            -parameters => {
                'filepath'     => '#filepath#',
                'bcftools_folder' => $self->o('bcftools_folder'),
                'bedtools_folder' => $self->o('bedtools_folder'),
                'r_folder' => $self->o('r_folder'),
                'r_scripts' => $self->o('r_scripts'),
		'length' => 10000,
                'height' => 11000,
                'chromName_cex' => 0.8,
                'genome' => $self->o('genome_alt_scaffolds_part1'),
                'work_dir' => $self->o('work_dir')
            },
	    -analysis_capacity => 5,
            -rc_name => '5Gb',
	    -flow_into => {
		1 => {'store_density_alt1_indels' => {'filename'=>'#density_plot_f#'}},
	    }
        },

	{   -logic_name => 'store_density_alt1_indels',
            -module        => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'type' => 'VCF_DENSITYPLOT',
                'final_dir' => $self->o('final_dir'),
                'newlayout' => $self->o('newlayout'),
                'add_date' => 'True',
                'extension' => 'indels.alt1.png'
            },
        },

	{   -logic_name => 'plot_density_alt_scaffolds2_indels',
            -module     => 'PyHive.VcfQC.PlotVariantDensity',
            -language   => 'python3',
            -parameters => {
                'filepath'     => '#filepath#',
                'bcftools_folder' => $self->o('bcftools_folder'),
                'bedtools_folder' => $self->o('bedtools_folder'),
                'r_folder' => $self->o('r_folder'),
                'r_scripts' => $self->o('r_scripts'),
                'length' => 10000,
                'height' => 11000,
                'chromName_cex' => 0.8,
                'genome' => $self->o('genome_alt_scaffolds_part2'),
                'work_dir' => $self->o('work_dir')
            },
	    -analysis_capacity => 5,
            -rc_name => '5Gb',
	    -flow_into => {
		1 => {'store_density_alt2_indels' => {'filename'=>'#density_plot_f#'}},
	    },
        },

	{   -logic_name => 'store_density_alt2_indels',
            -module        => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
                'username' => $self->o('username'),
                'hostname' => $self->o('hostname'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'type' => 'VCF_DENSITYPLOT',
                'final_dir' => $self->o('final_dir'),
                'newlayout' => $self->o('newlayout'),
                'add_date' => 'True',
                'extension' => 'indels.alt2.png'
            },
        },

	{   -logic_name => 'plot_density_alt_scaffolds3_indels',
            -module     => 'PyHive.VcfQC.PlotVariantDensity',
            -language   => 'python3',
            -parameters => {
                'filepath'     => '#filepath#',
                'bcftools_folder' => $self->o('bcftools_folder'),
                'bedtools_folder' => $self->o('bedtools_folder'),
                'r_folder' => $self->o('r_folder'),
                'r_scripts' => $self->o('r_scripts'),
                'length' => 10000,
                'height' => 11000,
                'chromName_cex' => 0.8,
                'genome' => $self->o('genome_alt_scaffolds_part3'),
                'work_dir' => $self->o('work_dir')
            },
            -analysis_capacity => 5,
            -rc_name => '5Gb',
	    -flow_into => {
		1 => {'store_density_alt3_indels' => {'filename'=>'#density_plot_f#'}},
	    },
        },

	{   -logic_name => 'store_density_alt3_indels',
            -module        => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'type' => 'VCF_DENSITYPLOT',
                'final_dir' => $self->o('final_dir'),
                'newlayout' => $self->o('newlayout'),
                'add_date' => 'True',
                'extension' => 'indels.alt3.png'
            },
        },

	{   -logic_name => 'plot_density_unloc_scaffolds_indels',
            -module     => 'PyHive.VcfQC.PlotVariantDensity',
            -language   => 'python3',
            -parameters => {
                'filepath'     => '#filepath#',
                'bcftools_folder' => $self->o('bcftools_folder'),
                'bedtools_folder' => $self->o('bedtools_folder'),
                'r_folder' => $self->o('r_folder'),
                'r_scripts' => $self->o('r_scripts'),
                'length' => 10000,
                'height' => 7000,
                'chromName_cex' => 0.8,
                'genome' => $self->o('genome_unlocalized_scaffolds'),
                'work_dir' => $self->o('work_dir')
            },
	    -analysis_capacity => 5,
            -rc_name => '5Gb',
	    -flow_into => {
		1 => {'store_density_unloc_scaffolds_indels' => {'filename'=>'#density_plot_f#'}},
	    }
        },
	
	{   -logic_name => 'store_density_unloc_scaffolds_indels',
            -module        => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'type' => 'VCF_DENSITYPLOT',
                'final_dir' => $self->o('final_dir'),
                'newlayout' => $self->o('newlayout'),
                'add_date' => 'True',
                'extension' => 'indels.unloc.png'
            },
        },

	{   -logic_name => 'plot_density_unplc_scaffolds1_indels',
            -module     => 'PyHive.VcfQC.PlotVariantDensity',
            -language   => 'python3',
            -parameters => {
                'filepath'     => '#filepath#',
                'bcftools_folder' => $self->o('bcftools_folder'),
                'bedtools_folder' => $self->o('bedtools_folder'),
                'r_folder' => $self->o('r_folder'),
                'r_scripts' => $self->o('r_scripts'),
		'length' => 10000,
                'height' => 7000,
                'chromName_cex' => 0.8,
                'genome' => $self->o('genome_unplaced_scaffolds_part1'),
                'work_dir' => $self->o('work_dir')
            },
	    -analysis_capacity => 5,
            -rc_name => '5Gb',
	    -flow_into => {
		1 => {'store_density_unplc_scaffolds1_indels' => {'filename'=>'#density_plot_f#'}},
	    },
        },

	{   -logic_name => 'store_density_unplc_scaffolds1_indels',
            -module        => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'type' => 'VCF_DENSITYPLOT',
                'final_dir' => $self->o('final_dir'),
                'newlayout' => $self->o('newlayout'),
                'add_date' => 'True',
                'extension' => 'indels.unplc1.png'
            }
        },

	{   -logic_name => 'plot_density_unplc_scaffolds2_indels',
            -module     => 'PyHive.VcfQC.PlotVariantDensity',
            -language   => 'python3',
            -parameters => {
                'filepath'     => '#filepath#',
                'bcftools_folder' => $self->o('bcftools_folder'),
                'bedtools_folder' => $self->o('bedtools_folder'),
                'r_folder' => $self->o('r_folder'),
                'r_scripts' => $self->o('r_scripts'),
                'length' => 10000,
                'height' => 7000,
                'chromName_cex' => 0.8,
                'genome' => $self->o('genome_unplaced_scaffolds_part2'),
                'work_dir' => $self->o('work_dir')
            },
            -analysis_capacity => 5,
            -rc_name => '5Gb',
	    -flow_into => {
		1 => {'store_density_unplc_scaffolds2_indels' => {'filename'=>'#density_plot_f#'}},
	     },
        },

	{   -logic_name => 'store_density_unplc_scaffolds2_indels',
            -module        => 'PyHive.File.StoreFile',
            -language   => 'python3',
            -parameters    => {
                'filename' => '#filename#',
                'hostname' => $self->o('hostname'),
                'username' => $self->o('username'),
                'port' => $self->o('port'),
                'db' => $self->o('db'),
                'pwd' => $self->o('pwd'),
                'type' => 'VCF_DENSITYPLOT',
                'final_dir' => $self->o('final_dir'),
                'newlayout' => $self->o('newlayout'),
                'add_date' => 'True',
                'extension' => 'indels.unplc2.png'
            }
        }
	];
}

1;

