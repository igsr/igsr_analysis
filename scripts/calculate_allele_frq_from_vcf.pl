#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use IPC::System::Simple qw(system);

my ($vcf,
	$sample_panel,
	$user_input_pop_string,
	$region,
	$out_dir,
  $user_out_file,
	$tabix,
	$no_tabix,
	$help,
  %sample_pop_hash,#keys are sample IDs and values the population for that sample
  %pops,#keys of this hash are the populations in the sample panel file
  @selected_pops,#processed list of selected pops
  @pop_cols,#record columns belonging to the selected populations
  $vcf_cmd,#command to retrieve VCF
);
	
&GetOptions( 
	'vcf=s'					=> \$vcf,
	'sample_panel=s'		=> \$sample_panel,
	'pop=s'					=> \$user_input_pop_string,
	'region=s'				=> \$region,
	'out_dir:s'				=> \$out_dir,
  'out_file:s'      => \$user_out_file,
	'tabix=s'				=> \$tabix,
	'no_tabix!'				=> \$no_tabix,
	'help!'					=> \$help,
);

if ($help) {
	 exec('perldoc', $0);
}

check_and_process_input();

#input and output files
open(VCF, $vcf_cmd." |") or die "Failed to open $vcf_cmd $!\n";
my $outfile;
if($user_out_file){
  $outfile = $user_out_file;
}
else{
  $outfile = "$out_dir/afc.".($region =~ s/:/./r).".proc$$.tsv";
}
open(OUT, ">", $outfile) or die "Cannot open $outfile.\n";

#print header
print OUT "#CHROM\tPOS\tID\tREF\tALT\t";
foreach my $sp (@selected_pops){
  print OUT $sp."_TOTAL_CNT\t".$sp."_ALT_CNT\t".$sp."_FRQ\t";
}
print OUT "ALL_TOTAL_CNT\tALL_ALT_CNT\tALL_FRQ\n";

#process VCF
my $line;
LINE: while(<VCF>){
  chomp;
  $line = $_;
  #ignore general header lines
  next LINE if($line =~ /^##/);
  #process column header line to note which samples are in which columns
  if($line =~ /^#[Cc]/){
    my @cols = split "\t", $line;
    #for each sample in header
    SAMPLE: for(my $i=9;$i<scalar(@cols);$i++){
      for(my $j=0;$j<scalar(@selected_pops);$j++){
        if($selected_pops[$j] eq $sample_pop_hash{$cols[$i]}){
          push @{$pop_cols[$j]}, $i;
          next SAMPLE;
        }
      }
    }
  }#process lines from VCF body
  else{
    my @body_fields = split "\t", $line;
    #print chr, pos, ID, ref and alt cols 
    for(my $i=0;$i<5;$i++){
      print OUT $body_fields[$i]."\t";
    }
    my $num_alleles = scalar(split ",", $body_fields[4]) + 1;
    #keep score for all pops combined
    my @all_allele_cnt = ();
    #start allele counts at 0
    for(my $a=0;$a<$num_alleles;$a++){
      $all_allele_cnt[$a] = 0;
    } 
    #process each population    
    for(my $i=0;$i<scalar(@selected_pops);$i++){
      my @pop_allele_cnt = ();
      #start pop specific allele counts at 0
      for(my $a=0;$a<$num_alleles;$a++){
        $pop_allele_cnt[$a] = 0;
      }
      #process the sample genotypes for this population
      foreach my $sample_index(@{$pop_cols[$i]}){
        my @alleles = split /[|\/]/, $body_fields[$sample_index];
        foreach my $allele (@alleles){
          $pop_allele_cnt[$allele]++;
        }
      }
      #handle output for the population, including adding to total counts
      my $pop_total = 0;
      for(my $a=0;$a<$num_alleles;$a++){
        $all_allele_cnt[$a] += $pop_allele_cnt[$a];
        $pop_total += $pop_allele_cnt[$a];
      }
      print OUT $pop_total."\t";
      print OUT join(',', @pop_allele_cnt[1..($num_alleles-1)])."\t";
      my @alt_freqs = ();
      foreach my $alt_cnt (@pop_allele_cnt[1..($num_alleles-1)]){
         push @alt_freqs, (sprintf "%.2f", $alt_cnt/$pop_total);
      }
      print OUT join(',', @alt_freqs)."\t";
    }
    #output for combined populations
    my $all_total = 0;
    for(my $a=0;$a<$num_alleles;$a++){
      $all_total += $all_allele_cnt[$a];
    }
    print OUT $all_total."\t";
    print OUT join(',', @all_allele_cnt[1..($num_alleles-1)])."\t";
    my @all_alt_freqs = ();
    foreach my $all_alt_cnt (@all_allele_cnt[1..($num_alleles-1)]){
      push @all_alt_freqs, (sprintf "%.2f", $all_alt_cnt/$all_total);
    }
    print OUT join(',', @all_alt_freqs)."\t";
    print OUT "\n";
  }
}


######
#SUBS
######
sub check_and_process_input{
  #check input options
  if($tabix and ! -x $tabix){
    die "The specified tabix, $tabix, is not an executable.\n";
  }
  #default to this installation if not specified at the command line
  $tabix = "/nfs/1000g-work/G1K/work/bin/tabix/tabix" if (!$tabix);

  if(! -e $vcf){
    die "VCF file $vcf does not exist.\n";
  }
  if(! -e $sample_panel){
    die "Sample panel $sample_panel does not exist.\n";
  }
  if($out_dir and ! -d $out_dir){
    die "The directory $out_dir does not exist.\n";
  }

  my $panel_hash = parse_sample_panel($sample_panel);

  #if user specified populations, parse and check that they exist in the panel
  if($user_input_pop_string){
    @selected_pops = split ',', $user_input_pop_string;
    if($selected_pops[0] eq "ALL"){
      @selected_pops = keys(%pops);
    }
    else{
      foreach my $p (@selected_pops){
        die "Population $p is not listed in $sample_panel.\n" if(! exists $pops{$p});
      }
    }
  }
  else{
    #assume user wants all populations
    @selected_pops = keys(%pops);
  }

  if($region){
    die "Region should be in format chr:start-end, i.e. 1:1-1000.\n" if($region !~ /\s*:[0-9]*-[0-9]*$/);
  }

  #put together command to get the VCF
  my $cmd;
  if($no_tabix){
    $vcf_cmd = "cat $vcf";
  }
  elsif($region){
    $vcf_cmd = "$tabix -h $vcf $region";
  }
  else{
    $vcf_cmd = "zcat $vcf";
  }
}

sub parse_sample_panel {
  my ($panel_file) = @_;

  my $fh;
  open($fh, "<", $panel_file);
  while (<$fh>) {
    chomp;
    my $sample_line = $_;
    my ($sample, $pop_in_panel) = split(/\t/, $sample_line);
    #add unless it looks like the header
    $sample_pop_hash{$sample} = $pop_in_panel unless $sample eq "sample";
    $pops{$pop_in_panel} = 1 unless $sample eq "sample";
  }
}

=pod

=head1 NAME

calculate_allele_frq_from_vcf_file.pl

=head1 SYNOPSIS

This script takes a VCF file, a matching sample panel file, a chromosomal region, population names, it then calculates population-wide allele 
frequency for sites within the chromosomal region defined.

When no population is specified, allele fequences will be calcuated for all populations in the VCF files, one at a time.

=head1 Dependency

	To run this script, you need to install the following software
	tabix: http://sourceforge.net/projects/samtools/files/tabix/
	vcftools: http://sourceforge.net/projects/vcftools/files/
	
=head1 Options

	-vcf			Input VCF file that contains genotype data for each samples; this file must be bgzipped and tabix indexed; if '-no_tabix' is used, the 
vcf file can be uncompressed and un-indexed.
	-sample_panel	A tab deliminated file lists mapping between sample and population (see example below) 
	-pop			Populations of interest; separated by ",".  This field can be null
	-region			chromosome region of interest, format is chr_number:start_end (optional). 
	-out_dir		Where temporary file and final output files should be written to. Default is current direcotry.
	-tabix			A path to tabix executable
	-vcftools_dir	A path to the vcftools base directory that contains vcftools executable and perl libraries
	-no_tabix		If the input VCF files is a pre-sliced VCF file containing a small number of sites, this option can be used so the vcf file doesn't have 
					to be tabix indexed or bgzipped.  This is to speed up run time for the web application.
	-help			Print this page when specified

=head1 EXAMPLE lines from a sample panel file. Only the first 2 columns are essential.

	HG00096	GBR	EUR
	HG00097	GBR	EUR
	HG00099	GBR	EUR
	HG00100	GBR	EUR

=head1 OUTPUT

The allele frequency of an user-specified population for sites within the user-specified chromosomal region is written to a file.  The headers of the 
output file are:

	CHR:		Chromosome
	POS:		Start position of the variant
	ID:		Identification of the variant
	REF:		Reference allele
	ALT:		Alternative allele
	TOTAL_CNT:	Total number of alleles in samples of the chosen population(s) 
	ALT_CNT:	Number of alternative alleles observed in samples of the chosen populations(s)
	FRQ:		Ratio of ALT_CNT to TOTAL_CNT


=head1 EXAMPLE

perl $ZHENG_RP/bin/calculate_allele_frq_from_vcf.pl \
-vcf /nfs/1000g-archive/vol1/ftp/phase1/analysis_results/integrated_call_sets/ALL.chr22.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz \
-out_dir . \
-sample_panel /nfs/1000g-archive/vol1/ftp/phase1/analysis_results/integrated_call_sets/integrated_call_samples.20101123.ALL.panel \
-region 22:17000000-17005000 \
-pop CEU,FIN \

OR (when the input VCF file is pre-sliced small size file)

perl $ZHENG_RP/bin/calculate_allele_frq_from_vcf.pl \
-vcf ALL.chr22_17000000_17005000.test.vcf \
-out_dir ~ \
-sample_panel /nfs/1000g-archive/vol1/ftp/phase1/analysis_results/integrated_call_sets/integrated_call_samples.20101123.ALL.panel \
-pop CEU,FIN \
-no_tabix 
