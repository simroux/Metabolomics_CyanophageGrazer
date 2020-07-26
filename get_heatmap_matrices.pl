#!/usr/bin/perl
use strict;
use autodie;
use Getopt::Long;
my $h='';
my $code='';
GetOptions ('help' => \$h, 'h' => \$h, 'c=s'=>\$code);
if ($h==1 || $code eq ""){ # If asked for help or did not set up any argument
	print "# Script to get an input matrix for pheatmap
# Arguments :
# -c : code of the csv files storing the FC (e.g. Media_ph_s or Pellet_cnorm_ph_s)\n";
	die "\n";
}

my $search=$code."_C*-vs-C_All_log*";
my @list=<Diff_detections_2/$search.csv>;

## Load compounds
my %check_std;
my $file_compound="Full_compound_list.tsv";
open my $tsv,"<",$file_compound;
while(<$tsv>){
      chomp($_);
      if ($_=~/^#/){next;}
      my @tab=split("\t",$_);
      if ($tab[2] eq "Yes"){
            print $tab[0]." is an internal standard\n";
            $check_std{$tab[0]}=1;
		my $new_name=$tab[0];
		$tab[0].="_MSMS";
            print $tab[0]." is an internal standard\n";
            $check_std{$tab[0]}=1;
		$new_name=~s/-/\./g;
		$check_std{$new_name}=1;
		$new_name.="_MSMS";
		$check_std{$new_name}=1;
      }
}
close $tsv;

## Load custom cutoffs for specific time points / conditions (based on internal standards)
my %custom_cutoff;
my $custom_cutoff_list="Special_cutoff.tsv";
open my $tsv,"<",$custom_cutoff_list;
while(<$tsv>){
	chomp($_);
	my @tab=split("\t",$_);
	if ($tab[3] ne ""){$custom_cutoff{$tab[0]}{$tab[1]}{$tab[2]}{"upper"}=$tab[3];}
	if ($tab[4] ne ""){$custom_cutoff{$tab[0]}{$tab[1]}{$tab[2]}{"lower"}=$tab[4];}
}
close $tsv;


my %ordered_sample=("T1"=>1,"T2"=>2,"T3"=>3,"T4"=>4,"T5"=>5,"T6"=>6,"All"=>7);

my %check_sample;
my %store;
my $tag=0;
foreach my $file (@list){
	## For each condition (CG, CP, CPG), we write a matrix including only compounds with at least 1 siginficant time point
      %store=();
      %check_sample=();
      $tag=0;
	print "Reading $file\n";
      $file=~/.*_([^_]+)-vs-C_All_(log\d+)/;
      my $treatment=$1;
      my $log_transfo=$2;
	open my $csv,"<",$file;
	while(<$csv>){
		chomp($_);
		$_=~s/\"//g;
		my @tab=split(",",$_);
		if ($tab[0] eq "Data"){next;}
		# my $s=$tab[4]."_T".$tab[3];
            my $s="T".$tab[2];
            if ($tab[2] eq "All"){$s=$tab[2];}
		$check_sample{$s}=1;
		print $tab[1]."\n";
		if ($tab[1]=~/^X\d/){
			$tab[1]=~s/^X//;
		}
		$store{$tab[1]}{$s}{"FC"}=$tab[4];
		$store{$tab[1]}{$s}{"FDR"}=$tab[6];
	}
	close $csv;
      my $out_file="Heatmap_input/".$code."_".$treatment."_".$log_transfo."_for_pheatmap.csv";
      my $out_file_forclust="Heatmap_input/".$code."_".$treatment."_".$log_transfo."_for_clustering.csv";
      my $out_file_std="Heatmap_input/".$code."_".$treatment."_".$log_transfo."_std_for_pheatmap.csv";
      my $out_file_std_forclust="Heatmap_input/".$code."_".$treatment."_".$log_transfo."_std_for_clustering.csv";
      print "Writing $out_file and $out_file_std along with $out_file_forclust and $out_file_std_forclust\n";
      my @t_s=sort {$ordered_sample{$a} <=> $ordered_sample{$b} or $a <=> $b} keys %check_sample;
      my $f_l=join(",",@t_s);
      open my $s1,">",$out_file;
      print $s1 "$f_l\n";
      open my $s1_c,">",$out_file_forclust;
      print $s1_c "$f_l\n";
      open my $s2,">",$out_file_std;
      print $s2 "$f_l\n";
      open my $s2_c,">",$out_file_std_forclust;
      print $s2_c "$f_l\n";
      $tag=0;
      foreach my $compound (sort keys %store){
      	my $line=$compound;
            my $line_c=$compound;
		my $line_std=$compound;
            $tag=0;
      	foreach my $samp (@t_s){
                  $line_c.=",".$store{$compound}{$samp}{"FC"};
                  if ($log_transfo eq "log2"){ ## Two different thresholds, one for log2, one for log10
                        if ($store{$compound}{$samp}{"FDR"}<=0.05 && abs($store{$compound}{$samp}{"FC"})>=1){}
                        else{$store{$compound}{$samp}{"FC"}="NA";}
                  }
                  elsif ($log_transfo eq "log10"){
                        if ($store{$compound}{$samp}{"FDR"}<=0.05 && abs($store{$compound}{$samp}{"FC"})>=0.1){} ## Lower threshold for log10 transformation
                        else{$store{$compound}{$samp}{"FC"}="NA";}
                  }
      		$line_std.=",".$store{$compound}{$samp}{"FC"};
			if ($log_transfo eq "log10"){
				# print $compound."\t".$code."\t".$treatment."\t".$samp."\t".$store{$compound}{$samp}{"FC"}."\n";
				if (defined($custom_cutoff{$code}{$treatment}{$samp}) &&  ($store{$compound}{$samp}{"FC"} ne "NA")){
					print "we have a custom cutoff for $code - $treatment - $samp -> ".$custom_cutoff{$code}{$treatment}{$samp}{"upper"}."\t".$custom_cutoff{$code}{$treatment}{$samp}{"lower"}."\n";
					# <STDIN>;
					if (defined($custom_cutoff{$code}{$treatment}{$samp}{"upper"}) && ($store{$compound}{$samp}{"FC"}>=0.1)){
						if (($custom_cutoff{$code}{$treatment}{$samp}{"upper"} > $store{$compound}{$samp}{"FC"})){
							print "$compound - $samp used to be significant (upper), but is not anymore in the new cutoff\n";
							$store{$compound}{$samp}{"FC"}="NA";
							# <STDIN>;
						}
					}
					if (defined($custom_cutoff{$code}{$treatment}{$samp}{"lower"}) && ($store{$compound}{$samp}{"FC"}<=-0.1)){
						if (($custom_cutoff{$code}{$treatment}{$samp}{"lower"} < $store{$compound}{$samp}{"FC"})){
							print "$compound - $samp used to be significant (lower), but is not anymore in the new cutoff\n";
							$store{$compound}{$samp}{"FC"}="NA";
						}
						else{
							print $store{$compound}{$samp}{"FC"}." is ok \n";
						}
						# <STDIN>;
					}
				}
			}
                  if ($log_transfo eq "log2"){ ## Two different thresholds, one for log2, one for log10
                        if ($store{$compound}{$samp}{"FDR"}<=0.05 && abs($store{$compound}{$samp}{"FC"})>=1){$tag=1;}
                  }
                  elsif ($log_transfo eq "log10"){
                        if ($store{$compound}{$samp}{"FDR"}<=0.05 && abs($store{$compound}{$samp}{"FC"})>=0.1){$tag=1;} ## Lower threshold for log10 transformation
                        else{$store{$compound}{$samp}{"FC"}="NA";}
                  }
      		$line.=",".$store{$compound}{$samp}{"FC"};
      	}
            if ($tag==1){
                  print $s1 $line."\n";
                  print $s1_c $line_c."\n";
            }
            if ($check_std{$compound}==1){
                  print $s2 $line_std."\n";
                  print $s2_c $line_c."\n";
                  ## This is the matrix including all internal standards, so we take all (and only) these internal standards
            }
		# print $compound."\t===".$check_std{$compound}."---\n";
      }
      close $s1;
      close $s1_c;
      close $s2;
      close $s2_c;
}
