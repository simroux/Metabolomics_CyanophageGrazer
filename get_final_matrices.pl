#!/usr/bin/perl
use strict;
use autodie;
use Getopt::Long;
use Text::ParseWords;
my $h='';
my $code='';
my $tag_ph=0;
GetOptions ('help' => \$h, 'h' => \$h, 'ph'=>\$tag_ph);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
	print "# Script to get the final matrices after removing outlier samples and correlated compounds
# Arguments :
# -ph if working from peak height (default: peak area)\n";
	die "\n";
}

my $type="Area";
my $code="pa";
if ($tag_ph==1){
	$type="Height";
	print "$type\n";
      $code="ph";
}
my @tab_type=("media","pellet_cnorm");

## Load clean list of compounds
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
		$tab[0].="_MSMS";
            print $tab[0]." is an internal standard\n";
            $check_std{$tab[0]}=1;
      }
}
close $tsv;

## Load samples that need to be removed
my $cor_samples="Correlated_samples.tsv";
my %remove_samples;
open my $tsv,"<",$cor_samples;
while(<$tsv>){
      chomp($_);
      if ($_=~/^#/){next;}
      my @tab=split("\t",$_);
      $remove_samples{$tab[0]}{$tab[1]}=1;
      $tab[0]=~s/_pa$//;
      $tab[0]=~s/_ph$//;
      $remove_samples{$tab[0]}{$tab[1]}=1;
      $tab[0]=~s/M/m/;
      $tab[0]=~s/P/p/;
      $remove_samples{$tab[0]}{$tab[1]}=1;
	if ($tab[0]=~/^pellet/){
		$tab[0]=~s/pellet/pellet_cnorm/;
	      $remove_samples{$tab[0]}{$tab[1]}=1;
		print "Removing $tab[1] from $tab[0]\n";
	}
}
close $tsv;

## Load compounds that need to be removed
my $corr_compounds="Correlated_compounds_decision_clean.tsv";
my %remove_compound;
open my $tsv,"<",$corr_compounds;
while(<$tsv>){
      chomp($_);
      if ($_=~/^#/){next;}
      my @tab=split("\t",$_);
      $remove_compound{$tab[0]}{$tab[2]}=1;
      $tab[0]=~s/M/m/;
      $tab[0]=~s/P/p/;
      $remove_compound{$tab[0]}{$tab[2]}=1;
	if ($tab[0]=~/^pellet/){
		$tab[0]=~s/pellet/pellet_cnorm/;
	      $remove_compound{$tab[0]}{$tab[2]}=1;
		print "Removing $tab[2] from $tab[0]\n";
	}
}
close $tsv;

my %seen;
my %store;
my %columns;
my %check_column;
foreach my $type (@tab_type){ ## For each type (Media and Pellet)
      my $in_file="Clean_data/All_".$type."_".$code.".csv";
      my $out_file="Selected_data/All_".$type."_".$code."_selected.csv";
      my $out_file_std="Selected_data/All_".$type."_".$code."_selected_std.csv";
      print "### $type - $in_file / $out_file / $out_file_std\n";
      ## Load clean data
      %store=();
      %columns=();
      %seen=();
      %check_column=();
      open my $csv,"<",$in_file;
      while(<$csv>){
            chomp($_);
            my @tab=split(",",$_);
            if ($tab[0] eq "Treatment"){
                  for (my $i=1;$i<=$#tab;$i++){
                        my $j=$i+1;
                        $columns{$j}=$tab[$i];
                  }
            }
            else{
                  if ($remove_samples{$type}{$tab[0]}==1){
                        print "We remove $tab[0]\n";
                        $seen{$tab[0]}=1;
                  }
                  else{
                        $store{$tab[0]}{"treatment"}=$tab[1];
                        for (my $i=2;$i<=$#tab;$i++){
                              my $name=$columns{$i};
                              if ($name=~/^\d/){$name="X".$name;}
                              $name=~s/-/./g;
                              if ($remove_compound{$type}{$name}==1){
                                    print "We remove $columns{$i}\n";
                                    $seen{$name}=1;
                              }
                              elsif($check_std{$columns{$i}}){
						print $columns{$i}." is an internal standard\n";
                                    $check_column{"std"}{$columns{$i}}=1;
                                    $check_column{"ok"}{$columns{$i}}=1;
                                    $store{$tab[0]}{$columns{$i}}=$tab[$i];
                              }
                              else{
                                    $check_column{"ok"}{$columns{$i}}=1;
                                    $store{$tab[0]}{$columns{$i}}=$tab[$i];
                              }
                        }
                  }
            }
      }
      close $csv;
      ## Write selected in a format usable for the next R script
      open my $s1,">",$out_file;
      my @t_compound=sort keys %{$check_column{"ok"}};
      print $s1 "Treatment,".join(",",@t_compound)."\n";
      foreach my $sample (sort keys %store){
            my $line=$sample.",".$store{$sample}{"treatment"};
            foreach my $compound (@t_compound){
                  $line.=",".$store{$sample}{$compound};
            }
            print $s1 $line."\n";
      }
      close $s1;
      ## Write standard
      open my $s1,">",$out_file_std;
      my @t_compound=sort keys %{$check_column{"std"}};
      print $s1 "Treatment,".join(",",@t_compound)."\n";
      foreach my $sample (sort keys %store){
            my $line=$sample.",".$store{$sample}{"treatment"};
            foreach my $compound (@t_compound){
                  $line.=",".$store{$sample}{$compound};
            }
            print $s1 $line."\n";
      }
      close $s1;
      ## check that we removed everything
      foreach my $sample (sort keys %{$remove_samples{$type}}){
            if ($seen{$sample}==1){print "All good\n";}
            else{
                  print "We did not remove $sample ??\n";
                  <STDIN>;
            }
      }
      foreach my $compound (sort keys %{$remove_compound{$type}}){
            if ($seen{$compound}==1){print "All good\n";}
            else{
                  print "We did not remove $compound ??\n";
                  <STDIN>;
            }
      }
}
