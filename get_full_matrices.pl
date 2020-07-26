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
	print "# Script to get the input matrices for the CpG Metabolomics analysis
# Arguments :
# -ph if working from peak height (default: peak area)\n";
	die "\n";
}

## Decide if we'll use peak area or height
my $type="Area";
if ($tag_ph==1){
	$type="Height";
	print "$type\n";
}


my $out_file_media="";
my $out_file_pellet="";
my $out_file_pellet_cnorm="";
my $out_file_blank="";
my @tab_files;

## Adjust input files based on the type of metric we use
my @tab_types=("Media","Media","Pellet","Pellet");
if ($type eq "Area"){
	@tab_files=("Raw_data/20170215_RL_Sullivan_media_pos/sheets/peak_area.tab","Raw_data/20170414_RL_Sullivan_media_neg_v4/sheets/peak_area.tab","Raw_data/20161213_Sullivan_pellet_pos/sheets/peak_area.tab","Raw_data/20170411_RL_Sullivan_pellet_NEG_v3-3/sheets/peak_area.tab");
	$out_file_media="Clean_data/All_media_pa.csv";
	$out_file_pellet="Clean_data/All_pellet_pa.csv";
	$out_file_pellet_cnorm="Clean_data/All_pellet_cnorm_pa.csv";
	$out_file_blank="Clean_data/All_blanks_and_others_pa.csv";
}
elsif($type eq "Height"){
	@tab_files=("Raw_data/20170215_RL_Sullivan_media_pos/sheets/peak_height.tab","Raw_data/20170414_RL_Sullivan_media_neg_v4/sheets/peak_height.tab","Raw_data/20161213_Sullivan_pellet_pos/sheets/peak_height.tab","Raw_data/20170411_RL_Sullivan_pellet_NEG_v3-3/sheets/peak_height.tab");
	$out_file_media="Clean_data/All_media_ph.csv";
	$out_file_pellet="Clean_data/All_pellet_ph.csv";
	$out_file_pellet_cnorm="Clean_data/All_pellet_cnorm_ph.csv";
	$out_file_blank="Clean_data/All_blanks_and_others_ph.csv";
}


## Load msms info
my %valid_msms;
my $file_compounds="Full_compound_list.tsv";
open my $tsv,"<",$file_compounds;
while(<$tsv>){
	chomp($_);
	if ($_=~/^#/){next;}
	my @tab=split("\t",$_);
	if ($tab[5] eq "Yes" || $tab[9] eq "Yes"){
		$valid_msms{"Media"}{$tab[0]}=1;
	}
	if ($tab[13] eq "Yes" || $tab[17] eq "Yes"){
		$valid_msms{"Pellet"}{$tab[0]}=1;
	}
}
close $tsv;

my %check_comp;
my %store_blank;
my %c_name;
my %store;
my %check;
my %s_to_t;
my $tag=0;

## For each file
for (my $i=0;$i<=$#tab_files;$i++){
	my $file_value=$tab_files[$i];
	my $type=$tab_types[$i];
	print "Reading $file_value\n";
	open my $csv,"<",$file_value;
	%c_name=();
	$tag=0;
	while(<$csv>){
		chomp($_);
		my @tab=split("\t",$_);
		if ($tab[0] eq "" || $tab[0] eq "group"){next;}
		if ($tag==0){
			### This is the first line, we take the list of samples from this line
			$tag=1;
			if ($tab[0] ne "file"){
				print "?! $file_value is weird -> $_\n";
				<STDIN>;
			}
			for (my $i=1;$i<=$#tab;$i++){
				my $sample=$tab[$i];
				if ($sample=~/.*_5uM.*/){ ## This is a blank sample
					$c_name{$i}="blank_".$sample;
					print "$tab[$i] - sample $sample ==> other $c_name{$i}\n";
				}
				elsif ($sample=~s/.*_Qex_(C_*P*_*G*_*\d+[^_]+)//){ ## We decode the sample name
					$sample=$1;
					$sample=~s/_//g;
					$c_name{$i}=$sample;
					$sample=~/(.*)[ABCD]$/;
					$s_to_t{$sample}=$1;
					print "$tab[$i] - sample $sample ==> $s_to_t{$sample}\n";
				}
				elsif ($sample=~/.*blank.*/i || $sample=~/.*QCmix.*/ || $sample=~/.*DUBA.*/){ ## This is a blank sample
					$c_name{$i}="blank_".$sample;
					print "$tab[$i] - sample $sample ==> blank $c_name{$i}\n";
				}
				elsif ($sample=~s/.*_HILIC_.*_(C\+*P*\+*G*\_*\d+[^_]+)//){ ## We decode the sample name
					$sample=$1;
					$sample=~s/\+//g;
					$sample=~s/_//g;
					$sample=~s/\.h\d//g;
# 					print "sample $tab[$i] -> $sample\n";
					$c_name{$i}=$sample;
					$sample=~/(.*)[ABCD]$/;
					$s_to_t{$sample}=$1;
					print "$tab[$i] - sample $sample ==> $s_to_t{$sample}\n";
				}
				elsif ($sample=~/.*M_*\d.*/i || $sample=~/.*CCMP3375.*/){
					$c_name{$i}="media";
					print "$tab[$i] - sample $sample ==> other $c_name{$i}\n";
				}
				else{
					print "Weird sample $tab[$i] ?\n";
					<STDIN>;
				}
			}
		}
		else{
			### This is not the first line, so this should correspond to the values we want for a given compound (in first column)
			my $compound=$tab[0];
			print "$compound\n";
			if ($compound=~/(.*)_positive_.*/){$compound=$1;}
			elsif ($compound=~/(.*)_negative_.*/){$compound=$1;}
			else{print "$compound ?? \n";<STDIN>;}
			my $tag_enough=0;
			for (my $i=1;$i<=$#tab;$i++){ ## We verify that we have enough relevant samples for this compound
				if ($tab[$i]>=10000){$tag_enough++;}
			}
			if ($tag_enough>=2){
				## This is an interesting compound with enough samples > 0, now we actually store the values
				for (my $i=1;$i<=$#tab;$i++){
					my $sample=$c_name{$i};
					if ($sample=~/^blank/ || $sample=~/^media/){
						$store_blank{$sample}{$compound}=$tab[$i];
					}
					elsif (defined($store{$sample}{$type}{$compound})){
						print "?! Already a value for $compound and sample $sample in $type ->  $store{$sample}{$type}{$compound} == ";
						## For compounds detected multiple times (e.g. positive and negative mode) we use the highest value
						if ($store{$sample}{$type}{$compound}<$tab[$i]){
							print "$tab[$i] vs $store{$sample}{$type}{$compound}, we replace\n";
							$store{$sample}{$type}{$compound}=$tab[$i];
						}
						print "\n";
					}
					else{
						$store{$sample}{$type}{$compound}=$tab[$i];
					}
					$check_comp{$type}{$compound}=1;
					$check_comp{"all"}{$compound}=1;
				}
			}
			else{
				print "############# $compound ==> WE DON'T WANT BECAUSE IT'S NEVER REALLY DETECTED\n";
				<STDIN>;
			}
		}
	}
	close $csv;
	print "Now we have processd $file_value\n";
}

print "#######################################################\n";
print "Everything read, now writing output\n";
## We first write the output file corresponding to the Media samples
print "Writing $out_file_media\n";
open my $s1,">",$out_file_media;
my $f_line="Treatment,";
foreach my $compound (sort keys %{$check_comp{"Media"}}){
	my $compound_name=$compound;
	if ($valid_msms{"Media"}{$compound}==1){
		$compound_name.="_MSMS";
		print "Compound $compound was detected confidently in Media\n";
	}
	$f_line.=$compound_name.",";
}
chop($f_line);
print $s1 "$f_line\n";
foreach my $sample (sort keys %store){
	my $line=$sample.",".$s_to_t{$sample};
	my $total=0;
	foreach my $compound (sort keys %{$check_comp{"Media"}}){$line.=",".$store{$sample}{"Media"}{$compound};$total+=$store{$sample}{"Media"}{$compound};}
	if ($total>0){print $s1 "$line\n";}
}
close $s1;

## Now we write the output file corresponding to the Pellet samples
print "Writing $out_file_pellet\n";
open my $s1,">",$out_file_pellet;
my $f_line="Treatment,";
foreach my $compound (sort keys %{$check_comp{"Pellet"}}){
	my $compound_name=$compound;
	if ($valid_msms{"Pellet"}{$compound}==1){
		$compound_name.="_MSMS";
		print "Compound $compound was detected confidently in Pellet\n";
	}
	$f_line.=$compound_name.",";
}
chop($f_line);
print $s1 "$f_line\n";
foreach my $sample (sort keys %store){
	my $line=$sample.",".$s_to_t{$sample};
	my $total=0;
	foreach my $compound (sort keys %{$check_comp{"Pellet"}}){$line.=",".$store{$sample}{"Pellet"}{$compound};$total+=$store{$sample}{"Pellet"}{$compound};}
	if ($total>0){print $s1 "$line\n";}
}
close $s1;

### And now we generate the output file for Pellet samples where values are normalized by the number of cells included in the pellet
print "#######################################################\n";
print "Make the cell normalization\n";
my $cell_count="Raw_data/CPG_cell_counts.tsv";
my %c_n;
open my $tsv,"<",$cell_count;
while(<$tsv>){
	chomp($_);
	my @tab=split("\t",$_);
	if ($tab[0] eq "Sample"){}
	else{
		$tab[0]=~s/_//g;
		$tab[0]=~s/\+//g;
		if ($tab[3] eq ""){$tab[3]=$tab[1]*40;}
		$c_n{$tab[0]}=$tab[3];
		print "$tab[0] -> $c_n{$tab[0]}\n";
	}
}
close $tsv;


open my $s1,">",$out_file_pellet_cnorm;
my @tab_comp=sort keys %{$check_comp{"Pellet"}}; # Listing all compounds detected in this type of sample
my $f_line="Treatment,";
foreach my $compound (@tab_comp){
	my $compound_name=$compound;
	if ($valid_msms{"Pellet"}{$compound}==1){
		$compound_name.="_MSMS";
		print "Compound $compound was detected confidently in Pellet\n";
	}
	$f_line.=$compound_name.",";
}
chop($f_line);
print $s1 $f_line."\n";
foreach my $sample  (sort keys %store){
	my $line=$sample.",".$s_to_t{$sample};
	my $total=0;
	if (!defined($c_n{$sample}) || $c_n{$sample}==0){
		print "no cell number for $sample ?\n";
		<STDIN>;
	}
	print "Normalizing sample $sample\n";
	foreach my $compound (@tab_comp){
		my $calc="";
		if (!defined($store{$sample}{"Pellet"}{$compound}) || $store{$sample}{"Pellet"}{$compound} eq ""){
			print "No info for $compound in $sample for Pellets\n";
		}
		else{
			$calc=$store{$sample}{"Pellet"}{$compound}/$c_n{$sample};
			$total++;
		}
		$line.=",".$calc;
	}
	print "$total compounds for sample $sample\n";
	if ($total>=0){print $s1 "$line\n";}
	else{print "$sample is removed because only $total compounds\n";}
}
close $s1;

### Finally, we generate a file with the values in the blank samples (for negative control)
print "Writing $out_file_blank\n";
open my $s1,">",$out_file_blank;
my $f_line="";
foreach my $compound (sort keys %{$check_comp{"all"}}){
	$f_line.=$compound.",";
}
chop($f_line);
print $s1 "$f_line\n";
foreach my $sample (sort keys %store_blank){
	my $line=$sample;
	foreach my $compound (sort keys %{$check_comp{"all"}}){$line.=",".$store_blank{$sample}{$compound};}
	print $s1 "$line\n";
}
close $s1;
