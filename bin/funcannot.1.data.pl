#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use FuhaoPerl5Lib::MiscKit qw/KeggJosnLoad/;
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 --input my.fa [Options]
Version: v20210127

Requirements:
    Modiles: Getopt::Long
             FuhaoPerl5Lib::MiscKit

Descriptions:
    Convert eggnog output to AnnotationForge datasets

Options:
    --help|-h
        Print this help/usage;
    --prefix|-p   <STR> Output prefix, default: MyOut

#  For eggnog files
    --input|-i    <Files> Eggnog output files ending with ".annotations"
                          Could be a list separated by comma or
                          specify multiple times or both
    --ignoreegg   <INT>   Ignore first INT comment lines, default:4
                          Will inore lines starting with '# ' or '#query_name'
    --rmeggtsnum          Trim transcript number ending with .\\d+\$
  
#  For custom GOs
    --gos|-g      <Files> 2-column data: 
                              column1: gene/transcript ID
                              column2: GOs sep by comma. example: GO:000000,GO000001
    --ignorego    <INT>   Ignore first INT comment lines, default:0
    --rmgotsnum           Trim transcript number ending with .\\d+\$
  
#  For custom KOs
    --kos|-k      <Files> 2-column data: 
                              column1: gene/transcript ID
                              column2: GOs sep by comma. example: K00001,K00002
    --ignoreko    <INT>   Ignore first INT comment lines, default:0
    --rmkotsnum           Trim transcript number ending with .\\d+\$
    --keggjson            KEGG Orthology file, downloaded from 
                          https://www.genome.jp/kegg-bin/show_brite?ko00001.keg

    --verbose
        Detailed output for trouble-shooting;
    --version|v!
        Print current SCRIPT version;

Example:
    perl $0 

Author:
    Fu-Hao Lu
    Professor, PhD
    State Key Labortory of Crop Stress Adaptation and Improvement
    College of Life Science
    Jinming Campus, Henan University
    Kaifeng 475004, P.R.China
    E-mail: lufuhao\@henu.edu.cn
EOH
###HELP ends#########################################################
die USAGE unless @ARGV;



###Receving parameter################################################
my ($help, $verbose, $debug, $version);
my ($input, $output);
my ($num_ignored_line1, $num_ignored_line2, $num_ignored_line3)=(4, 0, 0);
my ($transcript_num1, $transcript_num2, $transcript_num3)= (0,0,0);
my $prefix="MyOut";
my @eggnogfiles=();
my @gofiles=();
my @kofiles=();
my $kegg_json="";

GetOptions(
	"help|h!" => \$help,
	"input|i:s" => \@eggnogfiles,
	"gos|g:s" => \@gofiles,
	"kos|k:s" => \@kofiles,
	"prefix|p:s" => \$prefix,
	"ignoreegg:i" => \$num_ignored_line1,
	"ignorego:i" => \$num_ignored_line2,
	"ignoreko:i" => \$num_ignored_line3,
	"rmeggtsnum!" => \$transcript_num1,
	"rmgotsnum!" => \$transcript_num2,
	"rmkotsnum!" => \$transcript_num3,
	"keggjson:s" => \$kegg_json,
	"debug!" => \$debug,
	"verbose!" => \$verbose,
	"version|v!" => \$version) or die USAGE;
($help or $version) and die USAGE;



### Defaults ########################################################
$debug=0 unless (defined $debug);
$verbose=0 unless (defined $verbose);
my %gene2info=();
my %gene2go=();
my %gene2ko=();
my %gene2pathway=();



### input and output ################################################
@eggnogfiles=split(/,/,join(',',@eggnogfiles));
@gofiles=split(/,/,join(',',@gofiles));
@kofiles=split(/,/,join(',',@kofiles));
my $ko2pathway={};
if (scalar(@kofiles)>0) {
	if ($kegg_json ne "" and -s $kegg_json) {
		$ko2pathway=KeggJosnLoad($kegg_json);
	}
	else {
		print STDERR "Error: please provide a KEGG json file to convert KO to pathway\n";
		exit 100;
	}
}



### Main ############################################################
foreach my $x (@eggnogfiles) {
	my $linenum=0;
	open (INPUT, "<", $x) || die "Error: can not open Eggnog input file: $x\n";
	while (my $line=<INPUT>) {
		$linenum++;
		if ($linenum<=$num_ignored_line1) {
			next;
		}
		next if ($line=~/^#\s+/ or $line=~/^#query_name/);
		my @arr=();
		@arr=split(/\t/, $line);
		unless (scalar(@arr)==22) {
			print STDERR "Warnings: colnum !=22: $line in file $x, ignored\n";
			next;
		}
		if ($transcript_num1) {
			$arr[0]=~s/\.\d+$//;
		}
		#gene2info
		unless (exists $gene2info{$arr[0]} and exists ${$gene2info{$arr[0]}}{$arr[1]} and ${$gene2info{$arr[0]}}{$arr[1]}>$arr[3]) {
			$gene2info{$arr[0]}{$arr[1]}=$arr[3];
		}
		#%gene2go
		my @gos=();
		@gos=split(/,/, $arr[6]);
		foreach my $i (@gos) {
			$i=~s/^\s+//;$i=~s/\s+$//;
			if ($i=~/^GO:/i) {
#				$i=~s/^GO://i;
				$gene2go{$arr[0]}{$i}++;
			}
		}
		@gos=();
		#%gene2ko
		my @kos=();
		@kos=split(/,/, $arr[8]);
		foreach my $i (@kos) {
			$i=~s/^\s+//;$i=~s/\s+$//;
			if ($i=~/^ko:/i) {
				$i=~s/^ko://i;
				$gene2ko{$arr[0]}{$i}++;
			}
		}
		@kos=();
		#%gene2pathway
		my @pathway=();
		@pathway=split(/,/, $arr[9]);
		foreach my $i (@pathway) {
			$i=~s/^\s+//;$i=~s/\s+$//;
			if ($i=~/^ko/i) {
#				$i=~s/^ko://i;
				$gene2pathway{$arr[0]}{$i}++;
			}
		}
		@pathway=();
	}
	close INPUT;
}


my $target="uknGO000000001";
foreach my $y (@gofiles) {
	my $linenum=0;
	open (INPUT2, "<", $y) || die "Error: can not open Eggnog input file: $y\n";
	while (my $line=<INPUT2>) {
		$linenum++;
		if ($linenum<=$num_ignored_line2) {
			next;
		}
		my @arr=();
		@arr=split(/\t/, $line);
		unless (scalar(@arr)==2) {
			print STDERR "Warnings: colnum !=2: $line in file $y, ignored\n";
			next;
		}
		if ($transcript_num2) {
			$arr[0]=~s/\.\d+$//;
		}
		#gene2info
		unless (exists $gene2info{$arr[0]} and exists ${$gene2info{$arr[0]}}{$target}) {
			$gene2info{$arr[0]}{$target++}=0;
		}
		#%gene2go
		my @gos=();
		@gos=split(/,/, $arr[1]);
		foreach my $i (@gos) {
			$i=~s/^\s+//;$i=~s/\s+$//;
			if ($i=~/^GO:/i) {
#				$i=~s/^GO://i;
				$gene2go{$arr[0]}{$i}++;
			}
		}
		@gos=();
	}
	close INPUT2;
}

$target="uknKO000000001";
foreach my $z (@kofiles) {
	my $linenum=0;
	open (INPUT3, "<", $z) || die "Error: can not open KO input file: $z\n";
	while (my $line=<INPUT>) {
		$linenum++;
		if ($linenum<=$num_ignored_line3) {
			next;
		}
		my @arr=();
		@arr=split(/\t/, $line);
		next if (scalar(@arr)>=2 and defined $arr[1] and $arr[1] ne "");
		unless (scalar(@arr)==2) {
			print STDERR "Warnings: colnum !=2: $line in file $z, ignored\n";
			next;
		}
		if ($transcript_num3) {
			$arr[0]=~s/\.\d+$//;
		}
		#gene2info
		unless (exists $gene2info{$arr[0]} and exists ${$gene2info{$arr[0]}}{$target}) {
			$gene2info{$arr[0]}{$target++}=0;
		}
		#%gene2ko
		my @kos=();
		@kos=split(/,/, $arr[1]);
		foreach my $i (@kos) {
			$i=~s/^\s+//;$i=~s/\s+$//;
			if ($i=~/^K/i) {
				$gene2ko{$arr[0]}{$i}++;
				if (exists ${$ko2pathway}{$i}) {
					foreach my $ptw (keys %{${$ko2pathway}{$i}}) {
						$gene2pathway{$arr[0]}{$ptw}++;
					}
				}
			}
		}
		@kos=();
	}
	close INPUT3;
}



### write gene2info
my $file_gene2info="$prefix.gene_info.tab";
my $file_gene2go="$prefix.gene2go.tab";
my $file_gene2ko="$prefix.gene2ko.tab";
my $file_gene2pathway="$prefix.gene2pathway.tab";
open (INFO, " > ", $file_gene2info) || die "Error: can not write to geneinfo: $file_gene2info\n";
foreach my $a (sort keys(%gene2info)) {
	my $b=""; my $max=0;
	foreach my $c (sort keys(%{$gene2info{$a}})) {
		unless ($b=~/\S+/) {
			$b=$c; $max=$gene2info{$a}{$c};
		}
		elsif($max=$gene2info{$a}{$c}>$max) {
			$b=$c; $max=$gene2info{$a}{$c};
		}
	}
	print INFO $a, "\t", $b, "\n";
}
close INFO;
open (GO, " > ", $file_gene2go) || die "Error: can not write to gene2go: $file_gene2go\n";
foreach my $a (sort keys(%gene2go)) {
	foreach my $c (sort keys(%{$gene2go{$a}})) {
		print GO $a, "\t", $c, "\tIEA\n";
	}
	
}
close GO;
open (KO, " > ", $file_gene2ko) || die "Error: can not write to gene2ko: $file_gene2ko\n";
foreach my $a (sort keys(%gene2ko)) {
	foreach my $c (sort keys(%{$gene2ko{$a}})) {
		print KO $a, "\t", $c, "\n";
	}
}
close KO;
open (PATHWAY, " > ", $file_gene2pathway) || die "Error: can not write to file_gene2pathway: $file_gene2pathway\n";
foreach my $a (sort keys(%gene2ko)) {
	foreach my $c (sort keys(%{$gene2ko{$a}})) {
		print PATHWAY $a, "\t", $c, "\n";
	}
}
close PATHWAY;


#####################################################################
###                         sub functions                         ###
#####################################################################


1;
