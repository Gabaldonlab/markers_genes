#!/usr/bin/perl
#
# Ktreedist
# Calculation of the minimum branch length distance (K tree score) between phylogenetic trees
#
# Copyright (C) 2007 Victor Soria-Carrasco & Jose Castresana
# Institute of Molecular Biology of Barcelona (IBMB), CSIC, Jordi Girona 18, 08034 Barcelona, Spain
# vscagr@ibmb.csic.es (VSC) & jcvagr@ibmb.csic.es (JC)
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# Requirements:
#	Perl v.5.8.x or greater (previous version have not been tested)
#
# Usage:
#  Ktreedist.pl <options> -rt <reference tree file> -ct <comparison tree/s file>
#     Options:
#          -s [<filename>]  Output file for comparison tree/s after scaling
#          -t [<filename>]  Output file for table of results
#          -p [<filename>]  Output file for partitions table
#          -r  Output symmetric difference (Robinson-Foulds)
#          -n  Output number of partitions in the comparison tree/s
#          -a  Equivalent to all options\n"
#
# For example:
# Ktreedist.pl -rt example_ref.tree -ct example_comp.tree -a
# or, if you have not made  Ktreedist.pl executable in a UNIX system:
# perl Ktreedist.pl -rt example_ref.tree -ct example_comp.tree -a

use warnings;
use strict;


# Ktreedist	version 1.0
# Last modified 25 June 2007
my $version = "1.0";
my $monthyear = "June 2007";


# //////////////////////////////////////////////////////////////////////////////
# //////////////////////////////  HIDDEN OPTIONS  //////////////////////////////
# //////////////////////////////////////////////////////////////////////////////

# The variable hard_polytomies is "off" by default. By default, zero branch length partitions count in the calculation of the symmetric difference.
# With hard_polytomies "on", zero branch length partitions are collapsed and treated as polytomies for this calculation.
my $hard_polytomies = "off";


# Debug mode is "off" by default. With debug mode "on" you do not need to enter -rt and -ct for the reference tree and
# comparison tree/s files, respectively, but combinations of parameters to be entered may be less flexible.
my $debugmode = "off";


# Safe mode is "off" by default. With safe mode "on" you are warned in case of overwriting files. When safe mode is "on"
# a new visible option appears to control overwriting:
#          -f  Overwrite output files
my $safemode = "off";


# Number of decimal places for branch lenghts in scaled trees
my $precision = 10;



# //////////////////////////////////////////////////////////////////////////////
# ///////////////////////////////  PRELIMINARS  ////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////

print "\nKtreedist version $version -  $monthyear";
print "\nCalculation of the minimum branch length distance (K tree score) between trees\n\n";

# Get command line arguments
my ($reftree, $comptrees, $scaled, $table, $partitions, $robinson, $npartitions, $overwrite) = get_cmd_line(@ARGV);

# DEFAULT DEFINITIONS, PATHS & FILENAMES
#-------------------------------------------------------------------------------
# Get paths & filenames
my ($filename_reftree, $filename_reftree_woext, $path_reftree) = get_path_filename($reftree);
my ($filename_comptrees, $filename_comptrees_woext, $path_comptrees) = get_path_filename($comptrees);

#ÊThis is really necessary in case of no full path
$reftree = "$path_reftree$filename_reftree" if ($reftree ne "");
$comptrees = "$path_comptrees$filename_comptrees" if ($comptrees ne "");

my $default_scaled = "$path_comptrees$filename_comptrees.scaled";
my $default_table = "$path_comptrees$filename_comptrees.tab";
my $default_partitions = "$path_comptrees$filename_comptrees.part";

# EXTRAS
$scaled = $default_scaled if ($scaled eq "yes");
$table = $default_table if ($table eq "yes");
$partitions = $default_partitions if ($partitions eq "yes");

# Do initial checkings and get trees if there are no errors
print "Checkings:\n";
my ($reftr, $comptr, $nametr, $txbl, $root, $bl, $t0,  $tr1) = initial_check ($reftree, $comptrees, $scaled, $table, $partitions);
$reftree = ($$reftr);
$comptrees = ($$comptr);
my %nametrees = (%$nametr);
my $taxblock = ($$txbl);
my $rooted = ($$root);
my $bk_line = ($$bl);
my $tree0 = ($$t0);
my @trees1 = (@$tr1);


# //////////////////////////////////////////////////////////////////////////////
# ///////////////////////////////////  MAIN  ///////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////

if (scalar(@trees1) == 1){
	print "\nK tree score calculation:\n'$filename_reftree' vs '$filename_comptrees'\n";
}
else{
	print "\nK tree score calculation:\n'$filename_reftree' vs '$filename_comptrees' (".scalar(@trees1)." trees)\n";
}

my ($claA, $claB, $brl, $sp) = get_partitions_tree($tree0, "yes");
my @clados0A = (@$claA);
my @clados0B = (@$claB);
my %brlen0 = (%$brl);
my @species0 = (@$sp);
my @clados0 = (@clados0A, @clados0B);

my $ntree = 0;
my @tab = ();
my @part = ();
my @scaled_trees = ();

my $nameref = "";
if (exists $nametrees{$reftree}{0}){ $nameref = $nametrees{$reftree}{0}; }
else{ $nameref = $filename_reftree_woext; }
my $rtd = "";
if ($rooted eq "yes"){ $rtd = "rooted";}
else{ $rtd = "unrooted";}

my $n_clados0 = scalar(@clados0)/2;
if ($hard_polytomies eq "on"){
	print "Zero branch length partitions, if present, are collapsed.\n";
	my $min=0;
	foreach my $br (values %brlen0){
		if ($br == 0){
			$min++;
		}
	}
	$n_clados0 = $n_clados0 - $min/2;
	$min = 0;
}

print "Reference tree ($nameref) is $rtd and has ".scalar(@species0)." tips and ".$n_clados0." partitions.\n";

print "Number of partitions on comparison tree/s are listed below.\n" if ($npartitions eq "yes");

my $initial_lines = "Comparison_tree    K-score   Scale_factor";
$initial_lines = $initial_lines."  Symm_difference" if ($robinson eq "yes");
$initial_lines = $initial_lines."  N_partitions" if ($npartitions eq "yes");

$initial_lines = "\n  $initial_lines  \n-";
my $len = (length($initial_lines)-3);
for (my $i=0; $i < $len; $i++){ $initial_lines = "-".$initial_lines."-";}
$initial_lines =~ s/-$/\n/g;
print "$initial_lines";

foreach my $tree1 (@trees1){
	my ($cla, $clb, $brl) = get_partitions_tree($tree1, "yes");
	my @clados1A = (@$cla);
	my @clados1B = (@$clb);
	my @clados1 = (@clados1A, @clados1B);
	my %brlen1 = (%$brl);
	my $n_clados1 = scalar (@clados1)/2;

	if ($hard_polytomies eq "on"){
		my $min=0;
		foreach my $br (values %brlen1){
			if ($br == 0){
				$min++;
			}
		}
		$n_clados1 = $n_clados1 - $min/2;		
	}
	
	# Get matrix for each pair of trees
	my @matrix = get_matrix (\@clados0A, \@clados0B, \@clados1A, \@clados1B, \%brlen0, \%brlen1);
	# Matrix array structure
	# First index: Index of partition
	# Second index:
	# 0				-	partition (tips)
	# 1				-	symmetric partition (tips)
	# 2				-	branch length partition reference tree
	# 3				-	branch length partition comparison tree
	# 4				- 	branch length partition reference tree "NONE"
	# 5				- 	branch length partition comparison tree "NONE"
	# After scale factor
	# 6				- 	branch length partition comparison tree after scaling
	# 7				- 	branch length partition comparison tree after scaling "NONE"


	# Get K-value
	my $K_value = get_K_value(@matrix);

	# Get K-score
	my ($ks, $rf, $ma) = get_K_score (\@matrix, \$K_value);
	my $K_score = ($$ks);
	my $drf = ($$rf);
	@matrix = (@$ma);

	# Scale comparison tree & output if -s was invoked
	my ($tr, $nwbrl, $sct) = modify_brlen_tree(\$scaled, \$tree1, \%brlen1, \$K_value);
	$tree1 = ($$tr);
	%brlen1 = (%$nwbrl);
	push (@scaled_trees, ($$sct));

	# //// OUTPUT DISPLAY ////
	my $output = "  ".output_display(	\$comptrees, \$ntree, \%nametrees, 
									\$K_score, \$K_value, \$robinson, \$drf, \$n_clados1);
	
	print "$output";

	# //// GENERATE ARRAYS FOR OUTPUT FILES ////

	push (@tab, get_table ( \$reftree, \$comptrees, \%nametrees, 
							\$filename_comptrees_woext,	\$ntree, \$npartitions, 
							\$K_score, \$K_value, \$robinson, \$drf,
							\$n_clados0, \$n_clados1))
			if ($table ne "no");
	
	push (@part, get_parts_table ( \$comptrees, \%nametrees, \$ntree, \@matrix)) 
			if ($partitions ne "no");

	$ntree++;
}

my $final_line = "";
for (my $i=0; $i < $len; $i++){ $final_line = $final_line."-"; }
$final_line = $final_line."\n";
print $final_line;


# //////////////////////////////////////////////////////////////////////////////
# ///////////////////////////////  OUTPUT  FILES ///////////////////////////////
# //////////////////////////////////////////////////////////////////////////////


print "\nFile saving:" if (($table ne "no") || ($partitions ne "no") || $scaled ne "no");
print "\n";

# Output scaled trees
if ($scaled ne "no"){
	print "Generating output file for comparison tree/s after scaling... ";
	open (FILE, ">$scaled") or die ("Cannot write to $scaled");
		foreach my $st (@scaled_trees){
			print FILE "$st\n";
		}
	close (FILE);
	newick2nexus (\$scaled, \$taxblock, \%nametrees, \$comptrees) if (exists ($nametrees{$comptrees}{0}));
	unix2macdos($bk_line, $scaled) if ($bk_line =~ m/^MAC|DOS$/);
	print "Done\n";
}

# 	Output results in a tabulated text file
if ($table ne "no"){
	print "Generating table of results... ";
	open (FILE, ">$table") or die ("Cannot write to $table");
		my $line = "Trees\tK-score\tScale_factor";
		$line = $line."\tSymm_difference" if ($robinson eq "yes");
		$line = $line."\tN_partitions_ref_tree\tN_partitions_comp_tree" if ($npartitions eq "yes");
		print FILE "$line\n";
		foreach my $l (@tab){
			print FILE "$l\n";
		}
	close (FILE);
	unix2macdos($bk_line, $table) if ($bk_line =~ m/^MAC|DOS$/);
	print "Done\n";
}


# 	Output partitions table
if ($partitions ne "no"){
	print "Generating table of partitions... ";
	open (FILE, ">$partitions") or die ("Cannot write to $partitions");
		print FILE "File generated by Ktreedist version $version\n\n";
		print FILE "Reference file: '$reftree'\n";
		print FILE "Comparison file: '$comptrees'\n";
		print FILE "Zero branch-length partitions, if present, are collapsed.\n" if ($hard_polytomies eq "on");
		print FILE "\n";
		print FILE "Meaning of each column:\n";
		print FILE "\tBrl_ref_tree: Branch length of this partition on the reference tree.\n";
		print FILE "\tBrl_comp_tree: Branch length of this partition on the original comparison tree.\n";
		print FILE "\tBrl_comp_tree_K: Branch length of this partition on the comparison tree after scaling.\n";
		print FILE "\tPartition: Tip nodes that constitute this partition.\n\n";
		foreach my $l (@part){
			print FILE "$l";
		}
	close(FILE);
	unix2macdos($bk_line, $partitions) if ($bk_line =~ m/^MAC|DOS$/);
	print "Done\n";
}
print "\n";

exit();


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------- SUBROUTINES ----------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# /// Get paths & filenames (with and without extension) ///
sub get_path_filename{
	my $filename = shift(@_);
	my $path = $filename;
	$filename =~ s/.+[\/|\\]//g;
	my $filename_woext = $filename;
	$filename_woext =~ s/\..+//g;
	$path =~ s/\/*$filename$//g;
	$path = "./".$path if ($path !~ m/^[A-Z]*\:*[\/|\.\/|\\]/);
	$path = $path."/" if ($path !~ m/[\/|\\]$/);

	return($filename, $filename_woext, $path);
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# /// Print usage howto ///
sub usage{
	my $wrong_option = shift(@_);
	print "WARNING: Option '$wrong_option' not recognized.\n\n" if (defined ($wrong_option));
	print "Usage:\n";
	if ($debugmode eq "on"){
		print " Ktreedist.pl <reference tree file> <comparison tree/s file> [<options>]\n";
	}
	else{
		print " Ktreedist.pl -rt <reference tree file> -ct <comparison tree/s file> [<options>]\n";
	}
	print "    Options:\n";
	print "         -t [<filename>]  Output file for table of results\n";
	print "         -p [<filename>]  Output file for table of partitions\n";
	print "         -s [<filename>]  Output file for comparison tree/s after scaling\n";
	print "         -r  Output symmetric difference (Robinson-Foulds)\n";
	print "         -n  Output number of partitions in the comparison tree/s\n";
 	print "         -f  Overwrite output files\n" if ($safemode eq "on");
	print "         -a  Equivalent to all options\n";
	print "         \n";
	exit();
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# /// Get command line parameters ///
sub get_cmd_line{
	my $get = "";
	my $reference_tree = "";
	my $comp_trees = "";
	my $scaled = "no";
	my $table = "no";
	my $partitions = "no";
	my $npartitions = "no";
	my $robinson = "no";
	my $overwrite = "yes";
	$overwrite = "no" if ($safemode eq "on");

	my @cmd = @_;

	if ($debugmode eq "on"){
# 		print "\nCaution: Debug mode is on.\n\n";
		$reference_tree = $ARGV[0] if (exists ($ARGV[0]));
		$comp_trees = $ARGV[1] if (exists ($ARGV[1]));

		if (scalar(@ARGV) == 2){
			@cmd = ();
			push (@cmd, "-r");
			push (@cmd, "-n");
		}
	}

	for (my $i=0; $i < scalar (@cmd); $i++){
		if ($cmd[$i] =~ m/^-[a|A]$/){
			$scaled = "yes";
			$table = "yes";
			$partitions = "yes";
			$robinson = "yes";
			$npartitions = "yes";
			$overwrite = "yes";
			next;
		}
		if ($cmd[$i] =~ m/^-[h|H]$/){
			usage();
		}
		if ($cmd[$i] =~ m/^-[s|S]$/){
			$get = "scaled";
			$scaled = "yes";
			next;
		}
		elsif ($cmd[$i] =~ m/^-[t|T]$/){
			$get = "table";
			$table = "yes";
			next;
		}
		elsif ($cmd[$i] =~ m/^-[p|P]$/){
			$get = "partitions";
			$partitions = "yes";
			next;
		}
		elsif ($cmd[$i] =~ m/^-[r|R]$/){
			$robinson = "yes";
			next;
		}
		elsif ($cmd[$i] =~ m/^-[n|N]$/){
			$npartitions = "yes";
			next;
		}
		elsif ($cmd[$i] =~ m/^-[r|R][t|T]$/){
			$get = "reference";
			next;
		}
		elsif ($cmd[$i] =~ m/^-[c|C][t|T]$/){
			$get = "comp";
			next;
		}
		elsif ($cmd[$i] =~ m/^-[f|F]$/){
			$overwrite = "yes";
			next;
		}

		#*****

		if ($get eq "scaled" && $scaled eq "yes"){
			$scaled = $cmd[$i];
			next;
		}
		if ($get eq "table" && $table eq "yes"){
			$table = $cmd[$i];
			next;
		}
		if ($get eq "partitions" && $partitions eq "yes"){
			$partitions = $cmd[$i];
			next;
		}
		elsif ($i<2 && $debugmode eq "on") {next;}
		if ($get eq "reference" && $reference_tree eq ""){
			$reference_tree = $cmd[$i];
			next;
		}
		elsif ($i<2 && $debugmode eq "on") {next;}
		if ($get eq "comp" && $comp_trees eq ""){
			$comp_trees = $cmd[$i] ;
			next;
		}

		if ($cmd[$i] !~ m/^-[a|h|s|t|c|f]$/i){
			usage($cmd[$i]);
		}
	}
	usage () if (($reference_tree eq "") || ($comp_trees eq ""));
	return ($reference_tree, $comp_trees, $scaled, $table, $partitions, $robinson, $npartitions, $overwrite);
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# This subroutine checks several parameters before starting tree comparison:
#	- if input and output files exist (and can be writen)
#	- kind of line break
#	- input files format (nexus needs internal conversion to newick)
#	- if all trees are rooted or unrooted
#	- if trees have branch length information
#	- if trees have same tip names
sub initial_check{
	my $reftree = shift(@_);
	my $comptrees = shift(@_);
	my $scaled = shift(@_);
	my $table = shift(@_);
	my $partitions = shift(@_);
	my %nametrees = ();
	my $nexus = "";
	my $taxblock = "no";
	my $out = "no";

	my $rooted = "no";

	my $ERR = 0;
	my $ERROR = "";
	my $WARNING = "";

	# =====================================
	# Check existence of files
	# =====================================
	if (($comptrees eq "") || ($reftree eq "")){
		usage();
	}
	elsif ((!-e $reftree) || (!-e $comptrees)){
		if ((!-e $reftree) && (!-e $comptrees)){
			$ERROR = "     Files '$reftree' & '$comptrees' do not exist\n";
		}
		elsif (!-e $reftree){
			$ERROR = "     File '$reftree' does not exist\n";
		}
		else{
			$ERROR = "     File '$comptrees' does not exist\n";
		}
		die("ERROR:\n$ERROR\n");
	}

	# =====================================
	# =====================================

	# =====================================
	#Check line break
	# =====================================
	print "Line break...  ";
	my @lines_ref = break_line($reftree);
	my @lines_comp = break_line($comptrees);
	my $breakline_ref = shift(@lines_ref);
	my $breakline_comp = shift(@lines_comp);
	my $bk_line_line = "OK ($breakline_ref & $breakline_comp)";

	print "$bk_line_line\n";
	my $bk_line = $breakline_comp;

	# =====================================
	# =====================================

	# =====================================
	# Check output files
	# =====================================
	if (($scaled ne "no") || ($table ne "no") || ($partitions ne "no")){
		if ($overwrite eq "no"){
			my $printout = "";
			if (-e $scaled){
				$printout = "WARNING: File '$scaled' exists. Do you want to overwrite it? (y/n) ";
			}
			if (-e $table){
				if ($printout ne ""){
					$printout = "WARNING: Files '$scaled' & '$table' exist. Do you want to overwrite them? (y/n) ";
				}
				else{
					$printout = "WARNING: File '$table' exists. Do you want to overwrite it? (y/n) ";
				}
			}
			if (-e $partitions){
				if ($printout ne ""){
					if (grep(/them/,$printout)){
						$printout = "WARNING: Files '$scaled', '$table' & '$partitions' exist. Do you want to overwrite them? (y/n) ";
					}
					else{
						if (grep (/$table/, $printout)){
							$printout = "WARNING: Files '$table' & '$partitions' exist. Do you want to overwrite them? (y/n) ";
						}
						else{
							$printout = "WARNING: Files '$scaled' & '$partitions' exist. Do you want to overwrite them? (y/n) ";
						}
					}
				}
				else{
					$printout = "WARNING: File '$partitions' exists. Do you want to overwrite it? (y/n) ";
				}
			}
			if ($printout ne ""){
				my $answer = "";
				while ($answer !~ /y|Y|n|N/){
					print $printout;
					$answer = <STDIN>;
				}
				chomp($answer);
				if ($answer !~ m/^y|Y$/){
					die("ERROR:\n     Output files overwrite not allowed.\n\n");
				}
				else{
					unlink ("$scaled");
					unlink ("$table");
					unlink ("$partitions");
				}
			}
		}
		else{
			unlink ("$scaled");
			unlink ("$table");
			unlink ("$partitions");
		}
	}
	# =====================================
	# =====================================

	# =====================================
	# CHECK if files are nexus or newick
	# =====================================
	print "File format...  ";

	my @lines = ([@lines_ref], [@lines_comp]);
	my $fileformat_line = "";
	my $formaterr = "";
	my @reftrees = ();
	my @comptrees = ();
	for (my $l=0; $l < scalar(@lines); $l++){
		my $newick = "yes";
		my @aux = @{$lines[$l]};

	# Check newick
	# -------------------------------------
		my $save = "";
		foreach my $a (@aux){
			next if ($a =~ m/^[\n|\/|\#|\>]/);
 			$a =~ s/^\[.+\]//g if ($a =~ m/^\[.+\]/);
			$save = $save.$a;
			chomp($save);
		}
		my @trs = split (/;/, $save);
		foreach (@trs){
			$_ = $_.";";
			my @aux = split(/\,/, $_);
			my $n_intercom = scalar(@aux);
			@aux = split(/ *\: *[0-9]+\.*[0-9]*e*E*\-*[0-9]* *| *\, *| *\( *| *\) *[0-9]*\.*[0-9]* *| *\;/, $_);
			my $n_spp = 0;
			foreach my $a (@aux){
				$n_spp++ if ($a ne "");
			}
			@aux = split (/\(/,$_);
			my $o_par = scalar(@aux);
			@aux = split (/\)/,$_);
			my $c_par = scalar(@aux);
			if ((!/^\(/ && !/\;$/) || ($n_intercom != $n_spp) || ($o_par != $c_par)){
				$newick = "no";
				last;
			}
		}
		if ($newick eq "yes"){
			if ($fileformat_line eq ""){
				$fileformat_line = "newick";
			}
			else{
				$fileformat_line = $fileformat_line." & newick";
			}
			if ($l == 0){ @reftrees = @trs; }
			else { @comptrees = @trs; }

		}


	# Check nexus
	# -------------------------------------
		if ($newick eq "no"){
			my $nex = "no";
			my $confirmnex = "no";
			foreach (@aux){
				next if (/^[\n|\/|\>]|^[ |\t]+$/);
				if (/^[ |\t]*\#NEXUS/){
					$nex = "yes";
					next;
				}
				if ($nex eq "yes"){
					if (/TREE/i){
						if ($fileformat_line eq ""){
							$fileformat_line = "nexus";
						}
						else{
							$fileformat_line = $fileformat_line." & nexus";
						}
						$confirmnex = "yes";
						last;
					}
				}
			}

			if (($nex eq "no") || ($nex eq "yes" && $confirmnex eq "no")){
				$nexus = $nexus."no";
			}
			else{
				my ($nametr, $txbl, $trs) = nexus2newick(\@aux);
				my %auxname = (%$nametr);
				$taxblock = ($$txbl);
				if ($l == 0){ @reftrees = (@$trs); }
				else{ @comptrees = (@$trs); }
				$fileformat_line = $fileformat_line." with taxa block" if ($taxblock eq "yes");
				if ($l == 0){ $nametrees{$reftree} = \%auxname; }
				else { $nametrees{$comptrees} = \%auxname; }
				$nexus = $nexus."yes";
			}
		}
		else{ $nexus = $nexus."no"; }

		if (($newick eq "no") && ($nexus =~ m/no$/)){
			if ($formaterr ne ""){ $formaterr = $formaterr."\n          &\n          format of file "; }
			else{ $formaterr = "Format of file "; }
			if ($l == 0) { $formaterr = $formaterr." '$reftree' not recognized"; }
			else { $formaterr = $formaterr." '$comptrees' not recognized"; }
		}
	}

	if ($formaterr ne ""){
		print "\n";
		die ("     ERROR:\n          $formaterr.\n\n");
	}
	else{
		print "OK ($fileformat_line)\n";
	}
	# =====================================
	# =====================================

	# =====================================
	# Reading trees
	# =====================================
	my @rtrs = ();
	my @ctrs = ();
	my @auxtrees = (@reftrees, @comptrees);

	for (my $t=0; $t < scalar(@auxtrees); $t++){
		my @trees = ();
		my $line = "";
		chomp($auxtrees[$t]);
		next if ($auxtrees[$t] =~ m/^[\n|\/|\#|\>]/);
		$auxtrees[$t] =~ s /^\[.+\]//g;
		if ($auxtrees[$t] =~ m/\;$/){
			$line = $line.$auxtrees[$t];
			if ($t < scalar (@reftrees)){ push (@rtrs, $line); }
			else { push (@ctrs, $line); }
			$line = "";
		}
		else{
			$line = $line.$auxtrees[$t];
		}
	}
	# =====================================
	# =====================================

	# =====================================
	# Check if there is more than one reference tree
	# =====================================
	print "Number of reference trees...  ";

	if (scalar(@rtrs) > 1){
		$WARNING = "\n     WARNING:  There is more than one tree in the reference file. Only the first one will be used.";
		print "$WARNING\n";
	}
	else{
		print "OK\n";
	}
	my $t0 = shift(@rtrs); #ÊThe first reference tree will be used

	# =====================================
	# =====================================

	# =====================================
	# Check branch lengths
	# =====================================
	print "Checking branch lengths...  ";
	my @trs = ($t0, @ctrs);
	my $nobrlen = "no";
	my $wrongbrlen = "no";	
	my $indexes = "";
	my $indexes_wrbrlen = "";
	my $n_wo_brlen = 0;
	my $n_wr_wo_brlen = 0;	
	my $refbr = "yes";
	my $wrrefbr = "yes";
	for (my $i=0; $i < scalar(@trs); $i++){
		if ($i == 0){ # It's reftree
			if (!grep (/\:/, $trs[$i])){
				$refbr = "no";
				$nobrlen = "yes";
			}
		}
		else{ # It's comptrees
			if (!grep (/\:/, $trs[$i])){
				if ($nexus =~ m /yes$/){
					$indexes = $indexes.$nametrees{$comptrees}{$i-1}.", ";
				}
				else{
					$indexes = $indexes.($i).", ";
				}
				$nobrlen = "yes";
				$n_wo_brlen++;
			}
		}
		
	#ÊChecks branch lengths are right (not negative)
		my @aux = split (/\:|\)|\,/, $trs[$i]);
		foreach my $a (@aux){
			if ($a=~ m/^\-[0-9]+/){
				$wrongbrlen = "yes";
				if ($i == 0){					
					$wrrefbr = "no";
				}
				else{
					if ($nexus =~ m /yes$/){
						$indexes_wrbrlen = $indexes_wrbrlen.$nametrees{$comptrees}{$i-1}.", ";
					}
					else{
						$indexes_wrbrlen = $indexes_wrbrlen.($i).", ";
					}
					$n_wr_wo_brlen++;
				}
			}
		}		
	}
	$indexes =~ s/, $//g;
	$indexes_wrbrlen =~ s/, $//g;

	if ($nobrlen eq "yes"){
		my $brlenerr = "";
		if ($refbr eq "no") {
			if ($indexes eq ""){
				$brlenerr = "          Tree '$reftree' has no branch lengths.\n";
			}
			else{
				if ($n_wo_brlen > 1){
					$brlenerr = "          Tree '$reftree' has no branch lengths\n          &\n          Trees '$indexes' in file '$comptrees' have no branch lengths.\n";
				}
				else{
					$brlenerr = "          Tree '$reftree' has no branch lengths\n          &\n          Tree '$indexes' in file '$comptrees' has no branch lengths.\n";
				}
			}
		}
		else{
			if ($n_wo_brlen > 1){
				$brlenerr = "          Trees '$indexes' in file '$comptrees' have no branch lengths.\n";
			}
			else{
				$brlenerr = "          Tree '$indexes' in file '$comptrees' has no branch lengths.\n";
			}
		}

		$ERROR = $ERROR."\n     ERROR (branch lengths):\n$brlenerr\n";
		print "$ERROR";
		$ERR = 1;
	}
	else{

		$ERROR = "$ERROR\n\n" if ($ERROR ne "");
		my $brlenerr = "";	
		if ($wrongbrlen eq "yes"){
			if ($wrrefbr eq "no"){
				if ($indexes_wrbrlen eq ""){
					$brlenerr = "          Branch lengths of tree '$reftree' are incorrect (negative values).\n";
				}
				else{
					if ($n_wo_brlen > 1){
						$brlenerr = "          Branch lengths of trees '$reftree' are incorrect (negative values)\n          &\n          Branch lengths of trees '$indexes_wrbrlen' in file '$comptrees' are incorrect (negative values).\n";
					}
					else{
						$brlenerr = "          Branch lengths of tree '$reftree' are incorrect (negative values)\n          &\n          Branch lengths of tree '$indexes_wrbrlen' in file '$comptrees' are incorrect (negative values).\n";
					}
				}
			}
			else{
				if ($n_wr_wo_brlen > 1){
					$brlenerr = "          Branch lengths of trees '$indexes_wrbrlen' in file '$comptrees' are incorrect (negative values).\n";
				}
				else{
					$brlenerr = "          Branch lengths of tree '$indexes_wrbrlen' in file '$comptrees' are incorrect (negative values).\n";
				}			
			}
			if ($nobrlen eq "yes"){
				$ERROR = $ERROR."\n\n$brlenerr\n";
			}
			else{
				$ERROR = "\n     ERROR (branch lengths):\n$brlenerr\n";			
			}
			print "$ERROR";
			$ERR = 1;			
		}
		else{
			print "OK\n";
		}
	}	
	$ERROR = "";
	# =====================================
	# =====================================

	# =====================================
	# Check species - tips
	# =====================================
	print "Tips in trees...  ";
	my @species = ();
	foreach my $t (@trs){
		my @sps = ();
		my @auxsps = ();
		if (grep (/\:/, $t)){
			@auxsps = split (/ *\: *\-*[0-9]+\.*[0-9]*e*E*-*[0-9]* *| *\, *| *\( *| *\) *[0-9]*\.*[0-9]* *| *\; */, $t);
		}
		else{
			@auxsps = split (/ *\, *| *\( *| *\) *| *\; */, $t);
		}
		foreach my $s (@auxsps){
			my $auxs = $s;
			$auxs =~ s/\'|\"//g;
			$auxs =~ s/ /_/g;
			$t =~ s/$s/$auxs/g;
			push (@sps, $auxs) if ($auxs ne "");
		}
		push (@species, [@sps]);
		$t =~ s/ //g;
	}
	$t0 = $trs[0]; # Update reference tree
	@ctrs = @trs; # Update comparison trees
	shift(@ctrs);
	my @dif_n_sps = ();
	my @dif_sps = ();
	for (my $i=1; $i < scalar(@species); $i++){
		my $hits = 0;
		if (scalar(@{$species[0]}) != scalar (@{$species[$i]})){
			push (@dif_n_sps, $i);
		}
		else{
			for (my $j=0; $j < scalar(@{$species[0]}); $j++){
				if (grep (/^$species[0][$j]$/, @{$species[$i]})){
					$hits++;
				}
			}
			if ($hits != scalar(@{$species[0]})){
				push (@dif_sps, $i);
			}
		}
	}

	my $spperr = "";
	$indexes = "";
	my $indexes_spp = "";

	if (scalar (@dif_n_sps) > 0){
		foreach my $dns (@dif_n_sps){
			if (exists ($nametrees{$comptrees}{($dns-1)})){
				$indexes = $indexes.$nametrees{$comptrees}{($dns-1)}.", ";
			}
			else{
				$indexes = $indexes.($dns).", ";
			}
			$indexes_spp = $indexes_spp.scalar(@{$species[$dns]}).", ";
		}
		$indexes =~ s/, $//g;
		$indexes_spp =~ s/, $//g;
		my @aux_ind = split (/, /, $indexes_spp);
		$indexes_spp = "";
		for (my $i=0; $i < $#aux_ind; $i++){
			$indexes_spp = $indexes_spp.$aux_ind[$i].", ";
		}
		$indexes_spp =~ s/, $/ /g;
		$indexes_spp = $indexes_spp."and " if ($indexes_spp ne "");
		$indexes_spp = $indexes_spp."$aux_ind[$#aux_ind]";
		my $n_ind = scalar(@aux_ind);
		if ($n_ind > 1){
			$spperr = $spperr."          Reference tree has ".scalar(@{$species[0]})." tips\n          &\n          trees '$indexes' in '$comptrees' have ".$indexes_spp." tips, respectively.";
		}
		else{
			$spperr = $spperr."          Reference tree has ".scalar(@{$species[0]})." tips\n          &\n          tree '$indexes' in '$comptrees' has ".$indexes_spp." tips.";
		}
	}

	$indexes = "";
	if (scalar (@dif_sps) > 0){
		foreach my $ds (@dif_sps){
			if (exists ($nametrees{$comptrees}{($ds-1)})){
				$indexes = $indexes.$nametrees{$comptrees}{($ds-1)}.", ";
			}
			else{
				$indexes = $indexes.($ds).", ";
			}
		}
		$indexes =~ s/, $//g;
		my @aux_ind = split (/, /, $indexes);
		my $n_ind = scalar(@aux_ind);
		if ($n_ind > 1){
			$spperr = $spperr."\n          &\n" if ($spperr ne "");
			$spperr = $spperr."          The tips of trees '".$indexes."' in '$comptrees' are different to those of reference tree.";
		}
		else{
			$spperr = $spperr."\n          &\n" if ($spperr ne "");
			$spperr = $spperr."          The tips of tree '".$indexes."' in '$comptrees' are different to those of reference tree.";
		}
	}
	if ($spperr ne ""){
		$spperr = "     ERROR (species - tips):\n".$spperr;
		$ERROR = "$ERROR\n$spperr\n\n";
		print "$ERROR";
		$ERR = 1;
	}
	else{
		print "OK\n";
	}

	$ERROR = "";

	# =====================================
	# =====================================

	# =====================================
	# Check root
	# =====================================
	print "All trees rooted or all unrooted...  ";
	my $n_rooted = 0;
	my $n_trees = 0;
	my @treenumber = ();
	foreach my $t (@trs){
		my $aux = $t;
		$aux =~ s/\)[0-9]+\:[0-9]+\.*[0-9]*/\)/g;# for bootstrap
		$aux =~ s/^\(|\)\;$|\:[0-9]+\.*[0-9]*//g;
		my @aux2 = split (/(\(|\))/,$aux);
		my $opened = 0;
		my $closed = 0;
		my $chg = "";
		foreach my $a (@aux2){
			next if ($a eq "\n" || $a eq "");
			$opened ++ if ($a eq "\(");

			$closed ++ if ($a eq "\)");

			if ($opened != 0){
				$a = "\\(" if ($a eq "\(");
				$a = "\\)" if ($a eq "\)");
				$chg = $chg.$a;
			}
		if (($opened == $closed) && ($opened != 0) && ($chg ne $a)){
				$aux =~ s/$chg/X/g;
				$chg = "";
				$opened = 0;
				$closed = 0;
			}
		}
		my @cm = split (/\,/, $aux);
		if (scalar (@cm) == 2) {
			push (@treenumber, $n_trees."R");
			$n_rooted++;
		}
		else {
			push (@treenumber, $n_trees."U");
		}
		$n_trees++;
	}

	if (($n_rooted != scalar (@trs)) && ($n_rooted != 0)){
		my $about_trees = "";
		my $indexes = "";

		if ($treenumber[0] =~ /R$/){ # Reference tree is rooted
			$rooted = "yes";
			my $n_unrooted = 0;
			foreach my $tn (@treenumber){
				if (grep (/U/,$tn)){
					$tn =~ s/U//g;
					if ($nexus =~ m/yes$/){
						$indexes = $indexes.$nametrees{$comptrees}{$tn-1}.", ";
					}
					else{
						$indexes = $indexes.$tn.", ";
					}
					$n_unrooted++;
				}
			}
			shift (@treenumber);
			if ($n_unrooted == scalar(@treenumber)){
				if (scalar(@treenumber) == 1){
					$about_trees = "          Reference tree '$reftree' is rooted\n          &\n          '$comptrees' is unrooted";
				}
				else{
					$about_trees = "          Reference tree '$reftree' is rooted\n          &\n          all trees in '$comptrees' are unrooted";
				}
			}
			else{
				$indexes =~ s/\, $//g;
				my @aux = split(/,/,$indexes);
				my $n_ind = scalar(@aux);
				if ($n_ind > 1){
					$about_trees = "          Reference tree '$reftree' is rooted\n          &\n          trees '$indexes' in '$comptrees' are unrooted";
				}
				else{
					$about_trees = "          Reference tree '$reftree' is rooted\n          &\n          tree '$indexes' in '$comptrees' is unrooted";
				}
			}
		}
		else {
			$n_rooted = 0;
			foreach my $tn (@treenumber){
				if (grep (/R/,$tn)){
					$tn =~ s/R//g;
					if ($nexus =~ m/yes$/){
						$indexes = $indexes.$nametrees{$comptrees}{$tn-1}.", ";
					}
					else{
						$indexes = $indexes.$tn.", ";
					}
					$n_rooted++;
				}
			}
			shift (@treenumber);
			if ($n_rooted == scalar(@treenumber)){
				if (scalar(@treenumber) == 1){
					$about_trees = "          Reference tree '$reftree' is unrooted tree \n          &\n          '$comptrees' is rooted";
				}
				else{
					$about_trees = "          Reference tree '$reftree' is unrooted \n          &\n          all trees in '$comptrees' are rooted" ;
				}
			}
			else{
				$indexes =~ s/\, $//g;
				my @aux = split(/,/,$indexes);
				my $n_ind = scalar(@aux);
				if ($n_ind > 1){
					$about_trees = "          Reference tree '$reftree' is unrooted \n          &\n          trees '$indexes' in '$comptrees' are rooted";
				}
				else{
					$about_trees = "          Reference tree '$reftree' is unrooted \n          &\n          tree '$indexes' in '$comptrees' is rooted";
				}
			}
		}
		if ($ERROR ne ""){
			$ERROR = "$ERROR\n\nERROR (rooted - unrooted):\n$about_trees.\n\n     Trees should be all rooted or all unrooted.\n\n";
		}
		else{
			$ERROR = "\n     ERROR (rooted - unrooted):\n$about_trees.\n\n          Trees should be all rooted or all unrooted.\n\n";
		}
		print "$ERROR";
		$ERR = 1;
	}
	else{
		if ($treenumber[0] =~ /R$/){
			$rooted = "yes";
		}
		print "OK\n";
	}
	$ERROR = "";
	# =====================================
	# =====================================


	my $errormsg = "\n".'*********************************************'."\n".'  There were one or more ERRORS. See above.  '."\n".'*********************************************'."\n\n";

	die ("$errormsg") if ($ERR == 1);

	return (\$reftree, \$comptrees, \%nametrees, \$taxblock, \$rooted, \$bk_line, \$t0, \@ctrs);
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# /// Check kind of line break (and change it to UNIX if necessary) ///
sub break_line{
	my $file = shift(@_);
	my $break = "";
	my $br = "\n";
	open (FILE, "$file") or die ("Cannot open file $file");
		my @file = (<FILE>);
	close(FILE);
	foreach (@file){
		if (/\r\n/){ $break = 'DOS'; $br = "\r\n"; last; }
		if (/\r/)  { $break = 'MAC'; $br = "\r"; last; }
		if (/\n/)  { $break = 'UNIX'; last; }
	}
	my @lines = ();
	if ($break ne 'UNIX'){
		my $all_lines = "";
		foreach (@file){
			s/$br/\n/g;
			$all_lines = $all_lines.$_;
		}
		@lines = split (/\n/, $all_lines);
	}
	else{
		@lines = @file;
	}
	return($break, @lines);
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# /// Converts nexus tree to newick ///
sub nexus2newick{
	my ($a, $b) = @_;
	my @aux = (@$a);
	my @newicktrees = ();
	my %nametrees = ();
	my $taxblock = "no";
	my $treeblock = "no";
	if (grep (/translate/i, @aux)){ # Nexus with taxa block
		$taxblock = "yes";
		my %taxa = ();
		my $n = 0;
		my $continue = "no";
		foreach (@aux){

			if (/^Begin +TREES/i){
				$treeblock = "yes";
			}
			next if ($treeblock eq "no");
			last if (/End\;/ && $treeblock eq "yes");

			$continue = "yes" if (/Translate/i);
			next if ($continue eq "no");
			
			if (/^( *|\t*)[0-9]+/){
				s/^( +|\t+)//g;
				my @aux2 = split(/^([0-9]+)/, $_);
				my $number = $aux2[1];
				my $taxon = $aux2[2];
				$taxon =~ s/\t+| +|\,|\;|\n//g;
				$taxa{"$number"} = $taxon;
			}
			
			if (/^ *\t*[U|R]*TREE/i){
				my $tree = $_;
				my $nametree = $tree;
				$nametree =~ s/^ *\t*[U|R]*TREE *\t*|\[\&r*u*\]| *\t*(\[.+\])* *\t*\=.*|\'|\n//ig;
				$tree =~ s/.+=[ |\t|\[|\]|\&|a-z]*\(/\(/ig;
				if (grep (/\:/,$tree)){
					my @auxtree = split (/([0-9]+\:)/, $tree);
					$tree = "";
					foreach my $a (@auxtree){
						if ($a =~ m/^[0-9]+\:/){
							$a =~ s/\://g;
							$a = $taxa{$a};
							$a = $a.":";
						}
						$tree = $tree.$a;
					}
				}
				else{# No branch lengths
					my @auxtree = split (/([0-9]+)/, $tree);
					$tree = "";
					foreach my $a (@auxtree){
						$a = $taxa{$a} if ($a =~ m/^[0-9]+/);
						$tree = $tree.$a;
					}
				}
				push (@newicktrees, $tree);
				$nametrees{$n} = $nametree;
				$n++;
			}
		}
	}
	else{ # Nexus without taxa block
		my $n = 0;
		foreach (@aux){
			if (/^Begin +TREES/i){
				$treeblock = "yes";
			}
			next if ($treeblock eq "no");
			last if (/End\;/ && $treeblock eq "yes");

			if (/^ *\t*[U|R]*TREE/i){
				my $tree = $_;
				my $nametree = $tree;
				$nametree =~ s/^ *\t*[U|R]*TREE *\t*|\[\&r*u*\]| *\t*(\[.+\])* *\t*\=.*|\'|\n//ig;
				$tree =~ s/.+=[ |\t|\[|\]|\&|a-z]*\(/\(/ig;
				push (@newicktrees, $tree);
				$nametrees{$n} = $nametree;
				$n++;
			}
		}
	}

	return (\%nametrees, \$taxblock, \@newicktrees);
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# /// Converts a newick file to nexus ///
sub newick2nexus{
	my ($sc, $txbl, $nt, $cpt) = @_;
	my $scaled = ($$sc);
	my $taxblock = ($$txbl);
	my %nametrees = (%$nt);
	my $comptrees = ($$cpt);
	open (FILE, "$scaled") or die ("Cannot open $scaled");
		my @auxfile = (<FILE>);
	close (FILE);
	open (FILE, ">$scaled") or die ("Cannot write to $scaled");
		print FILE "#NEXUS\n\n";
		print FILE "\[\!\n";
		print FILE "File generated by Ktreedist version $version\n\n";
		print FILE "Reference file: '$reftree'\n";
		print FILE "Comparison file: '$comptrees'\n";
		print FILE "]\n\n";
		print FILE "Begin TREES;\n\n";
		my $ntips = 1;
		my %nametips = ();
		if ($taxblock eq "yes"){
			print FILE "\tTranslate\n";
			my @aux = ();
			foreach my $af (@auxfile){
				if ($af =~ m/^\(/){
					@aux = split (/\,|\(|\)|\;/, $af);
					last;
				}
			}
			my $ind = 1;
			my $taxblock_output = "";
			foreach my $a (@aux){
				if ($a =~ m/^.+\:/){
					$a =~ s/\:.+//g;
					$taxblock_output = $taxblock_output."\t\t$ind\t'$a',\n";
					$nametips{$ind} = $a;
					$ind++;
				}
			}
			$ntips = $ind;
			$taxblock_output =~ s/,\n$/\n\t;\n\n/g;
			print FILE $taxblock_output;
		}
		for (my $i=0; $i < scalar(@auxfile); $i++){
			if ($auxfile[$i] =~ m/^\(/){
				if ($taxblock eq "yes"){
					for (my $j=1; $j < $ntips; $j++){
						$auxfile[$i] =~ s/$nametips{$j}/$j/g;
					}
				}
				print FILE "TREE '$nametrees{$comptrees}{$i}' = $auxfile[$i]";
			}
			else{
				print FILE $auxfile[$i];
			}
		}
		print FILE "\nEnd;\n";
	close (FILE);
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# /// Get partitions ///
sub get_partitions_tree{
	my $tree = shift(@_);
	my $separated = shift(@_);
	$tree =~ s/\) *[0-9]*\.*[0-9]*/\)/g;
	my @aux = split(/(:|,|\)|\;|\()/, $tree);
	my @treesplit = ();
	foreach my $a (@aux){
		push (@treesplit, $a) if ($a ne "");
	}

	my %brlen = ();
	my %aux_brlen = ();
	my @cladosA = ();
	my @cladosB = ();
	my @species = ();
	my @opened = ();
	my $n_comas = 0;
	my $cont = -1;
	my $rooted = "no";

	# Get internal partitions
	for (my $i=0; $i < scalar(@treesplit); $i++) {
		$n_comas ++ if (($treesplit[$i] eq "\,") && ($#opened == 0));
		if ($treesplit[$i] eq "\("){
			$cont++;
			push (@opened, $cont);
		}
		if ($treesplit[$i] eq "\:"){
			if ($treesplit[$i-1] =~ m/\w+/){
				if ($treesplit[$i-2] ne "\)"){# it's a tip and not a bootstrap value
					push (@species, $treesplit[$i-1]);
					$aux_brlen{$species[$#species]} = $treesplit[$i+1];
					foreach my $o (@opened){
						push (@{$cladosA[$o]}, $treesplit[$i-1]);
					}
				}
			}
		}
		if ($treesplit[$i] eq "\)"){
			if ($treesplit[$i+1] ne ";"){
				$brlen{$cladosA[$opened[$#opened]]} = $treesplit[$i+2];
			}
			pop(@opened);
		}
	}
	$rooted = "yes" if ($n_comas == 1);

	my $index = $#cladosA;
	shift(@cladosA); # First partition is all the tree

	#ÊGet symmetric ones
	for (my $i=0; $i < scalar(@cladosA); $i++){
		foreach my $sp (@species){
			if (!grep(/^$sp$/, @{$cladosA[$i]})){
				push (@{$cladosB[$i]}, $sp);
			}
		}
		push (@{$cladosB[$i]}, "root") if ($rooted eq "yes");
		$brlen{$cladosB[$i]} = $brlen{$cladosA[$i]};
	}
	# Get tip partitions
	foreach my $sp (@species){
		$cladosA[$index][0] = $sp;
		$brlen{$cladosA[$#cladosA]} = $aux_brlen{$sp};
		foreach my $sp2 (@species){
			push (@{$cladosB[$index]}, $sp2) if ($sp2 ne $sp);
		}
		push (@{$cladosB[$index]}, "root") if ($rooted eq "yes");
		$brlen{$cladosB[$#cladosB]} = $aux_brlen{$sp};
		$index++;
	}

	my @clados = ();
	my %brlen_aux = ();

	# Join partitions when requested
	if ($separated eq "no"){
		for (my $i=0; $i < scalar (@cladosA); $i++){
			push (@clados, [@{$cladosA[$i]}]);
			$brlen_aux{$clados[$#clados]} = $brlen{$cladosA[$i]};
			push (@clados, [@{$cladosB[$i]}]);
			$brlen_aux{$clados[$#clados]} = $brlen{$cladosB[$i]};
		}
		%brlen = %brlen_aux;
		return (\@clados, \%brlen, \@species, \$rooted);
	}
	else{
		return (\@cladosA, \@cladosB, \%brlen, \@species, \$rooted);
	}
}

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#Ê//// Get matrix ////

sub get_matrix{
	my ($cl0A, $cl0B, $cl1A, $cl1B, $brl0, $brl1) = @_;
	my @clados0A = (@$cl0A);
	my @clados0B = (@$cl0B);
	my @clados1A = (@$cl1A);
	my @clados1B = (@$cl1B);
	my @clados1 = (@clados1A, @clados1B);
	my %brlen0 = (%$brl0);
	my %brlen1 = (%$brl1);

	my @matrix = ();
	my $m = 0;
	# Add partitions there are in tree 0
	for (my $i=0; $i < scalar(@clados0A); $i++){
		my $priorm = $m;
		for (my $j=0; $j < scalar(@clados1); $j++){
			my $n_equal = 0;
			if (scalar(@{$clados0A[$i]}) == scalar(@{$clados1[$j]})){
				foreach my $cl0 (@{$clados0A[$i]}){
					if (grep (/^$cl0$/, @{$clados1[$j]})){
						$n_equal++;
					}
				}

				if ($n_equal == scalar (@{$clados0A[$i]})){
					if (($hard_polytomies eq "off") || ($hard_polytomies eq "on" && $brlen0{$clados0A[$i]} != 0 && $brlen1{$clados1[$j]} != 0)){
						push (@{$matrix[$m]}, [@{$clados0A[$i]}]);
						push (@{$matrix[$m]}, [@{$clados0B[$i]}]);
						push (@{$matrix[$m]}, $brlen0{$clados0A[$i]});
						push (@{$matrix[$m]}, $brlen1{$clados1[$j]});
						push (@{$matrix[$m]}, $brlen0{$clados0A[$i]}); # "NONE" column
						push (@{$matrix[$m]}, $brlen1{$clados1[$j]}); # "NONE" column
						$m++;
					}
					last;
				}
			}
		}
		if ($priorm == $m){
			if ($hard_polytomies eq "off" || ($hard_polytomies eq "on" && $brlen0{$clados0A[$i]} != 0)){
				push (@{$matrix[$m]}, [@{$clados0A[$i]}]); 
				push (@{$matrix[$m]}, [@{$clados0B[$i]}]);
				push (@{$matrix[$m]}, $brlen0{$clados0A[$i]});
				push (@{$matrix[$m]}, 0);
				push (@{$matrix[$m]}, $brlen0{$clados0A[$i]}); # "NONE" column
				push (@{$matrix[$m]}, "NONE"); # "NONE" column
				$m++;
			}
		}
	}

	# Add partitions there are in tree 1 but not in tree 0
	for (my $i=0; $i < scalar(@clados1A); $i++){
		my $n_nonequal = 0;
		for (my $j=0; $j < scalar(@matrix); $j++){
			my $n_equal = 0;
			for (my $k=0; $k < 2; $k++){
				$n_equal = 0;
				if (scalar(@{$clados1A[$i]}) == scalar (@{$matrix[$j][$k]})){
					foreach my $cl1 (@{$clados1A[$i]}){
						if (grep (/^$cl1$/, @{$matrix[$j][$k]})){
							$n_equal++;
						}
					}
					last if ($n_equal == scalar(@{$clados1A[$i]}));
				}
			}
			if ($n_equal != scalar(@{$clados1A[$i]})){
				$n_nonequal ++;
			}
			else{
				last;
			}
		}
		if ($n_nonequal == scalar (@matrix)){
			if (($hard_polytomies eq "off") || ($hard_polytomies eq "on" && $brlen1{$clados1A[$i]} != 0)){
				push (@{$matrix[$m]}, [@{$clados1A[$i]}]);
				push (@{$matrix[$m]}, [@{$clados1B[$i]}]);
				push (@{$matrix[$m]}, 0);
				push (@{$matrix[$m]}, $brlen1{$clados1A[$i]});
 				push (@{$matrix[$m]}, "NONE"); # "NONE" column
				push (@{$matrix[$m]}, $brlen1{$clados1A[$i]}); # "NONE" column
				$m++;
			}
		}
	}
	return (@matrix);
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#Ê//// Calculate K-value ////

sub get_K_value{
	my @matrix = @_;
	my $kval = 0;

	# Numerator calculation
	my $l = 0;
	my $numerator = 0;
	for (my $i=0; $i < scalar(@matrix); $i++){
		$numerator = $numerator + $matrix[$i][2]*$matrix[$i][3];
	}

	#Denominator calculation
	my $denominator = 0;

	for (my $i=0; $i < scalar (@matrix); $i++){
		$denominator = $denominator+($matrix[$i][3]**2);
	}

	#K-value calculation
	$kval = $numerator/$denominator;

	return ($kval);
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#Ê//// Scale the tree by K-value ////

sub modify_brlen_tree{
	my ($fi, $tr, $brl, $Ka) = @_;
	my $file = ($$fi);
	my $tree = ($$tr);
	my %brlen = (%$brl);
	my $K = ($$Ka);

	my $scaled_tree = "";
	for (values %brlen){
		$_ = sprintf ("%.".$precision."f",($_*$K));
	}

	my @aux = split(/(:|,|\)|\;|\()/, $tree);
	my @treesplit = ();
	foreach my $a (@aux){
		push (@treesplit, $a) if ($a ne "");
	}
	$tree = "";
	for (my $i=0; $i < scalar(@treesplit); $i++) {
		$treesplit[$i+1] = sprintf ("%.".$precision."f",($treesplit[$i+1]*$K)) if ($treesplit[$i] eq "\:");
		$tree = $tree.$treesplit[$i];
	}

	$scaled_tree = $tree if ($file ne "no");  # output scaled tree

	return (\$tree, \%brlen, \$scaled_tree);
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#Ê//// Calculate K-score ////
sub get_K_score{
	my ($ma, $kv) = @_;
	my @matrix = (@$ma);
	my $K_value = ($$kv);
	my $BLD = 0;
	my $drf = 0;

	for (my $i=0; $i < scalar (@matrix); $i++){
		$matrix[$i][6] = $matrix[$i][3]*$K_value;
		$BLD = $BLD+($matrix[$i][2] - $matrix[$i][6])**2;
		if ($matrix[$i][4] eq "NONE" || $matrix[$i][5] eq "NONE"){
			$drf++;
			$matrix[$i][7] = "NONE" if ($matrix[$i][5] eq "NONE");
		}
	}
	$BLD = sqrt($BLD);
	return (\$BLD, \$drf, \@matrix);
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# Change line break from UNIX to DOS or MAC
sub unix2macdos{
	my $kind = shift(@_);
	my $file = shift(@_);
	my $bl = "\n";
	if ($kind eq "MAC"){ $bl = "\r"; }
	elsif ($kind eq "DOS"){	$bl = "\r\n"; }
	open (FILE, "$file") or die ("Cannot open $file");
		my @aux=(<FILE>);
	close(FILE);
	open (FILE, ">$file") or die ("Cannot write to $file");
	foreach (@aux){
		s/\n/$bl/g;
		print FILE $_;
	}
	close(FILE);
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#Ê/// Generate lines for display output ///
sub output_display{
	my ($ct, $nt, $namt, $ks, $kv, $r, $dr, $nc1) = @_;
	my $comptrees = ($$ct);
	my $ntree = ($$nt);
	my %nametrees = (%$namt);
	my $K_score = ($$ks);
	my $K_value = ($$kv);		
	my $robinson = ($$r);		
	my $drf = ($$dr);		
	my $n_clados1 = ($$nc1);
	
	# Comparison tree name	
	my $name = "";
	my $output = "";
	if (exists($nametrees{$comptrees}{$ntree})){ $name = $nametrees{$comptrees}{$ntree}; }
	else{ $name = "Tree ".($ntree+1); }
	$name = substr($name, 0, 15) if (length($name) > 15);
	my $pre_name_bl = "";
	my $post_name_bl = "";
	for (my $i=0; $i < int((15-length($name))/2); $i++){ $post_name_bl = $post_name_bl." "; }
	my $check = ((15-length($name))%2);
	if ($check != 0){
		for (my $i=0; $i < int(((15-length($name))/2)+1); $i++){ $pre_name_bl = $pre_name_bl." "; }
	}
	else{
		$pre_name_bl = $post_name_bl;
	}
	$output = $output."$name$pre_name_bl$post_name_bl";

	# K-score & K-value columns
	$output = $output.sprintf("    %.5f      %.5f ",$K_score, $K_value);

	# Symmetric difference score column
	if ($robinson eq "yes"){
		my $pre_rf_bl = "";
		my $post_rf_bl = "";
		for (my $i=0; $i < int((15-length($drf))/2); $i++){ $post_rf_bl = $post_rf_bl." "; }

		my $check = ((15-length($drf))%2);
		if ($check != 0){
			for (my $i=0; $i < int(((15-length($drf))/2)+1); $i++){ $pre_rf_bl = $pre_rf_bl." "; }
		}
		else{
			$pre_rf_bl = $post_rf_bl;
		}

		$output = $output.sprintf("  $pre_rf_bl%u$post_rf_bl", $drf);
	}

	# Number of partitions column
	if ($npartitions eq "yes"){
		my $parts = "$n_clados1";
		my $pre_parts_bl = "";
		my $post_parts_bl = "";
		for (my $i=0; $i < int((12-length($parts))/2); $i++){ $post_parts_bl = $post_parts_bl." "; }
		$check = ((12-length($parts))%2);
		if ($check != 0){
			for (my $i=0; $i < int(((12-length($parts))/2)+1); $i++){ $pre_parts_bl = $pre_parts_bl." "; }
		}
		else{
			$pre_parts_bl = $post_parts_bl;
		}
		$output = $output."  $pre_parts_bl$parts$post_parts_bl";
	}
	$output = $output."  \n";		
	
	return ($output);
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#Ê/// Generate lines for table file output ///
sub get_table {
	my ($rt, $ct, $namt, $fn, $nt, $np, $ks, $kv, $r, $dr, $nc0, $nc1) = @_;
	my $reftree = ($$rt);
	my $comptrees = ($$ct);
	my %nametrees = (%$namt);
	my $ntree = ($$nt);		
	my $filename_comptrees_woext = ($$fn);	
	my $npartitions = ($$np);			
	my $K_score = ($$ks);
	my $K_value = ($$kv);
	my $robinson = ($$r);
	my $drf = ($$dr);
	my $n_clados0 = ($$nc0);
	my $n_clados1 = ($$nc1);
	
	my $linea = "";
	if (exists ($nametrees{$reftree}{0})){
		$linea = "$nametrees{$reftree}{0} (ref)";
# 				$linea = $filename_reftree_woext."."."$nametrees{$reftree}{0} (ref)"; # For nexus, in case you want the filename to appear in the table
	}
	else{
		$linea = "$filename_reftree_woext (ref)";
	}
	if (exists ($nametrees{$comptrees}{$ntree})){
		$linea = $linea." - ".$nametrees{$comptrees}{$ntree};
#				$linea = $linea." - ".$filename_comptrees_woext.".".$nametrees{$comptrees}{$ntree}; # For nexus, in case you want the filename to appear in the table
	}
	else{
		$linea = $linea." - ".$filename_comptrees_woext." ("."Tree ".($ntree+1).")";
	}
	$linea = $linea.sprintf("\t%.5f\t%.5f", $K_score, $K_value);
	$linea = $linea.sprintf("\t%u",$drf) if ($robinson eq "yes");
	$linea = $linea."\t$n_clados0\t$n_clados1" if ($npartitions eq "yes");

	return ($linea);
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#Ê/// Generate lines for partitions file output ///
sub get_parts_table {
	my ($ct, $namt, $nt, $ma) = @_;

	my $comptrees = ($$ct);
	my %nametrees = (%$namt);
	my $ntree = ($$nt);		
	my @matrix = (@$ma);				
	
	my $whatree = "";
	if (exists ($nametrees{$comptrees}{$ntree})){
		$whatree = $nametrees{$comptrees}{$ntree};
	}
	else{
		$whatree = " Tree ".($ntree+1)." ";
	}
	my $chars = length($whatree);
	my $upanddown = "";
	for (my $i=0; $i < ($chars+2); $i++){ $upanddown = $upanddown."-"; }
	my $linea = "\n$upanddown\n\|$whatree\|\n$upanddown\n";
	$linea = $linea."Brl_ref_tree\tBrl_comp_tree\tBrl_comp_tree_K\tPartition\n";
	for (my $i=0; $i < scalar(@matrix); $i++){
		if ($matrix[$i][4] eq "NONE"){
			$linea = $linea."$matrix[$i][4]\t".sprintf ("%.5f\t%.5f\t", $matrix[$i][5],$matrix[$i][6]);
		}
		elsif ($matrix[$i][5] eq "NONE"){
			$linea = $linea.sprintf ("%.5f\t", $matrix[$i][4])."$matrix[$i][5]\t$matrix[$i][7]\t";
		}
		else{
			$linea = $linea.sprintf ("%.5f\t%.5f\t%.5f\t", $matrix[$i][2],$matrix[$i][3],$matrix[$i][6]);
		}

		$linea = $linea."(";
		my @part = ($matrix[$i][0]);
		for (my $h=0; $h < scalar (@part); $h++){
			for (my $q=0; $q < scalar (@{$part[$h]}); $q++){
 					$linea = $linea. "$part[$h][$q], ";
			}
		}
		$linea =~ s/, $//g;
		$linea = $linea.") / (";
		@part = ($matrix[$i][1]);
		for (my $h=0; $h < scalar (@part); $h++){
			for (my $q=0; $q < scalar (@{$part[$h]}); $q++){
 					$linea = $linea. "$part[$h][$q], ";
			}
		}
		$linea =~ s/, $/)\n/g;

	}

	return ($linea);		
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
	
