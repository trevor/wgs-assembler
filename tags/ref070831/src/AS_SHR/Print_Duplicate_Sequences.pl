#!/usr/local/bin/perl

##########################################################################
#
# This file is part of Celera Assembler, a software program that
# assembles whole-genome shotgun reads into contigs and scaffolds.
# Copyright (C) 1999-2004, The Venter Institute. All rights reserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received (LICENSE.txt) a copy of the GNU General Public
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
##########################################################################

# $Id: Print_Duplicate_Sequences.pl,v 1.1 2006-01-10 22:44:42 kli1000 Exp $

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_f);

getopts("f:");
my $usage = "usage: 
$0 
	-f <fasta file>

	This program will:
		Will read in a FASTA file, and print out id's that have the same sequence.
";

if(!(
	defined($opt_f))){
	die $usage;
}

###############################################################################
# Make sure files open before wasting any time doing anything

open(FASTA_FH, "<$opt_f") || die "Could not open $opt_f\n";

###############################################################################

my %seq_hash;

print STDERR "Processing FASTA file...\n";
my ($defline, $prev_defline, $sequence);
while(<FASTA_FH>){
	chomp;
	
	if(/^>/){
		$defline=$_;
		if($sequence ne ""){
			process_record($prev_defline, $sequence);
			$sequence="";
		}
		$prev_defline=$defline;
	}else{
		$sequence.=$_;
	}
}
process_record($prev_defline, $sequence);

close(FASTA_FH);

my $dup_counter=0;
foreach my $seq(keys %seq_hash){

	if(($#{$seq_hash{$seq}}+1)>=2){
		print $#{$seq_hash{$seq}}+1;
		print "\n";

	}


	if(($#{$seq_hash{$seq}}+1)>=10){
if(0){
		foreach my $id(@{$seq_hash{$seq}}){
			my ($rank, $x, $y)=split /_/, $id;
			print "$dup_counter,$x,$y\n";
		}
		$dup_counter++;
}
if(0){
		print (join ",", @{$seq_hash{$seq}});
		print "\n";
}

	}
}



print STDERR "Completed.\n";

###############################################################################


sub process_record{
	my $defline = shift;
	my $sequence = shift;

	my $id;
	if($defline=~/^>(\S+)/){
		$id=$1;
	}

	push @{$seq_hash{$sequence}}, $id;

}

#------------------------------------------------------------------------------
