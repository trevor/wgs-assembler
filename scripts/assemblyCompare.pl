#!/usr/local/bin/perl
# $Id: assemblyCompare.pl,v 1.1.2.1 2005-10-07 20:45:41 catmandew Exp $
#
# Program to compare two assemblies
#
#  Compare two qc files
#  Compare TAMPA results
#  Run nucmer & show-coords & show-snps
#  Compare show-coords & show-snps output
#
#   Written by Ian Dew
#

use Carp;
use strict;
use FileHandle;
use Getopt::Long;

my $MY_VERSION = " Version 1.01 (Build " . (qw/$Revision: 1.1.2.1 $/ )[1]. ")";
my $MY_APPLICATION = "assemblyCompare";

my $HELPTEXT = qq~
Compare two caqc-generated qc files

    assemblyCompare  [options]  <-r dir1>  <-q dir2>

    <-r dir1>   The 'reference' assembly directory

    <-q dir2>   The 'query' assembly directory
  
    options:
      -h               Print help.
  
      -v <level>       Set verbosity to level.

      -c               Do not compare qc files

      -t               Do not compare TAMPA results

      -m               Do not run mummer and compare results

$MY_VERSION

~;

my $QC = "qc";
my $TAMPA1 = "tampa1";
my $TAMPA2 = "tampa2";
my $MUMMER = "mummer";


######################################################################
# Parse the command line
######################################################################
my $helpRequested;
my $verboseLevel = 0;
my $dir1;
my $dir2;
my $dontQC;
my $dontTAMPA;
my $dontMummer;

GetOptions("r=s" => \$dir1,
           "q=s" => \$dir2,
           "c" => \$dontQC,
           "t" => \$dontTAMPA,
           "m" => \$dontMummer,
           "h|help" => \$helpRequested,
           "v|V|verbose:1" => \$verboseLevel
           ) or die $HELPTEXT;

if($helpRequested)
{
  print STDERR "Help requested:\n\n";
  print STDERR $HELPTEXT;
  exit 0;
}

if(!$dir1 || ! -d $dir1)
{
  print STDERR "Please specify a reference assembly directory\n\n";
  print STDERR $HELPTEXT;
  exit 1;
}

if(!$dir2 || ! -d $dir2)
{
  print STDERR "Please specify a query assembly directory\n\n";
  print STDERR $HELPTEXT;
  exit 1;
}

######################################################################
# Verify existence of required files
######################################################################
# get full lists
opendir(THISDIR, $dir2);
my @files2 = grep (-f, readdir(THISDIR));
closedir(THISDIR);

# need unambiguous qc file in each
# need unambiguous pair of tampa files in each
# need unambiguous pair of .scaffolds.fasta files in each
my %dir1Files = IdentifyRequiredFiles($dir1);
my %dir2Files = IdentifyRequiredFiles($dir2);


sub IdentifyRequiredFiles($)
{
  my $dir = shift;

  opendir(THISDIR, $dir);
  my @files = grep (-f, readdir(THISDIR));
  closedir(THISDIR);

  my %reqFiles;
  for(my $i = 0; $i <= $#files; $i++)
  {
    if(!$dontQC && $files[$i] =~ m/^(.*)\.qc$/)
    {
      if(defined($reqFiles{$QC}))
      {
        printf STDERR "Ambiguous - Multiple qc files in $dir\n";
        exit;
      }
      $reqFiles{$QC} = $files[$i];
      next;
    }
    if(!$dontTAMPA)
    {
      if($files[$i] =~ m/^(.*)\.intra\.summary\.tampa$/)
      {
        if(defined($reqFiles{$TAMPA1}))
        {
          printf STDERR "Ambiguous - Multiple intra-sequence TAMPA files in $dir\n";
          exit;
        }
        $reqFiles{$TAMPA1} = $files[$i];
        next;
      }
      if($files[$i] =~ m/^(.*)\.inter\.summary\.tampa$/)
  }
  return %reqFiles;
}
