#!/usr/bin/perl

#  Read in version 1 format fragment file, write version 2.
#
#  Assumes that ALL READS in the input are in the same library.  If
#  not, well, there is then no way to tell which library the unmated
#  reads are really from.  That library is assumed to be in a DST
#  record that appears BEFORE all reads.
#
#  If you give it "-v vecFile" it'll also populate the vector clear
#  range.


use strict;

my $vec;
my %clv;
my %clq;  #  Currently, we never have this info
my $lib;

my $noOBT = 0;
my $is454 = 0;
my $nft   = 0;

my $err = 0;
while (scalar(@ARGV) > 0) {
    my $arg = shift @ARGV;
    if      ($arg eq "-v") {
        $vec = shift @ARGV;
    } elsif ($arg eq "-noobt") {
        $noOBT = 1;
    } elsif ($arg eq "-is454") {
        $is454 = 1;
    } else {
        $err++;
    }
}
if ($err) {
    die "usage: $0 [-v vector-clear-file] [-noobt] [-is454] < old.frg > new.frg\n";
}

if (defined($vec)) {
    open(F, "< $vec") or die "Failed to open '$vec'\n";
    while (<F>) {
        s/^\s+//;
        s/\s$//;
        my @v = split '\s+', $_;
        $clv{$v[0]} = "$v[1],$v[2]";
    }
    close(F);
    print STDERR "Read vector info for ", scalar(keys %clv), " reads.\n";
}


sub readMultiLineDot {
    #$_ = <STDIN>;                #  Read and discard the tag
    my $save = $/;  $/ = ".\n";  #  Prepare to read the whole thing
    my $src = <STDIN>;           #  Read it
    $/ = $save;                  #  Reset our end of line marker
    $src =~ s/\s+//g;            #  Replace spaces and newlines
    chop $src;                   #  Get rid of the trailing .
    return($src);
}



print "{VER\n";
print "ver:2\n";
print "}\n";


while (!eof(STDIN)) {
    my $line = <STDIN>; chomp $line;

    if      ($line =~ m/^{BAT$/) {
        #  Discard BAT
        while ($line ne "}") {
            $line = <STDIN>; chomp $line;
        }
    } elsif ($line =~ m/^{ADT$/) {
        #  Discard ADT/ADL
        $line = <STDIN>; chomp $line;
        while ($line =~ m/^{ADL$/) {
            while ($line ne "}") {
                $line = <STDIN>; chomp $line;
            }
        }
        $line = <STDIN>; chomp $line;
        $line = <STDIN>; chomp $line;
    } elsif ($line =~ m/^{DST$/) {
        #  Convert a DST into a LIB.
        my $acc;
        my $mea;
        my $std;
        while ($line ne "}") {
            if ($line =~ m/^acc:(\d+)$/) {
                $acc = $1;
            }
            if ($line =~ m/^mea:(\d+\.\d+)$/) {
                $mea = $1;
            }
            if ($line =~ m/^std:(\d+\.\d+)$/) {
                $std = $1;
            }
            $line = <STDIN>; chomp $line;
        }
        print "{LIB\n";
        print "act:A\n";
        print "acc:$acc\n";
        print "ori:I\n";
        print "mea:$mea\n";
        print "std:$std\n";
        print "src:\n";
        print "convert-v1-to-v2\n";
        print ".\n";
        print "nft:2\n";
        print "fea:\n";
        print "is454=$is454\n";
        print "doNotOverlapTrim=$noOBT\n";
        print ".\n";
        print "}\n";
        $lib = $acc;
    } elsif ($line =~ m/^{FRG$/) {
        my $acc;
        my $src;
        my $seq;
        my $qlt;
        my $clr;
        while ($line ne "}") {

            if ($line =~ m/^acc:(\d+)$/) {
                $acc = $1;
            }
            if ($line =~ m/^clr:(\d+,\d+)$/) {
                $clr = $1;
            }

            if ($line =~ m/^src:$/) {
                $src = readMultiLineDot();
            }
            if ($line =~ m/^seq:$/) {
                $seq = readMultiLineDot();
            }
            if ($line =~ m/^qlt:$/) {
                $qlt = readMultiLineDot();
            }

            $line = <STDIN>; chomp $line;
        }
        print "{FRG\n";
        print "act:A\n";
        print "acc:$acc\n";
        print "rnd:1\n";
        print "sta:G\n";
        print "lib:$lib\n";
        print "pla:0\n";
        print "loc:0\n";
        print "src:\n";
        print ".\n";
        print "seq:\n";
        print "$seq\n";
        print ".\n";
        print "qlt:\n";
        print "$qlt\n";
        print ".\n";
        print "hps:\n";
        print ".\n";
        if (defined($clv{$acc})) {
            print "clv:$clv{$acc}\n";
        }
        if (defined($clq{$acc})) {
            print "clq:$clq{$acc}\n";
        }
        print "clr:$clr\n";
        print "}\n";
    } elsif ($line =~ m/^{LKG$/) {
        my $fg1;
        my $fg2;
        while ($line ne "}") {
            if ($line =~ m/^fg1:(\d+)$/) {
                $fg1 = $1;
            }
            if ($line =~ m/^fg2:(\d+)$/) {
                $fg2 = $1;
            }
            $line = <STDIN>; chomp $line;
        }
        print "{LKG\n";
        print "act:A\n";
        print "frg:$fg1\n";
        print "frg:$fg2\n";
        print "}\n";
    } else {
        die "Unsupported line '$_'\n";
    }
}