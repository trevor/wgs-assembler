#!/usr/bin/perl

use strict;

#  Convert BLASR default output to OBT overlaps:  aIID bIID [f|r] aBgn aEnd bBgn bEnd error

my %IDmap;

open(F, "< $ARGV[0]") or die "Failed to open '$ARGV[0]' for reading sequence names.\n";
while (<F>) {
    my ($uid, $iid, $name) = split '\s+', $_;

    $IDmap{$name} = $iid;
}
close(F);


while (<STDIN>) {
    my @v = split '\s+', $_;

    #  The first read seems to have a sub-read range appended.

    if ($v[0] =~ m/^(.*)\/\d+_\d+$/) {
        $v[0] = $1;
    }

    my $aiid  = $IDmap{$v[0]};
    my $biid  = $IDmap{$v[1]};
    my $fr    = ($v[8] == 0) ? "f" : "r";
    my $error = 100.0 - $v[3];

    ($v[9], $v[10]) = ($v[10], $v[9])  if ($v[8] == 1);

    die "No A iid found for '$v[0]'\n" if (!defined($aiid));
    die "No B iid found for '$v[1]'\n" if (!defined($biid));

    next if ($aiid == $biid);

    #  HACK!
    #$error = 0.01;

    print "$aiid\t$biid\t$fr\t$v[5]\t$v[6]\t$v[9]\t$v[10]\t$error\n";
}
