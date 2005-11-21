use strict;

sub unitigger {

    system("mkdir $wrk/4-unitigger") if (! -e "$wrk/4-unitigger");

    if (! -e "$wrk/4-unitigger/$asm.ofgList") {
        if (runCommand("ls -1 $wrk/0-preoverlap/*.ofg > $wrk/4-unitigger/$asm.ofgList")) {
            print STDERR "Failed to find the ofg's.\n";
            rename "$wrk/4-unitigger/$asm.ofgList", "$wrk/4-unitigger/$asm.ofgList.FAILED";
            exit(1);
        }
    }

    if (! -e "$wrk/4-unitigger/unitigger.success") {
        print STDERR "Starting b -- unitigger\n";

        my $cmd;
        $cmd  = "cd $wrk/4-unitigger && ";
        $cmd .= "$bin/unitigger ";
        $cmd .= " -B 250000 ";

        my $l = getGlobal("genomeSize");
        my $m = getGlobal("unitiggerEdges");
        my $n = getGlobal("unitiggerFragments");

        $cmd .= " $l " if defined($l);
        $cmd .= " $m " if defined($m);
        $cmd .= " $n " if defined($n);

        $cmd .= " -c -P -A 1 -d 1 -x 1 -z 10 -j 5 -U 1 -e 15 ";
        $cmd .= " -F $wrk/$asm.frgStore ";
        $cmd .= " -f ";
        $cmd .= " -o $wrk/4-unitigger/$asm.fgbStore ";
        $cmd .= " -L $wrk/4-unitigger/$asm.ofgList ";
        $cmd .= " -I $wrk/$asm.ovlStore ";
        $cmd .= " > $wrk/4-unitigger/unitigger.out ";
        $cmd .= " 2> $wrk/4-unitigger/unitigger.err ";

        if (runCommand($cmd)) {
            print STDERR "Failed to unitig.\n";
            exit(1);
        }

        touch("$wrk/4-unitigger/unitigger.success");
    }
}

1;
