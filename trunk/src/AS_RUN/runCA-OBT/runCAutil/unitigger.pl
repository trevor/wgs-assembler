use strict;

sub unitigger () {
    my $bin = getBinDirectory();


    if (0) {
        my $cmd = "$bin/removeMateOverlap -gkp $wrk/$asm.gkpStore -ovl $wrk/$asm.ovlStore";
        if (runCommand("$wrk", $cmd)) {
            caFailure("failed to remove mate overlaps", undef);
        }
    }


    if (! -e "$wrk/4-unitigger/unitigger.success") {

        #  Default to unitigger, unless the gkpStore says otherwise.
        #
        if (!defined(getGlobal("unitigger"))) {
            setGlobal("unitigger", "utg");

            open(F, "$bin/gatekeeper -dumplibraries $wrk/$asm.gkpStore |");
            while (<F>) {
                if (m/forceBOGunitigger=1/) {
                    setGlobal("unitigger", "bog");
                }
            }
            close(F);
        }

        system("mkdir $wrk/4-unitigger") if (! -e "$wrk/4-unitigger");

        my $l = getGlobal("utgGenomeSize");
        my $e = getGlobal("utgErrorRate");

        my $B = int($numFrags / getGlobal("cnsPartitions"));
        $B = getGlobal("cnsMinFrags") if ($B < getGlobal("cnsMinFrags"));

        my $unitigger = getGlobal("unitigger");

        my $cmd;

        if ($unitigger eq "bog") {
            my $bmd = getGlobal("bogBadMateDepth");

            $cmd  = "$bin/buildUnitigs ";
            $cmd .= " -O $wrk/$asm.ovlStore ";
            $cmd .= " -G $wrk/$asm.gkpStore ";
            $cmd .= " -T $wrk/$asm.tigStore ";
            $cmd .= " -B $B ";
            $cmd .= " -e $e ";
            $cmd .= " -s $l "   if (defined($l));
            $cmd .= " -b "      if (getGlobal("bogBreakAtIntersections") == 1);
            $cmd .= " -m $bmd " if (defined($bmd));
            $cmd .= " -o $wrk/4-unitigger/$asm ";
            $cmd .= " > $wrk/4-unitigger/unitigger.err 2>&1";
        } elsif ($unitigger eq "utg") {
            my $u = getGlobal("utgBubblePopping");

            $cmd  = "$bin/unitigger ";
            $cmd .= " -I $wrk/$asm.ovlStore ";
            $cmd .= " -F $wrk/$asm.gkpStore ";
            $cmd .= " -T $wrk/$asm.tigStore ";
            $cmd .= " -B $B ";
            $cmd .= " -e $e ";
            $cmd .= " -k " if (getGlobal("utgRecalibrateGAR") == 1);
            $cmd .= " -l $l " if defined($l);
            $cmd .= " -d 1 -x 1 -z 10 -j 5 -U $u ";
            $cmd .= " -o $wrk/4-unitigger/$asm ";
            $cmd .= " > $wrk/4-unitigger/unitigger.err 2>&1";
        } else {
            caFailure("unknown unitigger $unitigger; must be 'bog' or 'utg'", undef);
        }

        stopBefore("unitigger", $cmd);

        if (runCommand("$wrk/4-unitigger", $cmd)) {
            caFailure("failed to unitig", "$wrk/4-unitigger/unitigger.err");
        }

        touch("$wrk/4-unitigger/unitigger.success");
    }

    stopAfter("unitigger");
}

1;
