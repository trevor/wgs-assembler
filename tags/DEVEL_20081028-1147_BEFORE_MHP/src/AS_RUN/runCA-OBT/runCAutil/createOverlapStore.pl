use strict;

sub createOverlapStore {

    goto alldone if (-d "$wrk/$asm.ovlStore");

    if (runCommand($wrk, "find $wrk/1-overlapper -name \\*ovb.gz -print > $wrk/$asm.ovlStore.list")) {
        caFailure("Failed to generate a list of all the overlap files.\n");
    }

    my $bin = getBinDirectory();
    my $cmd;
    $cmd  = "$bin/overlapStore ";
    $cmd .= " -c $wrk/$asm.ovlStore.BUILDING ";
    $cmd .= " -g $wrk/$asm.gkpStore ";
    $cmd .= " -M " . getGlobal("ovlStoreMemory");
    $cmd .= " -L $wrk/$asm.ovlStore.list ";
    $cmd .= " > $wrk/$asm.ovlStore.err 2>&1";

    if (runCommand($wrk, $cmd)) {
        caFailure("Failed to create the overlap store.\n");
    }

    rename "$wrk/$asm.ovlStore.BUILDING", "$wrk/$asm.ovlStore";

    rmrf("$asm.ovlStore.list");
    rmrf("$asm.ovlStore.err");

  alldone:
    stopAfter("overlapper");
}

1;