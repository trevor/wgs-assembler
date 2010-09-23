use strict;

#  Don't do interleaved merging unless we are throwing stones.

sub CGW ($$$$$$) {
    my $thisDir     = shift @_;
    my $lastDir     = shift @_;
    my $cgiFile     = shift @_;
    my $stoneLevel  = shift @_;
    my $logickp     = shift @_;
    my $finalRun    = shift @_;

    return($thisDir) if (-e "$wrk/$thisDir/cgw.success");

    my $lastckp = findLastCheckpoint($lastDir)  if (defined($lastDir));
    my $ckp     = "-R $lastckp -N $logickp"  if (defined($lastckp) && defined($logickp));

    #  If there is a timing file here, assume we are restarting.  Not
    #  all restarts are possible, but we try hard to make it so.
    #
    if (-e "$wrk/$thisDir/$asm.timing") {
        my $restartckp = undef;

        open(F, "< $wrk/$thisDir/$asm.timing");
        while (<F>) {
            print STDERR $_;
            if (m/Writing.*ckp.(\d+)\s\(logical\s(.+)\)/) {
                $restartckp = "-R $1 -N $2";
            }
        }
        close(F);

        if (!defined($restartckp)) {
            print STDERR "Found an empty timing file, starting from the beginning: $ckp\n";
        } else {
            $ckp = $restartckp;
            print STDERR "Found a timing file, restarting: $ckp\n";
        }
    }

    system("mkdir $wrk/$thisDir")               if (! -d "$wrk/$thisDir");
    system("mkdir $wrk/$asm.SeqStore")          if (! -d "$wrk/$asm.SeqStore");

    if (!defined($cgiFile)) {
        open(F, "ls $wrk/5-consensus |");
        while (<F>) {
            chomp;
            if (m/cgi$/) {
                $cgiFile .= " $wrk/5-consensus/$_";
            }
        }
        close(F);
    } else {
        system("ln -s $cgiFile $wrk/$thisDir/$asm.cgi") if (! -e "$wrk/$thisDir/$asm.cgi");
        $cgiFile = "$wrk/$thisDir/$asm.cgi";
    }

    system("ln -s ../$asm.SeqStore  $wrk/$thisDir/$asm.SeqStore")     if (! -e "$wrk/$thisDir/$asm.SeqStore");

    system("ln -s ../$lastDir/$asm.ckp.$lastckp $wrk/$thisDir/$asm.ckp.$lastckp") if (defined($lastDir));

    if (-e "$wrk/$thisDir/cgw.out") {
        my $ckp = findLastCheckpoint($thisDir);
        my $ver = "00";
        while (-e "$wrk/$thisDir/cgw.out.$ver.ckp.$ckp") {
            $ver++;
        }
        rename "$wrk/$thisDir/cgw.out", "$wrk/$thisDir/cgw.out.$ver.ckp.$ckp"
    }

    my $sampleSize = getGlobal("cgwDistanceSampleSize");

    my $bin = getBinDirectory();
    my $cmd;
    my $astatLow = getGlobal("astatLowBound");
    my $astatHigh = getGlobal("astatHighBound");
    $cmd  = "$bin/cgw $ckp -j $astatLow -k $astatHigh -r 5 -s $stoneLevel ";
    $cmd .= " -S 0 "                               if (($finalRun == 0)   || (getGlobal("doResolveSurrogates") == 0));
    $cmd .= " -G "                                 if (($finalRun == 0)   && (getGlobal("cgwOutputIntermediate") == 0));
    $cmd .= " -z "                                 if (getGlobal("cgwDemoteRBP") == 1);
    $cmd .= " -p " . getGlobal("closurePlacement") if (defined(getGlobal("closurePlacement")));
    $cmd .= " -F "                                 if (getGlobal("toggleDoNotDemote") != 0);
    $cmd .= " -u $wrk/4-unitigger/$asm.unused.ovl" if (getGlobal("cgwUseUnitigOverlaps") != 0);
    $cmd .= " -m $sampleSize";
    $cmd .= " -g $wrk/$asm.gkpStore ";
    $cmd .= " -o $wrk/$thisDir/$asm ";
    $cmd .= " $cgiFile ";
    $cmd .= " > $wrk/$thisDir/cgw.out 2>&1";
    if (runCommand("$wrk/$thisDir", $cmd)) {
        caFailure("scaffolder failed", "$wrk/$thisDir/cgw.out");
    }


    open(F, "ls -1 $wrk/$thisDir |");
    while (<F>) {
        chomp;

        if (m/\.log$/) {
            system("mkdir $wrk/$thisDir/log")        if (! -d "$wrk/$thisDir/log");
            rename "$wrk/$thisDir/$_", "$wrk/$thisDir/log/$_";
        }

        if (m/\.analysis$/) {
            system("mkdir $wrk/$thisDir/analysis")   if (! -d "$wrk/$thisDir/analysis");
            rename "$wrk/$thisDir/$_", "$wrk/$thisDir/analysis/$_";
        }
    }
    close(F);


    if (getGlobal("cgwPurgeCheckpoints") != 0) {
        my $f = findFirstCheckpoint($thisDir);
        my $l = findLastCheckpoint($thisDir);

        while ($f < $l) {
            #print STDERR "Purging $wrk/$thisDir/$asm.ckp.$f\n";
            unlink "$wrk/$thisDir/$asm.ckp.$f";
            $f++;
        }
    }

    touch("$wrk/$thisDir/cgw.success");

    return $thisDir;
}


sub eCR ($$$) {
    my $thisDir = shift @_;
    my $lastDir = shift @_;
    my $iter    = shift @_;

    return $thisDir if (-e "$wrk/$thisDir/extendClearRanges.success");

    my $lastckp = findLastCheckpoint($lastDir);

    system("mkdir $wrk/$thisDir") if (! -d "$wrk/$thisDir");

    system("ln -s ../$lastDir/$asm.ckp.$lastckp $wrk/$thisDir/$asm.ckp.$lastckp")  if (! -e "$wrk/$thisDir/$asm.ckp.$lastckp");
    system("ln -s ../$asm.SeqStore              $wrk/$thisDir/$asm.SeqStore")      if (! -e "$wrk/$thisDir/$asm.SeqStore");

    #  Partition eCR.

    if (! -e "$wrk/$thisDir/extendClearRanges.partitionInfo") {
        my $bin = getBinDirectory();
        my $cmd;
        $cmd  = "$bin/extendClearRangesPartition ";
        $cmd .= " -g $wrk/$asm.gkpStore ";
        $cmd .= " -n $lastckp ";
        $cmd .= " -c $asm ";
        $cmd .= " -N 4 ";
        $cmd .= " -p $wrk/$thisDir/extendClearRanges.partitionInfo";
        $cmd .= "  > $wrk/$thisDir/extendClearRanges.partitionInfo.err 2>&1";

        if (runCommand("$wrk/$thisDir", $cmd)) {
            caFailure("extendClearRanges partitioning failed", "$wrk/$thisDir/extendClearRanges.partitionInfo.err");
        }
    }

    #  Read the partitioning info, create jobs.  No partitions?  No ECR jobs.

    my @jobs;

    open(P, "< $wrk/$thisDir/extendClearRanges.partitionInfo") or caFailure("failed to find extendClearRanges partitioning file $wrk/$thisDir/extendClearRanges.partitionInfo", undef);
    while (<P>) {
        #  Fields are: partitionNum BgnScf partitionNum EndScf NumFrags

        my @v = split '\s+', $_;

        my $curScaffold = substr("000000000$v[1]", -7);
        my $endScaffold = substr("000000000$v[3]", -7);

        my $j = "$wrk/$thisDir/extendClearRanges-scaffold.$curScaffold";

        if (! -e "$j.success") {
            my $bin = getBinDirectory();

            $lastckp = findLastCheckpoint($thisDir);

            open(F, "> $j.sh");
            print F "#!" . getGlobal("shell") . "\n\n";
            print F "\n";
            print F "AS_OVL_ERROR_RATE=", getGlobal("ovlErrorRate"), "\n";
            print F "AS_CNS_ERROR_RATE=", getGlobal("cnsErrorRate"), "\n";
            print F "AS_CGW_ERROR_RATE=", getGlobal("cgwErrorRate"), "\n";
            print F "export AS_OVL_ERROR_RATE AS_CNS_ERROR_RATE AS_CGW_ERROR_RATE\n";
            print F "\n";
            print F "$bin/extendClearRanges \\\n";
            print F " -g $wrk/$asm.gkpStore \\\n";
            print F " -n $lastckp \\\n";
            print F " -c $asm \\\n";
            print F " -b $curScaffold -e $endScaffold \\\n";
            print F " -i $iter \\\n";
            print F " > $j.err 2>&1\n";
            close(F);

            system("chmod +x $j.sh");

            push @jobs, "$j";
        }
    }
    close(F);

    #  Run jobs.

    foreach my $j (@jobs) {
        if (runCommand("$wrk/$thisDir", "$j.sh")) {
            caFailure("extendClearRanges failed", "$j.err");
        }
        touch("$j.success");
    }

    touch("$wrk/$thisDir/extendClearRanges.success");

    return $thisDir;
}


sub updateDistanceRecords ($) {
    my $thisDir = shift @_;

    return if (-e "$wrk/$thisDir/cgw.distupdate.success");

    #  Older versions needed to actually compute the updated
    #  distances.  Now, cgw outputs it!  Yay!

    my $bin = getBinDirectory();
    my $cmd;
    my $gkpErrorFile = "$wrk/$thisDir/gkp.distupdate.err";
    $cmd  = "$bin/gatekeeper ";
    $cmd .= " -a -o $wrk/$asm.gkpStore ";
    $cmd .= " -E $gkpErrorFile";
    $cmd .= " $wrk/$thisDir/stat/scaffold_final.distupdate.dst ";
    $cmd .= " $wrk/$thisDir/stat/contig_final.distupdate.dst ";
    $cmd .= " > $wrk/$thisDir/cgw.distupdate.err 2>&1";
    if (runCommand("$wrk/$thisDir", $cmd)) {
        caFailure("gatekeeper distance update failed", "$wrk/$thisDir/cgw.distupdate.err");
    }

    touch("$wrk/$thisDir/cgw.distupdate.success");
}


sub scaffolder ($) {
    my $cgiFile    = shift @_;
    my $lastDir    = undef;
    my $thisDir    = 0;
    my $stoneLevel = getGlobal("stoneLevel");

    goto alldone if (-e "$wrk/7-CGW/cgw.success");

    #  Do an initial CGW to update distances, then update the
    #  gatekeeper.  This initial run shouldn't be used for later
    #  CGW'ing.
    #
    if ((getGlobal("computeInsertSize") == 1) ||
        (getGlobal("computeInsertSize") == 0) && ($numFrags < 1000000)) {
        updateDistanceRecords(CGW("6-clonesize", undef, $cgiFile, $stoneLevel, undef, 0));
    }


    #  If we're not doing eCR, we just do a single scaffolder run, and
    #  get the heck outta here!  OK, we'll do resolveSurrogates(), maybe.
    #
    if (getGlobal("doExtendClearRanges") == 0) {
        $lastDir = CGW("7-$thisDir-CGW", $lastDir, $cgiFile, $stoneLevel, undef, 1);
        $thisDir++;
    } else {

        #  Do the initial CGW, making sure to not throw stones.
        #
        $lastDir = CGW("7-$thisDir-CGW", $lastDir, $cgiFile, 0, undef, 0);
        $thisDir++;

        #  Followed by at least one eCR
        #
        $lastDir = eCR("7-$thisDir-ECR", $lastDir, 1);
        $thisDir++;

        #  Iterate eCR: do another scaffolder still without stones,
        #  then another eCR.  Again, and again, until we get dizzy and
        #  fall over.
        #
        my $iterationMax = getGlobal("doExtendClearRanges") + 1;
        for (my $iteration = 2; $iteration < $iterationMax; $iteration++) {
            $lastDir = CGW("7-$thisDir-CGW", $lastDir, $cgiFile, 0, "ckp01-ABS", 0);
            $thisDir++;

            $lastDir = eCR("7-$thisDir-ECR", $lastDir, $iteration);
            $thisDir++;
        }

        #  Then another scaffolder, chucking stones into the big holes,
        #  filling in surrogates, and writing output.
        #
        $lastDir = CGW("7-$thisDir-CGW", $lastDir, $cgiFile, $stoneLevel, "ckp01-ABS", 1);
        $thisDir++;
    }


    #  And, finally, hold on, we're All Done!  Point to the correct output directory.
    #
    system("ln -s $lastDir $wrk/7-CGW") if (! -d "$wrk/7-CGW");

  alldone:
    stopAfter("scaffolder");
}


1;