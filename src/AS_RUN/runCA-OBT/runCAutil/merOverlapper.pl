use strict;

sub findOvermerryFailures ($$) {
    my $outDir   = shift @_;
    my $ovmJobs  = shift @_;
    my $failures = 0;

    for (my $i=1; $i<=$ovmJobs; $i++) {
        my $out = substr("0000" . $i, -4);
        if (! -e "$wrk/$outDir/seeds/$out.success") {
            $failures++;
        }
    }

    return $failures;
}


sub findOlapFromSeedsFailures ($$) {
    my $outDir   = shift @_;
    my $olpJobs  = shift @_;
    my $failures = 0;

    for (my $i=1; $i<=$olpJobs; $i++) {
        my $out = substr("0000" . $i, -4);
        if (! -e "$wrk/$outDir/olaps/$out.success") {
            $failures++;
        }
    }

    return $failures;
}


sub merOverlapper($) {
    my $isTrim = shift @_;

    return if (-d "$wrk/$asm.ovlStore");
    return if (-d "$wrk/$asm.obtStore") && ($isTrim eq "trim");

    caFailure("merOverlapper()-- Help!  I have no frags!\n") if ($numFrags == 0);
    caFailure("merOverlapper()-- I need to know if I'm trimming or assembling!\n") if (!defined($isTrim));

    my $outDir  = "1-overlapper";
    my $ovlOpt  = "";
    my $merSize = getGlobal("ovlMerSize");
    my $merComp = getGlobal("merCompression");
    my $merType = "ovl";

    if ($isTrim eq "trim") {
        $outDir  = "0-overlaptrim-overlap";
        $ovlOpt  = "-G";
        $merSize = getGlobal("obtMerSize");
        $merComp = getGlobal("merCompression");
        $merType = "obt";
    }

    #  To prevent infinite loops -- stop now if the overmerry script
    #  exists.  This will unfortunately make restarting from transient
    #  failures non-trivial.
    #
    caFailure("mer overlapper failed.") if (-e "$wrk/$outDir/overmerry.sh");

    system("mkdir $wrk/$outDir")       if (! -d "$wrk/$outDir");
    system("mkdir $wrk/$outDir/seeds") if (! -d "$wrk/$outDir/seeds");
    system("mkdir $wrk/$outDir/olaps") if (! -d "$wrk/$outDir/olaps");

    #  Make the 2-frgcorr directory (to hold the corrections
    #  output) and claim that fragment correction is all done.
    #  after this, the rest of the fragment/overlap correction
    #  pipeline Just Works.
    #
    if ($isTrim ne "trim") {
        system("mkdir $wrk/2-frgcorr") if (! -d "$wrk/2-frgcorr");
        touch("$wrk/2-frgcorr/jobsCreated.success");
    }


    my $ovmBatchSize = getGlobal("merOverlapperSeedBatchSize");
    my $ovmJobs      = int($numFrags / ($ovmBatchSize-1)) + 1;

    my $olpBatchSize = getGlobal("merOverlapperExtendBatchSize");
    my $olpJobs      = int($numFrags / ($olpBatchSize-1)) + 1;


    #  Create overmerry and olap-from-seeds jobs
    #
    if (! -e "$wrk/$outDir/overmerry.sh") {
        open(F, "> $wrk/$outDir/overmerry.sh") or caFailure("Can't open '$wrk/$outDir/overmerry.sh'\n");
        print F "#!/bin/sh\n";
        print F "\n";
        print F "jobid=\$SGE_TASK_ID\n";
        print F "if [ x\$jobid = x -o x\$jobid = xundefined ]; then\n";
        print F "  jobid=\$1\n";
        print F "fi\n";
        print F "if [ x\$jobid = x ]; then\n";
        print F "  echo Error: I need SGE_TASK_ID set, or a job index on the command line.\n";
        print F "  exit 1\n";
        print F "fi\n";
        print F "\n";
        print F "jobid=`printf %04d \$jobid`\n";
        print F "minid=`expr \$jobid \\* $ovmBatchSize - $ovmBatchSize + 1`\n";
        print F "maxid=`expr \$jobid \\* $ovmBatchSize`\n";
        print F "runid=\$\$\n";
        print F "\n";
        print F "if [ \$maxid -gt $numFrags ] ; then\n";
        print F "  maxid=$numFrags\n";
        print F "fi\n";
        print F "if [ \$minid -gt \$maxid ] ; then\n";
        print F "  echo Job partitioning error -- minid=\$minid maxid=\$maxid.\n";
        print F "  exit\n";
        print F "fi\n";
        print F "\n";
        print F "AS_OVL_ERROR_RATE=", getGlobal("ovlErrorRate"), "\n";
        print F "AS_CNS_ERROR_RATE=", getGlobal("cnsErrorRate"), "\n";
        print F "AS_CGW_ERROR_RATE=", getGlobal("cgwErrorRate"), "\n";
        print F "export AS_OVL_ERROR_RATE AS_CNS_ERROR_RATE AS_CGW_ERROR_RATE\n";
        print F "\n";
        print F "if [ ! -d $wrk/$outDir/seeds ]; then\n";
        print F "  mkdir $wrk/$outDir/seeds\n";
        print F "fi\n";
        print F "\n";
        print F "if [ -e $wrk/$outDir/seeds/\$jobid.success ]; then\n";
        print F "  echo Job previously completed successfully.\n";
        print F "  exit\n";
        print F "fi\n";

        print F getBinDirectoryShellCode();

        print F "\$bin/overmerry \\\n";
        print F " -g  $wrk/$asm.gkpStore \\\n";
        if ($ovmJobs > 1) {
            print F " -mc $wrk/0-mercounts/$asm-C-ms$merSize-cm$merComp \\\n";
            print F " -tb \$minid -te \$maxid \\\n";
        }
        print F " -m $merSize \\\n";
        print F " -c $merComp \\\n";
        print F " -t " . getGlobal("merOverlapperThreads") . "\\\n";
        print F " -o $wrk/$outDir/seeds/\$jobid.ovm \\\n";
        print F " > $wrk/$outDir/seeds/\$jobid.ovm.err 2>&1 \\\n";
        print F "&& \\\n";
        print F "touch $wrk/$outDir/seeds/\$jobid.success\n";
        close(F);

        system("chmod +x $wrk/$outDir/overmerry.sh");
    }

    if (! -e "$wrk/$outDir/olap-from-seeds.sh") {
        open(F, "> $wrk/$outDir/olap-from-seeds.sh") or caFailure("Can't open '$wrk/$outDir/olap-from-seeds.sh'\n");
        print F "#!/bin/sh\n";
        print F "\n";
        print F "jobid=\$SGE_TASK_ID\n";
        print F "if [ x\$jobid = x -o x\$jobid = xundefined ]; then\n";
        print F "  jobid=\$1\n";
        print F "fi\n";
        print F "if [ x\$jobid = x ]; then\n";
        print F "  echo Error: I need SGE_TASK_ID set, or a job index on the command line.\n";
        print F "  exit 1\n";
        print F "fi\n";
        print F "\n";
        print F "jobid=`printf %04d \$jobid`\n";
        print F "minid=`expr \$jobid \\* $olpBatchSize - $olpBatchSize + 1`\n";
        print F "maxid=`expr \$jobid \\* $olpBatchSize`\n";
        print F "runid=\$\$\n";
        print F "\n";
        print F "if [ \$maxid -gt $numFrags ] ; then\n";
        print F "  maxid=$numFrags\n";
        print F "fi\n";
        print F "if [ \$minid -gt \$maxid ] ; then\n";
        print F "  echo Job partitioning error -- minid=\$minid maxid=\$maxid.\n";
        print F "  exit\n";
        print F "fi\n";
        print F "\n";
        print F "AS_OVL_ERROR_RATE=", getGlobal("ovlErrorRate"), "\n";
        print F "AS_CNS_ERROR_RATE=", getGlobal("cnsErrorRate"), "\n";
        print F "AS_CGW_ERROR_RATE=", getGlobal("cgwErrorRate"), "\n";
        print F "export AS_OVL_ERROR_RATE AS_CNS_ERROR_RATE AS_CGW_ERROR_RATE\n";
        print F "\n";
        print F "if [ ! -d $wrk/$outDir/olaps ]; then\n";
        print F "  mkdir $wrk/$outDir/olaps\n";
        print F "fi\n";
        print F "\n";
        print F "if [ -e $wrk/$outDir/olaps/\$jobid.success ]; then\n";
        print F "  echo Job previously completed successfully.\n";
        print F "  exit\n";
        print F "fi\n";

        print F getBinDirectoryShellCode();

        print F "\$bin/olap-from-seeds \\\n";
        print F " -a -b \\\n";
        print F " -t 1 \\\n";
        print F " -S $wrk/$outDir/$asm.merStore \\\n";

        if ($isTrim eq "trim") {
            print F " -G \\\n";  #  Trim only
            print F " -o $wrk/$outDir/olaps/\$jobid.ovr \\\n";
            print F " $wrk/$asm.gkpStore \\\n";
            print F " \$minid \$maxid \\\n";
            print F " > $wrk/$outDir/olaps/$asm.\$jobid.ovb.err 2>&1 \\\n";
            print F " &&  \\\n";
            print F "\$bin/acceptableOBToverlap \\\n";
            print F " < $wrk/$outDir/olaps/\$jobid.ovr \\\n";
            print F " > $wrk/$outDir/olaps/\$jobid.ovb \\\n";
        } else {
            print F " -c $wrk/2-frgcorr/\$jobid.frgcorr \\\n";
            print F " -o $wrk/$outDir/olaps/\$jobid.ovb \\\n";
            print F " $wrk/$asm.gkpStore \\\n";
            print F " \$minid \$maxid \\\n";
            print F " > $wrk/$outDir/olaps/$asm.\$jobid.ovb.err 2>&1 \\\n";
            print F "&& \\\n";
            print F "touch $wrk/2-frgcorr/\$jobid.success \\\n";
        }

        print F "&& \\\n";
        print F "touch $wrk/$outDir/olaps/\$jobid.success\n";

        close(F);

        system("chmod +x $wrk/$outDir/olap-from-seeds.sh");
    }



    #  If we're partitioned, we need to run meryl.  We don't care
    #  about single-copy mers, only mers that occur in two or more
    #  frags.
    #
    #  Possibly should only compute for the $isTrim setting.
    #
    if ($ovmJobs > 1) {
        setGlobal("obtMerThreshold", 2);
        setGlobal("ovlMerThreshold", 2);
        meryl();
    }


    #  Submit to the grid (or tell the user to do it), or just run
    #  things here
    #
    if (findOvermerryFailures($outDir, $ovmJobs) > 0) {
        if (getGlobal("useGrid") && getGlobal("ovlOnGrid")) {
            my $sge        = getGlobal("sge");
            my $sgeOverlap = getGlobal("sgeMerOverlapSeed");

            my $SGE;
            $SGE  = "qsub $sge $sgeOverlap -N NAME \\\n";
            $SGE .= "  -t MINMAX \\\n";
            $SGE .= "  -j y -o $wrk/$outDir/seeds/\\\$TASK_ID.out \\\n";
            $SGE .= "  $wrk/$outDir/overmerry.sh\n";

            my $waitTag = submitBatchJobs("mer", $SGE, $ovmJobs, getGlobal("merOverlapperThreads"));
            submitScript($waitTag) if (runningOnGrid());
            exit(0);
        } else {
            for (my $i=1; $i<=$ovmJobs; $i++) {
                my $out = substr("0000" . $i, -4);
                #runCommand("$wrk/$outDir", "$wrk/$outDir/overmerry.sh $i > $wrk/$outDir/seeds/$out.out 2>&1");
                &scheduler::schedulerSubmit("sh $wrk/$outDir/overmerry.sh $i > $wrk/$outDir/seeds/$out.out 2>&1");
            }

            &scheduler::schedulerSetNumberOfProcesses(getGlobal("merOverlapperSeedConcurrency"));
            &scheduler::schedulerFinish();
        }
    }

    #  Make sure everything finished ok.
    #
    {
        my $f = findOvermerryFailures($outDir, $ovmJobs);
        caFailure("There were $f overmerry failures.\n") if ($f > 0);
    }

    if (! -e "$wrk/$outDir/$asm.merStore") {
        my $bin = getBinDirectory();
        my $cmd;
        $cmd  = "$bin/overlapStore";
        $cmd .= " -c $wrk/$outDir/$asm.merStore";
        $cmd .= " -g $wrk/$asm.gkpStore";
        $cmd .= " -M " . getGlobal("ovlStoreMemory");
        $cmd .= " $wrk/$outDir/seeds/*.ovm";
        $cmd .= " > $wrk/$outDir/$asm.merStore.err 2>&1";
        if (runCommand("$wrk/$outDir", $cmd)) {
            rename "$wrk/$outDir/$asm.merStore", "$wrk/$outDir/$asm.merStore.FAILED";
            caFailure("Failed.\n");
        }
    }



    #  Submit to the grid (or tell the user to do it), or just run
    #  things here
    #
    if (findOlapFromSeedsFailures($outDir, $olpJobs) > 0) {
        if (getGlobal("useGrid") && getGlobal("ovlOnGrid")) {
            my $sge        = getGlobal("sge");
            my $sgeOverlap = getGlobal("sgeMerOverlapExtend");

            my $SGE;
            $SGE  = "qsub $sge $sgeOverlap -N NAME \\\n";
            $SGE .= "  -t MINMAX \\\n";
            $SGE .= "  -j y -o $wrk/$outDir/olaps/\\\$TASK_ID.out \\\n";
            $SGE .= "  $wrk/$outDir/olap-from-seeds.sh\n";

            my $waitTag = submitBatchJobs("olp", $SGE, $olpJobs, 1);
            submitScript("$waitTag") if (runningOnGrid());
            exit(0);
        } else {
            for (my $i=1; $i<=$olpJobs; $i++) {
                my $out = substr("0000" . $i, -4);
                #runCommand("$wrk/$outDir", "$wrk/$outDir/olap-from-seeds.sh $i > $wrk/$outDir/olaps/$out.out 2>&1");
                &scheduler::schedulerSubmit("sh $wrk/$outDir/olap-from-seeds.sh $i > $wrk/$outDir/olaps/$out.out 2>&1");
            }

            &scheduler::schedulerSetNumberOfProcesses(getGlobal("merOverlapperExtendConcurrency"));
            &scheduler::schedulerFinish();
        }
    }

    #  Make sure everything finished ok.
    #
    {
        my $f = findOlapFromSeedsFailures($outDir, $olpJobs);
        caFailure("There were $f olap-from-seeds failures.\n") if ($f > 0);
    }
}
