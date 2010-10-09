use strict;

############################################################
#
#  Create overlap jobs
#
#  A better way to do this (maybe) is to make a job array that
#  contains the logic to decide which partition to work on.
#
sub createOverlapJobs {
    my $isTrim = shift @_;

    if (!defined($isTrim)) {
        die "createOverlapJobs()-- I need to know if I'm trimming or assembling!\n";
    }

    my $outDir = "1-overlapper";
    my $merDir = "0-preoverlap";
    my $ovlOpt = "";

    if ($isTrim eq "trim") {
        $outDir = "0-overlaptrim-overlap";
        $merDir = "0-overlaptrim-overlap";
        $ovlOpt = "-G";
    }

    system("mkdir $wrk/$outDir") if (! -d "$wrk/$outDir");

    return if (-e "$wrk/$outDir/jobsCreated.success");

    #  We make a giant job array for this -- we need to know hashBeg,
    #  hashEnd, refBeg and refEnd -- from that we compute batchName
    #  and jobName.
    #
    #  ovlopts.pl returns the batch name ($batchName), the job name
    #  ($jobName) and options to pass to overlap (-h $hashBeg-$hashEnd
    #  -r $refBeg-$refEnd).  From those, we can construct the command
    #  to run.
    #
    open(F, "> $wrk/$outDir/overlap.sh") or die "Can't open '$wrk/$outDir/overlap.sh'\n";
    print F "#!/bin/sh\n";
    print F "\n";
    print F "perl=perl\n";
    print F "if [ -e /usr/bin/perl ]; then\n";
    print F "  perl=/usr/bin/perl\n";
    print F "fi\n";
    print F "if [ -e /usr/local/bin/perl ]; then\n";
    print F "  perl=/usr/local/bin/perl\n";
    print F "fi\n";
    print F "\n";
    print F "bat=`\$perl $wrk/$outDir/ovlopts.pl bat \$SGE_TASK_ID`\n";
    print F "job=`\$perl $wrk/$outDir/ovlopts.pl job \$SGE_TASK_ID`\n";
    print F "opt=`\$perl $wrk/$outDir/ovlopts.pl opt \$SGE_TASK_ID`\n";
    print F "\n";
    print F "if [ ! -d $wrk/$outDir/\$bat ]; then\n";
    print F "  mkdir $wrk/$outDir/\$bat\n";
    print F "fi\n";
    print F "\n";
    print F "echo bat = \$bat\n";
    print F "echo job = \$job\n";
    print F "echo opt = \$opt\n";
    print F "\n";
    print F "echo out = $scratch/\$bat-\$job.ovl\n";
    print F "echo out = $wrk/$outDir/\$bat/\$job.ovl";
    print F "\n";
    print F "$processStats \\\n";
    print F "$gin/overlap -P $ovlOpt -M $ovlMemory -t $ovlThreads \\\n";
    print F "  \$opt \\\n";
    print F "  -k $wrk/0-preoverlap/$asm.nmers.fasta \\\n";
    print F "  -o $scratch/$asm.\$bat-\$job.$$.ovl \\\n";
    print F "  $wrk/$asm.frgStore \\\n";
    print F "&& \\\n";

    if ($isTrim eq "trim") {
        print F "$gin/overlap-convert -b $scratch/$asm.\$bat-\$job.$$.ovl $scratch/$asm.\$bat-\$job.$$.ovb \\\n";
        print F "&& \\\n";
        print F "mv $scratch/$asm.\$bat-\$job.$$.ovb \\\n";
        print F "   $wrk/$outDir/\$bat/\$job.ovb \\\n";
        #print F "&& \\\n";
        #print F "mv $scratch/$asm.\$bat-\$job.$$.ovl \\\n";
        #print F "   $wrk/$outDir/\$bat/\$job.ovl \\\n";
        print F "&& \\\n";
        print F "rm -f $scratch/$asm.\$bat-\$job.$$.ovl \\\n";
    } else {
        print F "mv $scratch/$asm.\$bat-\$job.$$.ovl \\\n";
        print F "   $wrk/$outDir/\$bat/\$job.ovl \\\n";
    }

    print F "&& \\\n";
    print F "touch $wrk/$outDir/\$bat/\$job.success\n";
    print F "exit 0\n";
    close(F);

    system("chmod +x $wrk/$outDir/overlap.sh");

    #  We segment the hash into $numFrags / $ovlHashBlockSize pieces,
    #  and the stream into $numFrags / $ovlRefBlockSize pieces.  Put
    #  all runs for the same hash into a subdirectory.

    my ($hashBeg, $hashEnd, $refBeg, $refEnd) = (1, 0, 1, 0);

    #  Saved for output to ovlopts.pl
    my @bat;
    my @job;
    my @opt;

    #  Number of jobs per batch directory
    #
    my $batchMax  = 200;
    my $batchSize = 0;
    my $batch     = 1;

    my $batchName = substr("0000000000" . $batch, -8);

    while ($hashBeg < $numFrags) {
        $hashEnd = $hashBeg + $ovlHashBlockSize - 1;
        $hashEnd = $numFrags if ($hashEnd > $numFrags);
        $refBeg = 0;
        $refEnd = 0;

        while ($refBeg < $hashEnd) {
            $refEnd = $refBeg + $ovlRefBlockSize - 1;
            $refEnd = $numFrags if ($refEnd > $numFrags);

            #print STDERR "hash: $hashBeg-$hashEnd ref: $refBeg-$refEnd\n";

            my $jobName;
            $jobName .= "h" . substr("0000000000" . $hashBeg, -8);
            $jobName .= "r" . substr("0000000000" . $refBeg, -8);

            push @bat, "$batchName";
            push @job, "$jobName";
            push @opt, "-h $hashBeg-$hashEnd  -r $refBeg-$refEnd";

            $refBeg = $refEnd + 1;

            $batchSize++;
            if ($batchSize >= $batchMax) {
                $batch++;
                $batchName = substr("0000000000" . $batch, -8);
                $batchSize = 0;
            }
        }

        $hashBeg = $hashEnd + 1;
    }

    open(SUB, "> $wrk/$outDir/ovlopts.pl") or die "Failed to open '$wrk/$outDir/ovlopts.pl'\n";
    print SUB "#!/usr/bin/perl\n";
    print SUB "use strict;\n";
    print SUB "my \@bat = (\n";  foreach my $b (@bat) { print SUB "\"$b\",\n"; }  print SUB ");\n";
    print SUB "my \@job = (\n";  foreach my $b (@job) { print SUB "\"$b\",\n"; }  print SUB ");\n";
    print SUB "my \@opt = (\n";  foreach my $b (@opt) { print SUB "\"$b\",\n"; }  print SUB ");\n";
    print SUB "my \$idx = int(\$ARGV[1]) - 1;\n";
    print SUB "if      (\$ARGV[0] eq \"bat\") {\n";
    print SUB "    print \"\$bat[\$idx]\";\n";
    print SUB "} elsif (\$ARGV[0] eq \"job\") {\n";
    print SUB "    print \"\$job[\$idx]\";\n";
    print SUB "} elsif (\$ARGV[0] eq \"opt\") {\n";
    print SUB "    print \"\$opt[\$idx]\";\n";
    print SUB "} else {\n";
    print SUB "    print STDOUT \"Got '\$ARGV[0]' and don't know what to do!\\n\";\n";
    print SUB "    print STDERR \"Got '\$ARGV[0]' and don't know what to do!\\n\";\n";
    print SUB "    die;\n";
    print SUB "}\n";
    print SUB "exit(0);\n";
    close(SUB);

    open(SUB, "> $wrk/$outDir/ovljobs.dat") or die "Failed to open '$wrk/$outDir/ovljobs.dat'\n";
    foreach my $b (@bat) { print SUB "$b "; }  print SUB "\n";
    foreach my $b (@job) { print SUB "$b "; }  print SUB "\n";
    close(SUB);


    my $jobs = scalar(@opt);

    #  Submit to the grid (or tell the user to do it), or just run
    #  things here
    #
    if ($useGrid) {
        my $SGE;
        $SGE .= "\n";
        $SGE .= "qsub -p 0 -hard -l num_proc=2 -r y -N ovl_${asm} \\\n";
        $SGE .= "  -t 1-$jobs \\\n";
        $SGE .= "  -o $wrk/$outDir/overlap.\\\$TASK_ID.out \\\n";
        $SGE .= "  -e $wrk/$outDir/overlap.\\\$TASK_ID.err \\\n";
        $SGE .= "  $wrk/$outDir/overlap.sh\n";

        pleaseExecute($SGE);
        touch("$wrk/$outDir/jobsCreated.success");
        exit(0);
    } else {
        for (my $i=1; $i<=$jobs; $i++) {
            my $out = substr("0000" . $i, -4);
            $ENV{'SGE_TASK_ID'} = $i;
            if (runCommand("$wrk/$outDir/overlap.sh > $wrk/$outDir/$out.out 2>&1")) {
                print STDERR "Failed job $i\n";
                exit(1);
            }
            delete $ENV{'SGE_TASK_ID'};
        }

        touch("$wrk/$outDir/jobsCreated.success");
    }
}

1;