#!/usr/bin/perl

use strict;
use Config;  #  for @signame

#  The only two globals:
#
my $wrkdir  = "/usr/local/projects/assembly_test/NIGHTLY";
my $wgscvs  = "/usr/local/projects/assembly_test/NIGHTLY/wgs-assembler-cvs";
my $kmersvn = "/usr/local/projects/assembly_test/NIGHTLY/kmer-svn";

#  Command line options are
#
#    'oper'   -- thing to do.
#
#    'ddir' -- properly formatted directory name to do it in.  op ==
#    checkout creates this directory.  If still supplied, it will
#    force a checkout from that date.
#
#  rsync
#     update the local repository in 'wgs-assembler-cvs'.
#
#  checkout date
#     Date must be of the form yyyy-mm-dd-hhmm.  This will checkout
#     the assembler current at that time.  If not present, the current
#     time is assumed.
#
#  build date
#     Builds the assembler checked out into 'date'.
#
#  assemble date specFile ....
#     Launches runCA on each specFile, running it in directory 'date'.
#     Each specFile must have a unique name; the assembly is named
#     after the specFile.
#
#  summarize date
#     Build a summary message for the directory 'date'.
#

{
    my $oper   = shift @ARGV;
    my $ddir   = shift @ARGV;

    my ($thisdate, $lastdate) = parseDate($ddir);

    if ($oper eq "rsync") {
        system("mkdir $wgscvs")  if (! -d "$wgscvs");
        system("mkdir $kmersvn") if (! -d "$kmersvn");

        system("cd $wgscvs  && rsync -av rsync://wgs-assembler.cvs.sourceforge.net/cvsroot/wgs-assembler/\* . > rsync.out 2>&1");
        system("cd $kmersvn && rsync -av kmer.svn.sourceforge.net::svn/kmer/\* . > rsync.out 2>&1");
    }


    if ($oper eq "checkout") {
        checkoutAndLogKmer($thisdate, $lastdate);
        checkoutAndLogCA($thisdate, $lastdate);
    }

    if ($oper eq "build") {
        buildKmer($thisdate);
        buildCA($thisdate);
    }

    if ($oper eq "assemble") {
        assemble($thisdate, $lastdate, @ARGV);
    }

    if (!defined($oper)) {
    }

    exit;
}


sub parseDate ($) {
    my $pathdate = shift @_;
    my $thisdate;
    my $lastdate;

    #  Until our Makefile is fixed, the ":" cannot be used, and so we need
    #  to adopt a slightly goofy timestamp name.

    if (defined($pathdate)) {
        if ($pathdate =~ m/^\d\d\d\d-\d\d-\d\d-\d\d\d\d$/) {
            $thisdate = $pathdate;
        } else {
            die "Malformed ddir '$pathdate'\n";
        }
    } else {
        my @v = localtime;
        $v[5] += 1900;
        $v[4]++;

        $v[5] = substr("0000$v[5]", -4);
        $v[4] = substr("0000$v[4]", -2);
        $v[3] = substr("0000$v[3]", -2);
        $v[2] = substr("0000$v[2]", -2);
        $v[1] = substr("0000$v[1]", -2);

        $thisdate = "$v[5]-$v[4]-$v[3]-$v[2]$v[1]";
    }

    open(F, "ls $wrkdir |");
    while (<F>) {
        chomp;
        if (m/^\d\d\d\d-\d\d-\d\d-\d\d\d\d$/) {
            if ($_ ne $thisdate) {
                $lastdate = $_;
            }
        }
    }
    close(F);

    print STDERR "Working on '$thisdate' -- last found is '$lastdate'.\n";

    return($thisdate, $lastdate);
}


sub checkoutAndLogKmer ($$) {
    my $thisdate = shift @_;
    my $lastdate = shift @_;

    print STDERR "Checking out $wrkdir/$thisdate (KMER)\n";

    if (-d "$wrkdir/$thisdate/wgs/kmer") {
        print STDERR "$wrkdir/$thisdate/wgs/kmer already exists.  Please remove to rerun.\n";
        return;
    }

    system("mkdir $wrkdir/$thisdate")     if (! -d "$wrkdir/$thisdate");
    system("mkdir $wrkdir/$thisdate/wgs") if (! -d "$wrkdir/$thisdate/wgs");

    #  Convert time-names to dates that svn can use.

    my $thisdatesvn;
    my $lastdatesvn;

    my $tz = "-05:00";

    if ($thisdate =~ m/(\d\d\d\d)-(\d\d)-(\d\d)-(\d\d)(\d\d)/) {
        $thisdatesvn = "$1-$2-$3T$4:$5$tz";
    } else {
    }

    if (defined($lastdate)) {
        if ($lastdate =~ m/(\d\d\d\d)-(\d\d)-(\d\d)-(\d\d)(\d\d)/) {
            $lastdatesvn = "$1-$2-$3T$4:$5$tz";
        } else {
        }
    }

    system("cd $wrkdir/$thisdate/wgs && svn co  -r \"{$thisdatesvn}\" file://$kmersvn/trunk kmer > kmer.checkout.err 2>&1");

    if ($lastdate ne "") {
        system("cd $wrkdir/$thisdate/wgs && svn log -r \"{$lastdatesvn}:{$thisdatesvn}\" kmer > kmer.updates");
    }

    print STDERR "$thisdate checked out!\n";
}

sub checkoutAndLogCA ($$) {
    my $thisdate = shift @_;
    my $lastdate = shift @_;

    print STDERR "Checking out $wrkdir/$thisdate (CA)\n";

    if (-d "$wrkdir/$thisdate/wgs/src") {
        print STDERR "$wrkdir/$thisdate/wgs/src already exists.  Please remove to rerun.\n";
        return;
    }

    system("mkdir $wrkdir/$thisdate")     if (! -d "$wrkdir/$thisdate");
    system("mkdir $wrkdir/$thisdate/wgs") if (! -d "$wrkdir/$thisdate/wgs");

    #  Convert time-names to dates that cvs can use.

    my $thisdatecvs;
    my $lastdatecvs;

    my $tz = "EDT";

    if ($thisdate =~ m/(\d\d\d\d)-(\d\d)-(\d\d)-(\d\d)(\d\d)/) {
        $thisdatecvs = "$1-$2-$3 $4:$5 $tz";
    } else {
    }

    if (defined($lastdate)) {
        if ($lastdate =~ m/(\d\d\d\d)-(\d\d)-(\d\d)-(\d\d)(\d\d)/) {
            $lastdatecvs = "$1-$2-$3 $4:$5 $tz";
        } else {
        }
    }

    #  Add -R to cvs options, for read-only repository; breaks on jcvi
    #  cvs; this seems to be a BSD extension?

    system("cd $wrkdir/$thisdate/wgs && cvs -r -d $wgscvs -z3 co  -N    -D '$thisdatecvs' src > src.checkout.err 2>&1");

    if ($lastdate ne "") {
        system("cd $wrkdir/$thisdate/wgs && cvs -r -d $wgscvs -z3 log -N -S -d '$lastdatecvs<$thisdatecvs' src > src.updates.raw");

        my $log;
        my $revs;
        my $date;
        my $auth;
        my $chng;
        my $file;
        my $inlog;
        my %logs;
        my %logdate;

        open(F, "< $wrkdir/$thisdate/wgs/src.updates.raw");
        while (<F>) {
            if (m/^==========.*========$/) {
                $inlog = 0;
                if (!defined($logs{$log})) {
                    $logs{$log}  = "$date\t$revs\t$chng\t$auth\t$file\n";
                } else {
                    $logs{$log} .= "$date\t$revs\t$chng\t$auth\t$file\n";
                }
                $logdate{"$date\0$log"}++;
                undef $log;
            }
            if (m/^----------------------------$/) {
                $inlog = 0;
                if (!defined($logs{$log})) {
                    $logs{$log}  = "$date\t$revs\t$chng\t$auth\t$file\n";
                } else {
                    $logs{$log} .= "$date\t$revs\t$chng\t$auth\t$file\n";
                }
                $logdate{"$date\0$log"}++;
                undef $log;
            }
            if ($inlog) {
                $log .= $_;
            }
            if (m/^RCS\s+file:\s+(.*),v/) {
                $file = $1;
                $file =~ s!$wgscvs/src/!!g;
            }
            if (m/^revision\s+(\d+.\d+)/) {
                $revs = $1;
            }
            if (m/^date:\s+(.*);\s+author:\s+(.*);\s+state.*lines:\s+(.*)\s*$/) {
                $date = $1;
                $auth = $2;
                $chng = $3;
                $inlog = 1;
            }
        }
        close(F);


        open(F, "> $wrkdir/$thisdate/wgs/src.updates");
        my @keys = sort keys %logdate;
        foreach my $l (@keys) {
            my ($d, $l) = split '\0', $l;
            
            if ((defined($logs{$l})) && (length($logs{$l} > 0))) {
                print F "----------------------------------------\n";
                print F "$logs{$l}\n";
                print F "$l\n";
                undef $logs{$l};
            }
        }
        close(F);
    }

    print STDERR "$thisdate checked out!\n";
}


sub buildKmer ($) {
    my $thisdate = shift @_;

    print STDERR "Building KMER $wrkdir/$thisdate\n";

    if (-e "$wrkdir/$thisdate/wgs/kmer/make.err") {
        print STDERR "$wrkdir/$thisdate/wgs/kmer was already built once (successuflly or not).  Please cleanup first.\n";
    } else {
        system("cd $wrkdir/$thisdate/wgs/kmer && sh configure.sh > configure.out 2> configure.err");
        system("cd $wrkdir/$thisdate/wgs/kmer && gmake install   > make.out.raw  2> make.err.raw");

        my %lines;

        open(F, "< $wrkdir/$thisdate/wgs/kmer/make.err.raw");
        open(G, "> $wrkdir/$thisdate/wgs/kmer/make.err");
        while (<F>) {
            chomp;
            next if (m/^ar:\s+creating/);
            if (!exists($lines{$_})) {
                $lines{$_}++;
                print G "$_\n";
            }
        }
        close(F);
        close(G);
    }
}


sub buildCA ($) {
    my $thisdate = shift @_;

    print STDERR "Building CA $wrkdir/$thisdate\n";

    if (-e "$wrkdir/$thisdate/wgs/src/make.err") {
        print STDERR "$wrkdir/$thisdate/wgs/src was already built once (successuflly or not).  Please cleanup first.\n";
    } else {
        system("cd $wrkdir/$thisdate/wgs/src && gmake > make.out.raw 2> make.err.raw");

        my %lines;

        open(F, "< $wrkdir/$thisdate/wgs/src/make.err.raw");
        open(G, "> $wrkdir/$thisdate/wgs/src/make.err");
        while (<F>) {
            chomp;
            next if (m/^ar:\s+creating/);
            if (!exists($lines{$_})) {
                $lines{$_}++;
                print G "$_\n";
            }
        }
        close(F);
        close(G);
    }
}



sub assemble ($$@) {
    my $thisdate  = shift @_;
    my $lastdate  = shift @_;
    my $assemblies;
    my $holds;

    my $arch;

    #  Ripped from runCA/util.pl
    my $syst = `uname -s`;    chomp $syst;  #  OS implementation
    my $arch = `uname -m`;    chomp $arch;  #  Hardware platform
    my $name = `uname -n`;    chomp $name;  #  Name of the system

    $arch = "amd64"  if ($arch eq "x86_64");
    $arch = "ppc"    if ($arch eq "Power Macintosh");
    #  end of rip.


    #  Figure out which are specfiles and which are email addresses.
    #
    my @spec;
    my @addresses;
    my $addresses;

    foreach my $arg (@_) {
        if ($arg =~ m/\@/) {
            push @addresses, $arg;
        } else {
            push @spec, $arg;
        }
    }

    $addresses = join ',', @addresses;

    open(F, "> $wrkdir/$thisdate/asm-done.sh");
    print F "#!/bin/sh\n";
    print F "\n";
    print F "perl $wrkdir/sanity-asm-done.pl JUNK $wrkdir \$1 $thisdate $lastdate ref\n";
    close(F);

    open(F, "> $wrkdir/$thisdate/all-done.sh");
    print F "#!/bin/sh\n";
    print F "\n";
    print F "perl $wrkdir/sanity-all-done.pl JUNK $wrkdir $thisdate $lastdate $addresses | /usr/sbin/sendmail -t -F CAtest\n";
    close(F);

    foreach my $s (@spec) {
        my ($n, $p) = split '\.', $s;

        print STDERR "----------------------------------------\n";
        print STDERR "Submitting assembly '$n'.\n";

        system("mkdir $wrkdir/$thisdate/$n") if (! -d "$wrkdir/$thisdate/$n");

        my $jl = "CAini_${n}_$$";  #  Name of the launcher
        my $jn = "CAtst_${n}_$$";  #  Name of the asm-done

        open(F, "> $wrkdir/$thisdate/$n/launch-assembly.sh");
        print F "#!/bin/sh\n";
        print F "\n";
        print F "#  runCA checks if this is set to decide if it is on the grid or not.  We want\n";
        print F "#  to pretend we are NOT on the grid, so runCA will resubmit itself immediately.\n";
        print F "#\n";
        print F "unset SGE_TASK_ID\n";
        print F "\n";
        print F "perl $wrkdir/$thisdate/wgs/$syst-$arch/bin/runCA \\\n";
        print F "  sgePropagateHold=$jn useGrid=1 scriptOnGrid=1\\\n";
        print F "  -p $n -d $wrkdir/$thisdate/$n -s $wrkdir/$n.$p\n";
        print F "\n";
        print F "#  Once that runCA finishes, we've updated the hold on $jn, and so can release\n";
        print F "#  our user hold on it.\n";
        print F "#\n";
        print F "qrls -h u $jn\n";
        print F "\n";
        print F "\n";
        close(F);

        #  A separate script is used to launch runCA, so that we can get a user hold on it.  This is messy.
        #  1)  Submit the launcher, holding it.
        #  2)  Submit the asm-done, holding it too.  Make it also hold_jid on the launcher.
        #  3)  Release the launcher.
        #  4)  The launcher submits runCA, with sgePropagateHold.  When that runCA finishes, the assembly
        #      is now running, AND it has updated the hold_jid on asm-done.
        #  5)  We can now release the asm-done.
        #
        #  Steps 1 and 2 are done here, step 3 is at the very end, steps 4 and 5 are done in the launch-assembly.sh above.

        system("cd $wrkdir/$thisdate/$n && qsub -b n -cwd -j y -o $wrkdir/$thisdate/$n/launch-assembly.err -l fast -h               -N $jl $wrkdir/$thisdate/$n/launch-assembly.sh");
        system("cd $wrkdir/$thisdate/$n && qsub -b n -cwd -j y -o $wrkdir/$thisdate/$n/asm-done.err        -l fast -h -hold_jid $jl -N $jn $wrkdir/$thisdate/asm-done.sh $n");

        if (defined($holds)) {
            $holds .= ",$jn";
        } else {
            $holds  = "$jn";
        }

        if (defined($assemblies)) {
            $assemblies .= ",$jl";
        } else {
            $assemblies  = "$jl";
        }
    }

    system("cd $wrkdir/$thisdate && qsub -b n -cwd -j y -o $wrkdir/$thisdate/all-done.err -l fast -hold_jid $holds -N CAfin_$thisdate $wrkdir/$thisdate/all-done.sh");

    #  Now, release everything to run.

    system("qrls -h u $assemblies");
}

exit(0);