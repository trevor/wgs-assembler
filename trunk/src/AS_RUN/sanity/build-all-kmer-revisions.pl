#!/usr/bin/perl



open(F, "< /work/NIGHTLY/kmer-svn/db/current");
my $latest = <F>;
chomp $latest;
close(F);


for (my $i=1917; $i<=$latest; $i++) {

    if (! -d "kmer$i") {
        print "Check out r$i\n";
        system("mkdir kmer$i");
        system("cd kmer$i && svn co -r $i file:///work/NIGHTLY/kmer-svn/trunk . > kmer-checkout.err 2>&1");
    }

    if (! -d "kmer$i/FreeBSD-amd64") {
        print "Compile r$i\n";
        system("cd kmer$i && gmake install > kmer-build.err 2>& 1 &");
    }
}
