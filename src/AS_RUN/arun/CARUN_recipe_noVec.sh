## Assembly Script.
## By Jason Miller, October 2006
## Modified by Marwan Oweis
## Modified by Sergey Koren
## Based on a run of Brian's script (runCA-OBT.pl) on GBAF genome.
## noVec ver. 2.0
## CA    ver. 4.0

## Files Required
##REQUIRED PREFIX.frg

echo $AS_OVL_ERROR_RATE
echo $AS_CNS_ERROR_RATE
echo $AS_CGW_ERROR_RATE

#> Print Version
$CA_BIN/version

#> Gatekeeper

mkdir -p $WORKDIR/0-preoverlap && cd $WORKDIR/0-preoverlap
$CA_BIN/gatekeeper -e 10000000 -Q -T $AS_GKP_BELIEVE_STATS \
   -o $WORKDIR/$PREFIX.gkpStore \
   $WORKDIR/$PREFIX.frg \
   > $WORKDIR/0-preoverlap/gatekeeper.out \
   2> $WORKDIR/0-preoverlap/gatekeeper.err

cd $WORKDIR/0-preoverlap
$CA_BIN/initialTrim -update -q 12  \
   -log $PREFIX.initialTrimLog  \
   -frg $WORKDIR/$PREFIX.gkpStore  \
   > initialTrim.err \
   2>&1

#> Meryl

cd $WORKDIR/0-preoverlap
$CA_BIN/meryl -m $KMER -n $NMER -K 10 \
   -s $WORKDIR/$PREFIX.gkpStore \
   -o $PREFIX.${KMER}mers${NMER}.ovl.fasta \
   > meryl.out \
   2> meryl.err

$CA_BIN/meryl -m $KMER -n $NMER_OBT -K 10 \
   -s $WORKDIR/$PREFIX.gkpStore \
   -o $PREFIX.${KMER}mers${NMER_OBT}.obt.fasta \
   > meryl.out \
   2> meryl.err

#> Overlapper, first run  
mkdir -p $WORKDIR/0-overlap && cd $WORKDIR/0-overlap
$CA_BIN/gatekeeper -L $WORKDIR/$PREFIX.gkpStore  >lastfraginstore.out
HIFRAGID=`sed 's/Last frag in store is iid = //'  < lastfraginstore.out `
# Required options: -M, -h, -r, -o. Do partial overlaps (-G).
$CA_BIN/overlap -G -M $MEMORY -t 2 -h 1-$HIFRAGID -r 0-$HIFRAGID \
  -k ${KMER} -k $WORKDIR/0-preoverlap/$PREFIX.${KMER}mers${NMER_OBT}.obt.fasta \
  -o $PREFIX.ovl \
  $WORKDIR/$PREFIX.gkpStore  

#> Convert overlaps from text to binary
cd $WORKDIR/0-overlap
$CA_BIN/acceptableOBToverlap < $PREFIX.ovl > $PREFIX.ovb 

#> Create overlap store
mkdir -p $WORKDIR/0-trim && cd $WORKDIR/0-trim
$CA_BIN/overlapStore \
  -c $WORKDIR/$PREFIX.obtStore -M 1024 -m $HIFRAGID \
  $WORKDIR/0-overlap/$PREFIX.ovb  \
  > $PREFIX.ovl.sorted 2> $PREFIX.ovl.sorted.err

#> Consolidate overlaps.
cd $WORKDIR/0-trim
$CA_BIN/consolidate \
  -ovs $WORKDIR/$PREFIX.obtStore > $PREFIX.ovl.consolidated

#> Overlap-based trimming
cd $WORKDIR/0-trim
$CA_BIN/merge-trimming \
  -log $PREFIX.mergeLog \
  -frg $WORKDIR/$PREFIX.gkpStore \
  -ovl $PREFIX.ovl.consolidated \
  > $PREFIX.ovl.consolidated.err 2>&1

#> Overlap-based chimer detection
cd $WORKDIR/0-trim
$CA_BIN/chimera  \
  -gkp $WORKDIR/$PREFIX.gkpStore  \
  -ovs $WORKDIR/$PREFIX.obtStore \
  -summary $PREFIX.chimera.summary  \
  -report $PREFIX.chimera.report  \
  > $PREFIX.chimera.err 2>&1

#> Overlapper, second run  
mkdir -p $WORKDIR/1-overlap && cd $WORKDIR/1-overlap
$CA_BIN/gatekeeper -L $WORKDIR/$PREFIX.gkpStore  >lastfraginstore.out
HIFRAGID=`sed 's/Last frag in store is iid = //'  < lastfraginstore.out `
# Required options: -M, -h, -r, -o. No partial overlaps (-G).
$CA_BIN/overlap -M $MEMORY -t 2 -h 1-$HIFRAGID -r 0-$HIFRAGID \
  -k ${KMER} -k $WORKDIR/0-preoverlap/$PREFIX.${KMER}mers${NMER}.ovl.fasta \
  -o $PREFIX.ovb \
  $WORKDIR/$PREFIX.gkpStore  

#> Create overlap store
cd $WORKDIR/1-overlap
$CA_BIN/overlapStore \
  -c $WORKDIR/$PREFIX.ovlStore -M 1024 -m $HIFRAGID \
  $WORKDIR/1-overlap/$PREFIX.ovb  \
  > grow-olap-store.out 2> grow-olap-store.err

#> Verifying overlap store exists (fails when there are no overlaps)
OVL_EXIST=`ls $WORKDIR/$PREFIX.ovlStore | wc -l`
CONTINUE=true
if [ $OVL_EXIST == 0 ]; then \
  CONTINUE=false; \
fi
$CONTINUE

#> Fragment correction  
mkdir -p $WORKDIR/2-frgcorr && cd $WORKDIR/2-frgcorr
# threads: -t 2
$CA_BIN/correct-frags \
  -S $WORKDIR/$PREFIX.ovlStore -o $PREFIX.frgcorr  \
  $WORKDIR/$PREFIX.gkpStore  1  $HIFRAGID  > correct-frags.err 2>&1 
$CA_BIN/cat-corrects  -o $PREFIX.corr  \
  $PREFIX.frgcorr > cat-corrects.out 2> cat-corrects.err

#> Overlap error correction

mkdir -p $WORKDIR/3-ovlcorr && cd $WORKDIR/3-ovlcorr 
$CA_BIN/correct-olaps \
  -S $WORKDIR/$PREFIX.ovlStore   -e $PREFIX.erate \
  $WORKDIR/$PREFIX.gkpStore  $WORKDIR/2-frgcorr/$PREFIX.corr  \
  1  $HIFRAGID  > correct-olaps.err 2>&1
# In single-partition mode, move file to simulate cat-erates
cd $WORKDIR/3-ovlcorr && mv $PREFIX.erate $PREFIX.erates

cd $WORKDIR/3-ovlcorr
# only update the erates when the previous step generated an erates files, otherwise this step will fail
# the overlapStore now handles updating the errates, through the -u option (which is not listed in the docs)
ERATES_EXIST=`ls $PREFIX.erates | wc -l`
if [ $ERATES_EXIST != 0 ] ; then \
   $CA_BIN/overlapStore \
      -u $WORKDIR/$PREFIX.ovlStore \
      $PREFIX.erates > update-erates.err 2>&1; \
fi

mkdir -p $WORKDIR/4-unitigger && cd $WORKDIR/4-unitigger
# Note Brian uses -B 75000. Code supports B=0 => single partition.
# The -e 15 parameter equates to 1.5% assumed sequencing error.
# This generates $PREFIX.cgb file.
cd $WORKDIR/4-unitigger 
$CA_BIN/unitigger -k -d 1 -x 1 -z 10 -j 5 $GENOMELENGTH -U $BUBBLE -e $ERATE \
  -F $WORKDIR/$PREFIX.gkpStore -o $PREFIX  \
  -I $WORKDIR/$PREFIX.ovlStore  \
  > unitigger.out  2> unitigger.err 
cp $WORKDIR/4-unitigger/$PREFIX.cga.0 $WORKDIR/.

#> Consensus on unitigs

# For single-partition case, we do not call partitionFragStore.
# Note Brian uses -S for frag store partition.
mkdir -p $WORKDIR/5-consensus && cd $WORKDIR/5-consensus
$CA_BIN/consensus -G -m -U   -o $PREFIX.cgi \
  $WORKDIR/$PREFIX.gkpStore  $WORKDIR/4-unitigger/$PREFIX.cgb  \
  > consensus.err 2>&1
mv $PREFIX.cgi $WORKDIR
ln -s $WORKDIR/$PREFIX.cgi .

#> Preliminary scaffolds, estimate mate distances, build SeqStore

mkdir -p $WORKDIR/6-distances && cd $WORKDIR/6-distances
# cgw options: -c for checkpoints, -j and -k for Astat cutoffs,
# -r for rocks, -s for stones, # -T to ignore transchunks, -P for text output.
# The -o creates $PROJECT.SeqStore directory.
# The -G means no intermediate results and the -S 0 means no resolveSurrogates
$CA_BIN/cgw  -c -j $ASTATLOW -k $ASTATHIGH -r 5 -s 2 -S 0 -G -z  \
  -m $AS_CGW_SAMPLE_SIZE \
  -g $WORKDIR/$PREFIX.gkpStore  \
  -o $PREFIX  $WORKDIR/$PREFIX.cgi  > cgw.out 2>&1
# SeqStore needed by subsequent cgw calls. Move it up one level.
mv $PREFIX.SeqStore $WORKDIR
ln -s $WORKDIR/$PREFIX.SeqStore .

#> Store distance estimates with gatekeeper
cd $WORKDIR/6-distances
$CA_BIN/gatekeeper  -a -D \
  -o $WORKDIR/$PREFIX.gkpStore  \
  $WORKDIR/6-distances/stat/scaffold_final.distupdate.dst \
  $WORKDIR/6-distances/stat/contig_final.distupdate.dst  \
  > gatekeeper.err 2>&1

#> Initial scaffolds

mkdir -p $WORKDIR/7-0-scaffold && cd $WORKDIR/7-0-scaffold
ln -s $WORKDIR/$PREFIX.SeqStore .
$CA_BIN/cgw  -c -j $ASTATLOW -k $ASTATHIGH -r 5 -s 0 \
  -m $AS_CGW_SAMPLE_SIZE \
  -S 0 -G -z \
  -g $WORKDIR/$PREFIX.gkpStore  \
  -o $PREFIX    $WORKDIR/$PREFIX.cgi  \
  > cgw.err 2>&1

# Begin : determine max checkpoint number
ls -1 $PREFIX.ckp.*   > pipe01
sed 's/.ckp./ /'      < pipe01 > pipe02 
awk '{print $2 ;}'    < pipe02 > pipe03
sort -nr              < pipe03 > pipe04
MAXCHKPT=`head -n 1 pipe04`
rm     pipe01 pipe02 pipe03 pipe04
# End : determine max checkpoint number

#> Initial extend clear ranges

mkdir -p $WORKDIR/7-1-extend && cd $WORKDIR/7-1-extend
ln -s $WORKDIR/$PREFIX.SeqStore .
ln -s $WORKDIR/7-0-scaffold/$PREFIX.ckp.$MAXCHKPT .
# Brian specifies a range of scaffolds (-b and -e) but we do not.
# Next command reads ckp files and SeqStore directory.
$CA_BIN/extendClearRanges  -c $PREFIX  -n $MAXCHKPT    \
  -g $WORKDIR/$PREFIX.gkpStore  \
  > extendClearRanges.err 2>&1
assert_exists $PREFIX.ckp.* > /dev/null
MAXCHKPT=`ls -1 $PREFIX.ckp.* | sed 's/.ckp./ /' | awk '{print $2 ;}' | sort -nr | head -n 1`

#> Secondary scaffolds

mkdir -p $WORKDIR/7-2-scaffold && cd $WORKDIR/7-2-scaffold 
ln -s $WORKDIR/$PREFIX.SeqStore .
ln -s $WORKDIR/7-1-extend/$PREFIX.ckp.$MAXCHKPT .
# The -R option says to start from previous checkpoint.
# We use -N 3 for all cgw runs except first and last.
$CA_BIN/cgw -y -N 3 -c -j $ASTATLOW -k $ASTATHIGH -r 5 -s 0 -P \
  -m $AS_CGW_SAMPLE_SIZE \
  -S 0 -G -z \
  -R $MAXCHKPT   \
  -g $WORKDIR/$PREFIX.gkpStore   \
  -o $PREFIX  $WORKDIR/$PREFIX.cgi  > cgw.err 2>&1

# Begin : determine max checkpoint number
ls -1 $PREFIX.ckp.*   > pipe01
sed 's/.ckp./ /'      < pipe01 > pipe02 
awk '{print $2 ;}'    < pipe02 > pipe03
sort -nr              < pipe03 > pipe04
MAXCHKPT=`head -n 1 pipe04`
rm     pipe01 pipe02 pipe03 pipe04
# End : determine max checkpoint number

#> Secondary extend clear ranges

mkdir -p $WORKDIR/7-3-extend && cd $WORKDIR/7-3-extend
ln -s $WORKDIR/$PREFIX.SeqStore .
ln -s $WORKDIR/7-2-scaffold/$PREFIX.ckp.$MAXCHKPT .
# Brian specifies a range of scaffolds (-b and -e) but we do not.
# Next command reads ckp files and SeqStore directory.
# Next command generates another checkpoint file, not used.
$CA_BIN/extendClearRanges  -c $PREFIX  -n $MAXCHKPT     \
  -g $WORKDIR/$PREFIX.gkpStore   \
  > extendClearRanges.err 2>&1

# Begin : determine max checkpoint number
ls -1 $PREFIX.ckp.*   > pipe01
sed 's/.ckp./ /'      < pipe01 > pipe02 
awk '{print $2 ;}'    < pipe02 > pipe03
sort -nr              < pipe03 > pipe04
MAXCHKPT=`head -n 1 pipe04`
rm     pipe01 pipe02 pipe03 pipe04
# End : determine max checkpoint number

#> Tertiary scaffolds

mkdir -p $WORKDIR/7-4-scaffold && cd $WORKDIR/7-4-scaffold 
ln -s $WORKDIR/$PREFIX.SeqStore .
ln -s $WORKDIR/7-3-extend/$PREFIX.ckp.$MAXCHKPT .
# The -R option says to start from previous checkpoint.
# We use -N 3 for all cgw runs except first and last.
$CA_BIN/cgw -y -N 3 -c -j $ASTATLOW -k $ASTATHIGH -r 5 -s 2 -z \
  -m $AS_CGW_SAMPLE_SIZE \
  -R $MAXCHKPT   \
  -g $WORKDIR/$PREFIX.gkpStore  \
  -o $PREFIX  $WORKDIR/$PREFIX.cgi  > cgw.err 2>&1

# Begin : determine max checkpoint number
ls -1 $PREFIX.ckp.*   > pipe01
sed 's/.ckp./ /'      < pipe01 > pipe02 
awk '{print $2 ;}'    < pipe02 > pipe03
sort -nr              < pipe03 > pipe04
MAXCHKPT=`head -n 1 pipe04`
rm     pipe01 pipe02 pipe03 pipe04
# End : determine max checkpoint number

#> Consensus on scaffolds

ln -s $WORKDIR/7-4-scaffold/ $WORKDIR/7-scaffold
mkdir -p $WORKDIR/8-consensus && cd $WORKDIR/8-consensus 
ln -s $WORKDIR/7-scaffold/$PREFIX.cgw_contigs .
# Brian uses -p and -S and -m to specify partition, but we do not.
$CA_BIN/consensus \
  -s $WORKDIR/$PREFIX.SeqStore  -V $MAXCHKPT    \
  -m -o $PREFIX.cns_contigs  $WORKDIR/$PREFIX.gkpStore   \
  $PREFIX.cgw_contigs > consensus.err 2>&1 

#> Terminator
#-B $EUIDBlockSize -n $EUIDNamespace
mkdir -p $WORKDIR/9-terminator && cd $WORKDIR/9-terminator 
assert_exists $WORKDIR/7-scaffold/$PREFIX.cgw  \
    $WORKDIR/8-consensus/$PREFIX.cns_contigs   \
    $WORKDIR/7-scaffold/$PREFIX.cgw_scaffolds > /dev/null
cat $WORKDIR/7-scaffold/$PREFIX.cgw  \
    $WORKDIR/8-consensus/$PREFIX.cns_contigs   \
    $WORKDIR/7-scaffold/$PREFIX.cgw_scaffolds  \
  | $CA_BIN/terminator -E $EUIDList -n $EUIDNamespace  \
  -g $WORKDIR/$PREFIX.gkpStore  \
  -o $PREFIX.asm  -m $PREFIX.map  > terminator.err 2>&1 
cd $WORKDIR
ln -s 9-terminator/$PREFIX.asm .

#> Generate QC report
cd $WORKDIR
$CA_BIN/caqc.pl $GENOMELENGTH -metrics $PREFIX.asm > caqc.err 2>&1

#> Generate scaffold FASTA files
cd $WORKDIR
$CA_BIN/asmProcessScaffolds_TER -f   \
  $PREFIX.scaffold.fasta   < $PREFIX.asm
ln -s $PREFIX.scaffold.fasta $PREFIX.scaffolds.fasta

#> Generate singleton FASTA files
cd $WORKDIR
$CA_BIN/dumpSingletons  \
  -g $WORKDIR/$PREFIX.gkpStore   \
  -c $WORKDIR/7-scaffold/$PREFIX -n $MAXCHKPT -U -S  \
  > $PREFIX.singleton.fasta   2> dumpSingletons.err

#> Generate contig FASTA files (includes degenerates)
cd $WORKDIR
$CA_BIN/asmOutputContigsFasta -d  < $PREFIX.asm > $PREFIX.fasta

