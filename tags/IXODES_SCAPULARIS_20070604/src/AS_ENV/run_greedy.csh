#!/usr/local/bin/tcsh -efx

# location of assembly bin directory (and extreme assembly helper scripts)
set asmBin=$1

set noFrgStore=1

#usage: readlist_easm.csh <bin directory> <list of read UIDs to use as seeds> <store prefix> [<cutoff > 0 and < 1>] [<list of read UIDs to extend with>]

set inputlist=$2
set inputfile=`echo $inputlist | awk -F/ '{print $NF}'`

set prefix=$3


## N.B. You probably want to have build the overlap store with the
## "metagenomic" setting of ERR_MODEL_IN_AS_GLOBAL_H -- i.e. set this
## to 10 or 15 or 20; above 25 not recommended.
## If you build your binaries this way, and then build an ovlStore
## you will have overlaps of up to 10%, 15% etc.

# maximal overlap error to follow; 0.15 = 15%
set cutoff=$4
if( $cutoff == "" ) then
  set cutoff=0.15
endif

set restrictToList=$5


set minasmlen=5000



### UID-to-IID mapping

if ( $noFrgStore == 1 ) then 
 ${asmBin}/gatekeeper -tabular -uid ${inputlist} -dumpfragments ${prefix}.gkpStore |\
   grep -v "UID" | awk '{print $2,$1}' \
 > ${inputfile}.iid2uid
 awk '{print $1}' ${inputfile}.iid2uid > ${inputfile}.iids
else
 ${asmBin}/fraguid2iid -s ${prefix}.frgStore -i ${inputlist}  > ${inputfile}.iids
endif


if($restrictToList != "") then
  cat ${restrictToList}  ${prefix}.IID2UIDwstatus.map |\
    awk 'NF==1{good[$1]=1;next}good[$2]==1{print $1,$2}' \
  > ${inputfile}.restrictediid2uid

  awk '{print $1}' ${inputfile}.restrictediid2uid > ${inputfile}.restrictediids
  set restriction="-i ${intputfile}.restrictediids"
else
  set restriction=""
endif

endif


### THE CORE EXECUTABLE ; output is the layout of fragments in "contigs"
if ( $noFrgStore == 1 ) then 
$asmBin/greedyFragmentTiling -e $cutoff -Q -N 40 -m 50 -R -P \
  -g ${prefix}.gkpStore \
  -o ${prefix}.ovlStore \
  -I ${inputfile}.iids \
  $restriction \
> ${inputfile}_e${cutoff}_m50_Q_N40.layout
else
$asmBin/greedyFragmentTiling -e $cutoff -Q -N 40 -m 50 -R -P \
  -f ${prefix}.frgStore \
  -g ${prefix}.gkpStore \
  -o ${prefix}.ovlStore \
  -I ${inputfile}.iids \
  $restriction \
> ${inputfile}_e${cutoff}_m50_Q_N40.layout
endif

### Convert the "layout" file into a set of IUM messages
cat  ${inputfile}_e${cutoff}_m50_Q_N40.layout |\
    ${asmBin}/greedy_layout_to_IUM \
  >  ${inputfile}_e${cutoff}_m50_Q_N40.cgb

### compute consensus sequence
if ( $noFrgStore == 1 ) then
  $asmBin/consensus -P -U -G ${prefix}.gkpStore ${inputfile}_e${cutoff}_m50_Q_N40.cgb
else
  $asmBin/consensus -P -U -G ${prefix}.frgStore ${inputfile}_e${cutoff}_m50_Q_N40.cgb
endif

### run terminator to assign UIDs
if ( $noFrgStore == 1 ) then
  cat ${inputfile}_e${cutoff}_m50_Q_N40.cgi | $asmBin/terminator -g ${prefix}.gkpStore -o ${inputfile}_e${cutoff}_m50_Q_N40.asm -m junk_mappings
else
   $asmBin/terminator -P -u -N -f ${prefix}.frgStore -g ${prefix}.gkpStore -i ${inputfile}_e${cutoff}_m50_Q_N40.cgi -o ${inputfile}_e${cutoff}_m50_Q_N40.asm -m junk_mappings
endif

\rm junk_mappings.*

### extract fasta sequence
cat ${inputfile}_e${cutoff}_m50_Q_N40.asm | ${asmBin}/utg2fasta ${inputfile}.iid2uid $minasmlen > ${inputfile}_e${cutoff}_m50_Q_N40.mfa

