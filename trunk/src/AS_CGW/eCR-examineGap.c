/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/


static const char *rcsid = "$Id: eCR-examineGap.c,v 1.28 2011-12-29 09:26:03 brianwalenz Exp $";
#include "eCR.h"

#include "GapWalkerREZ.h"  //  FindGapLength
#include "AS_UTL_reverseComplement.h"


// extern variables for controlling use of Local_Overlap_AS_forCNS

// initialized value is 12 -- no more than this many segments in the chain
extern int MaxGaps;

// init value is 200; this could be set to the amount you extend the clear
// range of seq b, plus 10 for good measure
extern int MaxBegGap;

// init value is 200; this could be set to the amount you extend the
// clear range of seq a, plus 10 for good measure
extern int MaxEndGap;

// initial value is 1000 (should have almost no effect) and defines
// the largest gap between segments in the chain
//
// Also: allowed size of gap within the alignment -- forcing
// relatively good alignments, compared to those allowed in
// bubble-smoothing where indel polymorphisms are expected
extern int MaxInteriorGap;

// boolean to cause the size of an "end gap" to be evaluated with
// regard to the clear range extension
extern int asymmetricEnds;


static int DefaultMaxBegGap;
static int DefaultMaxEndGap;
static int DefaultMaxGaps;
static int DefaultMaxInteriorGap;
static int DefaultAsymmetricEnds;

void
saveDefaultLocalAlignerVariables(void) {
  DefaultMaxBegGap      = MaxBegGap;
  DefaultMaxEndGap      = MaxEndGap;
  DefaultMaxGaps        = MaxGaps;
  DefaultMaxInteriorGap = MaxInteriorGap;
  DefaultAsymmetricEnds = asymmetricEnds;
}

void
restoreDefaultLocalAlignerVariables(void) {
  MaxBegGap      = DefaultMaxBegGap;
  MaxEndGap      = DefaultMaxEndGap;
  MaxGaps        = DefaultMaxGaps;
  MaxInteriorGap = DefaultMaxInteriorGap;
  asymmetricEnds = DefaultAsymmetricEnds;
}


int
examineGap(ContigT *lcontig, int lFragIid,
           ContigT *rcontig, int rFragIid,
           int gapNumber,
           int *ahang,
           int *olapLengthOut,
           int *bhang,
           int *currDiffs,
           int *lcontigBasesIntact,
           int *rcontigBasesIntact,
           int *closedGapDelta,
           int lBasesToNextFrag,
           int rBasesToNextFrag,
           int *leftFragFlapLength,
           int *rightFragFlapLength) {

  CIFragT *lFrag = NULL;
  CIFragT *rFrag = NULL;
  char lFragSeqBuffer[AS_READ_MAX_NORMAL_LEN+1];
  char rFragSeqBuffer[AS_READ_MAX_NORMAL_LEN+1];
  char lcompBuffer[AS_READ_MAX_NORMAL_LEN+CONTIG_BASES+1];
  char rcompBuffer[AS_READ_MAX_NORMAL_LEN+CONTIG_BASES+1];
  uint lclr_bgn=0, lclr_end=0;
  uint rclr_bgn=0, rclr_end=0;
  char *lSequence = NULL;
  char *rSequence = NULL;
  SequenceOrient lContigOrientation;
  SequenceOrient rContigOrientation;
  int lcontigBaseStart = 0;
  int lcontigBasesUsed = 0;
  int rcontigBasesUsed = 0;
  int lFragContigOverlapLength = 0;
  int rFragContigOverlapLength = 0;
  int i;

  static gkFragment fr;  //  static for performance only

  static VA_TYPE(char)           *lContigConsensus = NULL;
  static VA_TYPE(char)           *rContigConsensus = NULL;
  static VA_TYPE(char)           *lContigQuality = NULL;
  static VA_TYPE(char)           *rContigQuality = NULL;

  if (lContigConsensus == NULL) {
    lContigConsensus   = CreateVA_char(4096);
    rContigConsensus   = CreateVA_char(4096);
    lContigQuality     = CreateVA_char(4096);
    rContigQuality     = CreateVA_char(4096);
  }

#if 0
  if ((lFragIid == 746274) || (rFragIid == 1109314)) {
    debug.examineGapLV = 1;
    debug.examineGapFP = stderr;
  }
#endif

  // set some variables to control Local_Overlap_AS_forCNS
  //
  MaxGaps        = 5;
  //MaxBegGap    = 200;
  //MaxEndGap    = 200;
  MaxInteriorGap = 30;
  asymmetricEnds = TRUE;

  if (lcontig->offsetAEnd.mean < lcontig->offsetBEnd.mean)
    lContigOrientation.setIsForward();
  else
    lContigOrientation.setIsReverse();

  if (rcontig->offsetAEnd.mean < rcontig->offsetBEnd.mean)
    rContigOrientation.setIsForward();
  else
    rContigOrientation.setIsReverse();

  // Get the consensus sequences for both chunks from the Store

  GetConsensus(ScaffoldGraph->ContigGraph, lcontig->id, lContigConsensus, lContigQuality);
  GetConsensus(ScaffoldGraph->ContigGraph, rcontig->id, rContigConsensus, rContigQuality);

  lSequence = Getchar(lContigConsensus, 0);
  rSequence = Getchar(rContigConsensus, 0);

  // ----------------------> lContigOrientation == A_B
  //                  -----> frag is 5p->3p into gap, aligned with contig
  //
  // <---------------------- lContigOrientation == B_A
  //                  -----> frag is 5p->3p into gap, aligned opposite to contig
  //
  // the frag is oriented opposite to the contig in this case
  // flip contig sequence to its orientation in scaffold
  if (lContigOrientation.isReverse())
    reverseComplementSequence(lSequence, strlen(lSequence));

  // ----------------------> rContigOrientation == A_B
  // <-----                  frag is 5p->3p into gap, aligned opposite to contig
  //
  // <---------------------- rContigOrientation == B_A
  // <-----                  frag is 5p->3p into gap, aligned with contig
  //
  // the frag is oriented opposite to the contig in this case
  // flip contig sequence to its orientation in scaffold
  if (rContigOrientation.isReverse())
    reverseComplementSequence(rSequence, strlen(rSequence));

  if (lFragIid != -1) {
    lFrag = GetCIFragT(ScaffoldGraph->CIFrags, lFragIid);

    ScaffoldGraph->gkpStore->gkStore_getFragment(lFragIid, &fr, GKFRAGMENT_SEQ);
    fr.gkFragment_getClearRegion(lclr_bgn, lclr_end, AS_READ_CLEAR_ECR_0 + iterNumber - 1);
    assert(lclr_bgn < lclr_end);
    strcpy(lFragSeqBuffer, fr.gkFragment_getSequence());
  }

  if (rFragIid != -1) {
    rFrag = GetCIFragT(ScaffoldGraph->CIFrags, rFragIid);

    ScaffoldGraph->gkpStore->gkStore_getFragment(rFragIid, &fr, GKFRAGMENT_SEQ);
    fr.gkFragment_getClearRegion(rclr_bgn, rclr_end, AS_READ_CLEAR_ECR_0 + iterNumber - 1);
    assert(rclr_bgn < rclr_end);
    strcpy(rFragSeqBuffer, fr.gkFragment_getSequence());
  }

  //  Always, we want to flip the right frag.
  //
  if (rFragIid != -1) {
    int temp = 0;
    int len  = strlen(rFragSeqBuffer);

    reverseComplementSequence(rFragSeqBuffer, len);

    temp     = len - rclr_bgn;
    rclr_bgn = len - rclr_end;
    rclr_end = temp;
  }

  if (debug.examineGapLV > 0) {
    if (lFragIid != -1)
      fprintf(debug.examineGapFP, "lFrag:%d clr:%d,%d len: "F_SIZE_T"\n%s\n", lFragIid, lclr_bgn, lclr_end, strlen(lFragSeqBuffer), lFragSeqBuffer);
    if (rFragIid != -1)
      fprintf(debug.examineGapFP, "rFrag:%d clr:%d,%d len: "F_SIZE_T"\n%s\n", rFragIid, rclr_bgn, rclr_end, strlen(rFragSeqBuffer), rFragSeqBuffer);
  }

  ////////////////////////////////////////
  //  Create the left sequence

  // we use frag sequence from where the clear range ends to the end of the frag
  // ----------------------> lContigOrientation == A_B
  //               ----->    frag is 5p->3p into gap, aligned with contig
  // <---------------------- lContigOrientation == B_A
  //               ----->    frag is 5p->3p into gap, aligned opposite to contig

  lFragContigOverlapLength = 0;

  if (lFragIid != -1) {
    if (lContigOrientation.isForward())
      lFragContigOverlapLength = (int) (lcontig->bpLength.mean - lFrag->contigOffset3p.mean);
    else
      lFragContigOverlapLength = (int) (lFrag->contigOffset3p.mean);
  }

  // grab the last CONTIG_BASES bases of the lcontig consensus sequence
  lcontigBasesUsed = MIN(CONTIG_BASES - lFragContigOverlapLength,
                         (int) lcontig->bpLength.mean - lFragContigOverlapLength);

  // lcontigBaseStart is the base where we start using the consensus
  // sequence in lcompBuffer and thus also the number of bases from
  // the contig that are intact
  //
  lcontigBaseStart = strlen(lSequence) - lcontigBasesUsed - lFragContigOverlapLength;

  // grab the bases from the contig, ie, those not from the non-clear
  // range of the frag
  //
  for (i = 0; i < lcontigBasesUsed; i++)
    lcompBuffer[ i ] = lSequence[ lcontigBaseStart + i];
  lcompBuffer[ i ] = 0;

  if (debug.examineGapLV > 0) {
    fprintf(debug.examineGapFP, "lcompBuffer:  len:"F_SIZE_T"   lFragContigOverlapLength:%d lcontigBaseStart:%d lcontigBasesUsed:%d\n",
            strlen(lcompBuffer), lFragContigOverlapLength, lcontigBaseStart, lcontigBasesUsed);
  }

  // now tack on the 3p clr range extension to the bases of the contig
  // consensus sequence

  MaxEndGap = 100;

  if (lFragIid != -1) {

    // basesToNextFrag is the number of bases back to the first frag that gets us to 2x
    //
    MaxEndGap = strlen(lFragSeqBuffer) - lclr_end + lBasesToNextFrag + 20;  // 20 is slop

    for (i = lclr_end; i < strlen(lFragSeqBuffer); i++)
      lcompBuffer[ lcontigBasesUsed + i - lclr_end ] = lFragSeqBuffer[ i ];
    lcompBuffer[ lcontigBasesUsed + i - lclr_end ] = 0;
  }



  ////////////////////////////////////////
  //  Create the right sequence


  // we use frag sequence from where the clear range ends to the end of the frag
  // ----------------------> rContigOrientation == A_B
  //    <-----               frag is 5p->3p into gap, aligned opposite to contig

  // <---------------------- rContigOrientation == B_A
  //      <-----             frag is 5p->3p into gap, aligned with contig

  rFragContigOverlapLength = 0;

  if (rFragIid != -1) {
    if (rContigOrientation.isForward())
      rFragContigOverlapLength = (int) (rFrag->contigOffset3p.mean);
    else
      rFragContigOverlapLength = (int) (rcontig->bpLength.mean - rFrag->contigOffset3p.mean);
  }

  MaxBegGap = 100;

  if (rFragIid != -1) {

    // basesToNextFrag is the number of bases back to the first frag that gets us to 2x
    //
    MaxBegGap = rclr_bgn + rBasesToNextFrag + 20;  // 20 is slop

    // now if we have a right frag, grab the "5p" clr range extension
    // - remember the frag has been flipped

    for (i = 0; i < rclr_bgn; i++)
      rcompBuffer[ i ] = rFragSeqBuffer[ i ];
    rcompBuffer[i] = 0;
  }

  // grab the first CONTIG_BASES bases of the rcontig consensus
  // sequence the rcontig consensus has been flipped if necessary

  rcontigBasesUsed = MIN(CONTIG_BASES - rFragContigOverlapLength,
                         (int) rcontig->bpLength.mean - rFragContigOverlapLength);

  for (i = 0; i < rcontigBasesUsed; i++)
    rcompBuffer[ rclr_bgn + i ] = rSequence[ i + rFragContigOverlapLength ];
  rcompBuffer[ rclr_bgn + i ] = 0;

  if (debug.examineGapLV > 0) {
    fprintf(debug.examineGapFP, "> lcompBuffer gap %d (len: " F_SIZE_T "): \n%s\n", gapNumber, strlen(lcompBuffer), lcompBuffer);
    fprintf(debug.examineGapFP, "> rcompBuffer gap %d (len: " F_SIZE_T "): \n%s\n", gapNumber, strlen(rcompBuffer), rcompBuffer);
  }
  assert(strlen(lcompBuffer)>0);
  assert(strlen(rcompBuffer)>0);

  // now lcompBuffer and rcompBuffer hold the sequence of the fragments in the correct strand
  // now prepare for call to Local_Overlap_AS_forCNS

  ALNoverlap *overlap;
  int beg, end, opposite = FALSE;
  double erate, thresh, minlen;
  CompareOptions what;
  LengthT gapSize;

  beg    = -strlen (rcompBuffer);
  end    = strlen (lcompBuffer);
  assert((0.0 <= AS_CGW_ERROR_RATE) && (AS_CGW_ERROR_RATE <= AS_MAX_ERROR_RATE));
  erate  = AS_CGW_ERROR_RATE + 0.02;  //  Historically, erate == 0.12.
  thresh = CGW_DP_THRESH;
  minlen = CGW_DP_MINLEN + 10;        //  Historically, minlen == 30.  Not sure if we should hardcode 30 or offset from CGW_DP_MINLEN.
  what   = AS_FIND_LOCAL_ALIGN;

  if (debug.examineGapLV > 0)
    fprintf(debug.examineGapFP, "MaxBegGap: %d  MaxEndGap: %d", MaxBegGap, MaxEndGap);

  overlap = Local_Overlap_AS_forCNS(lcompBuffer,
                                    rcompBuffer,
                                    beg, end,   //  band
                                    -beg, end,  //  ahang,bhang not known
                                    opposite,
                                    erate,
                                    thresh,
                                    minlen,
                                    what);

  char *rcompBufferTrimmed = NULL;

  // not interested in overlaps with negative ahangs or bhangs
  //
  if ((overlap == NULL) || (overlap->begpos <= 0) || (overlap->endpos <= 0)) {
    if (debug.examineGapLV > 0)
      fprintf(debug.examineGapFP, "no overlap found between frags %d and %d\n", lFragIid, rFragIid);
    restoreDefaultLocalAlignerVariables();
    return 0;
  }


  if (debug.examineGapLV > 0) {
    fprintf(debug.examineGapFP, "initial ahang: %d, bhang:%d, length: %d, diffs: %d, diffs / length %%: %f\n",
            overlap->begpos, overlap->endpos, overlap->length, overlap->diffs,
            100.0 * overlap->diffs / overlap->length);
    PrintALNoverlap("", lcompBuffer, rcompBuffer, overlap);
    //debug.examineGapFP
  }

  /*
    \
    \
    \  left flap, right frag
    \
    -------------------
    -------------------
    \   right flap, left frag
    \
    \
  */

  // do flap trimming
  // left flap, right frag
  *rightFragFlapLength = 0;
  if (overlap->begpos == abs(overlap->trace[0]) - 1 && overlap->trace[0] < 0) {
    while (overlap->begpos == abs(overlap->trace[ *rightFragFlapLength ]) - 1 &&
           overlap->trace[ *rightFragFlapLength ] < 0)
      (*rightFragFlapLength)++;
  }

  // right flap, left frag
  {
    // first find the number of elements in trace
    int numTraceElements = 0, icnt = 0, rcompLength = 0;
    while (overlap->trace[ icnt++ ] != 0)
      numTraceElements++;

    *leftFragFlapLength = 0;
    if ((numTraceElements > 0) &&
        (overlap->endpos > 0) &&
        (overlap->trace[ numTraceElements - 1 ] > 0)) {
      icnt = numTraceElements - 1;
      rcompLength = strlen(rcompBuffer) - overlap->endpos + 1;
      while ((icnt >= 0) &&
             (overlap->trace[icnt] == rcompLength)) {
        (*leftFragFlapLength)++;
        icnt--;
      }
    }
  }

  lcompBuffer[ strlen(lcompBuffer) - *leftFragFlapLength ] = 0;
  rcompBufferTrimmed = &rcompBuffer[ *rightFragFlapLength ];

  // now do overlap again after trimming to make sure it is still
  // there, sometimes trimming makes them go away

  beg = -strlen(rcompBufferTrimmed);
  end = strlen(lcompBuffer);

  overlap = Local_Overlap_AS_forCNS(lcompBuffer,
                                    rcompBufferTrimmed,
                                    beg, end,
                                    -beg, end,
                                    opposite,
                                    erate,
                                    thresh,
                                    minlen,
                                    what);


  // not interested in overlaps with negative ahangs or bhangs
  //
  if ((overlap == NULL) || (overlap->begpos <= 0) || (overlap->endpos <= 0)) {
    if (debug.examineGapLV > 0)
      fprintf(debug.examineGapFP, "no overlap found between frags %d and %d (lost after flap trimming)\n", lFragIid, rFragIid);
    restoreDefaultLocalAlignerVariables();
    return 0;
  }


  if (debug.examineGapLV > 0) {
    fprintf(debug.examineGapFP, "post-flap trimming ahang: %d, bhang:%d, length: %d, diffs: %d, diffs / length %%: %f\n",
            overlap->begpos, overlap->endpos, overlap->length, overlap->diffs,
            100.0 * overlap->diffs / overlap->length);
    PrintALNoverlap("", lcompBuffer, rcompBufferTrimmed, overlap);
    //debug.examineGapFP
  }

  if (debug.examineGapLV > 0) {
    fprintf(debug.examineGapFP, "found overlap between frags %d and %d, length = %d\n",
            lFragIid, rFragIid, overlap->length);
    fprintf(debug.examineGapFP, "ahang + overlap->length: %d, strlen(lcompBuffer): "F_SIZE_T", diff: "F_SIZE_T"\n",
            overlap->begpos + overlap->length,
            strlen(lcompBuffer),
            overlap->begpos + overlap->length - strlen(lcompBuffer));
    fprintf(debug.examineGapFP, "overlap->length + bhang: %d, strlen(rcompBufferTrimmed): "F_SIZE_T", diff: "F_SIZE_T"\n",
            overlap->length + overlap->endpos,
            strlen(rcompBufferTrimmed),
            overlap->length + overlap->endpos - strlen(rcompBufferTrimmed));
  }

  *lcontigBasesIntact = lcontigBaseStart;
  *ahang              = MAX (0, overlap->begpos);
  *olapLengthOut      = overlap->length;
  *bhang              = MAX (0, overlap->endpos);
  *currDiffs          = overlap->diffs;
  *rcontigBasesIntact = MAX((int) rcontig->bpLength.mean - CONTIG_BASES, 0);

  // calculate how many bases would be changed if gap was closed
  // take the whole lcompBuffer, subtract lcontigBasesUsed and lFragContigOverlapLength from ahang
  // and add in the length of the overlap
  // take the whole rcompBuffer, subtract rcontigBasesUsed and rFragContigOverlapLength from bhang

  int baseChangeLeftContig  = overlap->begpos - lcontigBasesUsed - lFragContigOverlapLength;
  int baseChangeRightContig = overlap->endpos - rcontigBasesUsed - rFragContigOverlapLength;

  int basesAdded = baseChangeLeftContig + overlap->length + baseChangeRightContig;

  gapSize = FindGapLength(lcontig, rcontig, FALSE);

  totalContigsBaseChange += basesAdded;

  *closedGapDelta = basesAdded - (int) gapSize.mean;

  if (debug.examineGapLV > 0) {
    fprintf(debug.examineGapFP, "lcontigBasesIntact: %d\n", *lcontigBasesIntact);
    fprintf(debug.examineGapFP, "overlap->begpos:    %d\n", overlap->begpos);
    fprintf(debug.examineGapFP, "overlap->length:    %d\n", overlap->length);
    fprintf(debug.examineGapFP, "overlap->endpos:    %d\n", overlap->endpos);
    fprintf(debug.examineGapFP, "rcontigBasesIntact: %d\n", *rcontigBasesIntact);

    fprintf(debug.examineGapFP, "lcontig->bpLength.mean: %f (%d change)\n", lcontig->bpLength.mean, baseChangeLeftContig);
    fprintf(debug.examineGapFP, "rcontig->bpLength.mean: %f (%d change)\n", rcontig->bpLength.mean, baseChangeRightContig);

    fprintf(debug.examineGapFP, "would fill gap %d of size %d with %d bases, net change: %d\n", gapNumber, (int) gapSize.mean, basesAdded, basesAdded - (int) gapSize.mean);

    fprintf(debug.examineGapFP, "new contig size:        %f\n", lcontig->bpLength.mean + baseChangeLeftContig + rcontig->bpLength.mean + baseChangeRightContig + overlap->length);
    fprintf(debug.examineGapFP, "totalContigsBaseChange: %d\n", totalContigsBaseChange);
  }

  restoreDefaultLocalAlignerVariables();

  return 1;
}
