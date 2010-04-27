

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

static const char *rcsid = "$Id: AS_ALN_loverlapper.c,v 1.24 2010-02-17 01:32:57 brianwalenz Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "AS_global.h"
#include "AS_UTL_reverseComplement.h"
#include "CA_ALN_local.h"
#include "AS_ALN_aligners.h"

#define OKNAFFINE 1

#undef DEBUG_LOCALOVL
#undef DEBUG_LOCALOVL_VERBOSE
#undef DEBUG_TRACE_BUILDING
#undef DEBUG_TRACE_CONSTR
#undef DEBUG_FINAL_TRACE

//#define DEBUG_LOCALOVL
//#define DEBUG_LOCALOVL_VERBOSE
//#define DEBUG_TRACE_CONSTR
//#define DEBUG_FINAL_TRACE

#undef PRINT_SUMMARY
#undef VERBOSE_SUMMARY


//maximum number of matching segments that can be pieced together
int MaxGaps= 3;

//maximum allowed mismatch at end of overlap
int MaxBegGap= 200;

//maximum allowed mismatch at end of overlap
int MaxEndGap= 200;

//biggest gap internal to overlap allowed
int MaxInteriorGap=400;

//whether to treat the beginning of the b fragment and
//   the end of the a fragment as allowed to have more error
int asymmetricEnds=0;

//amount of mismatch at end of the overlap that can cause
//   an overlap to be rejected
int MaxFreeFlap=20;

// set useSizeToOrderBlocks to 1 to get block mismatches resolved such that
//    the smaller block comes first:
//
//       .........AA------.........
//       .........--BBBBBB........
//
// or set it to 0 to get resolution with the A block always first
int useSizeToOrderBlocks = 1;

//global variable holding largest block mismatch of last returned overlap
int max_indel_AS_ALN_LOCOLAP_GLOBAL;


int ENDGAPHACK=3;
int MINCORESEG=30;


#undef AFFINE_QUALITY   /* overlap diff and length reported in affine terms */


/* safely copy a substring into memory pointed to by a static pointer */
static void safe_suffix(char **Dest,int *DestLen,char *src,int start){
  int len=strlen(src+start);
  if(*DestLen==0||*DestLen<len+1){
    *DestLen=len+1;
    if(*Dest==NULL){
      *Dest=(char*)safe_malloc(sizeof(char)*(*DestLen));
    } else {
      *Dest=(char*)safe_realloc(*Dest,sizeof(char)*(*DestLen));
    }
  }
  strcpy(*Dest,src+start);
  assert(strlen(*Dest)==len);
}

/* print alignment of a "piece" -- one local alignment in the overlap chain*/
static void print_piece(Local_Overlap *O,int piece,char *aseq,char *bseq){
  int alen,blen,segdiff,spnt,epnt;
  static char *aseg,*bseg;
  static int aseglen=0,bseglen=0, *segtrace;

  alen=O->chain[piece].piece.aepos-O->chain[piece].piece.abpos;
  blen=O->chain[piece].piece.bepos-O->chain[piece].piece.bbpos;

  /* make sure there is (persistant) space for the strings */
  if(aseglen<alen+1){
    aseglen=2*(alen+1);
    if(aseg==NULL){
      aseg=(char*)safe_malloc(sizeof(char)*aseglen);
    } else {
      aseg=(char*)safe_realloc(aseg,sizeof(char)*aseglen);
    }
  }
  if(bseglen<blen+1){
    bseglen=2*(blen+1);
    if(bseg==NULL){
      bseg=(char*)safe_malloc(sizeof(char)*bseglen);
    } else {
      bseg=(char*)safe_realloc(bseg,sizeof(char)*bseglen);
    }
  }

  /* copy the segments */
  strncpy(aseg,aseq+O->chain[piece].piece.abpos,alen);
  aseg[alen]='\0';
  strncpy(bseg,bseq+O->chain[piece].piece.bbpos,blen);
  bseg[blen]='\0';

  assert(strlen(bseg)==blen);
  assert(strlen(aseg)==alen);

  /* guesstimate the required number of diagonals/edits to consider to
     get optimal alignment */
  segdiff=1+(int)((O->chain[piece].piece.aepos
                   -O->chain[piece].piece.abpos)
                  *1.5*O->chain[piece].piece.error);


  /* get trace for the segment from AS_ALN_OKNAlign */
  spnt=0;
  /* subtract from aseg, bseg because Gene likes to index from 1, not 0 */
#ifdef OKNAFFINE
  epnt=0;
  segtrace=AS_ALN_OKNAffine(aseg-1,alen,bseg-1,blen,&spnt,&epnt,segdiff);
#else
  segtrace=AS_ALN_OKNAlign(aseg-1,alen,bseg-1,blen,&spnt,segdiff);
#endif

#ifdef DEBUG_LOCALOVL
  { int p=0,pos=0,neg=0;
    while(segtrace[p]!=0){
      fprintf(stderr, "trace[%d] = %d\n",p,segtrace[p]);
      if(segtrace[p]>0){
	pos++;
      } else {
	neg++;
      }
      p++;
    }
    fprintf(stderr, "alen %d blen %d pos %d neg %d spnt %d epnt %d\n",alen,blen,pos,neg,spnt,epnt);
  }
#endif

  if ((alen > 0) && (blen > 0)) {
    if(epnt!=0){
      if(epnt<0){
        aseg[alen+epnt]='\0';
      } else {
        bseg[blen-epnt]='\0';
      }
    }
    if(spnt<0){
      int i=0;
      fprintf(stderr, "(B sequence on top (2nd place in src); spnt=%d!!!\n",spnt);
      while(segtrace[i]!=0){segtrace[i++]*=-1;}
      PrintAlign(stderr,-spnt,0,bseg,aseg,segtrace);
      i=0; while(segtrace[i]!=0){segtrace[i++]*=-1;}
    } else {
      PrintAlign(stderr,spnt,0,aseg,bseg,segtrace);
    }
  }
}


/* Create a trace to be interpreted as with DP_Compare_AS, but based
   on a Local_Overlap record.  A Local_Segment within the overlap
   will be aligned using AS_ALN_OKNAlign(), generating a subtrace.  Subtraces,
   with their indices appropriately adjusted, will be spliced together
   by an encoding of the gaps between segments; for now, we'll simply insert
   gaps as follows:

   A "gap" with x bases in A and y bases in B will become a section
   of the alignment x+y positions long, with the A fragment first
   and the B fragment second (with '=' indicating an aligned match):

   ====AAAAAAAAAA--------------======
   ====----------BBBBBBBBBBBBBB======

   Obviously, a more compact treatment is possible, but this makes clear
   the presumptive blocks involved in the mismatch; also, Karin says that it
   will make Consensus happy.  One slight change

   Assumptions:  (do they matter?)
   By the usual conventions, the ahang should be nonnegative, the
   bhang negative only if the ahang is positive, and both sequences
   should be in the forward orientation.

*/



static
int
compressChain(Local_Overlap *O, int i) {
  int I, j, k;

  //  Compress out any deleted segments.  Assumes the last piece is
  //  never deleted.

  for (j=0, k=0; k<=O->num_pieces; j++, k++) {

    //  Remember the new location of piece i.
    if (i == k)
      I = j;

    //  If j is deleted, shift k into it, and delete k
    if (O->chain[j].agap==0 &&
        O->chain[j].bgap==0 &&
        O->chain[j].piece.abpos==O->chain[j].piece.aepos &&
        O->chain[j].piece.bbpos==O->chain[j].piece.bepos) {

      //  Advance k to the next non-deleted piece
      while (O->chain[k].agap==0 &&
             O->chain[k].bgap==0 &&
             O->chain[k].piece.abpos==O->chain[k].piece.aepos &&
             O->chain[k].piece.bbpos==O->chain[k].piece.bepos)
        k++;

      if (j != k) {
#ifdef DEBUG_LOCALOVL
        fprintf(stderr, "COMPRESS piece %d into empty piece %d\n", k, j);
#endif

        O->chain[j] = O->chain[k];
        O->chain[k].agap=0;
        O->chain[k].bgap=0;
        O->chain[k].piece.abpos=O->chain[k].piece.aepos;
        O->chain[k].piece.bbpos=O->chain[k].piece.bepos;
      }
    }
  }

  if (O->num_pieces != j-1) {
#ifdef DEBUG_LOCALOVL
    fprintf(stderr, "COMPRESS num pieces from %d to %d (i changes from %d to %d)\n", O->num_pieces, j-1, i, I);
#endif
    O->num_pieces = j-1;
  }

  return(I);
}



int *AS_Local_Trace(Local_Overlap *O, char *aseq, char *bseq){
  static int *TraceBuffer=NULL;
  int i,j,k,segdiff,*segtrace;
  static int allocatedspace=0;
  int tracespace=0;
  static char *aseg=NULL,*bseg=NULL;

  static int aseglen=0,bseglen=0;
  int abeg=0,bbeg=0; /* begining of segment; overloaded */
  int tracep=0; /* index into TraceBuffer */
  int spnt=0; /* to pass to AS_ALN_OKNAlign */
#ifdef OKNAFFINE
  int epnt=0; /* to pass to AS_ALN_OKNAffine */
#endif
  int alen,blen;
  int ahang = 0;

#ifdef DEBUG_LOCALOVL
  fprintf(stderr, "================================================================================\n");
  fprintf(stderr, "Start of Trace construction!\n\n\n");
#endif

  assert(O->num_pieces>0);

  {
    int n=O->num_pieces;
    O->chain[n].piece.abpos = O->chain[n].agap + O->chain[n-1].piece.aepos;
    O->chain[n].piece.bbpos = O->chain[n].bgap + O->chain[n-1].piece.bepos;
  }

  /*Estimate length required to store trace*/
  tracespace=0;
  tracespace+=abs(O->begpos)+abs(O->endpos);
  for(i=0;i<=O->num_pieces;i++){
    tracespace+=MAX(O->chain[i].agap,1);
    tracespace+=MAX(O->chain[i].bgap,1);
    tracespace+=(int)((O->chain[i].piece.aepos
		       -O->chain[i].piece.abpos)
		      *1.5*O->chain[i].piece.error);
    tracespace+=1000;
  }

  /*(Re)allocate space for the trace as necessary;
    Note that this is persistent storage so ...
    ... it doesn't need to get allocated on every call
    ... it shouldn't get freed
    ... results stored here need to be copied elsewhere if they
    are to be saved
  */
  if(allocatedspace<tracespace){
    allocatedspace=2*tracespace;
    if(TraceBuffer==NULL){
      TraceBuffer=(int*)safe_malloc(sizeof(int)*allocatedspace);
    } else {
      TraceBuffer=(int*)safe_realloc(TraceBuffer,sizeof(int)*allocatedspace);
    }
  }

  /* for each Local_Overlap chain[i].piece,
     need to handle the gap at the beginning and
     (for all but the final piece) the associated aligned segment */
  for(i=0;i<=O->num_pieces;i++){

#ifdef DEBUG_LOCALOVL
    fprintf(stderr, "Top of main loop in AS_Local_Trace, segment %d alignment:\n",i);

    if(i<O->num_pieces)
      print_piece(O,i,aseq,bseq);

    for(j=0;j<=O->num_pieces;j++)
      fprintf(stderr, "chain[%d] agap=%d abpos=%d aepos=%d  bgap=%d bbpos=%d bepos=%d\n",
              j,
              O->chain[j].agap, O->chain[j].piece.abpos, O->chain[j].piece.aepos,
              O->chain[j].bgap, O->chain[j].piece.bbpos, O->chain[j].piece.bepos);
#endif

    i = compressChain(O, i);

    /* Building an alignment/trace under the usual assumptions does not allow
       a given position in one sequence to simultaneously align to two or more
       positions in the other sequence.  However, the Find_Local_Overlap()
       routine can chain together local alignment segments that overlap.

       In order to make the local overlaps compatible with everything else
       we do, we need to trim back the overlaps.

       Since we will "output" this segment at the end of the loop,
       we need to fix its overlap with the following segment in this
       cycle through the loop
    */

    k=i+1;
    while(k<O->num_pieces){

#ifdef DEBUG_LOCALOVL
      {
        int j;
        fprintf(stderr,"----------------------------------------\n");
        fprintf(stderr,"Checking for overlaps of %d to %d\n",i,k);
        //fprintf(stderr, "Currently, segment %d alignment:\n",k);
        //print_piece(O,k,aseq,bseq);
        for(j=0;j<=O->num_pieces;j++)
          fprintf(stderr, "chain[%d] agap=%d abpos=%d aepos=%d  bgap=%d bbpos=%d bepos=%d\n",
                  j,
                  O->chain[j].agap, O->chain[j].piece.abpos, O->chain[j].piece.aepos,
                  O->chain[j].bgap, O->chain[j].piece.bbpos, O->chain[j].piece.bepos);
      }
#endif

      if(O->chain[k].piece.abpos<O->chain[i].piece.aepos||
         O->chain[k].piece.bbpos<O->chain[i].piece.bepos){

#ifdef DEBUG_LOCALOVL
        fprintf(stderr,"Need to fix overlaps of %d to %d\n",i,k);
#endif

	/* handle possibility of the first segment being
           contained within the second;originally simply asserted
           against this; now, try to handle by deleting first segment */

	if(O->chain[i].piece.abpos>O->chain[k].piece.abpos||
	   O->chain[i].piece.bbpos>O->chain[k].piece.bbpos){

#ifdef DEBUG_LOCALOVL
	  fprintf(stderr,"Warning: deleting contained piece of local_overlap (#1): piece %d contained in piece %d\n", i, k);
#endif

	  O->chain[i].agap=0;
	  O->chain[i].bgap=0;
          O->chain[i].piece.abpos=O->chain[i].piece.aepos;
          O->chain[i].piece.bbpos=O->chain[i].piece.bepos;

          if (i == 0) {
            O->chain[k].agap=O->chain[k].piece.abpos - 1;
            O->chain[k].bgap=O->chain[k].piece.bbpos - 1;
          } else {
            O->chain[k].agap=O->chain[k].piece.abpos-O->chain[i-1].piece.aepos;
            O->chain[k].bgap=O->chain[k].piece.bbpos-O->chain[i-1].piece.bepos;
          }
          //  Nothing left to do since we deleted the 'master' piece.
          k = O->num_pieces;
          continue;

        } else if(O->chain[i].piece.aepos>O->chain[k].piece.aepos||
                  O->chain[i].piece.bepos>O->chain[k].piece.bepos){

          /* otherwise, check for 2nd piece contained within first */

          /* if the next piece is completely within current piece,
             effectively remove it */

#ifdef DEBUG_LOCALOVL
          fprintf(stderr,"Warning: deleting contained piece of local_overlap (#2): piece %d contains piece %d\n", i, k);
#endif

          O->chain[k].agap=0;
          O->chain[k].bgap=0;
          O->chain[k].piece.abpos=O->chain[k].piece.aepos;
          O->chain[k].piece.bbpos=O->chain[k].piece.bepos;

          if(k+1<=O->num_pieces){
            // Deleted piece k, so need to update gaps using first
            // previous non-deleted and first next non-deleted
            // piece.
            int fp=k, fn=k;

            while(O->chain[fp].agap==0 &&
                  O->chain[fp].bgap==0 &&
                  O->chain[fp].piece.abpos==O->chain[fp].piece.aepos &&
                  O->chain[fp].piece.bbpos==O->chain[fp].piece.bepos){
              fp--;
              assert(fp>=0);
            }
            while(O->chain[fn].agap==0 &&
                  O->chain[fn].bgap==0 &&
                  O->chain[fn].piece.abpos==O->chain[fn].piece.aepos &&
                  O->chain[fn].piece.bbpos==O->chain[fn].piece.bepos){
              fn++;
              assert(fn<=O->num_pieces);
            }

            O->chain[fn].agap = O->chain[fn].piece.abpos - O->chain[fp].piece.aepos;
            O->chain[fn].bgap = O->chain[fn].piece.bbpos - O->chain[fp].piece.bepos;

#ifdef DEBUG_LOCALOVL
            fprintf(stderr,"Resetting gaps (#1) for segment %d to point to end of segment %d -- GAPS for %d now %d,%d\n", fn, fp, fn, O->chain[fn].agap, O->chain[fn].bgap);
#endif
          }

          /* else, fix the overlap */
        } else {

#ifdef DEBUG_LOCALOVL
          fprintf(stderr, "fix_overlapping_pieces() i=%d k=%d\n", i, k);
#endif
          fix_overlapping_pieces(aseq,
                                 bseq,
                                 O,i,k);

          // if the second piece disappeared
          if(O->chain[k].piece.abpos==O->chain[k].piece.aepos||
             O->chain[k].piece.bbpos==O->chain[k].piece.bepos){

#ifdef DEBUG_LOCALOVL
            fprintf(stderr,"piece %d segment size now 0!\n", k);
#endif
            //  DELETE
            O->chain[k].agap=0;
            O->chain[k].bgap=0;
            O->chain[k].piece.abpos=O->chain[k].piece.aepos;
            O->chain[k].piece.bbpos=O->chain[k].piece.bepos;

            if(k+1<=O->num_pieces){
              // Deleted piece k (which isn't the last piece), so need
              // to update gaps using first previous non-deleted and
              // first next non-deleted piece.
              int fp=k, fn=k;

              while(O->chain[fp].agap==0 &&
                    O->chain[fp].bgap==0 &&
                    O->chain[fp].piece.abpos==O->chain[fp].piece.aepos &&
                    O->chain[fp].piece.bbpos==O->chain[fp].piece.bepos){
                fp--;
                assert(fp>=0);
              }
              while(O->chain[fn].agap==0 &&
                    O->chain[fn].bgap==0 &&
                    O->chain[fn].piece.abpos==O->chain[fn].piece.aepos &&
                    O->chain[fn].piece.bbpos==O->chain[fn].piece.bepos){
                fn++;
                assert(fn<=O->num_pieces);
              }

              O->chain[fn].agap = O->chain[fn].piece.abpos - O->chain[fp].piece.aepos;
              O->chain[fn].bgap = O->chain[fn].piece.bbpos - O->chain[fp].piece.bepos;

#ifdef DEBUG_LOCALOVL
              fprintf(stderr,"Resetting gaps (#2) for segment %d to point to end of segment %d -- GAPS for %d now %d,%d\n", fn, fp, fn, O->chain[fn].agap, O->chain[fn].bgap);
#endif
            }
          } else if (O->chain[i].piece.abpos==O->chain[i].piece.aepos||
                     O->chain[i].piece.bbpos==O->chain[i].piece.bepos){
            // if the first piece disappeared

#ifdef DEBUG_LOCALOVL
            fprintf(stderr,"piece %d segment size now 0!\n", i);
#endif

            assert(i < k);

            //  DELETE segment i
            O->chain[i].agap=0;
            O->chain[i].bgap=0;
            O->chain[i].piece.abpos=O->chain[i].piece.aepos;
            O->chain[i].piece.bbpos=O->chain[i].piece.bepos;

            // Deleted piece i, so need to update gaps of pieces
            // surrounding i, using first previous non-deleted and
            // first next non-deleted piece.
            {
              int fp=i, fn=i;

              while(fp > 0 &&
                    O->chain[fp].agap==0 &&
                    O->chain[fp].bgap==0 &&
                    O->chain[fp].piece.abpos==O->chain[fp].piece.aepos &&
                    O->chain[fp].piece.bbpos==O->chain[fp].piece.bepos){
                fp--;
                assert(fp>=0);
              }
              while(O->chain[fn].agap==0 &&
                    O->chain[fn].bgap==0 &&
                    O->chain[fn].piece.abpos==O->chain[fn].piece.aepos &&
                    O->chain[fn].piece.bbpos==O->chain[fn].piece.bepos){
                fn++;
                assert(fn<=O->num_pieces);
              }

              if (i == 0) {
                //  Deleted first piece.
                O->chain[fn].agap=O->chain[fn].piece.abpos - 1;
                O->chain[fn].bgap=O->chain[fn].piece.bbpos - 1;
              } else {
                O->chain[fn].agap=O->chain[fn].piece.abpos-O->chain[fp].piece.aepos;
                O->chain[fn].bgap=O->chain[fn].piece.bbpos-O->chain[fp].piece.bepos;
              }
            }
          }
        }  //  fix the overlap
      }  //  need to fix

      k++;
    }  //  over all pieces

#ifdef DEBUG_LOCALOVL
    fprintf(stderr, "After fixing, segment %d alignment:\n",i);

    if(i<O->num_pieces)
      print_piece(O,i,aseq,bseq);

    for(j=0;j<=O->num_pieces;j++)
      fprintf(stderr, "chain[%d] agap=%d abpos=%d aepos=%d  bgap=%d bbpos=%d bepos=%d\n",
              j,
              O->chain[j].agap, O->chain[j].piece.abpos, O->chain[j].piece.aepos,
              O->chain[j].bgap, O->chain[j].piece.bbpos, O->chain[j].piece.bepos);
#endif


    /* if conditions indicate the segment was deleted previously, skip! */
    if(O->chain[i].agap==0 &&
       O->chain[i].bgap==0 &&
       O->chain[i].piece.abpos==O->chain[i].piece.aepos &&
       O->chain[i].piece.bbpos==O->chain[i].piece.bepos) {
#ifdef DEBUG_LOCALOVL
      fprintf(stderr, "deleted piece %d, the one we were working on.\n", i);
#endif
      i--;  //  To do this piece again, once we compress.
      continue;
    }

    //
    //  COMPRESS AGAIN
    //
    i = compressChain(O, i);

    //  If we deleted the first piece, the begpos changes
    //O->begpos = O->chain[0].piece.abpos - O->chain[0].piece.bbpos;;

#ifdef DEBUG_LOCALOVL
    fprintf(stderr,"----------------------------------------\n");
    fprintf(stderr,"Processing piece %d out of %d\n",i, O->num_pieces);
#endif

    /* set up positions before which gaps are inserted to handle
       the gap portion of a chain piece */

    /* put gaps before beginning of aligned piece (but after the portion
       of aseq in the gap); location is relative to the beginning of the
       alignment (i.e., ignores ahang worth of positions)  */

    if (i < O->num_pieces)
      abeg=O->chain[i].piece.abpos;
    else
      abeg=O->chain[i-1].piece.aepos + O->chain[i].agap;

#ifdef DEBUG_LOCALOVL
    fprintf(stderr, "abeg=%d\n", abeg);
    fprintf(stderr, "endpos=%d agap=%d bgap=%d\n", O->endpos, O->chain[i].agap, O->chain[i].bgap);
#endif

#ifdef DEBUG_LOCALOVL
    fprintf(stderr, "abeg=%d\n", abeg);
    fprintf(stderr, "endpos=%d agap=%d bgap=%d\n", O->endpos, O->chain[i].agap, O->chain[i].bgap);
#endif

    /*handle boundary case to prevent gaps preceding the b sequence*/
    if((i==0) && (O->chain[i].bgap>0)) {
      assert(O->chain[i].agap>=0);

      if(O->begpos>=0){
        O->begpos = O->chain[i].piece.abpos - 1;
        assert(O->begpos>=0);
        O->chain[i].agap=0;

        if (O->chain[i].bgap != O->chain[i].piece.bbpos - 1) {
          fprintf(stderr, "ERROR!  bgap != bbpos-1\n");
          for(j=0;j<=O->num_pieces;j++)
            fprintf(stderr, "chain[%d] agap=%d abpos=%d aepos=%d  bgap=%d bbpos=%d bepos=%d\n",
                    j,
                    O->chain[j].agap, O->chain[j].piece.abpos, O->chain[j].piece.aepos,
                    O->chain[j].bgap, O->chain[j].piece.bbpos, O->chain[j].piece.bepos);
        }

        assert(O->chain[i].bgap == O->chain[i].piece.bbpos - 1);
      } else {
        O->begpos-=O->chain[i].bgap;
        O->chain[i].bgap=0;
      }
#ifdef DEBUG_LOCALOVL
      fprintf(stderr, "begpos=%d agap=%d bgap=%d\n", O->begpos, O->chain[i].agap, O->chain[i].bgap);
#endif
    }

    /* now prevent gaps at end of A sequence */
    //    if(i==O->num_pieces&&O->endpos>=0){
    //      O->endpos+=O->chain[i].bgap;
    //      O->chain[i].bgap=0;
    //    }

    /* now make sure that end mismatches are treated by tucking
       the shorter tail into a gap before the longer tail,
       or, ifdef FORCEPOSITIVEBHANG, by tucking the A tail into
       a gap before the B tail */

    if(i==O->num_pieces){
      if(O->endpos>=0){
	O->endpos+=O->chain[i].bgap;
	O->chain[i].bgap=0;
      }else {
	O->endpos-=O->chain[i].agap;
        abeg-=O->chain[i].agap;
        O->chain[i].agap=0;
      }
#ifdef DEBUG_LOCALOVL
      fprintf(stderr, "abeg=%d\n", abeg);
      fprintf(stderr, "endpos=%d agap=%d bgap=%d\n", O->endpos, O->chain[i].agap, O->chain[i].bgap);
#endif
    }

    /* put gaps before the portion of bseq in the gap; for the first
       piece, this means before position 0 */

    if (i==0) {
      ahang=O->begpos;
      bbeg=1-MIN(ahang,0);
    } else {
      bbeg=O->chain[i-1].piece.bepos;
    }

    /* now insert the right number of gaps! */

    ///////////////////////////////////////

    if(i==O->num_pieces){

      if(O->endpos<0){
	O->chain[i].agap+=-O->endpos;
	O->endpos=0;
      } else {
	O->chain[i].bgap+=O->endpos;
	O->endpos=0;
      }

      if(O->chain[i].agap <= O->chain[i].bgap){
	O->endpos=O->chain[i].bgap;
	O->chain[i].bgap=0;
      }else{
	O->endpos=-O->chain[i].agap;
	O->chain[i].agap=0;
      }

#ifdef DEBUG_LOCALOVL
      fprintf(stderr, "endpos=%d agap=%d bgap=%d\n", O->endpos, O->chain[i].agap, O->chain[i].bgap);
#endif
    }

    if(O->chain[i].agap <= O->chain[i].bgap || ! useSizeToOrderBlocks ){

      /* start by putting len(agap) gaps before the chunk of B in the gap */

      for(j=0;j<O->chain[i].agap;j++){
	TraceBuffer[tracep++]=bbeg;
#ifdef DEBUG_TRACE_BUILDING
        fprintf(stderr, "TraceBuffer1 %d = %d\n", tracep-1, TraceBuffer[tracep-1]);
#endif
      }


      /* then put len(bgap) gaps before the chunk of A in the gap */

      for(j=0;j<O->chain[i].bgap;j++){
	TraceBuffer[tracep++]=-abeg;
#ifdef DEBUG_TRACE_BUILDING
        fprintf(stderr, "TraceBuffer2 %d = %d\n", tracep-1, TraceBuffer[tracep-1]);
#endif
      }


    } else { // if the bgap is smaller,

      abeg-=O->chain[i].agap;
      bbeg+=O->chain[i].bgap;

      //      fprintf(stderr,"SWAPPING BLOCK MISMATCH GAPS\n");

      /* start by putting len(bgap) gaps before the chunk of A in the gap */

      for(j=0;j<O->chain[i].bgap   ;j++){
	TraceBuffer[tracep++]=-abeg;
#ifdef DEBUG_TRACE_BUILDING
        fprintf(stderr, "TraceBuffer3 %d = %d\n", tracep-1, TraceBuffer[tracep-1]);
#endif
      }

      /* then put len(agap) gaps before the chunk of B in the gap */

      for(j=0;j<O->chain[i].agap ;j++){
	TraceBuffer[tracep++]=bbeg;
#ifdef DEBUG_TRACE_BUILDING
        fprintf(stderr, "TraceBuffer4 %d = %d\n", tracep-1, TraceBuffer[tracep-1]);
#endif
      }

    }

    ///////////////////////////////////////

    /* if last piece, there is no aligned segment */

    if(i==O->num_pieces)break;


    /* set bbeg to beginning of aligned segment for piece */

    abeg=O->chain[i].piece.abpos ;
    bbeg=O->chain[i].piece.bbpos /* -MIN(ahang,0) ??*/;


    /* set lengths of segments */

    alen=O->chain[i].piece.aepos-abeg; /* check +1?? */
    blen=O->chain[i].piece.bepos-bbeg; /* check +1?? */


    /* create strings for just the parts of the sequences in the
       aligned segment */


    /* make sure there is (persistant) space for the strings */
    if(aseglen<alen+1){
      aseglen=2*(alen+1);
      if(aseg==NULL){
	aseg=(char*)safe_malloc(sizeof(char)*aseglen);
      } else {
	aseg=(char*)safe_realloc(aseg,sizeof(char)*aseglen);
      }
    }
    if(bseglen<blen+1){
      bseglen=2*(blen+1);
      if(bseg==NULL){
	bseg=(char*)safe_malloc(sizeof(char)*bseglen);
      } else {
	bseg=(char*)safe_realloc(bseg,sizeof(char)*bseglen);
      }
    }


    /* copy the segments */
    strncpy(aseg,aseq+abeg,alen);
    aseg[alen]='\0';
    strncpy(bseg,bseq+bbeg,blen);
    bseg[blen]='\0';

    assert(strlen(bseg)==blen);
    assert(strlen(aseg)==alen);

    /* guesstimate the required number of diagonals/edits to consider to
       get optimal alignment */
    segdiff=1+(int)((O->chain[i].piece.aepos
		     -O->chain[i].piece.abpos)
		    *1.5*O->chain[i].piece.error);


    /* get trace for the segment from AS_ALN_OKNAlign */
    spnt=0;
    /* subtract from aseg, bseg because Gene likes to index from 1, not 0 */
#ifdef OKNAFFINE
    epnt=0;
    segtrace=AS_ALN_OKNAffine(aseg-1,alen,bseg-1,blen,&spnt,&epnt,segdiff);

    if(epnt!=0){

      //  Find the next non deleted segment
      int ll = i+1;

      while (O->chain[ll].agap==0 &&
             O->chain[ll].bgap==0 &&
             O->chain[ll].piece.abpos == O->chain[ll].piece.aepos &&
             O->chain[ll].piece.bbpos == O->chain[ll].piece.bepos)
        ll++;

      assert(ll <= O->num_pieces);

      if(epnt>0){ /* throwing away some of B segment */
        O->chain[ll].bgap+=epnt;
        O->chain[i].piece.bepos-=epnt;
        assert(O->chain[i].piece.bbpos<=O->chain[i].piece.bepos);
      } else {
	O->chain[ll].agap-=epnt;
	O->chain[i].piece.aepos+=epnt;
	assert(O->chain[i].piece.abpos<=O->chain[i].piece.aepos);
      }
    }
#else
    segtrace=AS_ALN_OKNAlign(aseg-1,alen,bseg-1,blen,&spnt,segdiff);
#endif
    if(spnt>0){
      O->chain[i].agap+=spnt;
      O->chain[i].piece.abpos+=spnt;
    } else {
      O->chain[i].bgap-=spnt;
      O->chain[i].piece.bbpos-=spnt;
    }


#ifdef DEBUG_LOCALOVL_VERBOSE
    if(spnt<0){
      int i=0;
      fprintf(stderr, "(B sequence on top (2nd place in src); spnt=%d!!!\n",spnt);
      while(segtrace[i]!=0){segtrace[i++]*=-1;}
      PrintAlign(stderr,-spnt,0,bseg,aseg,segtrace);
      i=0; while(segtrace[i]!=0){segtrace[i++]*=-1;}
    } else {
      PrintAlign(stderr,spnt,0,aseg,bseg,segtrace);
    }
#endif

    assert(segtrace!=NULL);


#ifdef DEBUG_LOCALOVL_VERBOSE
    fprintf(stderr,"  Match A[%d,%d] to B[%d,%d]\n",
	    O->chain[i].piece.abpos,
	    O->chain[i].piece.aepos,
	    O->chain[i].piece.bbpos,
	    O->chain[i].piece.bepos);
#endif


    /* Now copy the segment trace into master trace, adjusting positions */
    j=0;

    //  Gaps for what we trimmed off
    if(spnt<0){
      int ctr;
      for(ctr=0;ctr<abs(spnt);ctr++){
	TraceBuffer[tracep++]=-abeg;
#ifdef DEBUG_TRACE_BUILDING
        fprintf(stderr, "TraceBuffer5 %d = %d\n", tracep-1, TraceBuffer[tracep-1]);
#endif
      }
    } else {
      int ctr;
      for(ctr=0;ctr<spnt;ctr++){
	TraceBuffer[tracep++]=bbeg;
#ifdef DEBUG_TRACE_BUILDING
        fprintf(stderr, "TraceBuffer6 %d = %d\n", tracep-1, TraceBuffer[tracep-1]);
#endif
      }
    }

    //  Then the actual segment

    while(segtrace[j]!=0){
      //fprintf(stderr," %d",segtrace[j]+(segtrace[j]<0 ? 1 +MAX(0,spnt) : -1-MAX(0,-spnt)));
      if(segtrace[j]<0){
	TraceBuffer[tracep++]=-abeg+segtrace[j++]+1 /* -MAX(0,spnt) ?? */;
#ifdef DEBUG_TRACE_BUILDING
        fprintf(stderr, "TraceBuffer7 %d = %d\n", tracep-1, TraceBuffer[tracep-1]);
#endif
      } else {
	TraceBuffer[tracep++]=bbeg+segtrace[j++]-1 /* +MAX(0,-spnt) ?? */;
#ifdef DEBUG_TRACE_BUILDING
        fprintf(stderr, "TraceBuffer8 %d = %d\n", tracep-1, TraceBuffer[tracep-1]);
#endif
      }
    }

#ifdef DEBUG_LOCALOVL_VERBOSE
    fprintf(stderr,"\n");
    if(spnt<0){
      int i=0;
      fprintf(stderr, "(B sequence on top (2nd place in src); spnt=%d!!!\n",spnt);
      while(segtrace[i]!=0){segtrace[i++]*=-1;}
      PrintAlign(stderr,-spnt,0,bseg,aseg,segtrace);
      i=0; while(segtrace[i]!=0){segtrace[i++]*=-1;}
    } else {
      PrintAlign(stderr,spnt,0,aseg,bseg,segtrace);
    }
#endif


    /* and back to the top of the loop for another overlap piece */

#ifdef DEBUG_LOCALOVL
    fprintf(stderr, "Done with chain loop for piece %d\n", i);
    for(j=0;j<=O->num_pieces;j++)
      fprintf(stderr, "chain[%d] agap=%d abpos=%d aepos=%d  bgap=%d bbpos=%d bepos=%d\n",
              j,
              O->chain[j].agap, O->chain[j].piece.abpos, O->chain[j].piece.aepos,
              O->chain[j].bgap, O->chain[j].piece.bbpos, O->chain[j].piece.bepos);
#endif
  }

  /* terminate the trace */
  TraceBuffer[tracep]=0;

  /* really non-robust coding on allocated(trace)space, so let's at least test
     for a violation! */
  assert(tracep<allocatedspace);


  //  IsInvalidAlign() attempts to catch a long-standing bug in the
  //  above code where the resulting alignment trace is invalid.  If
  //  decoded, the alignment would extend past the end of one of the
  //  fragments, and usually has a gigantic 'hole' in the middle.
  //
  //  Note that if it is invalid, PrintAlign() should fail an assert.

#ifdef DEBUG_TRACE_CONSTR
  fprintf(stderr, "End of Trace HERE IT IS!\n");

    for(j=0;j<=O->num_pieces;j++)
      fprintf(stderr, "chain[%d] agap=%d abpos=%d aepos=%d  bgap=%d bbpos=%d bepos=%d\n",
              j,
              O->chain[j].agap, O->chain[j].piece.abpos, O->chain[j].piece.aepos,
              O->chain[j].bgap, O->chain[j].piece.bbpos, O->chain[j].piece.bepos);
#endif

  int  invalidTrace = 0;

  if(O->begpos>=0){
    invalidTrace = IsInvalidAlign(O->begpos,O->endpos,aseq+1,bseq+1,TraceBuffer);
#ifdef DEBUG_TRACE_CONSTR
    if (!invalidTrace)
      PrintAlign(stderr,O->begpos,O->endpos,aseq+1,bseq+1,TraceBuffer);
#endif

  }else{
    for(i=0;i<tracep;i++)
      TraceBuffer[i]*=-1;

    invalidTrace = IsInvalidAlign(-O->begpos,-O->endpos,bseq+1,aseq+1,TraceBuffer);
#ifdef DEBUG_TRACE_CONSTR
    if (!invalidTrace)
      PrintAlign(stderr,-O->begpos,-O->endpos,bseq+1,aseq+1,TraceBuffer);
#endif

    for(i=0;i<tracep;i++)
      TraceBuffer[i]*=-1;
  }

#ifdef DEBUG_TRACE_CONSTR
  fprintf(stderr, "End of Trace construction in AS_Local_Trace()!\n");
  fprintf(stderr, "================================================================================\n\n\n");
#endif

  return((invalidTrace) ? NULL : TraceBuffer);
}


/* Given fragments a and b, find the best overlap between them using
   local alignment code.  The purpose is to permit "bubbles" in the fragment
   graph due to multi-base polymorphisms to be smoothed out, such that a
   unitig can be constructed.

   The method relies on Myers' local overlap code.  MORE DETAILS????

   The local overlap code returns the best overlap available as a chain of
   matching segments.

   The current function is concerned with three things:

   1) Providing a wrapper around the local alignment code so that users
   can carry out operations on the level of the "message" (IFM, OVL, etc).

   2) Determine whether the best possible overlap is good enough to constitute
   an overlap for the purposes of unitigging.  A preliminary set of criteria
   are given below.

   3) Combine bits of existing code to provide a full alignment specification
   for the local overlapper; one can debate whether large gaps should be
   encoded in an edit trace, but for the purposes of integration with the
   existing code base, the current function will aim to specify an
   alignment in terms of an ahang, a bhang and a trace, which will get
   converted into the core elements of an overlap message just as in
   DP_Compare_AS.

   PRELIMINARY DECISION CRITERIA:

   - if an overall %-match of less than erate with length MinLen[=40]
   exists, then the overlap should be allowed regardless of the
   following additional criteria
   - at least one segment of at least MinLen bases
   - if we treat an affine event as a single error (and the length involved
   to likewise be a single base) for the purpose of calculating an
   error rate, the error rate should be no more than erate.
   - no more than MaxGaps [=2?] unaligned regions ("gaps") in the overlap;
   this also entails no more than MaxGaps+1 matching segments

   Additional criteria that we might want to use:

   - no more than MaxEndGap [=20?] unaligned bases at either end; this
   is to prevent a classic branchpoint out of a repeat into distinct
   regions from being overlapped; on the other hand, the requirement
   that all fragments in a bubble go back together might mean that
   this filter is unnecessary.
   - only one end may end in a gap; having both ends
   be in polymorphic regions should be quite unlikely, and it is
   not currently our concern to handle crappy fragment ends (= spurs,
   handled by other means).


   AS_FIND_LOCAL_OVERLAP:
   Just find a good alignment to the boundary of the d.p. matrix.  From
   this extrapolate a rough overlap relationship without an alignment
   and return the result (if there is one).
   AS_FIND_LOCAL_ALIGN:
   For this option, further go to the trouble of computing the alignment
   and store it in the overlap message.
   AS_FIND_LOCAL_ALIGN_NO_TRACE
   A hack by SAK...don't compute delta encoding.

   As with DP_Compare_AS, the resulting overlap message occupies
   memory which is persistently owned (reused) by the function; it
   must be copied if it is to be retained beyond the given call.

*/


OverlapMesg *
Local_Overlap_AS(InternalFragMesg *a, InternalFragMesg *b,
                 int beg, int end,
                 int opposite,
                 double erate, double thresh, int minlen,
                 CompareOptions what, int *where) {
  int   alen,  blen;
  char *aseq, *bseq;
  int NumSegs=0;
  int seglen,sumseglen=0;
  int noCoreSeg=1;
  int i;
  int forcenext=0;

  assert((0.0 <= erate) && (erate <= 4 * AS_MAX_ERROR_RATE));

  static char *Ausable=NULL, *Busable=NULL;
  static int AuseLen=0, BuseLen=0;

  int coreseglen=MIN(MINCORESEG,minlen);

  double avgerror=0.;

  static OverlapMesg QVBuffer;

  Local_Segment *local_results=NULL;
  Local_Overlap *O=NULL;

  int *trace;

  assert(what == AS_FIND_LOCAL_ALIGN || what == AS_FIND_LOCAL_ALIGN_NO_TRACE
	 || what == AS_FIND_LOCAL_OVERLAP );

  aseq = a->sequence;  /* Setup sequence access */
  bseq = b->sequence;
  alen = strlen(aseq);
  blen = strlen(bseq);

  if (opposite)
    reverseComplement(b->sequence, b->quality, blen);


  /* now generate the substrings of the input sequences which contain the
     portion usable according to the user's begin/end parameters */
#define EXPAND_BAND 10 /* amount to pad ends of  user's begin..end range */
#define STRETCH 1.5  /* amount of "stretching" allowed in bounding
			overlap end */
#define BIGPAD 500 /* amount of extra slop allowed in bounding overlap end */


  // not really using the "suffix" part of next calls; just a safe string copy
  safe_suffix(&Ausable,&AuseLen,aseq,0);
  safe_suffix(&Busable,&BuseLen,bseq,0);

  for(i=0;i<MAX(0,beg-EXPAND_BAND);i++){
    Ausable[i]='N';
  }
  assert(alen>=i);
  for(i=0;i<MAX(0,-end-EXPAND_BAND);i++){
    Busable[i]='N';
  }
  assert(blen>=i);

  for(i=MAX(0,end)+(blen+MIN(0,end))*STRETCH+BIGPAD;i<alen;i++){
    Ausable[i]='N';
  }

  i=-MIN(0,beg)+(alen-MAX(0,beg))*STRETCH+BIGPAD;
  assert(i>=0);

  for(;i<blen;i++){
    Busable[i]='N';
  }

#ifdef DEBUG_LOCALOVL
  fprintf(stderr, "Local_Overlap_AS()--   BEGIN beg=%d end=%d opposite=%d\n", beg, end, opposite);
  fprintf(stderr, "a = %s\n", Ausable);
  fprintf(stderr, "b = %s\n", Busable);
#endif

  /* Notes on handling of reverse complement overlap searching:
     1) The implementation here uses Find_Local_Segments and
     Find_Local_Overlaps only in the forward direction; the interaction
     of these two routines on reverse orientation searches is non-obvious
     [even Gene agrees], so we avoid it.
     2) Thus, it is the responsibility of this routine to reverse the B
     sequence prior to running these routines.
  */



  local_results=Find_Local_Segments(Ausable,strlen(Ausable),
				    Busable,strlen(Busable),
				    LOCAL_FORW, 8, erate*2., &NumSegs);


#ifdef DEBUG_LOCALOVL
  fprintf(stderr,"Find_Local_Segments found %d.\n",NumSegs);
#ifdef DEBUG_LOCALOVL
  { int i;
    for(i=0;i<NumSegs;i++){
      fprintf(stderr, "Segment %d: (%d,%d)[----------](%d,%d), error %e\n",i,
              local_results[i].abpos,
              local_results[i].bbpos,
              local_results[i].aepos,
              local_results[i].bepos,
              local_results[i].error);
    }
  }
#endif
#endif

  if(NumSegs==0){
#ifdef DEBUG_LOCALOVL
    fprintf(stderr,"No segments found by Find_Local_Segments()\n");
#endif
    goto nooverlap;
  }

  O=Find_Local_Overlap(alen,blen,
		       0 /*comp==0 -> fwd orientation */,
		       0 /*nextbest==0 -> best overlap*/,
		       local_results,NumSegs,
		       minlen-2*ENDGAPHACK
                       ,1.);



#ifdef DEBUG_LOCALOVL
  if(O!=NULL){
    for(i=0;i<=O->num_pieces;i++){
      fprintf(stderr, "overlap segment %d: (%d,%d)[---------](%d,%d)\n",
              i,
              O->chain[i].piece.abpos,
              O->chain[i].piece.bbpos,
              O->chain[i].piece.aepos,
              O->chain[i].piece.bepos
              );
    }
  }
#endif




  if(O!=NULL){
  next_overlap:

#ifdef DEBUG_LOCALOVL
    fprintf(stderr, "\n\nNEXT OVERLAP:\n");
#endif

    if(forcenext){
      safe_free(O);
      forcenext=0;
      //      fprintf(stderr, "Looking at lower-scoring overlaps\n");
      O=Find_Local_Overlap(alen,blen,
                           0 /*comp==0 -> fwd orientation */,
                           1 /*nextbest==1 -> next best */,
                           local_results,NumSegs,
                           minlen-2*ENDGAPHACK
                           ,1.);
      if(O==NULL)goto nooverlap;

#ifdef DEBUG_LOCALOVL
      for(i=0;i<=O->num_pieces;i++){
        fprintf(stderr, "overlap segment %d: (%d,%d)[---------](%d,%d)\n",
                i,
                O->chain[i].piece.abpos,
                O->chain[i].piece.bbpos,
                O->chain[i].piece.aepos,
                O->chain[i].piece.bepos
                );
      }
#endif

    }
  }

#ifdef DEBUG_LOCALOVL
  if(O!=NULL){
    fprintf(stderr,"Overlap found.\n");
    Print_Local_Overlap(stderr,O,5);
    //Print_Local_Overlap_withAlign(stderr,O,Ausable,Busable);  //  WARNING: bogus traces kill this
  }
#endif


  if(O==NULL){
#ifdef DEBUG_LOCALOVL
    fprintf(stderr,"No overlap out of Find_Local_Overlap()\n");
#endif
    goto nooverlap;
  }


  if(O->diffs/(double)O->length >= erate){
    /* global version fails */

    if(asymmetricEnds){

      /* Test for Mike's clear-range-extended overlaps */

      /* For the beginning of the overlap, the overlap should be rejected if:
         - the first N bases of the B sequence are not aligned to A
         - N > MaxBegGap (which should be set to the number of bases
         the clear range was extended, plus a bit of slop)
         - there is a "flap" (unaligned bases) at the beginning of
         the alignment of sufficient size; this distinguishes cases
         which simply have a negative ahang from cases where there is
         a problematic mismatch */

      if(O->chain[0].bgap>MaxFreeFlap&&
	 O->chain[0].bgap+MAX(0,-O->begpos)>MaxBegGap){
#ifdef DEBUG_LOCALOVL
	fprintf(stderr,"OVERLAP REJECTED because: begin flap too big %d\n\n",
                O->chain[0].bgap+MAX(0,-O->begpos)>MaxBegGap);
#endif
	forcenext=1;
	goto next_overlap;
      }

      /* For the end of the overlap, as above but reverse orientation */

      if(O->chain[O->num_pieces].agap>MaxFreeFlap&&
	 O->chain[O->num_pieces].agap+MAX(0,-O->endpos)>MaxEndGap){

#ifdef DEBUG_LOCALOVL
	fprintf(stderr,"OVERLAP REJECTED because: end flap too big %d\n\n",
		O->chain[O->num_pieces].agap+MAX(0,-O->endpos)>MaxEndGap);
#endif

	forcenext=1;
	goto next_overlap;

      }


    } else {  //  asymmetricEnds == 0


      if(MIN(O->chain[0].agap,
	     O->chain[0].bgap)>MaxBegGap){
#ifdef DEBUG_LOCALOVL
	fprintf(stderr,"OVERLAP REJECTED because: begin gap too big %d\n\n",
		MIN(O->chain[0].agap,
		    O->chain[0].bgap));
#endif
	forcenext=1;
	goto next_overlap;
      }

      if(MIN(O->chain[O->num_pieces].agap,
	     O->chain[O->num_pieces].bgap)>MaxEndGap){
#ifdef DEBUG_LOCALOVL
	fprintf(stderr,"OVERLAP REJECTED because: end gap too big %d\n\n",
		MIN(O->chain[O->num_pieces].agap,
		    O->chain[O->num_pieces].bgap));
#endif

	forcenext=1;
	goto next_overlap;

      }


    }  //  if asymmetricEnds
  }

  sumseglen=0;
  avgerror=0.;


  for(i=0;i<O->num_pieces;i++){

    if(i!=0&&
       MAX(O->chain[i].agap,
	   O->chain[i].bgap)>MaxInteriorGap){
#ifdef DEBUG_LOCALOVL
      fprintf(stderr,"OVERLAP REJECTED because: interior gap %d too big %d\n\n",
	      i,
	      MAX(O->chain[i].agap,
		  O->chain[i].bgap));
#endif

      forcenext=1;
      goto next_overlap;

    }

    /* segments assumed to be plus-sense wrt A sequence */
    assert(O->chain[i].piece.aepos-
	   O->chain[i].piece.abpos>0);

    seglen=1+((O->chain[i].piece.aepos-
               O->chain[i].piece.abpos
               +abs(O->chain[i].piece.bepos-
                    O->chain[i].piece.bbpos))
              /2);

    if(seglen>=coreseglen)noCoreSeg=0;
    sumseglen+=seglen;
    avgerror+=O->chain[i].piece.error*seglen;
  }

  //Test for at least one segment bigger than minlen
  if(noCoreSeg){

    /* if none found so far, might need to allow for a segment with short
       end gaps due to a mismatch (given +1/-3 scoring, 1 mismatch kills
       2 matches, for example) */
    int last=O->num_pieces-1;
    int begseglen = 1+
      (O->chain[i].piece.aepos-
       O->chain[i].piece.abpos
       +abs(O->chain[i].piece.bepos-
	    O->chain[i].piece.bbpos)
       )/2
      +MIN(ENDGAPHACK,O->chain[0].agap);
    int endseglen = 1+
      (O->chain[last].piece.aepos-
       O->chain[last].piece.abpos
       +abs(O->chain[last].piece.bepos-
	    O->chain[last].piece.bbpos)
       )/2
      +MIN(ENDGAPHACK,O->chain[last+1].agap);
    if(last==0)begseglen+=MIN(ENDGAPHACK,O->chain[1].agap);
    if(last==0)endseglen+=MIN(ENDGAPHACK,O->chain[0].agap);

    /* if even allowing some slop at ends we can't find a big enough segment,
       quit! */
    if(MAX(begseglen,endseglen)<minlen){
#ifdef DEBUG_LOCALOVL
      fprintf(stderr,"OVERLAP REJECTED because no anchoring segment (min = %d)\n\n",
	      minlen);
#endif

      forcenext=1;
      goto next_overlap;
    }
  }


  QVBuffer.quality=avgerror;

  /* Now that you have finish point, first compute complete alignment
     by calling AS_ALN_OKNAlign or extrapolate start point as dictated by
     `what'.  Then build an overlap record modelling the overlap.   */

  trace = NULL;

  {
    int ahang, bhang;
    int i, j;

#ifdef DEBUG_LOCALOVL
    fprintf(stderr, "before AS_Local_Trace:\n");
#endif

    // coordinates from Find_Local routines will be one off from
    // those expected by the trace routines, so adjust them!
    //	O->begpos+= (O->begpos>=0 ? 1 : -1);
    //	O->endpos+= (O->endpos>=0 ? 1 : -1);
    for(i=0;i<O->num_pieces;i++){
      O->chain[i].piece.abpos++;
      O->chain[i].piece.bbpos++;
      O->chain[i].piece.aepos++;
      O->chain[i].piece.bepos++;
#ifdef DEBUG_LOCALOVL
      fprintf(stderr, "chain[%d] agap=%d abpos=%d aepos=%d bgap=%d bbpos=%d bepos=%d\n",
              i,
              O->chain[i].agap, O->chain[i].piece.abpos, O->chain[i].piece.aepos,
              O->chain[i].bgap, O->chain[i].piece.bbpos, O->chain[i].piece.bepos);
#endif
    }

    //AS_Local_Trace assumes string pointer one before start of string!
    trace = AS_Local_Trace(O,Ausable-1,Busable-1);

    //  Bail if the trace is messed up.  This catches a long-standing
    //  bug that BPW can't find.
    //
    if (trace == NULL) {
#ifdef DEBUG_LOCALOVL
      fprintf(stderr,"OVERLAP REJECTED because: end invalid trace\n\n");
#endif
      forcenext=1;
      goto next_overlap;
    }


    //TEST NUMBER OF GAPS AT SEGMENT BOUNDARIES
    {
      int final_num_pieces=0;
      for(i=0;i<O->num_pieces;i++){
        if(O->chain[i].piece.abpos < O->chain[i].piece.aepos) final_num_pieces++;
      }
      if(O->chain[0].bgap==0){final_num_pieces--;}
      if(O->chain[O->num_pieces].agap>0){final_num_pieces++;}
      if(final_num_pieces> MaxGaps){
#ifdef DEBUG_LOCALOVL
        fprintf(stderr,"OVERLAP REJECTED because: too many pieces %d\n\n",
                O->num_pieces-
                (O->chain[0].type==LOCAL_BOUNDARY||
                 O->chain[0].agap+
                 O->chain[0].bgap==2)+
                (O->chain[O->num_pieces].type!=LOCAL_BOUNDARY) );
#endif

        forcenext=1;
        goto next_overlap;
      }
    }

#ifdef DEBUG_LOCALOVL
    fprintf(stderr,"========================================Overlap after trimming BEGIN\n");

    //  Print_Local_Overlap expects unadjusted chain positions.  We
    //  unadjust, the readjust.

    for(i=0;i<O->num_pieces;i++){
      O->chain[i].piece.abpos--;
      O->chain[i].piece.bbpos--;
      O->chain[i].piece.aepos--;
      O->chain[i].piece.bepos--;

      fprintf(stderr, "overlap segment %d: (%d,%d)[---------](%d,%d) gap(%d,%d)\n",
              i,
              O->chain[i].piece.abpos, O->chain[i].piece.bbpos,
              O->chain[i].piece.aepos, O->chain[i].piece.bepos,
              O->chain[i].agap, O->chain[i].bgap);
    }

    fprintf(stderr, "overlap end gaps (%d,%d)\n",
            O->chain[O->num_pieces].agap,
            O->chain[O->num_pieces].bgap);

    Print_Local_Overlap(stderr,O,5);
    //Print_Local_Overlap_withAlign(stderr,O,Ausable,Busable);

    for(i=0;i<O->num_pieces;i++){
      O->chain[i].piece.abpos++;
      O->chain[i].piece.bbpos++;
      O->chain[i].piece.aepos++;
      O->chain[i].piece.bepos++;
    }

    fprintf(stderr,"========================================Overlap after trimming END\n");
#endif


    //  Can be changed in AS_Local_Trace
    ahang = O->begpos;
    bhang = O->endpos;

    //  FIX_END_GAPS
    {
      int i=0;
      int j=0;
      int changeahang=0;
      int changebhang=0;
      int fullLenA = strlen(a->sequence);
      int fullLenB = strlen(b->sequence);
#ifdef DEBUG_FINAL_TRACE
      char c;
      fprintf(stderr, "Local_Overlap_AS Final Trace (lens %d %d):\n",fullLenA,fullLenB);
#endif
      while(trace[i]!=0){
#ifdef DEBUG_FINAL_TRACE
       	c='*';
#endif

	if(trace[i]<-fullLenA){
	  changebhang++;
	} else if (trace[i]>fullLenB){
	  changebhang--;
	} else if (trace[i]==-1){
	  changeahang--;
	} else if (trace[i]==1){
	  changeahang++;
	} else {
#ifdef DEBUG_FINAL_TRACE
	  c=' ';
#endif
	  trace[j++]=trace[i];
	}
#ifdef DEBUG_FINAL_TRACE
	fprintf(stderr, " %c%d",c,trace[i]);
#endif
	i++;
      }
#ifdef DEBUG_FINAL_TRACE
      fprintf(stderr, "\n");
#endif
      trace[j]=0;
      ahang+=changeahang;
      bhang+=changebhang;
    }


    if( ahang > end|| ahang < beg ){
      forcenext=1;
#ifdef DEBUG_LOCALOVL
      fprintf(stderr,"OVERLAP REJECTED because: ahang problems ahang=%d, > end=%d, < beg=%d\n\n",
              ahang, end, beg);
#endif
      goto next_overlap;
    }


    if (ahang < 0 || (ahang == 0 && bhang > 0))
      { if (bhang >= 0)
          QVBuffer.overlap_type = AS_CONTAINMENT;
        else
          QVBuffer.overlap_type = AS_DOVETAIL;
        QVBuffer.ahg = -ahang;
        QVBuffer.bhg = -bhang;
        if (opposite)
          QVBuffer.orientation.setIsOuttie();
        else
          QVBuffer.orientation.setIsNormal();
        QVBuffer.aifrag = b->iaccession;
        QVBuffer.bifrag = a->iaccession;
        if (trace != NULL)
          { int i;
            for (i = 0; trace[i] != 0; i++)
              trace[i] = -trace[i];
          }
      }
    else
      { if (bhang <= 0)
          QVBuffer.overlap_type = AS_CONTAINMENT;
        else
          QVBuffer.overlap_type = AS_DOVETAIL;
        QVBuffer.ahg = ahang;
        QVBuffer.bhg = bhang;
        if (opposite)
          QVBuffer.orientation.setIsInnie();
        else
          QVBuffer.orientation.setIsNormal();
        QVBuffer.aifrag = a->iaccession;
        QVBuffer.bifrag = b->iaccession;
      }

    assert(trace!=NULL);

    QVBuffer.alignment_trace = trace;
    //QVBuffer.alignment_delta = 0L;

    *where = ahang;

#ifdef DEBUG_FINAL_TRACE
    {
      fprintf(stderr, "Local_Overlap_AS Final Final Trace:\n");
      int i;
      for (i = 0; trace[i] != 0; i++)
        fprintf(stderr, " %d",trace[i]);
      fprintf(stderr, "\n");
    }
#endif

  }


  //  Check the final quality


  {
    int alen,blen,del,sub,ins,affdel,affins,blockdel,blockins;
    double errRate,errRateAffine;
    int AFFINEBLOCKSIZE=4;

    //////////////////
    // TEST FINAL FLAPS
    //////////////////

    if(asymmetricEnds){

      /* Test for Mike's clear-range-extended overlaps */

      /* For the beginning of the overlap, the overlap should be rejected if:
         - the first N bases of the B sequence are not aligned to A
         - N > MaxBegGap (which should be set to the number of bases
         the clear range was extended, plus a bit of slop)
         - there is a "flap" (unaligned bases) at the beginning of
         the alignment of sufficient size; this distinguishes cases
         which simply have a negative ahang from cases where there is
         a problematic mismatch */

      // if there is a flap of sufficient size
      if(MAX(O->chain[0].agap,O->chain[0].bgap)>MaxFreeFlap&&
	 // and the the unaligned amount of B is too long ...
	 O->chain[0].bgap+MAX(0,-O->begpos)>MaxBegGap){
#ifdef DEBUG_LOCALOVL
	fprintf(stderr,"OVERLAP REJECTED because: begin flap too big %d\n\n",
                O->chain[0].bgap+MAX(0,-O->begpos)>MaxBegGap);
#endif
	forcenext=1;
	goto next_overlap;
      }

      /* For the end of the overlap, as above but reverse orientation */

      // if there is a flap of sufficient size
      if(MAX(O->chain[O->num_pieces].agap,O->chain[O->num_pieces].bgap)>MaxFreeFlap&&
	 // and the the unaligned amount of A is too long ...
	 O->chain[O->num_pieces].agap+MAX(0,-O->endpos)>MaxEndGap){

#ifdef DEBUG_LOCALOVL
	fprintf(stderr,"OVERLAP REJECTED because: end flap too big %d\n\n",
		O->chain[O->num_pieces].agap+MAX(0,-O->endpos)>MaxEndGap);
#endif

	forcenext=1;
	goto next_overlap;

      }


    } else {  //  asymmetricEnds == 0


      if(MIN(O->chain[0].agap+MAX(O->begpos,0),
	     O->chain[0].bgap+MAX(-O->begpos,0))>MaxBegGap){
#ifdef DEBUG_LOCALOVL
	fprintf(stderr,"OVERLAP REJECTED because: begin gap too big %d\n\n",
		MIN(O->chain[0].agap, O->chain[0].bgap));
#endif
	forcenext=1;
	goto next_overlap;
      }

      if(MIN(O->chain[O->num_pieces].agap+MAX(-O->endpos,0),
	     O->chain[O->num_pieces].bgap+MAX(O->endpos,0))>MaxEndGap){
#ifdef DEBUG_LOCALOVL
	fprintf(stderr,"OVERLAP REJECTED because: end gap too big %d\n\n",
		MIN(O->chain[O->num_pieces].agap, O->chain[O->num_pieces].bgap));
#endif

	forcenext=1;
	goto next_overlap;

      }


    }  //  if asymmetricEnds


    //////////

    if (opposite)
      reverseComplement(b->sequence, b->quality, strlen(b->sequence));

    Analyze_Affine_Overlap_AS(a,
                              b,
                              &QVBuffer,
                              AS_ANALYZE_ALL,
                              &alen,
                              &blen,
                              &del,
                              &sub,
                              &ins,
			      &affdel,
                              &affins,
                              &blockdel,
                              &blockins,
                              AFFINEBLOCKSIZE,
                              &max_indel_AS_ALN_LOCOLAP_GLOBAL);

    if (opposite)
      reverseComplement(b->sequence, b->quality, strlen(b->sequence));

    errRate       = (sub+   ins+   del) / (double)(alen+      ins);
    errRateAffine = (sub+affins+affdel) / (double)(alen-del+affins+affdel);


    if(errRateAffine>erate){
#ifdef DEBUG_LOCALOVL
      fprintf(stderr,"OVERLAP REJECTED because: affine-adjusted error rate %f > erate %f\n\n",
              errRateAffine,erate);
#endif
      forcenext=1;
      goto next_overlap;
    }
    //    fprintf(stderr, "avgerror came to %f\n",avgerror);
#ifdef AFFINE_QUALITY
    QVBuffer.quality = errRateAffine;
#endif


#ifdef DEBUG_FINAL_TRACE
    {
      fprintf(stderr, "Inner: Alen %d, Blen %d, del %d, sub %d, ins %d\n", alen,blen,del,sub,ins);
      fprintf(stderr, "affdel %d, affins %d, blockdel %d, blockins %d\n", affdel,affins,blockdel,blockins);
      fprintf(stderr, "Simple mismatch rate %f\n",errRate);
      fprintf(stderr, "Affine mismatch rate %f\n",errRateAffine);
      fprintf(stderr, "Local_Overlap_AS Final Final Final Trace:\n");
      int i;
      for (i = 0; trace[i] != 0; i++)
        fprintf(stderr, " %d",trace[i]);
      fprintf(stderr, "\n");
    }
#endif

  }


  if (opposite)
    reverseComplement(b->sequence, b->quality, strlen(b->sequence));

  safe_free(O);

  return (&QVBuffer);

 nooverlap:

  if (opposite)
    reverseComplement(b->sequence, b->quality, strlen(b->sequence));

  safe_free(O);

  return((OverlapMesg*)NULL);
}
