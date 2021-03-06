
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

/* NOTE: Everything a user needs to find and squash bubbles can be
   found in this file!  The other AS_CGB_Bubble_*.h files are for
   internal use only! */

#ifndef _AS_CGB_BUBBLE_H_
#define _AS_CGB_BUBBLE_H_

#include "AS_ALN_aligners.H"

static const char *rcsid__AS_CGB_BUBBLE_H_ = "$Id$";

/* For debugging only.  Set this to 1 only for test-size problems ... :-) */
#define AS_CGB_BUBBLE_VERBOSE 0

/* For traversal debugging only.  Produces too much output. Use only for
   teensy-weensy test-size problems ... */
#define AS_CGB_BUBBLE_VERY_VERBOSE 0

/* For outputting distribution bubble size distribution data to a
   file.  Set to 1 to create three files, "bubble.nfrags.celagram",
   "bubble.poppedfrags.celagram", and "bubble.sdiff.celagram", with the number
   of fragments in each bubble, the number of fragments in each popped
   bubble, and the largest block difference in an alignment within
   each bubble respectively. */
#define AS_CGB_BUBBLE_DIST_OUTPUT 1

/*
 *
 * YOUR ALL-IN-ONE, ONE-STOP, ONE-SIZE-FITS-ALL SUPERSTORE FOR BUBBLE SQUASHING
 *
 */

/* This method uses DAB and BOPPER to remove all possible bubbles from
   the graph.  It takes the fragments, edges, and chunks as input, as
   well as the handle of an output file into which to write the
   necessary overlaps, and a file into which to record log messages
   and statistics.

   VERY IMPORTANT NOTE:  All of the overlaps produced by bubble popping
   have NEGATIVE quality to indicate that they were generated by affine
   alignment as opposed to normal overlapping.  It is therefore important
   to take the ABS() of the quality before using it for comparison purposes. */
void
AS_CGB_Bubble_find_and_remove_bubbles
(gkStore *gkpStore,
 Tfragment *frags, Tedge *edges,
 TChunkMesg *chunks, TChunkFrag *cfrgs,
 float global_arrival_rate,
 const char * bubblename,
 const char * fileprefix);

/*
 *
 * The rest is for do-it-yourself types only ...
 *
 */

/*
 *
 * ROUTINES FOR BUBBLE FINDING
 *
 */

typedef struct AS_CGB_Bubble_List {
  IntFragment_ID start, end;
  int start_sx, end_sx;
  struct AS_CGB_Bubble_List *next;
} AS_CGB_Bubble_List;

typedef AS_CGB_Bubble_List * AS_CGB_Bubble_List_t;

/* This method simply returns a list of bubble's found by Dan's Automated
   Bubblefinder (DAB: A little DAB'll do ya!).  The start and end fragments
   are specified by fragment CGB Vertex ID (VID).

   More documentation on those odd numerical parameters forthcoming ... */
AS_CGB_Bubble_List_t
AS_CGB_Bubble_find_bubbles(Tfragment *frags, Tedge *edges, int sz, int age,
			   int max_outdegree);

/*
 *
 * ROUTINES FOR BUBBLE REMOVAL ("POPPING")
 *
 */

typedef struct BubblePopper * BubblePopper_t;


/* This method attempts to squash the given bubble with the BOPPER
   algorithm.  The bubble is specified by the start and end fragment
   CGB Vertex IDs (VIDs), as well as flags which indicate whether the
   bubble touches the suffix (=1) or the prefix (=0) of each fragment.

   The return value is either a pointer to an array of overlap messages,
   or NULL if the bubble could not be successfully squashed.
   Statistics about the bubble are accumulated in the BubblePopper,
   but the overlap array IS NOT persistent; it will be recycled with
   the next call to AS_CGB_Bubble_pop_bubble().  The number of
   overlaps found is placed in the num_olaps output parameter, which
   must be non-null.

   VERY IMPORTANT NOTE:  All of the overlaps produced by bubble popping
   have NEGATIVE quality to indicate that they were generated by affine
   alignment as opposed to normal overlapping.  It is therefore important
   to take the ABS() of the quality before using it for comparison purposes.

   Implemented in "AS_CGB_Bubble_PopperMethods.c" */
ALNoverlapFull *
AS_CGB_Bubble_pop_bubble(BubblePopper_t bp, IntFragment_ID start,
			 int start_sx, IntFragment_ID end,
			 int end_sx, int *num_olaps);

/*
 * Should place BubblePopper accessor functions (particularly for
 * statistics) here.  Any demand?
 */

/* Deallocates a bubble popper.

   Implemented in "AS_CGB_Bubble_Popper.c" */
void
AS_CGB_Bubble_Popper_destroy(BubblePopper_t bp);

#endif
