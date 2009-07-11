
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
#ifndef MULTIALIGNMENT_CNS_INCLUDE
#define MULTIALIGNMENT_CNS_INCLUDE

static const char *rcsid_MULTIALIGNMENT_CNS_INCLUDE = "$Id: MultiAlignment_CNS.h,v 1.52 2009-07-11 00:20:30 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_MSG_pmesg.h"
#include "AS_SDB_SequenceDB.h"
#include "MultiAlignStore_CNS.h"
#include "AS_ALN_aligners.h"


//  This is probably broken, or extremely inefficient, as of Nov 4 2007.
#undef PRINTUIDS

extern int DUMP_UNITIGS_IN_MULTIALIGNCONTIG;
extern int VERBOSE_MULTIALIGN_OUTPUT;

#define CNS_OPTIONS_SPLIT_ALLELES_DEFAULT  1
#define CNS_OPTIONS_MIN_ANCHOR_DEFAULT    11

typedef struct {
  int split_alleles;
  int smooth_win;
} CNS_Options;

typedef enum {
  CNS_QUIET       = (int)'Q', // quiet,  print nothing
  CNS_STATS_ONLY  = (int)'S', // print only 1-line statistic summary
  CNS_ALIGNMENT   = (int)'A', // print the multialignment, sans CNS
  CNS_CONSENSUS   = (int)'C', // print the multialignment, with CNS
  CNS_DOTS        = (int)'D', // print the multialignment, dot format
  CNS_NODOTS      = (int)'N', // print the multialignment, "nodot" format
  CNS_EDIT_SCORE  = (int)'E', // print the edit score column by column
  CNS_VIEW_UNITIG = (int)'U',  // show the unitigs in the contig alignment
  CNS_VERBOSE     = (int)'V'  // verbose pre-post refinment output
} CNS_PrintKey;   // determine the format for PrintAlignment

typedef enum {
  CNS_SMOOTH = 1, // only eliminate pairwise construction artifacts
  CNS_POLYX  = 2, // align poly-X regions
  CNS_INDEL  = 4  // push apart mushed block indels
}  CNS_RefineLevel;


MultiAlignT *MergeMultiAlignsFast_new(tSequenceDB *,
                                      gkStore *,
                                      VA_TYPE(IntElementPos) *,
                                      int,
                                      int,
                                      CNS_Options *opp);

MultiAlignT *ReplaceEndUnitigInContig(tSequenceDB *,
                                      gkStore * ,
                                      uint32,
                                      uint32,
                                      int,
                                      CNS_Options *opp);


int MultiAlignUnitig(IntUnitigMesg *,
                     gkStore *,
                     VA_TYPE(char) *,
                     VA_TYPE(char) *,
                     VA_TYPE(int32) *,
                     CNS_PrintKey,
                     CNS_Options *opp);

int MultiAlignContig(IntConConMesg *,
                     VA_TYPE(char) *,
                     VA_TYPE(char) *,
                     VA_TYPE(int32) *,
                     CNS_PrintKey ,
                     CNS_Options *opp);

//  Options to things in MultiAligment_CNS.c

extern int allow_neg_hang;

#endif
