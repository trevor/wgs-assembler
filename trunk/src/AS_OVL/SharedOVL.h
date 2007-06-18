
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

/*************************************************
* Module:  SharedOVL.h
* Description:
*   Shared declarations for overlap programs
* 
*    Programmer:  A. Delcher
*       Started:   15 February 2007
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: SharedOVL.h,v 1.3 2007-06-18 13:16:05 adelcher Exp $
 * $Revision: 1.3 $
*/


#ifndef  __SHAREDOVL_H_INCLUDED
#define  __SHAREDOVL_H_INCLUDED


#include "AS_OVL_delcher.h"


// Constants

#define  DEFAULT_CHAR_MATCH_VALUE  ERR_FRACTION_IN_AS_GLOBAL_H
  //  Default value to add for a match in scoring alignments.
  //  Corresponding error value is this value minus 1.0
  //  Using integer values didn't make alignments faster on DEC Alphas.
  //  An alignment with a matches and b mismatches scores
  //  (a + b) * Match_Value - b = a * Match_Value + b (Match_Value - 1.0)
  //     = a * Match_Value + b * Mismatch_Value
  //  Letting x = Match_Value
  //  a zero score occurs when 0 = ax + b(x-1) or
  //  a/b = (1-x)/x
  //  Defining p = b/(a+b) = Error_Rate gives
  //  1/p = (a+b)/b = 1 + a/b = 1 + (1-x)/x = 1/x
  //  whence p = x.

#define  FRAG_LEN_BITS            15
  //  Number of bits to store lengths and positions on fragments
#define  MIN_BRANCH_END_DIST      20
  // Branch points must be at least this many bases from the
  // end of the fragment to be reported
#if ERR_MODEL_IN_AS_GLOBAL_H > 6
  #define  MIN_BRANCH_TAIL_SLOPE  1.0
#else
  #define  MIN_BRANCH_TAIL_SLOPE  0.20
#endif
  // Branch point tails must fall off from the max by at least this rate
#define  MAX_ERRORS               (1 + (int) (MAX_ERROR_RATE * MAX_FRAG_LEN))
  // Most errors in any edit distance computation


// Type definitions

typedef  enum
  {
   DELETE, A_SUBST, C_SUBST, G_SUBST, T_SUBST,
   A_INSERT, C_INSERT, G_INSERT, T_INSERT, NO_VOTE,
   EXTENSION
  }  Vote_Value_t;

typedef  struct
  {
   unsigned  is_ID : 1;
   unsigned  keep_left : 1;     // set true if left overlap degree is low
   unsigned  keep_right : 1;    // set true if right overlap degree is low
   unsigned  iid : 28;
  }  Frag_ID_t;

typedef  struct
  {
   unsigned  is_ID : 1;
   unsigned  pos : 20;    // position in fragment
   unsigned  type : 11;
  }  Correction_t;

typedef  union
  {
   Frag_ID_t  frag;
   Correction_t  corr;
  }  Correction_Output_t;

typedef struct
  {
   unsigned  len : 12;
   unsigned  action : 2;   // 0,1,2,3 = insert,delete,substitute,noop resp.
   unsigned  ch : 2;       // 0,1,2,3 = a,c,g,t resp.
  }  Diff_Entry_t;

typedef struct
  {
   unsigned  a_lo : FRAG_LEN_BITS;
   unsigned  a_hi : FRAG_LEN_BITS;
   unsigned  b_lo : FRAG_LEN_BITS;
   unsigned  b_hi : FRAG_LEN_BITS;
   uint32  diff_len;
   Diff_Entry_t  * de;
  }  Sequence_Diff_t;


// Function prototypes

int  Fwd_Prefix_Edit_Dist
  (char A [], int m, char T [], int n, int Error_Limit,
   int * A_End, int * T_End, int * Match_To_End,
   double match_value, int * Delta, int * Delta_Len, int * * edit_array,
   int edit_match_limit [], int error_bound [], int doing_partial);
int  OVL_Max_int
  (int a, int b);
int  OVL_Min_int
  (int a, int b);
int  Rev_Prefix_Edit_Dist
  (char a_string [], int m, char t_string [], int n, int error_limit,
   int * a_end, int * t_end, int * leftover, int * match_to_end,
   double match_value, int * delta, int * delta_len, int * * edit_array,
   int edit_match_limit [], int error_bound [], int doing_partial);
void  Set_Fwd_Delta
  (int delta [], int * delta_len, int * * edit_array,
   int e, int d);
void  Set_Rev_Delta
  (int delta [], int * delta_len, int * * edit_array,
   int e, int d, int * leftover, int * t_end, int t_len);
int  Sign
  (int a);


#endif
