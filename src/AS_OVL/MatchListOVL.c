
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
/*********************************************************************
   Module:  AS_OVL
   Description:  Assembly Overlap Module--computes overlaps between
      pairs of DNA strings.
 *********************************************************************/

/* RCS info
 * $Id: MatchListOVL.c,v 1.3 2005-03-22 19:07:11 jason_miller Exp $
 * $Revision: 1.3 $
*/


// Component:
//
//   MatchListOVL.c
//
//   Last revised:  10 July 2001
//
// Description:
// 
//   These routines implement the data structure that holds matches
//   between fragments that could potentially overlap.  A match is
//   an exact-match substring between two fragments.
// 
// Design:
// 
// 
// Limitations:
// 
// 
// Status:
// 
// 
// Architecture and Dependencies:
// 
//   MatchListOVL.h   Header for this file
//   MatchListOVL.c   This file
//


#include  <stdlib.h>
#include  <stdio.h>
#include  <assert.h>
#include  <fcntl.h>
#include  <sys/types.h>
#include  <string.h>
#include  <dirent.h>
#include  <sys/stat.h>
#include  <unistd.h>

#include  "AS_OVL_delcher.h"
#include  "MatchListOVL.h"


int  Is_Empty
    (const Match_Set_t * matches)

//  Returns  TRUE  iff  matches  has no elements.

  {
   return  (matches -> first == 0);
  }

