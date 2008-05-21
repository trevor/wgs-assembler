
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
#ifndef AS_ALN_FORCNS_H
#define AS_ALN_FORCNS_H

#include "AS_ALN_aligners.h"

/* 

   Wrappers for finding fragment overlaps with moderate sized indels
   ("bubbles" in CGB that break up unitigging in the presence of
   moderate polymorphisms).

   The wrappers are around the routines Local_Overlap_AS  ([originally]
   defined in AS_ALN_loverlapper.c) and AS_ALN_affine_overlap ([originally]
   defined in AS_ALN_qvaligner.c).  See those functions for more details.

*/

Overlap *Local_Overlap_AS_forCNS(char *aseq, char *bseq,
                    int beg, int end, int opposite,
                    double erate, double thresh, int minlen,
				 CompareOptions what);


Overlap *Affine_Overlap_AS_forCNS(char *aseq, char *bseq,
                    int beg, int end, int opposite,
                    double erate, double thresh, int minlen,
				 CompareOptions what);

#endif
