
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
   CVS_ID:  $Id: AS_TER_terminator_funcs.h,v 1.2 2004-09-23 20:25:29 mcschatz Exp $
 *********************************************************************/
#ifndef AS_TER_TERMINATOR_FUNCS_H
#define AS_TER_TERMINATOR_FUNCS_H

#include "AS_global.h"

void output_snapshot(char* fragStoreName, char* bactigStoreName,
		     char* gkpStoreName,
		     char** inputFileList, int32 numInputFiles, 
		     char* outputFileName,
		     char* mapFileName,
		     OutputType output, 
		     int32 blockSize, int32 real, int32 quiet,
		     int32 random, uint64 uidStart,
		     int argc, char *argv[]);
/******************************************************************/
// Outputs the genome snapshot
// The parameter output specifies whether the out put should be
// ACII or binary. 
// blockSize is the number of UIDs requested in one request
// real is a boolean that indicates whether real or dummy UIDs
// should be assigned. If bactigStoreName is not NUll the terminator
// is in CA mode.
//
// Precondition: filePrefix and fragStoreName are not NULL and do 
// not point to an empty string.
// 
/*****************************************************************/

#endif





