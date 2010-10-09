
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
/* $Id: MPTypes.h,v 1.1.1.1 2004-04-14 13:52:04 catmandew Exp $ */
#ifndef MPTYPES_H
#define MPTYPES_H

#ifdef i386
#define _LU "%llu"
#define _LS "%ll"
typedef unsigned long long uint64;
typedef long long int64;
#else
#define _LU "%lu"
#define _LS "%l"
typedef unsigned long uint64;
typedef long int64;
#endif

typedef unsigned int uint32;
typedef int int32;
typedef float float32;

#define ID_TYPE uint64
#define ID_SCAN_FORMAT _LU
#define UNIT_TYPE int64
#define UNIT_SCAN_FORMAT  _LS

#define COINCIDENT_THRESHOLD  50  // bp
#define CONFIRMATION_THRESHOLD 2  // 2 agreeing unsatisfieds confirm each other
#define STDDEVS_THRESHOLD 3.0     // mean +/- N stddevs

typedef enum
{
  ICS_ALL,
  ICS_HIGHER,
  ICS_LOWER,
  ICS_NUM_SCOPES
} InterChromosomeScope;

#define INTER_C_DEFAULT_SCOPE  ICS_ALL

#endif //  MPTYPES_H