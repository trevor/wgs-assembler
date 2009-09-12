
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

#ifndef REPEAT_REZ_H
#define REPEAT_REZ_H

static const char *rcsid_REPEAT_REZ_H = "$Id: RepeatRez.h,v 1.7 2009-09-12 22:35:58 brianwalenz Exp $";

int  Fill_Gaps
    (char *, int, int redo_index);

int  Show_Reads_In_Gaps
    (char * prefix);

int  Hurl_Contained_Rocks
    (char * prefix, int level, int redo_index);

#define AGGRESSIVE_WALKING_STD_DEVS 5.0
#define CONSERVATIVE_WALKING_STD_DEVS 3.0

int  Walk_Gaps
    (char *, int, int startWalkFrom, double gapSizeStdDevs);

int  Throw_Stones
    (char *, int, int);

int  Toss_Contained_Stones
    (char * prefix, int level, int redo_index);

int Inter_Scaffold_Walking(void);

#endif
