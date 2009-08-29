
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

#ifndef INPUT_CGW_H
#define INPUT_CGW_H

static const char *rcsid_INPUT_CGW_H = "$Id: Input_CGW.h,v 1.10 2009-07-30 10:42:55 brianwalenz Exp $";

int ProcessInput(Global_CGW *data, int optind, int argc, char *argv[]);

void ProcessIUM_ScaffoldGraph(IntUnitigMesg *ium_mesg,
                              int32 length,
                              int sequenceOnly);

void  LoadDistData(void);

void LoadClosureReadData(void);

#endif
