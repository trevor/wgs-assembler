
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
/**************************************************************
 * AS_CNS/Utilities_CNS.h
 *
 * Utility functions for the CNS (consensus) subsystem 
 * of the Celera WGS assembler.
 *
 **************************************************************/
/*********************************************************************
 $Id: Utilities_CNS.h,v 1.1.1.1 2004-04-14 13:51:22 catmandew Exp $
 *********************************************************************/

#ifndef UTILITIES_CNS_INCLUDE
#define UTILITIES_CNS_INCLUDE

#include "PublicAPI_CNS.h"

void CleanExit(char *mesg, int lineno, int rc) ;

#endif

