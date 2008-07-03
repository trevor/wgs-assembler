
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
   Description:  Assembly Overlap Module.  Computes overlaps between
      pairs of DNA strings.
   Assumptions:  Input meets specifications in the ProtoIO documents
 *********************************************************************/
static char CM_ID[] = "$Id: AS_OVL_driver.c,v 1.6 2008-06-27 06:29:17 brianwalenz Exp $";

/* RCS info
 * $Id: AS_OVL_driver.c,v 1.6 2008-06-27 06:29:17 brianwalenz Exp $
 * $Revision: 1.6 $
*/


/******************************************************************************/
/* AS_OVL_driver
     This is the driver for the algorithmic overlap code in AS_OVL_overlap.
     Invoked by main.
*/


/* this is the frag version binary */
#undef CONTIG_OVERLAPPER_VERSION

#include "AS_OVL_driver_common.h"

