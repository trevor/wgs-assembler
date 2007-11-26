
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

// $Id: AS_PER_encodeSequenceQuality.h,v 1.5 2007-02-22 00:06:57 brianwalenz Exp $

#ifndef AS_PER_ENCODESEQUENCEQUALITY_H
#define AS_PER_ENCODESEQUENCEQUALITY_H

// Functions for encoding sequence and quality values one char per
// seq/quality pair.


void
encodeSequenceQuality(char *encoded,
                      char *sequence,
                      char *quality);

void
decodeSequenceQuality(char *encoded,
                      char *sequence,
                      char *quality);

#endif
