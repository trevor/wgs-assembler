
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
#ifndef OUTPUT_CGW_H
#define OUTPUT_CGW_H

void MarkContigEdges(void);
void OutputMateDists(ScaffoldGraphT *graph);
void OutputFrags(ScaffoldGraphT *graph);
void OutputUnitigs(ScaffoldGraphT *graph);
void OutputUnitigLinks(ScaffoldGraphT *graph);
void OutputContigs(ScaffoldGraphT *graph);
void OutputContigLinks(ScaffoldGraphT *graph, int outputOverlapOnlyContigEdges);
void OutputScaffolds(ScaffoldGraphT *graph);
void OutputScaffoldLinks(ScaffoldGraphT *graph);

void OutputUnitigLinksFromMultiAligns(void);
void OutputContigLinksFromMultiAligns(int outputOverlapOnlyContigEdges);
void OutputContigsFromMultiAligns(void);
void OutputUnitigsFromMultiAligns(void);

#endif
