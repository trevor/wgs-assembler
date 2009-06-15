
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, The Venter Institute. All rights reserved.
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

#ifndef INCLUDE_AS_BOG_CHUNKGRAPH
#define INCLUDE_AS_BOG_CHUNKGRAPH

static const char *rcsid_INCLUDE_AS_BOG_CHUNKGRAPH = "$Id: AS_BOG_ChunkGraph.hh,v 1.20 2009-06-15 07:01:37 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"

#define PROMISCUOUS

struct BestOverlapGraph;

struct ChunkGraph{
public:
  ChunkGraph(FragmentInfo *fi, BestOverlapGraph *bovlg);
  ~ChunkGraph(void) {
    delete [] _chunk_lengths;
  };

  uint32 nextFragByChunkLength(void) {
    static uint32 pos = 0;

    if (pos < _max_fragments)
      return _chunk_lengths[pos++].fragId;

    pos = 0;
    return(0);
  };

private:
  uint32 countFullWidth(BestOverlapGraph *BOG, uint32 *, uint32, uint32 );

  struct _chunk_length {
    uint32 fragId;
    uint32 cnt;

    bool operator<(_chunk_length const that) const {
      if (cnt == that.cnt)
        return(fragId < that.fragId);
      return(cnt > that.cnt);
    };
  };
  _chunk_length      *_chunk_lengths;

  uint32              _max_fragments;
};

#endif

