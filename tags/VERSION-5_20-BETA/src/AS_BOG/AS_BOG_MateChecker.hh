
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

#ifndef INCLUDE_AS_BOG_MATECHEKER
#define INCLUDE_AS_BOG_MATECHEKER

static const char *rcsid_INCLUDE_AS_BOG_MATECHEKER = "$Id: AS_BOG_MateChecker.hh,v 1.29 2008-10-08 22:02:54 brianwalenz Exp $";

#include "AS_BOG_UnitigGraph.hh"

extern "C" {
#include "AS_PER_gkpStore.h"
}

typedef std::map<iuid,iuid> IdMap;
typedef IdMap::iterator IdMapIter;
typedef IdMap::const_iterator IdMapConstIter;

typedef std::vector<int> DistanceList;
typedef DistanceList::const_iterator DistanceListCIter;

typedef std::map<iuid,DistanceList> LibraryDistances;
typedef LibraryDistances::const_iterator LibDistsConstIter;

static const SeqInterval NULL_SEQ_LOC = {0,0};

struct DistanceCompute {
  double stddev;
  double mean;
  double sumSquares;
  double sumDists;
  int numPairs;
  DistanceCompute() : stddev(0), mean(0), sumSquares(0), sumDists(0), numPairs(0) {}
};

typedef std::map<iuid,DistanceCompute> LibraryStats;

struct MateCounts {
  int badOtherTig;
  int otherTig;
  int total;
  int goodCircular;
  int good;
  int badOuttie;
  int badInnie;
  int badAntiNormal;
  int badNormal;
  MateCounts() : badOtherTig(0), otherTig(0), total(0), goodCircular(0), good(0),
                 badOuttie(0), badInnie(0), badAntiNormal(0), badNormal(0)
  {}

  MateCounts operator+= (MateCounts other) {
    badOtherTig   += other.badOtherTig;
    otherTig      += other.otherTig;
    total         += other.total;
    goodCircular  += other.goodCircular;
    good          += other.good;
    badOuttie     += other.badOuttie;
    badInnie      += other.badInnie;
    badAntiNormal += other.badAntiNormal;
    badNormal     += other.badNormal;
  };
};


struct MateChecker{
  MateChecker(FragmentInfo *fi);
  ~MateChecker();

  void checkUnitigGraph(UnitigGraph &, int badMateBreakThreshold);

private:

  // Checks size of mates internal to unitig
  LibraryStats* computeLibraryStats( Unitig* );

  // Compute good and bad coverage graphs for a unitig, returns split points
  UnitigBreakPoints* computeMateCoverage(Unitig *, BestOverlapGraph *, int badMateBreakThreshold);

  void moveContains(UnitigGraph&);
  void splitDiscontinuousUnitigs(UnitigGraph&);

  // Computes stddev and mate coverage over all unitigs
  void computeGlobalLibStats( UnitigGraph& );

private:
  LibraryDistances _dists; // all distances
  LibraryStats     _globalStats;
  FragmentInfo    *_fi;
};


struct MateLocationEntry {
  SeqInterval pos1;
  SeqInterval pos2;
  iuid        id1;
  iuid        id2;
  iuid        unitig1;
  iuid        unitig2; // in the future the table might be across unitigs
  bool        isBad;
};

static const MateLocationEntry NULL_MATE_ENTRY =
  {NULL_SEQ_LOC,NULL_SEQ_LOC,0,0,0,0};

typedef std::vector<MateLocationEntry> MateLocTable;
typedef MateLocTable::iterator MateLocIter;
typedef MateLocTable::const_iterator MateLocCIter;

class MateLocation {
public:

  MateLocation(FragmentInfo *fi) {
    goodGraph   = new std::vector<short>;
    badFwdGraph = new std::vector<short>;
    badRevGraph = new std::vector<short>;
    _fi         = fi;
  };
  ~MateLocation() {
    delete goodGraph;
    delete badFwdGraph;
    delete badRevGraph;
  };

  bool startEntry( iuid, iuid, SeqInterval);
  bool addMate( iuid, iuid, SeqInterval);
  bool hasFrag( iuid );
  MateLocationEntry getById( iuid );
  void sort();
  MateLocIter begin() { return _table.begin(); }
  MateLocIter end()   { return _table.end();   }

  void buildTable( Unitig *);
  MateCounts* buildHappinessGraphs( int tigLen, LibraryStats &);

  std::vector<short>* goodGraph;
  std::vector<short>* badFwdGraph;
  std::vector<short>* badRevGraph;

private:
  MateLocTable  _table;
  IdMap         _iidIndex;
  FragmentInfo *_fi;
};
inline bool operator==(SeqInterval a, SeqInterval b) {
  if (a.bgn == b.bgn && a.end == b.end ||
      a.bgn == b.end && a.end == b.bgn)
    return true;
  else
    return false;
};
inline bool operator<(SeqInterval a, SeqInterval b) {
  if ( isReverse(a) ) {
    if ( isReverse(b) ) return a.end < b.end;
    else                return a.end < b.bgn;
  } else {
    if ( isReverse(b) ) return a.bgn < b.end;
    else                return a.bgn < b.bgn;
  }
};
inline bool SeqInterval_less(SeqInterval a, SeqInterval b) {
  return a < b;
};
inline bool operator==(MateLocationEntry a, MateLocationEntry b) {
  if (a.pos1 == b.pos1 && a.pos2 == b.pos2)
    return true;
  else
    return false;
};
inline bool operator<(MateLocationEntry a, MateLocationEntry b) {
  if (a.pos1 < b.pos1)                     return true;
  if (a.pos1 == b.pos1 && a.pos2 < b.pos2) return true;
  else                                     return false;
};

#endif