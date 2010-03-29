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

/* RCS info
 * $Id: AS_BOG_MateChecker.hh,v 1.19 2007-12-05 23:46:57 brianwalenz Exp $
 * $Revision: 1.19 $
 */

#ifndef INCLUDE_AS_BOG_MATECHEKER
#define INCLUDE_AS_BOG_MATECHEKER

#include "AS_BOG_UnitigGraph.hh"

extern "C" {
#include "AS_PER_gkpStore.h"
}

namespace AS_BOG{
    typedef std::map<iuid,iuid> IdMap;
    typedef IdMap::iterator IdMapIter;
    typedef IdMap::const_iterator IdMapConstIter;
    typedef std::vector<int> DistanceList;
    typedef DistanceList::const_iterator DistanceListCIter;
    typedef std::map<iuid,DistanceList> LibraryDistances;
    typedef LibraryDistances::const_iterator LibDistsConstIter;

    struct MateInfo {
        iuid mate;
        iuid lib;
    };
    typedef std::map<iuid,MateInfo> MateMap;
    static const MateInfo NULL_MATE_INFO = {0,0};
    inline bool operator==(MateInfo a, MateInfo b) {
        if (a.mate == b.mate)
            return(a.lib == b.lib);
        return false;
    }
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
        friend std::ostream& operator<< (std::ostream&, MateCounts);

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

    ///////////////////////////////////////////////////////////////////////////

    struct MateChecker{
        ~MateChecker();

        iuid readStore(const char *);   // reads the gkpStore mate info into memory

        // returns the mate iid for the given iid, or zero if none
        MateInfo getMateInfo(iuid);

        // Checks size of mates internal to unitig
        LibraryStats* computeLibraryStats( Unitig* );

        // Compute good and bad coverage graphs for a unitig, returns split points
        FragmentEnds* computeMateCoverage( Unitig*, BestOverlapGraph *);

        // Make singleton unitigs out of unhappy contained reads
        void ejectUnhappyContains( UnitigGraph& );

        // Computes stddev and mate coverage over all unitigs
        void computeGlobalLibStats( UnitigGraph& );

        // Main entry point for running mate splitting
        void checkUnitigGraph( UnitigGraph& );

    private:
        MateMap _mates;
        LibraryDistances _dists; // all distances 
        LibraryStats _globalStats;
    };

    ///////////////////////////////////////////////////////////////////////////

    struct MateLocationEntry {
        SeqInterval pos1;
        SeqInterval pos2;
        iuid        id1;
        iuid        id2;
        iuid        unitig1;
        iuid        unitig2; // in the future the table might be across unitigs
        bool        isBad;
    };
    std::ostream& operator<< (std::ostream& os, MateLocationEntry&);

    static const MateLocationEntry NULL_MATE_ENTRY =
        {NULL_SEQ_LOC,NULL_SEQ_LOC,0,0,0,0};

    typedef std::vector<MateLocationEntry> MateLocTable;
    typedef MateLocTable::iterator MateLocIter;
    typedef MateLocTable::const_iterator MateLocCIter;

    class MateLocation {
    public:

        MateLocation(MateChecker* const check) : _checker(check) {
            goodGraph   = new std::vector<short>;
            badFwdGraph = new std::vector<short>;
            badRevGraph = new std::vector<short>;
        };
        ~MateLocation();

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

        MateLocTable _table;
        IdMap _iidIndex;
        MateChecker* _checker;
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

}
#endif