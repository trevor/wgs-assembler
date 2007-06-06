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
 * $Id: AS_BOG_MateChecker.cc,v 1.15 2007-06-06 19:56:26 eliv Exp $
 * $Revision: 1.15 $
*/

#include <math.h>
#include "AS_BOG_MateChecker.hh"

namespace AS_BOG{
    MateChecker::~MateChecker() {
    }

    void MateChecker::readStore(const char* gkpStorePath) {
        GateKeeperStore *gkpStore = openGateKeeperStore(gkpStorePath, FALSE);

        StreamHandle frags = openStream(gkpStore->frg, NULL, 0);
        resetStream(frags, STREAM_FROMSTART, STREAM_UNTILEND);

        GateKeeperFragmentRecord gkpf;
        iuid frgIID = 1;
        while(nextStream(frags, &gkpf)){
            MateInfo mi;
            mi.mate = gkpf.mateIID;
            mi.lib  = gkpf.libraryIID;
            _mates[ gkpf.readIID ] = mi;

            frgIID++;
        }
        closeStream(frags);
        closeGateKeeperStore(gkpStore);
    }

    ///////////////////////////////////////////////////////////////////////////

    MateInfo MateChecker::getMateInfo(iuid fragID)
    {
        if (_mates.find( fragID ) == _mates.end())
            return NULL_MATE_INFO;
        else
            return _mates[ fragID ];
    }

    ///////////////////////////////////////////////////////////////////////////

    LibraryStats* MateChecker::checkUnitig(Unitig* tig)
    {
        fprintf(stderr,"Check mates for tig %ld\n",tig->id());
        IdMap goodMates;
        iuid max = std::numeric_limits<iuid>::max(); // Sentinel value
        int numInTig = 0;
        int numWithMate = 0;
        LibraryStats *libs = new LibraryStats();
        std::vector<iuid> otherUnitig; // mates to another unitig
        IdMap otherTigHist; // count of mates to the other unitigs
        // Check each frag in unitig for it's mate
        DoveTailConstIter tigIter = tig->dovetail_path_ptr->begin();
        for(;tigIter != tig->dovetail_path_ptr->end(); tigIter++)
        {
            DoveTailNode frag = *tigIter;
            iuid fragId = frag.ident;
            MateInfo mateInfo = getMateInfo( fragId );
            if ( mateInfo.mate != NULL_FRAG_ID ) {
                numWithMate++;
                if ( tig->id() == Unitig::fragIn( mateInfo.mate ) ) {
                    numInTig++;
                    SeqInterval fragPos = frag.position;
                    if (goodMates.find( fragId ) != goodMates.end()) {
                        // already seen it's mate
                        if (isReverse(fragPos)) {
                            if ( goodMates[fragId] == max )
                                continue;
                            // 2nd frag seen is reverse, so good orientation
                            long mateDist = fragPos.bgn - goodMates[fragId];
                            goodMates[fragId] = mateDist;
                            // Record the good distance for later stddev recalculation
                            iuid distId = mateInfo.lib;
                            if (_dists.find( distId ) == _dists.end()) {
                                DistanceList newDL;
                                _dists[ distId ] = newDL;
                            } 
                            _dists[ distId ].push_back( mateDist );
                            
                            // Record the distance for unitig local stddev
                            if (libs->find( distId ) == libs->end()) {
                                DistanceCompute d;
                                d.numPairs   = 1;
                                d.sumDists   = mateDist;
                                d.stddev     = 0.0;
                                d.mean       = 0.0;
                                d.sumSquares = 0.0;
                                (*libs)[ distId ] = d;
                            } else {
                                (*libs)[distId].numPairs++;
                                (*libs)[distId].sumDists+=mateDist;
                            }
                        } else {
                            // 2nd frag seen is forward, so bad
                            goodMates[fragId] = max;
                        }
                    } else { // 1st of pair
                        if (isReverse(fragPos)) {
                            // 1st reversed, so bad
                            goodMates[mateInfo.mate] = max;
                        } else {
                            // 1st forward, so good store begin
                            goodMates[mateInfo.mate] = fragPos.bgn;
                        }
                    }
                } else {
                    // mate in other unitig, just create histogram currently
                    otherUnitig.push_back( fragId );
                    iuid otherTig = Unitig::fragIn( mateInfo.mate );
                    if (otherTigHist.find( otherTig ) == otherTigHist.end())
                        otherTigHist[ otherTig ] = 1;
                    else
                        otherTigHist[ otherTig ]++;
                }
            }
        }
        fprintf(stderr,"Num frags in tig %ld\n",tig->dovetail_path_ptr->size());
        fprintf(stderr,"Num frags with mate %d\n",numWithMate);
        fprintf(stderr,"Num with mate in unitig %d\n",numInTig);
        fprintf(stderr,"Num other unitig %d\n",otherUnitig.size());
        IdMapConstIter histIter = otherTigHist.begin();
        for(;histIter != otherTigHist.end(); histIter++) {
            iuid libId = histIter->first;
            iuid cnt   = histIter->second;
            fprintf(stderr,"Num mates to unitig %ld is %ld.\n",libId,cnt);
        }

        // Calculate the unitig local mean
        LibraryStats::iterator dcIter = libs->begin();
        for(; dcIter != libs->end(); dcIter++) {
            iuid lib = dcIter->first;
            DistanceCompute *dc = &(dcIter->second);
            dc->mean = dc->sumDists / dc->numPairs;
            fprintf(stderr,"Distance lib %ld has %ld pairs with mean dist %.1f\n",
                    lib, dc->numPairs, dc->mean );
        }
        // Sum of (x-mean)^2, not yet full stddev
        IdMapConstIter goodIter = goodMates.begin();
        for(;goodIter != goodMates.end(); goodIter++)
        {
            iuid fragId = goodIter->first;
            iuid mateDist = goodIter->second;
            if (mateDist == max)
                continue;
            iuid distId = getMateInfo( fragId ).lib;
            (*libs)[distId].sumSquares+=pow(mateDist-(*libs)[distId].mean,2);
        }
        // Calculate the real stddev
        dcIter = libs->begin();
        for(; dcIter != libs->end(); dcIter++) {
            iuid lib = dcIter->first;
            DistanceCompute *dc = &(dcIter->second);
            // really need to just collect all values and calculate in checkUnitigGraph()
            if (dc->numPairs == 1)
                dc->stddev = 0.0;
            else
                dc->stddev = sqrt( dc->sumSquares / (dc->numPairs-1) );
            fprintf(stderr,"Distance lib %ld has %ld pairs with stddev %.1f\n",
                    lib, dc->numPairs, dc->stddev );
        }
        return libs;
    }

    ///////////////////////////////////////////////////////////////////////////
    // main entry point into mate checking code
    void MateChecker::checkUnitigGraph( UnitigGraph& tigGraph )
    {
        LibraryStats globalStats;
        LibraryStats::iterator dcIter;
        _dists.clear(); // reset to seperate multiple Graphs
        UnitigsConstIter tigIter = tigGraph.unitigs->begin();
        for(; tigIter != tigGraph.unitigs->end(); tigIter++)
        {
            if (*tigIter == NULL ) 
                continue;
            LibraryStats* libs = checkUnitig(*tigIter);
            // Accumulate per unitig stats to compute global stddev's
            for(dcIter = libs->begin(); dcIter != libs->end(); dcIter++) {
                iuid lib = dcIter->first;
                DistanceCompute dc = dcIter->second;
                if (globalStats.find(lib) == globalStats.end() ) {
                    globalStats[ lib ] = dc;
                }
                else {
                    DistanceCompute *gdc = &(globalStats[ lib ]);
                    gdc->numPairs   += dc.numPairs;
                    gdc->sumDists   += dc.sumDists;
                    gdc->sumSquares += dc.sumSquares;
                }
            }
            delete libs; // Created in checkUnitig
        }
        // Calculate and output overall global mean
        for(dcIter= globalStats.begin(); dcIter != globalStats.end(); dcIter++){
            iuid lib = dcIter->first;
            DistanceCompute *dc = &(dcIter->second);
            dc->mean = dc->sumDists / dc->numPairs;
            fprintf(stderr,"Distance lib %ld has global %ld pairs with mean dist %.1f\n",
                    lib, dc->numPairs, dc->mean );
        }
        // Calculate and output overall global stddev
        for(dcIter= globalStats.begin(); dcIter != globalStats.end(); dcIter++)
        {
            iuid lib = dcIter->first;
            DistanceCompute *dc = &(dcIter->second);
            if (dc->numPairs == 1)
                dc->stddev = 0.0;
            else
                dc->stddev = sqrt( dc->sumSquares / (dc->numPairs-1) );
            fprintf(stderr,"Distance lib %ld has global %ld pairs with stddev %.1f\n",
                    lib, dc->numPairs, dc->stddev );
            // Now reset the global stats to zero so the real calculation below works
            dc->numPairs   = 0;
            dc->mean       = 0.0;
            dc->stddev     = 0.0;
            dc->sumSquares = 0.0;
            dc->sumDists   = 0.0;
        }
        // Tab delimited table header
        fprintf(stderr,
"DistLib\tnumDists\tmedian\t1/3rd\t2/3rd\tmaxDiff\tmin\tmax\tnumGood\tmean\tstddev\n");

        // Disregard outliers and recalculate global stddev
        LibDistsConstIter libDistIter = _dists.begin();
        for(; libDistIter != _dists.end(); libDistIter++)
        {
            iuid libId = libDistIter->first;
            DistanceList dl = libDistIter->second;
            sort(dl.begin(),dl.end());
            int size = dl.size();
            int median   = dl[ size / 2 ];
            int third    = dl[ size / 3 ];
            int twoThird = dl[ size * 2 / 3 ];
            int aproxStd = MAX( median - third, twoThird - median);
            int biggest  = median + aproxStd * 5;
            int smallest = median - aproxStd * 5;

            // now go through the distances and calculate the real stddev
            // including everything within 5 stddevs
            iuid numBad = 0;
            DistanceCompute *gdc = &(globalStats[ libId ]);
            DistanceListCIter dIter = dl.begin();
            for(;dIter != dl.end(); dIter++) {
                if (*dIter >= smallest && *dIter <= biggest ) {
                    gdc->numPairs++;
                    gdc->sumDists += *dIter;
                } else {
                    numBad++;
                } 
            }
            gdc->mean = gdc->sumDists / gdc->numPairs;
            // Compute sum of squares for stddev
            for(dIter = dl.begin(); dIter != dl.end(); dIter++)
            {
                if (*dIter >= smallest && *dIter <= biggest )
                    gdc->sumSquares += pow( *dIter - gdc->mean, 2);
            }
            if (gdc->numPairs > 1)
                gdc->stddev = sqrt( gdc->sumSquares / (gdc->numPairs-1) );

            fprintf(stderr, "%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%.1f\t%.1f\n",
                    libId, size, median, third, twoThird, aproxStd, smallest, biggest,
                            gdc->numPairs, gdc->mean, gdc->stddev );
        }
        UnitigVector* splits = new UnitigVector(); // save split unitigs
        UnitigsIter tigEditIter = tigGraph.unitigs->begin();
        for(; tigEditIter != tigGraph.unitigs->end(); tigEditIter++)
        {
            if (*tigEditIter == NULL ) 
                continue;

            FragmentEnds* breaks = computeMateCoverage( *tigEditIter, globalStats );
            tigGraph.accumulateSplitUnitigs( tigEditIter, breaks, splits );
            delete breaks;
        }
        fprintf(stderr,"Num mate based splits %d\n",splits->size());
        tigGraph.unitigs->insert( tigGraph.unitigs->end(), splits->begin(), splits->end());
        delete splits;
    }

    ///////////////////////////////////////////////////////////////////////////

    inline void incrRange( short graph[], short val, iuid n, iuid m ) {
        for(iuid i=n; i <=m ; i++)
            graph[i] += val;
    }
    ///////////////////////////////////////////////////////////////////////////

    IntervalList* findPeakBad( short* badGraph, int tigLen);

    static const bool MATE_3PRIME_END = true;
    FragmentEnds* MateChecker::computeMateCoverage( Unitig* tig, LibraryStats& globalStats )
    {
        int tigLen = tig->getLength();
        short goodGraph[tigLen]; 
        short  badFwdGraph[tigLen];
        short  badRevGraph[tigLen];
        memset( goodGraph, 0, tigLen * sizeof(short));
        memset(  badFwdGraph, 0, tigLen * sizeof(short));
        memset(  badRevGraph, 0, tigLen * sizeof(short));
        IdMap seenMates;
        MateLocation positions(this);
        // Build mate position table
        DoveTailConstIter tigIter = tig->dovetail_path_ptr->begin();
        for(;tigIter != tig->dovetail_path_ptr->end(); tigIter++)
        {
            DoveTailNode frag = *tigIter;
            iuid fragId = frag.ident;
            MateInfo mateInfo = getMateInfo( fragId );
            iuid mateId = mateInfo.mate;
            if ( mateInfo.mate != NULL_FRAG_ID )
            {
                if (positions.hasFrag( mateId ) )
                    positions.addMate( tig->id(), fragId, frag.position );
                else
                    positions.startEntry( tig->id(), fragId, frag.position );
            }
        }
        positions.sort();
        MateLocIter  posIter  = positions.begin();
        for(;        posIter != positions.end(); posIter++) {
            MateLocationEntry loc = *posIter;
            iuid fragId         =  loc.id1;
            MateInfo mateInfo   =  getMateInfo( fragId );
            iuid mateId         =  mateInfo.mate;
            iuid lib            =  mateInfo.lib;
            DistanceCompute *gdc = &(globalStats[ lib ]);
            int badMax = static_cast<int>(gdc->mean + 5 * gdc->stddev);
            int frgBgn = loc.pos1.bgn;
            int frgEnd = loc.pos1.end;
            if ( loc.unitig1 != loc.unitig2) {
                // mate in another tig, mark bad only if max range exceeded
                if (isReverse(loc.pos1)) {
                    if ( frgBgn > badMax ) {
                        if (MATE_3PRIME_END) frgEnd = loc.pos1.bgn;
                        incrRange( badRevGraph, -1, frgBgn - badMax, frgEnd );
                        posIter->isBad = true;
                        fprintf(stderr,"Bad mate %ld pos %ld %ld mate %ld lib %d\n",
                                fragId, frgBgn, loc.pos1.end, mateId, lib);
                    }
                } else {
                    if ( frgBgn + badMax < tigLen ) {
                        iuid beg = frgBgn;
                        if (MATE_3PRIME_END) beg = frgEnd;
                        incrRange( badFwdGraph, -1, beg, frgBgn + badMax );
                        posIter->isBad = true;
                        fprintf(stderr,"Bad mate %ld pos %ld %ld mate %ld lib %d\n",
                                fragId, frgBgn, frgEnd, mateId, lib);
                    }
                }
            } else {
                // both mates in this unitig
                int mateBgn =  loc.pos2.bgn;
                int mateEnd =  loc.pos2.end;
                if (isReverse( loc.pos1 )) {
                    // 1st reversed, so bad 
                    iuid beg = MAX( 0, frgBgn - badMax );
                    incrRange( badRevGraph, -1, beg, frgEnd);
                    posIter->isBad = true;
                    fprintf(stderr,"Bad mate %ld pos %ld %ld mate %ld lib %d\n",
                                fragId, frgBgn, frgEnd, mateId, lib);
                    iuid end;
                    if (isReverse( loc.pos2 )) {
                        // 2nd mate is reversed, so mark bad towards tig begin
                        beg = MAX( 0, mateBgn - badMax );
                        end = mateBgn;
                        if (MATE_3PRIME_END)
                            end = mateEnd;
                        incrRange( badRevGraph, -1, beg, end);
                    } else {
                        // 2nd mate is forward, so mark bad towards tig end
                        end = MIN( tigLen-1, mateBgn + badMax );
                        beg = mateBgn;
                        if (MATE_3PRIME_END)
                            beg = mateEnd;
                        incrRange( badFwdGraph, -1, beg, end);
                    }
                    fprintf(stderr,"Bad mate %ld pos %ld %ld mate %ld lib %d\n",
                                mateId, mateBgn, mateEnd, fragId, lib);

                } else {
                    // 1st forward
                    if (isReverse( loc.pos2 )) {
                        // 2nd reverse so good orient, check distance
                        uint16 mateLen = mateBgn - mateEnd;
                        int mateDist = mateBgn - frgBgn;  

                        int badMin = static_cast<int>(gdc->mean - 5 * gdc->stddev);
                        if (mateDist >= badMin && mateDist <= badMax)
                            incrRange(goodGraph,2, frgBgn, mateEnd);
                        else {
                            // both are bad, mate points towards tig begin
                            iuid beg = MAX(0, mateBgn - badMax); 
                            iuid end = mateBgn;
                            if (MATE_3PRIME_END)
                                end = mateEnd;

                            incrRange(badRevGraph, -1, beg, end);
                            posIter->isBad = true;
                            fprintf(stderr,"Bad mate %ld pos %ld %ld mate %ld lib %d\n",
                                    mateId, mateBgn, mateEnd, fragId, lib);

                            end = MIN( tigLen-1, frgBgn + badMax );
                            beg = frgBgn;
                            if (MATE_3PRIME_END)
                                beg = frgEnd;

                            incrRange(badFwdGraph,-1, beg, end);
                            fprintf(stderr,"Bad mate %ld pos %ld %ld mate %ld lib %d\n",
                                    fragId, frgBgn, frgEnd, mateId, lib);
                        }
                    } else {
                        // 1st and 2nd forward so both bad 
                        iuid end = MIN( tigLen-1, frgBgn + badMax );
                        iuid beg = frgBgn;
                        if (MATE_3PRIME_END)
                            beg = frgEnd;

                        incrRange(badFwdGraph,-1, beg, end);
                        posIter->isBad = true;

                        fprintf(stderr,"Bad mate %ld pos %ld %ld mate %ld lib %d\n",
                                fragId, frgBgn, frgEnd, mateId, lib);

                        // 2nd mate is forward, so mark bad towards tig end
                        end = MIN( tigLen-1, mateBgn + badMax );
                        beg = mateBgn;
                        if (MATE_3PRIME_END)
                            beg = mateEnd;

                        incrRange( badFwdGraph, -1, beg, end);
                        fprintf(stderr,"Bad mate %ld pos %ld %ld mate %ld lib %d\n",
                                mateId, mateBgn, mateEnd, fragId, lib);
                    }
                }
            }
        }
        for(posIter = positions.begin(); posIter != positions.end(); posIter++){
            MateLocationEntry loc = *posIter;
            std::cerr << loc << std::endl;
        }
        // do something with the good and bad graphs
        fprintf(stderr,"Per 300 bases good graph unitig %ld size %ld:\n",tig->id(),tigLen);
        long sum = 0;
        for(int i=0; i < tigLen; i++) {
            if (i > 1 && i % 300 == 0) {
                fprintf(stderr,"%d ", sum / 300);
                sum = 0;
            }
            sum += goodGraph[ i ];
        }
        IntervalList *fwdBads,*revBads;
        fprintf(stderr,"\nPer 300 bases bad fwd graph:\n");
        fwdBads = findPeakBad( badFwdGraph, tigLen );
        fprintf(stderr,"\nPer 300 bases bad rev graph:\n");
        revBads = findPeakBad( badRevGraph, tigLen );

        if (!fwdBads->empty()) {
            fprintf(stderr,"Non empty fwdBads, should handle same as revBads\n");
        }
        FragmentEnds* breaks = new FragmentEnds(); // return value
        posIter  = positions.begin();
        fprintf(stderr,"Num revBads is %d\n",revBads->size());
        IntervalList::const_iterator badIter = revBads->begin();
        bool inBad = false;
        for(; badIter != revBads->end(); badIter++) {
            SeqInterval bad = *badIter;
            fprintf(stderr,"Bad peak from %d to %d\n",bad.bgn,bad.end);
            for(;        posIter != positions.end(); posIter++) {
                MateLocationEntry loc = *posIter;
                if ( !inBad && bad < loc.pos1 ) {
                    fprintf(stderr,"First frg in peak bad range is %d\n", loc.id1);
                    inBad = true;
                }
                if ( inBad && loc.isBad ) {
                    fragment_end_type fragEndInTig = FIVE_PRIME;
                    if (isReverse( loc.pos1 ))
                        fragEndInTig = THREE_PRIME;
                    UnitigBreakPoint bp( loc.id1, fragEndInTig );
                    bp.position = loc.pos1;
                    bp.inSize = 100000;
                    bp.inFrags = 10;
                    breaks->push_back( bp );
                    fprintf(stderr,"First bad frg in peak bad range is %d\n", loc.id1);
                    inBad = false;
                    break;
                }
            }
        }
        delete fwdBads;
        delete revBads;
        return breaks;
    }

    ///////////////////////////////////////////////////////////////////////////

    IntervalList* findPeakBad( short* badGraph, int tigLen) {
        long sum = 0;
        for(int i=0; i < tigLen; i++) {
            if (i > 1 && i % 300 == 0) {
                fprintf(stderr,"%d ", sum / 300);
                sum = 0;
            }
            sum += badGraph[ i ];
        }
        fprintf(stderr,"\n");
        IntervalList* peakBads = new IntervalList();
        SeqInterval   peak = NULL_MATE_LOC;
        int badBegin, peakBad, lastBad;
        peakBad = lastBad = badBegin = 0;
        for(int i=0; i < tigLen; i++) {
            if( badGraph[ i ] < -3 ) {
                if (badBegin == 0)  // start bad region
                    badBegin = i;
                if(badGraph[i] < peakBad) {
                    peakBad   = badGraph[i];
                    peak.bgn = i;
                } else if (lastBad < 0 && lastBad == peakBad) {
                    peak.end = i-1;
                }
                lastBad = badGraph[i];
            } else {
                if (badBegin > 0) {  // end bad region
                    fprintf(stderr,"Bad mates >3 from %d to %d peak %d from %d to %d\n",
                            badBegin,i-1,peakBad,peak.bgn,peak.end);
                    peakBads->push_back( peak );
                    peakBad = lastBad = badBegin = 0;
                    peak = NULL_MATE_LOC;
                }
            }
        }
        return peakBads;
    }

    ///////////////////////////////////////////////////////////////////////////

    bool MateLocation::startEntry(iuid unitigID, iuid fragID, SeqInterval fragPos)
    {
        if ( _iidIndex.find( fragID) != _iidIndex.end() )
            return false; // Entry already exists, can't start new

        assert( fragID != NULL_FRAG_ID );
        MateLocationEntry entry;
        entry.id1     = fragID;
        entry.pos1    = fragPos;
        entry.unitig1 = unitigID; 
        entry.id2     = NULL_FRAG_ID;
        entry.pos2    = NULL_MATE_LOC;
        entry.unitig2 = NULL_FRAG_ID;; 
        entry.isBad   = false;

        _table.push_back( entry );
        _iidIndex[ fragID ] = _table.size()-1;
        return true;
    }

    ///////////////////////////////////////////////////////////////////////////

    bool MateLocation::addMate(iuid unitigId, iuid fragId, SeqInterval fragPos)
    {
        iuid mateId = _checker->getMateInfo( fragId ).mate;
        IdMapConstIter entryIndex = _iidIndex.find( mateId );
        if ( _iidIndex.find( fragId ) != _iidIndex.end() ||
                           entryIndex == _iidIndex.end() )
            return false; // Missing mate or already added

        iuid idx = entryIndex->second;
        _table[ idx ].id2      = fragId;
        _table[ idx ].pos2     = fragPos;
        _table[ idx ].unitig2  = unitigId;
        _iidIndex[ fragId ]    = idx;
        return true;
    }

    ///////////////////////////////////////////////////////////////////////////

    MateLocationEntry MateLocation::getById( iuid fragId ) {
        IdMapConstIter entryIndex = _iidIndex.find( fragId );
        if ( entryIndex == _iidIndex.end() )
            return NULL_MATE_ENTRY;
        else
            return _table[ entryIndex->second ];
    }

    ///////////////////////////////////////////////////////////////////////////

    bool MateLocation::hasFrag(iuid fragId)
    {
        if ( _iidIndex.find( fragId ) == _iidIndex.end() )
            return false;
        else
            return true;
    }

    ///////////////////////////////////////////////////////////////////////////

    void MateLocation::sort()
    {
        std::sort(begin(),end());
        MateLocCIter iter = begin();
        int i = 0;
        for(; iter != end(); iter++, i++) {
            MateLocationEntry entry = *iter;
            _iidIndex[ entry.id1 ] = i;
            _iidIndex[ entry.id2 ] = i;
        }
    }

    ///////////////////////////////////////////////////////////////////////////

    std::ostream& operator<< (std::ostream& os, MateLocationEntry& e)
    {
        int dist = e.pos2.bgn - e.pos1.bgn;
        os << e.pos1.bgn <<","<< e.pos1.end <<" -> "<< e.pos2.bgn <<","<< e.pos2.end
           << " "<< dist << " Frg " << e.id1 << " mate " << e.id2 << " isBad " << e.isBad ;
        return os;
    }

    ///////////////////////////////////////////////////////////////////////////
}
