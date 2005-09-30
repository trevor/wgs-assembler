
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
/*************************************************
* Module:  AS_BOG_BestOverlapGraph.hh
* Description:
*        Data structure to contain the best overlaps and containments
*        based on a defined metric.
* 
*    Programmer:  K. Li
*       Started:  20 July 2005
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: AS_BOG_BestOverlapGraph.hh,v 1.16 2005-09-30 14:17:31 eliv Exp $
 * $Revision: 1.16 $
*/

//  System include files

#ifndef INCLUDE_AS_BOG_BESTOVERLAPGRAPH
#define INCLUDE_AS_BOG_BESTOVERLAPGRAPH

#include <map>
#include "AS_BOG_Datatypes.hh"

extern "C" {
#include "OlapStoreOVL.h"
}

namespace AS_BOG{

    ///////////////////////////////////////////////////////////////////////

    struct BestEdgeOverlap{
        iuid frag_b_id;
        fragment_end_type bend;                
        int in_degree;
        float score;
    };

    struct BestFragmentOverlap{

        // Contains information on what a known fragment overlaps.
        // It is assumed that an index into an array of BestOverlap
        // will tell us what fragment has this best overlap

        BestEdgeOverlap five_prime; 
        BestEdgeOverlap three_prime;
    };

    ///////////////////////////////////////////////////////////////////////

    struct BestContainment{

        // Contains what kind of containment relationship exists between
        // fragment a and fragment b

        iuid container;
        float score;
        bool sameOrientation;
        int a_hang;
        int b_hang;
    };

    ///////////////////////////////////////////////////////////////////////

    struct BestOverlapGraph {

            // Constructor, parametrizing maximum number of overlaps
            BestOverlapGraph(int max_fragment_count);

            // Destructor
            ~BestOverlapGraph();

            // Interface to graph visitor
            // accept(BestOverlapGraphVisitor bog_vis){
            // bog_vis.visit(this);
            // }

            // Accessor Functions
            BestEdgeOverlap *getBestEdge(iuid frag_id, fragment_end_type which_end);
            void setBestEdge(const Long_Olap_Data_t& olap, float newScore);
            iuid getNumFragments() { return _num_fragments; }
            BestContainment *getBestContainer(iuid frag_id);

            // Graph building methods
            fragment_end_type AEnd(const Long_Olap_Data_t& olap);
            fragment_end_type BEnd(const Long_Olap_Data_t& olap);
            void scoreOverlap(const Long_Olap_Data_t& olap);
            short olapLength(const Long_Olap_Data_t& olap);
            bool checkForNextFrag(const Long_Olap_Data_t& olap);
            virtual float score( const Long_Olap_Data_t& olap) =0;
            void transitiveContainment();

            // FragStore related variables
        //These should be moved to protected
            static uint16 *fragLength;
            static uint16 fragLen( iuid );
            static ReadStructp fsread;
            static FragStoreHandle fragStoreHandle;
            static overlap_type getType(const Long_Olap_Data_t & olap);

            std::map<iuid, BestContainment> _best_containments;

        protected:
            iuid _num_fragments;
            iuid curFrag;
            int bestLength;
            BestFragmentOverlap* _best_overlaps;

    }; //BestOverlapGraph

    ///////////////////////////////////////////////////////////////////////////

    struct ErateScore : public BestOverlapGraph {
        ErateScore(int num) : BestOverlapGraph(num) {}
        float score( const Long_Olap_Data_t& olap);
    };

    struct LongestEdge : public BestOverlapGraph {
        LongestEdge(int num) : BestOverlapGraph(num) {}
        float score( const Long_Olap_Data_t& olap);
    };

    struct LongestHighIdent : public BestOverlapGraph {
        float mismatchCutoff;
        LongestHighIdent(int num, float maxMismatch)
            : BestOverlapGraph(num), mismatchCutoff(maxMismatch) {}
        float score( const Long_Olap_Data_t& olap);
    };

    ///////////////////////////////////////////////////////////////////////////

} //AS_BOG namespace


#endif //INCLUDE_AS_BOG_BESTOVERLAPGRAPH

