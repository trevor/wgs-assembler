
/**************************************************************************
 * Copyright (C) 2005, J Craig Venter Institute. All rights reserved.
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
 * Module:  FindContainedfragments.c
 * Description:
 *   Dump (to stdout) the contents of an overlap store
 * 
 *    Programmer:  E. Venter
 *       Started:  11 July 2005
 * 
 * Assumptions:
 * 
 * Notes:
 *
 *************************************************/

/* RCS info
 * $Id: FindContainedFrags.cc,v 1.32 2007-12-27 18:41:15 brianwalenz Exp $
 * $Revision: 1.32 $
 */

static const char CM_ID[] = "$Id: FindContainedFrags.cc,v 1.32 2007-12-27 18:41:15 brianwalenz Exp $";

//  System include files

#include<map>
#include<set>
#include<vector>
#include<cmath>   // for ceil and log10
#include<cstdlib> // for abs(int)
#include<iostream>

#include"AS_BOG_BestOverlapGraph.hh"

using std::cout;
using std::endl;
using std::map;
using std::set;
using AS_BOG::BestOverlapGraph;

//  Local include files
extern "C" {
#include "AS_OVS_overlapStore.h"
}

int  main
(int argc, char * argv [])

{
    OverlapStore  * my_store;

    // Get path/names of olap and frg stores from command line
    const char* OVL_Store_Path = argv[1];
    const char* FRG_Store_Path = argv[2];

    argc = AS_configure(argc, argv);

    // Open and initialize Overlap store
    my_store = AS_OVS_openOverlapStore(OVL_Store_Path);

    uint32  last = AS_OVS_lastFragInStore(my_store);

    // must be before creating scores, since it sets the num frags
    AS_BOG::BOG_Runner bogRunner(last);

    // Initialize our three different types of Best Overlap Graphs
    //   AS_BOG::ErateScore erScore;
    //   AS_BOG::LongestEdge lenScore;
    AS_BOG::BestOverlapGraph lenid25(0.025);
    AS_BOG::BestOverlapGraph lenid15(0.015);
    AS_BOG::BestOverlapGraph lenid10(0.010);

    // Put the three graphs into a vector, so we can step through them
    //   bogRunner.push_back(&erScore);
    //   bogRunner.push_back(&lenScore);
    bogRunner.push_back(&lenid25);
    bogRunner.push_back(&lenid15);
    bogRunner.push_back(&lenid10);

    // Go through the overlap stream, and populate the 3 overlap graphs
    bogRunner.processOverlapStream( my_store, FRG_Store_Path );

    // Free/clean up the frag/overlap store/stream handles
    AS_OVS_closeOverlapStore(my_store);

    // Compute the width (number of digits) in the largest IUID
    int pad = static_cast<int>(ceil( log10( last )));

    // For each IUID
    for(int i = 1; i <= last; i++) {
        AS_BOG::BestEdgeOverlap* five[bogRunner.metrics.size()];
        AS_BOG::BestEdgeOverlap* three[bogRunner.metrics.size()];

        // For each graph type
        for(int j = 0; j < bogRunner.metrics.size(); j++) { // output olap graph

            // Retrieve the best overlaps from the BOG for the current IUID
            five[j] = bogRunner.metrics[j]->getBestEdgeOverlap( i, AS_BOG::FIVE_PRIME );
            three[j] = bogRunner.metrics[j]->getBestEdgeOverlap( i, AS_BOG::THREE_PRIME );

            // Why isn't this just outside of the for j loop?
            if ( j == bogRunner.metrics.size()-1 ) {

                // Get the 5' overlap for the current IUID
                AS_IID    b0 = five[0]->frag_b_id;
                AS_BOG::BestEdgeOverlap* f1 = five[1];
                AS_BOG::BestEdgeOverlap* f2 = five[2];

                // Set up output formatting
                cout.flags(std::ios_base::left);
                cout.width(pad);

                // Output the fragments 5' ends in_degree, for each graph
                cout << i <<" 5' "<< five[0]->in_degree << " "
                     << f1->in_degree <<" "<< f2->in_degree;

                // Report whether best overlaps agree between graphs
                if (b0 == five[1]->frag_b_id && b0 == five[2]->frag_b_id) {
                    cout <<" best "; cout.width(pad * 3 + 2);
                    cout << b0 ;
                } else { 
                    cout<< " diff "          ; cout.width(pad);
                    cout<< b0 <<" "          ; cout.width(pad);
                    cout<< f1->frag_b_id<<" "; cout.width(pad);
                    cout<< f2->frag_b_id ;
                }

                // Get the 3' overlap for the current IUID
                b0 = three[0]->frag_b_id;
                f1 = three[1];
                f2 = three[2];

                // Output the fragments 3' ends in_degree, for each graph
                cout << " 3' "<< three[0]->in_degree << " "
                     << f1->in_degree <<" "<< f2->in_degree;

                // Report whether best overlaps agree between graphs
                if (b0 == three[1]->frag_b_id && b0 == three[2]->frag_b_id) {
                    cout <<" best "; cout.width(pad);
                    cout << b0 << endl;
                } else { 
                    cout<<" diff ";            cout.width(pad);
                    cout<< b0 <<" ";           cout.width(pad);
                    cout<< f1->frag_b_id<<" "; cout.width(pad);
                    cout<< f2->frag_b_id << endl;
                }

            }
        } // end for each metric
    } // end for each fragment

    // We should typedef this map for best_containments
    map<AS_IID,AS_BOG::BestContainment> c1 = bogRunner.metrics[0]->_best_containments;

    // Iterate through all the containees, this may miss containees that exists
    //   by other bogRunner.metrics/graph types
    for(map<AS_IID,AS_BOG::BestContainment>::const_iterator it = c1.begin();
        it != c1.end(); it++) {

        // First is containee IUID, Second is container IUID
        AS_IID    id = it->first;
        AS_BOG::BestContainment bst = it->second;
        AS_IID    cr1 = bst.container;
        int ahng1 = bst.a_hang;


	//  1: erScore
        cout << id << " c1 by "<< cr1 << " " << bst.contain_score << " " << ahng1 <<
            " " << bst.b_hang <<" sameOrient " << bst.sameOrientation << endl;

        //  2: lenScore
        AS_BOG::BestContainment* b2 = bogRunner.metrics[1]->getBestContainer(id);

       if ( b2 != NULL ) {
           AS_IID    cr2 = b2->container;
           int ahng2 = b2->a_hang;

           cout << id << " c2 by "<< cr2 << " " << b2->contain_score << " " << ahng2 <<
               " " << b2->b_hang <<" sameOrient " << b2->sameOrientation << endl;

           if ( cr1 == cr2 && ahng1 != ahng2 )
              cout<<"Error, hang mismatch: "<<cr1<<" "<<ahng1<<" "<<ahng2<<endl;
           cr1   = cr2;
           ahng1 = ahng2;
       }
       
        //  3: lenIdent
        AS_BOG::BestContainment* b3 = bogRunner.metrics[2]->getBestContainer(id);
        if ( b3 != NULL ) {
            AS_IID    cr3 = b3->container;
            int ahng3 = b3->a_hang;

            cout << id << " c3 by "<< cr3 << " " << b3->contain_score << " " << ahng3 <<
                " " << b3->b_hang << " sameOrient " << b3->sameOrientation << endl;

            if ( cr3 == cr1 && ahng3 != ahng1 )
                cout<<"Error, hang mismatch: "<<cr3<<" "<<ahng3<<" "<<ahng1<<endl;
        }

    }

    // Shouldn't these both be n a  destructor in BOG?
    delete[] BestOverlapGraph::fragLength;

    return  0;
}
