
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
static char *rcsid = "$Id: ChunkOverlap_CGW.c,v 1.37 2009-01-08 14:39:46 brianwalenz Exp $";

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "AS_UTL_Hash.h"
#include "AS_UTL_Var.h"
#include "AS_CGW_dataTypes.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"    // For DeleteCIOverlapEdge
#include "AS_ALN_aligners.h"
#include "CommonREZ.h"
#include "UtilsREZ.h"

// this is the initial range we use to compute
// overlaps with a different min/max range
#define BEGENDSLOP 10

// This will enable the use of Local_Overlap_AS_forCNS as a fallback
// when DP_Compare fails.  It will slow down cgw by over an order of
// magnitude (cgw on some microbes was taking more than 12 hours).
//
// In p.cnpt3, enabling this found 3 more overlaps.
//
//      3 only found with Local_Overlap
//     20 only found with DP_Compare
//    618 found by both
//   2711 found by none
//   3352 tries
//
#undef USE_LOCAL_OVERLAP_AS_FALLBACK



/* ChunkOverlap_CGW provides tools for invoking Gene's dpalign tool to compute
   overlaps between chunks.  Such overlaps are first 'collected' and then 'computed'.
   Then, they may be 'looked up'.

   This sparse 'database' of relationships between chunks (both overlaps and lack of overlaps) is
   stored in a symbol table based on AS_UTL_Hash to facilitate quick lookups and insertions.
   Key to this database is a canonical overlap representation.  The canonical storage for a
   potential overlap is as follows:
   if(orientation is symettric -- AB_BA or BA_AB )
   then the overlap is stored with the lower numbered chunk first
   if(orientation is antinormal -- BA_BA)
   then the overlap is stored as an normal overlap with the chunks reversed

   The idea of the collection phase, is to determine the maximal overlap range between
   two chunks that is worth computing.  So, during the construction of the raw
   mate edges for the extended chunk graph, any mate link that implies a possible overlap is
   collected for later evaluation.  If multiple overlaps are collected for the same chunk pair
   and orientation, the maximal overlap interval is collected.

   Once a set of inter-chunk relationships have been collected, Gene's dpalign tool is invoked
   on the consensus sequence for each chunk, and the results are stored in the database.  As
   an option, these overlap edges can also be added to the extended chunk graph.  The implied
   overlaps that are detected are thus added to the set of raw edge mates prior to merging.

   Once the extended chunk graph has been merged, we look for potential overlaps between
   unique chunks that are implied by the confirmed edges are checked.  These overlaps are NOT
   added tot he extended chunk graph, but are simply stored in the database for later
   retrieval in the scaffold construction code.

   Saul Kravitz
   April 1999


   Additions made by Knut Reinert
   09/10/99

   We make use of the above hash table to compute and store quality values with the overlaps
   Initially no meaningful quality value is stored in the ChunkOverlapCheckT::overlap
   We add bits to the ChunkOverlapCheckT struct to indicate whether there was a quality
   computation, and which function computed the quality value.
   Whenever a new method for computing the quality of an overlap should be tested,
   one should augment ChunkOverlapCheckT by the appropriate bit and add a value
   to the enum QualityFuncsT.

*/





//  Hash table support functions
//
static
int CanOlapCmp(uint64 cO1, uint64 cO2){
  ChunkOverlapSpecT *c1 = (ChunkOverlapSpecT *)(INTPTR)cO1;
  ChunkOverlapSpecT *c2 = (ChunkOverlapSpecT *)(INTPTR)cO2;
  int diff;

  diff = c1->cidA - c2->cidA;
  if(diff)
    return diff;

  diff = c1->cidB - c2->cidB;
  if(diff)
    return diff;

  if(c1->orientation == c2->orientation)
    return 0;
  else return (c1 - c2);
}

static
int CanOlapHash(uint64 cO, uint32 length){
  return  Hash_AS((uint8 *)(INTPTR)cO, length, 37);
}


//external
ChunkOverlapperT *CreateChunkOverlapper(void){
  ChunkOverlapperT *chunkOverlapper = (ChunkOverlapperT *)safe_malloc(sizeof(ChunkOverlapperT));
  chunkOverlapper->hashTable = CreateGenericHashTable_AS(CanOlapHash, CanOlapCmp);
  chunkOverlapper->ChunkOverlaps = AllocateHeap_AS(sizeof(ChunkOverlapCheckT));
  return chunkOverlapper;
}

//external
void DestroyChunkOverlapper(ChunkOverlapperT *chunkOverlapper){
  DeleteHashTable_AS(chunkOverlapper->hashTable);
  FreeHeap_AS(chunkOverlapper->ChunkOverlaps);
}


//external
int InsertChunkOverlap(ChunkOverlapperT *chunkOverlapper,
                       ChunkOverlapCheckT *olap){
  ChunkOverlapCheckT *nolap = (ChunkOverlapCheckT *)GetHeapItem_AS(chunkOverlapper->ChunkOverlaps);
  *nolap = *olap;
  assert(nolap->overlap==0||nolap->errorRate >= 0.0);
  return(InsertInHashTable_AS(chunkOverlapper->hashTable,
                              (uint64)(INTPTR)&nolap->spec,
                              sizeof(olap->spec),
                              (uint64)(INTPTR)nolap,
                              0));
}


//external
void  SaveChunkOverlapperToStream(ChunkOverlapperT *chunkOverlapper, FILE *stream){

  /* We will save the ChunkOverlaps to a file, and rebuild the heap after
     we read the file.  In other words, we don't worry about hashTable persistence. */

  HashTable_Iterator_AS iterator;
  uint64 key, value;
  uint32 valuetype;
  int64 numOverlaps = 0;

  // Iterate over all hashtable elements, just to count them

  InitializeHashTable_Iterator_AS(chunkOverlapper->hashTable, &iterator);

  while(NextHashTable_Iterator_AS(&iterator, &key, &value, &valuetype)){
    numOverlaps++;
  }

  AS_UTL_safeWrite(stream, &numOverlaps, "SaveChunkOverlapperToStream", sizeof(int64), 1);

  // Iterate over all hashtable elements, writing

  InitializeHashTable_Iterator_AS(chunkOverlapper->hashTable, &iterator);

  while(NextHashTable_Iterator_AS(&iterator, &key, &value, &valuetype)){
    ChunkOverlapCheckT *olap = (ChunkOverlapCheckT*)(INTPTR)value;

    AS_UTL_safeWrite(stream, olap, "SaveChunkOverlapperToStream", sizeof(ChunkOverlapCheckT), 1);
  }
}


//external
ChunkOverlapperT *  LoadChunkOverlapperFromStream(FILE *stream){
  int64 numOverlaps;
  int status;
  int64 overlap;
  ChunkOverlapCheckT olap = {0};
  ChunkOverlapperT *chunkOverlapper;

  // open the chunkStore at chunkStorepath
  status = AS_UTL_safeRead(stream, &numOverlaps, "LoadChunkOverlapperFromStream", sizeof(int64), 1);
  assert(status == 1);

  chunkOverlapper = (ChunkOverlapperT *)safe_malloc(sizeof(ChunkOverlapperT));
  chunkOverlapper->hashTable = CreateGenericHashTable_AS(CanOlapHash, CanOlapCmp);
  chunkOverlapper->ChunkOverlaps = AllocateHeap_AS(sizeof(ChunkOverlapCheckT));

  for(overlap = 0; overlap < numOverlaps; overlap++){
    status = AS_UTL_safeRead(stream, &olap, "LoadChunkOverlapperFromStream", sizeof(ChunkOverlapCheckT), 1);
    assert(status == 1);
    assert(olap.errorRate > 0.0);
    InsertChunkOverlap(chunkOverlapper, &olap);
  }

  return chunkOverlapper;
}




/************************************************************************
 * A canonical overlap hs the following properties
 *       if(orientation is symmetric -- AB_BA or BA_AB )
 *       	then the overlap is stored with the lower numbered chunk first
 *	if(orientation is antinormal -- BA_BA)
 *	        then the overlap is stored as an normal overlap with the chunks reversed
 *************************************************************************/

//external
int InitCanonicalOverlapSpec(CDS_CID_t cidA, CDS_CID_t cidB,
                             ChunkOrientationType orientation,
                             ChunkOverlapSpecT *spec){
  int canonical = TRUE;

  switch(orientation){

    case BA_BA:
      spec->orientation = AB_AB;
      spec->cidA = cidB;
      spec->cidB = cidA;
      canonical = FALSE;
      break;
    case AB_AB:
      spec->orientation = orientation;
      spec->cidA = cidA;
      spec->cidB = cidB;
      break;
    case BA_AB:
    case AB_BA:
      spec->orientation = orientation;
      if(cidA > cidB){
	spec->cidA = cidB;
	spec->cidB = cidA;
	canonical = FALSE;
      }else{
	spec->cidA = cidA;
	spec->cidB = cidB;
      }
      break;
    default:
      assert(0);
      break;
  }

  return canonical;
}



/* Given a graph edge, create an overlap in the hashtable and mark it as computed */
//external
void CreateChunkOverlapFromEdge(GraphCGW_T *graph,
                                EdgeCGW_T *edge,
                                int bayesian){
  ChunkOverlapCheckT olap = {0};
  double delta = sqrt(edge->distance.variance) * 3.0;
  assert((0.0 <= AS_CGW_ERROR_RATE) && (AS_CGW_ERROR_RATE <= AS_MAX_ERROR_RATE));
  InitCanonicalOverlapSpec(edge->idA, edge->idB, edge->orient, &olap.spec);
  olap.computed = TRUE;
  olap.overlap = -edge->distance.mean;
  olap.minOverlap = -edge->distance.mean - delta;
  olap.maxOverlap = -edge->distance.mean + delta;;
  olap.fromCGB = FALSE;
  olap.cgbMinOverlap = 0;
  olap.cgbMaxOverlap = 0;
  olap.hasBayesianQuality = bayesian;
  olap.errorRate = AS_CGW_ERROR_RATE;
  olap.quality = edge->quality;
  olap.ahg = olap.bhg = 0;
  olap.min_offset = olap.max_offset = 0;
  InsertChunkOverlap(graph->overlapper, &olap);
}



/* Given a graph edge, create an overlap in the hashtable */
//external
void FillChunkOverlapWithEdge(EdgeCGW_T *edge, ChunkOverlapCheckT *olap){
  double delta = sqrt(edge->distance.variance) * 3.0;
  assert((0.0 <= AS_CGW_ERROR_RATE) && (AS_CGW_ERROR_RATE <= AS_MAX_ERROR_RATE));
  InitCanonicalOverlapSpec(edge->idA, edge->idB, edge->orient, &olap->spec);
  olap->computed = FALSE;
  olap->overlap = -edge->distance.mean;

  // might be unsafe for big variances after tandem mark propagation
  olap->minOverlap = (CDS_COORD_t) -edge->distance.mean - delta;
  olap->maxOverlap = (CDS_COORD_t) -edge->distance.mean + delta;

  olap->minOverlap = (CDS_COORD_t) -edge->distance.mean - 20;
  olap->maxOverlap = (CDS_COORD_t) -edge->distance.mean + 20;
  olap->fromCGB = FALSE;
  olap->cgbMinOverlap = 0;
  olap->cgbMaxOverlap = 0;
  olap->hasBayesianQuality = FALSE;
  olap->errorRate = AS_CGW_ERROR_RATE;
  olap->quality = edge->quality;
  olap->ahg = olap->bhg = 0;
  olap->min_offset = olap->max_offset = 0;
}



/* Given a graph edge, create an overlap in the hashtable */
static
void FillChunkOverlapWithUOM(ChunkOverlapCheckT *olap, UnitigOverlapMesg *uom_mesg){
  assert((0.0 <= AS_CGW_ERROR_RATE) && (AS_CGW_ERROR_RATE <= AS_MAX_ERROR_RATE));
  InitCanonicalOverlapSpec(uom_mesg->chunk1, uom_mesg->chunk2, uom_mesg->orient, &olap->spec);
  olap->computed = FALSE;
  olap->overlap = uom_mesg->best_overlap_length;
  olap->minOverlap = uom_mesg->min_overlap_length;
  olap->maxOverlap = uom_mesg->max_overlap_length;
  olap->hasBayesianQuality = !ScaffoldGraph->alignOverlaps;
  olap->errorRate = AS_CGW_ERROR_RATE;
  olap->quality = uom_mesg->quality;
  olap->ahg = olap->bhg = 0;
  olap->min_offset = olap->max_offset = 0;
}


//external
ChunkOverlapCheckT *LookupCanonicalOverlap(ChunkOverlapperT *chunkOverlapper,
                                           ChunkOverlapSpecT *spec){
  return (ChunkOverlapCheckT *)(INTPTR)LookupValueInHashTable_AS(chunkOverlapper->hashTable, (uint64)(INTPTR)spec, sizeof(*spec));
}


// Returns FALSE if none found
// Returns TRUE and overwrites *olap if found
//
//external
int LookupOverlap(GraphCGW_T *graph,
		  CDS_CID_t cidA, CDS_CID_t cidB,
		  ChunkOrientationType orientation,
		  ChunkOverlapCheckT *olap){
  ChunkOverlapperT *chunkOverlapper = graph->overlapper;
  int isCanonical ;
  ChunkOverlapSpecT spec;
  ChunkOverlapCheckT *lookup;
  isCanonical = InitCanonicalOverlapSpec(cidA, cidB, orientation, &spec);
  lookup = LookupCanonicalOverlap(chunkOverlapper, &spec);
  if(!lookup)  // We didn't find anything
    return FALSE;

  *olap = *lookup;

  if(isCanonical){  // we're done
    return TRUE;
  }
  // Otherwise we have to fix up the retrieved canonical overlap for the non-canonical query
  //
  olap->spec.orientation = orientation;
  olap->spec.cidA = cidA;
  olap->spec.cidB = cidB;
  if(olap->BContainsA || olap->AContainsB)
    {
      int swap;
      NodeCGW_T *a = GetGraphNode(graph, cidA);
      NodeCGW_T *b = GetGraphNode(graph, cidB);

      swap = olap->BContainsA;
      olap->BContainsA = olap->AContainsB;
      olap->AContainsB = swap;
      /// NEW!!!!
      if(olap->AContainsB){

	switch(orientation){
          case AB_AB:
          case AB_BA:
            olap->overlap = b->bpLength.mean + olap->bhg;
            break;
          case BA_AB:
          case BA_BA:
            olap->overlap = b->bpLength.mean - olap->ahg;
            break;
          default:
            assert(0);
            break;
	}
      }else if(olap->BContainsA){
	switch(orientation){
          case AB_AB:
          case AB_BA:
            olap->overlap = a->bpLength.mean - olap->bhg;
            break;
          case BA_AB:
          case BA_BA:
            olap->overlap = a->bpLength.mean + olap->ahg;
            break;
          default:
            assert(0);
            break;
	}
      }
      /// END NEW!
    }

  return TRUE;

}




static
int DeleteChunkOverlap(ChunkOverlapperT *chunkOverlapper,
                       ChunkOverlapCheckT *olap){
  return DeleteFromHashTable_AS(chunkOverlapper->hashTable, (uint64)(INTPTR)&olap->spec, sizeof(olap->spec));
}



/* Insert a computed overlap as a CIEdgeT into the Scaffold Graph. */

//external
CDS_CID_t InsertComputedOverlapEdge(GraphCGW_T *graph,
                                    ChunkOverlapCheckT *olap){
  CDS_CID_t eid;
  int fudge;
  int isDoveTail = FALSE;
  LengthT overlap;
  ChunkOrientationType orient = olap->spec.orientation;
  EdgeCGW_T *existing = FindGraphOverlapEdge(graph, olap->spec.cidA, olap->spec.cidB, orient);

  overlap.mean   = -olap->overlap;
  overlap.variance = MAX(1.0, ComputeFudgeVariance((double)olap->overlap));
  fudge = sqrt(overlap.variance);

  isDoveTail = !(olap->AContainsB || olap->BContainsA);


  // If there is an existing overlap edge, don't insert this one!
  if(existing){
    double diff = abs(existing->distance.mean - overlap.mean);
    if(diff < 5.0){ // this is the same edge
      CDS_CID_t eid = GetVAIndex_EdgeCGW_T(graph->edges, existing);
      return eid;
    }
  }

  eid = AddGraphEdge(graph,
                     olap->spec.cidA,
                     olap->spec.cidB,
                     NULLINDEX, NULLINDEX, // frags
                     NULLINDEX,  // dist
                     overlap,
                     olap->quality,
                     fudge,
                     orient,
                     FALSE, // inducedByUnknownOrientation
                     isDoveTail,  // isOverlap
                     olap->AContainsB,                 // isAContainsB
                     olap->BContainsA,                 // isBContainsA
                     FALSE,                        // isTransChunk
                     FALSE, FALSE,  // extremalA and extremalB
                     UNKNOWN_EDGE_STATUS,
                     FALSE,
                     TRUE /* insert*/ );

  return eid;
}






//external
void CollectChunkOverlap(GraphCGW_T *graph,
                         CDS_CID_t cidA, CDS_CID_t cidB,
                         ChunkOrientationType orientation,
                         float   meanOverlap, float   deltaOverlap,
                         float   quality, int bayesian,
                         int fromCGB,
			 int verbose){
  ChunkOverlapperT *chunkOverlapper = graph->overlapper;
  ChunkOverlapCheckT canOlap={0}, *olap;
  CDS_COORD_t delta;
  CDS_COORD_t minOverlap,maxOverlap;

  delta = (CDS_COORD_t)(3.0 * deltaOverlap);
  delta = MAX(delta, 10);
  minOverlap = MAX(meanOverlap - delta, 0);
  maxOverlap = meanOverlap + delta;

  if(maxOverlap < 0){
    // There is no chance that these overlap!
    return;
  }

  // Create a canonical representation of the overlap
  InitCanonicalOverlapSpec(cidA, cidB, orientation, &canOlap.spec);

  // Lookup to see if we've already stored such an overlap
  olap = LookupCanonicalOverlap(chunkOverlapper, &canOlap.spec);
  if(!olap){
    assert((0.0 <= AS_CGW_ERROR_RATE) && (AS_CGW_ERROR_RATE <= AS_MAX_ERROR_RATE));
    canOlap.computed = FALSE;
    canOlap.overlap = FALSE;
    canOlap.quality = 1.0;
    canOlap.hasBayesianQuality = FALSE;
    canOlap.minOverlap = minOverlap;
    canOlap.maxOverlap = maxOverlap;
    canOlap.fromCGB = fromCGB;
    if(fromCGB && bayesian){
      canOlap.computed = TRUE;
      canOlap.quality = quality;
      canOlap.hasBayesianQuality = TRUE;
      canOlap.overlap = (canOlap.minOverlap + canOlap.maxOverlap)/2;
    }
    canOlap.cgbMinOverlap = minOverlap;
    canOlap.cgbMaxOverlap = maxOverlap;
    canOlap.errorRate = AS_CGW_ERROR_RATE;

    // Add it to the symbol table
    if(InsertChunkOverlap(chunkOverlapper, &canOlap) != HASH_SUCCESS)
      assert(0);

  }else{ // we found an overlap
    // We found one.  So now we need to update the existing record so that
    // it reflects the maximal interval we're interested in overlapping.
    if(olap->computed){
      // If we've already computed this one, only recompute if the range is expanded
      if(olap->overlap == 0 &&
	 (minOverlap < olap->minOverlap ||
          maxOverlap > olap->maxOverlap)){
	olap->computed = FALSE; // Recompute
	olap->hasBayesianQuality = FALSE;
	olap->minOverlap = MIN(minOverlap, olap->minOverlap);
	olap->maxOverlap = MAX(maxOverlap, olap->maxOverlap);
      }


    }else{
      if(fromCGB){  // a CGB overlap clobbers whatever is there
        if(!olap->fromCGB){ // Not from the chunk graph builder
          olap->cgbMinOverlap = minOverlap;
          olap->cgbMaxOverlap = maxOverlap;
          olap->overlap = (minOverlap + maxOverlap)/2;
          olap->fromCGB = fromCGB;
        }else{              // From the chunk graph builder
          if(quality < olap->quality){
            olap->cgbMinOverlap = minOverlap;
            olap->cgbMaxOverlap = maxOverlap;
          }
        }
	olap->quality = quality;
	olap->hasBayesianQuality = bayesian;
	if(bayesian){
	  olap->overlap = (olap->cgbMinOverlap + olap->cgbMaxOverlap)/2;
	  olap->computed = TRUE;
	}

      }
      olap->minOverlap = MIN(minOverlap, olap->minOverlap);
      olap->maxOverlap = MAX(maxOverlap, olap->maxOverlap);
    }
  }
}






//external
Overlap* OverlapSequences( char *seq1, char *seq2,
                           ChunkOrientationType orientation,
                           CDS_COORD_t min_ahang, CDS_COORD_t max_ahang,
                           double erate, double thresh, CDS_COORD_t minlen)
{
  Overlap *dp_omesg = NULL;
  Overlap *lo_omesg = NULL;
  int flip = 0;

  //  This function takes two sequences, their orientation and an
  //  assumed minimum and maximum ahang for which it checks.

  // if the orientation is BA_AB or BA_BA, we need to reverse
  // complement the first contig
  if (orientation == BA_AB || orientation == BA_BA)
    Complement_Seq( seq1 );

  // if the orientation is AB_BA or BA_BA, we need to set the flip
  // variable for the second contig
  if (orientation == AB_BA || orientation == BA_BA)
    flip = 1;

  // min_ahang and end are essentially bounds on the a-hang

  dp_omesg = DP_Compare(seq1, seq2,
                        min_ahang, max_ahang,
                        strlen(seq1), strlen(seq2),
                        flip,
                        erate, thresh, minlen,
                        AS_FIND_ALIGN);

  if ((dp_omesg != NULL) && (dp_omesg->length <= minlen))
    dp_omesg = NULL;

  //if (dp_omesg != NULL)
  //  fprintf(stderr, "OverlapSequences()-- Found overlap with DP_Compare   begpos=%d endpos=%d length=%d diffs=%d comp=%d/%d\n",
  //          dp_omesg->begpos, dp_omesg->endpos, dp_omesg->length, dp_omesg->diffs, dp_omesg->comp, flip);

#ifdef USE_LOCAL_OVERLAP_AS_FALLBACK
  if (!dp_omesg)
    lo_omesg = Local_Overlap_AS_forCNS(seq1, seq2,
                                       min_ahang, max_ahang,
                                       strlen(seq1), strlen(seq2),
                                       flip,
                                       erate, thresh, minlen,
                                       AS_FIND_LOCAL_OVERLAP);

  if ((lo_omesg != NULL) && (lo_omesg->length <= minlen))
    lo_omesg = NULL;

  //if (lo_omesg != NULL)
  //  fprintf(stderr, "OverlapSequences()-- Found overlap with Local_Overlap   begpos=%d endpos=%d length=%d diffs=%d comp=%d/%d\n",
  //          lo_omesg->begpos, lo_omesg->endpos, lo_omesg->length, lo_omesg->diffs, lo_omesg->comp, flip);
#endif

  if (orientation == BA_AB || orientation == BA_BA)
    Complement_Seq( seq1 );

  // omesg->begpos is the a-hang, omesg->endpos is the b-hang

  return((dp_omesg) ? dp_omesg : lo_omesg);
}



static VA_TYPE(char) *consensusA = NULL;
static VA_TYPE(char) *consensusB = NULL;
static VA_TYPE(char) *qualityA = NULL;
static VA_TYPE(char) *qualityB = NULL;



static
void ComputeCanonicalOverlap_new(GraphCGW_T *graph,
                                 ChunkOverlapCheckT *canOlap)
{
  CDS_COORD_t lengthA, lengthB;
  ChunkOverlapperT *chunkOverlapper = graph->overlapper;
  Overlap * tempOlap1;
  ChunkOverlapSpecT inSpec;

  if(consensusA == NULL)
    {
      consensusA = CreateVA_char(2048);
      consensusB = CreateVA_char(2048);
      qualityA = CreateVA_char(2048);
      qualityB = CreateVA_char(2048);
    }

  // Save the input spec
  inSpec = canOlap->spec;

  canOlap->BContainsA = FALSE;
  canOlap->AContainsB = FALSE;
  canOlap->computed = TRUE;
  canOlap->overlap = FALSE;
  canOlap->ahg = canOlap->bhg = 0;

  if(canOlap->maxOverlap < 0) // no point doing the expensive part if there can be no overlap
    return;


  // Get the consensus sequences for both chunks from the ChunkStore
  lengthA = GetConsensus(graph, canOlap->spec.cidA, consensusA, qualityA);
  lengthB = GetConsensus(graph, canOlap->spec.cidB, consensusB, qualityB);

  if(canOlap->minOverlap > (lengthA+lengthB-CGW_DP_MINLEN)) // no point doing the expensive part if there can be no overlap
    return;

  // Return value is length of chunk sequence/quality
  // overlap 'em
  {
    char *seq1, *seq2;
    CDS_COORD_t min_ahang, max_ahang;
    seq1   = Getchar(consensusA, 0);
    seq2   = Getchar(consensusB, 0);

    min_ahang = lengthA - canOlap->maxOverlap;
    max_ahang = lengthA - canOlap->minOverlap;

    // tempOlap1 is a static down inside of DP_Compare, don't free it
    tempOlap1 =
      OverlapSequences(seq1, seq2, canOlap->spec.orientation,
                       min_ahang, max_ahang,
                       canOlap->errorRate,
                       CGW_DP_THRESH, CGW_DP_MINLEN);

    if (tempOlap1 ) {     // Found one....

      if( tempOlap1->begpos < 0 && tempOlap1->endpos > 0) // ahang is neg and bhang is pos
        canOlap->BContainsA = TRUE;
      else if( tempOlap1->begpos > 0 && tempOlap1->endpos < 0) // ahang is pos and bhang is neg
        canOlap->AContainsB = TRUE;

      //	    Print_Overlap_AS(GlobalData->stderrc,&AFR,&BFR,O);
      canOlap->computed = TRUE;
      canOlap->ahg = tempOlap1->begpos;
      canOlap->bhg = tempOlap1->endpos;

      //  Make the overlap field be the number of bases from the tail of
      //  the A sequence to the beginning of the B sequence
      canOlap -> overlap = tempOlap1 -> length;
      if  (canOlap -> ahg < 0)
        canOlap -> overlap -= canOlap -> ahg;
      if  (canOlap -> bhg < 0)
        canOlap -> overlap -= canOlap -> bhg;

      // fields are no longer used in DP_Compare (as opposed to DP_Compare_AS)
      canOlap->quality = 0.0;
      // canOlap->min_offset = O->min_offset;
      // canOlap->max_offset = O->max_offset;

      // here we set the suspicious flag based on what the overlap is
      // if the sequences have slid (e.g., an AB_AB has become a BA_BA)
      // then we change the orientation and set the suspicious flag

      // not dealing with containments here - they can go in either way?

      if (canOlap->ahg < 0 && canOlap->bhg < 0){

        // Try to delete the overlap from the hashtable.  It may or
        // may not be there.  We will (re)insert it later.  If we
        // didn't do this, the hashtable would be corrupted, since the
        // chains in the buckets are ordered by (ida,idb,orientation),
        // so we can't screw around with these, without removing the
        // entry and reinserting it.

        DeleteChunkOverlap(chunkOverlapper, canOlap);

        canOlap->suspicious = TRUE;
        canOlap -> overlap = tempOlap1 -> length;
        switch  (canOlap -> spec . orientation)
          {
            case  AB_AB :
              // we want to switch to a non-canonical orientation
              // ...canOlap->spec.orientation = BA_BA; but since
              // we can't, we switch order of operands instead
              canOlap -> spec. cidA = inSpec . cidB;
              canOlap -> spec. cidB = inSpec . cidA;
              canOlap -> ahg = - canOlap -> ahg;
              canOlap -> bhg = - canOlap -> bhg;
              break;

            case  AB_BA :
              canOlap -> spec. orientation = BA_AB;
              canOlap -> ahg = - tempOlap1 -> endpos;
              canOlap -> bhg = - tempOlap1 -> begpos;
              break;

            case  BA_AB :
              canOlap -> spec. orientation = AB_BA;
              canOlap -> ahg = - tempOlap1 -> endpos;
              canOlap -> bhg = - tempOlap1 -> begpos;
              break;

            default :
              fprintf (GlobalData->stderrc, "Non_canonical orientation = %c\n",
                       canOlap -> spec . orientation);
              assert (FALSE);
          }

        fprintf(GlobalData->stderrc,">>> Fixing up suspicious overlap (" F_CID "," F_CID ",%c) (ahg:" F_COORD " bhg:" F_COORD ") to (" F_CID "," F_CID ",%c) (ahg:" F_COORD " bhg:" F_COORD ") len: " F_COORD "\n",
                inSpec.cidA, inSpec.cidB,
                inSpec.orientation,
                tempOlap1->begpos, tempOlap1->endpos,
                canOlap->spec.cidA, canOlap->spec.cidB,
                canOlap->spec.orientation,
                canOlap->ahg, canOlap->bhg,
                canOlap->overlap);

        // Add it to the symbol table
        InsertChunkOverlap(chunkOverlapper, canOlap);
      }
    }
  }
}



static
int checkChunkOverlapCheckT(ChunkOverlapCheckT *co1,
                            CDS_COORD_t minOverlap,
                            CDS_COORD_t maxOverlap,
                            float errorRate)
{
  if( co1->errorRate != errorRate )
    return FALSE;
  if( minOverlap >= co1->minOverlap && maxOverlap <= co1->maxOverlap &&
      ( minOverlap <= co1->overlap && maxOverlap >= co1->overlap ))
    return TRUE;

  return FALSE;
}



//external
ChunkOverlapCheckT OverlapChunks( GraphCGW_T *graph,
                                  CDS_CID_t cidA, CDS_CID_t cidB,
                                  ChunkOrientationType orientation,
                                  CDS_COORD_t minOverlap,
                                  CDS_COORD_t maxOverlap,
                                  float errorRate,
                                  int insertGraphEdges)
{
  /* this function takes two chunks cidA and cidB, their orientation
     and an assumed minimum and maximum overlap for which it checks.
     It then tries to lookup whether the overlap was computed earlier
     or, if not, it computes the overlap and inserts the result in the
     hash table only if there was not such symbol there before.
  */

  int recompute = FALSE;
  int insert    = FALSE;
  int isCanonical;
  // If we init olap the return value is stored
  // here indicating whether the orientation of the two chunks
  // was canonical in Saul's definition or not.

  ChunkOverlapCheckT *lookup;
  // This pointer holds the return value of LookupCanonicalOverlap

  ChunkOverlapCheckT olap = {0};
  // This is a temporary variable. The return value will be in lookup
  // or the pointer returned by the lookup following the insert


  isCanonical = InitCanonicalOverlapSpec(cidA,cidB,orientation,&olap.spec);
  // initalize olap with the chunk IDs and their orientation and record
  // whether the orientation was already canonical or not

  lookup = LookupCanonicalOverlap(graph->overlapper,&olap.spec);
  // lookup whether the overlap had already been computed


  if( lookup != NULL ){
    olap = *lookup;
    if( checkChunkOverlapCheckT(lookup,minOverlap,maxOverlap,errorRate) == FALSE )
      recompute = TRUE;
    insert = FALSE;
  }
  else
    {
      recompute = TRUE;
      insert = insertGraphEdges;
    }

  if( recompute == TRUE ){
    // compute new overlap and store it into the table
    // If it has an overlap add an edge mate to the CG
    // and return TRUE
    olap.computed      = FALSE;
    // olap.overlap       = 0;
    olap.overlap       = (minOverlap + maxOverlap) / 2;
    olap.minOverlap    = minOverlap;
    olap.maxOverlap    = maxOverlap;
    olap.fromCGB       = FALSE;
    olap.cgbMinOverlap = minOverlap;
    olap.cgbMaxOverlap = maxOverlap;
    olap.errorRate     = errorRate;
    olap.suspicious = FALSE;

    ComputeCanonicalOverlap_new(graph, &olap);

    if(insert)
      { // Insert new entries in hashtable, if requested
        //
        int suspicious = olap.suspicious;
        if(olap.suspicious){
          olap.suspicious = FALSE;
          lookup = LookupCanonicalOverlap(graph->overlapper,&olap.spec); // see if it is already in the table
          if(!lookup){
            if(InsertChunkOverlap(graph->overlapper, &olap) != HASH_SUCCESS) {
              fprintf(stderr, "Failed to insert overlap into hash table.\n");
              assert(0);
            }
          }
        }else{
          if(InsertChunkOverlap(graph->overlapper, &olap) != HASH_SUCCESS) {
            fprintf(stderr, "Failed to insert overlap into hash table.\n");
            assert(0);
          }
        }
        lookup = LookupCanonicalOverlap(graph->overlapper,&olap.spec);
        assert(lookup != NULL);
        // ComputeCanonicalOverlap does not return the olap, so we look it up again.

        olap = *lookup;
        olap.suspicious = suspicious;
      }


    if(olap.overlap && insert){ // Insert all non-zero overlaps, if requested
      InsertComputedOverlapEdge(graph, &olap);
    }
    /* Make sure the orientation of the edge we return is IDENTICAL with the spec returned */
    // if the input was not canonical we set the cid's and orientation
    // back to the input value (see als LookupOverlap)
    if( !olap.suspicious  && ! isCanonical ){
      int swap;

      olap.spec.orientation = orientation;

      // If this is non-canonical, swap things back
      olap.spec.cidA = cidA;
      olap.spec.cidB = cidB;
      swap = olap.BContainsA;
      olap.BContainsA = olap.AContainsB;
      olap.AContainsB = swap;
      swap = olap.ahg;
      olap.ahg = olap.bhg;
      olap.bhg = swap;
    }
  }

  if(olap.overlap==0){olap.quality=0;}
  return olap;
}




//external
Overlap* OverlapContigs(NodeCGW_T *contig1, NodeCGW_T *contig2,
                        ChunkOrientationType *overlapOrientation,
                        CDS_COORD_t minAhang, CDS_COORD_t maxAhang,
                        int computeAhang)
{
  Overlap *tempOlap1;
  char *seq1, *seq2;
  double erate, thresh;
  CDS_COORD_t minlen;

static VA_TYPE(char) *consensus1 = NULL;
static VA_TYPE(char) *consensus2 = NULL;
static VA_TYPE(char) *quality1 = NULL;
static VA_TYPE(char) *quality2 = NULL;

  assert((0.0 <= AS_CGW_ERROR_RATE) && (AS_CGW_ERROR_RATE <= AS_MAX_ERROR_RATE));
  erate = AS_CGW_ERROR_RATE;
  thresh = CGW_DP_THRESH;
  minlen = CGW_DP_MINLEN;

  // if computeAhang is TRUE, allow a lot of slide in ahang
  if (computeAhang == TRUE)
    {
      minAhang = - (CDS_COORD_t) contig2->bpLength.mean;
      maxAhang = (CDS_COORD_t) contig1->bpLength.mean;
      minlen -= 3;  // we subtract 3 because of an asymmetry in DP_COMPARE re AB_BA and BA_AB
    }

  if(consensus1 == NULL)
    {
      consensus1 = CreateVA_char(2048);
      consensus2 = CreateVA_char(2048);
      quality1 = CreateVA_char(2048);
      quality2 = CreateVA_char(2048);
    }else{
      ResetVA_char(consensus1);
      ResetVA_char(consensus2);
      ResetVA_char(quality1);
      ResetVA_char(quality2);
    }
  // Get the consensus sequences for both chunks from the Store
  GetConsensus(ScaffoldGraph->RezGraph, contig1->id, consensus1, quality1);
  GetConsensus(ScaffoldGraph->RezGraph, contig2->id, consensus2, quality2);

  seq1 = Getchar(consensus1, 0);
  seq2 = Getchar(consensus2, 0);

  // tempOlap1 is a static down inside of DP_Compare
  tempOlap1 =
    OverlapSequences( seq1, seq2,
                      *overlapOrientation, minAhang, maxAhang,
                      erate, thresh, minlen);

  return(tempOlap1);
}












// This is the top level routine that computes all new potential overlaps.
//
#define NUM_SECTIONS (5)

//external
void ComputeOverlaps(GraphCGW_T *graph, int addEdgeMates,
                     int recomputeCGBOverlaps)
{
  int i = 0;
  HashTable_Iterator_AS iterator;
  uint64 key, value;
  uint32 valuetype;
  int sectionOuter, sectionOuterMin, sectionOuterMax;
  int sectionInner, sectionInnerMin, sectionInnerMax;
  int numOverlaps = 0;

  // Iterate over all hashtable elements, computing overlaps
  for ( sectionOuter = 0; sectionOuter < NUM_SECTIONS; sectionOuter++)
    {
      sectionOuterMin = sectionOuter * (GetNumGraphNodes(graph)) / NUM_SECTIONS;
      sectionOuterMax = (sectionOuter + 1) * (GetNumGraphNodes(graph)) / NUM_SECTIONS;

      for ( sectionInner = 0; sectionInner < NUM_SECTIONS; sectionInner++)
	{
	  sectionInnerMin = sectionInner * (GetNumGraphNodes(graph)) / NUM_SECTIONS;
	  sectionInnerMax = (sectionInner + 1) * (GetNumGraphNodes(graph)) / NUM_SECTIONS;

	  fprintf(GlobalData->stderrc,"ComputeOverlaps section (o %d,i %d) outer:[%d,%d) inner:[%d,%d)\n",
                  sectionOuter,  sectionInner,
                  sectionOuterMin, sectionOuterMax,
                  sectionInnerMin, sectionInnerMax);

	  InitializeHashTable_Iterator_AS(graph->overlapper->hashTable, &iterator);

	  while(NextHashTable_Iterator_AS(&iterator, &key, &value, &valuetype))
            {
              ChunkOverlapCheckT *olap = (ChunkOverlapCheckT*)(INTPTR)value;

              assert(key == value);

              {
                int inSection = FALSE;
                CDS_CID_t smaller, bigger;

                smaller = MIN( olap->spec.cidA, olap->spec.cidB);
                bigger = MAX( olap->spec.cidA, olap->spec.cidB);

                inSection = (smaller < sectionOuterMax && smaller >= sectionOuterMin) &&
                  (bigger < sectionInnerMax && bigger >= sectionInnerMin);

                // Only do the overlaps where the larger of the ids is within range.
                // The overlaps order of appearance is sorted by the smaller of the indices.
                if(olap->computed || !inSection)
                  continue;
              }


              if(!olap->computed &&  /* We haven't computed this one, and it isn't from cgb, or recomputeCGBOverlaps is set */
                 (!olap->fromCGB || recomputeCGBOverlaps)){
                ChunkOrientationType orientation = olap->spec.orientation;

                // set errRate to old value
                assert((0.0 <= AS_CGW_ERROR_RATE) && (AS_CGW_ERROR_RATE <= AS_MAX_ERROR_RATE));
                olap->errorRate = AS_CGW_ERROR_RATE;

                // first we trust that overlap
                olap->suspicious = FALSE;

                if(olap->maxOverlap < 0){ // Dummy!  Who put this overlap in the table?  No overlap is possible.....SAK
                  olap->overlap = 0;
                  continue;
                }
                if((++i % 100000) == 0){
                  fprintf(GlobalData->stderrc,
                          "* ComputeOverlaps %d  (" F_CID "," F_CID ",%c)\n",
                          i, olap->spec.cidA, olap->spec.cidB,
                          olap->spec.orientation);
                }

                numOverlaps++;

                {
		  ChunkOverlapSpecT inSpec;
		  inSpec = olap->spec;
		  ComputeCanonicalOverlap_new(graph, olap);
		  if(olap->suspicious)
                    {
                      int lengthA, lengthB;
                      lengthA = GetConsensus(graph, olap->spec.cidA, consensusA, qualityA);
                      lengthB = GetConsensus(graph, olap->spec.cidB, consensusB, qualityB);

                      fprintf(GlobalData->stderrc,"* CO: SUSPICIOUS Overlap found! Looked for (" F_CID "," F_CID ",%c)[" F_COORD "," F_COORD "]"
                              "found (" F_CID "," F_CID ",%c) " F_COORD "; contig lengths as found (%d,%d)\n",
                              inSpec.cidA, inSpec.cidB, orientation, olap->minOverlap, olap->maxOverlap,
                              olap->spec.cidA, olap->spec.cidB, olap->spec.orientation, olap->overlap,
                              lengthA,lengthB);
                    }
                }

                if(addEdgeMates && !olap->fromCGB && olap->overlap)
                  InsertComputedOverlapEdge(graph, olap);
              }
            }
	}
    }
}
