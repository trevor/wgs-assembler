
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

#ifndef SCAFFOLD_GRAPH_H
#define SCAFFOLD_GRAPH_H

static const char *rcsid_SCAFFOLD_GRAPH_H = "$Id: ScaffoldGraph_CGW.h,v 1.46 2011-01-25 21:26:56 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_MSG_pmesg.h"
#include "AS_PER_gkpStore.h"
#include "InputDataTypes_CGW.h"
#include "GraphCGW_T.h"
#include "Globals_CGW.h"
#include "MultiAlign.h"
#include "MultiAlignStore.h"
#include "AS_OVS_overlapStore.h"
#include "AS_CGW_dataTypes.h"

#define CGW_MISSED_OVERLAP CGW_DP_MINLEN /* size the overlapper may have missed */
#define MAX_OVERLAP_SLOP_CGW 10

// controls whether or not we assert in CheckCIScaffoldTLength, ScaffoldSanity
//#define STRICT_SCAFFOLD_CHECKING

typedef EdgeCGW_T CIEdgeT;
typedef CIEdgeT SEdgeT;

typedef NodeCGW_T ChunkInstanceT;
typedef ChunkInstanceT ContigT;
typedef ChunkInstanceT CIScaffoldT;



void ClearChunkInstance(ChunkInstanceT *ci);

void InitializeChunkInstance(ChunkInstanceT *ci, ChunkInstanceType type);

void InitializeContig(ContigT *contig, ChunkInstanceType type);

void InitializeScaffold(CIScaffoldT *scaffold, ChunkInstanceType type);

int isDeadCIScaffoldT(CIScaffoldT *scaffold);

int RepeatRez(int repeatRezLevel, char *name);

VA_DEF(ChunkInstanceT)
VA_DEF(CIScaffoldT)
VA_DEF(ContigT)
VA_DEF(CIEdgeT)
VA_DEF(SEdgeT)

typedef struct{
  VA_TYPE(CIFragT)        *CIFrags;
  VA_TYPE(DistT)          *Dists;
  VA_TYPE(ChunkInstanceT) *ChunkInstances;  // CIs and Contigs
  VA_TYPE(ContigT)        *Contigs;
  VA_TYPE(CIScaffoldT)    *CIScaffolds;
  VA_TYPE(CIEdgeT)        *CIEdges;
  VA_TYPE(CIEdgeT)        *ContigEdges;
  VA_TYPE(SEdgeT)         *SEdges;
  char                    name[256];
  int32                   checkPointIteration; // Index of next checkpoint
  int32                   numContigs;  // Number of contigs...they may be interspersed
  int32                   numOriginalCIs;
  int32                   numLiveCIs; // Number of currently instantiated CIs, including splits
  int32                   numDiscriminatorUniqueCIs;
  int32                   numLiveScaffolds;
  GraphCGW_T             *CIGraph;
  GraphCGW_T             *ContigGraph;
  GraphCGW_T             *ScaffoldGraph;
  ChunkOverlapperT       *ChunkOverlaps;
  gkStore                *gkpStore;
  MultiAlignStore        *tigStore;
  OverlapStore           *frgOvlStore;
}ScaffoldGraphT;



/* Constructor */
ScaffoldGraphT *CreateScaffoldGraph(char *name);
void InsertRepeatCIsInScaffolds(ScaffoldGraphT *sgraph);
void BuildCIEdges(ScaffoldGraphT *graph);
void BuildSEdgesForScaffold(ScaffoldGraphT * graph,
                            CIScaffoldT * scaffold,
                            int canonicalOnly,
                            int includeNegativeEdges);
void BuildSEdges(ScaffoldGraphT *graph, int canonicalOnly, int includeNegativeEdges);

/* Destructor */
void DestroyScaffoldGraph(ScaffoldGraphT *sgraph);

/* Maintentance */


/*****************************************************************************
	Operations on Chunk Instances
****************************************************************************/

void DumpChunkInstances(FILE *stream, ScaffoldGraphT *graph, int confirmedOnly,
			int scaffoldedOnly, int uniqueToUniqueOnly,
			int verbose);
void DumpChunkInstance(FILE *stream, ScaffoldGraphT *graph,
                       ChunkInstanceT *chunk,
		       int confirmedOnly, int scaffoldedOnly,
		       int uniqueToUniqueOnly, int verbose);




// Iterators
//
//	Iterate over all Edges incident no a particular CI
//      Iterate over all Edges incident on a pair (i,j)


/****************************************************************************
	Operations on Scaffolds
***************************************************************************/


void InsertCIInScaffold(ScaffoldGraphT *sgraph, CDS_CID_t ci, CDS_CID_t sid,
                        LengthT aEndOffset, LengthT bEndOffset,
                        int AEndToBend, int contigNow);

void  MarkCIElementsForScaffoldMembership(ChunkInstanceT *chunkInstance,
                                          CDS_CID_t scaffoldID);


/*
  RemoveCIFromScaffold
  Remove chunk instance ci from scaffold sid
  Returns 0 if successful.
*/
int RemoveCIFromScaffold(ScaffoldGraphT *sgraph, CIScaffoldT *scaffold,
                         ChunkInstanceT *CI, int adjustPositions);

int ContigContainment(CIScaffoldT *scaffold, NodeCGW_T *prevCI,
                      NodeCGW_T *thisCI, EdgeCGW_T *overlapEdge,
                      int tryHarder);

int CleanupAScaffold(ScaffoldGraphT *graph, CIScaffoldT *scaffold,
                     int lookForSmallOverlaps,
                     int32 maxContigsInMerge,
                     int deleteUnmergedSurrogates);
int CleanupScaffolds(ScaffoldGraphT *graph, int lookForSmallOverlaps,
                     int32 maxContigsInMerge,
                     int deleteUnmergedSurrogates);

/* Try desperately to overcome failues in merge multialigns by
   performing the merges iteratively on subsets of the original failed merge.*/
int CleanupFailedMergesInScaffolds(ScaffoldGraphT *sgraph);

/*   Discard short (2kb or less) contigs that contain ONLY surrogates and
     should have been contigged, except for the fact the Merge Multi aligns
     failed. We delete these from their scaffold, delete the associated
     surrogate contigs and CIs, and adjust surrogate instance counts, as
     required */
int  DeleteAllSurrogateContigsFromFailedMerges(CIScaffoldT *scaffold,
                                               NodeCGW_T *contig);




// Build contigs for a single scaffold
int ContigAScaffold(ScaffoldGraphT *graph, CDS_CID_t sid);

int  CreateAContigInScaffold(CIScaffoldT *scaffold,
                             VA_TYPE(IntElementPos) *ContigPositions,
                             LengthT offsetAEnd, LengthT offsetBEnd);


void CheckContigs(void);

void
dumpContigInfo(ChunkInstanceT *contig);

void
GetContigPositionInScaffold(ChunkInstanceT *contig, int *left_end, int *right_end,
                            int *contigScaffoldOrientation);

void
GetFragmentPositionInScaffold(CIFragT *frag,
                              int *left_end, int *right_end,
                              int *fragmentScaffoldOrientation);

void DumpContig(FILE *stream, ScaffoldGraphT *graph, ContigT *contig, int raw);

void DumpContigInScfContext(FILE *stream, ScaffoldGraphT *graph,
                            ContigT *contig, int raw);

void DumpContigs(FILE *stream, ScaffoldGraphT *graph, int raw);

#ifndef USE_EARLY_CONTIGGING
/* Convert all of the scaffolds from scaffolds of CIs to scaffolds of Contigs.
   This involves creating the appropriate contigs and inserting them in the
   scaffolds, in lieu of the CIs.  Following the construction of the contigs,
   the contig edges are built and merged. */
int BuildContigs(ScaffoldGraphT *graph);
#endif


// Create a contig for each CI
// Create all of the merged/raw CIEdges for the graph of Contigs
int BuildInitialContigs(ScaffoldGraphT *graph);

#if 1
/* Construct the Contig Edges from the CIEdges */
int BuildContigEdges(ScaffoldGraphT *graph);
#endif


/*
  Iterators

  Iterate over all Scaffold Edges incident no a particular scaffold
  Iterate over all Scaffold Edges incident on a pair (i,j)
*/

/****************************************************************************
	Operations on Scaffold Graph
***************************************************************************/

/*
  Rebuild Scaffold Edges

  Using the CIEdges from all scaffolded chunkInstances, rebuild
  and merge the Scaffold Edges.  We can probably do this incrementally
  at some point.
*/
int RebuildScaffoldEdges(ScaffoldGraphT *sgraph);

/*
  AssignFragmentsToCIs

  Assigns a list of FragInfoT records to each ChunkInstance that is not a
  contig. We expect this to be peformed prior to microhet seperation and/or
  prior to output.
*/
int AssignFragmentsToCIs(ScaffoldGraphT *sgraph);

void PrintCIScaffoldHeader(FILE *stream, ScaffoldGraphT *graph,
                           CIScaffoldT *scaffold);
void DumpCIScaffold(FILE *stream, ScaffoldGraphT *graph,
                    CIScaffoldT *scaffold, int raw);
void DumpCIScaffolds(FILE *stream, ScaffoldGraphT *graph, int raw);
void DumpACIScaffold(FILE *stream, ScaffoldGraphT *graph,
                     CIScaffoldT *scaffold, int raw);
void DumpACIScaffoldNew(FILE *stream, ScaffoldGraphT *graph,
                        CIScaffoldT *scaffold, int raw);
void PrintContigEdgeInScfContext(FILE *fp, GraphCGW_T *graph,
                                 char *label, EdgeCGW_T *edge,
                                 CDS_CID_t cid);

typedef enum {
  RECOMPUTE_OK = 0,
  RECOMPUTE_SINGULAR = 1,
  RECOMPUTE_LAPACK = 2,
  RECOMPUTE_NO_GAPS = 3,
  RECOMPUTE_FAILED_REORDER_NEEDED = 4,
  RECOMPUTE_NOT_ENOUGH_CLONES = 5,
  RECOMPUTE_CONTIGGED_CONTAINMENTS = 6,
  RECOMPUTE_FAILED_CONTIG_DELETED  = 7
}RecomputeOffsetsStatus;

/*
  RecomputeOffsetsInScaffold

  Recomputes the positions of the CIs in a scaffold using a least
  square error approach
  Arguments:
  allowOrderChanges -- if TRUE, reordering the CIs in a scaffold can occur
  if FALSE, reordering will cause a return value
  of RECOMPUTE_FAILED_REORDER_NEEDED
  Preconditions:
  scaffold should be internally connected (IsScaffoldInternallyConnected)
  Return Values:
  RECOMPUTE_OK    CI positions are updated and scaffold's least
  square error measure and number of least square
  clones are set.
*/
RecomputeOffsetsStatus RecomputeOffsetsInScaffold(ScaffoldGraphT *sgraph,
                                                  CDS_CID_t scaffoldID,
                                                  int allowOrderChanges,
                                                  int forceNonOverlaps,
                                                  int verbose);


int IsScaffold2EdgeConnected(ScaffoldGraphT *graph, CIScaffoldT *scaffold);

/*
  IsScaffoldInternallyConnected

  Determines whether the scaffold is connected by edges marked TRUSTED
  and TENTATIVELY_TRUSTED. This is a necessary condition for boths sanity
  and successful recomputation of positions of Scaffold CI positions.  Also
  interesting to evaluate this after MarkInternalCIEdgeStatus. edgeTypes
  defines the set of edges used.  LeastSquares uses ALL_TRUSTED_EDGES,
  other manipulations use ALL_EDGES.
  Returns TRUE if connected, FALSE if not connected.
*/
int IsScaffoldInternallyConnected(ScaffoldGraphT *graph,
                                  CIScaffoldT *scaffold, int32 edgeTypes);


// New test code to partly substitute for the status given by
// MarkInternalEdgeStatus, for to help handle slightly messier cases!
//
int IsInternalEdgeStatusVaguelyOK(EdgeCGW_T *edge,CDS_CID_t thisCIid);



// Uses isScaffoldInternallyConnected to decide if a scaffold is
// connected.  If not, splits the scaffold into its components.
//
int CheckScaffoldConnectivityAndSplit(ScaffoldGraphT *graph,
                                      CDS_CID_t sid,
                                      int32 edgeTypes, int verbose);


void PrintCIEdgeT(FILE *fp, ScaffoldGraphT *graph,
                  char *label, CIEdgeT *edge, CDS_CID_t cid);
void PrintSEdgeT(FILE *fp, ScaffoldGraphT *graph,
                 char *label, SEdgeT *edge, CDS_CID_t sid);



/* Build Scaffolds containing only Unique CIs based on the Unique CI
   subset of the enhanced unitig graph.
*/
void BuildUniqueCIScaffolds(ScaffoldGraphT *graph,
                            int markShakyBifurcations, int verbose);



// Cleans up a scaffold where the first contig has been placed at a negative coordinate.
void  CheckLSScaffoldWierdnesses(char *string, ScaffoldGraphT *graph,
                                 CIScaffoldT *scaffold);


/* Recompute mean and variance gap estimates based on clone mate pairs
   within a scaffold for all scaffolds.
*/

void LeastSquaresGapEstimates(ScaffoldGraphT *graph, int markEdges,
			      int useGuides, int forceNonOverlaps,
			      int checkConnectivity, int verbose);

/***** Celamy *****/
void DumpCelamyColors(FILE *file);
void DumpCelamyMateColors(FILE *file);
void DumpCelamyFragColors(FILE *file);

void CelamyScaffold(FILE *fout, CIScaffoldT *scaffold, int64 scaffoldAEndCoord, int64 scaffoldBEndCoord);
void CelamyAssembly(char *name);


void FindScaffoldComponents(ScaffoldGraphT *graph, int findPaths);
int MergeScaffoldPaths(ScaffoldGraphT *sgraph);


/* Check that all trusted edges are intra-scaffold,
   generates output to log file
*/
void CheckAllTrustedEdges(ScaffoldGraphT * sgraph);

/* Checks that all trusted edges incident on cid are intra-scaffold */
void CheckTrustedEdges(ScaffoldGraphT * sgraph, CDS_CID_t cid);

/* Check means and variance in all real scaffolds */
void CheckCIScaffoldTs(ScaffoldGraphT *sgraph);

/* Check means and variance in a single */
void CheckCIScaffoldT(ScaffoldGraphT *sgraph, CIScaffoldT *scaffold);

/* Scans the scaffold and sets its length */
void SetCIScaffoldTLength(ScaffoldGraphT *sgraph,
                          CIScaffoldT *ciScaffold, int verbose);

/* Scans all scaffolds and sets their lengths */
void SetCIScaffoldTLengths(ScaffoldGraphT *sgraph, int verbose);

void CheckCIScaffoldTLengths(ScaffoldGraphT *sgraph);
void CheckCIScaffoldTLength(ScaffoldGraphT *sgraph, CIScaffoldT *scaffold);

/* Checks that all edges (save those marked untrusted) incident on
   cid are intra-scaffold
*/
int CheckAllEdges(ScaffoldGraphT *sgraph, CDS_CID_t sid, CDS_CID_t cid);

void CheckpointScaffoldGraph(const char *logicalname, const char *location);
void LoadScaffoldGraphFromCheckpoint(char *name, int32 checkPointNum, int writable);

void ReportMemorySize(ScaffoldGraphT *graph, FILE *stream);


void AddInferredScaffoldEdges(ScaffoldGraphT *graph,  int verbose);


/* Insert a computed overlap as a CIEdgeT into the Scaffold Graph. */
int InsertComputedOverlapCIEdgeGraph(ScaffoldGraphT *graph,
                                     ChunkOverlapCheckT *olap);


int GetConsensus(GraphCGW_T *graph, CDS_CID_t CIindex,
                 VA_TYPE(char) *consensus, VA_TYPE(char) *quality);
// Return value is length of unitig or contig  sequence/quality (-1 if failure)
// Consensus and quality are COPIED into the VAs

int GetCoverageStat(ChunkInstanceT *CI);
int GetNumInstances(ChunkInstanceT *CI);


void CheckScaffoldOrder(CIScaffoldT *scaffold, ScaffoldGraphT *graph);

void ScaffoldSanity(CIScaffoldT *scaffold, ScaffoldGraphT *graph);

void DumpScaffoldGraph(ScaffoldGraphT *graph);

void SetCIScaffoldIds(ChunkInstanceT *CI, CDS_CID_t scaffoldID);
void SetContigScaffoldIds(ContigT *contig, CDS_CID_t scaffoldID);


/*
  MarkInternalEdgeStatus

  Using the CI positions, determines the trustworthiness of all CI edges
  that are internal to the scaffold, and marks them with a (new) status.
  This routine only changes the status of edges which are not masked by the
  doNotChange mask which should be the bitwise OR of CIEdgeStatus enums
  the user does not want changed.
  Edges with variance > maxVariance will be marked with
  LARGE_VARIANCE_EDGE_STATUS
  As an additional side effect, it also updates the scaffolds count of
  internal edges and confirmed internal edges.
  If operateOnMerged == TRUE, checks merged edges
  else                        checks raw edges
*/
void MarkInternalEdgeStatus(ScaffoldGraphT *graph,
                            CIScaffoldT *scaffold,
                            float pairwiseChiSquaredThreshhold,
                            float maxVariance,
                            int markTrusted,
                            int markUntrusted,
                            int doNotChange,
                            int operateOnMerged);





/* *********************************************************************** */
/* Add a fixed amount to the offsetAEnd and offsetBEnd starting from a given
   CI to the end of the Scaffold                                */
/* *********************************************************************** */

void AddDeltaToScaffoldOffsets(ScaffoldGraphT *graph,
                               CDS_CID_t scaffoldIndex,
                               CDS_CID_t indexOfCI,
                               int aEndToBEnd,
                               int verbose,
                               LengthT delta);

void AddDeltaToScaffoldOffsets(ScaffoldGraphT *graph,
                               CDS_CID_t scaffoldIndex,
                               CDS_CID_t indexOfCI,
                               int aEndToBEnd,
                               int verbose,
                               LengthT delta,
                               uint32 mark);


void RebuildScaffolds(ScaffoldGraphT *ScaffoldGraph,
                      int markShakyBifurcations);
void TidyUpScaffolds(ScaffoldGraphT *ScaffoldGraph);

void CheckAllContigFragments(void);

/* Globals */
extern ScaffoldGraphT *ScaffoldGraph;

void FixupLengthsScaffoldTs(ScaffoldGraphT *sgraph);
void FixupLengthScaffoldT(ScaffoldGraphT *sgraph, CIScaffoldT *scaffold);

EdgeCGW_T *FindOverlapEdgeChiSquare(ScaffoldGraphT *graph,
                                    NodeCGW_T *sourceCI, CDS_CID_t targetId,
                                    PairOrient edgeOrient,
                                    double inferredMean,
                                    double inferredVariance,
                                    float *chiSquaredValue,
                                    float chiSquareThreshold, int *alternate,
                                    int verbose);

void CheckInternalEdgeStatus(ScaffoldGraphT *graph, CIScaffoldT *scaffold,
                             float pairwiseChiSquaredThreshhold,
                             float maxVariance,
                             int doNotChange, int verbose);

int CheckForContigs(ScaffoldGraphT *sgraph,
                    CDS_CID_t cid, CDS_CID_t sid,
                    LengthT offsetAEnd, LengthT offsetBEnd);

/* DemoteSmallSingletonScaffolds
   We want to demote the contigs/unitigs in small singleton scaffolds so
   that they can be candidates for stone/rock throwing.
*/
void DemoteSmallSingletonScaffolds(void);

void PrintSEdgesForScaffold(ScaffoldGraphT * graph,
                            CIScaffoldT * scaffold);

void ShatterScaffoldsConnectedByLowWeight(FILE *stream, ScaffoldGraphT *graph, uint32 minWeight, int verbose);

#include "CIScaffoldT_Merge_CGW.h"

#endif
