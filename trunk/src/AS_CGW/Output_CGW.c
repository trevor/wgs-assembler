
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
static char *rcsid = "$Id: Output_CGW.c,v 1.38 2009-02-02 13:51:14 brianwalenz Exp $";

#include <assert.h>
#include <math.h>
#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "UtilsREZ.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "Globals_CGW.h"
#include "ChiSquareTest_CGW.h"

VA_DEF(IntMate_Pairs);

static VA_TYPE(IntMate_Pairs) *JumpList = NULL;

/*********/

/* This Routine outputs the MateDist messages and also set the
   outMateStat field for the IntAugFrag messages */

void OutputMateDists(ScaffoldGraphT *graph){
  int                   i;
  GenericMesg		pmesg;
  IntMateDistMesg	imd;
  DistT			*dptr;

  pmesg.m = &imd;
  pmesg.t = MESG_IMD;

  assert(graph->doRezOnContigs);

  for(i = 1; i < GetNumDistTs(graph->Dists); i++){
    dptr = GetDistT(graph->Dists, i);

    //  Believe whatever estimate is here.  We used to reset to zero
    //  and the input (except we had already munged the input stddev)
    //  if there were 30 or fewer samples.

    imd.refines     = i;
    imd.mean        = dptr->mu;
    imd.stddev      = dptr->sigma;
    imd.min         = CDS_COORD_MIN;
    imd.max         = CDS_COORD_MAX;
    imd.num_buckets = 0;
    imd.histogram   = NULL;

    // the histogram does not get stored in a checkpoint
    // if the current run of CGW did not have enough samples to recompute the histogram, we have to live without it 
    if (dptr->histogram != NULL) {
      imd.min         = dptr->min;
      imd.max         = dptr->max;
      imd.num_buckets = dptr->bnum;
      imd.histogram   = dptr->histogram;
    }

    if (GlobalData->cgwfp)
      WriteProtoMesg_AS(GlobalData->cgwfp,&pmesg);

    safe_free(dptr->histogram);
    dptr->histogram  = NULL;
    dptr->numSamples = 0;
    dptr->bnum       = 0;
  }
}


/* Must be called after OutputMateDists */
void OutputFrags(ScaffoldGraphT *graph){
  CDS_CID_t		i;

  int goodMates = 0;
  int ctenUntrust=0;
  int ctenTrust  =0;
  int ctrusted   =0;
  int cuntrust   =0;
  int clongvar   =0;
  int cunknown   =0;
  int cinvalid   =0;
  int cwrongScf  =0;

  // Output fragments in iid order
  //
  for(i = 0; i < GetNumInfoByIIDs(graph->iidToFragIndex); i++) {
    CIFragT          *cifrag = NULL;
    InfoByIID        *info   = GetInfoByIID(graph->iidToFragIndex, i);
    GenericMesg       pmesg;
    IntAugFragMesg    iaf;

    if(!info->set)
      continue;

    cifrag = GetCIFragT(graph->CIFrags, info->fragIndex);

    assert(cifrag->iid == i);

    switch(cifrag->flags.bits.edgeStatus){
      case INVALID_EDGE_STATUS:             cinvalid++;    break;
      case TRUSTED_EDGE_STATUS:             ctrusted++;    break;
      case TENTATIVE_TRUSTED_EDGE_STATUS:   ctenTrust++;   break;
      case UNTRUSTED_EDGE_STATUS:           cuntrust++;    break;
      case TENTATIVE_UNTRUSTED_EDGE_STATUS: ctenUntrust++; break;
      case LARGE_VARIANCE_EDGE_STATUS:      clongvar++;    break;
      case INTER_SCAFFOLD_EDGE_STATUS:      cwrongScf++;   break;
      case UNKNOWN_EDGE_STATUS:             cunknown++;    break;
    }

    //  Terminator sets the final fragment clear range based on the
    //  fragStore.

    iaf.iaccession     = cifrag->iid;
    iaf.type           = (FragType)cifrag->type;
    iaf.chaff          = cifrag->flags.bits.isChaff;
    iaf.mate_status    = cifrag->flags.bits.mateDetail;
    iaf.clear_rng.bgn  = -1;
    iaf.clear_rng.end  = -1;

    pmesg.m = &iaf;
    pmesg.t = MESG_IAF;

    if (GlobalData->cgwfp)
      WriteProtoMesg_AS(GlobalData->cgwfp,&pmesg);
  }

  // Output mates in iid order
  //
  for(i = 0; i < GetNumInfoByIIDs(graph->iidToFragIndex); i++) {
    CIFragT            *cif1 = NULL, *cif2 = NULL;
    InfoByIID          *inf1 = NULL, *inf2 = NULL;
    GenericMesg         pmesg;
    IntAugMatePairMesg  iam;

    inf1 = GetInfoByIID(graph->iidToFragIndex, i);

    if(!inf1->set)
      continue;

    cif1 = GetCIFragT(graph->CIFrags, inf1->fragIndex);

    if (cif1->mateOf < 1)
      continue;

    cif2 = GetCIFragT(graph->CIFrags, cif1->mateOf);

    if (cif1->iid > cif2->iid)
      continue;

    inf2 = GetInfoByIID(graph->iidToFragIndex, cif2->iid);

    if(!inf2->set)
      continue;

    assert(inf1->fragIndex == cif2->mateOf);
    assert(inf2->fragIndex == cif1->mateOf);

    assert(cif1->flags.bits.edgeStatus == cif2->flags.bits.edgeStatus);
    assert(cif1->flags.bits.mateDetail == cif2->flags.bits.mateDetail);

    iam.fragment1   = cif1->iid;
    iam.fragment2   = cif2->iid;
    iam.mate_status = cif1->flags.bits.mateDetail;

    pmesg.m = &iam;
    pmesg.t = MESG_IAM;

    if (GlobalData->cgwfp)
      WriteProtoMesg_AS(GlobalData->cgwfp,&pmesg);
  }
}


/* This routine not only outputs the PCM messages, but also sets up
   some values used for ICL & IMD messages.  Thus it must be called
   before OutputConigLinks and OutputMateDists */

void MarkContigEdges(void){
  CIScaffoldT *scaffold;
  GraphNodeIterator scaffolds;

  assert(ScaffoldGraph->doRezOnContigs);

  fprintf(GlobalData->stderrc,"* MarkContigEdges\n");

  // Mark the trustedness of the intra-scaffold, inter-contig edges

  InitGraphNodeIterator(&scaffolds, ScaffoldGraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL){
    if(scaffold->type != REAL_SCAFFOLD)
      continue;
    MarkInternalEdgeStatus(ScaffoldGraph, scaffold, PAIRWISECHI2THRESHOLD_CGW,
                           100000000000.0, TRUE, TRUE, 0, FALSE);
  }

  InitGraphNodeIterator(&scaffolds, ScaffoldGraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL){
    ContigT *contig;
    CIScaffoldTIterator Contigs;

    InitCIScaffoldTIterator(ScaffoldGraph, scaffold,TRUE, FALSE, &Contigs);
    while((contig = NextCIScaffoldTIterator(&Contigs)) != NULL){
      GraphEdgeIterator edges;
      EdgeCGW_T *edge;

      InitGraphEdgeIterator(ScaffoldGraph->ContigGraph, contig->id, ALL_END, ALL_EDGES, GRAPH_EDGE_RAW_ONLY , &edges);
      while((edge = NextGraphEdgeIterator(&edges)) != NULL){
        ContigT *mcontig;

        assert(edge->flags.bits.isRaw);
        if((edge->idA != contig->id) || isSingletonOverlapEdge(edge))
          continue;

        mcontig = GetGraphNode(ScaffoldGraph->ContigGraph, edge->idB);

        if(contig->scaffoldID != mcontig->scaffoldID)
          SetEdgeStatus(ScaffoldGraph->ContigGraph, edge, INTER_SCAFFOLD_EDGE_STATUS);

        PropagateEdgeStatusToFrag(ScaffoldGraph->ContigGraph, edge);
      }
    }
  }

  //  This is needed (??) to update the edge status that we screwed up
  //  above.  The huge int disables any update to the distances.
  //
  //ComputeMatePairStatisticsRestricted(CONTIG_OPERATIONS, 2147483647, "MarkContigEdges");
}



/****************************************************************************/
void OutputContigsFromMultiAligns(void){
  GenericMesg		pmesg;
  IntConConMesg		icm_mesg;
  IntUnitigPos		*uptr;
  GraphCGW_T *graph = ScaffoldGraph->ContigGraph;
  GraphNodeIterator     nodes;
  ContigT		*ctg;
  int32 ubufSize = 100;

  pmesg.m = &icm_mesg;
  pmesg.t = MESG_ICM;

  icm_mesg.unitigs = (IntUnitigPos *) safe_malloc(ubufSize*sizeof(IntUnitigPos));

  InitGraphNodeIterator(&nodes, graph, GRAPH_NODE_DEFAULT);
  /* 1st get min and max values */
  while((ctg = NextGraphNodeIterator(&nodes)) != NULL){
    AS_IID    i;

    if(ctg->flags.bits.isChaff){
      //      fprintf(GlobalData->stderrc,"* # Contig " F_CID " is CHAFF\n", ctg->id);
      continue;
    }

    {
      CIScaffoldT *scaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, ctg->scaffoldID);
      AS_IID    numFrag;
      AS_IID    numUnitig;
      IntMultiPos *mp;
      IntUnitigPos *up;

      MultiAlignT *ma = loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, ctg->id, FALSE);

      numFrag = GetNumIntMultiPoss(ma->f_list);
      mp = GetIntMultiPos(ma->f_list,0);
      numUnitig = GetNumIntUnitigPoss(ma->u_list);
      up = GetIntUnitigPos(ma->u_list,0);

      if(numUnitig >= ubufSize){
        ubufSize = numUnitig * 2;
        icm_mesg.unitigs = (IntUnitigPos *) safe_realloc(icm_mesg.unitigs, ubufSize*sizeof(IntUnitigPos));
      }
      uptr = icm_mesg.unitigs;
      for(i = 0; i < numUnitig; i++){
        IntUnitigPos *iup = up + i;
        NodeCGW_T *unitig = GetGraphNode(ScaffoldGraph->CIGraph, iup->ident);
        if(unitig->type == DISCRIMINATORUNIQUECHUNK_CGW){
          uptr[i].type = AS_UNIQUE_UNITIG;
        }else{
          if(unitig->scaffoldID != NULLINDEX){
            if(!unitig->flags.bits.isSurrogate){
              uptr[i].type = AS_ROCK_UNITIG;
            }else  if(unitig->flags.bits.isStoneSurrogate){
              uptr[i].type = AS_STONE_UNITIG;
            }else{
              uptr[i].type = AS_PEBBLE_UNITIG;
            }
          }else{
            uptr[i].type = AS_SINGLE_UNITIG;
          }
        }
        uptr[i].position = iup->position;
        uptr[i].delta_length = iup->delta_length;
        uptr[i].delta = iup->delta;
        if(unitig->type == RESOLVEDREPEATCHUNK_CGW){
          iup->ident = unitig->info.CI.baseID; // map back to the parent of this instance
        }
        uptr[i].ident = iup->ident;
      }
      ctg->outputID = ctg->id ;  // cid++;
      icm_mesg.placed = (scaffold && (scaffold->type == REAL_SCAFFOLD)?AS_PLACED:AS_UNPLACED);
      icm_mesg.iaccession = ctg->outputID;
      icm_mesg.forced = 0;
      icm_mesg.num_pieces = numFrag;
      icm_mesg.num_vars = GetNumIntMultiVars(ma->v_list); // affects .cgw/ICM
      //    icm_mesg.num_vars = 60;   // affects .cgw/ICM
      icm_mesg.pieces = mp;
      icm_mesg.num_unitigs = numUnitig;
      icm_mesg.length = GetMultiAlignLength(ma);
      if(icm_mesg.num_unitigs > 1){
        icm_mesg.consensus = ""; // Getchar(ma->consensus,0);
        icm_mesg.quality = ""; // Getchar(ma->quality,0);
      }else{
        icm_mesg.consensus = Getchar(ma->consensus,0);
        icm_mesg.quality = Getchar(ma->quality,0);
      }

      if(icm_mesg.num_unitigs > 1){
        assert(ctg->scaffoldID != NULLINDEX);
        if (GlobalData->ctgfp)
          WriteProtoMesg_AS(GlobalData->ctgfp,&pmesg);
      }else{
        if(ctg->scaffoldID == NULLINDEX) {// contig is not placed
          NodeCGW_T *unitig = GetGraphNode(ScaffoldGraph->CIGraph, ctg->info.Contig.AEndCI);

          assert(unitig != NULL);
          if(unitig->info.CI.numInstances == 0){ // If this unitig has been placed as a surrogate, don't output contig
            if (GlobalData->ctgfp)
              WriteProtoMesg_AS(GlobalData->ctgfp,&pmesg); // write the contig
          }else{
            // do nothing. The unitig in this contig appears as a surrogate elsewhere in the assembly
          }
        }else{ // Contig is placed
          if (GlobalData->ctgfp)
            WriteProtoMesg_AS(GlobalData->ctgfp,&pmesg); // write the contig
        }
      }
    }

    //clearCacheSequenceDB(sgraph->sequenceDB);
  }
  safe_free(icm_mesg.unitigs);
}

static int SurrogatedSingleUnitigContig( NodeCGW_T* contig)
{
  if(contig->info.Contig.numCI > 1)
    {
      return 0;
    }
  else
    {
      if(contig->scaffoldID == NULLINDEX) // contig is not placed
	{
	  NodeCGW_T *unitig = GetGraphNode(ScaffoldGraph->CIGraph, contig->info.Contig.AEndCI);

	  assert(unitig != NULL);

	  if(unitig->info.CI.numInstances == 0) // this unitig has not been placed as a surrogate
            {
              return 0;
            }
	  else
            {
              return 1;  // The unitig in this contig appears as a surrogate elsewhere in the assembly
            }
	}
      else
	{ // Contig is placed
	  return 0;
	}
    }
}



void OutputContigLinks(ScaffoldGraphT *graph, int outputOverlapOnlyContigEdges)
{
  IntContigLinkMesg		clm;
  GenericMesg			pmesg;
  GraphNodeIterator nodes;
  ContigT *ctg;
  pmesg.m = &clm;
  pmesg.t = MESG_ICL;

  if(JumpList == NULL)
    JumpList = CreateVA_IntMate_Pairs(256);

  InitGraphNodeIterator(&nodes, ScaffoldGraph->ContigGraph, GRAPH_NODE_DEFAULT);
  while((ctg = NextGraphNodeIterator(&nodes)) != NULL){
    ContigT		*mate;
    GraphEdgeIterator	edges;
    CIEdgeT		*edge, *redge;
    CIFragT		*frag;
    int 		edgeTotal;
    int 		edgeCount;	// This var used for sanity checks
    IntMate_Pairs	imp;
    // MateStatType	mstat;

    if(ctg->flags.bits.isChaff)
      continue;

    if (SurrogatedSingleUnitigContig( ctg ))
      continue;

    clm.contig1 = ctg->outputID;
    InitGraphEdgeIterator(ScaffoldGraph->ContigGraph, ctg->id, ALL_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);
    while((edge = NextGraphEdgeIterator(&edges)) != NULL){

      if (edge->idA != ctg->id)
	continue;

      ResetVA_IntMate_Pairs(JumpList);

      mate = GetGraphNode(ScaffoldGraph->ContigGraph, edge->idB);
      if(mate->flags.bits.isChaff)
	continue;

      if (SurrogatedSingleUnitigContig( mate ))
        continue;

      clm.contig2 = mate->outputID;

      /* Don't need to map orientation, always using canonical orientation*/
      clm.orientation = edge->orient;
      if(!isOverlapEdge(edge)){
	clm.overlap_type = AS_NO_OVERLAP;
      }else {
	clm.overlap_type = AS_OVERLAP;
      }


      switch(GetEdgeStatus(edge)){
        case LARGE_VARIANCE_EDGE_STATUS:
        case UNKNOWN_EDGE_STATUS:
        case INTER_SCAFFOLD_EDGE_STATUS:
          clm.status = AS_UNKNOWN_IN_ASSEMBLY;
          break;
        case TENTATIVE_TRUSTED_EDGE_STATUS:
        case TRUSTED_EDGE_STATUS:
          clm.status = AS_IN_ASSEMBLY;
          break;
        case TENTATIVE_UNTRUSTED_EDGE_STATUS:
        case UNTRUSTED_EDGE_STATUS:
          clm.status = AS_BAD;
          break;
        default:
          assert(0 /* Invalid edge status */);
      }

      clm.is_possible_chimera = edge->flags.bits.isPossibleChimera;
      clm.includes_guide = FALSE;
      clm.mean_distance = edge->distance.mean;
      clm.std_deviation = sqrt(edge->distance.variance);
      edgeTotal = clm.num_contributing = edge->edgesContributing;
      if (clm.overlap_type != AS_NO_OVERLAP)
	--edgeTotal;
      if (!edgeTotal && !outputOverlapOnlyContigEdges)
	continue;	// don't output pure overlap edges


      if (edge->flags.bits.isRaw) {
	assert(edgeTotal <= 1);		// sanity check
	if(edgeTotal == 1){
	  frag = GetCIFragT(ScaffoldGraph->CIFrags, edge->fragA);
	  //	  frag->outMateStat = mstat;
	  frag->flags.bits.edgeStatus = GetEdgeStatus(edge);
	  imp.in1 = frag->iid;
	  frag = GetCIFragT(ScaffoldGraph->CIFrags, edge->fragB);
	  frag->flags.bits.edgeStatus = GetEdgeStatus(edge);
	  //	  frag->outMateStat = mstat;
	  imp.in2 = frag->iid;
	}else{
	  imp.in1 = imp.in2 = 0;
	}
	if(isOverlapEdge(edge)){
	  assert(outputOverlapOnlyContigEdges);
	  imp.type = 'X';
	}else{
          imp.type = AS_MATE;
	  AppendIntMate_Pairs(JumpList, &imp);
	}
      }
      else {
	redge = edge;

	assert(redge->nextRawEdge != NULLINDEX); // must have >= 1 raw edge
	edgeCount = 0;

	while (redge->nextRawEdge != NULLINDEX) {
	  redge = GetGraphEdge(ScaffoldGraph->ContigGraph,redge->nextRawEdge);
	  if (isOverlapEdge(redge))
	    continue;		// overlap edges don't count
	  ++edgeCount;
	  frag = GetCIFragT(ScaffoldGraph->CIFrags, redge->fragA);
	  imp.in1 = frag->iid;
	  frag->flags.bits.edgeStatus = GetEdgeStatus(edge);
	  frag = GetCIFragT(ScaffoldGraph->CIFrags, redge->fragB);
	  frag->flags.bits.edgeStatus = GetEdgeStatus(edge);
	  imp.in2 = frag->iid;
          assert(!isOverlapEdge(redge));
          imp.type = AS_MATE;
          AppendIntMate_Pairs(JumpList, &imp);
	}
	assert(GetNumIntMate_Pairss(JumpList) == edgeTotal);
	assert(edgeCount == edgeTotal);
      }		// if (edge . . .
      clm.jump_list = GetIntMate_Pairs(JumpList,0);

      if (GlobalData->scffp)
        WriteProtoMesg_AS(GlobalData->scffp,&pmesg);
    }	// while (edge . . .
  }	// for (i . . .
}


static void OutputScaffoldLink(ScaffoldGraphT * graph,
                               CIScaffoldT * scaffold,
                               CIEdgeT * edge)
{
  InternalScaffoldLinkMesg slm;
  GenericMesg pmesg;
  CIScaffoldT *mate;
  CIEdgeT *redge;
  CIFragT *frag;
  IntMate_Pairs	imp;
  int edgeTotal = 0;
  int edgeCount = 0; // This var used for sanity checks

  pmesg.m = &slm;
  pmesg.t = MESG_ISL;
  slm.iscaffold1 = scaffold->id;

  if(JumpList == NULL)
    JumpList = CreateVA_IntMate_Pairs(256);
  else
    ResetVA_IntMate_Pairs(JumpList);

  mate = GetGraphNode(ScaffoldGraph->ScaffoldGraph, edge->idB);

  slm.iscaffold2 = mate->id;

  /* Don't need to map orientation, always using canonical orientation*/
  slm.orientation = edge->orient;
  assert(!isOverlapEdge(edge));

  slm.includes_guide = FALSE;
  slm.mean_distance = edge->distance.mean;
  slm.std_deviation = sqrt(edge->distance.variance);
  edgeTotal = slm.num_contributing = edge->edgesContributing;

  redge = edge;
#if 1
  if(edgeTotal < 2)
    return;
#endif
  if (edge->flags.bits.isRaw) {
    assert(edgeTotal <= 1);		// sanity check
    if(edgeTotal == 1){
      frag = GetCIFragT(ScaffoldGraph->CIFrags, edge->fragA);
      frag->flags.bits.edgeStatus = GetEdgeStatus(edge);
      imp.in1 = frag->iid;
      frag = GetCIFragT(ScaffoldGraph->CIFrags, edge->fragB);
      frag->flags.bits.edgeStatus = GetEdgeStatus(edge);
      imp.in2 = frag->iid;
    }else{
      imp.in1 = imp.in2 = 0;
    }
    imp.type = AS_MATE;
    AppendIntMate_Pairs(JumpList, &imp);
    edgeCount = 1;
  }else{

    assert(!edge->flags.bits.isRaw);

    assert(redge->nextRawEdge != NULLINDEX); // must have >= 1 raw edge

    edgeCount = 0;

    while (redge->nextRawEdge != NULLINDEX) {
      redge = GetGraphEdge(ScaffoldGraph->ScaffoldGraph,redge->nextRawEdge);
      assert(!isOverlapEdge(redge));
      ++edgeCount;
      frag = GetCIFragT(ScaffoldGraph->CIFrags, redge->fragA);
      imp.in1 = frag->iid;
      frag->flags.bits.edgeStatus = GetEdgeStatus(edge);
      frag = GetCIFragT(ScaffoldGraph->CIFrags, redge->fragB);
      frag->flags.bits.edgeStatus = GetEdgeStatus(edge);
      imp.in2 = frag->iid;
      imp.type = AS_MATE;
      AppendIntMate_Pairs(JumpList, &imp);
    }
  }
  assert(GetNumIntMate_Pairss(JumpList) == edgeTotal);
  assert(edgeCount == edgeTotal);
  slm.jump_list = GetIntMate_Pairs(JumpList,0);
  if (GlobalData->scffp)
    WriteProtoMesg_AS(GlobalData->scffp,&pmesg);
}


static void OutputScaffoldLinksForScaffold(ScaffoldGraphT * graph,
                                           CIScaffoldT * scaffold)
{
  GraphEdgeIterator	edges;
  CIEdgeT		*edge;

  InitGraphEdgeIterator(graph->ScaffoldGraph, scaffold->id,
                        ALL_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);
  while((edge = NextGraphEdgeIterator(&edges)) != NULL)
    {
      if (edge->idA != scaffold->id)
        continue;
      OutputScaffoldLink(graph, scaffold, edge);
    }
}


void OutputScaffoldLinks(ScaffoldGraphT *graph)
{
  GraphNodeIterator nodes;
  CIScaffoldT *scaffold;

  fprintf(GlobalData->stderrc,"* OutputScaffoldLinks *\n");
  if(JumpList == NULL){
    fprintf(GlobalData->stderrc,"* Creating JumpList *\n");
    JumpList = CreateVA_IntMate_Pairs(256);
  }

  InitGraphNodeIterator(&nodes, graph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while((scaffold = NextGraphNodeIterator(&nodes)) != NULL)
    {
      OutputScaffoldLinksForScaffold(graph, scaffold);
    }
}

/********************************************************************************/
void OutputUnitigsFromMultiAligns(void){
  GenericMesg			pmesg;
  ContigT			*ci;
  IntUnitigMesg			ium_mesg;
  GraphNodeIterator nodes;
  int numCIs = (int) GetNumGraphNodes(ScaffoldGraph->CIGraph);
  CDS_CID_t cid = 0;

  pmesg.m = &ium_mesg;
  pmesg.t = MESG_IUM;

  InitGraphNodeIterator(&nodes, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while((ci = NextGraphNodeIterator(&nodes)) != NULL){
    UnitigStatus   status;

    assert(ci->id>=0 && ci->id< numCIs);

    if(ci->flags.bits.isChaff){
      //      fprintf(GlobalData->stderrc,"* # Unitig " F_CID " is CHAFF\n", ci->id);
      continue;
    }

    switch(ci->type){
      case DISCRIMINATORUNIQUECHUNK_CGW:
        status = AS_UNIQUE;
        //	fprintf(GlobalData->stderrc,"* Unitig " F_CID " is UNIQUE: DISCRIMINATOR output " F_CID " \n",ci->id, cid);
        break;
      case UNIQUECHUNK_CGW:
        status = AS_UNIQUE;
        //	fprintf(GlobalData->stderrc,"* Unitig " F_CID " is UNIQUE: output " F_CID " \n",ci->id, cid);
        break;
      case UNRESOLVEDCHUNK_CGW:
        if(ci->info.CI.numInstances > 0){
          assert(!ci->flags.bits.isUnique);
          status = AS_SEP;
          //	fprintf(GlobalData->stderrc,"* Unitig " F_CID " has %d instances--- output " F_CID " SEP\n",ci->id, ci->info.CI.numInstances,cid);
        }else{
          if(ci->scaffoldID != NULLINDEX){
            //	  fprintf(GlobalData->stderrc,"* Unitig " F_CID " has %d instances--- output " F_CID " UNIQUE\n",ci->id, ci->info.CI.numInstances,cid);
            status = AS_UNIQUE;
          }else{
            //	  fprintf(GlobalData->stderrc,"* Unitig " F_CID " has %d instances--- output " F_CID " NOTREZ\n",ci->id, ci->info.CI.numInstances,cid);
            status = AS_NOTREZ;
          }
        }
        break;
      case RESOLVEDREPEATCHUNK_CGW:
        /* SKIP THESE */
        //      fprintf(GlobalData->stderrc,"* Skipping unitig " F_CID " --- RESOLVEDREPEAT\n",ci->id);
        continue;
      default:
        assert(0);
    }

    {
      MultiAlignT *ma = loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, ci->id, TRUE);
      AS_IID    numFrag = GetNumIntMultiPoss(ma->f_list);
      assert (ci->type != CONTIG_CGW);

      ci->outputID = cid++;
      //assert(ci->outputID == ci->id); // TRUE FOR UNITIGS UNTIL WE SPLIT

      ium_mesg.iaccession = ci->id;
      ium_mesg.coverage_stat = ci->info.CI.coverageStat;
      ium_mesg.microhet_prob = ci->info.CI.microhetProb;
      ium_mesg.status = status;
      ium_mesg.unique_rept = ci->info.CI.forceUniqueRepeat;
      ium_mesg.length = GetMultiAlignLength(ma);
      ium_mesg.consensus = Getchar(ma->consensus,0);
      ium_mesg.quality = Getchar(ma->quality,0);
      ium_mesg.forced = 0;
      ium_mesg.num_frags = GetNumIntMultiPoss(ma->f_list);
      ium_mesg.f_list = GetIntMultiPos(ma->f_list,0);

      if (GlobalData->cgwfp)
        WriteProtoMesg_AS(GlobalData->cgwfp,&pmesg);  //  write the unitig
    }

    //clearCacheSequenceDB(sgraph->sequenceDB);
  }	// while NextGraphNode
}



void OutputUnitigLinksFromMultiAligns(void){
  IntUnitigLinkMesg		ulm;
  GenericMesg			pmesg;
  GraphNodeIterator nodes;
  ChunkInstanceT *ci;
  pmesg.m = &ulm;
  pmesg.t = MESG_IUL;

  fprintf(GlobalData->stderrc,"* OutputUnitigLinksFromMultiAligns *\n");

  if(JumpList == NULL){
    fprintf(GlobalData->stderrc,"* Creating JumpList *\n");
    JumpList = CreateVA_IntMate_Pairs(256);
    AssertPtr(JumpList);
  }

  InitGraphNodeIterator(&nodes, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while((ci = NextGraphNodeIterator(&nodes)) != NULL){
    ContigT		*mate;
    GraphEdgeIterator	edges;
    CIEdgeT		*edge, *redge;
    CIFragT		*frag;
    int 		edgeTotal;
    int 		edgeCount;	// This var used for sanity checks
    IntMate_Pairs     imp;

    AssertPtr(ci);
    assert (ci->type != CONTIG_CGW);

    // We skip these...
    if(ci->type == RESOLVEDREPEATCHUNK_CGW)
      continue;
    if(ci->flags.bits.isChaff)
      continue;

    if(ci->id % 50000 == 0)
      fprintf(GlobalData->stderrc,"* Outputing links incident on unitig " F_CID "\n", ci->id);

    ulm.unitig1 = ci->id;
    InitGraphEdgeIterator(ScaffoldGraph->CIGraph, ci->id, ALL_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);
    while((edge = NextGraphEdgeIterator(&edges)) != NULL){
      AssertPtr(edge);
      ResetVA_IntMate_Pairs(JumpList);
      if (edge->idA != ci->id ||
          edge->flags.bits.isInferred ||
          edge->flags.bits.isInferredRemoved ||
          edge->flags.bits.isMarkedForDeletion)
        continue;
      ulm.unitig2 = edge->idB;
      /* Don't need to map orientation, always using canonical orientation*/
      ulm.orientation = edge->orient;
      if(!isOverlapEdge(edge)){
        ulm.overlap_type = AS_NO_OVERLAP;
      }else {
        ulm.overlap_type = AS_OVERLAP;
      }

      ulm.is_possible_chimera = edge->flags.bits.isPossibleChimera;
      ulm.includes_guide = FALSE;
      ulm.mean_distance = edge->distance.mean;
      ulm.std_deviation = sqrt(edge->distance.variance);
      edgeTotal = ulm.num_contributing = edge->edgesContributing;

      if (ulm.overlap_type != AS_NO_OVERLAP)
        --edgeTotal;
      if (!edgeTotal)
        continue;	// don't output pure overlap edges

      mate = GetGraphNode(ScaffoldGraph->CIGraph, edge->idB);
      if(mate->flags.bits.isChaff)
        continue;

      {
        int numBad = 0;
        int numGood = 0;
        int numUnknown = 0;
        CIFragT *fragA, *fragB;
        // Look through the fragment pairs in this edge.  If any of the fragments are
        // marked BAD ==> bad
        // Else, if any are marked good ==> good
        // Else, mark it unknown
        if(edge->flags.bits.isRaw){
          redge = edge;
        }else{
          redge = GetGraphEdge(ScaffoldGraph->CIGraph, edge->nextRawEdge);
        }
        assert(redge && edge);

        for(; redge != NULL; redge = GetGraphEdge(ScaffoldGraph->CIGraph, redge->nextRawEdge)){
          AssertPtr(redge);
          if(isOverlapEdge(redge))
            continue;
          fragA = GetCIFragT(ScaffoldGraph->CIFrags, redge->fragA);
          fragB = GetCIFragT(ScaffoldGraph->CIFrags, redge->fragB);

          assert(fragA && fragB);
          assert(fragA->flags.bits.edgeStatus== fragB->flags.bits.edgeStatus);

          if(fragA->flags.bits.edgeStatus == UNTRUSTED_EDGE_STATUS ||
             fragA->flags.bits.edgeStatus == TENTATIVE_UNTRUSTED_EDGE_STATUS){
            numBad++;
          }else if(fragA->flags.bits.edgeStatus == TRUSTED_EDGE_STATUS ||
                   fragA->flags.bits.edgeStatus == TENTATIVE_TRUSTED_EDGE_STATUS){
            numGood++;
          } else{
            numUnknown++;
          }
        }

        if(numBad > 0){
          ulm.status = AS_BAD;
        }else if(numGood > 0){
          ulm.status = AS_IN_ASSEMBLY;
        }else ulm.status = AS_UNKNOWN_IN_ASSEMBLY;

      }

      ResetVA_IntMate_Pairs(JumpList);

      if (edge->flags.bits.isRaw) {
        assert(edgeTotal == 1);		// sanity check
        edgeCount = 1;
        frag = GetCIFragT(ScaffoldGraph->CIFrags, edge->fragA);
        AssertPtr(frag);
        imp.in1 = frag->iid;
        frag = GetCIFragT(ScaffoldGraph->CIFrags, edge->fragB);
        imp.in2 = frag->iid;
        AssertPtr(frag);
        assert(!isOverlapEdge(edge));
        imp.type = AS_MATE;
        AppendIntMate_Pairs(JumpList,&imp);
      } else { // not raw
        redge = edge;

        assert(redge->nextRawEdge != NULLINDEX); // must have >= 1 raw edge
        edgeCount = 0;
        assert(edgeTotal > 0);

        while(redge->nextRawEdge != NULLINDEX) {
          redge = GetGraphEdge(ScaffoldGraph->CIGraph,redge->nextRawEdge);
          if (isOverlapEdge(redge))
            continue;		// overlap edges don't count
          ++edgeCount;
          frag = GetCIFragT(ScaffoldGraph->CIFrags, redge->fragA);
          AssertPtr(frag);
          imp.in1 = frag->iid;
          frag = GetCIFragT(ScaffoldGraph->CIFrags, redge->fragB);
          AssertPtr(frag);
          imp.in2 = frag->iid;
          assert(!isOverlapEdge(redge));
          imp.type = AS_MATE;

          AppendIntMate_Pairs(JumpList,&imp);
        }
      }
      assert(GetNumIntMate_Pairss(JumpList) == edgeTotal);

      if(edgeCount != edgeTotal){
        fprintf(GlobalData->stderrc,"* edgeCount = %d edgeTotal = %d\n",
                edgeCount, edgeTotal);
        PrintGraphEdge(GlobalData->stderrc,ScaffoldGraph->CIGraph," ", edge, edge->idA);
        fflush(GlobalData->stderrc);

        redge = edge;
        assert(redge->nextRawEdge != NULLINDEX); // must have >= 1 raw edge

        while(redge->nextRawEdge != NULLINDEX) {
          redge = GetGraphEdge(ScaffoldGraph->CIGraph,redge->nextRawEdge);
          PrintGraphEdge(GlobalData->stderrc,ScaffoldGraph->CIGraph," ", redge, redge->idA);
        }
        assert(edgeCount == edgeTotal);
      }
      ulm.jump_list = GetIntMate_Pairs(JumpList,0);

      if (GlobalData->cgwfp)
        WriteProtoMesg_AS(GlobalData->cgwfp,&pmesg);  //  write the unitig link
    }
  }
}



//#include "obsolete/output_unitig_links"


void OutputScaffolds(ScaffoldGraphT *graph)
{
  AS_IID			sid, pairCount;
  IntScaffoldMesg		ism;
  int				buffSize=2048;
  GenericMesg			pmesg;
  IntContigPairs		*cptr;
  CIScaffoldT			*scaf;
  CIScaffoldTIterator		Contigs;
  ChunkInstanceT		*curr, *last;
  GraphNodeIterator             scaffolds;
  AS_IID    cnt = 0;

  pmesg.m = &ism;
  pmesg.t = MESG_ISF;

  ism.contig_pairs = (IntContigPairs *) safe_malloc(sizeof(IntContigPairs)*buffSize);
  assert(ism.contig_pairs != NULL);

  InitGraphNodeIterator(&scaffolds, graph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while((scaf = NextGraphNodeIterator(&scaffolds)) != NULL){
    ChunkOrient orientLast, orientCurr;
    sid = scaf->outputID = cnt++;

    ism.iaccession = scaf->id;
    ism.num_contig_pairs = scaf->info.Scaffold.numElements-1;
    if(scaf->type != REAL_SCAFFOLD)
      continue;
    assert(ism.num_contig_pairs >= 0);

    if (ism.num_contig_pairs > buffSize) {
      buffSize = ism.num_contig_pairs * 2;
      ism.contig_pairs = (IntContigPairs *) safe_realloc(ism.contig_pairs, sizeof(IntContigPairs)*buffSize);
      assert(ism.contig_pairs != NULL);
    }
    cptr = ism.contig_pairs;

    InitCIScaffoldTIterator(graph, scaf, TRUE, FALSE, &Contigs);
    last = NextCIScaffoldTIterator(&Contigs);
    orientLast = (last->offsetAEnd.mean < last->offsetBEnd.mean) ? A_B : B_A;

    assert(last->scaffoldID == scaf->id);

    if(ism.num_contig_pairs == 0){
      cptr->contig1 = last->outputID;
      cptr->contig2 = last->outputID;
      cptr->mean = 0.0;
      cptr->stddev = 0.0;
      cptr->orient = AB_AB; // got to put something
    }else{
      pairCount = 0;	// only used for sanity check
      while((curr = NextCIScaffoldTIterator(&Contigs)) != NULL){
        assert(pairCount < ism.num_contig_pairs);

        assert(curr->scaffoldID == scaf->id);

        cptr->contig1 = last->outputID;
        cptr->contig2 = curr->outputID;
        orientCurr = (curr->offsetAEnd.mean < curr->offsetBEnd.mean) ? A_B : B_A;
        if(orientLast == A_B){
          if(orientCurr == A_B){
            cptr->mean = curr->offsetAEnd.mean - last->offsetBEnd.mean;
            cptr->stddev = sqrt(curr->offsetAEnd.variance -
                                last->offsetBEnd.variance);
            cptr->orient = AB_AB;
          }else{//orientCurr == B_A
            cptr->mean = curr->offsetBEnd.mean - last->offsetBEnd.mean;
            cptr->stddev = sqrt(curr->offsetBEnd.variance -
                                last->offsetBEnd.variance);
            cptr->orient = AB_BA;
          }
        }else{//orientLast == B_A
          if(orientCurr == A_B){
            cptr->mean = curr->offsetAEnd.mean - last->offsetAEnd.mean;
            cptr->stddev = sqrt(curr->offsetAEnd.variance -
                                last->offsetAEnd.variance);
            cptr->orient = BA_AB;
          }else{//orientCurr == B_A
            cptr->mean = curr->offsetBEnd.mean - last->offsetAEnd.mean;
            cptr->stddev = sqrt(curr->offsetBEnd.variance -
                                last->offsetAEnd.variance);
            cptr->orient = BA_BA;
          }
        }
        last = curr;
        orientLast = orientCurr;
        ++cptr;
        ++pairCount;
      }		// while (curr . . .
    }
    if (GlobalData->scffp)
      WriteProtoMesg_AS(GlobalData->scffp,&pmesg);
  }		// for (sid=0; . . .
  safe_free(ism.contig_pairs);
  return;
}

//#include "obsolete/mean-and-variance-according-to-edge"
