
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
#define FILTER_EDGES
static char CM_ID[] = "$Id: Input_CGW.c,v 1.45 2007-08-28 22:50:10 brianwalenz Exp $";

/*   THIS FILE CONTAINS ALL PROTO/IO INPUT ROUTINES */


//#define DEBUG 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_version.h"
#include "AS_CGW_dataTypes.h"
#include "AS_PER_gkpStore.h"
#include "ScaffoldGraph_CGW.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "Output_CGW.h"
#include "Input_CGW.h"


static int32 DiscriminatorUniques = 0;
static int32 ShortDiscriminatorUniques = 0;

static int32 totalReadFrags = 0;
static int32 inUniqueReadFrags = 0;
static int32 inRepeatReadFrags = 0;
static int32 inTeenyUnitigReadFrags = 0;
static int32 inSingletonUnitigReadFrags = 0;
static int32 onEndReadFrags = 0;

static int32 TouchesContained = 0;
static int32 TransChunk = 0;
static int32 Containment = 0;
static int32 DoveTail = 0;
static int32 Tandem = 0;
static int32 BetweenContained = 0;
static int32 ContainStack = 0;
static int32 BadQuality = 0;


static
void
ProcessFrags(void) {
  CDS_CID_t i;
  int32 unmatedFrags = 0;
  GateKeeperFragmentRecord gkf;
  int                      err = 0;

  //  Do one pass through, reading from the gatekeeper store to fill
  //  out the cifrag info.

  for(i = 0; i < GetNumInfoByIIDs(ScaffoldGraph->iidToFragIndex); i++){
    InfoByIID *ciinfo = GetInfoByIID(ScaffoldGraph->iidToFragIndex, i);
    CIFragT   *cifrag = GetCIFragT(ScaffoldGraph->CIFrags, ciinfo->fragIndex);

    if(!ciinfo->set)
      continue;

    assert(cifrag->iid == i);  //  If !set, this fails.

    getGateKeeperFragment(ScaffoldGraph->gkpStore, i, &gkf);

    if (gkf.mateIID != 0) {
      InfoByIID *miinfo = GetInfoByIID(ScaffoldGraph->iidToFragIndex, gkf.mateIID);

      if (miinfo && miinfo->set) {
        cifrag->mateOf   = miinfo->fragIndex;
        cifrag->dist     = gkf.libraryIID;
        if (gkf.orientation == AS_READ_ORIENT_INNIE)
          cifrag->flags.bits.innieMate = TRUE;
        cifrag->flags.bits.linkType   = AS_MATE;
        cifrag->flags.bits.edgeStatus = UNKNOWN_EDGE_STATUS;
        cifrag->flags.bits.hasMate    = TRUE;
      } else {
        fprintf(stderr, "ProcessFrags()-- WARNING!  fragiid=%d,index=%d mateiid=%d,index=%d -- MATE DOESN'T EXIST!\n",
                i, ciinfo->fragIndex, gkf.mateIID, miinfo->fragIndex);
        //  This is not a critical failure, but does indicate
        //  something amiss with either the store or the unitigs.
        //
        err++;
      }

      //fprintf(stderr, "Frag: iid=%d,index=%d mateiid=%d,index=%d\n", i, ciinfo->fragIndex, gkf.mateIID, miinfo->fragIndex);
    }

    if (cifrag->flags.bits.hasMate == FALSE)
      unmatedFrags++;
  }

  if (err)
    assert(err == 0);

  fprintf(stderr,"* Unmated fragments %d\n", unmatedFrags);
  fprintf(stderr,"* Total IIDs:       %d\n", (int)GetNumInfoByIIDs(ScaffoldGraph->iidToFragIndex));

  //  Do a second pass to clean up mates.  We don't need to access the
  //  gatekeeper here, since everything is already populated.  If we
  //  put this stuff in the first loop, we'll be randomly accessing
  //  the gatekeeper.

  for(i = 0; i < GetNumInfoByIIDs(ScaffoldGraph->iidToFragIndex); i++){
    InfoByIID *ciinfo = GetInfoByIID(ScaffoldGraph->iidToFragIndex, i);
    CIFragT   *cifrag = GetCIFragT(ScaffoldGraph->CIFrags, ciinfo->fragIndex);
    InfoByIID *miinfo = NULL;
    CIFragT   *mifrag = NULL;

    //  this frag not used, or no mate
    if ((!ciinfo->set) || (cifrag->mateOf == NULLINDEX))
      continue;

    mifrag = GetCIFragT(ScaffoldGraph->CIFrags, cifrag->mateOf);
    miinfo = GetInfoByIID(ScaffoldGraph->iidToFragIndex, mifrag->iid);

    if ((mifrag == NULL) || (mifrag->mateOf == NULLINDEX)) {
      //  We set up links to a dead frag, clean up...

      cifrag->mateOf   = NULLINDEX;
      cifrag->dist     = NULLINDEX;
      cifrag->flags.bits.linkType   = AS_UNKNOWN;
      cifrag->flags.bits.edgeStatus = INVALID_EDGE_STATUS;
      cifrag->flags.bits.hasMate    = FALSE;

      if (mifrag) {
        mifrag->mateOf   = NULLINDEX;
        mifrag->dist     = NULLINDEX;
        mifrag->flags.bits.linkType   = AS_UNKNOWN;
        mifrag->flags.bits.edgeStatus = INVALID_EDGE_STATUS;
        cifrag->flags.bits.hasMate    = FALSE;
      }
    } else {
      //  Both guys are alive, and we're mated.  Throw some asserts

      if (ciinfo->set == 0) {
        fprintf(stderr, "ERROR: cifrag iid=%d ciinfo->set == 0; cifrag not in the assembly\n");
      }
      if (miinfo->set == 0) {
        fprintf(stderr, "ERROR: mifrag iid=%d miinfo->set == 0; mifrag not in the assembly\n");
      }
      if (cifrag->dist   != mifrag->dist) {
        fprintf(stderr, "ERROR: cifrag iid=%d mifrag iid=%d -- cifrag->dist=%d != mifrag->dist=%d; libraries not the same\n",
                cifrag->iid, mifrag->iid, cifrag->dist, mifrag->dist);
      }
      if (cifrag->mateOf != miinfo->fragIndex) {
        fprintf(stderr, "ERROR: cifrag iid=%d mifrag iid==%d -- cifrag->mateOf=%d != miinfo->fragIndex=%d; messed up mate index/iid\n",
                cifrag->iid, mifrag->iid, cifrag->mateOf, miinfo->fragIndex);
      }
      if (mifrag->mateOf != ciinfo->fragIndex) {
        fprintf(stderr, "ERROR: mifrag iid=%d cifrag iid==%d -- mifrag->mateOf=%d != ciinfo->fragIndex=%d; messed up mate index/iid\n",
                mifrag->iid, cifrag->iid, mifrag->mateOf, ciinfo->fragIndex);
      }

      assert(ciinfo->set);
      assert(miinfo->set);
      assert(cifrag->dist   == mifrag->dist);
      assert(cifrag->mateOf == miinfo->fragIndex);
      assert(mifrag->mateOf == ciinfo->fragIndex);
    }

  }  //  for each frag

  fprintf(stderr,"* Found %d unmated frags (%g %%)\n",
          unmatedFrags,
          100.0 * unmatedFrags / GetNumCIFragTs(ScaffoldGraph->CIFrags));
}






int ProcessInput(Global_CGW *data, int optind, int argc, char *argv[]){
  GenericMesg   *pmesg;
  FILE *infp;
  int i,j = 0;
  int32 numIUM = 0;

  StartTimerT(&GlobalData->InputTimer);

  for(i = optind; i < argc; i++){
    infp = fopen(argv[i],"r");

    while ((EOF != ReadProtoMesg_AS(infp, &pmesg))) {
      if (pmesg->t == MESG_IUM) {
        IntUnitigMesg *ium_mesg = (IntUnitigMesg *)pmesg->m;
        MultiAlignT   *ma       = CreateMultiAlignTFromIUM(ium_mesg, GetNumCIFragTs(ScaffoldGraph->CIFrags),FALSE);
        CDS_COORD_t   length    = GetMultiAlignUngappedLength(ma);

        InsertMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, ium_mesg->iaccession, TRUE, ma, TRUE);
        DuplicateEntryInSequenceDB(ScaffoldGraph->sequenceDB, ium_mesg->iaccession, TRUE, ium_mesg->iaccession, FALSE, FALSE);
        ProcessIUM_ScaffoldGraph(ium_mesg, length, FALSE);
        UnloadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, ium_mesg->iaccession, TRUE);

        numIUM++;
      }
    }
    fclose(infp);
  }
  
  ScaffoldGraph->numLiveCIs = ScaffoldGraph->numOriginalCIs = GetNumGraphNodes(ScaffoldGraph->CIGraph);

  fprintf(stderr,"* cgw read the following messages:\n");
  fprintf(stderr,"\tIUM:%d  (max IUM acc = %d) with %d fragments\n",
          numIUM,
          (int)GetNumGraphNodes(ScaffoldGraph->CIGraph),
	  (int)GetNumCIFragTs(ScaffoldGraph->CIFrags));
  
  ProcessFrags();

  fprintf(stderr,"* Total Long Discriminator Uniques : %d   Short Uniques: %d\n",
	  DiscriminatorUniques, ShortDiscriminatorUniques);
  fprintf(stderr,"* Total Reads:%d in discriminator unique:%d in other:%d ; in teeny: %d in singles:%d on ends:%d\n",
	  totalReadFrags, inUniqueReadFrags, inRepeatReadFrags, inTeenyUnitigReadFrags, inSingletonUnitigReadFrags,
	  onEndReadFrags);

  StopTimerT(&GlobalData->InputTimer);

  return(0);
}





void ProcessIUM_ScaffoldGraph(IntUnitigMesg *ium_mesg, CDS_COORD_t length, int sequenceOnly){
  CDS_CID_t cfr;
  CDS_COORD_t simLength;
  ChunkInstanceT CI = {0};

  CI.id = ium_mesg->iaccession;
  CI.bpLength.mean = length;
  CI.bpLength.variance = MAX(1.0,ComputeFudgeVariance(CI.bpLength.mean));
  CI.edgeHead = NULLINDEX;
  CI.microhetScore = NULLINDEX;
  CI.setID = NULLINDEX;
  CI.scaffoldID = NULLINDEX;
  CI.indexInScaffold = NULLINDEX;
  CI.prevScaffoldID = NULLINDEX;
  CI.numEssentialA = 0;
  CI.numEssentialB = 0;
  CI.essentialEdgeA = NULLINDEX;
  CI.essentialEdgeB = NULLINDEX;
  CI.smoothExpectedCID = NULLINDEX;
  CI.BEndNext = CI.AEndNext = NULLINDEX;
  CI.info.CI.headOfFragments = GetNumCIFragTs(ScaffoldGraph->CIFrags);
  CI.info.CI.numFragments = ium_mesg->num_frags;
  CI.info.CI.coverageStat = (ium_mesg->coverage_stat < -1000.0? -1000:ium_mesg->coverage_stat);
  CI.info.CI.contigID = NULLINDEX;
  CI.info.CI.numInstances = 0;
  CI.info.CI.instances.in_line.instance1 = 0;
  CI.info.CI.instances.in_line.instance2 = 0;
  CI.info.CI.instances.va = NULL;
  CI.info.CI.source = NULLINDEX;
  CI.flags.all = 0;
  CI.offsetAEnd.mean = 0.0;
  CI.offsetAEnd.variance = 0.0;
  CI.offsetBEnd = CI.bpLength;

#ifdef AS_ENABLE_SOURCE
  if(ium_mesg->source){
    char *c = ium_mesg->source;
    CI.info.CI.source = GetNumchars(ScaffoldGraph->SourceFields);
    while(*c != '\0'){
      Appendchar(ScaffoldGraph->SourceFields,c++);
    }
    Appendchar(ScaffoldGraph->SourceFields,c);
  }else{
    CI.info.CI.source = NULLINDEX;
  }
#endif

  // Collect the microhetScore if available
  CI.microhetScore = 1.01;
#ifdef AS_ENABLE_SOURCE
  {
    char *mhp = strstr(ium_mesg->source,"mhp:");
    if(mhp)
      CI.microhetScore = atof(mhp+4);
      //fprintf(stderr,"* %s\n*  mhp:%g found *\n", ium_mesg->source, CI.microhetScore);
  }
#endif

  // See if this is a repeat, or we can pin it down to an interval
  {
    char *interval;
    char *type;
    int result;
    //	  fprintf(stderr,"* source = %s\n", ium_mesg->source);
	  
    CI.flags.bits.cgbType = XX_CGBTYPE;
    CI.aEndCoord = CI.bEndCoord = -1;
    simLength = CI.bpLength.mean;

    // See if this is a repeat, or we can pin it down to an interval
#ifdef AS_ENABLE_SOURCE
    type = strstr(ium_mesg->source,"gen> ");
    if(type){
      type += 5;
      if(!strncmp(type,"uu",2) || !strncmp(type,"@@",2)){
	CI.flags.bits.cgbType = (unsigned int)UU_CGBTYPE;
      }else if(!strncmp(type,"ru",2)){
	CI.flags.bits.cgbType = (unsigned int)RU_CGBTYPE;
      }else if(!strncmp(type,"rr",2)){
	CI.flags.bits.cgbType = (unsigned int)RR_CGBTYPE;
      }else if(!strncmp(type,"ur",2)){
	CI.flags.bits.cgbType = (unsigned int)UR_CGBTYPE;
      }

      if((interval = strstr(ium_mesg->source,"["))){
	//	    fprintf(stderr,"* interval = %s\n", interval);
	result = sscanf(interval + 1," " F_COORD "," F_COORD,
			&CI.aEndCoord, &CI.bEndCoord);
	simLength = abs(CI.aEndCoord - CI.bEndCoord);
      }else{
	CI.aEndCoord = CI.bEndCoord = -1;
	simLength = CI.bpLength.mean;
      }
    }
#endif
  }

  if(ium_mesg->coverage_stat >= GlobalData->cgbUniqueCutoff){
    if(length < CGW_MIN_DISCRIMINATOR_UNIQUE_LENGTH ||
       ium_mesg->num_frags < CGW_MIN_READS_IN_UNIQUE){
      ShortDiscriminatorUniques++;
    }else{
      DiscriminatorUniques++;
    }
  }

  
  {
    int isUnique = FALSE;
    if(ium_mesg->coverage_stat >= GlobalData->cgbUniqueCutoff &&
       length >= CGW_MIN_DISCRIMINATOR_UNIQUE_LENGTH &&
       ium_mesg->num_frags >= CGW_MIN_READS_IN_UNIQUE){
      // microhetScore is actually the probability of the sequence
      // being UNIQUE, based on microhet considerations.
      // Falling below threshhold makes something a repeat.
      if( CI.microhetScore < GlobalData->cgbMicrohetProb){
	if(ium_mesg->coverage_stat < GlobalData->cgbApplyMicrohetCutoff){
	  fprintf(stderr,"* CI " F_CID " with astat: %g classified as repeat based on microhet unique prob of %g < %g\n",
		  CI.id, ium_mesg->coverage_stat, CI.microhetScore, GlobalData->cgbMicrohetProb);
	  isUnique = FALSE;
	  if(CI.flags.bits.cgbType == XX_CGBTYPE)
	    CI.flags.bits.cgbType = RR_CGBTYPE;
	  CI.type = UNRESOLVEDCHUNK_CGW;
	}else{
	  isUnique = TRUE;
	  fprintf(stderr,"* WARNING: CI " F_CID " with coverage %g WOULD HAVE BEEN classified as repeat based on microhet unique prob of %g < %g\n",
		  CI.id, ium_mesg->coverage_stat, CI.microhetScore, GlobalData->cgbMicrohetProb);
	}
      }else{
	isUnique = TRUE;
      }
    }else{
      isUnique = FALSE;
    }

    if(isUnique){
      ScaffoldGraph->numDiscriminatorUniqueCIs++;
      CI.flags.bits.isUnique = 1;
      CI.type = DISCRIMINATORUNIQUECHUNK_CGW;
      if(CI.flags.bits.cgbType == XX_CGBTYPE)
	CI.flags.bits.cgbType = UU_CGBTYPE;
    }else{
      CI.flags.bits.isUnique = 0;
      CI.type = UNRESOLVEDCHUNK_CGW;
      if(CI.flags.bits.cgbType == XX_CGBTYPE)
	CI.flags.bits.cgbType = RR_CGBTYPE;
    }
  }

  CI.flags.bits.smoothSeenAlready = FALSE;
  CI.flags.bits.isCI = TRUE;
  CI.flags.bits.isChaff = FALSE;


  if( ! sequenceOnly ) {
      CDS_CID_t extremalA = NULLINDEX;
      CDS_CID_t extremalB = NULLINDEX;
      CDS_COORD_t minOffset = CDS_COORD_MAX;
      CDS_COORD_t maxOffset = CDS_COORD_MIN;
	  
      /* Determine extremal fragments so we can label the fragments */

      for(cfr = 0; cfr < ium_mesg->num_frags; cfr++){
	IntMultiPos *cfr_mesg = ium_mesg->f_list + cfr;
	CDS_COORD_t end = MAX( cfr_mesg->position.end, cfr_mesg->position.bgn);
	CDS_COORD_t beg = MIN( cfr_mesg->position.end, cfr_mesg->position.bgn);
	    
	if(minOffset > beg){
	  minOffset = beg;
	  extremalA = cfr;
	}
	if(maxOffset < end){
	  maxOffset = end;
	  extremalB = cfr;
	}
      }
	  
	  
      for(cfr = 0; cfr < ium_mesg->num_frags; cfr++){
	CIFragT        cifrag;
	InfoByIID  info, *old_info;
	CDS_CID_t fragid = GetNumCIFragTs(ScaffoldGraph->CIFrags);
	IntMultiPos *cfr_mesg = ium_mesg->f_list + cfr;

	cifrag.iid      = cfr_mesg->ident;
        cifrag.mateOf   = NULLINDEX;
        cifrag.dist     = 0;
	cifrag.cid      = ium_mesg->iaccession;
	cifrag.CIid     = ium_mesg->iaccession;

	// These get set in UpdateNodeFragments, called below
        cifrag.offset5p.mean      = 0.0;
        cifrag.offset5p.variance  = 0.0;
        cifrag.offset3p.mean      = 0.0;
        cifrag.offset3p.variance  = 0.0;

	cifrag.contigID                 = NULLINDEX;
        cifrag.contigOffset5p.mean      = 0.0;
        cifrag.contigOffset5p.variance  = 0.0;
        cifrag.contigOffset3p.mean      = 0.0;
        cifrag.contigOffset3p.variance  = 0.0;

        cifrag.type      = cfr_mesg->type;
        cifrag.label     = AS_SINGLETON;

	cifrag.flags.bits.hasInternalOnlyCILinks     = FALSE; // set in CreateCIEdge
	cifrag.flags.bits.hasInternalOnlyContigLinks = FALSE; // set in CreateCIEdge
	cifrag.flags.bits.isPlaced                   = FALSE;
	cifrag.flags.bits.isSingleton                = FALSE;
	cifrag.flags.bits.isChaff                    = FALSE;
        cifrag.flags.bits.innieMate                  = FALSE;
        cifrag.flags.bits.hasMate                    = FALSE;
        cifrag.flags.bits.linkType                   = AS_UNKNOWN;
	cifrag.flags.bits.edgeStatus                 = INVALID_EDGE_STATUS;
        cifrag.flags.bits.mateDetail                 = UNASSIGNED_MATE;

        //  Singleton chunks are chaff; singleton frags are chaff
        //  unless proven otherwise
        //
        if (ium_mesg->num_frags < 2) {
          CI.flags.bits.isChaff         = TRUE;
          cifrag.flags.bits.isSingleton = TRUE;
          cifrag.flags.bits.isChaff     = TRUE;
	}

	info.fragIndex   = fragid;
	info.set         = TRUE;

	// Check to see if we've already seen this fragment by IID!!!!
	old_info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, cifrag.iid);
	if(old_info && old_info->set){
	  CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags, old_info->fragIndex);
	  fprintf(stderr,"*** FATAL ERROR:  Fragment with IID " F_CID " appears more than once with id " F_CID " and " F_CID "\n",
                  cifrag.iid, old_info->fragIndex, fragid);
	  fprintf(stderr,"***               First appearance was in unitig " F_CID ", currently found in unitig " F_CID "\n",
		  frag->cid, cifrag.cid);
	  exit(1);
	}

	SetInfoByIID(ScaffoldGraph->iidToFragIndex, cifrag.iid, &info);
	    
	// Collect read stats
        if (AS_FA_READ(cfr_mesg->type)) { 
	  totalReadFrags++;

	  if(CI.flags.bits.isUnique)
	    inUniqueReadFrags++;
	  else
	    inRepeatReadFrags++;

	  if(ium_mesg->num_frags <= 2)
            inTeenyUnitigReadFrags++;

          if(ium_mesg->num_frags < 2)
	    inSingletonUnitigReadFrags++;

	  if(cfr == extremalA || cfr == extremalB)
	    onEndReadFrags++;
	}

        //else if(AS_FA_SHREDDED(cfr_mesg->type)){ 
	//  CI.flags.bits.includesFinishedBacFragments = TRUE; 
        //}
	   
        AppendCIFragT(ScaffoldGraph->CIFrags, &cifrag);
      }
    }

  // Insert the Chunk Instance
  SetChunkInstanceT(ScaffoldGraph->CIGraph->nodes, CI.id, &CI);

  // Mark all frags as being members of this CI, and set their offsets within
  // the CI
  if( ! sequenceOnly )
    UpdateNodeFragments(ScaffoldGraph->CIGraph,CI.id, CI.type == DISCRIMINATORUNIQUECHUNK_CGW, TRUE ); // mark unitigs and contigs
}




void
LoadDistData(void) {
  int32 numDists = getNumGateKeeperLibraries(ScaffoldGraph->gkpStore);
  CDS_CID_t i;
  
  for(i = 1; i <= numDists; i++){
    DistT dist;
    GateKeeperLibraryRecord  *gkpl = getGateKeeperLibrary(ScaffoldGraph->gkpStore, i);

    dist.mu             = gkpl->mean;
    dist.sigma          = gkpl->stddev;
    dist.numSamples     = 0;
    dist.min            = CDS_COORD_MAX;
    dist.max            = CDS_COORD_MIN;
    dist.bnum           = 0;
    dist.bsize          = 0;
    dist.histogram      = NULL;
    dist.lower          = dist.mu - CGW_CUTOFF * dist.sigma;
    dist.upper          = dist.mu + CGW_CUTOFF * dist.sigma;
    dist.numReferences  = 0;
    dist.numBad         = 0;

    fprintf(GlobalData->stderrc,"* Loaded dist "F_UID","F_CID" (%g +/- %g)\n",
            gkpl->libraryUID, i, dist.mu, dist.sigma);

    SetDistT(ScaffoldGraph->Dists, i, &dist);
  }
}