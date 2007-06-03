
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
static const char CM_ID[] = "$Id: AS_CGW_main.c,v 1.40 2007-05-22 17:04:25 granger_sutton Exp $";


static const char *usage = 
"usage: %s [options] -g <GatekeeperStoreName> -o <OutputPath> <InputCGB.ext>\n"
"\n"
"   [-a]           align overlaps  (default)\n"
"   [-b]           don't ignore UOM between contained     (default)\n"
"   [-c]           Generate checkpoints\n"
"   [-d]           DumpGraphs\n"
"   [-e <thresh>]  Microhet score probability cutoff\n"
"   [-g]           gkp Store path (required)\n"
"   [-h]           Don't fail on merge alignment failure\n"
"   [-i <thresh>]  Set max coverage stat for microhet determination of non-uniqueness (default -1)\n"
"   [-j <thresh>]  Set min coverage stat for definite uniqueness\n"
"   [-k <thresh>]  Set max coverage stat for possible uniqueness\n"
"   [-l <maxdegree> ]\n"
"   [-m <minSamplesForOverride>]   For ComputeMatePairStatisticsRestricted, default is 1000\n"
"   [-n] <scaffnum>  Starting scaffold number for stones\n"
"   [-o]           Output Name (required)\n"
"   [-q <cutoff>]  Transquality cutoff\n"
"   [-r <repeatRezlevel>]\n"
"   [-s <stoneLevel> ]\n"
"   [-t]           Don't ignore UOM Transchunks\n"
"   [-T]           ignore UOM transchunks\n"
"   [-u]           don't ignore UOM containments     (default)\n"
"   [-U]           ignore UOM containments\n"
"   [-v]           verbose\n"
"   [-w <walkLevel> ]\n"
"   [-x]           Dump stats on checkpoint load\n"
"   [-y]           Turn off Check for -R < -N, for restarting small assemblies\n"
"   [-z]           Turn on Check for Repeat Branch Pattern (demotes some unique unitigs to repeat)\n"
"\n"
"   [-A]           don't align overlaps (quality values reflect bayesian)\n"
"   [-B]           ignore UOM between contained\n"
"   [-C]           Don't cleanup scaffolds\n"
"   [-D <debugLevel>]\n"
"   [-E]           output overlap only contig edges\n"
"   [-G]           Don't generate output (cgw or cam)\n"
"   [-H]           fail on merge alignment failure   (default)\n"
"   [-I]           ignore chaff unitigs\n"
"   [-J]           annotate celamy output with Bactig content (default = false)\n"
"   [-M]           don't do interleaved scaffold merging\n"
"   [-N <ckp>]     restart from checkpoint location ckp\n"
"   [-O [icC]]     i = immediate, c = continue, C = celamy output only\n"
"   [-P]                            Proto Output\n"
"   [-R <ckp>]     restart from checkpoint file ckp\n"
"   [-S]           Walk scaffolds smallest to biggest (biggest first is default)\n"
"   [-W <startWalkFromScaffold> ]\n"
"   [-X <Estimated number of nodes>]\n"
"   [-Y <Estimated number of edges>]\n"
"   [-Z]           Don't demote singleton scaffolds\n"
"   [-4]           Allow forced fragments inside consensus ... useful for high-error runs\n"
"\n"
"CGBInputFiles: The file with new IUM,OUM, etc records to process.\n"
"\n"
"Checkpoint File: File named <outputName>.ckp.n\n"
"\n"
"The Chunk Graph Walker processes the input file, reporting errors and\n"
"warnings.  Any message that causes an error warning is output\n"
"to stderr, along with the associated warning/error messages.\n"
"\n"
"Copious log info is sent to <inputFile>.cgwlog\n"
"Celamy output is sent to <inputFile>.cam\n"
"\n"
"Opens ALL [<InputFileName>.<ext>]* to read input\n"
"\n"
"Writes diagnostic output to <OutputPath>.cgwlog\n"
"Writes multiAlignments to <OutputPath>.SeqStore\n"
"Writes output to <OutputPath>.cgw\n";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include <math.h>
#include <time.h>

#ifdef X86_GCC_LINUX
#include <fpu_control.h>
#endif

#define DEBUG_MERGE_SCAF FALSE

#undef INSTRUMENT_CGW
#undef CHECK_CONTIG_ORDERS
#undef CHECK_CONTIG_ORDERS_INCREMENTAL

#ifdef CHECK_CONTIG_ORDERS_INCREMENTAL
#define POPULATE_COC_HASHTABLE 1
#else
#define POPULATE_COC_HASHTABLE 0
#endif

//  If -1, do not test or fix edges.  If 0, test but do not fix.  If
//  1, test and fix edges.
//
#define FIX_CONTIG_EDGES -1

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "UtilsREZ.h"
#include "AS_UTL_timer.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "Output_CGW.h"
#include "SplitChunks_CGW.h"
#include "GreedyOverlapREZ.h"
#include "CommonREZ.h"
#include "RepeatRez.h"
#include "FbacREZ.h"
#include "Stats_CGW.h"
#include "AS_ALN_forcns.h"
#include "Instrument_CGW.h"
#include "AS_CGW_EdgeDiagnostics.h"
#include "Checkpoints_CGW.h"

extern int allow_forced_frags;

void DemoteUnitigsWithRBP(FILE *stream, GraphCGW_T *graph);
void RemoveSurrogateDuplicates(void);

FILE *  File_Open (const char * Filename, const char * Mode, int exitOnFailure);


int main(int argc, char *argv[]){
  int checkRepeatBranchPattern = FALSE;
  int preMergeRezLevel = -1;
  int generateOutput = 1;
  int annotateUnitigs = 0;
  int doInterleavedScaffoldMerging = 1;
  int debugLevel = 0;
  int repeatRezLevel = 0;
  int stoneLevel = 0;
  int startScaffoldWalkFrom = NULLINDEX;
  int starting_stone_scaffold = 0;
  int walkScaffoldsBiggestFirst = TRUE;
  int walkLevel = 0;
  int numNodes = 1024;
  int numEdges = 1024;
  int minSamplesForOverride = 1000;
  int ignoreUOMs = 1;
  int ignoreUOMContainments = 0; // -u -- turn off with -U
  int ignoreUOMContainStack = 1; //
  int ignoreUOMBetweenContained = 1; // -b turn off with -B
  int ignoreUOMTranschunk = 1; // -t turn off with -T
  int ignoreChaffUnitigs = 0;
  int failOn_NoOverlapFound = 0; 
  int geneOutput = 1;  // output simulated coordinates
  int checkPoint = 0;
  int dumpGraphs = 0;
  int outputOverlapOnlyContigEdges = 0;
  int doRezOnContigs = 1;
  int immediateOutput = 0;
  int camFileOnly = 0;
  int performCleanupScaffolds = 1;
  int alignOverlaps = 1; // if true (set to true by -a, false -A) recompute all overlaps supplied by cgb
  int dumpStatsOnRestart = FALSE;
  int32 restartFromCheckpoint = NULLINDEX;
  int32 restartFromLogicalCheckpoint = NULLINDEX;
  int demoteSingletonScaffolds = TRUE;
  Global_CGW *data;
  FILE *infp = NULL;
  float transQualityCutoff = 0.1; // quality cutoff for TransChunkEdges
  float cgbMicrohetProb = 1.e-05;      // scores less than this are considered repeats
  float cgbApplyMicrohetCutoff = -1; // This basically turns it off, unless enabled
  float cgbUniqueCutoff = CGB_UNIQUE_CUTOFF; 
  float cgbDefinitelyUniqueCutoff = CGB_UNIQUE_CUTOFF; 
  int maxDegree = 3; // maximum edges to keep for 'nonUnique' nodes
  int maxDegreeUnique = 30; // maximum edges to keep for 'Unique' nodes
  int setGatekeeperStore = 0;
  char *outputPath = NULL;
  int dumpScaffoldSnapshots = 0;
  int checkpointChecker = 1;
 
#ifdef X86_GCC_LINUX
  /*
  ** Set the x86 FPU control word to force double
  ** precision rounding rather than `extended'
  ** precision rounding. This causes base
  ** calls and quality values on x86 GCC-Linux
  ** (tested on RedHat Linux) machines to be
  ** identical to those on IEEE conforming UNIX
  ** machines.
  */
  fpu_control_t fpu_cw;

  fpu_cw = ( _FPU_DEFAULT & ~_FPU_EXTENDED ) | _FPU_DOUBLE;

  _FPU_SETCW( fpu_cw );
#endif
 
  fprintf( stderr, "Version: %s",CM_ID);
#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
  ContigOrientChecker * coc;
  coc = CreateContigOrientChecker();
  assert(coc != NULL);
#endif
  
  GlobalData  = data = CreateGlobal_CGW();
  GlobalData->stderrc = stderr;
  GlobalData->aligner=DP_Compare;

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    char  * p;

    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,
                                    "abcde:f:g:hi:j:k:l:m:n:o:p:q:r:s:tuvw:xyz"
                                    "ABCD:EFGHIJK:L:N:MO:PQR:STUV:W:X:Y:Z"
				    "4" )) != EOF)){
      switch(ch) {
        case 'O':
          fprintf(GlobalData->stderrc,"* Immediate output specified: arg %s\n", optarg);
          switch(*optarg){
            case 'C':
              immediateOutput = 3;
              camFileOnly = 1;
              generateOutput = 0;
              break;
            case 'i':
              immediateOutput = 2;
              break;
            case 'c':
              immediateOutput = 1;
              break;
            case 'm': //MergeScaffolds no cgw file
              immediateOutput = -1;
              camFileOnly = 1;
              generateOutput = 0;
              break;
            default:
              break;
          }
        case 'a':
          alignOverlaps = 1;
          break;
        case 'A':
          alignOverlaps = 0;
          break;
        case 'C':
          performCleanupScaffolds = 0;
          fprintf(GlobalData->stderrc,"* -C ==> No Cleanup Scaffolds!\n");
          fflush(stderr);
          break;
        case 'g':
          strcpy( data->Gatekeeper_Store_Name, argv[optind - 1]);
          setGatekeeperStore = 1;
          break;	  
        case 'h':
          failOn_NoOverlapFound = 0;
          break;
        case 'G':
          generateOutput = 0;
          break;
        case 'H':
          failOn_NoOverlapFound = 1;
          break;
        case 'I':
          ignoreChaffUnitigs = 1;
          break;
        case 'J':
          annotateUnitigs = 1;
          break;
        case 'M':
          doInterleavedScaffoldMerging = 0;
          break;
        case 'v':
          ignoreUOMs = 0;
          break;
        case 'V':
          ignoreUOMs = 1;
          break;
        case 'Q':
          fprintf(stderr,"* -Q  dumpScaffoldSnapshots == 1\n");
          dumpScaffoldSnapshots = 1;;
          break;
        case 'u':
          ignoreUOMContainments = 0;
          break;
        case 'U':
          ignoreUOMContainments = 1;
          break;
        case 't':
          ignoreUOMTranschunk = 0;
          break;
        case 'T':
          ignoreUOMTranschunk = 1;
          break;
        case 'b':
          ignoreUOMBetweenContained = 0;
          break;
        case 'B':
          ignoreUOMBetweenContained = 1;
          break;
        case 'K':
          switch(optarg[0]) {
            case 'D':
              GlobalData->aligner=DP_Compare;
              break;
            case 'L':
              GlobalData->aligner=Local_Overlap_AS_forCNS;
              break;
            case 'A':
              GlobalData->aligner=Affine_Overlap_AS_forCNS;
              break;
            default:
              GlobalData->aligner=DP_Compare;
              fprintf(stderr,"Unrecognized value for option -%c (%s)",optopt,optarg);
          } 
          break;
        case 'q':
          transQualityCutoff = (float)atof(optarg);
          fprintf(GlobalData->stderrc,"* transQualityCutoff set to %f\n", transQualityCutoff);
          break;
        case 'i':
          cgbApplyMicrohetCutoff = (float)atof(optarg);
          fprintf(GlobalData->stderrc,"* cgbApplyMicrohetCutoff set to %f\n", cgbApplyMicrohetCutoff);
          break;
        case 'j':
          cgbUniqueCutoff = (float)atof(optarg);
          fprintf(GlobalData->stderrc,"* cgbUniqueCutoff set to %f\n", cgbUniqueCutoff);
          break;
        case 'k':
          cgbDefinitelyUniqueCutoff = (float)atof(optarg);
          fprintf(GlobalData->stderrc,"* cgbDefinitelyUniqueCutoff set to %f\n", cgbDefinitelyUniqueCutoff);
          break;
        case 'l':
          maxDegree = atoi(optarg);
          fprintf(GlobalData->stderrc,"* maxDegree set to %d\n", maxDegree);
          break;
        case 'm':
          minSamplesForOverride = atoi(optarg);
          fprintf(GlobalData->stderrc,"* minSamplesForOverride set to %d\n", minSamplesForOverride);
          break;
        case 'n':
          starting_stone_scaffold = atoi(optarg);
          fprintf (GlobalData -> stderrc,"* Starting Scaffold for Stones set to %d\n",
                   starting_stone_scaffold);
          break;
        case 'L':
          maxDegreeUnique = atoi(optarg);
          fprintf(GlobalData->stderrc,"* maxDegreeUnique set to %d\n", maxDegreeUnique);
          break;
        case 'X':
          numNodes = atoi(optarg);
          fprintf(GlobalData->stderrc,"* numNodes set to %d\n", numNodes);
          break;
        case 'y':
          checkpointChecker = 0;
          break;
        case 'Y':
          numEdges = atoi(optarg);
          fprintf(GlobalData->stderrc,"* numEdges set to %d\n", numEdges);
          break;

        case 'r':
          repeatRezLevel = atoi(optarg);
          fprintf(GlobalData->stderrc,"* repeatRezLevel set to %d\n", repeatRezLevel);
          break;
        case 's':
          stoneLevel = atoi(optarg);
          fprintf(GlobalData->stderrc,"* stoneLevel set to %d\n", stoneLevel);
          break;
        case 'S':
          walkScaffoldsBiggestFirst = FALSE;
          fprintf(GlobalData->stderrc,"* walkScaffoldsBiggestFirst set to FALSE\n");
          break;
        case 'w':
          walkLevel = atoi(optarg);
          fprintf(GlobalData->stderrc,"* walkLevel set to %d\n", walkLevel);
          break;
        case 'W':
          startScaffoldWalkFrom = atoi(optarg);
          fprintf(GlobalData->stderrc,"* startScaffoldWalkFrom scaffold %d\n", startScaffoldWalkFrom);
          break;
        case 'D':
          debugLevel = atoi(optarg);
          fprintf(GlobalData->stderrc,"* debugLevel set to %d\n", debugLevel);
          break;
        case 'o':
          outputPath = strdup(optarg);
          break;
        case 'p':
          preMergeRezLevel = atoi(optarg);
          if(preMergeRezLevel < 0){
            fprintf(GlobalData->stderrc,"* pre-first-scaffold-merge repeatRezLevel must be >= 0\n");
            fprintf(stderr,"* pre-first-scaffold-merge repeatRezLevel must be >= 0\n");
            exit(-99);
          }
          fprintf(GlobalData->stderrc,"* pre-first-scaffold-merge repeatRezLevel set to %d\n", preMergeRezLevel);
          break;
        case 'c':
          checkPoint = 1;
          break;
        case 'd':
          dumpGraphs = 1;
          fprintf(GlobalData->stderrc,"* dumpGraphs set to 1\n");
          break;
        case 'e':
          cgbMicrohetProb = atof(optarg);
          fprintf(GlobalData->stderrc,"* Microhet probability threshhold set to %g\n",cgbMicrohetProb);
          break;
        case 'E':
          outputOverlapOnlyContigEdges = 1;
          fprintf(GlobalData->stderrc,"* Output overlap only  contig edges set to TRUE\n");
          break;
        case 'P':
          fprintf(stderr, "-P default.\n");
          break;
        case 'R':
          restartFromCheckpoint = atoi(optarg);
          fprintf(GlobalData->stderrc,"* Restarting from checkpoint %d\n", restartFromCheckpoint);
          break;
        case 'N':
          restartFromLogicalCheckpoint = atoi(optarg);
          fprintf(GlobalData->stderrc,"* Restarting from logical checkpoint %d\n", restartFromLogicalCheckpoint);
          break;
        case 'x':
          dumpStatsOnRestart = TRUE;
          fprintf(GlobalData->stderrc,"* Dumping Stats on restart\n");
          break;
        case 'Z':
          fprintf(GlobalData->stderrc,"* DON'T demote singleton scaffolds\n");
          demoteSingletonScaffolds = FALSE;
          break;
        case 'z':
          fprintf(GlobalData->stderrc,"* Check for Repeat Branch Pattern\n");
          checkRepeatBranchPattern = TRUE;
          break;

        case '4':
          allow_forced_frags=1;
          break;

        case '?':
          fprintf(GlobalData->stderrc,"Unrecognized option -%c",optopt);
        default :
          errflg++;
      }
    }

    if ((argc - optind < 1 ) || (setGatekeeperStore == 0) || (outputPath == NULL)) {
      fprintf(GlobalData->stderrc,"* argc = %d optind = %d setGatekeeperStore = %d outputPath = %s\n",
              argc, optind,setGatekeeperStore, outputPath);
      fprintf(GlobalData->stderrc, usage, argv[0]);
      exit (EXIT_FAILURE);
    }


    if(cgbDefinitelyUniqueCutoff < cgbUniqueCutoff){
      cgbDefinitelyUniqueCutoff = cgbUniqueCutoff;
      fprintf(GlobalData->stderrc,"* cgbDefinitelyUniqueCutoff set to %f\n", cgbUniqueCutoff);
    }
    data->annotateUnitigs = annotateUnitigs;
    data->doInterleavedScaffoldMerging = doInterleavedScaffoldMerging;
    data->dumpScaffoldSnapshots = dumpScaffoldSnapshots;
    data->maxSequencedbSize = MAX_SEQUENCEDB_SIZE;
    data->maxSequencedbCacheSize = MAX_SEQUENCEDB_CACHE_SIZE;
    data->maxDegreeUnique = maxDegreeUnique;
    data->maxDegree = maxDegree;
    data->transQualityCutoff = transQualityCutoff;
    data->cgbUniqueCutoff = cgbUniqueCutoff;
    if(preMergeRezLevel >= 0){
      data->repeatRezLevel = preMergeRezLevel;
    } else {
      data->repeatRezLevel = repeatRezLevel;
    }
    data->stoneLevel     = stoneLevel;
    data->walkLevel      = walkLevel;
    data->debugLevel     = debugLevel;
    data->failOn_NoOverlapFound = failOn_NoOverlapFound;
    data->walkScaffoldsBiggestFirst = walkScaffoldsBiggestFirst;
    data->ignoreChaffUnitigs = ignoreChaffUnitigs;
    data->performCleanupScaffolds = performCleanupScaffolds;
    data->cgbUniqueCutoff = cgbUniqueCutoff;
    data->cgbApplyMicrohetCutoff = cgbApplyMicrohetCutoff;
    data->cgbDefinitelyUniqueCutoff = cgbDefinitelyUniqueCutoff;
    data->cgbMicrohetProb = cgbMicrohetProb;
    data -> starting_stone_scaffold = starting_stone_scaffold;

    if(optind < argc) 
      {
        char  filepath[2048];

	strcpy(data->Input_File_Name, argv[optind]);
	infp = File_Open (data->Input_File_Name, "r", TRUE);     // frg file
	AssertPtr(infp);

	strcpy(data->File_Name_Prefix, outputPath);

        if (generateOutput) {
          sprintf(filepath, "%s.cgw", outputPath);
          data->cgwfp = File_Open(filepath, "w", TRUE);

          sprintf(filepath, "%s.cgw_contigs", outputPath);
          data->ctgfp = File_Open(filepath, "w", TRUE);

          sprintf(filepath, "%s.cgw_scaffolds", outputPath);
          data->scffp = File_Open(filepath, "w", TRUE);
        }
        
        sprintf(filepath, "%s.timing", outputPath);
        data->timefp = File_Open(filepath, "a", TRUE); // timing file

	{
	  time_t t;
	  t = time(0);
	  fprintf(data->timefp,"\n\n>>>>*************************************************************************<<<<\n");
	  fprintf(data->timefp,"* Restarting from checkpoint %d at time %s\n", restartFromCheckpoint,ctime(&t));
	}

	//	optind++;
      }
    /* End of command line parsing */
  }

  if(!doRezOnContigs && (walkLevel > 0 || stoneLevel > 0)){
    fprintf(GlobalData->stderrc,"#### Stone Throwing and Gap Walking can only be run in conjunction with the -C option ####\n");
    fprintf(GlobalData->stderrc,".....exiting.....\n");
    exit(1);
  }

  if(immediateOutput > 0 &&
     restartFromCheckpoint <= CHECKPOINT_AFTER_BUILDING_CIGRAPH){
    fprintf(GlobalData->stderrc,"* Can't generate output unless -R n, n > %d is specified...ignored\n",
	    CHECKPOINT_AFTER_BUILDING_CIGRAPH);
    immediateOutput = 0;
  }

  fprintf(GlobalData->stderrc,"* Scaffolding will be performed on %s\n",
	  (doRezOnContigs?"Contigs":"CIs"));

  if(ignoreUOMs)
    fprintf(GlobalData->stderrc,"* UOMs will be ignored *\n");
  if(ignoreUOMContainments)
    fprintf(GlobalData->stderrc,"* UOMs specifying containment will be ignored *\n");
  if(ignoreUOMBetweenContained)
    fprintf(GlobalData->stderrc,"* UOMs specifying BetweenContained will be ignored *\n");
  if(ignoreUOMTranschunk)
    fprintf(GlobalData->stderrc,"* UOMs specifying Transchunk will be ignored *\n");
  if(ignoreUOMContainStack)
    fprintf(GlobalData->stderrc,"* UOMs specifying ContainStack will be ignored *\n");
  if(alignOverlaps)
    fprintf(GlobalData->stderrc,"* UOMs will be recomputed within cgw *\n");
  else
    fprintf(GlobalData->stderrc,"* UOMs in .cgb file contain quality value alignments and will NOT be recomputed *\n");

  ProcessInputADT(data, infp, argc, argv); // Handle ADT

  data->saveCheckPoints = checkPoint;
  data->outputCalculatedOffsets = (geneOutput == 0);

  if(checkpointChecker && restartFromCheckpoint < restartFromLogicalCheckpoint &&
     restartFromLogicalCheckpoint !=  CHECKPOINT_BEFORE_FINAL_CLEANUP){
    fprintf(GlobalData->stderrc,"* Logical Checkpoint (%d)  MUST be no greater than restart checkpoint (%d)\n", restartFromLogicalCheckpoint, restartFromCheckpoint);
    exit(1);
  }

  if(restartFromLogicalCheckpoint == NULLINDEX){
    restartFromLogicalCheckpoint = restartFromCheckpoint;
  }


  if(restartFromLogicalCheckpoint < CHECKPOINT_AFTER_READING_INPUT){
    // create the checkpoint from scratch
    ScaffoldGraph = CreateScaffoldGraph(doRezOnContigs, data->File_Name_Prefix, numNodes, numEdges);
    
    // If we are ingoring these, just return, don't add any edges
    ScaffoldGraph->ignoreUOMs = ignoreUOMs;
    ScaffoldGraph->ignoreUOMContains = ignoreUOMContainments;
    ScaffoldGraph->ignoreUOMBetweenContained = ignoreUOMBetweenContained;
    ScaffoldGraph->ignoreUOMTranschunk = ignoreUOMTranschunk;
    ScaffoldGraph->ignoreUOMContainStack = ignoreUOMContainStack;

    ScaffoldGraph->alignOverlaps = alignOverlaps;
    ScaffoldGraph->doRezOnContigs = doRezOnContigs;
    ScaffoldGraph->checkPointIteration = 0;

    ProcessInput(data, optind, argc, argv);    // Handle the rest of the first file

    LoadDistData();

    CheckpointScaffoldGraph(ScaffoldGraph, CHECKPOINT_AFTER_READING_INPUT);

  } else {
    // Load the checkpoint from a file
    ScaffoldGraph = LoadScaffoldGraphFromCheckpoint(data->File_Name_Prefix,restartFromCheckpoint, !immediateOutput);

    if(dumpStatsOnRestart &&
       restartFromCheckpoint >= CHECKPOINT_AFTER_BUILDING_AND_CLEANING_SCAFFOLDS){
      char tmpBuffer[2048];

      fprintf(GlobalData->stderrc," Dumping Stats from Checkpoint ");
      sprintf(tmpBuffer,"LoadCKP%d",restartFromCheckpoint);
      GeneratePlacedContigGraphStats(tmpBuffer,0);
      GenerateScaffoldGraphStats(tmpBuffer,0);
    }
  }

  if(restartFromLogicalCheckpoint < CHECKPOINT_AFTER_BUILDING_CIGRAPH) {

    // Split chimeric unitigs
    if( restartFromLogicalCheckpoint < CHECKPOINT_AFTER_UNITIG_SPLITTING ){
      fprintf(GlobalData->stderrc,"* Splitting chimeric input unitigs\n");
      fflush(GlobalData->stderrc);
      SplitInputUnitigs(ScaffoldGraph);
      fprintf(GlobalData->timefp," Dumping checkpoint %d after SplitInputUnitigs\n",
              ScaffoldGraph->checkPointIteration);	  
      CheckpointScaffoldGraph(ScaffoldGraph, CHECKPOINT_AFTER_UNITIG_SPLITTING);
    }

    //    StatsMultiAlignStore(ScaffoldGraph->CIGraph->maStore, stderr);
    if(restartFromLogicalCheckpoint < CHECKPOINT_AFTER_BUILDING_EDGES){
      fprintf(GlobalData->stderrc,"* Calling ComputeMatePairStatisticsRestricted (UNITIG_OPERATIONS)\n");
      fflush(stderr);
      ComputeMatePairStatisticsRestricted( UNITIG_OPERATIONS, minSamplesForOverride /* update distance estimates */, 
                                           "unitig_initial");

      fprintf(data->stderrc,"** Before BUILDCIEDGES **\n");
      if(GlobalData->debugLevel > 0){
        DumpChunkInstances(data->stderrc, ScaffoldGraph, FALSE, FALSE, FALSE, FALSE);
      }

      // Allocate space for edges
      ReallocGraphEdges(ScaffoldGraph->CIGraph, numEdges);

      //    BuildCIEdges(ScaffoldGraph);
      BuildGraphEdgesDirectly(ScaffoldGraph->CIGraph);
      if(checkPoint){
        fprintf(GlobalData->stderrc," Dumping checkpoint %d after buildgraphedgesdirectly\n",
                ScaffoldGraph->checkPointIteration); 
        fprintf(GlobalData->timefp," Dumping checkpoint %d after buildgraphedgesdirectly\n",
                ScaffoldGraph->checkPointIteration);
        CheckpointScaffoldGraph(ScaffoldGraph, CHECKPOINT_AFTER_BUILDING_EDGES);
      }
    }  //  checkpoint < CHECKPOINT_AFTER_BUILDING_EDGES

    // Compute all overlaps implied by mate links between pairs of unique unitigs
    ComputeOverlaps( ScaffoldGraph->CIGraph, TRUE, alignOverlaps);  

    CheckEdgesAgainstOverlapper(ScaffoldGraph->CIGraph);
#if 0
    PropagateTandemMarks(ScaffoldGraph->CIGraph);
    if(GlobalData->debugLevel > 0){
      CheckEdgesAgainstOverlapper(ScaffoldGraph->CIGraph);
    }
#endif
    MergeAllGraphEdges(ScaffoldGraph->CIGraph, FALSE);
    CheckEdgesAgainstOverlapper(ScaffoldGraph->CIGraph);

    /* Mark some Unitigs/Chunks/CIs as repeats based on overlaps GRANGER 2/2/07 */

    if(checkRepeatBranchPattern){
      DemoteUnitigsWithRBP(GlobalData->stderrc, ScaffoldGraph->CIGraph);
    }

    /* At this Point we've constructed the CIGraph */

    /*    if(dumpGraphs)
	  DumpScaffoldGraph(ScaffoldGraph);*/

    if(GlobalData->debugLevel > 0){
      DumpChunkInstances(data->stderrc, ScaffoldGraph, FALSE, FALSE, FALSE, FALSE);
      CheckEdgesAgainstOverlapper(ScaffoldGraph->CIGraph);
    }

    /* Insert the guide edges into the graph AFTER merging */
    //    InsertGuideEdges(ScaffoldGraph->CIGraph); 
    {      
      /* Looks like we never keep BOTH the CIGraph overlaps and the ContigGraph overlaps */
      ChunkOverlapperT *tmp = ScaffoldGraph->ContigGraph->overlapper;
      ScaffoldGraph->ContigGraph->overlapper = ScaffoldGraph->CIGraph->overlapper;
      ScaffoldGraph->CIGraph->overlapper = tmp;
      ScaffoldGraph->RezGraph = ScaffoldGraph->ContigGraph;
    }
    /* Construct the Contigs and Contig Edges from the Chunk Instances */
    ClearCacheSequenceDB(ScaffoldGraph->sequenceDB, TRUE);
    BuildInitialContigs(ScaffoldGraph);
    if(GlobalData->debugLevel){
      CheckEdgesAgainstOverlapper(ScaffoldGraph->ContigGraph);
    }
    // Clear both sequence caches...they will reload as necessary
    fprintf(GlobalData->stderrc,"* Start Flushing sequenceDB caches!\n");
    ClearCacheSequenceDB(ScaffoldGraph->sequenceDB, TRUE);
    ClearCacheSequenceDB(ScaffoldGraph->sequenceDB, FALSE);
    fprintf(GlobalData->stderrc,"* Done Flushing sequenceDB caches!\n");


    if(checkPoint){
      fprintf(GlobalData->timefp," Dumping checkpoint %d after BuildInitialContigs\n",
              ScaffoldGraph->checkPointIteration); 
      CheckpointScaffoldGraph(ScaffoldGraph, CHECKPOINT_AFTER_BUILDING_SCAFFOLDS);
    }

    CheckCIScaffoldTs(ScaffoldGraph);
    if(GlobalData->debugLevel > 0){
      CheckEdgesAgainstOverlapper(ScaffoldGraph->RezGraph);
    }

  }  //  checkpoint < CHECKPOINT_AFTER_BUILDING_CIGRAPH


  // Build scaffolds and do rocks

  if(immediateOutput == 0 && 
     ((restartFromLogicalCheckpoint < CHECKPOINT_AFTER_BUILDING_AND_CLEANING_SCAFFOLDS) &&
      data->repeatRezLevel > 0)){
    int skipInitialScaffolds = 0;

    //#include obsolete/rat_lbac_reactivation

    if(GlobalData->debugLevel > 0)
      DumpContigs(data->stderrc,ScaffoldGraph, FALSE);

    ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);

    skipInitialScaffolds = (restartFromLogicalCheckpoint >= CHECKPOINT_AFTER_BUILDING_SCAFFOLDS);

    fprintf(GlobalData->stderrc,"**** Running BuildScaffoldsFromFirstPriniciples ****\n");
    fprintf(GlobalData->stderrc,"**** with skipInitialScaffolds = %d\n", skipInitialScaffolds);

    BuildScaffoldsFromFirstPriniciples(ScaffoldGraph, skipInitialScaffolds);  // rocks is called inside of here

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
    ResetContigOrientChecker(coc);
    AddAllScaffoldsToContigOrientChecker(ScaffoldGraph, coc);
#endif
    
    if(GlobalData->debugLevel > 0)
      DumpCIScaffolds(GlobalData->stderrc,ScaffoldGraph, FALSE);

    //  The checkpoints are written in BuildScaffoldsFromFirstPriniciples(), either for
    //  AFTER_BUILDING_SCAFFOLDS or AFTER_BUILDING_AND_CLEANING_SCAFFOLDS, depending on
    //  (GlobalData->repeatRezLevel > 0).

  }

  //#include obsolete/compute_new_dist_estimates

  // Conservative external gap walking

  if((immediateOutput == -1) || (immediateOutput <= 1 &&
                                 (restartFromLogicalCheckpoint <= CHECKPOINT_BEFORE_CONSERVATIVE_WALKING))){
    if( GlobalData->walkLevel > 0 )
      {
        fprintf(GlobalData->stderrc,"**** Running Conservative Walking level %d ****\n",GlobalData->walkLevel);
        if(GlobalData->debugLevel > 0){
          DumpCIScaffolds(GlobalData->stderrc,ScaffoldGraph, FALSE);
          CheckEdgesAgainstOverlapper(ScaffoldGraph->RezGraph);
        }
        CheckCIScaffoldTs(ScaffoldGraph);
        StartTimerT(&data->GapWalkerTimer);
        Walk_Gaps(GlobalData, GlobalData->File_Name_Prefix, GlobalData->walkLevel, startScaffoldWalkFrom, 
                  CONSERVATIVE_WALKING_STD_DEVS);
        StopTimerT(&data->GapWalkerTimer);
        ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);
        CheckCIScaffoldTs(ScaffoldGraph);


        fprintf(GlobalData->stderrc,"**** Finished Conservative Walking level %d ****\n",GlobalData->walkLevel);
        if(GlobalData->debugLevel > 0)
          DumpCIScaffolds(GlobalData->stderrc,ScaffoldGraph, FALSE);

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
        fprintf(stderr,
                "---Checking contig orders after Walk_Gaps (1)\n\n");
        CheckAllContigOrientationsInAllScaffolds(ScaffoldGraph, coc, POPULATE_COC_HASHTABLE);
#endif
#ifdef CHECK_CONTIG_ORDERS_INCREMENTAL
        ResetContigOrientChecker(coc);
        AddAllScaffoldsToContigOrientChecker(ScaffoldGraph, coc);
#endif

        fprintf(data->timefp,"Checkpoint %d written after Conservative Walking\n",ScaffoldGraph->checkPointIteration);
        CheckpointScaffoldGraph(ScaffoldGraph, CHECKPOINT_BEFORE_CONSERVATIVE_WALKING+1);
      }
  }

  if((immediateOutput == -1) || (immediateOutput <= 1 &&
                                 (restartFromLogicalCheckpoint <= CHECKPOINT_BEFORE_1ST_SCAFF_MERGE))){
    fprintf(GlobalData->stderrc,"*** CS before SM ***\n");
    fflush(GlobalData->stderrc);

    if(GlobalData->debugLevel > 0){
      fprintf(GlobalData->stderrc,"**** BEFORE cleanupScaffolds ***\n");
      DumpCIScaffolds(GlobalData->stderrc, ScaffoldGraph, FALSE);
    }
    CleanupScaffolds(ScaffoldGraph,FALSE, NULLINDEX, FALSE);
  
    if(GlobalData->debugLevel > 0){
      fprintf(GlobalData->stderrc,"**** AFTER cleanupScaffolds ***\n");
      DumpCIScaffolds(GlobalData->stderrc, ScaffoldGraph, FALSE);
    }

    CheckCIScaffoldTs(ScaffoldGraph);

    if(GlobalData->dumpScaffoldSnapshots){
      DumpScaffoldSnapshot("PreScafMerge");
    }
    /* First we try to merge Scaffolds agressively */
    MergeScaffoldsAggressive(ScaffoldGraph, CHECKPOINT_BEFORE_1ST_SCAFF_MERGE, DEBUG_MERGE_SCAF);

    ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);
    CleanupScaffolds(ScaffoldGraph, FALSE, NULLINDEX, FALSE);

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
    fprintf(stderr,
            "---Checking contig orders after MergeScaffoldsAggressive (1)\n\n");
    CheckAllContigOrientationsInAllScaffolds(ScaffoldGraph, coc, POPULATE_COC_HASHTABLE);
#endif
#ifdef CHECK_CONTIG_ORDERS_INCREMENTAL
    ResetContigOrientChecker(coc);
    AddAllScaffoldsToContigOrientChecker(ScaffoldGraph, coc);
#endif

    fprintf(data->timefp,"Checkpoint %d written after 1st Scaffold Merge\n",ScaffoldGraph->checkPointIteration);
    CheckpointScaffoldGraph(ScaffoldGraph, CHECKPOINT_BEFORE_1ST_SCAFF_MERGE+1);
  } // No immediate output

  CheckScaffoldGraphCache(ScaffoldGraph);

  /*
    now that we are done with initial scaffold merge, we want to use the
    standard/default repeatRezLevel. Up to now, the value of preMergeRezLevel
    was in use if set on the command line
  */
  data->repeatRezLevel = repeatRezLevel;
  
  // Convert single-contig scaffolds that are marginally unique back
  // to unplaced contigs so they might be placed as stones
  //
  if (demoteSingletonScaffolds)
    DemoteSmallSingletonScaffolds();

  /* Now we throw stones */
  if  (((immediateOutput == -1) ||
        ((immediateOutput <= 1) &&
         (restartFromLogicalCheckpoint <= CHECKPOINT_BEFORE_STONES))) &&
       (GlobalData->stoneLevel > 0)) {
    fprintf(GlobalData->stderrc, "**** Running Stone Throwing level %d ****\n",
            GlobalData->stoneLevel);

    CheckCIScaffoldTs(ScaffoldGraph);

    StartTimerT(&data->StoneThrowingTimer);
    Throw_Stones(GlobalData->File_Name_Prefix, GlobalData->stoneLevel, FALSE);
    StopTimerT(&data->StoneThrowingTimer);

    CheckCIScaffoldTs(ScaffoldGraph);

    //  Throw_Stones does this now.
    //Force_Increasing_Variances();

    fprintf(GlobalData->stderrc, "**** Finished Stone Throwing level %d ****\n",
            GlobalData->stoneLevel);

    ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);

    //  Throw_Stones does this now.
    //CleanupScaffolds(ScaffoldGraph, FALSE, NULLINDEX, FALSE);

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
    fprintf(stderr, "---Checking contig orders after Throw_Stones\n\n");
    CheckAllContigOrientationsInAllScaffolds(ScaffoldGraph, coc, POPULATE_COC_HASHTABLE);
#endif

#ifdef CHECK_CONTIG_ORDERS_INCREMENTAL
    ResetContigOrientChecker(coc);
    AddAllScaffoldsToContigOrientChecker(ScaffoldGraph, coc);
#endif

    fprintf(data->timefp, "Checkpoint %d written after Stone Throwing and CleanupScaffolds\n",
            ScaffoldGraph->checkPointIteration);
    CheckpointScaffoldGraph(ScaffoldGraph, CHECKPOINT_BEFORE_STONES+1);

    GenerateLinkStats(ScaffoldGraph->CIGraph, "Stones", 0);
    GeneratePlacedContigGraphStats("Stones", 0);
    GenerateLinkStats(ScaffoldGraph->ContigGraph, "Stones", 0);
    GenerateScaffoldGraphStats("Stones", 0);
  }


  // More aggressive external gap walking

  if((immediateOutput == -1) || (immediateOutput <= 1 &&
                                 (restartFromLogicalCheckpoint <= CHECKPOINT_BEFORE_AGGRESSIVE_WALKING))){

#ifdef INSTRUMENT_CGW
    MateInstrumenter mi_before;
    ScaffoldGraphInstrumenter * sg_inst;
    
    sg_inst =
      CreateScaffoldGraphInstrumenter(ScaffoldGraph, INST_OPT_ALL_MATES);
#endif
  
    if( GlobalData->walkLevel > 0 ){
      fprintf(GlobalData->stderrc,"**** Running Aggressive Walking level %d ****\n",GlobalData->walkLevel);
      if(GlobalData->debugLevel > 0){
	DumpCIScaffolds(GlobalData->stderrc,ScaffoldGraph, FALSE);
	CheckEdgesAgainstOverlapper(ScaffoldGraph->RezGraph);
      }
      CheckCIScaffoldTs(ScaffoldGraph);
      
#ifdef INSTRUMENT_CGW
      InstrumentScaffoldGraph(ScaffoldGraph, sg_inst,
                              30000, CDS_UINT32_MAX,
                              InstrumenterVerbose1, GlobalData->stderrc);
      GetMateInstrumenterFromScaffoldGraphInstrumenter(&mi_before, sg_inst);
#endif
        
      StartTimerT(&data->GapWalkerTimer);
      Walk_Gaps(GlobalData, GlobalData->File_Name_Prefix, GlobalData->walkLevel, startScaffoldWalkFrom, 
		AGGRESSIVE_WALKING_STD_DEVS);
      StopTimerT(&data->GapWalkerTimer);

#ifdef INSTRUMENT_CGW
      InstrumentScaffoldGraph(ScaffoldGraph, sg_inst,
                              30000, CDS_UINT32_MAX,
                              InstrumenterVerbose1, GlobalData->stderrc);
      CompareMateInstrumenters(&(mi_before),
                               &(sg_inst->scaffold.mates),
                               InstrumenterVerbose1,
                               GlobalData->stderrc);
      DestroyScaffoldGraphInstrumenter(sg_inst);
#endif

      ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);
      CheckCIScaffoldTs(ScaffoldGraph);


      fprintf(GlobalData->stderrc,"**** Finished Aggressive Walking level %d ****\n",GlobalData->walkLevel);
      if(GlobalData->debugLevel > 0)
	DumpCIScaffolds(GlobalData->stderrc,ScaffoldGraph, FALSE);

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
      fprintf(stderr,
              "---Checking contig orders after Walk_Gaps (2)\n\n");
      CheckAllContigOrientationsInAllScaffolds(ScaffoldGraph, coc, POPULATE_COC_HASHTABLE);
#endif
#ifdef CHECK_CONTIG_ORDERS_INCREMENTAL
      ResetContigOrientChecker(coc);
      AddAllScaffoldsToContigOrientChecker(ScaffoldGraph, coc);
#endif
      fprintf (data -> timefp, "Checkpoint %d written after Aggressive Walking\n",
               ScaffoldGraph->checkPointIteration);
      CheckpointScaffoldGraph(ScaffoldGraph, CHECKPOINT_BEFORE_AGGRESSIVE_WALKING+1);
    }
  }


  if(immediateOutput == 0 && 
     (restartFromLogicalCheckpoint <= CHECKPOINT_BEFORE_2ND_SCAFF_MERGE)){
    fprintf(GlobalData->stderrc,"* Before Final MergeScaffoldsAggressive\n");
    CheckCIScaffoldTs(ScaffoldGraph);

#if 0
    if(GlobalData->dumpScaffoldSnapshots){
      DumpScaffoldSnapshot("PreScafMerge");
    }
#endif

    MergeScaffoldsAggressive(ScaffoldGraph, CHECKPOINT_BEFORE_2ND_SCAFF_MERGE, DEBUG_MERGE_SCAF);

    ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);
    CleanupScaffolds(ScaffoldGraph, FALSE, NULLINDEX, FALSE);

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
    fprintf(stderr,
            "---Checking contig orders after MergeScaffoldsAggressive (2)\n\n");
    CheckAllContigOrientationsInAllScaffolds(ScaffoldGraph, coc, POPULATE_COC_HASHTABLE);
#endif
#ifdef CHECK_CONTIG_ORDERS_INCREMENTAL
    ResetContigOrientChecker(coc);
    AddAllScaffoldsToContigOrientChecker(ScaffoldGraph, coc);
#endif

    fprintf(GlobalData->stderrc,"Checkpoint %d written after 2nd Aggressive Scaffold Merge\n",ScaffoldGraph->checkPointIteration);
    fprintf(data->timefp,"Checkpoint %d written after 2nd Aggressive Scaffold Merge\n",ScaffoldGraph->checkPointIteration);
    CheckpointScaffoldGraph(ScaffoldGraph, CHECKPOINT_BEFORE_2ND_SCAFF_MERGE+1);
  }  // No immediate output


  CheckScaffoldGraphCache(ScaffoldGraph);

  if  (immediateOutput == 0
       && (restartFromLogicalCheckpoint <= CHECKPOINT_BEFORE_FINAL_ROCKS)
       && GlobalData -> repeatRezLevel > 0)
    {
      const int  MAX_EXTRA_ROCKS_ITERS = 5;
      int  iter = 0, extra_rocks;

      fprintf (GlobalData -> stderrc,
               "** Running Final Round of RepeatRez (Rocks) **\n");
      do
        {
          extra_rocks = Fill_Gaps (GlobalData, GlobalData -> File_Name_Prefix,
                                   GlobalData -> repeatRezLevel, iter);
          fprintf (GlobalData -> stderrc, "Threw additional %d rocks on iter %d\n",
                   extra_rocks, iter);

          if  (extra_rocks > 0)
            CheckScaffoldGraphCache (ScaffoldGraph);
        }  while  (extra_rocks > 1 && iter < MAX_EXTRA_ROCKS_ITERS);

      if  (checkPoint)
        {
          fprintf (GlobalData -> stderrc,
                   "Checkpoint %d written after Final Rocks\n",
                   ScaffoldGraph -> checkPointIteration);
          fprintf (data->timefp,
                   "Checkpoint %d written after Final Rocks\n",
                   ScaffoldGraph -> checkPointIteration);
          CheckpointScaffoldGraph (ScaffoldGraph, CHECKPOINT_BEFORE_FINAL_ROCKS+1);
        }
    }


  if  (immediateOutput == 0
       && (restartFromLogicalCheckpoint <= CHECKPOINT_BEFORE_PARTIAL_STONES)
       && GlobalData -> stoneLevel > 0)
    {
      int  partial_stones;

      fprintf (GlobalData -> stderrc, "* Starting Partial_Stones\n");

      CheckCIScaffoldTs (ScaffoldGraph);

      StartTimerT (& data -> StoneThrowingTimer);
      partial_stones = Throw_Stones (GlobalData -> File_Name_Prefix,
                                     GlobalData -> stoneLevel,
                                     TRUE);
      StopTimerT (& data -> StoneThrowingTimer);

      CheckCIScaffoldTs (ScaffoldGraph);
      ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);

      //  We used to CleanupScaffolds() here, but Throw_Stones now
      //  does that inline.

      ClearCacheSequenceDB (ScaffoldGraph -> sequenceDB, FALSE);

      fprintf (stderr, "Threw %d partial stones\n", partial_stones);
#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
      fprintf(stderr,
              "---Checking contig orders after partial_stones\n\n");
      CheckAllContigOrientationsInAllScaffolds(ScaffoldGraph, coc, POPULATE_COC_HASHTABLE);
#endif
#ifdef CHECK_CONTIG_ORDERS_INCREMENTAL
      ResetContigOrientChecker(coc);
      AddAllScaffoldsToContigOrientChecker(ScaffoldGraph, coc);
#endif
      fprintf (GlobalData -> stderrc,
               "Checkpoint %d written after Partial Stones\n",
               ScaffoldGraph->checkPointIteration);
      fprintf (data -> timefp,
               "Checkpoint %d written after Partial Stones\n",
               ScaffoldGraph->checkPointIteration);
      CheckpointScaffoldGraph (ScaffoldGraph, CHECKPOINT_BEFORE_PARTIAL_STONES+1);

      GenerateLinkStats (ScaffoldGraph -> CIGraph, "PStones", 0);
      GeneratePlacedContigGraphStats ("PStones", 0);
      GenerateLinkStats(ScaffoldGraph -> ContigGraph, "PStones", 0);
      GenerateScaffoldGraphStats ("PStones", 0);
    }


  if  (immediateOutput == 0
       && (restartFromLogicalCheckpoint <= CHECKPOINT_BEFORE_FINAL_CONTAINED_STONES)
       && GlobalData -> stoneLevel > 0)
    {
      int  contained_stones;

      fprintf (GlobalData -> stderrc, "* Starting Final Contained Stones\n");

      CheckCIScaffoldTs (ScaffoldGraph);
      StartTimerT (& data -> StoneThrowingTimer);
      contained_stones
        = Toss_Contained_Stones (GlobalData -> File_Name_Prefix,
                                 GlobalData -> stoneLevel, 0);
      StopTimerT (& data -> StoneThrowingTimer);
      ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);
      CheckCIScaffoldTs (ScaffoldGraph);
      fprintf (GlobalData -> stderrc,
               "**** Finished Final Contained Stones level %d ****\n",
               GlobalData -> stoneLevel);
      
      CleanupScaffolds (ScaffoldGraph, FALSE, NULLINDEX, FALSE);
      ClearCacheSequenceDB (ScaffoldGraph -> sequenceDB, FALSE);

      fprintf (stderr, "Threw %d contained stones\n", contained_stones);

      // Remove copies of surrogates which are placed multiple times in the same place in a contig

      RemoveSurrogateDuplicates();

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
      fprintf(stderr,
              "---Checking contig orders after contained_stones\n\n");
      CheckAllContigOrientationsInAllScaffolds(ScaffoldGraph, coc, POPULATE_COC_HASHTABLE);
#endif
#ifdef CHECK_CONTIG_ORDERS_INCREMENTAL
      ResetContigOrientChecker(coc);
      AddAllScaffoldsToContigOrientChecker(ScaffoldGraph, coc);
#endif

      fprintf(data -> timefp, "Checkpoint %d written after Final Contained Stones\n",
              ScaffoldGraph->checkPointIteration);
      CheckpointScaffoldGraph(ScaffoldGraph, CHECKPOINT_BEFORE_FINAL_CONTAINED_STONES+1);

      GenerateLinkStats (ScaffoldGraph -> CIGraph, "CStones", 0);
      GeneratePlacedContigGraphStats ("CStones", 0);
      GenerateLinkStats(ScaffoldGraph -> ContigGraph, "CStones", 0);
      GenerateScaffoldGraphStats ("CStones", 0);
    }


  if(immediateOutput == 0 && GlobalData->walkLevel > 0 &&
     (restartFromLogicalCheckpoint <= CHECKPOINT_BEFORE_INTER_SCAFFOLD_WALKING))
    {
      int mergedScaffolds = TRUE, round = 0;
      char middleName[32];

      fprintf(GlobalData->stderrc,"* Before Inter_Scaffold_Walking\n");
	
      while ( mergedScaffolds && round < 10)
	{
	  mergedScaffolds = Inter_Scaffold_Walking();
	  if (mergedScaffolds)
            {
              sprintf( middleName, "rnd_%d", round++);
              localeCam( middleName );
            }
	}
#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
      fprintf(stderr,
              "---Checking contig orders after Inter_Scaffold_Walking\n\n");
      CheckAllContigOrientationsInAllScaffolds(ScaffoldGraph, coc, POPULATE_COC_HASHTABLE);
#endif
#ifdef CHECK_CONTIG_ORDERS_INCREMENTAL
      ResetContigOrientChecker(coc);
      AddAllScaffoldsToContigOrientChecker(ScaffoldGraph, coc);
#endif
      ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);
      fprintf(data->timefp,"Checkpoint %d written after Inter_Scaffold_Walking\n", 
              ScaffoldGraph->checkPointIteration);
      CheckpointScaffoldGraph(ScaffoldGraph, CHECKPOINT_BEFORE_INTER_SCAFFOLD_WALKING+1);
    }


  if(immediateOutput == 0 && 
     (restartFromLogicalCheckpoint <= CHECKPOINT_BEFORE_FINAL_CLEANUP)){
    fprintf(GlobalData->stderrc,"* Before CleanupFailedMerges\n");


    // Try to cleanup failed merges, and if we do, generate a checkpoint
    if(CleanupFailedMergesInScaffolds(ScaffoldGraph)){
      fprintf(data->timefp,"Checkpoint %d written after CleanupFailedMergesInScaffolds\n",ScaffoldGraph->checkPointIteration);
      CheckpointScaffoldGraph(ScaffoldGraph, -1);

      ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);
      
      // This call deletes surrogate-only contigs that failed to merge
      if( CleanupScaffolds(ScaffoldGraph, FALSE, NULLINDEX, TRUE)){
        fprintf(data->timefp,"Checkpoint %d written after DeleteAllSurrogateContigsFromFailedMerges\n",ScaffoldGraph->checkPointIteration);
        CheckpointScaffoldGraph(ScaffoldGraph, CHECKPOINT_BEFORE_FINAL_CLEANUP+1);
	
#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
        fprintf(stderr,
                "---Checking contig orders after final cleanup\n\n");
        CheckAllContigOrientationsInAllScaffolds(ScaffoldGraph, coc, POPULATE_COC_HASHTABLE);
#endif
        ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);
      }
    }
  }


  Show_Reads_In_Gaps (GlobalData -> File_Name_Prefix);

  // now recompute mate pair statistics, once on scaffolds, once on contigs, with 
  // the results on contigs being the ones that are output in OutputMateDists
  if (immediateOutput == 0) {
    fprintf(GlobalData->stderrc,"* Calling ComputeMatePairStatisticsRestricted (SCAFFOLD_OPERATIONS)\n");
    ComputeMatePairStatisticsRestricted( SCAFFOLD_OPERATIONS, minSamplesForOverride, "scaffold_final");

    fprintf(GlobalData->stderrc,"* Calling ComputeMatePairStatisticsRestricted (CONTIG_OPERATIONS)\n");
    ComputeMatePairStatisticsRestricted( CONTIG_OPERATIONS, minSamplesForOverride, "contig_final");
  }

  fprintf(GlobalData->stderrc,"* Output cam files *\n");
  GenerateCIGraph_U_Stats();
  GenerateLinkStats(ScaffoldGraph->CIGraph,"final",0);
  GeneratePlacedContigGraphStats("final",0);
  GenerateLinkStats(ScaffoldGraph->ContigGraph,"final",0);
  GenerateScaffoldGraphStats("final",0);
  GenerateSurrogateStats("final");
  //  GenerateContigAlignmentStats("final");

  CheckScaffoldGraphCache(ScaffoldGraph);

  FixupLengthsScaffoldTs(ScaffoldGraph);
  MarkMisplacedContigs();

  if(camFileOnly || generateOutput){
    CelamyCIScaffolds(data->File_Name_Prefix,ScaffoldGraph);
    CelamyAssembly(data->File_Name_Prefix);
  }

  /************* Output ***************/

  if(GlobalData->debugLevel > 0)
    CheckSmallScaffoldGaps(ScaffoldGraph);

  if(generateOutput){
    StartTimerT(&data->OutputTimer);

    MarkContigEdges();
    OutputMateDists(ScaffoldGraph);

    ComputeMatePairDetailedStatus();
    OutputFrags(ScaffoldGraph);

    // We always have multiAlignments for Unitigs
    OutputUnitigsFromMultiAligns();
    OutputUnitigLinksFromMultiAligns();

    if(GlobalData->debugLevel > 0){
      DumpContigs(data->stderrc,ScaffoldGraph, FALSE);
      DumpCIScaffolds(GlobalData->stderrc,ScaffoldGraph, FALSE);
    }

    if (!ScaffoldGraph->doRezOnContigs)
      assert(0);

    OutputContigsFromMultiAligns();
    OutputContigLinks(ScaffoldGraph, outputOverlapOnlyContigEdges);

    OutputScaffolds(ScaffoldGraph);
    OutputScaffoldLinks(ScaffoldGraph);

    StopTimerT(&data->OutputTimer);
    fprintf(GlobalData->stderrc,"* CGW Output took %g seconds\n",
	    TotalTimerT(&data->OutputTimer, NULL));
  }


  {
    long cycles;
    fprintf(GlobalData->stderrc,"* Time in Chunk Selection %g seconds (%ld calls)\n",
	    TotalTimerT(&data->ChooseChunksTimer, &cycles), cycles);
    fprintf(GlobalData->stderrc,"* Time in Consistency Check %g seconds (%ld calls)\n",
	    TotalTimerT(&data->ConsistencyCheckTimer, &cycles), cycles);
    fprintf(GlobalData->stderrc,"* Time in Update %g seconds (%ld calls)\n",
	    TotalTimerT(&data->UpdateTimer, &cycles), cycles);
    fprintf(GlobalData->stderrc,"* Time in RecomputeOffsets %g seconds (%ld calls)\n",
	    TotalTimerT(&data->RecomputeOffsetsTimer, &cycles), cycles);
    fprintf(GlobalData->stderrc,"* Time in MergeScaffolds %g seconds (%ld calls)\n",
	    TotalTimerT(&data->MergeScaffoldsTimer, &cycles), cycles);
    fprintf(GlobalData->stderrc,"* Time in Gap Fill %g seconds (%ld calls)\n",
	    TotalTimerT(&data->GapFillTimer, &cycles), cycles);
    fprintf(GlobalData->stderrc,"* Time in Gap Walking %g seconds (%ld calls)\n",
	    TotalTimerT(&data->GapWalkerTimer, &cycles), cycles);
    fprintf(GlobalData->stderrc,"* Time in Update after Gap Walking %g seconds (%ld calls)\n",
	    TotalTimerT(&data->WalkUpdateTimer, &cycles), cycles);
    fprintf(GlobalData->stderrc,"* Time in Stone Throwing %g seconds (%ld calls)\n",
	    TotalTimerT(&data->StoneThrowingTimer, &cycles), cycles);
    fprintf(GlobalData->stderrc,"* Time in Consensus %g seconds (%ld calls)\n",
	    TotalTimerT(&data->ConsensusTimer, &cycles), cycles);
  }

  DestroyScaffoldGraph(ScaffoldGraph);

  if(restartFromCheckpoint == NULLINDEX)
    DeleteGlobal_CGW(data);

  fprintf(stderr,"* Bye *\n");

  exit(0);
}
