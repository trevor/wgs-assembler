
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
/*********************************************************************
   Module:       Consensus_CNS.c
   Description:  the consensus module 
                 processes IUM or ICM messages, forming multialignment
                 and generating consensus sequence/quality
   Assumptions:  
                 
 *********************************************************************/

static const char CM_ID[] = "$Id: Consensus_CNS.c,v 1.53 2007-08-03 20:45:03 brianwalenz Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include <math.h>
#include <time.h>

#ifdef X86_GCC_LINUX
#include <fpu_control.h>
#endif

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_genericStore.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_fileIO.h"
#include "UtilsREZ.h"
#include "AS_UTL_version.h"
#include "AS_SDB_SequenceDBPartition.h"
#include "AS_ALN_forcns.h"

#include "MultiAlignStore_CNS.h"
#include "MultiAlignment_CNS.h"
#include "Globals_CNS.h"
#include "PublicAPI_CNS.h"

#define MAX_NUM_UNITIG_FAILURES 100
#define MAX_NUM_CONTIG_FAILURES 100

extern int NumColumnsInUnitigs;
extern int NumRunsOfGapsInUnitigReads;
extern int NumGapsInUnitigs;
extern int NumColumnsInContigs;
extern int NumRunsOfGapsInContigReads;
extern int NumGapsInContigs;
extern int NumAAMismatches;
extern int NumVARRecords;
extern int NumVARStringsWithFlankingGaps;
extern int NumUnitigRetrySuccess;

//  Multialignment_CNS.c options
//
extern int DUMP_UNITIGS_IN_MULTIALIGNCONTIG;
extern int VERBOSE_MULTIALIGN_OUTPUT;
extern int FORCE_UNITIG_ABUT;
extern int clear_range_to_use;


float CNS_SEQUENCING_ERROR_EST = .02; // Used to calculate '-' probability
float CNS_SNP_RATE   = 0.0003; // Used to calculate BIAS
int   CNS_HAPLOTYPES = 1;   // Used to calculate BIAS
int   CNS_USE_PUBLIC = 0;   // Used to direct basecalling to include public data
int   CNS_CALL_PUBLIC = 0;   // Used to direct basecalling to favor public data

int IntUnitigPositionCmpLeft( const IntUnitigPos *l, const IntUnitigPos *m) {
  int ltmp,mtmp;
  ltmp = (l->position.bgn<l->position.end)?l->position.bgn:l->position.end;
  mtmp = (m->position.bgn<m->position.end)?m->position.bgn:m->position.end;
  if (ltmp == mtmp) return 0;
  return (ltmp > mtmp ) ? 1 : -1;
}



int HandleDir(char * filePathAndName, char *fileName) {
  // Make sure that the file directory path of the filename exists so
  // that the file can be created.  The output is the cleaned
  // "fileName".
   char *suffix;
   char *DirName;
   char *FileName;
   suffix = strrchr(filePathAndName,(int)'/');
   if ( suffix != NULL ) {
      *suffix = '\0';
      DirName = filePathAndName; 
      if ( DirName != NULL )
        AS_UTL_mkdir(DirName);
      *suffix = '/';
      FileName = filePathAndName;
    } else {
      FileName = filePathAndName;
      DirName = NULL;
    }
    strcpy(fileName,FileName);
    return 1;
}

int32 GetUngappedSequenceLength(char *seq) {
  int32 ungappedLength=0;
  char *c;
  for(c = seq;
      *c != '\0';
      c++){

    if(*c != '-')
      ungappedLength++;
  }
  return ungappedLength;
}

static void
help_message(int argc, char *argv[])
{
    fprintf(stderr,"  Usage:\n\n"
    "  %s [-P] [-v level] [-I] [-a [DLA]] [-X expert_options] GateKeeperStoreDir [CGWStream]\n"
    "\n Standard option flags:\n"
    "    -P           Force ASCII .cns output \n"
    "    -v [0-4]     Verbose:  0 = verbose off \n"
    "                           1 = horizontal multi-alignment print in .clg\n"
    "                           2 = 'dots'     multi-alignment print in .clg\n"
    "                           3 = like 2, but dots are replaced with whitespace\n"
    "                           4 = like 1, but with unitigs in  multi-alignment print in .clg\n"
    "    -K           don't split alleles when calling consensus\n"
    "    -N           don't output variation record to .cns file\n"
    "    -w win_size  specify the size of the 'smoothing window' that will be used in consensus calling\n"
    "                 If two SNPs are located win_size or less bases apart one from another,\n"
    "                 then they will be treated as one block\n"
    "    -T secs      time threshold which, if exceeded, should trigger clean exit.\n"
    "    -S partition Use gkpStorePartition partition\n"
    "    -m           Load gkpStorePartition into memory (default reads from disk)\n"
    "    -U           Unitigs ONLY\n"
    "    -I           IUM message alignment, write to .cgi file (instead of .cns)\n"
    "                 (will also process contigs if they exist in the input file\n"
    "                 (if -I is not specified, ICMs will be processed with the assumption\n"
    "                  that IUMs were processed in a prior call.\n"
    "    -r s,e       Process only message within range [s,e)\n"
    "                 if s==0, header messages are  passed through\n"
    "                 if e==last, trailing messages are  passed through\n"
    "                 (facilitates multiprocessing in batch mode\n"
    "    -a [DLA]     Specify aligner to use should DP_Compare fail\n"
    "                 L = Local_Aligner (default)\n"
    "                 D = standard DP_Compare (will cause failed overlaps to terminate the run)\n"
    "                 A = Affine_Aligner\n"
    "    -d int       Depth of Celera coverage below which to include external data in basecalling\n"
    "                    0 (default) indicates that external data should always be used\n"
    "                    1 yields the traditional behavior, which uses external only in absence of Celera\n"
    "                  > 1 will include publice data is the Celera depth falls below the given value\n"
    "    -X           Allow 'expert' options (following)\n"
    "\n Expert option flags:\n"
    "    -D opt       Enable debugging option 'opt'.  One of 'dumpunitigs', 'verbosemultialign',\n"
    "                    and 'forceunitigabut'.  (-X not needed).\n"
    "    -R %%d       Restart from the given ICM/IUM by internal id, appending to output file\n"
    "    -i           Realign IUM messages (while processing .cgw file)\n"
    "    -q string    Override default quality call parameters\n"
    "                    string is colon separated list of the form '%%f:%%d:%%f'\n"
    "                    where first field is estimated sequencing error rate (default: .015)\n"
    "                         second field is number of sequenced haplotypes (default: 1)\n"
    "                          third field is estimated SNP rate (default: 1/1000)\n"
    "    -e #%%d      Extract only a single ICM/IUM by internal id\n"
    "    -e idfile    Extract list of ICM/IUMs by internal ids provided in idfile\n"
    "\n Arguments:\n"
    "  GateKeeperStoreDir   path to previously created Fragment Store\n"
    "  [InputStream]        previously created .cgw/.cgb file (if not specified, stdin)\n\n"
    "\n Output:\n"
    "   Creates a .cns file by default (or appends to .cns if -R is specified)\n"
    "   -I sends output to a .cgi (post-unitigging consensus) file instead.\n"
    "   -o <filename>     Overrides default output filename\n"
    "   -o -              Overrides default output filename, sending output to stdout\n"
    "   -l <filename>     Overrides default log filename\n"
    "   -l -              Overrides default log filename, sending output to stderr\n", 
    argv[0]);
    exit(1);
}

static void
OutputScores(int NumColumnsInUnitigs,        int NumRunsOfGapsInUnitigReads, 
             int NumGapsInUnitigs,           int NumColumnsInContigs, 
             int NumRunsOfGapsInContigReads, int NumGapsInContigs,
             int NumAAMismatches,          
             int NumVARRecords,              int NumVARStringsWithFlankingGaps,
             int NumUnitigRetrySuccess)
{
     fprintf(stderr, "\nNumColumnsInUnitigs      = %d\n", NumColumnsInUnitigs);
     fprintf(stderr, "NumGapsInUnitigs           = %d\n", NumGapsInUnitigs);
     fprintf(stderr, "NumRunsOfGapsInUnitigReads = %d\n", 
         NumRunsOfGapsInUnitigReads);
     fprintf(stderr, "NumColumnsInContigs        = %d\n", NumColumnsInContigs);
     fprintf(stderr, "NumGapsInContigs           = %d\n", NumGapsInContigs);
     fprintf(stderr, "NumRunsOfGapsInContigReads = %d\n", 
         NumRunsOfGapsInContigReads);
     fprintf(stderr, "NumAAMismatches            = %d\n", NumAAMismatches);
     fprintf(stderr, "NumVARRecords              = %d\n", NumVARRecords);
     fprintf(stderr, "NumVARStringsWithFlankingGaps = %d\n", 
         NumVARStringsWithFlankingGaps);
     fprintf(stderr, "NumUnitigRetrySuccess      = %d\n", NumUnitigRetrySuccess);
}

static void
writeFailure(char *OutputFileName, int std_output, GenericMesg *pmesg) {

  if (std_output) {
    fprintf(stderr, "------------------------------------------------------------\n");
    WriteProtoMesg_AS(stderr, pmesg);
    fprintf(stderr, "------------------------------------------------------------\n");
  } else {
    FILE *fout;
    char  fname[1000];

    sprintf(fname, "%s.failed", OutputFileName);

    errno = 0;
    fout = fopen(fname, "a");
    if (errno) {
      fprintf(stderr, "Failed to open '%s' for storing the failed messge: %s\n", fname, strerror(errno));
      fprintf(stderr, "------------------------------------------------------------\n");
      WriteProtoMesg_AS(stderr, pmesg);
      fprintf(stderr, "------------------------------------------------------------\n");
    } else {
      WriteProtoMesg_AS(fout,pmesg); // pass through the Unitig message and continue
      fclose(fout);
    }
  }
}



int main (int argc, char *argv[]) 
{
    char InputFileName[FILENAME_MAX];
    char OutputNameBuffer[FILENAME_MAX];
    char LogNameBuffer[FILENAME_MAX];
    char CamFileName[FILENAME_MAX];
    char SeqStoreFileName[FILENAME_MAX];
    char *sublist_file = NULL;
    FILE *cgwin;
    FILE *cam;
    FILE *sublist;
    FILE *failout;
    int sdb_version=-1;
    int sdb_partition=-1;
    int64  i;

    HashTable_AS  *tig_iids = NULL;
    HashTable_AS  *tig_iids_found = NULL;

    /**************** Process Command Line Arguments *********************/
    /* Parse the argument list using "man 3 getopt". */
    int align_ium=0;
    int no_contigs=0;
    int cgbout=0;
    int printcns=1;
    int extract=-1;
    int continue_at=0;
    int beyond=1;
    int process_sublist=0;
    int expert=0;
    int output_override=0;
    int log_override=0;
    int std_input=0;
    int std_error_log=0;
    int input_lengths=0;
    int output_lengths=0;
    int range=0;
    int partition=0;
    int in_memory=0;
    int do_rez=0;
    int noop=0;
    CNS_Options options = { CNS_OPTIONS_SPLIT_ALLELES_DEFAULT,
                            CNS_OPTIONS_SMOOTH_WIN_DEFAULT,
                            CNS_OPTIONS_MAX_NUM_ALLELES };
    static int num_unitig_failures = 0;
    static int num_contig_failures = 0;

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

    Overlap *(*COMPARE_FUNC)(COMPARE_ARGS)=Local_Overlap_AS_forCNS;
    SeqInterval tig_range;
    CNS_PrintKey printwhat=CNS_STATS_ONLY;
    int ch,errflg=0,illegal_use=0,help_flag=0,iflags=0;
    int num_of_threads = 0;
    time_t time_limit = 0, tp1 = 0;

    fprintf(stderr,"Version: %s\n",CM_ID);

    USE_SDB=0;
    USE_SDB_PART=0;

    optarg = NULL;
    terminate_cond = 1;
    partitioned=0;
    allow_forced_frags=0;
    allow_neg_hang=0;
    allow_neg_hang_retry=0;
    ALIGNMENT_CONTEXT=AS_CONSENSUS;
   
    NumColumnsInUnitigs = 0;
    NumRunsOfGapsInUnitigReads = 0;
    NumGapsInUnitigs = 0;
    NumColumnsInContigs = 0;
    NumRunsOfGapsInContigReads = 0;
    NumGapsInContigs = 0;
    NumAAMismatches = 0;
    NumUnitigRetrySuccess = 0;

    argc = AS_configure(argc, argv);

    while ( !errflg && 
           ( (ch = getopt(argc, argv, 
                 "a:d:e:fghil:mno:p:q:r:s:t:v:w:D:GIKM:NO:PR:S:T:UV:X")) != EOF))
    {
        switch(ch) {
        case 'n':
          noop = 1;
          iflags++;
          break;
        case 'f':
          allow_forced_frags = 1;
          iflags++;
          break;
        case 'g':
          allow_neg_hang = 1;
          iflags++;
          break;
        case 'G':
          allow_neg_hang_retry = 1;
          iflags++;
          break;
        case 'P':
          fprintf(stderr, "-P is depricated; protoIO is default.\n");
          iflags++;
          break;
        case 'K':
          options.split_alleles = 0;                   
          iflags++;
          break;
        case 'v':
          switch( atoi(optarg) ) {
          case 0:
            fprintf(stderr,"Command line switch %c turned off\n",ch); 
            break;
          case 1:
            printwhat = CNS_CONSENSUS;
            break;
          case 2:
            printwhat = CNS_DOTS;
            break;
          case 3:
            printwhat = CNS_NODOTS;
            break;
          case 4:
            printwhat = CNS_VIEW_UNITIG;
            break;
          case 5:
            printwhat = CNS_VERBOSE;
            break;
          default:
            fprintf(stderr,"Command line switch %c %d not supported; ignoring...\n",
                    ch,atoi(optarg)); 
          }
          iflags++;
          iflags++;
          break;
        case 'o':
          output_override=1;
          iflags++;
          if ( optarg[0] == '-' ) {
            std_output=1;
          } else {
            strcpy(OutputNameBuffer, optarg);
          }
          iflags++;
          break;
        case 'r':
          range=1;
          sscanf(optarg,"%d,%d", &tig_range.bgn, &tig_range.end);
          iflags++;
          iflags++;
          break;
        case 'S':
          partitioned = 1;
          partition = atoi(optarg);
          iflags++;
          iflags++;
          break;
        case 'w':
          options.smooth_win = atoi(optarg);
          iflags++;
          iflags++;
          break;
        case 'M':
          options.max_num_alleles = atoi(optarg);
          iflags++;
          iflags++;
          break;
        case 'm':
          in_memory = 1;
          iflags++;
          break;
        case 'I':
          align_ium = 1;
          cgbout = 1;
          iflags++;
          break;
        case 'U':
          no_contigs = 1;
          align_ium = 1;
          cgbout = 1;
          iflags++;
          break;
        case 's':
          USE_SDB = 1;
          strcpy(SeqStoreFileName, optarg);
          iflags++;
          iflags++;
          break;
        case 'p':
          USE_SDB_PART = 1;
          sdb_partition = atoi(optarg);
          iflags++;
          iflags++;
          break;
        case 'V':
          sdb_version = atoi(optarg);
          iflags++;
          iflags++;
          break;
        case 'X':
          expert = 1;
          iflags++;
          break;

        case 'D':
          if        (strcmp(optarg, "dumpunitigs") == 0) {
            DUMP_UNITIGS_IN_MULTIALIGNCONTIG = 1;
          } else if (strcmp(optarg, "verbosemultialign") == 0) {
            VERBOSE_MULTIALIGN_OUTPUT = 1;
          } else if (strcmp(optarg, "forceunitigabut") == 0) {
            FORCE_UNITIG_ABUT = 1;
          } else {
            fprintf(stderr, "Unrecognized option '%s' to -D.\n", optarg);
          }
          iflags++;
          iflags++;
          break;

        case 'R':  // restart (formerly 'continue')
          if ( ! expert ) {
             fprintf(stderr,"Command line switch %c requires -X; try adding -X...\n",
                  ch); 
             illegal_use = 1;
          } else {
            continue_at = atoi(optarg);
            beyond=0;
          }
          iflags++;
          iflags++;
          break;
        case 'i':
          if ( ! expert ) {
             fprintf(stderr,"Command line switch %c requires -X; try adding -X...\n",
                  ch); 
             illegal_use = 1;
          } else {
            align_ium = 1;
           }
          iflags++;
          break;
        case 'd':
          {
            CNS_USE_PUBLIC = atoi(optarg);
          }
          iflags++;
          iflags++;
          break;
        case 'q':
          if ( ! expert ) {
             fprintf(stderr,"Command line switch %c requires -X; try adding -X...\n",
                  ch); 
             illegal_use = 1;
          } else {
            sscanf(optarg,"%f:%d:%f",&CNS_SEQUENCING_ERROR_EST,&CNS_HAPLOTYPES,
                &CNS_SNP_RATE);
            if (!(CNS_SEQUENCING_ERROR_EST > 0) || 
                CNS_SEQUENCING_ERROR_EST > .10 ) 
            {
              fprintf(stderr,"ERROR: Sequencing error estimate (-q flag) should be "
                  "within (0,.10) (%4f was specified\n",
                  CNS_SEQUENCING_ERROR_EST);
              illegal_use = 1;
            }
            if (CNS_HAPLOTYPES < 1) {
              fprintf(stderr,"ERROR: Haplotypes sampled (-h flag) must be > 0 "
                             "(%d was specified\n",CNS_HAPLOTYPES);
              illegal_use = 1;
            }
            if ((CNS_SNP_RATE < 0) || CNS_SNP_RATE > .10 ) {
              fprintf(stderr,
                  "ERROR: SNP rate estimate (-s flag) should be within [0,.10) "
                  "(%4f was specified\n",CNS_SNP_RATE);
              illegal_use = 1;
            }
          }
          iflags++;
          iflags++;
          break;
        case 'e':
          if ( ! expert ) {
             fprintf(stderr,
                 "Command line switch %c requires -X; try adding -X...\n",
                  ch); 
              illegal_use = 1;
          } else {
            if ( optarg[0] == '#' ) {
              extract = atoi(&optarg[1]);
            } else {
              process_sublist = 1;
              extract = -2; // special value to indicate that extrating a sublist, rather than ind.
              sublist_file = optarg;
            }
          }
          iflags++;
          iflags++;
          break;
        case 'l':
          log_override=1;
          if ( optarg[0] == '-' ) {
            std_error_log=1;
          } else {
            strcpy(LogNameBuffer, optarg);
          }
          iflags++;
          iflags++;
          break;
        case 't':
          num_of_threads = atoi(optarg);
          iflags++;
          iflags++;
          break;
        case 'a':
          switch(optarg[0]) {
          case 'D':
            COMPARE_FUNC=DP_Compare;
            break;
          case 'L':
            COMPARE_FUNC=Local_Overlap_AS_forCNS;
            break;
          case 'A':
            COMPARE_FUNC=Affine_Overlap_AS_forCNS;
            break;
          default:
            fprintf(stderr,"Unrecognized value for option -%c (%s)",optopt,optarg);
            help_flag=1;
          } 
          iflags++;
          iflags++;
          break;
        case 'T':
          time_limit = atoi(optarg);
          iflags++;
          iflags++;
          fprintf(stderr,"The time limit is = " F_TIME_T " seconds.\n",
                  time_limit);
          // A time limit for creating profiling runs.
          break;
        case '?':
        case 'h':
          help_flag = 1;
          break;
        default :
          {
          help_flag = 1;
          fprintf(stderr,"Unrecognized option -%c",optopt);
          }
        }
    }
    if ( (argc - iflags) == 2) { std_input = 1; iflags--; }
    if ( (argc - iflags) != 3 ) help_flag = 1;
    if (help_flag) 
        help_message(argc, argv);
   
    if ( illegal_use ) {
        fprintf(stderr,"\n consensus -h provides usage information.\n");
        exit(1);
    }

    /****************          Open Gatekeeper Store             ***********/

    gkpStore = openGateKeeperStore(argv[optind++], FALSE);

    if (partitioned)
      loadGateKeeperPartition(gkpStore, partition);
    else if (in_memory)
      loadGateKeeperStorePartial(gkpStore, 0, 0, FRAG_S_QLT);

    /****************      Initialize reusable stores          ***********/
    sequenceStore = NULL;
    qualityStore = NULL;
    beadStore = NULL;
    columnStore = NULL;
    manodeStore = NULL;

    {
      /**************** Determine which messages to process and  ***********/ 
      /****************          what files to create            ***********/
      char *suffix=NULL;
      char *extract_id=NULL;
      int eid_len;
       
      if ( ! std_input ) 
      {
        suffix = strrchr(argv[optind],(int)'.');
        fprintf(stderr,"Input file is %s suffix is %s\n",argv[optind], suffix);
        strcpy(InputFileName,argv[optind]);
        cgwin = fopen(InputFileName,"r");
        if (cgwin == NULL ) {
          fprintf(stderr,"Could not open %s for CGW input.\n",InputFileName);
          CleanExit("",__LINE__,1);
        }
        if(suffix) *suffix = '\0';
      } else {
        cgwin = stdin;
      }

      switch(extract) 
      {
        case -2:
          fprintf(stderr,"Extracting ");
          eid_len = fprintf(stderr,"_%s_sublist",(align_ium)?"IUM":"ICM");
          extract_id  = (char *) safe_malloc(eid_len+1);
          sprintf(extract_id,"_%s_sublist",(align_ium)?"IUM":"ICM");
          fprintf(stderr,"\n");
          break;
        case -1:
          break;
          default:
          fprintf(stderr,"Extracting ");
          eid_len = fprintf(stderr,"_%s_%d",(align_ium)?"IUM":"ICM",extract);
          extract_id  = (char *) safe_malloc(eid_len+1);
          sprintf(extract_id,"_%s_%d",(align_ium)?"IUM":"ICM",extract);
          fprintf(stderr,"\n");
      }


      if ( !std_output) 
      {
        if (output_override) 
        {
           HandleDir(OutputNameBuffer,OutputFileName);
        } else if (cgbout == 1) {
          sprintf(OutputFileName,"%s%s.cgi",argv[optind],
                (extract != -1)?extract_id:"");
        } 
        else 
        {
          sprintf(OutputFileName,"%s%s.cns",argv[optind],
                (extract != -1)?extract_id:"");
        }
        sprintf(OutputFileNameTmp,"%s_tmp",OutputFileName);
        fprintf(stderr,"Output temporary file name is %s \n",OutputFileNameTmp);
        fprintf(stderr,"Output final file name is %s \n",OutputFileName);
        if (continue_at > 0) 
        {
          cnsout = fopen(OutputFileNameTmp, "a");     // append to existing cns file
        } else {
          cnsout = fopen(OutputFileNameTmp, "w");     // write new cns file
        }
        if (cnsout == NULL ) {
          fprintf(stderr,"Failure to create output temporary file %s\n", 
              OutputFileNameTmp);
          CleanExit("",__LINE__,1);
        }
      }   
      else 
      {
        cnsout = stdout;
      }


      if ( ! std_input ) {
         sprintf(LogFileName,"%s%s.clg",OutputFileNameTmp,
              (extract != -1)?extract_id:"");
      } else {
         sprintf(LogFileName,"cns_%d_%s.clg",getpid(),
              (extract != -1)?extract_id:"");
      }
      if ( ! std_error_log ) 
      {
        if (log_override ) {
          HandleDir(LogNameBuffer,LogFileName);
        }
        if (extract == -1 ) 
        {
          if ( ! std_input ) {
            sprintf(LogFileName,"%s%s.clg",OutputFileName,
              (extract != -1)?extract_id:"");
          } else {
            sprintf(LogFileName,"cns_%d_%s.clg",getpid(),
              (extract != -1)?extract_id:"");
          }
          fprintf(stderr,"Creating log file %s\n", LogFileName);
       
          if (continue_at > 0) {
            cnslog = fopen(LogFileName,"a");             // append to existing log file for cns
            fprintf(stderr,"Opened logfile %s\n", LogFileName);
          } else {
            cnslog = fopen(LogFileName,"w");             // start new log file for cns
            fprintf(stderr,"Opened logfile %s\n", LogFileName);
          }
          if (cnslog == NULL ) {
            fprintf(stderr,"Failure to create log file %s\n", LogFileName);
            CleanExit("",__LINE__,1);
          }
        }
      } 
      else 
      {
        cnslog = stderr;
      }
      if ( cnslog == NULL ) { 
        cnslog = fopen(LogFileName,"w");             // start new log file for cns
        fprintf(stderr,"Opened logfile %s\n", LogFileName);
      }
      if ( cnslog == NULL ) { 
        cnslog = stderr;   // write log to stderr
      }

      if (process_sublist) 
      {
        char   string[1000];
        int    num_uids;
        sublist = fopen(sublist_file,"r");
        if( sublist == NULL )
        {
          fprintf( stderr, "Failed to open list file %s for reading.\n", 
            sublist_file );
          CleanExit("",__LINE__,1);
        }
        num_uids = 0;
        while( fgets( string, 1000, sublist ) )
        {
          num_uids++;
        }
        rewind( sublist );
        tig_iids       = CreateScalarHashTable_AS( num_uids );
        tig_iids_found = CreateScalarHashTable_AS( num_uids );
        if( tig_iids == NULL || tig_iids_found == NULL ) {
            return EXIT_FAILURE;
        }
        for( i = 0; i < num_uids; i++ )
        {
          fgets( string, 1000, sublist );
          InsertInHashTable_AS(tig_iids, STR_TO_UID(string, NULL, 10), 0, 0, 0);
        }
        fclose( sublist );
      }
      safe_free(extract_id);
    }


    /**************** Prepare Unitig Store ****************************/
    if ( USE_SDB ) {
      if ( USE_SDB_PART ) {
        sequenceDB_part = openSequenceDBPartition(SeqStoreFileName, sdb_version, 
            sdb_partition);
      } else {
        sequenceDB = OpenSequenceDB(SeqStoreFileName, FALSE, sdb_version);
      }
    } 
    else 
    {
      unitigStore = CreateMultiAlignStoreT(0);
    }
#if 0 
    /* this may be introduced to save i/o time at contigging stage */
    sprintf(MAStoreFileName,"%s.uma",argv[optind]);
    if (cgbout == 1) {
      unitigStore = CreateMultiAlignStoreT();
    }  else {
      umain = fopen(MAStoreFileName,"r");
      unitigStore = LoadMultiAlignStoreTFromStream(umain);
    }
#endif

    /**************** Loop on Input Messages **************************/
    tp1 = time(NULL);
    {
      int contig_count=0,unitig_count=0;
      GenericMesg *pmesg;  
      GenericMesg tmesg;  
      IntConConMesg *pcontig;
      IntUnitigMesg *iunitig;
      MultiAlignT *ma;
      VA_TYPE(int32) *deltas=CreateVA_int32(1);
      VA_TYPE(char) *sequence=CreateVA_char(200000);
      VA_TYPE(char) *quality=CreateVA_char(200000);
      time_t t;
      t = time(0);
      fprintf(stderr,"# Consensus $Revision: 1.53 $ processing. Started %s\n",
        ctime(&t));
      InitializeAlphTable();
      if ( ! align_ium && USE_SDB && extract > -1 ) 
      {
        clear_range_to_use = AS_READ_CLEAR_LATEST;
        IntConConMesg ctmp;
        if ( USE_SDB_PART ) {
          ma = loadFromSequenceDBPartition(sequenceDB_part, extract);
        } else {
          ma =  LoadMultiAlignTFromSequenceDB(sequenceDB, extract, FALSE);
        }
        ctmp.iaccession = extract;
        ctmp.forced      = 0;
        ctmp.length      = GetNumchars(ma->consensus);
        ctmp.num_pieces  = GetNumIntMultiPoss(ma->f_list);
        ctmp.pieces      = GetIntMultiPos(ma->f_list,0);
        ctmp.num_unitigs = GetNumIntUnitigPoss(ma->u_list);
        ctmp.unitigs     = GetIntUnitigPos(ma->u_list,0);
        ctmp.placed      = AS_PLACED;
#if 0
        fprintf(stderr, "Before MultiAlignContig #%d: ctmp.num_vars = %d\n", ctmp.iaccession, ctmp.num_vars);
#endif
        ctmp.num_vars    = GetNumIntMultiVars(ma->v_list);
        if (ctmp.num_vars == 0)
        {
            ctmp.num_vars = 1;
            ctmp.v_list = safe_malloc(sizeof(IntMultiVar));
        }
        if (MultiAlignContig(&ctmp, sequence, quality, deltas, printwhat,
            COMPARE_FUNC, &options) != EXIT_SUCCESS)
        {
            num_contig_failures++;
            if (num_contig_failures <= MAX_NUM_CONTIG_FAILURES)
            {
              fprintf(stderr,"MultiAlignContig failed for contig %d\n", ctmp.iaccession);
              writeFailure(OutputFileName, std_output, pmesg);
            }
            else
            {
                CleanExit("MultiAlignContig failed  more than MAX_NUM_CONTIG_FAILURES times.Exit."
                    ,__LINE__,1);
            }
        } 
#if 0
        fprintf(stderr, "After  MultiAlignContig #%d: ctmp.num_vars = %d\n", ctmp.iaccession, ctmp.num_vars);
#endif
	tmesg.t = MESG_ICM; 
	tmesg.m = &ctmp; 
	WriteProtoMesg_AS(cnsout,&tmesg); 
	fflush(cnsout);

        if ( printwhat != CNS_STATS_ONLY && cnslog != NULL )
        {
           MultiAlignT *ma1 = CreateMultiAlignTFromICM(&ctmp,-1,0);
           PrintMultiAlignT(cnslog,ma1,gkpStore, 
                            1, 0, clear_range_to_use);
           fflush(cnslog);
           DeleteVA_char(ma1->consensus);
           DeleteVA_char(ma1->quality);
           DeleteVA_int32(ma1->delta);
           DeleteVA_IntMultiPos(ma1->f_list);
           DeleteVA_int32(ma1->udelta);
           DeleteVA_IntUnitigPos(ma1->u_list);
           DeleteVA_IntMultiVar(ma1->v_list);
           ctmp.num_vars = 0;
        }

        OutputScores(NumColumnsInUnitigs, NumRunsOfGapsInUnitigReads, 
                     NumGapsInUnitigs, NumColumnsInContigs, 
                     NumRunsOfGapsInContigReads, NumGapsInContigs,
                     NumAAMismatches, 
                     NumVARRecords, NumVARStringsWithFlankingGaps,
                     NumUnitigRetrySuccess);

        if (num_contig_failures > 0) {
          fprintf(stderr, "%d contigs failed.\n", num_contig_failures);
          exit(1);
        }

        exit(0); 
      }  //  end of "if ( ! align_ium && USE_SDB && extract > -1 )"

      while ( (ReadProtoMesg_AS(cgwin,&pmesg) != EOF)
            ) 
      { 
        switch(pmesg->t)
        {
          case MESG_IUM:
          {
            iunitig = (IntUnitigMesg *)(pmesg->m);
            clear_range_to_use = AS_READ_CLEAR_OBT;
            if (align_ium) 
            {
       
              // Form alignment of IUM
              if (extract > -1 && iunitig->iaccession != extract) break;
              if (!beyond && continue_at > 0 && iunitig->iaccession < continue_at) 
              {
                if ( iunitig->iaccession == continue_at -1) beyond=1;
                break;
              } 
              else if (range && iunitig->iaccession < tig_range.bgn) 
              {
                break;
              } 
              else {
                beyond=1;
              }
              if (process_sublist && ExistsInHashTable_AS(tig_iids, iunitig->iaccession, 0))
              {
                fprintf(stderr,"Processing IUM %d from sublist\n",iunitig->iaccession);
                InsertInHashTable_AS(tig_iids_found, iunitig->iaccession, 0, 0, 0);
              } 
              else if (process_sublist) 
              {
                // pass through the Unitig message only if extract == -1
                if (extract == -1) WriteProtoMesg_AS(cnsout,pmesg); 
                break;
              }
              if (extract != -1 ) {
                int eid_len;
                char *extract_id;
                fprintf(stderr,"Extracting ");
                eid_len = fprintf(stderr,"_IUM_%d",iunitig->iaccession);
                extract_id  = (char *) safe_malloc(eid_len+1);
                sprintf(extract_id,"_%s_%d",(align_ium)?"IUM":"ICM",
                    iunitig->iaccession);
                fprintf(stderr,"\n");
                sprintf(LogFileName,"%s%s.clg",OutputFileName,
                    (extract != -1)?extract_id:"");
                if (continue_at > 0) {
                   if ( cnslog == NULL ) {
                      cnslog = fopen(LogFileName,"a");   // append to existing log file for cns
                      fprintf(stderr,"Opened logfile %s\n", LogFileName);
                   }
                } else {
                   if ( cnslog == NULL ){
                      cnslog = fopen(LogFileName,"w");    // start new log file for cns
                      fprintf(stderr,"Opened logfile %s\n", LogFileName);
                   }
                }
                sprintf(CamFileName,"%s%s.cns.cam",argv[optind],
                  (extract != -1)?extract_id:"");
                cam = fopen(CamFileName,"w");         // cam file
                safe_free(extract_id);
              }
              if (range && 
                  (iunitig->iaccession < tig_range.bgn || 
                   iunitig->iaccession > tig_range.end )) 
              {
                if (iunitig->iaccession > tig_range.end ) 
                    exit(0); 
                break;
              }

              {
                int   unitigfail = 0;

                unitigfail = MultiAlignUnitig(iunitig, gkpStore, sequence,
                                              quality, deltas, printwhat, do_rez, COMPARE_FUNC, &options);

                if ((unitigfail == EXIT_FAILURE) &&
                    (allow_neg_hang_retry) &&
                    (allow_neg_hang == 0))
                  {
                    allow_neg_hang = 1;
                    unitigfail = MultiAlignUnitig(iunitig, gkpStore, sequence,
                                                  quality, deltas, printwhat, do_rez, COMPARE_FUNC, &options);
                    allow_neg_hang = 0;
                    if (unitigfail != EXIT_FAILURE)
                      NumUnitigRetrySuccess++;
                  }

                if (unitigfail == EXIT_FAILURE) {
                  num_unitig_failures++;
                  if (num_unitig_failures <= MAX_NUM_UNITIG_FAILURES)
                    { 
                      fprintf(stderr,"MultiAlignUnitig failed for unitig %d\n", iunitig->iaccession);
                      writeFailure(OutputFileName, std_output, pmesg);
                    }
                  else
                    {
                      CleanExit("MultiAlignUnitig failed  more than MAX_NUM_UNITIG_FAILURES times.Exit."
                                ,__LINE__,1);
                    } 
                }
              }

              // Create a MultiAlignT from the MANode
            }

            if ( ! no_contigs & ! USE_SDB ) {
                ma = CreateMultiAlignTFromIUM(iunitig,-1,0);
                SetMultiAlignInStore(unitigStore,ma->id,ma);
            }
            input_lengths+=GetUngappedSequenceLength(iunitig->consensus);

            if (range == 1 &&
               ((align_ium &&
                 iunitig->iaccession >= tig_range.bgn &&
                 iunitig->iaccession <= tig_range.end) ||
                (!align_ium && tig_range.bgn == 0)) ) 
            {
              WriteProtoMesg_AS(cnsout,pmesg); // pass through the Unitig message
            } 
            else if (extract == -1 && beyond) 
            {
              if (CNS_HAPLOTYPES == 1) {
                 // RemoveSNPs(iunitig);
              }
              WriteProtoMesg_AS(cnsout,pmesg); // pass through the Unitig message
            } 
            else if (align_ium && (process_sublist || 
                       iunitig->iaccession == extract)) 
            {
        //camview(cam,iunitig->iaccession,iunitig->f_list,iunitig->num_frags,NULL,0,gkpStore);
        //fclose(cam);
        //fclose(cnslog);
              if (CNS_HAPLOTYPES == 1) {
           // RemoveSNPs(iunitig);
              }
              WriteProtoMesg_AS(cnsout,pmesg); // pass through the Unitig message and continue
              if (iunitig->iaccession == extract) 
                 exit(0);
            }
            unitig_count++;

            break;
          }
          case MESG_ICM:
          {
            pcontig = (IntConConMesg *)(pmesg->m);
            clear_range_to_use = AS_READ_CLEAR_LATEST;
            if (extract > -1 && pcontig->iaccession != extract) break;
            if (extract != -1 && align_ium) break;
            if (!beyond && continue_at > 0 && pcontig->iaccession < continue_at ) 
            {
              break;
            } else { 
              beyond=1;
            }
            if (process_sublist && ExistsInHashTable_AS(tig_iids, iunitig->iaccession, 0))
            {
              fprintf(stderr,"Processing ICM %d from sublist\n", pcontig->iaccession);
              InsertInHashTable_AS(tig_iids_found, pcontig->iaccession, 0, 0, 0);
            } else if ( range == 1 && 
                       ( (!align_ium) && 
                        ( (pcontig->iaccession<tig_range.bgn) || 
                          (pcontig->iaccession>tig_range.end))))
            {
              if ( pcontig->iaccession > tig_range.end ) 
                  exit(0);
              break;
            } 
            else if (process_sublist)
            {
              // pass through the Contig message and continue
              if (extract == -1) WriteProtoMesg_AS(cnsout,pmesg); 
              break;
            }
        //qsort(pcontig->unitigs, pcontig->num_unitigs, sizeof(IntUnitigPos),
        //      (int (*)(const void *,const void *))IntUnitigPositionCmpLeft);
            if ( ! noop > 0 ) {
             
//             if (pcontig->num_vars == 0)
//              {
//                  pcontig->num_vars = 1;
//                  pcontig->v_list = safe_malloc(sizeof(IntMultiVar));
//              }
                pcontig->num_vars == 0;
                pcontig->v_list == NULL;
                if (MultiAlignContig(pcontig, sequence, quality, deltas, printwhat,
                    COMPARE_FUNC, &options) != EXIT_SUCCESS)
                {
                    num_contig_failures++;
                    if (num_contig_failures <= MAX_NUM_CONTIG_FAILURES)
                    {
                      fprintf(stderr,"MultiAlignContig failed for contig %d\n", pcontig->iaccession);
                      writeFailure(OutputFileName, std_output, pmesg);
                      break;
                    }
                    else
                    {
                        CleanExit("MultiAlignContig failed  more than MAX_NUM_CONTIG_FAILURES times.Exit."
                            ,__LINE__,1);
                    }
                }
            }
            if ( printwhat == CNS_CONSENSUS && cnslog != NULL && 
                 pcontig->num_pieces > 0)
            {
                ma = CreateMultiAlignTFromICM(pcontig,-1,0);
                PrintMultiAlignT(cnslog,ma,gkpStore, 1, 0, clear_range_to_use);
            }
            output_lengths+=GetUngappedSequenceLength(pcontig->consensus);
            pmesg->t = MESG_ICM; 
            pmesg->m = pcontig; 
            if ( range == 1 && 
                ( (!align_ium) && 
                 ( (pcontig->iaccession>=tig_range.bgn) && 
                   (pcontig->iaccession<=tig_range.end))))
            {
              WriteProtoMesg_AS(cnsout,pmesg); // pass through the Contig message
            } else if (extract == -1) {
              WriteProtoMesg_AS(cnsout,pmesg);
            } else if ( pcontig->iaccession == extract) {
              //camview(cam,pcontig->iaccession,pcontig->pieces,pcontig->num_pieces,pcontig->unitigs,
              //        pcontig->num_unitigs,gkpStore);
              WriteProtoMesg_AS(cnsout,pmesg);
              OutputScores(NumColumnsInUnitigs, NumRunsOfGapsInUnitigReads,
                           NumGapsInUnitigs, NumColumnsInContigs,
                           NumRunsOfGapsInContigReads, NumGapsInContigs,
                           NumAAMismatches, 
                           NumVARRecords, NumVARStringsWithFlankingGaps,
                           NumUnitigRetrySuccess);
                  exit(0);
            }
            if (pcontig->v_list != NULL) 
            {
                int i;
                for (i=0; i<pcontig->num_vars; i++)
                {
                    safe_free(pcontig->v_list[i].nr_conf_alleles);
                    safe_free(pcontig->v_list[i].weights);
                    safe_free(pcontig->v_list[i].var_seq);
                    safe_free(pcontig->v_list[i].conf_read_iids);
                }
                safe_free(pcontig->v_list);
            }
            pcontig->num_vars = 0;
            fflush(cnsout);
            fflush(cnslog);
            contig_count++;
            break;
          } 

          case MESG_ADT:
          {
            if (beyond || ( range == 1 && tig_range.bgn==0) ) 
            {
              AuditMesg *adt_mesg;

              adt_mesg = (AuditMesg *)(pmesg->m);
              pmesg->t = MESG_ADT;

#if 0
            {
              AuditLine auditLine;
              AppendAuditLine_AS(adt_mesg, &auditLine, t,
                                 "Consensus", "$Revision: 1.53 $","(empty)");
            }
#endif
              VersionStampADT(adt_mesg,argc,argv);
              WriteProtoMesg_AS(cnsout,pmesg);
            }
          }
          break;
          default:
          { }
    /*    Passing through any other messages  */
          if ((range &&
               ((tig_range.bgn == 0 && unitig_count+contig_count == 0)
                || (align_ium && tig_range.end == unitig_count)
                || (!align_ium && tig_range.end == contig_count)))
              || (beyond && extract == -1)) {
            WriteProtoMesg_AS(cnsout,pmesg);
          }
        }
        fflush(cnsout);
        fflush(cnslog);
      }

      t = time(0);
      fprintf(stderr,"# Consensus $Revision: 1.53 $ Finished %s\n",ctime(&t));
      if (printcns) 
      {
        int unitig_length = (unitig_count>0)? (int) input_lengths/unitig_count: 0; 
        int contig_length = (contig_count>0)? (int) output_lengths/contig_count: 0;
         
        fprintf(stderr,"\nProcessed %d unitigs and %d contigs.\n",
          unitig_count,contig_count);
        fprintf(stderr,"\nAverage unitig length: %d  Effective coverage: %.2f\n",
           unitig_length, EffectiveCoverage(unitig_length));
        fprintf(stderr,"\nAverage contig length: %d  Effective coverage: %.2f\n",
          contig_length,EffectiveCoverage(contig_length));
        fprintf(cnslog,"\nProcessed %d unitigs and %d contigs.\n",
          unitig_count,contig_count);
          fprintf(cnslog,"\nAverage unitig length: %d  Effective coverage: %.2f\n",
         unitig_length, EffectiveCoverage(unitig_length));
        fprintf(cnslog,"\nAverage contig length: %d  Effective coverage: %.2f\n",
          contig_length,EffectiveCoverage(contig_length));
#ifdef CNS_TIME_OVERLAPS
        fprintf(stderr,"%d olaps computed in %7.1f sec\n",OverlapCount,
            (double) OlapTime/ CLOCKS_PER_SEC);
#endif
      }
    }


    if (cgbout == 1) {  /* This may be used later to bypass proto i/o for MultialignTs */
//    umaout = fopen(MAStoreFileName,"w");
//    SaveMultiAlignStoreTToStream(unitigStore,umaout,0);
//    fclose(umaout);
    }
    if ( unitigStore ) DeleteMultiAlignStoreT(unitigStore);
    {
      fclose(cnsout);
      if ( ! std_output ) {
        int ierr = rename( OutputFileNameTmp, OutputFileName );
        if(ierr != 0) {
          perror("ERROR: failure moving the cnsout file to final file name.\n");
          num_contig_failures++;
        }
      }
    }

    OutputScores(NumColumnsInUnitigs, NumRunsOfGapsInUnitigReads,
                 NumGapsInUnitigs, NumColumnsInContigs,
                 NumRunsOfGapsInContigReads, NumGapsInContigs,
                 NumAAMismatches, 
                 NumVARRecords, NumVARStringsWithFlankingGaps,
                 NumUnitigRetrySuccess);


    if (num_unitig_failures)
        fprintf(stderr, "\nTotal number of unitig failures= %d\n", num_unitig_failures);
    if (num_contig_failures)
        fprintf(stderr, "\nTotal number of contig failures= %d\n", num_contig_failures);
    if (num_unitig_failures || num_contig_failures)
        return EXIT_FAILURE;
    else
        return EXIT_SUCCESS;
}