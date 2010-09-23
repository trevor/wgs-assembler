
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

static char CM_ID[] = "$Id: Consensus_CNS.c,v 1.1.1.1 2004-04-14 13:51:12 catmandew Exp $";

// Operating System includes:
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

// Celera Assembler includes:
#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
#include "AS_PER_genericStore.h"
#include "AS_UTL_Var.h"
#include "UtilsREZ.h"
#include "AS_UTL_ID_store.h"
#include "PrimitiveVA.h"
#include "PrimitiveVA_MSG.h"
#include "AS_UTL_version.h"
#include "AS_SDB_SequenceDBPartition.h"
#include "AS_ALN_forcns.h"

// Consensus includes:
#include "MultiAlignStore_CNS.h"
#include "MultiAlignment_CNS.h"
#include "Globals_CNS.h"
#include "PublicAPI_CNS.h"

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
   mode_t mode = S_IRWXU | S_IRWXG | S_IROTH;
   char *suffix;
   char *DirName;
   char *FileName;
   DIR *Dir;
   suffix = strrchr(filePathAndName,(int)'/');
   if ( suffix != NULL ) {
      *suffix = '\0';
      DirName = filePathAndName; 
      if ( DirName != NULL ) {
        Dir = opendir(DirName);
        if ( Dir == NULL ) {
          if(mkdir(DirName,mode)){
            fprintf(stderr,"Failure to create directory %s\n", DirName);
            exit(1);
          }
        }
      }
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

int main (int argc, char *argv[]) {
  USE_SDB=0;
  USE_SDB_PART=0;

    MesgReader   reader;
    MesgWriter   writer;
    int binary_io;
    char InputFileName[FILENAME_MAX];
    char OutputNameBuffer[FILENAME_MAX];
#if 0
    char MAStoreFileName[FILENAME_MAX];
    FILE *umaout;
    FILE *umain;
#endif
    char LogNameBuffer[FILENAME_MAX];
    char CamFileName[FILENAME_MAX];
    char BactigStoreFileName[FILENAME_MAX];
    char SeqStoreFileName[FILENAME_MAX];
    char *sublist_file = NULL;
    FILE *cgwin;
    FILE *cam;
    FILE *sublist;
    int sdb_version=-1;
    int sdb_partition=-1;
    ID_Arrayp  tig_iids = NULL;
    ID_Arrayp  tig_iids_found = NULL;
    cds_int64  this_id;

    /**************** Process Command Line Arguments *********************/
    /* Parse the argument list using "man 3 getopt". */
    int align_ium=0;
    int no_contigs=0;
    int cgbout=0;
    int printcns=1;
    int extract=-1;
    int continue_at=0;
    int beyond=1;
    int bactigs=0;
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
    Overlap *(*COMPARE_FUNC)(COMPARE_ARGS)=Local_Overlap_AS_forCNS;
    SeqInterval tig_range;
    CNS_PrintKey printwhat=CNS_STATS_ONLY;
    int ch,errflg=0,illegal_use=0,help_flag=0,iflags=0;
    OutputType output=AS_BINARY_OUTPUT;
    int num_of_threads = 0;
    time_t time_limit = 0, tp1 = 0;

    optarg = NULL;
    debug_out = 1;
    terminate_cond = 1;
    partitioned=0;
    allow_forced_frags=0;
    allow_neg_hang=0;
    ALIGNMENT_CONTEXT=AS_CONSENSUS;
    
    while (!errflg && ((ch = getopt(argc, argv, "gnfhPv:d:O:o:r:S:mIUXAD:ECR:iGp:q:e:l:t:T:zs:V:a:")) != EOF)) {
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
        case 'P':
          output = AS_PROTO_OUTPUT;
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
        case 'O':
          bactigs = 1;
          terminate_cond = 0;
          debug_out = 1;
          strcpy(BactigStoreFileName, optarg);
          ALIGNMENT_CONTEXT = AS_OVERLAY;
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
        case 'z':
          //do_rez = 1;
          iflags++;
          break;
        case 'X':
          expert = 1;
          iflags++;
          break;
        case 'D':
          if ( ! expert ) {
             fprintf(stderr,"Command line switch %c requires -X; try adding -X...\n",
                  ch); 
             illegal_use = 1;
          } else {
            debug_out = atoi(optarg);
            iflags++;
            iflags++;
          }
          break;
        case 'A':
        case 'E':
        case 'C':
          if ( ! expert ) {
             fprintf(stderr,"Command line switch %c requires -X; try adding -X...\n",
                  ch); 
             illegal_use = 1;
          } else {
            fprintf(stderr,"Command line switch %c <%d> not supported; ignoring...\n",
                  ch, atoi(optarg)); 
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
        case 'G':
          if ( ! expert ) {
             fprintf(stderr,"Command line switch %c requires -X; try adding -X...\n",
                  ch); 
             illegal_use = 1;
          } else {
            CNS_CALL_PUBLIC = 1;
          }
          iflags++;
          break;
        case 'q':
          if ( ! expert ) {
             fprintf(stderr,"Command line switch %c requires -X; try adding -X...\n",
                  ch); 
             illegal_use = 1;
          } else {
            sscanf(optarg,"%f:%d:%f",&CNS_SEQUENCING_ERROR_EST,&CNS_HAPLOTYPES,&CNS_SNP_RATE);
            if (!(CNS_SEQUENCING_ERROR_EST > 0) || CNS_SEQUENCING_ERROR_EST > .10 ) {
              fprintf(stderr,"ERROR: Sequencing error estimate (-q flag) should be "
                            "within (0,.10) (%4f was specified\n",CNS_SEQUENCING_ERROR_EST);
              illegal_use = 1;
            }
            if (CNS_HAPLOTYPES < 1) {
              fprintf(stderr,"ERROR: Haplotypes sampled (-h flag) must be > 0 "
                             "(%d was specified\n",CNS_HAPLOTYPES);
              illegal_use = 1;
            }
            if ((CNS_SNP_RATE < 0) || CNS_SNP_RATE > .10 ) {
              fprintf(stderr,"ERROR: SNP rate estimate (-s flag) should be within [0,.10) "
                             "(%4f was specified\n",CNS_SNP_RATE);
              illegal_use = 1;
            }
          }
          iflags++;
          iflags++;
          break;
        case 'e':
          if ( ! expert ) {
             fprintf(stderr,"Command line switch %c requires -X; try adding -X...\n",
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
    if (help_flag) {
        fprintf(stderr,"  Usage:\n\n");
        fprintf(stderr,"  %s [-P] [-v level] [-O BactigStoreDir] [-I] [-a [DLA]] [-X expert_options] FragStoreDir [CGWStream]\n",argv[0]);
        fprintf(stderr,"\n Standard option flags:\n");
        fprintf(stderr,"    -P           Force ASCII .cns output \n");
        fprintf(stderr,"    -v [0-4]     Verbose:  0 = verbose off \n");
        fprintf(stderr,"                           1 = horizontal multi-alignment print in .clg\n");
        fprintf(stderr,"                           2 = 'dots'     multi-alignment print in .clg\n");
        fprintf(stderr,"                           3 = like 2, but dots are replaced with whitespace\n");
        fprintf(stderr,"                           4 = like 1, but with unitigs in  multi-alignment print in .clg\n");
        fprintf(stderr,"    -O BtigStore Overlay Assembler mode: BactigStoreDir argument required\n");
        fprintf(stderr,"    -T secs      time threshold which, if exceeded, should trigger clean exit.\n");
        fprintf(stderr,"    -S partition Use fragStorePartition partition\n");
        fprintf(stderr,"    -m           Load fragStorePartition into memory (default reads from disk)\n");
        fprintf(stderr,"    -U           Unitigs ONLY\n");
        fprintf(stderr,"    -I           IUM message alignment, write to .cgi file (instead of .cns)\n");
        fprintf(stderr,"                 (will also process contigs if they exist in the input file\n");
        fprintf(stderr,"                 (if -I is not specified, ICMs will be processed with the assumption\n");
        fprintf(stderr,"                  that IUMs were processed in a prior call.\n");
        fprintf(stderr,"    -r s,e       Process only message within range [s,e)\n");
        fprintf(stderr,"                 if s==0, header messages are  passed through\n");
        fprintf(stderr,"                 if e==last, trailing messages are  passed through\n");
        fprintf(stderr,"                 (facilitates multiprocessing in batch mode\n");
        fprintf(stderr,"    -a [DLA]     Specify aligner to use should DP_Compare fail\n");
        fprintf(stderr,"                 L = Local_Aligner (default)\n");
        fprintf(stderr,"                 D = standard DP_Compare (will cause failed overlaps to terminate the run)\n");
        fprintf(stderr,"                 A = Affine_Aligner\n");
        fprintf(stderr,"    -d int       Depth of Celera coverage below which to include external data in basecalling\n");
        fprintf(stderr,"                    0 (default) indicates that external data should always be used\n");
        fprintf(stderr,"                    1 yields the traditional behavior, which uses external only in absence of Celera\n");
        fprintf(stderr,"                  > 1 will include publice data is the Celera depth falls below the given value\n");
        fprintf(stderr,"    -X           Allow 'expert' options (following)\n");
        fprintf(stderr,"\n Expert option flags:\n");
        fprintf(stderr,"    -R %%d        Restart from the given ICM/IUM by internal id, appending to output file\n");
        fprintf(stderr,"    -i           Realign IUM messages (while processing .cgw file)\n");
        fprintf(stderr,"    -q string    Override default quality call parameters\n");
        fprintf(stderr,"                    string is colon separated list of the form '%%f:%%d:%%f'\n");
        fprintf(stderr,"                    where first field is estimated sequencing error rate (default: .015)\n");
        fprintf(stderr,"                         second field is number of sequenced haplotypes (default: 1)\n");
        fprintf(stderr,"                          third field is estimated SNP rate (default: 1/1000)\n");
        fprintf(stderr,"    -e #%%d       Extract only a single ICM/IUM by internal id\n");
        fprintf(stderr,"    -e idfile    Extract list of ICM/IUMs by internal ids provided in idfile\n");
        fprintf(stderr,"\n Arguments:\n");
        fprintf(stderr,"   FragStoreDir      path to previously created Fragment Store\n");
        fprintf(stderr,"  [InputStream]      previously created .cgw/.cgb file (if not specified, stdin)\n\n");
        fprintf(stderr,"\n Output:\n");
        fprintf(stderr,"   Creates a .cns file by default (or appends to .cns if -R is specified)\n");
        fprintf(stderr,"   -I sends output to a .cgi (post-unitigging consensus) file instead.\n");
        fprintf(stderr,"   -o <filename>     Overrides default output filename\n");
        fprintf(stderr,"   -o -              Overrides default output filename, sending output to stdout\n");
        fprintf(stderr,"   -l <filename>     Overrides default log filename\n");
        fprintf(stderr,"   -l -              Overrides default log filename, sending output to stderr\n");
        exit(1);
    }
    if ( illegal_use ) {
        fprintf(stderr,"\n consensus -h provides usage information.\n");
        exit(1);
    }

    /****************          Open Fragment Store             ***********/
    if ( partitioned ) {
      global_fragStorePartition = openFragStorePartition(argv[optind++],partition,in_memory);
      global_fragStore = NULLFRAGSTOREHANDLE;
    } else {
      if ( in_memory ) {
        global_fragStore = loadFragStore(argv[optind++]);
      } else {
        global_fragStore = openFragStore(argv[optind++], "rb");
      }
      global_fragStorePartition = NULL;
    }
    if (bactigs) {
      global_bactigStore = openFragStore(BactigStoreFileName, "rb");
    }

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
       
    if ( ! std_input ) {
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
    binary_io = (InputFileType_AS( cgwin ) == ReadProtoMesg_AS) ? 0 : 1;
    reader = (binary_io == 1) ? ReadBinaryMesg_AS : ReadProtoMesg_AS;



    switch(extract) {
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
    if ( !std_output) {
      if (output_override) {
         HandleDir(OutputNameBuffer,OutputFileName);
      } else if (cgbout == 1) {
        sprintf(OutputFileName,"%s%s.cgi",argv[optind],(extract != -1)?extract_id:"");
      } else {
        sprintf(OutputFileName,"%s%s.cns",argv[optind],(extract != -1)?extract_id:"");
      }
      sprintf(OutputFileNameTmp,"%s_tmp",OutputFileName);
      fprintf(stderr,"Output temporary file name is %s \n",OutputFileNameTmp);
      fprintf(stderr,"Output final file name is %s \n",OutputFileName);
      if (continue_at > 0) {
        cnsout = fopen(OutputFileNameTmp, "a");     // append to existing cns file
      } else {
        cnsout = fopen(OutputFileNameTmp, "w");     // write new cns file
      }
      if (cnsout == NULL ) {
        fprintf(stderr,"Failure to create output temporary file %s\n", OutputFileNameTmp);
        CleanExit("",__LINE__,1);
      }
    } else {
      cnsout = stdout;
    }
    if ( ! std_input ) {
       sprintf(LogFileName,"%s%s.clg",OutputFileNameTmp,(extract != -1)?extract_id:"");
    } else {
       sprintf(LogFileName,"cns_%d_%s.clg",getpid(),(extract != -1)?extract_id:"");
    }
    if ( ! std_error_log ) {
      if (log_override ) {
        HandleDir(LogNameBuffer,LogFileName);
      }
      if (extract == -1 ) {
        if ( ! std_input ) {
          sprintf(LogFileName,"%s%s.clg",OutputFileName,(extract != -1)?extract_id:"");
        } else {
          sprintf(LogFileName,"cns_%d_%s.clg",getpid(),(extract != -1)?extract_id:"");
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
    } else {
      cnslog = stderr;
    }
    if ( cnslog == NULL ) { 
      cnslog = fopen(LogFileName,"w");             // start new log file for cns
      fprintf(stderr,"Opened logfile %s\n", LogFileName);
    }
    if ( cnslog == NULL ) { 
      cnslog = stderr;   // write log to stderr
    }
    if ( output == AS_PROTO_OUTPUT ) {
      writer = WriteProtoMesg_AS;
    } else {
      writer = (binary_io == 1) ? WriteBinaryMesg_AS : WriteProtoMesg_AS;
    }
    if (process_sublist) {
      char   string[1000];
      int    num_uids;
      sublist = fopen(sublist_file,"r");
      if( sublist == NULL )
      {
        fprintf( stderr, "Failed to open list file %s for reading.\n", sublist_file );
        CleanExit("",__LINE__,1);
      }
      num_uids = 0;
      while( fgets( string, 1000, sublist ) )
      {
        num_uids++;
      }
      rewind( sublist );
      tig_iids = AllocateID_Array( num_uids );
      tig_iids_found = AllocateID_Array( num_uids );
      if( tig_iids == NULL || tig_iids_found == NULL ) return 1;
      for( this_id = 0; this_id < num_uids - 1; this_id++ )
      {
        fgets( string, 1000, sublist );
        AppendToID_Array( tig_iids, STR_TO_UID(string, NULL, 10), 0 );
      }
      fgets( string, 1000, sublist );
      AppendToID_Array( tig_iids, STR_TO_UID( string, NULL, 10), 1 );
  
      fclose( sublist );
    }
    if (extract_id) free(extract_id);
    }

    /**************** Prepare Unitig Store ****************************/
    if ( USE_SDB ) {
      if ( USE_SDB_PART ) {
        sequenceDB_part = openSequenceDBPartition(SeqStoreFileName, sdb_version, sdb_partition);
      } else {
        sequenceDB = OpenSequenceDB(SeqStoreFileName, FALSE, sdb_version);
      }
    } else {
      unitigStore = CreateMultiAlignStoreT(0);
      if (bactigs ) {
        bactig_delta_length = CreateVA_int32(getLastElemFragStore(global_bactigStore));
        bactig_deltas = CreateVA_PtrT(getLastElemFragStore(global_bactigStore));
      }
    }
#if 0 
    /* this may be introduced to save i/o time at contigging stage */
    sprintf(MAStoreFileName,"%s.uma",argv[optind]);
    if (cgbout == 1 || bactigs) {
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
    fprintf(stderr,"# Consensus $Revision: 1.1.1.1 $ processing. Started %s\n",ctime(&t));
    InitializeAlphTable();
    if ( ! align_ium && USE_SDB && extract > -1 ) {
       IntConConMesg ctmp;
       if ( USE_SDB_PART ) {
         ma = loadFromSequenceDBPartition(sequenceDB_part, extract);
       } else {
         ma =  LoadMultiAlignTFromSequenceDB(sequenceDB, extract, FALSE);
       }
       ctmp.iaccession = extract;
       ctmp.forced =0;
       ctmp.length = GetNumchars(ma->consensus);
       ctmp.num_pieces = GetNumIntMultiPoss(ma->f_list);
       ctmp.pieces = GetIntMultiPos(ma->f_list,0);
       ctmp.num_unitigs = GetNumIntUnitigPoss(ma->u_list);
       ctmp.unitigs= GetIntUnitigPos(ma->u_list,0);
       MultiAlignContig(&ctmp, sequence, quality, deltas, printwhat,COMPARE_FUNC);
       if ( printwhat != CNS_STATS_ONLY && cnslog != NULL ){
          ma = CreateMultiAlignTFromICM(&ctmp,-1,0);
          PrintMultiAlignT(cnslog,ma,global_fragStore,global_fragStorePartition, global_bactigStore,
                           1,0,READSTRUCT_LATEST);
         tmesg.t = MESG_ICM; 
         tmesg.m = &ctmp; 
         writer(cnslog,&tmesg); 
         fflush(cnslog);
       }
       exit(0); 
    }
    while ( (reader(cgwin,&pmesg) != EOF)
#if 1
            &&((time_limit == 0) || ((time(NULL)-tp1) <= time_limit))
#endif            
            ) { 
    switch(pmesg->t){
    case MESG_IUM:
      {
      iunitig = (IntUnitigMesg *)(pmesg->m);
      if (align_ium) {
        // Form alignment of IUM
        if (extract > -1 && iunitig->iaccession != extract) break;
        if (!beyond && continue_at > 0 && iunitig->iaccession < continue_at) {
          if ( iunitig->iaccession == continue_at -1) beyond=1;
          break;
        } else if (range && iunitig->iaccession < tig_range.bgn) {
          break;
        } else {
          beyond=1;
        }
        if( process_sublist && (this_id = FindID_ArrayID( tig_iids, iunitig->iaccession)) > -1 )
        {
          fprintf(stderr,"Processing IUM %d from sublist\n",iunitig->iaccession);
          AppendToID_Array( tig_iids_found, iunitig->iaccession, 1 );
        } else if (process_sublist) {
          // pass through the Unitig message only if extract == -1
          if (extract == -1) writer(cnsout,pmesg); 
          break;
        }
        if (extract != -1 ) {
          int eid_len;
          char *extract_id;
          fprintf(stderr,"Extracting ");
          eid_len = fprintf(stderr,"_IUM_%d",iunitig->iaccession);
          extract_id  = (char *) safe_malloc(eid_len+1);
          sprintf(extract_id,"_%s_%d",(align_ium)?"IUM":"ICM",iunitig->iaccession);
          fprintf(stderr,"\n");
          sprintf(LogFileName,"%s%s.clg",OutputFileName,(extract != -1)?extract_id:"");
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
          sprintf(CamFileName,"%s%s.cns.cam",argv[optind],(extract != -1)?extract_id:"");
          cam = fopen(CamFileName,"w");         // cam file
          free(extract_id);
        }
        if (range && (iunitig->iaccession < tig_range.bgn || iunitig->iaccession > tig_range.end )) {
          if (iunitig->iaccession > tig_range.end ) exit(0); 
           break;
        }
        if (MultiAlignUnitig(iunitig,global_fragStore,sequence,quality,deltas,printwhat,do_rez,COMPARE_FUNC)==-1 ) {
            fprintf(stderr,"MultiAlignUnitig failed for unitig %d\n",iunitig->iaccession);
            //break;  // un-comment this line to allow failures... intended for diagnosis only.
            assert(FALSE);
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
       (!align_ium && tig_range.bgn == 0)) ) {
        writer(cnsout,pmesg); // pass through the Unitig message
      } else if (extract == -1 && beyond) {
        if (CNS_HAPLOTYPES == 1) {
           // RemoveSNPs(iunitig);
        }
        writer(cnsout,pmesg); // pass through the Unitig message
      } else if (align_ium && (process_sublist || iunitig->iaccession == extract)) {
        //camview(cam,iunitig->iaccession,iunitig->f_list,iunitig->num_frags,NULL,0,global_fragStore);
        //fclose(cam);
        //fclose(cnslog);
        if (CNS_HAPLOTYPES == 1) {
           // RemoveSNPs(iunitig);
        }
        writer(cnsout,pmesg); // pass through the Unitig message and continue
        if (iunitig->iaccession == extract) exit(0);
      }
      unitig_count++;
      break;
      }
    case MESG_ICM:
      {
      pcontig = (IntConConMesg *)(pmesg->m);
      if (extract > -1 && pcontig->iaccession != extract) break;
      if (extract != -1 && align_ium) break;
      if (!beyond && continue_at > 0 && pcontig->iaccession < continue_at ) {
        break;
      } else { 
        beyond=1;
      }
      if( process_sublist && (this_id = FindID_ArrayID( tig_iids, pcontig->iaccession)) > -1 )
      {
        fprintf(stderr,"Processing ICM %d from sublist\n",pcontig->iaccession);
        AppendToID_Array( tig_iids_found, pcontig->iaccession, 1 );
      } else if (range == 1 && ( (!align_ium) && ( (pcontig->iaccession<tig_range.bgn) || (pcontig->iaccession>tig_range.end)))){
        if ( pcontig->iaccession > tig_range.end ) exit(0);
        break;
      } else if (process_sublist){
        // pass through the Contig message and continue
        if (extract == -1) writer(cnsout,pmesg); 
        break;
      }
      //qsort(pcontig->unitigs, pcontig->num_unitigs, sizeof(IntUnitigPos),
      //      (int (*)(const void *,const void *))IntUnitigPositionCmpLeft);
      if ( ! noop > 0 ) {
        MultiAlignContig(pcontig, sequence, quality, deltas, printwhat,COMPARE_FUNC);
      }
      if ( printwhat == CNS_CONSENSUS && cnslog != NULL && pcontig->num_pieces > 0){
          ma = CreateMultiAlignTFromICM(pcontig,-1,0);
          PrintMultiAlignT(cnslog,ma,global_fragStore,global_fragStorePartition, global_bactigStore,
                           1,0,READSTRUCT_LATEST);
      }
      output_lengths+=GetUngappedSequenceLength(pcontig->consensus);
      pmesg->t = MESG_ICM; 
      pmesg->m = pcontig; 
      if (range == 1 && ( (!align_ium) && ( (pcontig->iaccession>=tig_range.bgn) && (pcontig->iaccession<=tig_range.end)))){
        writer(cnsout,pmesg); // pass through the Contig message
      } else if (extract == -1) {
        writer(cnsout,pmesg);
      } else if ( pcontig->iaccession == extract) {
        //camview(cam,pcontig->iaccession,pcontig->pieces,pcontig->num_pieces,pcontig->unitigs,
        //        pcontig->num_unitigs,global_fragStore);
        writer(cnsout,pmesg);
        exit(0);
      }
      fflush(cnsout);
      fflush(cnslog);
      contig_count++;
      break;
      } 
    case MESG_ADT:
      {
        if (beyond || ( range == 1 && tig_range.bgn==0) ) {
        AuditMesg *adt_mesg;

        adt_mesg = (AuditMesg *)(pmesg->m);
        pmesg->t = MESG_ADT;

#if 0
        {
          AuditLine auditLine;
          AppendAuditLine_AS(adt_mesg, &auditLine, t,
                             "Consensus", "$Revision: 1.1.1.1 $","(empty)");
        }
#endif
        VersionStampADT(adt_mesg,argc,argv);
        writer(cnsout,pmesg);
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
        writer(cnsout,pmesg);
      }
    }
    fflush(cnsout);
    fflush(cnslog);
  }
  t = time(0);
  fprintf(stderr,"# Consensus $Revision: 1.1.1.1 $ Finished %s\n",ctime(&t));
  if (printcns) {
    int unitig_length = (unitig_count>0)? (int) input_lengths/unitig_count: 0; 
    int contig_length = (contig_count>0)? (int) output_lengths/contig_count: 0;
         
    fprintf(stderr,"\nProcessed %d unitigs and %d contigs.\n",unitig_count,contig_count);
    fprintf(stderr,"\nAverage unitig length: %d  Effective coverage: %.2f\n",
         unitig_length, EffectiveCoverage(unitig_length));
    fprintf(stderr,"\nAverage contig length: %d  Effective coverage: %.2f\n",contig_length,EffectiveCoverage(contig_length));
    fprintf(cnslog,"\nProcessed %d unitigs and %d contigs.\n",unitig_count,contig_count);
    fprintf(cnslog,"\nAverage unitig length: %d  Effective coverage: %.2f\n",
         unitig_length, EffectiveCoverage(unitig_length));
    fprintf(cnslog,"\nAverage contig length: %d  Effective coverage: %.2f\n",contig_length,EffectiveCoverage(contig_length));
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
    char buffer[2*FILENAME_MAX+4];
    int ierr;
    fclose(cnsout);
    if ( ! std_output ) {
      sprintf(buffer,"mv %s %s",OutputFileNameTmp,OutputFileName);
      ierr = system(buffer);
      if(ierr) {
        fprintf(stderr,
                "ERROR: failure moving the cnsout file to final file name.\n");}
      }
   }
   return 0;
}