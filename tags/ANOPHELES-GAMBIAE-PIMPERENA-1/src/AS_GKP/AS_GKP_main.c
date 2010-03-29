
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

static char CM_ID[] = "$Id: AS_GKP_main.c,v 1.39 2007-05-16 11:48:37 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_PER_genericStore.h"
#include "AS_PER_gkpStore.h"
#include "AS_UTL_PHash.h"
#include "AS_UTL_version.h"
#include "AS_MSG_pmesg.h"
#include "AS_GKP_include.h"

GateKeeperStore *gkpStore;

static
void
usage(char *filename) {
  fprintf(stderr, "usage1: %s -o gkpStore [append/create options] <input.frg> <input.frg> ...\n", filename);
  fprintf(stderr, "usage2: %s -P partitionfile gkpStore\n", filename);
  fprintf(stderr, "usage3: %s [dump-options] gkpStore\n", filename);
  fprintf(stderr, "\n");
  fprintf(stderr, "The first usage will append to or create a GateKeeper store:\n");
  fprintf(stderr, "  -a                     append to existing tore\n");
  fprintf(stderr, "  -e <errorThreshhold>   set error threshhold\n");
  fprintf(stderr, "  -o <gkpStore>          append to or create gkpStore\n");
  fprintf(stderr, "  -G                     gatekeeper for assembler Grande (default)\n");
  fprintf(stderr, "  -T                     gatekeeper for assembler Grande with Overlap Based Trimming\n");
  fprintf(stderr, "  -Q                     don't check quality-based data quality\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  -v <vector-info>       load vector clear ranges into each read.\n");
  fprintf(stderr, "                         MUST be done on an existing, complete store.\n");
  fprintf(stderr, "                         example: -a -v vectorfile -o that.gkpStore\n");
  fprintf(stderr, "                         format: 'UID vec-clr-begin vec-clr-end'\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "The second usage will partition an existing store, allowing\n");
  fprintf(stderr, "the entire store partition to be loaded into memory.\n");
  fprintf(stderr, "  -P <partitionfile>     a list of (partition fragiid)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "The third usage will dump the contents of a GateKeeper store.\n");
  fprintf(stderr, "  [selection of what objects to dump]\n");
  fprintf(stderr, "  -b <begin-iid>         dump starting at this batch, library or read (1)\n");
  fprintf(stderr, "  -e <ending-iid>        dump stopping after this iid (1)\n");
  fprintf(stderr, "  -uid <uid-file>        dump only objects listed in 'uid-file' (1)\n");
  fprintf(stderr, "  -iid <iid-file>        dump only objects listed in 'iid-file' (1)\n");
  fprintf(stderr, "  -randommated <lib> <n> pick n mates (2n frags) at random from library lib,\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  [how to dump it]\n");
  fprintf(stderr, "  -tabular               dump info, batches, libraries or fragments in a tabular\n");
  fprintf(stderr, "                         format (for -dumpinfo, -dumpbatch, -dumplibraries,\n");
  fprintf(stderr, "                         and -dumpfragments, ignores -withsequence and -clear)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  [format of dump]\n");
  fprintf(stderr, "  -dumpinfo              print information on the store\n");
  fprintf(stderr, "    -lastfragiid         just print the last IID in the store\n");
  fprintf(stderr, "  -dumpbatch             dump all batch records\n");
  fprintf(stderr, "  -dumplibraries         dump all library records\n");
  fprintf(stderr, "  -dumpfragments         dump fragment info, no sequence\n");
  fprintf(stderr, "    -withsequence          ...and include sequence\n");
  fprintf(stderr, "    -clear <clr>           ...in clear range <clr>, default=UNTRIM\n");
  fprintf(stderr, "  -dumpfasta[seq|qlt]    dump fragment sequence or quality, as fasta format\n");
  fprintf(stderr, "    -allreads              ...all reads, regardless of deletion status\n");
  fprintf(stderr, "    -decoded               ...quality as integers ('20 21 19')\n");
  fprintf(stderr, "    -clear <clr>           ...in clear range <clr>, default=LATEST\n");
  fprintf(stderr, "  -dumpfrg               extract LIB, FRG and LKG messages\n");
  fprintf(stderr, "    -donotfixmates         ...only extract the fragments given, do not add in\n");
  fprintf(stderr, "                              missing mated reads\n");
  fprintf(stderr, "    -clear <clr>           ...use clear range <clr>, default=ORIG\n");
  fprintf(stderr, "    -format2               ...extract using frg format version 2\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  (1) - must have a -dump option, e.g., -uid file -tabular -dumpfragments some.gkpStore\n");
  fprintf(stderr, "\n");
}

#define DUMP_NOTHING     0
#define DUMP_INFO        1
#define DUMP_BATCHES     2
#define DUMP_LIBRARIES   3
#define DUMP_FRAGMENTS   4
#define DUMP_FASTA       5
#define DUMP_FRG         6
#define DUMP_LASTFRG     7

int
main(int argc, char **argv) {

  //  Options for everyone.  Everyone needs a GateKeeper!
  //
  char            *gkpStoreName    = NULL;

  //  Options for appending or creating:
  //
  int              append          = 0;
  int              outputExists    = 0;
  int              check_qvs       = 1;
  char            *vectorClearFile = NULL;
  int              assembler       = AS_ASSEMBLER_GRANDE;
  int              nerrs           = 0;   // Number of errors in current run
  int              maxerrs         = 1;   // Number of errors allowed before we punt
  int              firstFileArg    = 0;

  //  Options for partitioning
  //
  char            *partitionFile = NULL;

  //  Options for dumping:
  //
  CDS_IID_t        begIID            = 0;
  CDS_IID_t        endIID            = 2000000000;  //  I hope I never see an assembly with 2 billion IIDs!
  char            *uidFileName       = NULL;
  char            *iidFileName       = NULL;
  int              dump              = DUMP_NOTHING;
  int              dumpTabular       = 0;
  int              dumpWithSequence  = 0;
  int              dumpFastaAllReads = 0;
  int              dumpClear         = AS_READ_CLEAR_UNTRIM;
  int              dumpFRGClear      = AS_READ_CLEAR_ORIG;
  int              dumpFastaClear    = AS_READ_CLEAR_LATEST;
  int              dumpFastaQuality  = 0;
  int              doNotFixMates     = 0;
  int              dumpFormat        = 1;
  int              dumpRandMateLib   = 0;
  int              dumpRandMateNum   = 0;
  char            *iidToDump         = NULL;



  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-a") == 0) {
      append = 1;
    } else if (strcmp(argv[arg], "-b") == 0) {
      begIID = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-e") == 0) {
      maxerrs = atoi(argv[++arg]);
      endIID  = maxerrs;
    } else if (strcmp(argv[arg], "-h") == 0) {
      err++;
    } else if (strcmp(argv[arg], "-o") == 0) {
      gkpStoreName = argv[++arg];
    } else if (strcmp(argv[arg], "-v") == 0) {
      vectorClearFile = argv[++arg];
    } else if (strcmp(argv[arg], "-G") == 0) {
      assembler = AS_ASSEMBLER_GRANDE;
    } else if (strcmp(argv[arg], "-T") == 0) {
      assembler = AS_ASSEMBLER_OBT;
    } else if (strcmp(argv[arg], "-Q") == 0) {
      check_qvs = 0;

    } else if (strcmp(argv[arg], "-P") == 0) {
      partitionFile = argv[++arg];

      //  Except for -b and -e, all the dump options are below.

    } else if (strcmp(argv[arg], "-tabular") == 0) {
      dumpTabular = 1;
    } else if (strcmp(argv[arg], "-uid") == 0) {
      uidFileName = argv[++arg];
    } else if (strcmp(argv[arg], "-iid") == 0) {
      iidFileName = argv[++arg];
    } else if (strcmp(argv[arg], "-randommated") == 0) {
      dumpRandMateLib = atoi(argv[++arg]);
      dumpRandMateNum = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-dumpinfo") == 0) {
      dump = DUMP_INFO;
    } else if ((strcmp(argv[arg], "-lastfragiid") == 0) ||
               (strcmp(argv[arg], "-L") == 0)) {
      //  -L for compatability
      dump = DUMP_LASTFRG;
    } else if (strcmp(argv[arg], "-dumpbatch") == 0) {
      dump = DUMP_BATCHES;
    } else if (strcmp(argv[arg], "-dumplibraries") == 0) {
      dump = DUMP_LIBRARIES;
    } else if (strcmp(argv[arg], "-dumpfragments") == 0) {
      dump = DUMP_FRAGMENTS;
    } else if (strcmp(argv[arg], "-withsequence") == 0) {
      dumpWithSequence = 1;
    } else if (strcmp(argv[arg], "-clear") == 0) {
      dumpClear      = AS_PER_decodeClearRangeLabel(argv[++arg]);
      dumpFRGClear   = dumpClear;
      dumpFastaClear = dumpClear;
    } else if (strcmp(argv[arg], "-format2") == 0) {
      dumpFormat = 2;
    } else if (strcmp(argv[arg], "-dumpfastaseq") == 0) {
      dump = DUMP_FASTA;
      dumpFastaQuality = 0;
    } else if (strcmp(argv[arg], "-dumpfastaqlt") == 0) {
      dump = DUMP_FASTA;
      dumpFastaQuality = 1;
    } else if (strcmp(argv[arg], "-allreads") == 0) {
      dumpFastaAllReads = 1;
    } else if (strcmp(argv[arg], "-decoded") == 0) {
      dumpFastaQuality = 2;
    } else if (strcmp(argv[arg], "-dumpfrg") == 0) {
      dump = DUMP_FRG;
    } else if (strcmp(argv[arg], "-donotfixmates") == 0) {
      doNotFixMates = 1;

      //  End of dump options

    } else if (strcmp(argv[arg], "--") == 0) {
      firstFileArg = arg++;
      arg = argc;
    } else if (argv[arg][0] == '-') {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    } else {
      firstFileArg = arg;

      if ((dump != DUMP_NOTHING) || (partitionFile))
        gkpStoreName = argv[arg];

      arg = argc;
    }

    arg++;
  }
  if ((err) || (gkpStoreName == NULL)) {
    usage(argv[0]);
    exit(1);
  }
  
  outputExists = testOpenGateKeeperStore(gkpStoreName, FALSE);

  //  Construct an IID map of objects we care about.
  //
  if (uidFileName || iidFileName) {
    GateKeeperStore *gkp      = openGateKeeperStore(gkpStoreName, FALSE);
    CDS_IID_t        lastElem = getLastElemFragStore(gkp) + 1;
    FILE            *F        = NULL;
    char             L[1024];

    if (gkp == NULL) {
      fprintf(stderr, "Failed to open %s\n", gkpStoreName);
      exit(1);
    }

    iidToDump = (char *)safe_calloc(lastElem, sizeof(char));

    if (iidFileName) {
      errno = 0;
      if (strcmp(iidFileName, "-") == 0)
        F = stdin;
      else
        F = fopen(iidFileName, "r");
      if (errno) {
        fprintf(stderr, "%s: Couldn't open -iid file '%s': %s\n", argv[0], iidFileName, strerror(errno));
        exit(1);
      }
      fgets(L, 1024, F);
      while (!feof(F)) {
        CDS_UID_t      iid = STR_TO_UID(L, 0L, 10);
        if (iid >= lastElem)
          fprintf(stderr, "%s: IID "F_UID" too big, ignored.\n", argv[0], iid);
        else
          iidToDump[iid]++;
        fgets(L, 1024, F);
      }
      if (F != stdin)
        fclose(F);
    }

    if (uidFileName) {
      errno = 0;
      if (strcmp(uidFileName, "-") == 0)
        F = stdin;
      else
        F = fopen(uidFileName, "r");
      if (errno) {
        fprintf(stderr, "%s: Couldn't open -uid file '%s': %s\n", argv[0], uidFileName, strerror(errno));
        exit(1);
      }
      fgets(L, 1024, F);
      while (!feof(F)) {
        CDS_UID_t      uid = STR_TO_UID(L, 0L, 10);
        PHashValue_AS  value;

        if (HASH_FAILURE == getGatekeeperUIDtoIID(gkp, uid, &value)) {
          fprintf(stderr, "%s: UID "F_UID" doesn't exist, ignored.\n", argv[0], uid);
        } else {
          if (value.IID >= lastElem)
            fprintf(stderr, "%s: UID "F_UID" is IID "F_IID", and that's too big, ignored.\n", argv[0], uid, value.IID);
          else
            iidToDump[value.IID]++;
        }
        fgets(L, 1024, F);
      }
      if (F != stdin)
        fclose(F);
    }

    closeGateKeeperStore(gkp);
  }

  if (dumpRandMateNum > 0) {
    GateKeeperStore *gkp        = openGateKeeperStore(gkpStoreName, FALSE);
    CDS_IID_t        lastElem   = getLastElemFragStore(gkp) + 1;

    //  No way to know (currently) how many reads are in a library,
    //  without actually counting it.  We just allocate enough space
    //  to hold every read.

    uint32          *candidatesA   = (uint32 *)safe_malloc(lastElem * sizeof(uint32));
    uint32          *candidatesB   = (uint32 *)safe_malloc(lastElem * sizeof(uint32));
    uint32           candidatesLen = 0;

    if (gkp == NULL) {
      fprintf(stderr, "Failed to open %s\n", gkpStoreName);
      exit(1);
    }

    iidToDump = (char *)safe_calloc(lastElem, sizeof(char));

    fragRecord   *fr = new_fragRecord();
    FragStream   *fs = openFragStream(gkp, FRAG_S_INF);
    StoreStat     stat;

    int           i;

    statsStore(gkp->frg, &stat);

    if (begIID < stat.firstElem)
      begIID = stat.firstElem;
    if (stat.lastElem < endIID)
      endIID = stat.lastElem;

    resetFragStream(fs, begIID, endIID);

    //  Scan the whole fragstore, looking for mated reads in the
    //  correct library, and save the lesser of the two reads.
    //
    while (nextFragStream(fs, fr)) {
      if ((getFragRecordLibraryIID(fr) == dumpRandMateLib) &&
          (getFragRecordMateIID(fr) > 0) &&
          (getFragRecordIID(fr) < getFragRecordMateIID(fr))) {
        candidatesA[candidatesLen] = getFragRecordIID(fr);
        candidatesB[candidatesLen] = getFragRecordMateIID(fr);
        candidatesLen++;
      }
    }

    //  Now pick N reads from our list of candidates, and let the
    //  other guys fill in the mates.
    //
    srand48(time(NULL));
    for (i=0; (i<dumpRandMateNum) && (candidatesLen > 0); i++) {
      int  x = lrand48() % candidatesLen;
      iidToDump[candidatesA[x]] = 1;
      iidToDump[candidatesB[x]] = 1;
      candidatesLen--;
      candidatesA[x] = candidatesA[candidatesLen];
      candidatesB[x] = candidatesB[candidatesLen];
    }

    closeFragStream(fs);
    closeGateKeeperStore(gkp);

    safe_free(candidatesA);
    safe_free(candidatesB);
  }




  if (dump != DUMP_NOTHING) {
    if (outputExists == 0) {
      fprintf(stderr,"* Gatekeeper Store %s doesn't exist, but you want to dump it.  Exit.\n", gkpStoreName);
      exit(1);
    }

    switch (dump) {
      case DUMP_INFO:
        dumpGateKeeperInfo(gkpStoreName);
        break;
      case DUMP_BATCHES:
        dumpGateKeeperBatches(gkpStoreName, begIID, endIID, iidToDump, dumpTabular);
        break;
      case DUMP_LIBRARIES:
        dumpGateKeeperLibraries(gkpStoreName, begIID, endIID, iidToDump, dumpTabular);
        break;
      case DUMP_FRAGMENTS:
        dumpGateKeeperFragments(gkpStoreName, begIID, endIID, iidToDump,
                                dumpWithSequence,
                                dumpClear,
                                dumpTabular);
        break;
      case DUMP_FASTA:
        dumpGateKeeperAsFasta(gkpStoreName, begIID, endIID, iidToDump,
                              dumpFastaAllReads,
                              dumpFastaClear,
                              dumpFastaQuality);
        break;
      case DUMP_FRG:
        dumpGateKeeperAsFRG(gkpStoreName, dumpFormat, begIID, endIID, iidToDump,
                            doNotFixMates,
                            dumpFRGClear);
        break;
      case DUMP_LASTFRG:
        {
          GateKeeperStore *gkp = openGateKeeperStore(gkpStoreName, FALSE);
          fprintf(stdout, "Last frag in store is iid = "F_S64"\n", getLastElemFragStore(gkp));
          closeGateKeeperStore(gkp);
        }
        break;
      default:
        break;
    }

    exit(0);
  }






  if (partitionFile) {
    Build_Partition(gkpStoreName, partitionFile, FRAG_S_ALL);
    exit(0);
  }





  //  This sure is sloppy (and slow) but only used on format 1.
  //  Load in any vector info.
  //
  if (vectorClearFile) {
    char          line[256];
    int           nlines  = 0;
    int           nupdate = 0;

    errno = 0;
    FILE   *v = fopen(vectorClearFile, "r");
    if (errno) {
      fprintf(stderr, "%s: couldn't open '%s' to read vector clear ranges: %s\n",
              argv[0], vectorClearFile, strerror(errno));
      exit(1);
    }

    gkpStore = openGateKeeperStore(gkpStoreName, TRUE);

    fgets(line, 256, v);
    while (!feof(v)) {
      PHashValue_AS  value;
      char          *pine = line;

      CDS_UID_t uid = STR_TO_UID(pine, &pine, 10);
      int       l   = strtol(pine, &pine, 10);
      int       r   = strtol(pine, &pine, 10);

      if (uid == 0)
        fprintf(stderr, "unexpected line: %s", line);

      if (HASH_FAILURE != getGatekeeperUIDtoIID(gkpStore, uid, &value)) {
        fragRecord    fr;
        uint32        iid = value.IID;

        getFrag(gkpStore, iid, &fr, FRAG_S_INF);

        //  Assume they are base-based.
        fr.gkfr.hasVectorClear = 1;
        if (l < r) {
          fr.gkfr.clearBeg[AS_READ_CLEAR_VEC] = l - 1;
          fr.gkfr.clearEnd[AS_READ_CLEAR_VEC] = r;
        } else {
          fr.gkfr.clearBeg[AS_READ_CLEAR_VEC] = r - 1;
          fr.gkfr.clearEnd[AS_READ_CLEAR_VEC] = l;
        }

        nupdate++;
        setFrag(gkpStore, iid, &fr);
      }

      nlines++;
      fgets(line, 256, v);
    }

    closeGateKeeperStore(gkpStore);

    fclose(v);

    fprintf(stderr, "in %d lines, updated %d fragments.\n", nlines, nupdate);

    exit(0);
  }



   

  if ((append == 0) && (outputExists == 1)) {
    fprintf(stderr,"* Gatekeeper Store %s exists and append flag not supplied.  Exit.\n", gkpStoreName);
    exit(1);
  }
  if ((append == 1) && (outputExists == 0)) {
    //  Silently switch over to create
    append = 0;
  }
  if ((append == 0) && (outputExists == 1)) {
    fprintf(stderr,"* Gatekeeper Store %s exists, but not told to append.  Exit.\n", gkpStoreName);
    exit(1);
  }

  if (append)
    gkpStore = openGateKeeperStore(gkpStoreName, TRUE);
  else
    gkpStore = createGateKeeperStore(gkpStoreName);

  for (; firstFileArg < argc; firstFileArg++) {
    FILE            *inFile            = NULL;
    GenericMesg     *pmesg             = NULL;
    int              fileIsCompressed  = 0;

    if        (strcmp(argv[firstFileArg] + strlen(argv[firstFileArg]) - 3, ".gz") == 0) {
      char  cmd[1024];
      sprintf(cmd, "gzip -dc %s", argv[firstFileArg]);
      errno = 0;
      inFile = popen(cmd, "r");
      fileIsCompressed = 1;
    } else if (strcmp(argv[firstFileArg] + strlen(argv[firstFileArg]) - 4, ".bz2") == 0) {
      char  cmd[1024];
      sprintf(cmd, "bzip2 -dc %s", argv[firstFileArg]);
      errno = 0;
      inFile = popen(cmd, "r");
      fileIsCompressed = 1;
    } else {
      errno = 0;
      inFile = fopen(argv[firstFileArg], "r");
    }

    if (errno)
      fprintf(stderr, "%s: failed to open input '%s': %s\n", argv[0], argv[firstFileArg], strerror(errno)), exit(1);
    if (inFile == NULL)
      fprintf(stderr, "%s: failed to open input '%s': (returned null pointer)\n", argv[0], argv[firstFileArg]), exit(1);

    while (EOF != ReadProtoMesg_AS(inFile, &pmesg)) {
      int success = GATEKEEPER_SUCCESS;

      if (pmesg->t == MESG_BAT) {
        success = Check_BatchMesg((BatchMesg *)pmesg->m);
      } else if (pmesg->t == MESG_DST) {
        success = Check_DistanceMesg((DistanceMesg *)pmesg->m);
      } else if (pmesg->t == MESG_LIB) {
        success = Check_LibraryMesg((LibraryMesg *)pmesg->m);
      } else if (pmesg->t == MESG_FRG) {
        success = Check_FragMesg((FragMesg *)pmesg->m, check_qvs, assembler);
      } else if (pmesg->t == MESG_LKG) {
        success = Check_LinkMesg((LinkMesg *)pmesg->m);
      } else if (pmesg->t == MESG_VER) {
        //  Ignore
      } else {
        fprintf(stderr,"GKP Error: Unknown message with type %s.\n", MessageTypeName[pmesg->t]);
        success = GATEKEEPER_FAILURE;
      }

      if (success != GATEKEEPER_SUCCESS) {
        fprintf(stderr,"GKP Error: at line %d:\n", GetProtoLineNum_AS());
        WriteProtoMesg_AS(stderr,pmesg);
        if (pmesg->t == MESG_BAT) {
          fprintf(stderr, "GKP Error: Invalid BAT message, can't continue.\n");
          return(GATEKEEPER_FAILURE);
        }
        nerrs++;
      }

      if (nerrs >= maxerrs) {
        fprintf(stderr, "GKP Error: Too many errors (%d), can't continue.\n", nerrs);
        return(GATEKEEPER_FAILURE);
      }
    }

    if (fileIsCompressed) {
      errno = 0;
      pclose(inFile);
      if (errno)
        fprintf(stderr, "%s: WARNING!  Failed to close '%s': %s\n", argv[0], argv[firstFileArg], strerror(errno));
    } else {
      fclose(inFile);
    }
  }

  closeGateKeeperStore(gkpStore);

  fprintf(stderr, "GKP finished with %d errors.\n", nerrs);

  exit(0);
}