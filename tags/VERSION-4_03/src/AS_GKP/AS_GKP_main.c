
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

static char const *rcsid = "$Id: AS_GKP_main.c,v 1.55 2007-08-31 21:06:16 brianwalenz Exp $";

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
#include "AS_UTL_version.h"
#include "AS_UTL_fileIO.h"
#include "AS_MSG_pmesg.h"
#include "AS_GKP_include.h"

GateKeeperStore *gkpStore  = NULL;
FILE            *errorFP   = NULL;

static
void
usage(char *filename, int longhelp) {
  fprintf(stdout, "usage1: %s -o gkpStore [append/create options] <input.frg> <input.frg> ...\n", filename);
  fprintf(stdout, "usage2: %s -P partitionfile gkpStore\n", filename);
  fprintf(stdout, "usage3: %s [dump-options] gkpStore\n", filename);
  fprintf(stdout, "\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "The first usage will append to or create a GateKeeper store:\n");
  fprintf(stdout, "  -a                     append to existing tore\n");
  fprintf(stdout, "  -e <errorThreshhold>   set error threshhold\n");
  fprintf(stdout, "  -o <gkpStore>          append to or create gkpStore\n");
  fprintf(stdout, "  -G                     gatekeeper for assembler Grande (default)\n");
  fprintf(stdout, "  -T                     gatekeeper for assembler Grande with Overlap Based Trimming\n");
  fprintf(stdout, "  -Q                     don't check quality-based data quality\n");
  fprintf(stdout, "  -E <error.frg>         write errors to this file\n");
  fprintf(stdout, "  -D                     do NOT limit insert size stddev to 10% of the mean\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "  -v <vector-info>       load vector clear ranges into each read.\n");
  fprintf(stdout, "                         MUST be done on an existing, complete store.\n");
  fprintf(stdout, "                         example: -a -v vectorfile -o that.gkpStore\n");
  fprintf(stdout, "                         format: 'UID vec-clr-begin vec-clr-end'\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "The second usage will partition an existing store, allowing\n");
  fprintf(stdout, "the entire store partition to be loaded into memory.\n");
  fprintf(stdout, "  -P <partitionfile>     a list of (partition fragiid)\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "The third usage will dump the contents of a GateKeeper store.\n");
  fprintf(stdout, "  [selection of what objects to dump]\n");
  fprintf(stdout, "  -b <begin-iid>         dump starting at this batch, library or read (1)\n");
  fprintf(stdout, "  -e <ending-iid>        dump stopping after this iid (1)\n");
  fprintf(stdout, "  -uid <uid-file>        dump only objects listed in 'uid-file' (1)\n");
  fprintf(stdout, "  -iid <iid-file>        dump only objects listed in 'iid-file' (1)\n");
  fprintf(stdout, "  -randommated <lib> <n> pick n mates (2n frags) at random from library lib,\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "  [how to dump it]\n");
  fprintf(stdout, "  -tabular               dump info, batches, libraries or fragments in a tabular\n");
  fprintf(stdout, "                         format (for -dumpinfo, -dumpbatch, -dumplibraries,\n");
  fprintf(stdout, "                         and -dumpfragments, ignores -withsequence and -clear)\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "  [format of dump]\n");
  fprintf(stdout, "  -dumpinfo              print information on the store\n");
  fprintf(stdout, "    -lastfragiid         just print the last IID in the store\n");
  fprintf(stdout, "  -dumpbatch             dump all batch records\n");
  fprintf(stdout, "  -dumplibraries         dump all library records\n");
  fprintf(stdout, "  -dumpfragments         dump fragment info, no sequence\n");
  fprintf(stdout, "    -withsequence          ...and include sequence\n");
  fprintf(stdout, "    -clear <clr>           ...in clear range <clr>, default=UNTRIM\n");
  fprintf(stdout, "  -dumpfasta[seq|qlt]    dump fragment sequence or quality, as fasta format\n");
  fprintf(stdout, "    -allreads              ...all reads, regardless of deletion status\n");
  fprintf(stdout, "    -decoded               ...quality as integers ('20 21 19')\n");
  fprintf(stdout, "    -clear <clr>           ...in clear range <clr>, default=LATEST\n");
  fprintf(stdout, "  -dumpfrg               extract LIB, FRG and LKG messages\n");
  fprintf(stdout, "    -donotfixmates         ...only extract the fragments given, do not add in\n");
  fprintf(stdout, "                              missing mated reads\n");
  fprintf(stdout, "    -clear <clr>           ...use clear range <clr>, default=ORIG\n");
  fprintf(stdout, "    -format2               ...extract using frg format version 2\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "  (1) - must have a -dump option, e.g., -uid file -tabular -dumpfragments some.gkpStore\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "\n");
  if (longhelp == 0) {
    fprintf(stdout, "Use '-h' to get a discussion of what gatekeeper is.\n");
    fprintf(stdout, "\n");
  } else {
    fprintf(stdout, "The Gatekeeper ensures that data entering the assembly system meets\n");
    fprintf(stdout, "the data specification (see GateKeeper design document).  It is also\n");
    fprintf(stdout, "used for examining and partitioning the assembler data store.\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "Each input message is checked for semantic consistency as described in\n");
    fprintf(stdout, "the defining document for that stage.  Messages containing a UID are\n");
    fprintf(stdout, "converted to a UID,IID pair -- the assembler modules require\n");
    fprintf(stdout, "consecutive IID beginning at 1 for efficient indexing of internal and\n");
    fprintf(stdout, "disk-based data structures; gatekeeper performs this task.  Finally,\n");
    fprintf(stdout, "each message is inserted into the assembly data store.\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "The GateKeeper succeeds if it consumes its entire input with less than\n");
    fprintf(stdout, "a specified number of errors (the -e option).  Upon successful exit,\n");
    fprintf(stdout, "the store reflects all of the records that were successfully read.\n");
    fprintf(stdout, "Unsuccessful records are reported to stderr, along with a brief\n");
    fprintf(stdout, "explanation of the problem.\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "If unsuccessful, the store is partially updated.\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "Resoure Requirements\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "The key gatekeeper data structures are in-memory copies of its store.\n");
    fprintf(stdout, "This store should scale linearly with the number of fragments.\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "No formal benchmarking of the gatekeeper has been performed to date.\n");
    fprintf(stdout, "However, each LKG message requires four random disk accesses -- two to\n");
    fprintf(stdout, "read the linked fragment records, and two two write the updated\n");
    fprintf(stdout, "fragment records.  This can cause problems when gatekeeper is run over\n");
    fprintf(stdout, "low-performance or heavily used NFS mount points.\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "3. General Pre Conditions\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "a) Each input UID must be unique.  A new message with a duplicate UID\n");
    fprintf(stdout, "will be rejected as an error.\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "b) Any object referred to in a message must be defined:\n");
    fprintf(stdout, "def-before-ref.  Likewise, before an object can be deleted, all\n");
    fprintf(stdout, "references to it must be removed: unref-before-undef.\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "c) The input specification is defined elsewhere.\n");
    fprintf(stdout, "\n");
  }
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
  int              append             = 0;
  int              outputExists       = 0;
  int              check_qvs          = 1;
  char            *vectorClearFile    = NULL;
  int              assembler          = AS_ASSEMBLER_GRANDE;
  int              nerrs              = 0;   // Number of errors in current run
  int              maxerrs            = 1;   // Number of errors allowed before we punt
  int              firstFileArg       = 0;
  char            *errorFile          = NULL;
  int              believeInputStdDev = 0;

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

  gkpStore = NULL;
  errorFP  = stdout;

  int arg = 1;
  int err = 0;
  int hlp = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-a") == 0) {
      append = 1;
    } else if (strcmp(argv[arg], "-b") == 0) {
      begIID = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-e") == 0) {
      maxerrs = atoi(argv[++arg]);
      endIID  = maxerrs;
    } else if (strcmp(argv[arg], "-h") == 0) {
      hlp++;
      err++;
    } else if (strcmp(argv[arg], "-o") == 0) {
      gkpStoreName = argv[++arg];
    } else if (strcmp(argv[arg], "-v") == 0) {
      vectorClearFile = argv[++arg];
      firstFileArg    = 1;  // gets us around the input file sanity check, unused otherwise
    } else if (strcmp(argv[arg], "-G") == 0) {
      assembler = AS_ASSEMBLER_GRANDE;
    } else if (strcmp(argv[arg], "-T") == 0) {
      assembler = AS_ASSEMBLER_OBT;
    } else if (strcmp(argv[arg], "-Q") == 0) {
      check_qvs = 0;
    } else if (strcmp(argv[arg], "-E") == 0) {
      errorFile = argv[++arg];
    } else if (strcmp(argv[arg], "-D") == 0) {
      believeInputStdDev = 1;

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

      //  End of dump options, SECRET OPTIONS below

    } else if (strcmp(argv[arg], "--rebuildmap") == 0) {
      if ((arg + 2) == argc) {
         fprintf(stderr, "Must input a gkpStore and map file name for map re-generation.\n");
         exit(1);
      }
      
      rebuildMap(argv[arg+1], argv[arg+2]);
      exit(0);
    } else if (strcmp(argv[arg], "--rearrange") == 0) {
      //  Takes three args:  UID order file, oldStore, newStore
      //
      rearrangeStore(argv[arg+1], argv[arg+2], argv[arg+3]);
      exit(0);

    } else if ((strcmp(argv[arg], "--edit") == 0) ||
               (strcmp(argv[arg], "--testedit") == 0)) {
      //  Takes two args:  edit file, gkpStore
      //
      editStore(argv[arg+1], argv[arg+2], (strcmp(argv[arg], "--edit") == 0));
      exit(0);

      //  End of SECRET options

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

  if ((err) || (gkpStoreName == NULL) || (firstFileArg == 0)) {
    usage(argv[0], hlp);
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
        CDS_IID_t      iid = getGatekeeperUIDtoIID(gkp, uid, NULL);

        if (iid == 0)
          fprintf(stderr, "%s: UID "F_UID" doesn't exist, ignored.\n", argv[0], uid);
        else if (iid >= lastElem)
          fprintf(stderr, "%s: UID "F_UID" is IID "F_IID", and that's too big, ignored.\n", argv[0], uid, iid);
        else
          iidToDump[iid]++;

        fgets(L, 1024, F);
      }
      if (F != stdin)
        fclose(F);
    }

    closeGateKeeperStore(gkp);
  }

  if (vectorClearFile) {
    updateVectorClear(vectorClearFile, gkpStoreName);
    exit(0);
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

  if (errorFile) {
    errno = 0;
    errorFP = fopen(errorFile, "w");
    if (errno) {
      fprintf(stderr, "%s: cannot open error file '%s': %s\n", errorFile, strerror(errno));
      exit(1);
    }
  }

  for (; firstFileArg < argc; firstFileArg++) {
    FILE            *inFile            = NULL;
    GenericMesg     *pmesg             = NULL;
    int              fileIsCompressed  = 0;

    fprintf(stderr, "Starting file '%s' at line %d.\n", argv[firstFileArg], GetProtoLineNum_AS());

    AS_MSG_setFormatVersion(1);

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

      if (pmesg->t == MESG_ADT) {
        //  Ignore
      } else if (pmesg->t == MESG_BAT) {
        success = Check_BatchMesg((BatchMesg *)pmesg->m);
      } else if (pmesg->t == MESG_DST) {
        success = Check_DistanceMesg((DistanceMesg *)pmesg->m, believeInputStdDev);
      } else if (pmesg->t == MESG_LIB) {
        success = Check_LibraryMesg((LibraryMesg *)pmesg->m, believeInputStdDev);
      } else if (pmesg->t == MESG_FRG) {
        success = Check_FragMesg((FragMesg *)pmesg->m, check_qvs, assembler);
      } else if (pmesg->t == MESG_LKG) {
        success = Check_LinkMesg((LinkMesg *)pmesg->m);
      } else if (pmesg->t == MESG_VER) {
        //  Ignore
      } else {
        fprintf(errorFP,"# GKP Error: Unknown message with type %s.\n", MessageTypeName[pmesg->t]);
        success = GATEKEEPER_FAILURE;
      }

      if (success != GATEKEEPER_SUCCESS) {
        fprintf(errorFP,"# GKP Error: at line %d:\n", GetProtoLineNum_AS());
        WriteProtoMesg_AS(errorFP,pmesg);
        nerrs++;
      }

      if (nerrs >= maxerrs) {
        fprintf(errorFP, "# GKP Error: Too many errors (%d), can't continue.\n", nerrs);
        goto done;
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

 done:
  closeGateKeeperStore(gkpStore);

  if (errorFile)
    fclose(errorFP);

  fprintf(stderr, "GKP finished with %d errors.\n", nerrs);

  return(0);
}