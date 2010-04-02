
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute. All rights reserved.
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

const char *mainid = "$Id: overlapStore.c,v 1.25 2010-01-29 07:15:25 brianwalenz Exp $";

#include "overlapStore.h"
#include "AS_OVS_overlap.h"   //  Just to know the sizes of structs
#include "AS_PER_gkpStore.h"  //  Just to know clear region labels

#include <ctype.h>
#include <unistd.h>  //  sysconf()

int
main(int argc, char **argv) {
  uint32          operation   = OP_NONE;
  char           *storeName   = NULL;
  char           *gkpName     = NULL;
  uint32          clearRegion = AS_READ_CLEAR_ERROR;
  uint32          dumpBinary  = FALSE;
  double          dumpERate   = 100.0;
  uint32          dumpType    = 0;
  uint32          bgnIID      = 0;
  uint32          endIID      = 1000000000;
  uint32          qryIID      = 0;
  uint64          memoryLimit = 512 * 1024 * 1024;
  uint32          nThreads    = 4;
  uint32          doFilterOBT = 0;
  uint32          fileListLen = 0;
  uint32          fileListMax = 10 * 1024;  //  If you run more than 10,000 overlapper jobs, you'll die.
  char          **fileList    = (char **)safe_malloc(sizeof(char *) * fileListMax);
  Ovl_Skip_Type_t ovlSkipOpt  = ALL;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {

    if        (strcmp(argv[arg], "-c") == 0) {
      if (storeName)
        fprintf(stderr, "ERROR: only one of -c, -m, -d, -q, -s, -S or -u may be supplied.\n"), err++;
      storeName   = argv[++arg];
      operation   = OP_BUILD;

    } else if (strcmp(argv[arg], "-m") == 0) {
      if (storeName)
        fprintf(stderr, "ERROR: only one of -c, -m, -d, -q, -s, -S or -u may be supplied.\n"), err++;
      storeName   = argv[++arg];
      operation   = OP_MERGE;

    } else if (strcmp(argv[arg], "-d") == 0) {
      if (storeName)
        fprintf(stderr, "ERROR: only one of -c, -m, -d, -q, -s, -S or -u may be supplied.\n"), err++;
      storeName   = argv[++arg];
      operation   = OP_DUMP;

    } else if (strcmp(argv[arg], "-p") == 0) {
      if (storeName)
        fprintf(stderr, "ERROR: only one of -c, -m, -d, -q, -s, -S or -u may be supplied.\n"), err++;
      qryIID      = atoi(argv[++arg]);
      storeName   = argv[++arg];
      gkpName     = argv[++arg];
      clearRegion = gkStore_decodeClearRegionLabel(argv[++arg]);
      operation   = OP_DUMP_PICTURE;

    } else if (strcmp(argv[arg], "-s") == 0) {
      if (storeName)
        fprintf(stderr, "ERROR: only one of -c, -m, -d, -q, -s, -S or -u may be supplied.\n"), err++;
      storeName   = argv[++arg];
      operation   = OP_STATS_DUMP;

    } else if (strcmp(argv[arg], "-S") == 0) {
      if (storeName)
        fprintf(stderr, "ERROR: only one of -c, -m, -d, -q, -s, -S or -u may be supplied.\n"), err++;
      storeName   = argv[++arg];
      operation   = OP_STATS_REBUILD;

    } else if (strcmp(argv[arg], "-u") == 0) {
      if (storeName)
        fprintf(stderr, "ERROR: only one of -c, -m, -d, -q, -s, -S or -u may be supplied.\n"), err++;
      storeName   = argv[++arg];
      operation   = OP_UPDATE_ERATES;

    } else if (strcmp(argv[arg], "-t") == 0) {
      nThreads    = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-E") == 0) {
      dumpERate = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-d5") == 0) {
      dumpType |= DUMP_5p;

    } else if (strcmp(argv[arg], "-d3") == 0) {
      dumpType |= DUMP_3p;

    } else if (strcmp(argv[arg], "-dC") == 0) {
      dumpType |= DUMP_CONTAINS;

    } else if (strcmp(argv[arg], "-dc") == 0) {
      dumpType |= DUMP_CONTAINED;

    } else if (strcmp(argv[arg], "-B") == 0) {
      dumpBinary = TRUE;

    } else if (strcmp(argv[arg], "-b") == 0) {
      bgnIID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      endIID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-q") == 0) {
      if (storeName)
        fprintf(stderr, "ERROR: only one of -c, -m, -d, -q, -s, -S or -u may be supplied.\n"), err++;
      bgnIID    = atoi(argv[++arg]);
      endIID    = bgnIID;
      qryIID    = atoi(argv[++arg]);
      storeName = argv[++arg];
      operation = OP_DUMP;

    } else if (strcmp(argv[arg], "-O") == 0) {
      doFilterOBT++;

    } else if (strcmp(argv[arg], "-M") == 0) {
      memoryLimit  = atoi(argv[++arg]);  //  convert first, then multiply so we don't
      memoryLimit *= 1024 * 1024;        //  overflow whatever type atoi() is.

    } else if (strcmp(argv[arg], "-g") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-L") == 0) {
      char *line;

      //  The next arg is a file with the list of files to use
      errno = 0;
      FILE *F = fopen(argv[++arg], "r");
      if (errno)
        fprintf(stderr, "Can't open '%s': %s\n", argv[arg], strerror(errno)), exit(1);

      line = (char *)safe_malloc(sizeof(char) * FILENAME_MAX);
      fgets(line, FILENAME_MAX, F);
      while (!feof(F)) {
        chomp(line);
        fileList[fileListLen++] = line;
        if (fileListLen >= fileListMax)
          fprintf(stderr, "Too many input files, increase fileListMax.\n"), exit(1);
        line = (char *)safe_malloc(sizeof(char) * FILENAME_MAX);
        fgets(line, FILENAME_MAX, F);
      }
      safe_free(line);
      fclose(F);
    
    } else if (strcmp(argv[arg], "-i") == 0) {      
      switch (atoi(argv[++arg])) {
         case 0:
            ovlSkipOpt = NONE;
            break;
         case 1:
            ovlSkipOpt = ALL;
            break;
         case 2:
            ovlSkipOpt = INTERNAL;
            break;
         default:
            fprintf(stderr, "%s: unknown overlap ignore option '%s'. Must be one of 0, 1, or 2.\n", argv[0], argv[arg]);
            err++;
            break;
      }
    
    } else if ((argv[arg][0] == '-') && (argv[arg][1] != 0)) {
      fprintf(stderr, "%s: unknown option '%s'.\n", argv[0], argv[arg]);
      err++;

    } else {
      //
      //  Assume it's an input file
      //
      fileList[fileListLen++] = argv[arg];
      if (fileListLen >= fileListMax)
        fprintf(stderr, "Too many input files, increase fileListMax.\n"), exit(1);
    }
    arg++;
  }
  if ((operation == OP_NONE) || (storeName == NULL) || (err)) {
    fprintf(stderr, "usage: %s -c storeName [-M x (MB)] [-t threads] [-g gkpStore] [-L list-of-ovl-files] ovl-file ...\n", argv[0]);
    fprintf(stderr, "       %s -m storeName mergeName\n", argv[0]);
    fprintf(stderr, "       %s -d storeName [-B] [-E erate] [-b beginIID] [-e endIID]\n", argv[0]);
    fprintf(stderr, "       %s -s storeName\n", argv[0]);
    fprintf(stderr, "       %s -S storeName -g gkpStore\n", argv[0]);
    fprintf(stderr, "       %s -q iid iid storeName\n", argv[0]);
    fprintf(stderr, "       %s -p iid storeName gkpStore clr\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -c  create a new store, fails if the store exists\n");
    fprintf(stderr, "  -m  merge store mergeName into store storeName\n");
    fprintf(stderr, "  -d  dump a store\n");
    fprintf(stderr, "  -s  dump statistics about a store\n");
    fprintf(stderr, "  -S  recompute overlap statistics\n");
    fprintf(stderr, "  -q  report the a,b overlap, if it exists.\n");
    fprintf(stderr, "  -p  dump a picture of overlaps to fragment 'iid', using clear region 'clr'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "CREATION\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -O           Filter overlaps for OBT.\n");
    fprintf(stderr, "  -M x         Use 'x'MB memory for sorting overlaps.\n");
    fprintf(stderr, "  -t t         Use 't' threads for sorting overlaps.\n");
    fprintf(stderr, "  -L f         Read overlaps from files listed in 'f'.\n");
    fprintf(stderr, "  -I F         Ignore the ovls for the reads listed in 'f'.\n");
    fprintf(stderr, "  -i x         Where x is either 0, 1, or 2. Requires -I option. 1 (default) - Delete all overlaps. 2 - Delete only the overlaps between reads listed in 'f'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "MERGING\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -m storeName mergeName   Merge the store 'mergeName' into 'storeName'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "DUMPING\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -B                Dump the store as binary, suitable for input to create a new store.\n");
    fprintf(stderr, "  -E erate          Dump only overlaps <= erate error.\n");
    fprintf(stderr, "  -d5               Dump only overlaps off the 5' end of the A frag.\n");
    fprintf(stderr, "  -d3               Dump only overlaps off the 3' end of the A frag.\n");
    fprintf(stderr, "  -dC               Dump only overlaps that are contained in the A frag (B contained in A).\n");
    fprintf(stderr, "  -dc               Dump only overlaps that are containing the A frag (A contained in B).\n");
    fprintf(stderr, "  -b beginIID       Start dumping at 'beginIID'.\n");
    fprintf(stderr, "  -e endIID         Stop dumping after 'endIID'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "DUMPING PICTURES\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -p iid storeName gkpStore clr\n");
    fprintf(stderr, "                    clr is usually OBTINITIAL for obtStore.\n");
    fprintf(stderr, "                    clr is usually OBTCHIMERA for ovlStore when OBT is used.\n");
    fprintf(stderr, "                    clr is usually CLR        for ovlStore when OBT is not used.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Limited to %d open files.\n", sysconf(_SC_OPEN_MAX));
    fprintf(stderr, "\n");
    fprintf(stderr, "OVSoverlap     %d bytes\n", sizeof(OVSoverlap));
    fprintf(stderr, "OVSoverlapINT  %d bytes\n", sizeof(OVSoverlapINT));
    fprintf(stderr, "OVSoverlapDAT  %d bytes\n", sizeof(OVSoverlapDAT));
    fprintf(stderr, "               %d %d %d\n", sizeof(struct OVSoverlapOVL), sizeof(struct OVSoverlapMER), sizeof(struct OVSoverlapOBT));
    fprintf(stderr, "AS_OVS_NWORDS  %d\n", AS_OVS_NWORDS);
    exit(1);
  }
  if ((fileListLen == 0) && (operation == OP_BUILD)) {
    fprintf(stderr, "No input files?\n");
    exit(1);
  }
  if (dumpType == 0)
    dumpType = DUMP_5p | DUMP_3p | DUMP_CONTAINED | DUMP_CONTAINS;

  switch (operation) {
    case OP_BUILD:
      buildStore(storeName, gkpName, memoryLimit, nThreads, doFilterOBT, fileListLen, fileList, ovlSkipOpt);
      break;
    case OP_MERGE:
      mergeStore(storeName, fileList[0]);
      break;
    case OP_DUMP:
      dumpStore(storeName, dumpBinary, dumpERate, dumpType, bgnIID, endIID, qryIID);
      break;
    case OP_DUMP_PICTURE:
      dumpPicture(storeName, gkpName, clearRegion, dumpERate, dumpType, qryIID);
      break;
    case OP_STATS_DUMP:
      dumpStats(storeName);
      break;
    case OP_STATS_REBUILD:
      rebuildStats(storeName, gkpName);
      break;
    case OP_UPDATE_ERATES:
      updateErates(storeName, fileList[0]);
      break;
    default:
      break;
  }


  exit(0);
}
