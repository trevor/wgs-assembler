
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

static const char *rcsid = "$Id: overlapStore_build.c,v 1.25 2009-07-06 20:03:41 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <assert.h>
#include <time.h>
#include <sys/stat.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "AS_UTL_qsort_mt.h"
#include "AS_OBT_acceptableOverlap.h"
#include "AS_OVS_overlap.h"
#include "AS_OVS_overlapFile.h"
#include "AS_OVS_overlapStore.h"

#include "overlapStore.h"

int
OVSoverlap_sort(const void *a, const void *b) {
  OVSoverlap const *A = (OVSoverlap const *)a;
  OVSoverlap const *B = (OVSoverlap const *)b;
  if (A->a_iid   < B->a_iid)    return(-1);
  if (A->a_iid   > B->a_iid)    return(1);
  if (A->b_iid   < B->b_iid)    return(-1);
  if (A->b_iid   > B->b_iid)    return(1);
  if (A->dat.dat < B->dat.dat)  return(-1);
  if (A->dat.dat > B->dat.dat)  return(1);
  return(0);
}


static
off_t
sizeOfFile(const char *path) {
  struct stat  s;
  int          r;
  off_t        size = 0;

  errno = 0;
  r = stat(path, &s);
  if (errno) {
    fprintf(stderr, "Failed to stat() file '%s': %s\n", path, strerror(errno));
    exit(1);
  }

  //  gzipped files contain a file contents list, which we can
  //  use to get the uncompressed size.
  //
  //  gzip -l <file>
  //  compressed        uncompressed  ratio uncompressed_name
  //       14444               71680  79.9% up.tar
  //
  //  bzipped files have no contents and we just guess.

  if        (strcasecmp(path+strlen(path)-3, ".gz") == 0) {
    char   cmd[256];
    FILE  *F;

    sprintf(cmd, "gzip -l %s", path);
    F = popen(cmd, "r");
    fscanf(F, " %*s %*s %*s %*s ");
    fscanf(F, " %*d %lld %*s %*s ", &size);
    pclose(F);
  } else if (strcasecmp(path+strlen(path)-4, ".bz2") == 0) {
    size = s.st_size * 14 / 10;
  } else {
    size = s.st_size;
  }

  return(size);
}


void
writeToDumpFile(OVSoverlap          *overlap,
                BinaryOverlapFile  **dumpFile,
                uint32               dumpFileMax,
                uint64              *dumpLength,
                uint32               iidPerBucket,
                char                *storeName) {

  int df = overlap->a_iid / iidPerBucket;

  if (df >= dumpFileMax) {
    fprintf(stderr, "Too many bucket files.  Increase memory size, set the correct\n");
    fprintf(stderr, "number of fragment IIDs in the input, or split your overlaps\n");
    fprintf(stderr, "into multiple batches.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "The runCA option for increasing the memory size is ovlStoreMemory.\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  if (dumpFile[df] == NULL) {
    char name[FILENAME_MAX];
    sprintf(name, "%s/tmp.sort.%03d", storeName, df);
    dumpFile[df]   = AS_OVS_createBinaryOverlapFile(name, FALSE);
    dumpLength[df] = 0;
  }

  AS_OVS_writeOverlap(dumpFile[df], overlap);
  dumpLength[df]++;
}





void
buildStore(
      char *storeName, 
      char *gkpName, 
      uint64 memoryLimit, 
      uint32 nThreads, 
      uint32 doFilterOBT, 
      uint32 fileListLen, 
      char **fileList, 
      char *ovlSkipName, 
      Ovl_Skip_Type_t ovlSkipOpt) {

  if (gkpName == NULL) {
    fprintf(stderr, "overlapStore: The '-g gkpName' parameter is required.\n");
    exit(1);
  }

  //  We create the store early, allowing it to fail if it already
  //  exists, or just cannot be created.
  //
  OverlapStore    *storeFile = AS_OVS_createOverlapStore(storeName, TRUE);

  storeFile->gkp = new gkStore(gkpName, FALSE, FALSE);

  //  Decide on some sizes.  We need to decide on how many IID's to
  //  put in each bucket.  Except for running out of file descriptors
  //  (an OS limit), there isn't much of a penalty for having lots of
  //  buckets -- our BinaryOverlapFile buffers writes, and, in fact,
  //  we could open/close the file each time if things get too bad.
  //
  //  The 2x multiplier isn't really true -- MER overlaps don't need
  //  to be flipped, and so mer overlaps count the true number.
  //  Maybe.
  //
  uint32  i                   = 0;
  uint64  numOverlaps         = 0;
  uint64  overlapsPerIID      = 0;
  uint64  iidPerBucket        = 0;
  uint64  maxIID              = storeFile->gkp->gkStore_getNumFragments() + 1;

  if (fileList[i][0] != '-') {
    for (i=0; i<fileListLen; i++) {
      uint64  no = sizeOfFile(fileList[i]);
      if (no == 0)
        fprintf(stderr, "No overlaps found (or file not found) in '%s'.\n", fileList[i]);
      
      numOverlaps += 2 * no / sizeof(OVSoverlap);
    }

    assert(numOverlaps > 0);

    //  Small datasets die with the default maxIID; reset it to not die.
    if (maxIID > numOverlaps)
      maxIID = numOverlaps;

    uint64  overlapsPerBucket   = memoryLimit / sizeof(OVSoverlap);
    uint64  overlapsPerIID      = numOverlaps / maxIID;

    iidPerBucket                = overlapsPerBucket / overlapsPerIID;

    fprintf(stderr, "For %.3f million overlaps, in "F_U64"MB memory, I'll put "F_U64" IID's (approximately "F_U64" overlaps) per bucket.\n",
            numOverlaps / 1000000.0,
            memoryLimit / (uint64)1048576,
            iidPerBucket,
            overlapsPerBucket);
  } else {
    iidPerBucket = maxIID / 100;
    fprintf(stderr, "Using stdin; cannot size buckets properly.  Using 100 buckets, "F_U64" IIDs per bucket.  Good luck!\n",
            iidPerBucket);
  }

  int                      dumpFileMax = sysconf(_SC_OPEN_MAX);
  BinaryOverlapFile      **dumpFile    = (BinaryOverlapFile **)safe_calloc(sizeof(BinaryOverlapFile *), dumpFileMax);
  uint64                  *dumpLength  = (uint64 *)safe_calloc(sizeof(uint64), dumpFileMax);

  // record the UIDs we should not process in the hashtable if we are supplied a file listing them  
  HashTable_AS *readIIDsToSkip = CreateScalarHashTable_AS();
  if (doFilterOBT == 0 && ovlSkipName != NULL) {
     int line_len = ( 16 * 1024 * 1024);
     char *currLine = (char *)safe_malloc(sizeof(char)*line_len);
     char fileName[1024];
     errno = 0;
     FILE *file = fopen(ovlSkipName, "r");
     if (errno) {
        errno = 0;
     } else {
        while (fgets(currLine, line_len-1, file) != NULL) {
           AS_UID read = AS_UID_lookup(currLine, NULL);
           InsertInHashTable_AS(readIIDsToSkip, storeFile->gkp->gkStore_getUIDtoIID(read, NULL), 0, 1, 0);
         }
         fclose(file);
      }
      safe_free(currLine);
  }

  for (i=0; i<fileListLen; i++) {
    BinaryOverlapFile  *inputFile;
    OVSoverlap          fovrlap;
    OVSoverlap          rovrlap;
    int                 df;

    fprintf(stderr, "bucketizing %s\n", fileList[i]);

    inputFile = AS_OVS_openBinaryOverlapFile(fileList[i], FALSE);

    while (AS_OVS_readOverlap(inputFile, &fovrlap)) {

      //  If filtering for OBT, skip the crap.
      if ((doFilterOBT == 1) && (AS_OBT_acceptableOverlap(fovrlap) == 0))
        continue;

      //  If filtering for OBTs dedup, skip the good
      if ((doFilterOBT == 2) && (AS_OBT_acceptableOverlap(fovrlap) == 1))
        continue;

      if (doFilterOBT == 0) {
         int firstIgnore = ExistsInHashTable_AS(readIIDsToSkip, fovrlap.a_iid, 0);
         int secondIgnore = ExistsInHashTable_AS(readIIDsToSkip, fovrlap.b_iid, 0);
         
         // option means don't ignore them at all
         if (ovlSkipOpt == NONE) {
         }
         // option means don't overlap them at all
         else if (ovlSkipOpt == ALL && ((firstIgnore == TRUE || secondIgnore == TRUE))) {
            continue;
         }
         // option means let them overlap other reads but not each other
         else if (ovlSkipOpt == INTERNAL && ((firstIgnore == TRUE && secondIgnore == TRUE))) {
            continue;
         }
      }

      writeToDumpFile(&fovrlap, dumpFile, dumpFileMax, dumpLength, iidPerBucket, storeName);

      //  flip the overlap -- copy all the dat, then fix whatever
      //  needs to change for the flip.

      switch (fovrlap.dat.ovl.type) {
        case AS_OVS_TYPE_OVL:
          rovrlap.a_iid = fovrlap.b_iid;
          rovrlap.b_iid = fovrlap.a_iid;
          rovrlap.dat   = fovrlap.dat;
          if (fovrlap.dat.ovl.flipped) {
            rovrlap.dat.ovl.a_hang = fovrlap.dat.ovl.b_hang;
            rovrlap.dat.ovl.b_hang = fovrlap.dat.ovl.a_hang;
          } else {
            rovrlap.dat.ovl.a_hang = -fovrlap.dat.ovl.a_hang;
            rovrlap.dat.ovl.b_hang = -fovrlap.dat.ovl.b_hang;
          }

          writeToDumpFile(&rovrlap, dumpFile, dumpFileMax, dumpLength, iidPerBucket, storeName);
          break;
        case AS_OVS_TYPE_OBT:
          rovrlap.a_iid = fovrlap.b_iid;
          rovrlap.b_iid = fovrlap.a_iid;
          rovrlap.dat   = fovrlap.dat;
          if (fovrlap.dat.obt.fwd) {
            rovrlap.dat.obt.a_beg = fovrlap.dat.obt.b_beg;
            rovrlap.dat.obt.a_end = fovrlap.dat.obt.b_end;
            rovrlap.dat.obt.b_beg = fovrlap.dat.obt.a_beg;
            rovrlap.dat.obt.b_end = fovrlap.dat.obt.a_end;
          } else {
            rovrlap.dat.obt.a_beg = fovrlap.dat.obt.b_end;
            rovrlap.dat.obt.a_end = fovrlap.dat.obt.b_beg;
            rovrlap.dat.obt.b_beg = fovrlap.dat.obt.a_end;
            rovrlap.dat.obt.b_end = fovrlap.dat.obt.a_beg;
          }

          writeToDumpFile(&rovrlap, dumpFile, dumpFileMax, dumpLength, iidPerBucket, storeName);
          break;
        case AS_OVS_TYPE_MER:
          //  Not needed; MER outputs both overlaps
          break;
        default:
          assert(0);
          break;
      }
    }

    AS_OVS_closeBinaryOverlapFile(inputFile);

#if 0
    //  This was VERY useful for creating a VERY big store.  As each
    //  file was finished, we could (manually) delete them.  (Yes, I
    //  had an offline backup of what I deleted).
    {
      char newname[1024];
      sprintf(newname, "%s.loaded", fileList[i]);
      rename(fileList[i], newname);
    }
#endif
  }
  DeleteHashTable_AS(readIIDsToSkip);

  for (i=0; i<dumpFileMax; i++)
    AS_OVS_closeBinaryOverlapFile(dumpFile[i]);

  fprintf(stderr, "bucketizing DONE!\n");

  //
  //  Read each bucket, sort it, and dump it to the store
  //

  uint64 dumpLengthMax = 0;
  for (i=0; i<dumpFileMax; i++)
    if (dumpLengthMax < dumpLength[i])
      dumpLengthMax = dumpLength[i];

  OVSoverlap         *overlapsort = NULL;
  overlapsort = (OVSoverlap *)safe_malloc(sizeof(OVSoverlap) * dumpLengthMax);

  time_t  beginTime = time(NULL);

  for (i=0; i<dumpFileMax; i++) {
    char                name[FILENAME_MAX];
    int                 x;
    BinaryOverlapFile  *bof = NULL;

    if (dumpLength[i] == 0)
      continue;

    //  We're vastly more efficient if we skip the AS_OVS interface
    //  and just suck in the whole file directly....BUT....we can't do
    //  that because the AS_OVS interface is rearranging the data to
    //  make sure the store is cross-platform compatible.

    sprintf(name, "%s/tmp.sort.%03d", storeName, i);
    fprintf(stderr, "reading %s (%d)\n", name, time(NULL) - beginTime);

    bof = AS_OVS_openBinaryOverlapFile(name, FALSE);
    x   = 0;
    while (AS_OVS_readOverlap(bof, overlapsort + x))
      x++;
    AS_OVS_closeBinaryOverlapFile(bof);

    assert(x == dumpLength[i]);
    assert(x <= dumpLengthMax);

    //  There's no real advantage to saving this file until after we
    //  write it out.  If we crash anywhere during the build, we are
    //  forced to restart from scratch.  I'll argue that removing it
    //  early helps us to not crash from running out of disk space.
    //
    unlink(name);

    fprintf(stderr, "sorting %s (%d)\n", name, time(NULL) - beginTime);
    qsort_mt(overlapsort, dumpLength[i], sizeof(OVSoverlap), OVSoverlap_sort, nThreads, 16 * 1024 * 1024);

    fprintf(stderr, "writing %s (%d)\n", name, time(NULL) - beginTime);
    for (x=0; x<dumpLength[i]; x++)
      AS_OVS_writeOverlapToStore(storeFile, overlapsort + x);
  }

  AS_OVS_closeOverlapStore(storeFile);

  safe_free(overlapsort);

  //  And we have a store.
  //
  exit(0);
}
