
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

static char CM_ID[] = "$Id: overlapStore_build.c,v 1.3 2007-03-09 07:29:23 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "AS_OVS_overlap.h"
#include "AS_OVS_overlapFile.h"
#include "AS_OVS_overlapStore.h"

#include "overlapStore.h"


int
OVSoverlap_sort(const void *a, const void *b) {
  OVSoverlap const *A = (OVSoverlap const *)a;
  OVSoverlap const *B = (OVSoverlap const *)b;
  if (A->a_iid < B->a_iid)  return(-1);
  if (A->a_iid > B->a_iid)  return(1);
  if (A->b_iid < B->b_iid)  return(-1);
  if (A->b_iid > B->b_iid)  return(1);
  return(0);
}


#define DONT_FLIP_INPUT_OVERLAPS


static
off_t
sizeOfFile(const char *path) {
  struct stat  s;
  int          r;

  errno = 0;
  r = stat(path, &s);
  if (errno) {
    fprintf(stderr, "Failed to stat() file '%s': %s\n", path, strerror(errno));
    exit(1);
  }

  return(s.st_size);
}


void
writeToDumpFile(OVSoverlap          *overlap,
                BinaryOverlapFile  **dumpFile,
                uint32               dumpFileMax,
                uint32              *dumpLength,
                uint32               iidPerBucket,
                char                *storeName) {

  int df = overlap->a_iid / iidPerBucket;

  if (df >= dumpFileMax) {
    fprintf(stderr, "Too many bucket files.  Increase memory size, set the correct\n");
    fprintf(stderr, "number of fragment IIDs in the input, or split your overlaps\n");
    fprintf(stderr, "into multiple batches.\n");
    exit(1);
  }
  
  if (dumpFile[df] == NULL) {
    char name[FILENAME_MAX];
    sprintf(name, "%s/tmp.sort.%03d", storeName, df);
    dumpFile[df]   = AS_OVS_createBinaryOverlapFile(name, FALSE, TRUE);
    dumpLength[df] = 0;
  }

  AS_OVS_writeOverlap(dumpFile[df], overlap);
  dumpLength[df]++;
}





void
buildStore(char *storeName, uint64 memoryLimit, uint64 maxIID, uint32 fileListLen, char **fileList) {
  int   i;


  //  We create the store early, allowing it to fail if it already
  //  exists, or just cannot be created.
  //
  OverlapStore  *storeFile = AS_OVS_createOverlapStore(storeName, TRUE);


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
  uint64  numOverlaps    = 0;
  for (i=0; i<fileListLen; i++)
    numOverlaps += 2 * sizeOfFile(fileList[i]) / sizeof(OVSoverlap);

  uint64  overlapsPerBucket   = memoryLimit / sizeof(OVSoverlap);
  uint64  overlapsPerIID      = numOverlaps / maxIID;
  uint64  iidPerBucket        = overlapsPerBucket / overlapsPerIID;

  fprintf(stderr, "For %.3f million overlaps, in "F_U64"MB memory, I'll put "F_U64" IID's (approximately "F_U64" overlaps) per bucket.\n",
          numOverlaps / 1000000.0,
          memoryLimit / 1048576,
          iidPerBucket,
          overlapsPerBucket);

  int                      dumpFileMax = sysconf(_SC_OPEN_MAX);
  BinaryOverlapFile      **dumpFile    = (BinaryOverlapFile **)safe_calloc(sizeof(BinaryOverlapFile *), dumpFileMax);
  uint32                  *dumpLength  = (uint32 *)safe_calloc(sizeof(uint32), dumpFileMax);

  for (i=0; i<fileListLen; i++) {
    BinaryOverlapFile  *inputFile;
    OVSoverlap          fovrlap;
    OVSoverlap          rovrlap;
    int                 df;

    fprintf(stderr, "bucketizing %s\n", fileList[i]);

    inputFile = AS_OVS_createBinaryOverlapFile(fileList[i], FALSE, FALSE);

    while (AS_OVS_readOverlap(inputFile, &fovrlap)) {

      writeToDumpFile(&fovrlap, dumpFile, dumpFileMax, dumpLength, iidPerBucket, storeName);

      //  flip the overlap -- copy all the dat, then fix whatever
      //  needs to change for the flip.

#ifndef DONT_FLIP_INPUT_OVERLAPS

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

#endif

    }

    AS_OVS_closeBinaryOverlapFile(inputFile);
  }

  for (i=0; i<dumpFileMax; i++)
    AS_OVS_closeBinaryOverlapFile(dumpFile[i]);

  fprintf(stderr, "bucketizing DONE!\n");

  //
  //  Read each bucket, sort it, and dump it to the store
  //

  for (i=0; i<dumpFileMax; i++) {
    FILE   *mergeFile;
    char    name[FILENAME_MAX];
    int     x;

    if (dumpLength[i] == 0)
      continue;

    OVSoverlap *overlapsort = (OVSoverlap *)safe_malloc(sizeof(OVSoverlap) * dumpLength[i]);

    sprintf(name, "%s/tmp.sort.%03d", storeName, i);
    errno = 0;
    mergeFile = fopen(name, "rb");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for reading: %s\n", name, strerror(errno)), exit(1);

    //  We're vastly more efficient if we skip the AS_OVS interface
    //  and just suck in the whole file directly.

    AS_UTL_safeRead(mergeFile, overlapsort, "sortOverlap", sizeof(OVSoverlap), dumpLength[i]);
    if (errno)
      fprintf(stderr, "Failed to read "F_U32" overlaps from '%s': %s\n", dumpLength[i], name, strerror(errno)), exit(1);

    fclose(mergeFile);

    fprintf(stderr, "sorting %s\n", name);

    qsort(overlapsort, dumpLength[i], sizeof(OVSoverlap), OVSoverlap_sort);

    for (x=0; x<dumpLength[i]; x++)
      AS_OVS_writeOverlapToStore(storeFile, overlapsort + x);

    safe_free(overlapsort);
  }

  AS_OVS_closeOverlapStore(storeFile);

  //
  //  And we have a store.  Clean up all the temporary files.
  //

  for (i=0; i<dumpFileMax; i++) {
    char name[FILENAME_MAX];
    sprintf(name, "%s/tmp.sort.%03d", storeName, i);
    unlink(name);
  }

  exit(0);
}
