
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

/* 	$Id: AS_PER_gkpStore.h,v 1.31 2007-04-23 22:21:59 brianwalenz Exp $	 */

#ifndef AS_PER_GKPFRGSTORE_H
#define AS_PER_GKPFRGSTORE_H

#include <time.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_PER_genericStore.h"
#include "AS_UTL_PHash.h"
#include "AS_UTL_Hash.h"

#define NULL_LINK 0

#define AS_PER_NAME_LEN      256
#define AS_PER_COMMENT_LEN   256

//  The following counts represent the number of records of each type
//  PRIOR to processing this batch.  To get the range of batch i, take
//  the difference between batch i+1 and batch i.
//
typedef struct {
  CDS_UID_t      batchUID;
  char           name[AS_PER_NAME_LEN];
  char           comment[AS_PER_COMMENT_LEN];
  uint64         created;

  unsigned int   deleted:1;
  unsigned int   spare:31;

  int32          numFragments;
  int32          numLibraries;
  int32          numLibraries_s;
} GateKeeperBatchRecord;

#define AS_READ_ORIENT_UNKNOWN    0x00
#define AS_READ_ORIENT_INNIE      0x01
#define AS_READ_ORIENT_OUTTIE     0x02
#define AS_READ_ORIENT_NORMAL     0x03
#define AS_READ_ORIENT_ANTINORMAL 0x04

static const char *AS_READ_ORIENT_NAMES[8] = {
  "U", "I", "O", "N", "A", "X", "X", "X"
};

typedef struct {
  CDS_UID_t      libraryUID;

  char           name[AS_PER_NAME_LEN];
  char           comment[AS_PER_COMMENT_LEN];
  uint64         created;

  //  Features: you can add boolean flags and small-value types to the
  //  64-bit-wide bit-vector immediately below.  And just in case you
  //  need A LOT of space, you've got ONE-HUNDRED-AND-TWENTY-EIGHT
  //  bits!!  (OK, minus 5, that are currently used).
  //
  uint64         spare2:64;
  uint64         spare1:59;
  uint64         orientation:3;
  uint64         redefined:1;
  uint64         deleted:1;

  double         mean;
  double         stddev;
  
  CDS_IID_t      prevInstanceID;    // Previous definitions are linked by this reference
  CDS_IID_t      prevID;            // If redefined == TRUE, the original ID of this

  uint16         birthBatch;        // This entry is valid
  uint16         deathBatch;        // [birthBatch, deatchBatch)

  //  Features: you can add more complicated data to the structure
  //  below.  It is in a union, of size 16KB, so that lots of data can
  //  be added without breaking backwards compatibility -- as long as
  //  your features don't break if they are all zero-valued (hint, use
  //  a boolean flag to indicate the data is valid) you can add a new
  //  feature without breaking older stores.
  //
  //  The idea here is that we allow the on-disk data (the structure
  //  in the union below) to be extended to accomodate new features,
  //  and also extend the AS_MSG library message with some variable,
  //  freeform, annotations.  This should allow both stores and
  //  message format to be flexible enough to extend to new features.
  //
#define AS_PER_LIBRARY_FEATURE_DATA_SIZE   16 * 1024
  union {
    char         bytes[AS_PER_LIBRARY_FEATURE_DATA_SIZE];

    struct {
      //  Example 1; indicate that the reads in this library should go
      //  between some specific clone.
      //
      //CDS_UID_t   cloneRestrictLeft;
      //CDS_UID_t   cloneRestrictRight;

      //  Example 2: just a collection of whatever.  This would be accessed
      //  as gkpl->features.data.whatever.a
      //
      //struct {
      //  uint32      a;
      //  uint32      b;
      //  double      c;
      //  char        m[256];
      //} whatever;

    } data;
  } features;

} GateKeeperLibraryRecord;


#define AS_READ_STATUS_G   0x00
#define AS_READ_STATUS_B   0x01
#define AS_READ_STATUS_U   0x02
#define AS_READ_STATUS_W   0x03
#define AS_READ_STATUS_X   0x04
#define AS_READ_STATUS_V   0x05
#define AS_READ_STATUS_E   0x06
#define AS_READ_STATUS_I   0x07
#define AS_READ_STATUS_R   0x08

static const char *AS_READ_STATUS_NAMES[9] = {
  "G", "B", "U", "W", "X", "V", "E", "I", "R"
};


// Clients must specify which clear range to use.
//
// The default 'get' function returns the latest clear range.  Latest
// is defined as the highest number here.
//
// UNTRIM is a special case that returns untrimmed sequence (0, len);
// see getFragRecordClearRegion() for example.
//
// NUM _must_ be the number of real clear ranges we store, not the
// number of real and virtual (like UNTRIM) ranges.
//
#define AS_READ_CLEAR_ORIG     0  //  read only
#define AS_READ_CLEAR_QLT      1  //  read only
#define AS_READ_CLEAR_VEC      2  //  read only
#define AS_READ_CLEAR_OBTINI   3
#define AS_READ_CLEAR_OBT      4
#define AS_READ_CLEAR_UTG      5  //  future use
#define AS_READ_CLEAR_ECR1     6
#define AS_READ_CLEAR_ECR2     7
#define AS_READ_CLEAR_UNTRIM   8  //  read only, virtual
#define AS_READ_CLEAR_LATEST   (AS_READ_CLEAR_ECR2)

//  These are private to AS_PER.  Please don't use!
#define AS_READ_CLEAR_NUMREAL  (AS_READ_CLEAR_LATEST + 1)
#define AS_READ_CLEAR_NUMVIRT  1

static const char *AS_READ_CLEAR_NAMES[AS_READ_CLEAR_NUMREAL + AS_READ_CLEAR_NUMVIRT] = {
  "ORIG", "QLT", "VEC", "OBTINI", "OBT", "UTG", "ECR1", "ECR2", "UNTRIM"
};

static
uint32
AS_PER_decodeClearRangeLabel(const char *label) {
  uint32 clr = AS_READ_CLEAR_LATEST;

  for (clr=0; clr<AS_READ_CLEAR_NUMREAL + AS_READ_CLEAR_NUMVIRT; clr++)
    if (strcasecmp(label, AS_READ_CLEAR_NAMES[clr]) == 0)
      return(clr);

  fprintf(stderr, "AS_PER_decodeClearRangeLabel()-- unknown clear range label '%s'\n", label);
  fprintf(stderr, "AS_PER_decodeClearRangeLabel()-- valid labels are:\n");
  fprintf(stderr, "AS_PER_decodeClearRangeLabel()--");
  for (clr=0; clr<AS_READ_CLEAR_NUMREAL + AS_READ_CLEAR_NUMVIRT; clr++)
    fprintf(stderr, " %s", AS_READ_CLEAR_NAMES[clr]);
  fprintf(stderr, "\n");

  exit(1);
  return(AS_READ_CLEAR_LATEST);
}


typedef struct{
  CDS_UID_t        readUID;
 
  CDS_IID_t        readIID;
  CDS_IID_t        mateIID;

  CDS_UID_t        plateUID;

  CDS_IID_t        libraryIID;
  uint32           deleted:1;
  uint32           nonrandom:1;
  uint32           status:4;
  uint32           orientation:3;  //  copied from the library
  uint32           plateLocation:8;
  uint32           pad1:15;

  uint64           seqLen:12;
  uint64           hpsLen:12;
  uint64           srcLen:12;
  uint64           birthBatch:10;  //  This entry is valid
  uint64           deathBatch:10;  //  [birthBatch, deatchBatch)
  uint64           pad2:8;

  VLSTRING_SIZE_T  clearBeg[AS_READ_CLEAR_NUMREAL];
  VLSTRING_SIZE_T  clearEnd[AS_READ_CLEAR_NUMREAL];

  uint64           seqOffset;
  uint64           qltOffset;
  uint64           hpsOffset;
  uint64           srcOffset;
} GateKeeperFragmentRecord;





////////////////////////////////////////////////////////////

//  The fragRecord is usually how one should access fragments in the
//  gatekeeper store.  It gets you the fragment info above, and
//  sequence, quality, hps data and source string.

#define MAX_SEQ_LENGTH (AS_READ_MAX_LEN + 1)
#define MAX_HPS_LENGTH (AS_READ_MAX_LEN + 1)
#define MAX_SRC_LENGTH (512 + sizeof(int32) + sizeof(int64))

#define FRAG_S_INF 0x00
#define FRAG_S_SEQ 0x01
#define FRAG_S_QLT 0x02
#define FRAG_S_HPS 0x04
#define FRAG_S_SRC 0x08
#define FRAG_S_ALL 0x0f

typedef struct {
  GateKeeperFragmentRecord   gkfr;
  uint32                     hasSEQ:1;
  uint32                     hasQLT:1;
  uint32                     hasHPS:1;
  uint32                     hasSRC:1;
  char                       seq[MAX_SEQ_LENGTH];
  char                       qlt[MAX_SEQ_LENGTH];
  char                       hps[MAX_HPS_LENGTH];
  char                       src[MAX_SRC_LENGTH];
} fragRecord;


fragRecord *new_fragRecord(void);
void        del_fragRecord(fragRecord *fr);
void        clr_fragRecord(fragRecord *fr);

void        setFragRecordClearRegion(fragRecord *fr, 
                                     uint32 start,
                                     uint32 end,
                                     uint32 which);

static
CDS_UID_t   getFragRecordUID(fragRecord *fr) {
  return(fr->gkfr.readUID);
};

static
CDS_IID_t   getFragRecordIID(fragRecord *fr) {
  return(fr->gkfr.readIID);
};

static
CDS_IID_t   getFragRecordMateIID(fragRecord *fr) {
  return(fr->gkfr.mateIID);
};

static
CDS_IID_t   getFragRecordLibraryIID(fragRecord *fr) {
  return(fr->gkfr.libraryIID);
};

static
int         getFragRecordIsDeleted(fragRecord *fr) {
  return(fr->gkfr.deleted);
};

static
int         getFragRecordIsNonRandom(fragRecord *fr) {
  return(fr->gkfr.nonrandom);
};

void        getFragRecordClearRegion(fragRecord *fr, uint32 *start, uint32 *end, uint32 flags);
uint32      getFragRecordClearRegionBegin(fragRecord *fr, uint32 flags);
uint32      getFragRecordClearRegionEnd  (fragRecord *fr, uint32 flags);


static
int         getFragRecordSequenceLength(fragRecord *fr) {
  return(fr->gkfr.seqLen);
}
static
int         getFragRecordQualityLength(fragRecord *fr) {
  return(fr->gkfr.seqLen);
}
static
int         getFragRecordHPSLength(fragRecord *fr) {
  return(fr->gkfr.hpsLen);
}
static
int         getFragRecordSourceLength(fragRecord *fr) {
  return(fr->gkfr.srcLen);
}


static
char       *getFragRecordSequence(fragRecord *fr) {
  assert(fr->hasSEQ);
  return(fr->seq);
}
static
char       *getFragRecordQuality(fragRecord *fr) {
  assert(fr->hasQLT);
  return(fr->qlt);
}
static
char       *getFragRecordHPS(fragRecord *fr) {
  assert(fr->hasHPS);
  return(fr->hps);
}
static
char       *getFragRecordSource(fragRecord *fr) {
  assert(fr->hasSRC);
  return(fr->src);
}

////////////////////////////////////////////////////////////

#define INDEXSTORE_DEF_EXTEND(type)\
static int deleteAndMark ## type ## Store(type ## Store fs, int index, int batchID){\
  type ## Record dr;\
  getIndexStore(fs,index,&dr); \
  dr.deleted = TRUE;\
  dr.deathBatch = batchID;\
  setIndexStore(fs,index,&dr);\
  return(0);\
}

#define INDEXSTORE_DEF(type)\
typedef StoreHandle type ## Store;\
static int commit ## type ## Store(type ## Store sh){\
  return commitStore(sh);\
}\
static type ## Store reset ## type ## Store(type ## Store sh, int firstID){\
  return resetIndexStore(sh, firstID);\
}\
static int close ## type ## Store(type ## Store sh){\
  return closeStore(sh);\
}\
static int delete ## type ## Store(type ## Store fs, int index){\
  type ## Record dr;\
  getIndexStore(fs,index,&dr); \
  dr.deleted = TRUE;\
  setIndexStore(fs,index,&dr);\
  return(0);\
}\
static int get ## type ## Store(type ## Store fs, int index, type ## Record *dr){\
  return getIndexStore(fs,index,dr); \
}\
static int set ## type ## Store(type ## Store fs, int index, type ## Record *dr){\
  return setIndexStore(fs,index,dr); \
}\
static type ## Store create ## type ## Store(char *StorePath, char *ext, int firstID){\
  type ## Store s = createIndexStore(StorePath,ext, sizeof(type ## Record), 1, firstID);\
  return s;\
}\
static type ## Store open ## type ## Store(char *StorePath, char *rw){\
  return openStore(StorePath, rw);\
}\
static int append ## type ## Store(type ## Store store, type ## Record *element){\
  return appendIndexStore(store,element);\
}\
static int32 getNum ## type ## s(type ## Store store){\
  StoreStat stat;\
  statsStore(store, &stat);\
  return(stat.lastElem);\
}

INDEXSTORE_DEF(GateKeeperBatch);
INDEXSTORE_DEF(GateKeeperFragment);
INDEXSTORE_DEF_EXTEND(GateKeeperFragment);
INDEXSTORE_DEF(GateKeeperLibrary);
INDEXSTORE_DEF_EXTEND(GateKeeperLibrary);

#define NUM_GKP_FILES 10

// 1  is gatekeeper store info
// 2  is batches
// 3  is fragments
// 4  is libraries
// 5  is shadow libraries
// 6  is sequence
// 7  is quality
// 8  is homopolymer spacing and etc
// 9  is source
// 10 is gkp.phash

typedef struct {
  uint64    gkpMagic;
  uint64    gkpVersion;
  uint32    gkpBatchRecordSize;
  uint32    gkpLibraryRecordSize;
  uint32    gkpFragmentRecordSize;
} GateKeeperStoreInfo;

#define UID_NAMESPACE_AS 'U'

typedef struct {
  char                     storePath[FILENAME_MAX];

  GateKeeperStoreInfo      gkp;

  GateKeeperBatchStore     bat;
  GateKeeperFragmentStore  frg;
  GateKeeperLibraryStore   lib;    
  GateKeeperLibraryStore   lis;

  StoreHandle              seq;
  StoreHandle              qlt;
  StoreHandle              hps;

  StoreHandle              src;

  //  This is now a _private_ member.  It is not allocated when the
  //  gatekeeper is initially loaded -- because on big assemblies, it
  //  blows the address space of grid-based jobs like overlap,
  //  consensus, etc.
  //
  PHashTable_AS           *phs_private;

  //  The rest are for a partitioned fragment store.
  //
  //  Notice that we do not load the 'seq' store; we load the 'qlt'
  //  store which contains both the sequence and quality values.
  //
  //  We load all frg and qlt in this partition into memory, and
  //  optionally load hps and src.  The map converts an iid (global
  //  fragment iid) into a pointer to the correct frg record, which we
  //  can then use to grab the sequence/quality/hps/src.  (whatever
  //  builds the partitions needs to reset the offsets in the
  //  partitioned frg to be the correct offsets in these stores).
  //
  int32                    partnum;

  StoreHandle              partfrg;
  StoreHandle              partqlt;
  StoreHandle              parthps;
  StoreHandle              partsrc;

  HashTable_AS            *partmap;
} GateKeeperStore;



//  The only public accessor for the persistent hash in the
//  gatekeeper.
//
static
int
getGatekeeperUIDtoIID(GateKeeperStore *gkp, CDS_UID_t key, PHashValue_AS *value) {

  if (gkp->phs_private == NULL) {
    char name[FILENAME_MAX];

    sprintf(name,"%s/phs", gkp->storePath);
    gkp->phs_private = OpenReadOnlyPHashTable_AS(name);
    if (gkp->phs_private == NULL) {
      fprintf(stderr,"**** Failed to open GateKeeper Persistent HashTable...\n");
      assert(0);
    }

    fprintf(stderr, "Loaded GateKeeperStore Persistent Hash Table.\n");
  }

  return(LookupInPHashTable_AS(gkp->phs_private,
                               UID_NAMESPACE_AS,
                               key,
                               value));
}


////////////////////////////////////////////////////////////////////////////////

GateKeeperStore *createGateKeeperStore(const char *path);
int              testOpenGateKeeperStore(const char *path, int writable);
GateKeeperStore *openGateKeeperStore(const char *path, int writable);
void             closeGateKeeperStore(GateKeeperStore *gkpStore);

GateKeeperStore *createGateKeeperPartition(const char *path, uint32 partnum);

void             loadGateKeeperPartition(GateKeeperStore *gkp, uint32 partnum);

void             loadGateKeeperStorePartial(GateKeeperStore *gkpStore,
                                            int64 firstElem,
                                            int64 lastElem,
                                            int flags);

void    clearGateKeeperBatchRecord(GateKeeperBatchRecord *g);
void    clearGateKeeperLibraryRecord(GateKeeperLibraryRecord *g);
void    clearGateKeeperFragmentRecord(GateKeeperFragmentRecord *g);

void    getFrag(GateKeeperStore *gkp, int64 iid, fragRecord *fr, int32 flags);
void    setFrag(GateKeeperStore *gkp, int64 iid, fragRecord *fr);
void    delFrag(GateKeeperStore *gkp, int64 iid);

#define getFirstElemFragStore(GKP)  getFirstElemStore((GKP)->frg)
#define getLastElemFragStore(GKP)   getLastElemStore((GKP)->frg)


////////////////////////////////////////////////////////////////////////////////

typedef struct {
  GateKeeperStore   *gkp;
  StreamHandle       frg;
  StreamHandle       seq;
  StreamHandle       qlt;
  StreamHandle       hps;
  StreamHandle       src;

  char              *frgBuffer;
  char              *seqBuffer;
  char              *qltBuffer;
  char              *hpsBuffer;
  char              *srcBuffer;

  int                flags;
} FragStream;

FragStream      *openFragStream(GateKeeperStore *gkp, int flags);
void             resetFragStream(FragStream *fs, int64 startIndex, int64 endIndex);
void             closeFragStream(FragStream *fs);

int64            getStartIndexFragStream(FragStream *fs);

int              nextFragStream(FragStream *fs, fragRecord *fr);




#endif
