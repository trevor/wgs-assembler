
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
/* 	$Id: AS_PER_gkpStore.h,v 1.14 2007-02-09 21:17:40 brianwalenz Exp $	 */
#ifndef AS_PER_GKPFRGSTORE_H
#define AS_PER_GKPFRGSTORE_H
/*************************************************************************
 Module:  AS_PER_gkpfrgStore
 Description:
    A thin layer on top of the IndexStore supporing the storage and
 retrieval of records used by the gatekeeper records.
    The idea is to provide easier to use shortcuts for the common
 operations, and let the other operations be accessed through the
 generic Index Store API.

 Assumptions:
    Nothing special beyond genericStore.rtf

 Document:
      GenericStore.rtf

 *************************************************************************/

#include <sys/types.h>
#include <time.h>

#include "AS_MSG_pmesg.h"
#include "AS_PER_genericStore.h"
#include "AS_UTL_PHash.h"

#define NULL_LINK 0

// The following counts represent the number of records of each type
// PRIOR to processing this batch.  To get the range of batch i, take
// the difference between batch i+1 and batch i.
//
typedef struct{
  CDS_UID_t  UID;
  char           name[256];
  time_t         created;
#ifndef __x86_64__ // 8 byte time_t on x86_64, so pad elsewhere
  uint32         padtime_t;
#endif
  char           comment[256];
  unsigned int   deleted:1;
  unsigned int   spare:31;
  int32          numFragments;
  int32          numDistances;
  int32          num_s_Distances; // shadowed for redefintions
  int32          numLinks;
}GateKeeperBatchRecord;

typedef struct{
  unsigned int   deleted:1;
  unsigned int   type:8;  /* From AS_MSG_pmesg.h FragType */
  unsigned int   numLinks:8; /* Number of LIVE links */
  uint           spare:15;
  CDS_IID_t      linkHead;           /* Index into Link Table */
  CDS_UID_t      readUID;            /* Accession ID of this read */
  uint16         birthBatch;         /* This entry is valid */
  uint16         deathBatch;         /* [birthBatch, deatchBatch) */
}GateKeeperFragmentRecord;

// One for each distance record
typedef struct{
  CDS_UID_t      UID;
  unsigned int   deleted:1;
  unsigned int   redefined:1;
  unsigned int   spare:30;
  CDS_IID_t      prevInstanceID; // Previous definitions are linked by this reference
  CDS_IID_t      prevID;         // If redefined == TRUE, the original ID of this
  float32        mean;
  float32        stddev;
  uint16         birthBatch;         /* This entry is valid */
  uint16         deathBatch;         /* [birthBatch, deatchBatch) */
}GateKeeperDistanceRecord;

#define AS_GKP_UNKNOWN 0
#define AS_GKP_INNIE 1
#define AS_GKP_OUTTIE 2
#define AS_GKP_NORMAL 3
#define AS_GKP_ANTINORMAL 4

typedef struct{
  unsigned int   deleted:1;
  unsigned int   type:8;  
  unsigned int   orientation:3;
  unsigned int   spare:23;
  CDS_IID_t      distance;   // iid of distance

  CDS_IID_t      frag1;      // iid of frag1
  CDS_IID_t      frag2;      // iid of frag2

  CDS_IID_t      frag1Next;
  CDS_IID_t      frag2Next;

  uint16         birthBatch;         /* This entry is valid */
  uint16         deathBatch;         /* [birthBatch, deatchBatch) */
  uint32         padTo8byteWords;
}GateKeeperLinkRecord;

static char getLinkOrientation(GateKeeperLinkRecord *gkpl){
  switch(gkpl->orientation){
  case AS_GKP_UNKNOWN:
    return '?';
  case AS_GKP_INNIE:
    return 'I';
  case AS_GKP_OUTTIE:
    return 'O';
  case AS_GKP_NORMAL:
    return 'N';
  case AS_GKP_ANTINORMAL:
    return 'A';
  default:
    return '-';
  }
  return '-';
}

// Stores


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



#define NUM_GKP_FILES 6

// 1. for gkp.bat
INDEXSTORE_DEF(GateKeeperBatch)

// 2. for gkp.frg
INDEXSTORE_DEF(GateKeeperFragment)
INDEXSTORE_DEF_EXTEND(GateKeeperFragment)

// 3. for gkp.lnk
INDEXSTORE_DEF(GateKeeperLink)
INDEXSTORE_DEF_EXTEND(GateKeeperLink)

// 4 & 5. for gkp.dst (distance definitions) & gkp.s_dst (distance redefinitions)
INDEXSTORE_DEF(GateKeeperDistance)
INDEXSTORE_DEF_EXTEND(GateKeeperDistance)

// 6. is gkp.phash


/***********************************************************************************
 * GateKeeperLinkIterator
 * Description:
 *     An iterator for GateKeepLinkRecords.
 *     Starting from a given entry in the gateKeeperLink store, the iterator
 *     supports traversing all of the gatekeeperLinks for a given fragment (followFrag).
 *     It maintains two points, one to the next link and one to the current link.
 *     This is useful for doing deletions from the singly-linked list.
 ***********************************************************************************/

typedef struct {
  GateKeeperLinkStore store;
  CDS_IID_t prevLinkRecord;      /* actually, the current record */
  CDS_IID_t linkRecord;          /* the next record */
  CDS_IID_t followFrag;          /* The fragment we're interested in */
}GateKeeperLinkRecordIterator;

/***********************************************************************************
 * Function: CreateGateKeeperLinkRecordIterator
 * Description:
 *     Create an iterator for GateKeeperLinkRecords
 *
 * Inputs:
 *     store      Handle of GateKeeperLink store
 *     startFromLink  The index of the link in the store to start from
 *     followFrag     The id of the fragment whose links we're traversing
 * I/O
 *     iterator      * GateKeeperLinkIterator that we're initializing    
 *
 * Return Value:
 *     Zero if success.
 ***********************************************************************************/
int CreateGateKeeperLinkRecordIterator(GateKeeperLinkStore store, CDS_IID_t startFromLink, 
				       CDS_IID_t followFrag, GateKeeperLinkRecordIterator *iterator);


/***********************************************************************************
 * Function: CreateGateKeeperLinkRecordFromFragmentIterator
 * Description:
 *     Create an iterator for GateKeeperLinkRecords
 *
 * Inputs:
 *     store      Handle of GateKeeperLink store
 *     followFrag     The id of the fragment whose links we're traversing
 * I/O
 *     iterator      * GateKeeperLinkIterator that we're initializing    
 *
 * Return Value:
 *     Zero if success.
 ***********************************************************************************/
int CreateGateKeeperLinkRecordFromFragmentIterator(GateKeeperLinkStore store,  CDS_IID_t followFrag, 
						   GateKeeperLinkRecordIterator *iterator);

/***********************************************************************************
 * Function: NextGateKeeperLinkRecordIterator
 * Description:
 *     Get the next GateKeeperLinkRecord from the store
 *
 * Inputs:
 *     iterator      * GateKeeperLinkIterator that we're initializing    
 * I/O
 *     link          * GateKeeperLinkRecord that we're retrieving
 *
 * Return Value:
 *     Non-Zero if success. Zero if we're done.
 ***********************************************************************************/
int NextGateKeeperLinkRecordIterator(GateKeeperLinkRecordIterator *iterator, GateKeeperLinkRecord *link);



/***********************************************************************************
 * Function: findLink
 * Description:
 *     Searches the GateKeeperLinkStore gkplStore from the record with index linkHead
 *     for links that match link.  If found, returns the index of the matching
 *     link.
 ***********************************************************************************/
int findLink(GateKeeperLinkStore store, 
	          CDS_IID_t frag,
		  CDS_IID_t linkHead, 
		  GateKeeperLinkRecord *searchlink,
		  GateKeeperLinkRecord *foundlink);

/***********************************************************************************
 * Function: unlinkLink_GKP
 * Description:
 *     Marks the link with index deleteLinkindex as deleted, and unlinks it
 *     it from other links.  This may involve modifying the gkpStore as well
 *     as the gkplStore.
 ***********************************************************************************/
int unlinkLink_GKP(GateKeeperLinkStore gkplStore, 
		 GateKeeperFragmentStore     gkpStore, 
		 CDS_IID_t frag1,
		 CDS_IID_t frag2,
 	         GateKeeperFragmentRecord *gkf1, 
		 GateKeeperFragmentRecord *gkf2,
		 GateKeeperLinkRecord *newLink,
		 int deleteLinkIndex);


/***********************************************************************************
 * Function: linkLink_GKP
 * Description:
 *     Inserts newLink at an appropriate place in the gkplStore, possibly modifying
 *     gkpStore in the process.
 ***********************************************************************************/
int linkLink_GKP(GateKeeperLinkStore gkplStore, 
		 GateKeeperFragmentStore     gkpStore, 
		 GateKeeperLinkRecord *newLink,
		 CDS_IID_t frag1,
		 CDS_IID_t frag2,
 	         GateKeeperFragmentRecord *gkf1, 
		 GateKeeperFragmentRecord *gkf2);




/* GateKeeperStore */
typedef struct {
  char storePath[FILENAME_MAX];

  PHashTable_AS           *hashTable;
  GateKeeperBatchStore     batStore;
  GateKeeperFragmentStore  frgStore;
  GateKeeperLinkStore      lnkStore;
  GateKeeperDistanceStore  dstStore;    
  GateKeeperDistanceStore  s_dstStore;    // Store for Distances that have been redefined
} GateKeeperStore;

int  CreateGateKeeperStore(GateKeeperStore *gkpStore);
int  OpenGateKeeperStore(GateKeeperStore *gkpStore);
int  OpenReadOnlyGateKeeperStore(GateKeeperStore *gkpStore);
int  CopyGateKeeperStoreFiles(GateKeeperStore *gkpStore, char *path);
int  RemoveGateKeeperStoreFiles(GateKeeperStore *gkpStore);
int  TestOpenGateKeeperStore(GateKeeperStore *gkpStore);
int  TestOpenReadOnlyGateKeeperStore(GateKeeperStore *gkpStore);
void InitGateKeeperStore(GateKeeperStore *gkpStore, const char *path);
void CloseGateKeeperStore(GateKeeperStore *gkpStore);
int  UpgradeGateKeeperStore(GateKeeperStore *gkpStore);

#endif