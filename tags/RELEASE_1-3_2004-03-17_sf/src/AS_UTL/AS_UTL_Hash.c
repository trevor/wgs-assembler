
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
static char CM_ID[] = "$Id: AS_UTL_Hash.c,v 1.1.1.1 2004-04-14 13:53:47 catmandew Exp $";

#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "AS_global.h"
#include "AS_UTL_HashCommon.h"
#include "AS_UTL_Hash.h"

//#define DEBUGHASH
//#define DEBUG_HASH

/* Return to free list */
static void FreeHashNode_AS(HashTable_AS *table, HashNode_AS *node);

/* Get from Free List */
static HashNode_AS *AllocHashNode_AS(HashTable_AS *table, 
				     void *key, uint32 keyLength, void *value);

static int InsertNodeInHashBucket(HashTable_AS *table, HashNode_AS *newnode);

static int Int64HashFn_AS(const void *pointerToInt64, int length){
  assert(length == sizeof(uint64));
  return Hash_AS((uint8 *)pointerToInt64, sizeof(uint64), 37);
}
static int Int64HashCmpFn_AS(const void *pint64a, const void *pint64b){
  int64 tmp;
  int flag;
  tmp = *(uint64 *)pint64a - *(uint64 *)pint64b;
  flag = (tmp > 0 ? 1 : (tmp < 0 ? -1 : 0 ));
  return flag;
}

int Int32HashFn_AS(const void *pointerToInt32, int length){
  assert(length == sizeof(uint32));
  return Hash_AS((uint8 *)pointerToInt32, sizeof(uint32), 37);
}
static int Int32HashCmpFn_AS(const void *pint32a, const void *pint32b){
  int32 tmp;
  int flag;
  tmp = *(int32 *)pint32a - *(int32 *)pint32b;
  flag = (tmp > 0 ? 1 : (tmp < 0 ? -1 : 0 ));
  return flag;
}

int StringHashFn_AS(const void *pointerToString, int length){
  assert(length == strlen((char *)pointerToString));
  return Hash_AS((uint8 *)pointerToString, length, 37);
}

/********************* HashNodePool *********************/

HashTable_AS *CreateHashTable_uint64_AS(int numItemsToHash){
  return CreateHashTable_AS(numItemsToHash, 
			    Int64HashFn_AS , Int64HashCmpFn_AS);
}

HashTable_AS *CreateHashTable_int32_AS(int numItemsToHash){
  return CreateHashTable_AS(numItemsToHash, 
			    Int32HashFn_AS , Int32HashCmpFn_AS);
}


int ReallocHashTable_AS(HashTable_AS *htable){
  /* Increase the size to the next power of two */
  int oldNumNodes = htable->numNodesAllocated;
  int oldNumBuckets = htable->numBuckets;

  /* Allocate and clear the new buckets */
    htable->freeList = NULL;
    htable->numBuckets = oldNumBuckets * 2;
    htable->numNodesAllocated = oldNumNodes *2;
    htable->hashmask = ((htable->hashmask + 1)<<1) - 1;
    free(htable->buckets);
    htable->buckets = (HashNode_AS **)malloc
      (htable->numBuckets * sizeof(HashNode_AS *));
    //    htable->allocated = (HashNode_AS *)realloc
    //      (htable->allocated, htable->numNodesAllocated * sizeof(HashNode_AS));
#ifdef DEBUG_HASH
    fprintf(stderr,"***Reallocating %d buckets and %d nodes, mask = 0x%x\n",
	    htable->numBuckets, htable->numNodesAllocated, htable->hashmask);
#endif
    {
      int i;
      HashNode_AS *node;
      HEAP_ITERATOR(HashNode_AS) iterator;

      InitHashNode_ASHeapIterator(htable->allocated, &iterator);

      for(i = 0; i < htable->numBuckets; i++){
	htable->buckets[i] = NULL;
      }

      while(NULL != (node = NextHashNode_ASHeapIterator(&iterator))){
	if(node->isFree || (node->key == NULL) || (node->keyLength == 0))
	  continue;
	node->next = NULL;
	InsertNodeInHashBucket(htable, node);
      }

    }
    return(HASH_SUCCESS);

}

HashTable_AS *CreateHashTable_AS(int numItemsToHash, 
				 HashFn_AS hashfn, HashCmpFn_AS cmpfn){
  HashTable_AS *table= (HashTable_AS *)malloc(sizeof(HashTable_AS));

  table->collisions = 0;
  table->numNodes = 0;
  table->numNodesAllocated = 0;
  table->compare = cmpfn;
  table->hash = hashfn;
  table->freeList = NULL;
  assert(numItemsToHash > 0);
  {
    int logsize = ceil_log2(numItemsToHash * 2);
    int i;
    table->numBuckets = hashsize(logsize);
    table->hashmask = hashmask(logsize);
//  Commenting out this diagnostic. -- KAR
    //  fprintf(stderr,"***Allocating %d buckets for %d elements, mask = 0x%x\n",
    //	    table->numBuckets, numItemsToHash, table->hashmask);

    table->buckets = (HashNode_AS **)malloc(table->numBuckets 
					    * sizeof(HashNode_AS *));
    for(i = 0; i < table->numBuckets; i++){
      table->buckets[i] = NULL;
    }
  }	
  {
    int logsize = ceil_log2(numItemsToHash);
    table->numNodesAllocated = hashsize(logsize);

    table->allocated = AllocateHashNode_ASHeap(table->numNodesAllocated);
  }
  assert(table->freeList == NULL);
  return table;
}


void ResetHashTable_AS(HashTable_AS *table){
  int i;

  table->collisions = 0;
  table->numNodes = 0;
  table->freeList = NULL;
    for(i = 0; i < table->numBuckets; i++){
      table->buckets[i] = NULL;
    }
    WipeHashNode_ASHeap(table->allocated);
}

void DeleteHashTable_AS(HashTable_AS *table){
  free(table->buckets);
  FreeHashNode_ASHeap(table->allocated);
  free(table);
}


int InsertNodeInHashBucket(HashTable_AS *table, HashNode_AS *newnode){
  int32 hashkey = (*table->hash)(newnode->key, newnode->keyLength);
  int bucket = hashkey & table->hashmask ; 
  HashNode_AS *node, *prevnode;
  int comparison;

#ifdef DEBUGHASH
  fprintf(stderr,
	  "InsertNodeInHashBucket key 0x%d keyLength %d  "
	  "value 0x%x hashkey = %d bucket = %d\n",
	  newnode->key, newnode->keyLength, newnode->value, hashkey, bucket);
#endif
  node = table->buckets[bucket];
  /* If bucket is empty, insert at head of bucket */
  if(!node){   /* Should usually be this way */
    table->buckets[bucket] = newnode;
    return HASH_SUCCESS;
  }
#ifdef DEBUG_HASH
    fprintf(stderr,"*** Collision on insert in bucket %d for key %ld\n",
	    bucket, newnode->key);
#endif
    table->collisions++;

  comparison = (*table->compare)(node->key,newnode->key);
  if(comparison == 0){
    /* element already present */
    return(HASH_FAILURE);
  }
  if(comparison < 0){   /* Insert at head of bucket */
#ifdef DEBUG_HASH
    fprintf(stderr,"*** Inserted befire key %ld\n",
	    table->buckets[bucket]->key);
#endif
    newnode->next = table->buckets[bucket];
    table->buckets[bucket] = newnode;
    return HASH_SUCCESS;
  }

  prevnode = node;
  node = node->next;
  while(node){
    comparison = (*table->compare)(node->key,newnode->key);
    if(comparison == 0){
      /* element already present */
      return(HASH_FAILURE);
    }
    if(comparison < 0){
#ifdef DEBUG_HASH
    fprintf(stderr,"*** Inserted after key %ld and before %ld\n",
	    prevnode->key, node->key);
#endif
      newnode->next = node;
      prevnode->next = newnode;
      return(HASH_SUCCESS);
    }
    prevnode = node;
    node = node->next;
  }
  /* if we got here, we're appending on the end of the list */
#ifdef DEBUG_HASH
    fprintf(stderr,"*** Inserted after key %ld\n",    prevnode->key);
#endif
  newnode->next = NULL;
  prevnode->next = newnode;
  return(HASH_SUCCESS);


}

/* Insert in HashTable */
/* Returns NULL if already in Hash Table */
int InsertInHashTable_AS(HashTable_AS *table, 
			 void *key, int keyLengthInBytes, void *value){
  
  HashNode_AS * newnode = AllocHashNode_AS(table, 
					   key, keyLengthInBytes, value);
  assert(NULL != newnode);
  
  return(InsertNodeInHashBucket(table, newnode));
}

void *LookupInHashTable_AS(HashTable_AS *table, 
			   void *key, int keyLengthInBytes){
  int32 hashkey = (*table->hash)(key, keyLengthInBytes);
  int bucket = hashkey & table->hashmask ; 
  HashNode_AS *node;

#ifdef DEBUG_HASH
    fprintf(stderr,"Lookup of key %x hashkey %x bucket %d\n", key, hashkey, bucket);
#endif
  assert(bucket >= 0 && bucket <= table->numBuckets);

  node = table->buckets[bucket];
  /* If bucket is empty, return */
  if(!node)
    return NULL;

  while(node){
    int comparison = (*table->compare)(node->key,key);
#ifdef DEBUG_HASH
    fprintf(stderr,"\tComparing with key %ld, comparison = %d\n",
	    node->key, comparison);
#endif
    if(comparison == 0){
#ifdef DEBUG_HASH
      fprintf(stderr,"\tFound it\n");
#endif
      return node->value;
    }else if(comparison < 0){
#ifdef DEBUG_HASH
      fprintf(stderr,"\tGave Up\n");
#endif
      return NULL; /* The nodes are ordered by key */
    }

    node = node->next;
  }
  return NULL;
}

int DeleteFromHashTable_AS(HashTable_AS *table, void *key, 
			   int keyLengthInBytes){


  int32 hashkey = (*table->hash)(key, keyLengthInBytes);
  int bucket = hashkey & table->hashmask ; 
  HashNode_AS *node, *previousNode;

  node = table->buckets[bucket];
  /* If bucket is empty, return */
  if(!node)
    return HASH_FAILURE;  /* Failure, didn't find it */

  previousNode = NULL;

  while(node){
    int comparison = (*table->compare)(node->key, key);
    if(comparison == 0){  /* Found it */
      if(previousNode){
	previousNode->next = node->next;
      }else{
	table->buckets[bucket] = node->next;
      }
      FreeHashNode_AS(table, node);
      return HASH_SUCCESS;
    }else if(comparison < 0){
      return HASH_FAILURE; /* The nodes are ordered by key */
    }

    previousNode = node;
    node = node->next;
  }
  return HASH_FAILURE;


}


/* Return to free list */
void FreeHashNode_AS(HashTable_AS *table, HashNode_AS *node){
  node->next = table->freeList;
  table->freeList = node;
  table->numNodes--;
  node->isFree = TRUE;
}

/* Get from Free List */
HashNode_AS *AllocHashNode_AS(HashTable_AS *table, 
			      void *key, uint32 keyLength, void *value){
  HashNode_AS *node;

  if(table->freeList){
    node = table->freeList;
    table->freeList = node->next;
    table->numNodes++;
  }else{
    if(table->numNodes++ >= table->numNodesAllocated)
      ReallocHashTable_AS(table);
    node = GetHashNode_ASHeapItem(table->allocated);
  }
  node->isFree = FALSE;
  node->next = NULL;
  node->keyLength = keyLength;
  node->key = key;
  node->value = value;
#ifdef DEBUG_HASH
  fprintf(stderr,"* Allocated node %x with key %x and length %u\n",
	  node, node->key, node->keyLength);
#endif
  return node;
}

/***********************************************************************************/


int InitializeHashTable_Iterator_AS(HashTable_AS *table, HashTable_Iterator_AS *iterator){

  iterator->table = table;
  InitHashNode_ASHeapIterator(table->allocated, &(iterator->iterator));
  return HASH_SUCCESS;
}

/***********************************************************************************
 * Function: NextHashTable_Iterator_AS
 * Description:
 *     Retrieves next element of HashTable using iterator
 *
 * Inputs:
 *     iterator      * to a HashTable_Iterator_AS
 * Input/Outputs
 *     value         * to a HashValue_AS *  (buffer is updated)
 *     key           * to a void *         (buffer is updated)
 *
 * Return Value:
 *     HASH_SUCCESS if successful (value is valid)
 *     HASH_FAILURE if failure
 ***********************************************************************************/
int NextHashTable_Iterator_AS(HashTable_Iterator_AS *iterator, 
			       void **key, 
			      void **value){
   HashNode_AS *node = NextHashNode_ASHeapIterator(&iterator->iterator);

   while(node && node->isFree)
     node = NextHashNode_ASHeapIterator(&iterator->iterator);

   if(node == NULL){
     *key = NULL;
     *value = NULL;
     return HASH_FAILURE;
   }
   *key = node->key;
   *value = node->value;
     return HASH_SUCCESS;
 }

