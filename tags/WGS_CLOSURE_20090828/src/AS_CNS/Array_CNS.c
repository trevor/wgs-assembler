
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

static const char *rcsid = "$Id: Array_CNS.c,v 1.22 2009-06-10 18:05:13 brianwalenz Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_PER_gkpStore.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_reverseComplement.h"
#include "UtilsREZ.h"
#include "MultiAlignStore_CNS.h"
#include "MultiAlignment_CNS.h"
#include "Array_CNS.h"

typedef struct LaneNode {
  IntMultiPos *read;
  int read_length;
  char *sequence;
  char *quality;
  struct LaneNode *prev;
  struct LaneNode *next;
} LaneNode;

typedef struct Lane {
  struct LaneNode *first;
  struct LaneNode *last;
  int lastcol;
} Lane;

VA_DEF(Lane)



static
LaneNode *
createLaneNode(IntMultiPos *read) {
  LaneNode *t = (LaneNode *)safe_malloc(sizeof(LaneNode));
  t->read = read;
  t->read_length = 0;
  t->sequence=NULL;
  t->quality=NULL;
  t->prev=NULL;
  t->next=NULL;
  return t;
}

static
int
 freeLaneNode(LaneNode *node) {
  safe_free(node->sequence);
  safe_free(node->quality);
  safe_free(node);
  return 1;
}

static
int
 PushLaneNode(LaneNode *new_lane_node, Lane *lane) {
  int leftpos = (new_lane_node->read->position.bgn<new_lane_node->read->position.end) ?
    new_lane_node->read->position.bgn : new_lane_node->read->position.end;
  if (leftpos < lane->lastcol+3) return 0;
  if ( lane->last == NULL ) {
    lane->first = new_lane_node;
    lane->last = new_lane_node;
    new_lane_node->next=NULL;
  } else {
    new_lane_node->prev = lane->last;
    lane->last = new_lane_node;
    if (new_lane_node->prev != NULL) new_lane_node->prev->next = new_lane_node;
  }
  lane->lastcol = leftpos +
    new_lane_node->read->delta_length +
    new_lane_node->read_length;
  return 1;
}

static
void
ClearLane(Lane *lane) {
  lane->first=NULL;
  lane->last=NULL;
  lane->lastcol = -3;
}

static
void
 FreeLane(Lane *lane) {
  if (lane->first) {
    // free nodes in lane
    LaneNode *node=lane->first;
    LaneNode *next;
    while (node) {
      next = node->next;
      freeLaneNode(node);
      node = next;
    }
  }
}

static
int
 IntMultiPositionCmp( const IntMultiPos *l, const IntMultiPos *m) {
  int ltmp,mtmp;
  ltmp = (l->position.bgn<l->position.end)?l->position.bgn:l->position.end;
  mtmp = (m->position.bgn<m->position.end)?m->position.bgn:m->position.end;
  if (ltmp == mtmp) return 0;
  return (ltmp > mtmp ) ? 1 : -1;
}



int
 IMP2Array(IntMultiPos *all_frags,
              int num_pieces,
              int length,
              gkStore *frag_store,
              int *depth,
              char ***array,
              int ***id_array,
              int ***ori_array,
              int show_cel_status,
              uint32 clrrng_flag) {
  char **multia = NULL;
  int **ia = NULL;
  int **oa = NULL;
  int i,lane_depth;
  uint clr_bgn, clr_end;
  FragMesg frag;
  Lane null_lane;
  int next_lane; Lane *free_lane,*lane; LaneNode *new_mlp; Lane space;
  int rc;
  char seq[AS_READ_MAX_LEN+1];
  char qv[AS_READ_MAX_LEN+1];
  gkFragment fsread;

  VA_TYPE(Lane) *Packed;
  lane_depth = 32;
  Packed = (VA_TYPE(Lane) *) CreateVA_Lane(lane_depth);
  frag.action = AS_ADD;
  frag.source = NULL;
  null_lane.first=NULL;
  null_lane.last=NULL;
  null_lane.lastcol=-3;
  for (i=0;i<lane_depth;i++){
    SetLane(Packed,i,&null_lane);
  }
  // Sort the fragments by leftmost position within contig
  qsort(all_frags, num_pieces, sizeof(IntMultiPos),
        (int (*)(const void *,const void *))IntMultiPositionCmp);
  next_lane=0;
  for (i=0;i<num_pieces;i++) {
    new_mlp = createLaneNode(&all_frags[i]);

    frag_store->gkStore_getFragment(all_frags[i].ident,&fsread,GKFRAGMENT_QLT);
    fsread.gkFragment_getClearRegion(clr_bgn, clr_end, clrrng_flag);
    new_mlp->read_length = fsread.gkFragment_getSequenceLength();

    strcpy(seq, fsread.gkFragment_getSequence());
    strcpy(qv,  fsread.gkFragment_getQuality());

    frag.eaccession = fsread.gkFragment_getReadUID();

    // All this frag stuff is defined in case it becomes important to print the frag
    frag.iaccession = all_frags[i].ident;
    frag.sequence = seq;
    frag.quality = qv;
    frag.type = all_frags[i].type;
    frag.clear_rng.bgn = clr_bgn;
    frag.clear_rng.end = clr_end;

    new_mlp->read_length = clr_end-clr_bgn;
    new_mlp->sequence = (char *)safe_malloc(sizeof(char)*(new_mlp->read_length+1));
    new_mlp->quality = (char *)safe_malloc(sizeof(char)*(new_mlp->read_length+1));

    seq[clr_end] = '\0';
    qv[clr_end] = '\0';

    strcpy(new_mlp->sequence,&seq[clr_bgn]);
    strcpy(new_mlp->quality,&qv[clr_bgn]);
    if (new_mlp->read->position.bgn > new_mlp->read->position.end) {
      reverseComplement(new_mlp->sequence, new_mlp->quality, strlen(new_mlp->sequence));
    }
    for (next_lane=0;next_lane<lane_depth;next_lane++) {
      free_lane = GetLane(Packed,next_lane);
      if (PushLaneNode(new_mlp,free_lane)) break;
    }
    if (next_lane==lane_depth) {  // an additional lane is needed
      ClearLane(&space);
      PushLaneNode(new_mlp,&space);
      SetLane(Packed, next_lane, &space);
      lane_depth++;
    }
  }
  {
    IntMultiPos *read;
    int col,cols;
    char *srow,*qrow;
    char laneformat[40];

    lane_depth = GetNumLanes(Packed);
    *depth =-1; // initializing is a good idea, obviously
    for (i=0;i<lane_depth;i++) {
      lane = GetLane(Packed,i);
      if (lane->first == NULL) {
        *depth = i;
        break;
      }
    }
    // if all lanes are occupied, then depth will still be -1 here,
    //  and should be set to lane_depth
    if ( *depth == -1 ) *depth = lane_depth;
    if ( *depth <= lane_depth ) {
      rc = 1;

      multia = (char **)safe_malloc(2*(*depth)*sizeof(char *));

      ia = (int **)safe_malloc((*depth)*sizeof(int *));
      oa = (int **)safe_malloc((*depth)*sizeof(int *));

      sprintf(laneformat,"%%%ds",length);
      {int j;
        for (i=0;i<(*depth);i++) {
          ia[i] = (int *) safe_malloc( length*sizeof(int));
          oa[i] = (int *) safe_malloc( length*sizeof(int));
          for (j=0;j<length;j++) {
            ia[i][j] = 0;
            oa[i][j] = 0;
          }
        }
      }
      for (i=0;i<2*(*depth);i++) {
        multia[i] = (char *) safe_malloc((length+1)*sizeof(char));
        sprintf(multia[i],laneformat," ");
        *(multia[i]+length) = '\0';
      }
      for (i=0;i<(*depth);i++) {
        int j,lastcol,firstcol,seglen;
        srow = multia[2*i];
        qrow = multia[2*i+1];
        lane = GetLane(Packed,i);
        lastcol = 0;
        if (lane->first == NULL) {
          *depth = i;
          break;
        }
        for (new_mlp=lane->first; new_mlp!=NULL; new_mlp=new_mlp->next) {
          read = new_mlp->read;
          firstcol = read->position.bgn;
          if (firstcol>read->position.end) firstcol = read->position.end;
          col = firstcol;
          cols = 0;
          for (j=0;j<read->delta_length;j++) {
            seglen = read->delta[j] - ( (j>0)?read->delta[j-1]:0 ) ;
            memcpy(srow+col,&new_mlp->sequence[cols],seglen);
            memcpy(qrow+col,&new_mlp->quality[cols],seglen);
            col+=seglen;
            srow[col] = '-';
            qrow[col] = '-';
            col++;
            cols+=seglen;
          }
          memcpy(srow+col,&new_mlp->sequence[cols],new_mlp->read_length-cols);
          memcpy(qrow+col,&new_mlp->quality[cols],new_mlp->read_length-cols);
        }
        // now, set the ids
        for (new_mlp=lane->first; new_mlp!=NULL; new_mlp=new_mlp->next) {
          int lastcol;
          int orient=0;
          read = new_mlp->read;
          firstcol = read->position.bgn;
          if (read->position.bgn>read->position.end) {
            orient = -1;
            firstcol = read->position.end;
            lastcol = read->position.bgn;
          } else {
            orient = +1;
            firstcol = read->position.bgn;
            lastcol = read->position.end;
          }
          for (col=firstcol;col<lastcol;col++) {
            if ( show_cel_status ) {
              ia[i][col] = AS_FA_READ(read->type) ? 1:0;
            } else {
              ia[i][col]=read->ident;
            }
            oa[i][col] = orient;
          }
        }
      }
    } else {
      rc = 0;
    }
  }
  { Lane *this_lane;
    for (i=0;i<lane_depth;i++){
      this_lane = GetLane(Packed,i);
      FreeLane(this_lane);
    }
    Delete_VA(Packed);
  }
  *array = multia;
  *id_array = ia;
  *ori_array = oa;
  return rc;
}
