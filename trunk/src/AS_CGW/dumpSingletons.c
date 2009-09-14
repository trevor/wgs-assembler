
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

const char *mainid = "$Id: dumpSingletons.c,v 1.30 2009-09-14 13:28:45 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"

#include "AS_UTL_fasta.h"
#include "AS_UTL_reverseComplement.h"

#include "SYS_UIDclient.h"



AS_UID
getFragmentClear(int    iid,
                 int    reversecomplement,
                 char  *toprint) {

  static gkFragment  fs;
  unsigned int  clr_bgn, clr_end;

  ScaffoldGraph->gkpStore->gkStore_getFragment(iid, &fs, GKFRAGMENT_SEQ);

  fs.gkFragment_getClearRegion(clr_bgn, clr_end);

  strcpy(toprint, fs.gkFragment_getSequence() + clr_bgn);
  toprint[clr_end - clr_bgn] = 0;

  if (reversecomplement)
    reverseComplementSequence(toprint, clr_end - clr_bgn);

  return(fs.gkFragment_getReadUID());
}



int
main( int argc, char **argv) {
  int          ckptNum           = NULLINDEX;
  int          makeMiniScaffolds = 1;
  uint64       uidStart          = 1230000;
  UIDserver   *uids              = NULL;

  GlobalData = new Globals_CGW();

  argc = AS_configure(argc, argv);

  int err=0;
  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-p") == 0) {
      ckptNum = GlobalData->setPrefix(argv[++arg]);
    } else if (strcmp(argv[arg], "-c") == 0) {
      strcpy(GlobalData->outputPrefix, argv[++arg]);
    } else if (strcmp(argv[arg], "-g") == 0) {
      strcpy(GlobalData->gkpStoreName, argv[++arg]);
    } else if (strcmp(argv[arg], "-n") == 0) {
      ckptNum = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-U") == 0) {
      uidStart = 0;
    } else if (strcmp(argv[arg], "-S") == 0) {
      makeMiniScaffolds = 0;
    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err = 1;
    }
    arg++;
  }

  if ((GlobalData->outputPrefix[0]      == 0) ||
      (GlobalData->gkpStoreName[0] == 0)) {
    fprintf(stderr, "usage: %s [[-p prefix] | [-c name -g gkpstore -n ckptNum]] [-U] [-S]\n", argv[0]);
    fprintf(stderr, "  -p      Attempt to locate the last checkpoint in directory 7-CGW.\n");
    fprintf(stderr, "  -c      Look for checkpoints in 'name'\n");
    fprintf(stderr, "  -g      Path to gkpStore\n");
    fprintf(stderr, "  -n      Checkpoint number to load\n");
    fprintf(stderr, "  -U      Use real UIDs for miniscaffolds, otherwise, UIDs start at 1230000\n");
    fprintf(stderr, "  -S      Do NOT make mini scaffolds.\n");
    exit(1);
  }

  uids = UIDserverInitialize(256, uidStart);

  char *toprint = (char *)safe_malloc(sizeof(char) * (AS_READ_MAX_LEN + 51 + AS_READ_MAX_LEN + 2));

  LoadScaffoldGraphFromCheckpoint(GlobalData->outputPrefix, ckptNum, FALSE);

  int ifrag;
  for (ifrag=0; ifrag < GetNumVA_CIFragT(ScaffoldGraph->CIFrags); ifrag++) {
    CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags, ifrag);
    CIFragT *mate = NULL;

    assert(frag->cid != NULLINDEX);
    assert((frag->flags.bits.hasMate == 0) || (frag->mate_iid != NULLINDEX));

    //  Fix for missing mates -- OBT used to not delete mate links, leaving
    //  dangling mates.  Somebody else seems to be doing this too.
    //
    if (frag->flags.bits.hasMate) {
      mate = GetCIFragT(ScaffoldGraph->CIFrags, frag->mate_iid);
      if (mate == NULL)
        frag->flags.bits.hasMate = 0;
    }

    //  If this fragment is not chaff, we have nothing to do here.
    //
    if (GetGraphNode(ScaffoldGraph->CIGraph,frag->cid)->flags.bits.isChaff == 0)
      continue;

    //  Print a singleton if there is no mate, the mate isn't chaff,
    //  or we were told to not make miniscaffolds.
    //
    if ((mate == NULL) ||
        (mate->flags.bits.isChaff == 0) ||
        (makeMiniScaffolds == 0)) {
      AS_UID  fUID = getFragmentClear(frag->read_iid, 0, toprint);

      AS_UTL_writeFastA(stdout,
                        toprint, strlen(toprint),
                        ">%s /type=singleton\n", AS_UID_toString(fUID));

    } else if ((mate != NULL) &&
               (mate->flags.bits.isChaff == 1) &&
               (makeMiniScaffolds == 1) &&
               (frag->read_iid < mate->read_iid)) {

      //  make sure the following chain of Ns is divisible by three;
      //  the exact length is arbitrary but Doug Rusch points out that
      //  by making it divisible by 3, we can get lucky and maintain
      //  the phase of a protein ...  which helps in the
      //  auto-annotation of environmental samples

      AS_UID  fUID = getFragmentClear(frag->read_iid, 0, toprint);

      strcat(toprint, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");

      AS_UID  mUID = getFragmentClear(mate->read_iid, 1, toprint + strlen(toprint));

      AS_UTL_writeFastA(stdout,
                        toprint, strlen(toprint),
                        ">"F_U64" /type=mini_scaffold /frgs=(%s,%s)\n",
                        getUID(uids),
                        AS_UID_toString(fUID),
                        AS_UID_toString(mUID));
    }
  }

  delete GlobalData;

  exit(0);
}
