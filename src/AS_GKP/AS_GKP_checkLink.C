
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

static char const *rcsid = "$Id$";

#include <stdio.h>
#include <stdlib.h>

#include "AS_global.H"
#include "AS_GKP_include.H"
#include "AS_PER_gkpStore.H"

int
Check_LinkMesg(LinkMesg *lkg_mesg) {
  AS_IID                     frag1IID;
  AS_IID                     frag2IID;

  if (lkg_mesg->action == AS_ADD)
    gkpStore->inf.lkgInput++;

  if (lkg_mesg->action == AS_IGNORE)
    return 0;

  //  Check that the fragments are different
  //
  if (AS_UID_compare(lkg_mesg->frag1, lkg_mesg->frag2) == 0) {
    AS_GKP_reportError(AS_GKP_LKG_SELF_LINK, 0,
                       AS_UID_toString(lkg_mesg->frag1));
    if (lkg_mesg->action == AS_ADD)
      gkpStore->inf.lkgErrors++;
    return(1);
  }

  //  Check that it's a currently supported link type
  //
  if (lkg_mesg->type.isMatePair() == false) {
    AS_GKP_reportError(AS_GKP_LKG_UNSUPPORTED_TYPE, 0,
                       lkg_mesg->type.toLetter(),
                       AS_UID_toString(lkg_mesg->frag1),
                       AS_UID_toString(lkg_mesg->frag2), AS_UID_toString(lkg_mesg->distance));
    if (lkg_mesg->action == AS_ADD)
      gkpStore->inf.lkgErrors++;
    return(1);
  }


  //  Check that the two fragments are alive and well
  //
  frag1IID = gkpStore->gkStore_getUIDtoIID(lkg_mesg->frag1, NULL);
  if (frag1IID == 0) {
    AS_GKP_reportError(AS_GKP_LKG_FRG_DOESNT_EXIST, 0,
                       AS_UID_toString(lkg_mesg->frag1));
    if (lkg_mesg->action == AS_ADD)
      gkpStore->inf.lkgErrors++;
    return(1);
  }

  frag2IID = gkpStore->gkStore_getUIDtoIID(lkg_mesg->frag2, NULL);
  if (frag2IID == 0) {
    AS_GKP_reportError(AS_GKP_LKG_FRG_DOESNT_EXIST, 0,
                       AS_UID_toString(lkg_mesg->frag2));
    if (lkg_mesg->action == AS_ADD)
      gkpStore->inf.lkgErrors++;
    return(1);
  }

  //fprintf(stderr, "lkg_mesg: %s,%d <-> %s,%d\n",
  //        AS_UID_toString(lkg_mesg->frag1), frag1IID,
  //        AS_UID_toString(lkg_mesg->frag2), frag2IID);

  //  Now grab the reads, we'll need them soon enough anyway.
  //
  gkpStore->gkStore_getFragment(frag1IID, gkFrag1, GKFRAGMENT_INF);
  gkpStore->gkStore_getFragment(frag2IID, gkFrag2, GKFRAGMENT_INF);


  //  Make sure they're not deleted
  //
  if (gkFrag1->gkFragment_getIsDeleted()) {
    AS_GKP_reportError(AS_GKP_LKG_FRG_DELETED, 0,
                       AS_UID_toString(lkg_mesg->frag1));
    if (lkg_mesg->action == AS_ADD)
      gkpStore->inf.lkgErrors++;
    return(1);
  }

  if (gkFrag2->gkFragment_getIsDeleted()) {
    AS_GKP_reportError(AS_GKP_LKG_FRG_DELETED, 0,
                       AS_UID_toString(lkg_mesg->frag2));
    if (lkg_mesg->action == AS_ADD)
      gkpStore->inf.lkgErrors++;
    return(1);
  }

  //  Make sure they're not already mated
  //
  if (lkg_mesg->action == AS_ADD) {
    int err = 0;

    if ((gkFrag1->gkFragment_getMateIID() > 0) && (gkFrag1->gkFragment_getMateIID() != frag2IID)) {
      AS_GKP_reportError(AS_GKP_LKG_ALREADY_MATED, 0,
                         AS_UID_toString(gkFrag1->gkFragment_getReadUID()), gkFrag1->gkFragment_getReadIID(),
                         gkFrag1->gkFragment_getMateIID(),
                         AS_UID_toString(lkg_mesg->frag2), frag2IID);
      err++;
    }
    if ((gkFrag2->gkFragment_getMateIID() > 0) && (gkFrag2->gkFragment_getMateIID() != frag1IID)) {
      AS_GKP_reportError(AS_GKP_LKG_ALREADY_MATED, 0,
                         AS_UID_toString(gkFrag2->gkFragment_getReadUID()), gkFrag2->gkFragment_getReadIID(),
                         gkFrag2->gkFragment_getMateIID(),
                         AS_UID_toString(lkg_mesg->frag1), frag1IID);
      err++;
    }
    if (err) {
      gkpStore->inf.lkgErrors++;
      return(1);
    }

    //  If only one of the frags has the mate relationship, and the
    //  other has mateIID == 0, well, this should be difficult to do,
    //  and your store is hosed anyway.  Something else should check
    //  the store for sanity, not gatekeeper.


    //  And if they ARE mated, just quit because their mates are
    //  correct.
    //
    if ((gkFrag1->gkFragment_getMateIID() > 0) &&
        (gkFrag2->gkFragment_getMateIID() > 0))
      return(0);
  }


  //  Version 1 encodes the library in the mate, not the read.  We
  //  need to check that the library (from a distance record) is
  //  there, and get the library IID to set in the reads.
  //
  if (AS_UID_isDefined(lkg_mesg->distance) == TRUE) {
    gkFrag1->gkFragment_setLibraryIID(gkpStore->gkStore_getUIDtoIID(lkg_mesg->distance, NULL));
    gkFrag2->gkFragment_setLibraryIID(gkFrag1->gkFragment_getLibraryIID());

    if (gkFrag1->gkFragment_getLibraryIID() == 0) {
      AS_GKP_reportError(AS_GKP_LKG_LIB_DOESNT_EXIST, 0,
                         AS_UID_toString(lkg_mesg->distance));
      if (lkg_mesg->action == AS_ADD)
        gkpStore->inf.lkgErrors++;
      return(1);
    }

    gkFrag1->gkFragment_setOrientation(AS_READ_ORIENT_UNKNOWN);
    gkFrag2->gkFragment_setOrientation(AS_READ_ORIENT_UNKNOWN);

    if        (lkg_mesg->link_orient.isInnie()) {
      gkFrag1->gkFragment_setOrientation(AS_READ_ORIENT_INNIE);
      gkFrag2->gkFragment_setOrientation(AS_READ_ORIENT_INNIE);
    } else if (lkg_mesg->link_orient.isOuttie()) {
      gkFrag1->gkFragment_setOrientation(AS_READ_ORIENT_OUTTIE);
      gkFrag2->gkFragment_setOrientation(AS_READ_ORIENT_OUTTIE);
    } else if (lkg_mesg->link_orient.isNormal()) {
      gkFrag1->gkFragment_setOrientation(AS_READ_ORIENT_NORMAL);
      gkFrag2->gkFragment_setOrientation(AS_READ_ORIENT_NORMAL);
    } else if (lkg_mesg->link_orient.isAnti()) {
      gkFrag1->gkFragment_setOrientation(AS_READ_ORIENT_ANTINORMAL);
      gkFrag2->gkFragment_setOrientation(AS_READ_ORIENT_ANTINORMAL);
    } else {
      //  This really cannot happen, unless someone adds a new orientation to PairOrient.
      assert(0);
    }
  }

  //  Now make absolutely sure the two reads are in the same library.
  //
  if (gkFrag1->gkFragment_getLibraryIID() != gkFrag2->gkFragment_getLibraryIID()) {
    AS_GKP_reportError(AS_GKP_LKG_DIFFERENT_LIB, 0,
                       frag1IID, gkFrag1->gkFragment_getLibraryIID(),
                       frag2IID, gkFrag2->gkFragment_getLibraryIID());
    if (lkg_mesg->action == AS_ADD)
      gkpStore->inf.lkgErrors++;
    return(1);
  }

  //  Make sure that the library is expecting mated reads.
  //
  //  Case 1: Library has unknown (aka unset) orientation.  Version 2 should be setting this
  //  explicitly, and version 1 defaults to INNIE.
  //
  //  Case 2: Reads differ from library.  Version 2 sets the orientation when the fragment is added.
  //  Version 1 sets the orientation a few blocks of code above.
  //
  uint64  orient = gkpStore->gkStore_getLibrary(gkFrag1->gkFragment_getLibraryIID())->orientation;

  if (orient == AS_READ_ORIENT_UNKNOWN) {
    AS_GKP_reportError(AS_GKP_LKG_UNMATED_LIB, 0,
                       frag1IID, frag2IID, gkFrag1->gkFragment_getLibraryIID());
    if (lkg_mesg->action == AS_ADD)
      gkpStore->inf.lkgErrors++;
    return(1);
  }

  if ((orient != gkFrag1->gkFragment_getOrientation()) ||
      (orient != gkFrag2->gkFragment_getOrientation())) {
    AS_GKP_reportError(AS_GKP_LKG_DIFFERENT_ORIENT, 0,
                       frag1IID,gkFrag1->gkFragment_getOrientation(),
                       frag2IID, gkFrag1->gkFragment_getOrientation(),
                       gkFrag1->gkFragment_getLibraryIID(), orient);
    if (lkg_mesg->action == AS_ADD)
      gkpStore->inf.lkgErrors++;
    return(1);
  }

  //  And commit the change.

  if (lkg_mesg->action == AS_ADD) {
    gkFrag1->gkFragment_setMateIID(frag2IID);
    gkFrag2->gkFragment_setMateIID(frag1IID);

    gkpStore->gkStore_setFragment(gkFrag1);
    gkpStore->gkStore_setFragment(gkFrag2);

    //fprintf(stderr, "ADD link from %s,%d (%d) to %s,%d (%d)\n",
    //        AS_UID_toString(gkFrag1->gkFragment_getReadUID()), gkFrag1->gkFragment_getReadIID(), frag1IID,
    //        AS_UID_toString(gkFrag2->gkFragment_getReadUID()), gkFrag2->gkFragment_getReadIID(), frag2IID);


    gkpStore->inf.lkgLoaded++;
  } else if (lkg_mesg->action == AS_DELETE) {
    gkFrag1->gkFragment_setMateIID(0);
    gkFrag2->gkFragment_setMateIID(0);

    gkpStore->gkStore_setFragment(gkFrag1);
    gkpStore->gkStore_setFragment(gkFrag2);
  } else {
    AS_GKP_reportError(AS_GKP_LKG_UNKNOWN_ACTION, 0);
    return 1;
  }

  return 0;
}
