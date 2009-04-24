
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, The Venter Institute. All rights reserved.
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

static const char *rcsid = "$Id: AS_BOG_UnitigGraph.cc,v 1.111 2009-04-24 14:13:21 skoren Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_UnitigGraph.hh"
#include "AS_BOG_BestOverlapGraph.hh"

#undef max

#define VERBOSEBUILD
#define VERBOSEBREAK


void
UnitigGraph::checkUnitigMembership(void) {
  int nutg = 0;
  int nfrg = 0;

  fprintf(stderr, "checkUnitigMembership()--  numfrags=%d\n", _fi->numFragments());

  iuid   *inUnitig = new iuid [_fi->numFragments()+1];

  for (int i=0; i<_fi->numFragments()+1; i++)
    inUnitig[i] = 987654321;

  for (int  ti=0; ti<unitigs->size(); ti++) {
    Unitig  *utg = (*unitigs)[ti];

    if (utg) {
      nutg++;
      for (DoveTailIter it=utg->dovetail_path_ptr->begin(); it != utg->dovetail_path_ptr->end(); it++) {
        if (it->ident > _fi->numFragments())
          fprintf(stderr, "HUH?  ident=%d numfrags=%d\n", it->ident, _fi->numFragments());
        inUnitig[it->ident] = utg->id();
        nfrg++;
      }
    }
  }

  int lost = 0;
  int found = 0;

  for (int i=0; i<_fi->numFragments()+1; i++) {
    if (_fi->fragmentLength(i) > 0) {
      if (inUnitig[i] == 0) {
        fprintf(stderr, "ERROR frag %d is in unitig 0!\n", i);
      } else if (inUnitig[i] != 987654321) {
        found++;
      } else {
        fprintf(stderr, "ERROR frag %d disappeared!\n", i);
        lost++;
      }
    }
  }

  fprintf(stderr, "checkUnitigMembership()-- nutg=%d nfrg=%d lost=%d found=%d\n", nutg, nfrg, lost, found);

  assert(lost == 0);
};


//  For every unitig, report the best overlaps contained in the
//  unitig, and all overlaps contained in the unitig.
void
UnitigGraph::reportOverlapsUsed(const char *filename) {

  return;

  FILE *F = fopen(filename, "w");

  if (F == NULL)
    return;

  for (int  ti=0; ti<unitigs->size(); ti++) {
    Unitig  *utg = (*unitigs)[ti];

    if (utg == NULL)
      continue;

    for (int fi=0; fi<utg->dovetail_path_ptr->size(); fi++) {
      DoveTailNode  *frg = &(*utg->dovetail_path_ptr)[fi];

      //  Where is our best overlap?  Contained or dovetail?

      BestEdgeOverlap *bestedge5 = bog_ptr->getBestEdgeOverlap(frg->ident, FIVE_PRIME);
      BestEdgeOverlap *bestedge3 = bog_ptr->getBestEdgeOverlap(frg->ident, THREE_PRIME);

      int              bestident5 = 0;
      int              bestident3 = 0;

      if (bestedge5)
        bestident5 = bestedge5->frag_b_id;

      if (bestedge3)
        bestident3 = bestedge3->frag_b_id;

      //  Now search ahead, reporting any overlap to any fragment.
      //
      for (int oi=fi+1; oi<utg->dovetail_path_ptr->size(); oi++) {
        DoveTailNode  *ooo = &(*utg->dovetail_path_ptr)[oi];

        int frgbgn = MIN(frg->position.bgn, frg->position.end);
        int frgend = MAX(frg->position.bgn, frg->position.end);

        int ooobgn = MIN(ooo->position.bgn, ooo->position.end);
        int oooend = MAX(ooo->position.bgn, ooo->position.end);

        if ((frgbgn <= ooobgn) && (ooobgn + 40 < frgend)) {
          BestContainment *bestcont  = bog_ptr->getBestContainer(ooo->ident);

          int              bestident  = 0;
          if (bestcont)
            bestident = bestcont->container;

          bool isBest = ((frg->ident == bestident) ||
                         (ooo->ident == bestident5) ||
                         (ooo->ident == bestident3));

          fprintf(F, "%d\t%d%s\n", frg->ident, ooo->ident, (isBest) ? ((bestident) ? "\tbc" : "\tbe") : "");
        }

        if (frgend < ooobgn)
          break;
      }
    }
  }

  fclose(F);
}


UnitigGraph::UnitigGraph(FragmentInfo *fi, BestOverlapGraph *bp) {
  unitigs = new UnitigVector;
  _fi     = fi;
  bog_ptr = bp;
}

UnitigGraph::~UnitigGraph() {
  UnitigVector::iterator utg_itr;
  for(utg_itr=unitigs->begin(); utg_itr!=unitigs->end(); utg_itr++)
    delete *utg_itr;
  delete unitigs;
}


void UnitigGraph::build(ChunkGraph *cg_ptr, bool unitigIntersectBreaking, char *output_prefix) {
  bool verbose = false;

  // Initialize where we've been to nowhere; "Do not retraverse list"
  Unitig::resetFragUnitigMap( _fi->numFragments() );

  // Invert the containment map to key by container, instead of containee
  ContainerMap      cMap;

  for (iuid c=0; c<=_fi->numFragments(); c++)
    if (bog_ptr->isContained(c))
      cMap[bog_ptr->getBestContainer(c)->container].push_back(c);

  // Step through all the fragments

  fprintf(stderr, "==> BUILDING UNITIGS from %d fragments.\n", _fi->numFragments());

  iuid frag_idx;
  while( frag_idx = cg_ptr->nextFragByChunkLength() ) {

    if (_fi->fragmentLength(frag_idx) == 0)
      //  Deleted fragment.
      continue;

    if (Unitig::fragIn(frag_idx) != 0)
      //  Fragment already placed.
      continue;

    if (bog_ptr->isContained(frag_idx) == true)
      //  A contained fragment, will deal with later
      continue;

    //  Make new unitig, then follow edges off the 5' end.

    Unitig *utg = new Unitig(verbose);

    unitigs->push_back(utg);

    populateUnitig(utg,
                   frag_idx, FIVE_PRIME,
                   NULL_FRAG_ID, NULL,
                   verbose);

    //  Attempt to go off the 3' end

    BestEdgeOverlap   *tpBest   = bog_ptr->getBestEdgeOverlap(frag_idx, THREE_PRIME);

    if (tpBest->frag_b_id == NULL_FRAG_ID)
      //  No next fragment
      continue;

    if (Unitig::fragIn(tpBest->frag_b_id) == utg->id())
      //  Already in this unitig; deal with circularity later
      continue;

    if (Unitig::fragIn(tpBest->frag_b_id)) {
      //  Save the outgoing intersection edge for later splitting.
      fprintf(stderr,"unitigIntersect: unitig %5d frag %7d -> unitig %5d frag %7d\n",
              utg->id(), frag_idx, Unitig::fragIn(tpBest->frag_b_id), tpBest->frag_b_id);
      unitigIntersect[tpBest->frag_b_id].push_back(frag_idx);
      continue;
    }

    //  Otherwise, we can extend the unitig, walking off the other end
    //  of the unitig.

    utg->reverseComplement();

    fprintf(stderr, "continue unitig %d (length = %d, lastfraglen = %d) with frag %d (hang %d %d)\n",
            utg->id(), utg->getLength(), _fi->fragmentLength(frag_idx),
            tpBest->frag_b_id, tpBest->ahang, -tpBest->bhang);

    populateUnitig(utg,
                   tpBest->frag_b_id, (tpBest->bend == THREE_PRIME) ? FIVE_PRIME : THREE_PRIME,
                   frag_idx, tpBest,
                   verbose);
  }

  reportOverlapsUsed("overlaps.afterbuild");

  fprintf(stderr, "==> BUILDING UNITIGS catching missed fragments.\n");

  // Pick up frags missed above, possibly circular unitigs
  for(frag_idx=1; frag_idx<=_fi->numFragments(); frag_idx++){

    if (_fi->fragmentLength(frag_idx) == 0)
      continue;

    if (Unitig::fragIn(frag_idx) > 0)
      continue;

    iuid fp_dst_frag_id = bog_ptr->getBestEdgeOverlap(frag_idx, FIVE_PRIME) ->frag_b_id;
    iuid tp_dst_frag_id = bog_ptr->getBestEdgeOverlap(frag_idx, THREE_PRIME)->frag_b_id;

    if (bog_ptr->isContained(frag_idx) == 0)
      fprintf(stderr, "frag %d missed; fp_dst_frag_id=%d tp_dst_frag_id=%d contained=%d\n",
              frag_idx, fp_dst_frag_id, tp_dst_frag_id, bog_ptr->isContained(frag_idx));

    if ((fp_dst_frag_id != NULL_FRAG_ID) &&
        (tp_dst_frag_id != NULL_FRAG_ID)) {

      Unitig *utg = new Unitig(true);

      populateUnitig(utg, frag_idx, FIVE_PRIME, NULL_FRAG_ID, NULL, verbose);

      unitigs->push_back(utg);
    } else {
      // Should both be null or neither null otherwise main loop
      // failed to find it
      assert(fp_dst_frag_id == NULL_FRAG_ID);
      assert(tp_dst_frag_id == NULL_FRAG_ID);
    }
  }

  reportOverlapsUsed("overlaps.afterbuild2");

  fprintf(stderr, "==> BREAKING UNITIGS.\n");

  if (unitigIntersectBreaking)
    breakUnitigs(cMap, output_prefix);

  reportOverlapsUsed("overlaps.afterbreak");

  fprintf(stderr, "==> PLACING CONTAINED FRAGMENTS\n");

  for (int  ti=0; ti<unitigs->size(); ti++) {
    Unitig  *thisUnitig = (*unitigs)[ti];
    if (thisUnitig)
      thisUnitig->recomputeFragmentPositions(cMap, bog_ptr);
  }

  reportOverlapsUsed("overlaps.aftercontains");

  ////////////////////////////////////////
  //
  //  This is a huge hack to get around a bug somewhere before us.  We
  //  seem to not be placing all fragments.  So, run through all
  //  fragments, and for anything not placed, toss it into a new
  //  unitig.  Let scaffolder figure out where to put it.
  //
  //  Notice the similarity between this code and checkUnitigMembership().
  //
  {
    fprintf(stderr, "==> SEARCHING FOR ZOMBIES\n");

    iuid   *inUnitig   = new iuid [_fi->numFragments()+1];
    int     numZombies = 0;

    //  Sorry, poor fella that has more than 987,654,321 unitigs.  I
    //  can't imagine the rest of the pipeline would run though.
    //
    //  Mark fragments as dead.
    //
    for (int i=0; i<_fi->numFragments()+1; i++)
      inUnitig[i] = 987654321;

    //  ZZZzzzzaapaapppp!  IT'S ALIVE!
    //
    for (int  ti=0; ti<unitigs->size(); ti++) {
      Unitig  *utg = (*unitigs)[ti];

      if (utg)
        for (DoveTailIter it=utg->dovetail_path_ptr->begin(); it != utg->dovetail_path_ptr->end(); it++)
          inUnitig[it->ident] = utg->id();
    }

    //  Anything still dead?
    //
    for (int i=0; i<_fi->numFragments()+1; i++) {
      if (_fi->fragmentLength(i) > 0) {
        if (inUnitig[i] == 0) {
          //  We'll catch this error inna second in checkUnitigMembership().
        } else if (inUnitig[i] != 987654321) {
          //  We'll count this inna second there too.
        } else {
          //  Ha!  Gotcha!  You're now a resurrected brain eating
          //  zomibie?  Some day we'll figure out how to put you in
          //  properly.  For now, enjoy the ride.

          Unitig *utg = new Unitig(false);

          DoveTailNode frag;

          frag.type         = AS_READ;
          frag.ident        = i;
          frag.contained    = 0;
          frag.delta_length = 0;
          frag.delta        = NULL;

          frag.position.bgn = 0;
          frag.position.end = _fi->fragmentLength(i);

          utg->addFrag(frag, 0, true);
          unitigs->push_back(utg);

          numZombies++;
        }
      }
    }

    if (numZombies > 0)
      fprintf(stderr, "RESURRECTED %d ZOMBIE FRAGMENTS.\n", numZombies);

    delete inUnitig;
  }

  checkUnitigMembership();
}


void
UnitigGraph::setParentAndHang(ChunkGraph *cg) {

  for (int  ti=0; ti<unitigs->size(); ti++) {
    Unitig        *utg = (*unitigs)[ti];
    DoveTailNode  *frg;

    if (utg == NULL)
      continue;

    if (utg->dovetail_path_ptr->size() == 0)
      continue;

    //  First fragment has no parent or hangs set.
    frg = &(*utg->dovetail_path_ptr)[0];

    frg->parent       = 0;
    frg->ahang        = 0;
    frg->bhang        = 0;

    //  All other fragments can have a parent, if the best overlap is
    //  in our unitig.
    //
    for (int fi=1; fi<utg->dovetail_path_ptr->size(); fi++) {
      frg = &(*utg->dovetail_path_ptr)[fi];

      //  Where is our best overlap?  Contained or dovetail?

      BestEdgeOverlap *bestedge5 = bog_ptr->getBestEdgeOverlap(frg->ident, FIVE_PRIME);
      BestEdgeOverlap *bestedge3 = bog_ptr->getBestEdgeOverlap(frg->ident, THREE_PRIME);
      BestContainment *bestcont  = bog_ptr->getBestContainer(frg->ident);

      if (frg->contained) {
        //  Fragment is contained.  Don't use the edges off the ends.
        bestedge5 = NULL;
        bestedge3 = NULL;
      } else {
        //  Fragment is dovetail.  Don't use the best containment.
        bestcont == NULL;
      }

      if (frg->position.bgn < frg->position.end) {
        //  Fragment is forward.  Looking for the best edge off of our 5' end.
        bestedge3 = NULL;
      } else {
        //  Fragment is reverse.  Looking for the best edge off of our 3' end.
        bestedge5 = NULL;
      }


      //  Initialize to 'unset'.
      frg->parent       = 0;
      frg->ahang        = 0;
      frg->bhang        = 0;


      //  Now search for the correct overlap
      //
      for (int oi=fi-1; oi>=0; oi--) {
        DoveTailNode  *ooo = &(*utg->dovetail_path_ptr)[oi];

        if (bestcont) {
          //  Frag is contained.

          if (ooo->ident == bestcont->container) {
            //  And the container is here too!
            frg->parent = bestcont->container;

            //  The hangs assume the container is forward; adjust if not so.
            if (ooo->position.bgn < ooo->position.end) {
              frg->ahang  = bestcont->a_hang;
              frg->bhang  = bestcont->b_hang;
            } else {
              frg->ahang  = -bestcont->b_hang;
              frg->bhang  = -bestcont->a_hang;
            }
          }

        } else {
          //  Frag is not contained.
          if ((bestedge5) && (ooo->ident == bestedge5->frag_b_id)) {
            //  Our fragment is forward.  Reverse the overlap.
            //
            //  reality  A -----------      bestedge  B    -------->
            //                B --------->            A -------
            //  
            frg->parent = bestedge5->frag_b_id;
            frg->ahang  = -bestedge5->ahang;
            frg->bhang  = -bestedge5->bhang;

            assert(frg->ahang >= 0);
          }
          if ((bestedge3) && (ooo->ident == bestedge3->frag_b_id)) {
            //  Our fragment is backward.  Reverse the overlap.
            //
            //  reality  A ---------       bestedge  B ------->
            //                 B <-------            A    --------->
            //
            frg->parent = bestedge3->frag_b_id;
            frg->ahang  = bestedge3->bhang;
            frg->bhang  = bestedge3->ahang;

            assert(frg->ahang >= 0);
          }
        }
      }

      //  If we didn't find our best, pick any thing that overlaps us.
      //  This should only happen for contained fragments.

#if 0
      //  Unfortunately, we didn't save the overlap info, so we cannot
      //  set the hangs.

      if (frg->parent == 0) {
        for (int oi=fi-1; (frg->parent == 0) && (oi >= 0); oi--) {
          if (tigGraph.bog_ptr->containHaveEdgeTo(frg->ident, oi)) {
            frg->parent = oi;
            frg->ahang  = 0;
            frg->bhang  = 0;
          }
        }
      }
#endif

    }
  }
}


void
UnitigGraph::populateUnitig(Unitig           *unitig,
                            iuid              firstFragID,
                            fragment_end_type walkEnd,
                            iuid              lastID,
                            BestEdgeOverlap  *lastEdge,
                            bool              verbose){
  iuid              fragID   = firstFragID;
  iuid              nextID   = 0;
  int               fragBgn  = 0;
  int               fragEnd  = 0;
  BestEdgeOverlap  *nextEdge = NULL;
  DoveTailNode     frag;

  fprintf(stderr, "populateUnitig()--  STARTS for id %d\n", unitig->id());

  if (fragID == NULL_FRAG_ID)
    return;

  if (Unitig::fragIn(fragID) != 0)
    return;

  memset(&frag, 0, sizeof(DoveTailNode));

  //  Decide on a placement for the end of this fragment.  Easy if
  //  there are no existing fragments, just the length.
  //
  if (unitig->getLength() <= 0) {
    fragBgn = 0;
    fragEnd = _fi->fragmentLength(fragID);
    fprintf(stderr, "new unitig; new frag at %d - %d\n", fragBgn, fragEnd);
  }

  //  With existing fragments, we need to grab the last fragment, get
  //  the edge that got us to this fragment, and use the bhang there
  //  to place the end.
  //
  if (unitig->getLength() > 0) {
    int ahang = (lastEdge->ahang > 0) ? lastEdge->ahang : -lastEdge->bhang;
    int bhang = (lastEdge->ahang > 0) ? lastEdge->bhang : -lastEdge->ahang;

    DoveTailNode last = unitig->dovetail_path_ptr->back();

    assert(last.ident == lastID);

    int bgn = (last.position.bgn < last.position.end) ? last.position.bgn : last.position.end;
    int end = (last.position.bgn < last.position.end) ? last.position.end : last.position.bgn;

    //  Deep coverage (and lots of contains) result in placement based
    //  on hangs being wrong.  We don't know the length of the
    //  overlap, and any indel in the overlap change the amount
    //  hanging.  Reset the end to be based on the length of the fragment.
    //
    //  This is also reset in AS_BOG_Unitig.c::placeContains().

    fragBgn = bgn + ahang;
    fragEnd = end + bhang;

#warning not knowing the overlap length really hurts.
    if (fragBgn < fragEnd)
      fragEnd = fragBgn + _fi->fragmentLength(fragID);
    else
      fragBgn = fragEnd + _fi->fragmentLength(fragID);

    fprintf(stderr, "existing unitig; new frag at %d - %d\n", fragBgn, fragEnd);
  }

  do {

    //  Grab the edge to the next fragment.

    nextEdge = bog_ptr->getBestEdgeOverlap(fragID, walkEnd);
    nextEdge->print(stderr);

    //  We only work with dovetails.

    assert(((nextEdge->ahang >= 0) && (nextEdge->bhang >= 0)) ||
           ((nextEdge->ahang <= 0) && (nextEdge->bhang <= 0)));

    //  Our ahang and bhang are relative to fragID.

    int ahang = (nextEdge->ahang > 0) ? nextEdge->ahang : -nextEdge->bhang;
    int bhang = (nextEdge->ahang > 0) ? nextEdge->bhang : -nextEdge->ahang;

    frag.type         = AS_READ;
    frag.ident        = fragID;
    frag.contained    = 0;
    frag.parent       = nextEdge->frag_b_id;
    frag.sourceInt    = -1;
    frag.ahang        = ahang;
    frag.bhang        = bhang;
    frag.position.bgn = (walkEnd == FIVE_PRIME) ? fragEnd : fragBgn;
    frag.position.end = (walkEnd == FIVE_PRIME) ? fragBgn : fragEnd;
    frag.delta_length = 0;
    frag.delta        = NULL;

    unitig->addFrag(frag, 0, verbose);

    //  Find the start of the next fragment.  See above comment.

    fragBgn += ahang;
    fragEnd += bhang;

#warning not knowing the overlap length really hurts.
    if (fragBgn < fragEnd)
      fragEnd = fragBgn + _fi->fragmentLength(nextEdge->frag_b_id);
    else
      fragBgn = fragEnd + _fi->fragmentLength(nextEdge->frag_b_id);

    // Move to next fragment in the chain

    lastID   = fragID;
    fragID   = nextEdge->frag_b_id;
    walkEnd  = (nextEdge->bend == FIVE_PRIME) ? THREE_PRIME : FIVE_PRIME;

    //  But stop if we're in a circle

    if (fragID == firstFragID)
      fragID = NULL_FRAG_ID;

  } while ((fragID != NULL_FRAG_ID) &&
           (Unitig::fragIn(fragID) == 0));

  fprintf(stderr, "populateUnitig()--  FINISHES for id %d (fragID=%d in unitig %d)\n", unitig->id(), fragID, Unitig::fragIn(fragID));

  //  Save this outgoing intersection point for future use
  //  (another circular unitig test here too)

  if ((lastID != NULL_FRAG_ID) &&
      (fragID != NULL_FRAG_ID) &&
      (Unitig::fragIn(fragID) != 0)) {
    fprintf(stderr,"unitigIntersect: unitig %5d frag %7d -> unitig %5d frag %7d\n",
            unitig->id(), lastID, Unitig::fragIn(fragID), fragID);
    unitigIntersect[fragID].push_back(lastID);
  }
}

void UnitigGraph::breakUnitigs(ContainerMap &cMap, char *output_prefix) {
  FILE *breakFile;

  {
    char name[FILENAME_MAX];
    sprintf(name, "%s.breaks.ovl", output_prefix);

    errno = 0;
    breakFile = fopen(name, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s' to write unitig breaks.\n", name), exit(1);
  }

  //  Debug output
  for (int  ti=0; ti<unitigs->size(); ti++) {
    Unitig        *tig   = (*unitigs)[ti];
    DoveTailNode   first = tig->dovetail_path_ptr->front();
    iuid           prev  = 0;
    DoveTailNode   last  = tig->getLastBackboneNode(prev);

    if ( prev == 0 )
      continue; // skip singletons

    iuid id1 = first.ident;
    iuid id2 = last.ident;

    assert(prev != id2);
    assert(last.contained == 0);

    BestEdgeOverlap *firstBest = bog_ptr->getBestEdgeOverlap(id1, (first.position.bgn < first.position.end) ? FIVE_PRIME  : THREE_PRIME);
    BestEdgeOverlap *lastBest  = bog_ptr->getBestEdgeOverlap(id2, (last.position.bgn  < last.position.end)  ? THREE_PRIME : FIVE_PRIME);

    // Skip non bubbles, anchored ends
    if (lastBest->frag_b_id == 0 || firstBest->frag_b_id == 0)
      continue;

    if (firstBest->frag_b_id == lastBest->frag_b_id)
      fprintf(stderr, "Bubble %d to %d len %d (self bubble %d to %d)\n", firstBest->frag_b_id, lastBest->frag_b_id, tig->getLength(), id1, id2);
    else
      fprintf(stderr, "Bubble %d to %d len %d\n", firstBest->frag_b_id, lastBest->frag_b_id, tig->getLength());
  }


  //  Stop when we've seen all current unitigs.  Replace tiMax
  //  in the for loop below with unitigs->size() to recursively
  //  split unitigs.
  int  tiMax = unitigs->size();

  for (int  ti=0; ti<tiMax; ti++) {
    Unitig             *tig = (*unitigs)[ti];
    UnitigBreakPoints   breaks;
    DoveTailNode        lastBackbone;

    int                 numFragsInUnitig = 0;
    int                 fragCount        = 0;
    int                 fragIdx;

    tig->sort();

    //  Count the number of fragments in this unitig, including
    //  yet-to-be-placed contained fragments.
    for (fragIdx=0; fragIdx<tig->dovetail_path_ptr->size(); fragIdx++) {
      DoveTailNode  *f = &(*tig->dovetail_path_ptr)[fragIdx];

      if (cMap.find(f->ident) != cMap.end())
        numFragsInUnitig += cMap[f->ident].size();

      numFragsInUnitig++;

      //  Contained fragments should not be placed yet.
      assert(f->contained == 0);
    }


    for (fragIdx=0; fragIdx<tig->dovetail_path_ptr->size(); fragIdx++) {
      DoveTailNode  *f = &(*tig->dovetail_path_ptr)[fragIdx];

      fragCount++;
      if (cMap.find(f->ident) != cMap.end())
        fragCount += cMap[f->ident].size();

#if 1
      //  First fragment?
      if (fragIdx == 0) {
        fragment_end_type  dtEnd = (isReverse(f->position)) ? THREE_PRIME : FIVE_PRIME;
        BestEdgeOverlap   *bEdge = bog_ptr->getBestEdgeOverlap(f->ident, dtEnd);

#ifdef VERBOSEBREAK
        if ((bEdge) && (bEdge->frag_b_id > 0))
          fprintf(stderr,"unitig %d %c' frag %d points to unitig %d frag %d\n",
                  tig->id(), (dtEnd == THREE_PRIME) ? '3' : '5', f->ident,
                  Unitig::fragIn(bEdge->frag_b_id), bEdge->frag_b_id);
#endif
      }

      //  Last fragment?
      if (fragIdx + 1 == tig->dovetail_path_ptr->size()) {
        fragment_end_type  dtEnd = (isReverse(f->position)) ? FIVE_PRIME : THREE_PRIME;
        BestEdgeOverlap   *bEdge = bog_ptr->getBestEdgeOverlap(f->ident, dtEnd);

#ifdef VERBOSEBREAK
        if ((bEdge) && (bEdge->frag_b_id > 0))
          fprintf(stderr,"unitig %d %c' frag %d points to unitig %d frag %d\n",
                  tig->id(), (dtEnd == THREE_PRIME) ? '3' : '5', f->ident,
                  Unitig::fragIn(bEdge->frag_b_id), bEdge->frag_b_id);
#endif
      }
#endif

      FragmentEdgeList::const_iterator edge_itr = unitigIntersect.find(f->ident);

      if (edge_itr != unitigIntersect.end()) {

        // We have a set of best edges incoming from other unitigs
        for (FragmentList::const_iterator fragItr = edge_itr->second.begin();
             fragItr != edge_itr->second.end();
             fragItr++) {
          iuid inFrag = *fragItr;

          // check if it's incoming frag's 5' best edge.  If not, it must be the 3' edge.
          fragment_end_type bestEnd  = FIVE_PRIME;
          BestEdgeOverlap  *bestEdge = bog_ptr->getBestEdgeOverlap(inFrag, bestEnd);

          if (bestEdge->frag_b_id != f->ident) {
            bestEnd  = THREE_PRIME;
            bestEdge = bog_ptr->getBestEdgeOverlap(inFrag, bestEnd);
            assert(bestEdge->frag_b_id == f->ident);
          }

          int pos = (bestEdge->bend == FIVE_PRIME) ? f->position.bgn : f->position.end;

#warning DANGEROUS assume unitig is at id-1 in vector
          Unitig *inTig = (*unitigs)[Unitig::fragIn(inFrag)-1];

          if (inTig) {
            UnitigBreakPoint breakPoint(f->ident, bestEdge->bend);

            breakPoint.fragPos     = f->position;
            breakPoint.fragsBefore = fragCount;
            breakPoint.fragsAfter  = numFragsInUnitig - fragCount;
            breakPoint.inSize      = inTig->getLength();
            breakPoint.inFrags     = inTig->getNumFrags();

            breaks.push_back(breakPoint);
          } else {
            //fprintf(stderr, "  NullBreak tig %5d at frag %7d\n", tig->id(), f->ident);
          }

#ifdef VERBOSEBREAK
          fprintf(stderr, "unitig %d (%d frags, len %d) frag %d end %c' into unitig %d frag %d end %c' pos %d\n",
                  Unitig::fragIn(inFrag),
                  inTig->getNumFrags(),
                  inTig->getLength(),
                  inFrag,
                  (bestEnd == FIVE_PRIME) ? '5' : '3',
                  tig->id(),
                  f->ident,
                  (bestEdge->bend == FIVE_PRIME) ? '5' : '3',
                  pos);
#endif
          {
            GenericMesg  pmesg;
            OverlapMesg  omesg;

            omesg.aifrag          = inFrag;
            omesg.bifrag          = f->ident;
            omesg.ahg             = bestEdge->ahang;
            omesg.bhg             = bestEdge->bhang;
            omesg.orientation     = AS_UNKNOWN;
            omesg.overlap_type    = AS_DOVETAIL;
            omesg.quality         = 0.0;
            omesg.min_offset      = 0;
            omesg.max_offset      = 0;
            omesg.polymorph_ct    = 0;
            omesg.alignment_trace = NULL;
#ifdef AS_MSG_USE_OVL_DELTA
            omesg.alignment_delta = NULL;
#endif

            if ((bestEnd == FIVE_PRIME) && (bestEdge->bend == FIVE_PRIME))
              omesg.orientation = AS_OUTTIE;
            if ((bestEnd == FIVE_PRIME) && (bestEdge->bend == THREE_PRIME))
              omesg.orientation = AS_ANTI;
            if ((bestEnd == THREE_PRIME) && (bestEdge->bend == FIVE_PRIME))
              omesg.orientation = AS_NORMAL;
            if ((bestEnd == THREE_PRIME) && (bestEdge->bend == THREE_PRIME))
              omesg.orientation = AS_INNIE;

            pmesg.t = MESG_OVL;
            pmesg.m = &omesg;

            WriteProtoMesg_AS(breakFile, &pmesg);
          }
        }
      }
    }

    if (breaks.empty() == false) {
      DoveTailNode  *f = &tig->dovetail_path_ptr->back();

      //  create a final fake bp for the last frag so we
      //  have a reference point to the end of the tig for
      //  filtering the breakpoints.  fakeEnd, for searching.
      //
      UnitigBreakPoint breakPoint(f->ident, (isReverse(f->position)) ? FIVE_PRIME : THREE_PRIME);

      //  +1 below seems like a bug, but it reproduces what was here
      //  before.
#warning possible bug

      breakPoint.fragPos     = f->position;
      breakPoint.fragsBefore = fragCount + 1;
      breakPoint.fragsAfter  = 0;
      breakPoint.inSize      = std::numeric_limits<int>::max();
      breakPoint.inFrags     = 0;

      breaks.push_back(breakPoint);

      UnitigVector* newUs = breakUnitigAt(cMap, tig, breaks);

      if (newUs != NULL) {
        delete tig;
        (*unitigs)[ti] = NULL;
        unitigs->insert(unitigs->end(), newUs->begin(), newUs->end());
      }

      delete newUs;
    }
  }

  fclose(breakFile);
}




float UnitigGraph::getGlobalArrivalRate(long total_random_frags_in_genome, long genome_size){

  float _globalArrivalRate;

  //  If the genome size has not been specified, estimate the GAR.
  if(genome_size == 0){

    float total_rho=0, avg_rho;
    float total_arrival_frags=0;
    size_t rho_gt_10000 = 0;

    // Go through all the unitigs to sum rho and unitig arrival frags
    UnitigVector::const_iterator iter;
    for(
        iter=unitigs->begin();
        iter!=unitigs->end();
        iter++){

      if (*iter == NULL)
        continue;

      avg_rho = (*iter)->getAvgRho(_fi);
      total_rho += avg_rho;
      if (avg_rho > 10000.0)
        rho_gt_10000 += (size_t)avg_rho / 10000;

      float unitig_random_frags = (*iter)->getNumRandomFrags();
      if (--unitig_random_frags < 0)
        unitig_random_frags = 0;

      total_arrival_frags += unitig_random_frags;
      (*iter)->setLocalArrivalRate( unitig_random_frags / avg_rho );
    }
    // Estimate GAR
    _globalArrivalRate = (total_rho > 0) ? (total_arrival_frags / total_rho): 0;

    std::cerr << "Calculated Global Arrival rate " << _globalArrivalRate <<
      std::endl;
    // Now recalculate based on big unitigs, copied from AS_CGB/AS_CGB_cgb.c
    if (rho_gt_10000 * 20000 > total_rho) {
      float min_10_local_arrival_rate          = _globalArrivalRate;
      float median_local_arrival_rate          = _globalArrivalRate;
      float max_local_arrival_rate             = _globalArrivalRate;
      float recalibrated_fragment_arrival_rate = _globalArrivalRate;
      size_t num_arrival_rates=0;
      int median_index;
      std::vector<float> arrival_rate_array(rho_gt_10000);

      for( iter=unitigs->begin(); iter!=unitigs->end(); iter++) {
        if (*iter == NULL)
          continue;
        avg_rho = (*iter)->getAvgRho(_fi);
        if (avg_rho > 10000.0) {
          const int num_10000 = (size_t)avg_rho / 10000;
          const float local_arrival_rate =
            (*iter)->getNumRandomFrags() / avg_rho;
          assert(num_10000 > 0);
          int i;
          for(i=0;i<num_10000;i++){
            assert(i < rho_gt_10000);
            arrival_rate_array[i] = local_arrival_rate;
            num_arrival_rates++;
          }
          if(num_arrival_rates > 0){
            float tmp_fragment_arrival_rate, max_diff_arrival_rate;
            float prev_arrival_rate, cur_arrival_rate, diff_arrival_rate;
            int max_diff_index;
            std::sort(arrival_rate_array.begin(),arrival_rate_array.end());
            min_10_local_arrival_rate = arrival_rate_array[num_arrival_rates / 10];
            median_index = (num_arrival_rates * 5) / 10;
            median_local_arrival_rate = arrival_rate_array[median_index];
            max_local_arrival_rate = arrival_rate_array[num_arrival_rates-1];
            recalibrated_fragment_arrival_rate =
              arrival_rate_array[(num_arrival_rates * 19) / 20];
            prev_arrival_rate = min_10_local_arrival_rate;
            max_diff_arrival_rate = 0.0;
            for(i=num_arrival_rates / 10;i<median_index;i++){
              cur_arrival_rate = arrival_rate_array[i];
              diff_arrival_rate = cur_arrival_rate - prev_arrival_rate;
              prev_arrival_rate = cur_arrival_rate;
              if(diff_arrival_rate > max_diff_arrival_rate){
                max_diff_arrival_rate = diff_arrival_rate;
              }
            }
            max_diff_arrival_rate *= 2.0;
            max_diff_index = num_arrival_rates - 1;
            for(i=median_index;i<num_arrival_rates;i++){
              cur_arrival_rate = arrival_rate_array[i];
              diff_arrival_rate = cur_arrival_rate - prev_arrival_rate;
              prev_arrival_rate = cur_arrival_rate;
              if(diff_arrival_rate > max_diff_arrival_rate){
                max_diff_arrival_rate = diff_arrival_rate;
                max_diff_index = i - 1;
                break;
              }
            }
            max_diff_arrival_rate = arrival_rate_array[max_diff_index];
            tmp_fragment_arrival_rate =MIN(min_10_local_arrival_rate * 2.0,
                                           median_local_arrival_rate * 1.25);
            if(tmp_fragment_arrival_rate < recalibrated_fragment_arrival_rate){
              recalibrated_fragment_arrival_rate = tmp_fragment_arrival_rate;
            }
            if(max_diff_arrival_rate < recalibrated_fragment_arrival_rate){
              recalibrated_fragment_arrival_rate = max_diff_arrival_rate;
            }
          }
          if(recalibrated_fragment_arrival_rate > _globalArrivalRate){
            _globalArrivalRate = recalibrated_fragment_arrival_rate;
            std::cerr <<
              "Used recalibrated global_fragment_arrival_rate="
                      << _globalArrivalRate << std::endl
                      << "Used recalibrated global_fragment_arrival_distance="
                      <<((_globalArrivalRate > 0.) ? 1./(_globalArrivalRate) : 0.)
                      << std::endl
                      << "Chunk arrival rates sorted at 1/100s"
                      << std::endl ;
            for(i=0;i<100;i++) {
              std::cerr<<arrival_rate_array[((num_arrival_rates * i) / 100)]
                       << std::endl;
            }
            std::cerr << max_local_arrival_rate << std::endl;
          }
        }
      }
    }

   std::cerr << 
      "Computed genome_size="
      << (_globalArrivalRate > 0.f ? ((float)total_random_frags_in_genome / _globalArrivalRate) : 0.f) 
      << std::endl;
  }else{
    // Compute actual GAR
    _globalArrivalRate = (float)total_random_frags_in_genome / (float)genome_size;
  }

  return(_globalArrivalRate);

}






UnitigBreakPoint UnitigGraph::selectSmall(ContainerMap &cMap,
                                          const Unitig *tig,
                                          const UnitigBreakPoints &smalls,
                                          const UnitigBreakPoint& big,
                                          int   &lastBPCoord,
                                          int   &lastBPFragNum) {
  UnitigBreakPoint selection;

  double difference = 0.0;
  int  rFrgs   = big.fragsBefore;
  bool bRev    = isReverse(big.fragPos);
  int bContain = 0;

  int right    = (big.fragEnd.fragEnd() == FIVE_PRIME) ? big.fragPos.bgn : big.fragPos.end;

  if (cMap.find(big.fragEnd.fragId()) != cMap.end())
    bContain = cMap[big.fragEnd.fragId()].size();

  if (((bRev == true)  && (big.fragEnd == THREE_PRIME)) ||
      ((bRev == false) && (big.fragEnd == FIVE_PRIME)))
    rFrgs -= 1 + bContain;

  for(UnitigBreakPoints::const_iterator sIter = smalls.begin(); sIter != smalls.end(); sIter++) {
    UnitigBreakPoint small = *sIter;

    if (small.fragEnd.fragId() == big.fragEnd.fragId())
      continue;

    bool rev     = isReverse(small.fragPos);
    int sContain = 0;
    int lFrgs    = small.fragsBefore;
    iuid sid     = small.fragEnd.fragId();

    if (cMap.find(sid) != cMap.end())
      sContain = cMap[sid].size();

    // left side of the frag in the unitig, don't count it
    if (((rev == true)  && (small.fragEnd == THREE_PRIME)) ||
        ((rev == false) && (small.fragEnd == FIVE_PRIME)))
      lFrgs -= 1 + sContain;

    if (rFrgs - lFrgs == 1)
      continue;

    // use middle of frag instead of end, to get some overlap
    double bp    = (small.fragPos.end + small.fragPos.bgn) / 2.0;

    double lRate = (lFrgs - lastBPFragNum) / (bp - lastBPCoord);
    double rRate = (rFrgs - lFrgs) / (right - bp);
    double ratio = (lRate > rRate) ? lRate / rRate : rRate / lRate;
    double diff  = fabs( lRate - rRate);

    if ((ratio > 1.8) && (diff > difference)) {
      //fprintf(stderr, "Break frg %7d b %4d l %4d pos b %5d e %5.0f lRate %.4f\n", sid, lastBPFragNum, lFrgs, lastBPCoord, bp, lRate );
      //fprintf(stderr, "     diff %4d r %4d pos %5d rRate %.4f ratio %.2f to frag %7d\n", rFrgs - lFrgs, rFrgs, right, rRate, ratio, big.fragEnd.fragId());
      //fprintf(stderr,"     select frg %d for break on arrival rate diff %.4f\n", sid, diff);
      difference = diff;
      selection = small;
    }
  }

  lastBPCoord   = right;
  lastBPFragNum = rFrgs;

  return selection;
}



void UnitigGraph::filterBreakPoints(ContainerMap &cMap,
                                    Unitig *tig,
                                    UnitigBreakPoints &breaks) {

  int lastBPCoord = 0;
  int lastBPFragNum = 0;

  // Check for sentinal at end, not added by MateChecker
  UnitigBreakPoint fakeEnd = breaks.back();
  if (fakeEnd.inSize == std::numeric_limits<int>::max())
    breaks.pop_back();

  UnitigBreakPoints smallBPs;
  UnitigBreakPoints newBPs;

  for(UnitigBreakPoints::iterator iter = breaks.begin(); iter != breaks.end(); iter++) {
    UnitigBreakPoint nextBP = *iter;

    if ((nextBP.inFrags > UnitigGraph::MIN_BREAK_FRAGS) && (nextBP.inSize > UnitigGraph::MIN_BREAK_LENGTH)) {

      // big one, compare against smalls -- if we haven't
      // seen this break point already.
      //
      if (newBPs.empty() || (nextBP.fragEnd != newBPs.back().fragEnd)) {

        if (smallBPs.empty() == false) {
          //  smalls exist, select one.
          UnitigBreakPoint small = selectSmall(cMap, tig, smallBPs, nextBP, lastBPCoord, lastBPFragNum);
          if (small.fragEnd.fragId() > 0)
            newBPs.push_back(small);
          smallBPs.clear();

        } else {
          //  No smalls.  Update state to move past the current big breakpoint.
          lastBPCoord   = (nextBP.fragEnd.fragEnd() == FIVE_PRIME) ? nextBP.fragPos.bgn : nextBP.fragPos.end;
          lastBPFragNum = nextBP.fragsBefore;

          int bContain = 0;

          if (cMap.find(nextBP.fragEnd.fragId()) != cMap.end())
            bContain = cMap[nextBP.fragEnd.fragId()].size();

          bool bRev = isReverse(nextBP.fragPos);

          if (((bRev == true)  && (nextBP.fragEnd == THREE_PRIME)) ||
              ((bRev == false) && (nextBP.fragEnd == FIVE_PRIME)))
            lastBPFragNum -= 1 + bContain;
        }

        //  Haven't seen this break, so add it.
        newBPs.push_back( nextBP );
      }
    } else {
      //  Not a big breaker.  Save the small breaker if we've not seen it yet.
      if (newBPs.empty() || ((nextBP.fragEnd != newBPs.back().fragEnd) &&
                             (nextBP.fragEnd != smallBPs.back().fragEnd)))
        smallBPs.push_back(nextBP);
    }
  }

  //  If we've got small ones saved, select one.

  if (smallBPs.empty() == false) {
    UnitigBreakPoint small = selectSmall(cMap, tig, smallBPs, fakeEnd, lastBPCoord, lastBPFragNum);
    if (small.fragEnd.fragId() > 0)
      newBPs.push_back(small);
    smallBPs.clear();
  }

  breaks = newBPs;
}

// Doesn't handle contained frags yet
UnitigVector* UnitigGraph::breakUnitigAt(ContainerMap &cMap,
                                         Unitig *tig,
                                         UnitigBreakPoints &breaks) {

  if (breaks.empty())
    return NULL;

  // remove small break points
  filterBreakPoints(cMap, tig, breaks);

  // if we filtered all the breaks out return an empty list
  if (breaks.empty())
    return NULL;

  //  This is to catch a bug in....something before us.
  //  Sometimes the breakpoint list contains nonsense that tells
  //  us to split a unitig before the first fragment, or after
  //  the last.
  //
  //  That code is pretty dense (and big) so instead we'll make
  //  one pass through the list of breaks and see if we'd do any
  //  splitting.
  //
  //  It's a huge code duplication of the main loop.
  //
  {
    UnitigBreakPoints breakstmp = breaks;

    UnitigBreakPoint  breakPoint = breakstmp.front();
    breakstmp.pop_front();

    UnitigBreakPoint nextBP;

    int   newUnitigsConstructed = 0;
    int   newUnitigExists = 0;

    if (!breakstmp.empty())
      nextBP = breakstmp.front();

    for(DoveTailIter dtIter = tig->dovetail_path_ptr->begin(); dtIter != tig->dovetail_path_ptr->end(); dtIter++) {
      DoveTailNode frg = *dtIter;

      bool bothEnds = false;
      while ( !breakstmp.empty() && nextBP.fragEnd.fragId() == breakPoint.fragEnd.fragId() ) {
        if (nextBP.fragEnd.fragEnd() != breakPoint.fragEnd.fragEnd())
          bothEnds = true;
        breakstmp.pop_front();
        if (!breakstmp.empty())
          nextBP = breakstmp.front();
      }

      bool reverse = isReverse(frg.position);

      if (breakPoint.fragEnd.fragId() == frg.ident) {

        if (bothEnds) {
          newUnitigsConstructed++;
          newUnitigExists = 0;
        }

        else if (breakPoint.fragEnd.fragEnd() ==  FIVE_PRIME && !reverse ||
                 breakPoint.fragEnd.fragEnd() == THREE_PRIME && reverse) {
          newUnitigsConstructed++;
          newUnitigExists = 1;
        }

        else if (breakPoint.fragEnd.fragEnd() ==  FIVE_PRIME && reverse ||
                 breakPoint.fragEnd.fragEnd() == THREE_PRIME && !reverse ) {
          if (newUnitigExists == 0)
            newUnitigsConstructed++;
          newUnitigExists = 0;
        } else {
          assert(0);
        }

        if (breakstmp.empty()) {
          breakPoint.fragEnd = FragmentEnd();
        } else {
          breakPoint = breakstmp.front();
          breakstmp.pop_front();
          if (!breakstmp.empty())
            nextBP = breakstmp.front();
        }
      } else {
        if (newUnitigExists == 0) {
          newUnitigsConstructed++;
          newUnitigExists = 1;
        }
      }
    }

    if (newUnitigsConstructed < 2) {
      fprintf(stderr, "SPLITTING BUG DETECTED!  Adjusting for it.\n");
      return NULL;
    }

  }  //  end of explicit split test




  UnitigVector *splits = new UnitigVector();
  Unitig       *newTig = NULL;
  int           offset = 0;

  //  Emit a log of unitig creation and fragment moves
  //
  bool          verbose = false;

  UnitigBreakPoint breakPoint = breaks.front();
  breaks.pop_front();

  UnitigBreakPoint nextBP;

  if (!breaks.empty())
    nextBP = breaks.front();

  for(DoveTailIter dtIter = tig->dovetail_path_ptr->begin(); dtIter != tig->dovetail_path_ptr->end(); dtIter++) {
    DoveTailNode frg = *dtIter;

    // reduce multiple breaks at the same fragment end down to one
    bool bothEnds = false;
    while ( !breaks.empty() && nextBP.fragEnd.fragId() == breakPoint.fragEnd.fragId() ) {
      if (nextBP.fragEnd.fragEnd() != breakPoint.fragEnd.fragEnd())
        bothEnds = true;
      breaks.pop_front();
      if (!breaks.empty())
        nextBP = breaks.front();
    }

    bool reverse = isReverse(frg.position);

    if (breakPoint.fragEnd.fragId() == frg.ident) {

      //  There is method to this madness.  The three if
      //  blocks below have the same two basic code blocks
      //  in various orders.
      //
      //  o One code block makes a new unitig.
      //
      //  o The other adds a fragment to it, possibly
      //    remembering the beginning offset.
      //
      //  Finally, we delay making a new unitig until we are
      //  absolutely sure we're going to put fragments into
      //  it.

      //  Notes on fixing the break/join bug -- the bug is when we have a big breaker
      //  come in and chop off a small guy at the end.  We should be joining the two
      //  large unitigs.
      //
      //  To fix this, BPW was going to store enough stuff in the breakPoint so that
      //  the case can be detected here.  We need to know how much is before/after
      //  (depends on which end of the unitig we're at) our incoming point, (e.g.,
      //  fragsBefore and fragsAfter), and the fragments involved (fragEnd and inEnd).
      //
      //  When detected, push fragEnd and inEnd onto a list of joins.  After all breaks
      //  are finished, use the list of joins to get the unitigs to join and, well,
      //  join them.

      if (bothEnds) {
        //
        //  Break at both ends, create a singleton for this fragment.
        //
        if (verbose)
          fprintf(stderr,"  Break tig %d at both ends of %d num %d\n", tig->id(), breakPoint.fragEnd.fragId(), breakPoint.fragsBefore);

        newTig = new Unitig(verbose);  //  always make a new unitig, we put a frag in it now
        splits->push_back( newTig );

        if (newTig->dovetail_path_ptr->empty())
          offset = reverse ? -frg.position.end : -frg.position.bgn;
        newTig->addFrag( frg, offset, verbose );

        newTig = NULL;  //  delay until we need to make it
      }

      else if (breakPoint.fragEnd.fragEnd() ==  FIVE_PRIME && !reverse ||
               breakPoint.fragEnd.fragEnd() == THREE_PRIME && reverse) {
        //
        //  Break at left end of frg, frg starts new tig
        //
        if (verbose)
          fprintf(stderr,"  Break tig %d before %d num %d\n", tig->id(), breakPoint.fragEnd.fragId(), breakPoint.fragsBefore);

        newTig = new Unitig(verbose);  //  always make a new unitig, we put a frag in it now
        splits->push_back( newTig );

        if (newTig->dovetail_path_ptr->empty())
          offset = reverse ? -frg.position.end : -frg.position.bgn;
        newTig->addFrag( frg, offset, verbose );
      }

      else if (breakPoint.fragEnd.fragEnd() ==  FIVE_PRIME && reverse ||
               breakPoint.fragEnd.fragEnd() == THREE_PRIME && !reverse ) {
        //
        //  Break at right end of frg, frg goes in existing tig, then make new unitig
        //

        if (newTig == NULL) {  //  delayed creation?
          newTig = new Unitig(verbose);
          splits->push_back( newTig );
        }

        if (newTig->dovetail_path_ptr->empty())
          offset = reverse ? -frg.position.end : -frg.position.bgn;
        newTig->addFrag( frg, offset, verbose );

        if (verbose)
          fprintf(stderr,"  Break tig %d after %d num %d\n", tig->id(), breakPoint.fragEnd.fragId(), breakPoint.fragsBefore);

        //  Delay making a new unitig until we need it.
        newTig = NULL;
      } else {
        // logically impossible!
        assert(0);
      }


      // Done breaking, continue adding remaining frags to newTig
      if (breaks.empty()) {
        breakPoint.fragEnd = FragmentEnd();
      } else {
        breakPoint = breaks.front();
        breaks.pop_front();
        if (!breaks.empty())
          nextBP = breaks.front();
      }
    } else {
      //
      //  frag is not a unitig break point fragment, add to current new tig
      //

      if (newTig == NULL) {  //  delayed creation?
        newTig = new Unitig(verbose);
        splits->push_back( newTig );
      }
      if (newTig->dovetail_path_ptr->empty())
        offset = reverse ? -frg.position.end : -frg.position.bgn;
      newTig->addFrag( frg, offset, verbose );
    }
  }

  return splits;
}



void UnitigGraph::writeIUMtoFile(char *fileprefix, int fragment_count_target){
  int         fragment_count         = 0;
  int         file_count             = 1;
  char        filename[FILENAME_MAX] = {0};
  int         iumiid                 = 0;
  GenericMesg mesg;

  int nf = 0;
  int nu = 0;

  // Open up the initial output file

  sprintf(filename, "%s_%03d.cgb", fileprefix, file_count++);
  FILE *file = fopen(filename,"w");
  assert(NULL != file);

  sprintf(filename, "%s.iidmap", fileprefix);
  FILE *iidm = fopen(filename,"w");
  assert(NULL != iidm);

  // Step through all the unitigs

  checkUnitigMembership();

  for (int  ti=0; ti<unitigs->size(); ti++) {
    Unitig  *utg = (*unitigs)[ti];

    if (utg == NULL) {
      //fprintf(stderr, "unitig %d is null, skip.\n", ti);
      continue;
    }

    if (utg->getNumFrags() == 0) {
      fprintf(stderr, "unitig %d HAS NO FRAGS?\n", ti);
      continue;
    }

    IntUnitigMesg *ium_mesg_ptr = new IntUnitigMesg;

    assert(utg->getLength() > 0);
    assert(utg->getNumFrags() == utg->dovetail_path_ptr->size());

    ium_mesg_ptr->iaccession    = iumiid++;
    ium_mesg_ptr->coverage_stat = utg->getCovStat(_fi);
    ium_mesg_ptr->microhet_prob = 1.01;
    ium_mesg_ptr->status        = AS_UNASSIGNED;
    ium_mesg_ptr->unique_rept   = AS_FORCED_NONE;
    ium_mesg_ptr->length        = utg->getLength();
    ium_mesg_ptr->consensus     = "";
    ium_mesg_ptr->quality       = "";
    ium_mesg_ptr->forced        = 0;
    ium_mesg_ptr->num_frags     = utg->getNumFrags();
    ium_mesg_ptr->f_list        = &(utg->dovetail_path_ptr->front());

    fprintf(iidm, "Unitig %d == IUM %d (with %d frags)\n", utg->id(), ium_mesg_ptr->iaccession, utg->getNumFrags());

    fragment_count += ium_mesg_ptr->num_frags;

    if ((fragment_count_target >= 0) &&
        (fragment_count >= fragment_count_target)) {
      fclose(file);
      sprintf(filename, "%s_%03d.cgb", fileprefix, file_count++);
      file = fopen(filename,"w");
      assert(NULL != file);
      fragment_count = ium_mesg_ptr->num_frags;
    }

    mesg.m = ium_mesg_ptr;
    mesg.t = MESG_IUM;

    nf += ium_mesg_ptr->num_frags;
    nu += 1;

    WriteProtoMesg_AS(file, &mesg);

    delete ium_mesg_ptr;
  }

  fclose(file);
  fclose(iidm);
}


//  For every unitig, report the best overlaps contained in the
//  unitig, and all overlaps contained in the unitig.
void
UnitigGraph::writeOVLtoFile(char *fileprefix) {
  char         filename[FILENAME_MAX] = {0};
  GenericMesg  pmesg;
  OverlapMesg  omesg;

  sprintf(filename, "%s.unused.ovl", fileprefix);
  FILE *file = fopen(filename, "w");
  assert(file != NULL);

  for (int  ti=0; ti<unitigs->size(); ti++) {
    Unitig  *utg = (*unitigs)[ti];

    if (utg == NULL)
      continue;

    for (int fi=0; fi<utg->dovetail_path_ptr->size(); fi++) {
      DoveTailNode  *frg = &(*utg->dovetail_path_ptr)[fi];

      //  Where is our best overlap?  Contained or dovetail?

      BestEdgeOverlap *bestedge5 = bog_ptr->getBestEdgeOverlap(frg->ident, FIVE_PRIME);
      BestEdgeOverlap *bestedge3 = bog_ptr->getBestEdgeOverlap(frg->ident, THREE_PRIME);

      int              bestident5 = 0;
      int              bestident3 = 0;

      if (bestedge5) {
        bestident5 = bestedge5->frag_b_id;

        if ((bestident5 > 0) && (utg->fragIn(bestident5) != utg->id())) {
          omesg.aifrag          = frg->ident;
          omesg.bifrag          = bestident5;
          omesg.ahg             = bestedge5->ahang;
          omesg.bhg             = bestedge5->bhang;
          omesg.orientation     = AS_UNKNOWN;
          omesg.overlap_type    = AS_DOVETAIL;
          omesg.quality         = 0.0;
          omesg.min_offset      = 0;
          omesg.max_offset      = 0;
          omesg.polymorph_ct    = 0;
          omesg.alignment_trace = NULL;
#ifdef AS_MSG_USE_OVL_DELTA
          omesg.alignment_delta = NULL;
#endif

          //  This overlap is off of the 5' end of this fragment.
          if (bestedge5->bend == FIVE_PRIME)
            omesg.orientation = AS_OUTTIE;
          if (bestedge5->bend == THREE_PRIME)
            omesg.orientation = AS_ANTI;

          pmesg.t = MESG_OVL;
          pmesg.m = &omesg;

          WriteProtoMesg_AS(file, &pmesg);
        }
      }

      if (bestedge3) {
        bestident3 = bestedge3->frag_b_id;

        if ((bestident3 > 0) && (utg->fragIn(bestident3) != utg->id())) {
          omesg.aifrag          = frg->ident;
          omesg.bifrag          = bestident3;
          omesg.ahg             = bestedge3->ahang;
          omesg.bhg             = bestedge3->bhang;
          omesg.orientation     = AS_UNKNOWN;
          omesg.overlap_type    = AS_DOVETAIL;
          omesg.quality         = 0.0;
          omesg.min_offset      = 0;
          omesg.max_offset      = 0;
          omesg.polymorph_ct    = 0;
          omesg.alignment_trace = NULL;
#ifdef AS_MSG_USE_OVL_DELTA
          omesg.alignment_delta = NULL;
#endif

          //  This overlap is off of the 3' end of this fragment.
          if (bestedge3->bend == FIVE_PRIME)
            omesg.orientation = AS_NORMAL;
          if (bestedge3->bend == THREE_PRIME)
            omesg.orientation = AS_INNIE;

          pmesg.t = MESG_OVL;
          pmesg.m = &omesg;

          WriteProtoMesg_AS(file, &pmesg);
        }
      }
    }
  }

  fclose(file);
}
