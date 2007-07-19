
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

static char CM_ID[] = "$Id: AS_FGB_main.c,v 1.16 2007-07-19 09:50:32 brianwalenz Exp $";

#include "AS_UTL_version.h"  
#include "AS_CGB_all.h"

int (*compare_edge_function)(const void *a, const void *b);

//  AS_CGB_fga.c
void fragment_graph_analysis(IntFragment_ID max_frag_iid,
                             Tfragment frags[],
                             Tedge     edges[],
                             FILE      *ffga);


//  AS_CGB_fgb.c
void transitive_edge_marking(TStateGlobals * gstate,
                             THeapGlobals  * heapva,
                             Tfragment     * frags,
                             Tedge         * edges,
                             TIntEdge_ID   * next_edge_obj,
                             const int walk_depth,
                             const int cutoff_fragment_end_degree,
                             const int work_limit_per_candidate_edge,
                             const IntFragment_ID iv_start);

//  AS_CGB_edgemate.c
void append_the_edge_mates (Tfragment frags[],
                            Tedge edges[],
                            TIntEdge_ID * next_edge_obj);



void
input_messages_from_a_file(FILE       *fovl,
                           Tfragment  frags[],
                           Tedge      edges[],
                           FragmentHashObject *afr_to_avx,
                           TIntEdge_ID       *next_edge_obj,
                           const int dvt_double_sided_threshold_fragment_end_degree,
                           const int con_double_sided_threshold_fragment_end_degree,
                           const int intrude_with_non_blessed_overlaps_flag,
                           const uint32 overlap_error_threshold);

void
process_gkp_store_for_fragments(char *gkpStoreName,
                                Tfragment   *frags,
                                Tedge       *edges,
                                FragmentHashObject *afr_to_avx,
                                IntFragment_ID    *min_frag_iid,
                                IntFragment_ID    *max_frag_iid);

void
process_ovl_store(char * OVL_Store_Path,
                  Tfragment  frags[],
                  Tedge      edges[],
                  FragmentHashObject *afr_to_avx,
                  TIntEdge_ID         next_edge_obj[],
                  const int dvt_double_sided_threshold_fragment_end_degree,
                  const int con_double_sided_threshold_fragment_end_degree,
                  const int intrude_with_non_blessed_overlaps_flag,
                  const uint32 overlap_error_threshold);
  



static void output_mesgs(const Tfragment frags[],
                         const Tedge     edges[],
                         FILE *fcgb) {

  // Output the OVL messages:

  const IntEdge_ID  nedge = GetNumEdges(edges);
  IntEdge_ID ie;
  for(ie=0;ie<nedge;ie++){
    OverlapMesg ovl_mesg;

    const IntFragment_ID avx = get_avx_edge(edges,ie);
    const int asx = get_asx_edge(edges,ie);
    const int ahg = get_ahg_edge(edges,ie);

    const IntFragment_ID bvx = get_bvx_edge(edges,ie);
    const int bsx = get_bsx_edge(edges,ie);
    const int bhg = get_bhg_edge(edges,ie);

    const Tnes    nes = get_nes_edge(edges,ie);
    const uint32  qua = get_qua_edge(edges,ie);

    const IntFragment_ID aid = get_iid_fragment(frags,avx);
    const IntFragment_ID bid = get_iid_fragment(frags,bvx);
    // Assembler internal Fragment ids.

    signed char delta[1] = {AS_ENDOF_DELTA_CODE};
      
    ovl_mesg.aifrag = aid;
    ovl_mesg.bifrag = bid;

    ovl_mesg.ahg = ahg;
    ovl_mesg.bhg = bhg;
    ovl_mesg.min_offset = ahg;
    ovl_mesg.max_offset = ahg;

    ovl_mesg.orientation =
      ( asx ?
        ( bsx ? AS_INNIE : AS_NORMAL ) :
        ( bsx ? AS_ANTI  : AS_OUTTIE ) );

    switch(nes) {
      case AS_CGB_DOVETAIL_EDGE:
      case AS_CGB_INTERCHUNK_EDGE:
      case AS_CGB_INTRACHUNK_EDGE:
      case AS_CGB_TOUCHES_CONTAINED_EDGE:
      case AS_CGB_BETWEEN_CONTAINED_EDGE:
        ovl_mesg.overlap_type = AS_DOVETAIL; break;
      case AS_CGB_REMOVED_BY_TRANSITIVITY_DVT:
        ovl_mesg.overlap_type = AS_DOVETAIL_CHORD; break;
      case AS_CGB_CONTAINED_EDGE:
        ovl_mesg.overlap_type = AS_CONTAINMENT; break;
      case AS_CGB_REMOVED_BY_TRANSITIVITY_CON:
        ovl_mesg.overlap_type = AS_CONTAINMENT_CHORD; break;
      default:
        fprintf(stderr,"Unexpected overlap edge type: nes=%d\n", nes);
        assert(FALSE);
    }

    ovl_mesg.quality = AS_OVS_decodeQuality(qua);
    ovl_mesg.polymorph_ct = 0;
    ovl_mesg.delta = delta;

    if(
       ((is_a_dvt_simple(ahg,bhg))&&(aid < bid))
       // Output only one dovetail overlap edge record per overlap.
       ||
       ((!is_a_dvt_simple(ahg,bhg))&&(asx))
       // Output only one containment overlap edge record per overlap.
       // Choose the NORMAL and INNIE orientations.
       ){
      GenericMesg   pmesg;
      pmesg.t = MESG_OVL;
      pmesg.m = &ovl_mesg;
      WriteProtoMesg_AS(fcgb,&pmesg);
    }
  }
}





static void process_one_ovl_file(const char Batch_File_Name[],
                                 TStateGlobals  * gstate,
                                 THeapGlobals   * heapva,
                                 FragmentHashObject * afr_to_avx,
                                 const int dvt_double_sided_threshold_fragment_end_degree,
                                 const int con_double_sided_threshold_fragment_end_degree,
                                 const int intrude_with_non_blessed_overlaps_flag,
                                 const uint32 overlap_error_threshold,
                                 const int check_point_level) {
  
  if( check_point_level == 0 ) {
    FILE *fovl = fopen(Batch_File_Name,"r");
    if(NULL == fovl){
      fprintf(stderr,"* Can not open input file %s\n",Batch_File_Name);
      exit(1);
    }
    
    input_messages_from_a_file(fovl,
                               heapva->frags,
                               heapva->edges,
                               afr_to_avx, 
                               heapva->next_edge_obj,
                               dvt_double_sided_threshold_fragment_end_degree,
                               con_double_sided_threshold_fragment_end_degree,
                               intrude_with_non_blessed_overlaps_flag,
                               overlap_error_threshold);
    
    fclose(fovl);
  }
}




static void delete_duplicate_edges
(/* Input Only */
 Tfragment frags[], 
 /* Modify */
 Tedge edges[],
 TIntEdge_ID *next_edge_obj
 )
{
  /* 
     This routine marks duplicated edges for deletion.
     It is only effective when the edges are fully sorted.
  */
  //const IntFragment_ID nfrag = GetNumFragments(frags);
  const IntEdge_ID nedge = GetNumEdges(edges);
  IntEdge_ID ie1=0;

  IntEdge_ID count=0;

  for(ie1=1; ie1 < nedge; ie1++) {
    const IntEdge_ID ie0 = ie1 - 1;
    Aedge edge0 = *GetVA_Aedge(edges,ie0);
    Aedge edge1 = *GetVA_Aedge(edges,ie1);
    edge0.blessed = FALSE;
    edge1.blessed = FALSE;
    edge0.reflected = FALSE;
    edge1.reflected = FALSE;
    {
      const int icompare = (*compare_edge_function)(&edge1,&edge0);
      
      assert(icompare >= 0);
      
      if( icompare == 0 ) {
        count++;

        if( is_a_dvt_edge(edges,ie1) == TRUE )
          set_nes_edge(edges,ie1,AS_CGB_REMOVED_BY_DUPLICATE_DVT);
        else
          set_nes_edge(edges,ie1,AS_CGB_REMOVED_BY_DUPLICATE_CON);
      }
    }
  }

  fprintf(stderr, " " F_IID " duplicate edges found.\n",count);

  if(count > 0)
    ResetToRangeVA_IntEdge_ID(next_edge_obj,0);
}



int compare_edge_weak(const void * const aa, const void * const bb) 
{
  // This comparison function is used with ANSI qsort() for sorting
  // the edges to form contiguous segments.
  
  // The lesser edge is the one we keep in the Reaper.
  int icom;
  Aedge *a = (Aedge *)aa;
  Aedge *b = (Aedge *)bb;

  icom = ((a->avx) - (b->avx));
  if( icom == 0 ) {
    icom = ((a->asx) - (b->asx));
    //if( icom == 0 ) {
    //icom = (b->blessed - a->blessed);
      if( icom == 0 ) {
        icom = ((a->ahg) - (b->ahg));
        // Favor the minimum ahg.
        if( icom == 0 ) {
          icom = ((b->bhg) - (a->bhg));
          // Favor the maximum bhg.
          
          if( icom == 0 ) {
            // The following is unnecessary, but useful for the binary
            // search in the adjaceny lists and regression output.
            icom = ((a->bvx) - (b->bvx));
            if( icom == 0 )
              icom = ((a->bsx) - (b->bsx));
          }
        } // End of regression stuff.
      }
      //}
  } 
  return icom ;
}


int compare_edge_strong(const void * const aa, const void * const bb) 
{
  // This comparison function is used with ANSI qsort() for sorting
  // the edges to form contiguous segments.
  
  // The lesser edge is the one we keep in the Reaper.
  int icom;
  Aedge *a = (Aedge *)aa;
  Aedge *b = (Aedge *)bb;

  icom = ((a->avx) - (b->avx));
  if( icom == 0 ) {
    icom = ((a->asx) - (b->asx));
    if( icom == 0 ) {
      icom = (b->blessed - a->blessed);
      if( icom == 0 ) {
        icom = ((a->ahg) - (b->ahg));
        // Favor the minimum ahg.
        if( icom == 0 ) {
          icom = ((b->bhg) - (a->bhg));
          // Favor the maximum bhg.
          
          if( icom == 0 ) {
            // The following is unnecessary, but useful for the binary
            // search in the adjaceny lists and regression output.
            icom = ((a->bvx) - (b->bvx));
            if( icom == 0 ) {
              icom = ((a->bsx) - (b->bsx));
              if( icom == 0 )
                icom = (a->reflected - b->reflected);
            }
          }
        } // End of regression stuff.
      }
    }
  }
  return icom ;
}



int main_fgb(TStateGlobals * gstate,
             THeapGlobals  * heapva,
             UnitiggerGlobals * rg) {
  int status = 0;
  int ierr = 0;

  /* The Store/Checkpoint */
  FragmentHashObject *afr_to_avx = NULL;

  int did_processing_phase_2 = FALSE;

  compare_edge_function = compare_edge_strong;

  
  // Re-hash the fragment IID to fragment VID mapping using the
  // fragments in the store.
  assert(afr_to_avx == NULL);
  afr_to_avx = build_FragmentHash(heapva->frags, rg->as_cgb_max_frag_iid );
  assert(afr_to_avx != NULL);

  // WARNING do we need to rebuild the next_edge_obj here?

  if((rg->frag_store) && (rg->frag_store[0] != '\0'))
    process_gkp_store_for_fragments(rg->frag_store,
                                    heapva->frags,
                                    heapva->edges,
                                    afr_to_avx,
                                    &(gstate->min_frag_iid),
                                    &(gstate->max_frag_iid));

  if (rg->ovl_files_list_fname != NULL)
    process_one_ovl_file(rg->ovl_files_list_fname,
                         gstate,
                         heapva,
                         afr_to_avx,
                         rg->dvt_double_sided_threshold_fragment_end_degree,
                         rg->con_double_sided_threshold_fragment_end_degree,
                         rg->intrude_with_non_blessed_overlaps_flag,
                         rg->overlap_error_threshold,
                         rg->check_point_level);
  


  // Process the blessed ovl file.
  if(rg->blessed_overlaps_input_filename){
    assert(0 == GetNumEdges(heapva->edges));

    process_one_ovl_file(rg->blessed_overlaps_input_filename,
                         gstate,
                         heapva,
                         afr_to_avx,
                         rg->dvt_double_sided_threshold_fragment_end_degree,
                         rg->con_double_sided_threshold_fragment_end_degree,
                         rg->intrude_with_non_blessed_overlaps_flag,
                         AS_OVS_encodeQuality(1.0),
                         rg->check_point_level);

    {
      const IntFragment_ID nfrag = GetNumFragments(heapva->frags);
      const IntEdge_ID nedge = GetNumEdges(heapva->edges);
      IntEdge_ID ie;
      for( ie=0; ie < nedge; ie ++) {
        const IntFragment_ID avx = get_avx_edge(heapva->edges,ie);
        const IntFragment_ID asx = get_asx_edge(heapva->edges,ie);
        // assert(nedge > ie);
        assert(nfrag > avx);
        set_blessed_edge(heapva->edges, ie, TRUE);
        set_blessed_vertex(heapva->frags, avx, asx, TRUE);
      }
    }
  }
  

  // Process the bubble smoothing ovl file.
  if(rg->bubble_overlaps_filename)
    process_one_ovl_file(rg->bubble_overlaps_filename,
                         gstate,
                         heapva,
                         afr_to_avx,
                         rg->dvt_double_sided_threshold_fragment_end_degree,
                         rg->con_double_sided_threshold_fragment_end_degree,
                         rg->intrude_with_non_blessed_overlaps_flag,
                         rg->overlap_error_threshold,
                         rg->check_point_level);


  if ((rg->OVL_Store_Path) && (rg->OVL_Store_Path[0] != '\0'))
    process_ovl_store(rg->OVL_Store_Path,
                      heapva->frags,
                      heapva->edges,
                      afr_to_avx,
                      heapva->next_edge_obj,
                      rg->dvt_double_sided_threshold_fragment_end_degree,
                      rg->con_double_sided_threshold_fragment_end_degree,
                      rg->intrude_with_non_blessed_overlaps_flag,
                      rg->overlap_error_threshold);

  destroy_FragmentHash(afr_to_avx);

  /////////////////////////////////////////////////////////////////

  //count_fragment_and_edge_labels ( heapva->frags, heapva->edges, "Before separate_fragments_as_solo_hanging_thru");

  separate_fragments_as_solo_hanging_thru(heapva->frags,heapva->edges);

  compare_edge_function = compare_edge_weak;

  // Now do not distinguish the blessed overlaps.

  reorder_edges( heapva->frags, heapva->edges, heapva->next_edge_obj);
  //count_fragment_and_edge_labels( heapva->frags, heapva->edges, "RISM_reorder_edges");
  //check_symmetry_of_the_edge_mates( heapva->frags, heapva->edges, heapva->next_edge_obj);
    
  append_the_edge_mates( heapva->frags, heapva->edges, heapva->next_edge_obj);

  // Currently the spur finding and micro-bubble smoothing code
  // needs the dovetail edge mates to exist.

  //count_fragment_and_edge_labels( heapva->frags, heapva->edges, "RISM_append_the_edge_mates");
  //check_symmetry_of_the_edge_mates( heapva->frags, heapva->edges, heapva->next_edge_obj);

  delete_duplicate_edges( heapva->frags, heapva->edges, heapva->next_edge_obj);
  //count_fragment_and_edge_labels( heapva->frags, heapva->edges, "RISM_delete_duplicate_edges");
  //check_symmetry_of_the_edge_mates( heapva->frags, heapva->edges, heapva->next_edge_obj);

  // Currently the spur finding code needs the dovetail edge mates to
  // exist.

  if( (GetNumFragments(heapva->frags) > 0) &&
      (GetNumEdges(heapva->edges) > 0) &&
      (rg->dechord_the_graph) &&
      (!did_processing_phase_2) 
      ) /* The FGB computations need a populated graph store. */ {

    transitive_edge_marking
      ( gstate, heapva, heapva->frags, heapva->edges, heapva->next_edge_obj,
        rg->walk_depth,
        rg->cutoff_fragment_end_degree,
        rg->work_limit_per_candidate_edge,
        rg->iv_start);
  }
  
  identify_early_spur_fragments( heapva->frags, heapva->edges);
  
  //count_fragment_and_edge_labels ( heapva->frags, heapva->edges, "Before contained_fragment_marking_frc");
  contained_fragment_marking_frc( heapva->frags, heapva->edges);
  //count_fragment_and_edge_labels ( heapva->frags, heapva->edges, "After contained_fragment_marking_frc");
  

  char  thePath[FILENAME_MAX];
  sprintf(thePath,"%s.fga.ckp",rg->Output_Graph_Store_Prefix);
  FILE *ffga = fopen(thePath,"w");
    
  fragment_graph_analysis(gstate->max_frag_iid,
                          heapva->frags,
                          heapva->edges,
                          ffga);
  fclose(ffga);
    

  if( rg->create_dump_file ) {
    FILE *folp = NULL;
    char strtmp[FILENAME_MAX];
    
    fprintf(stderr,"Opening dump file to write a batch of "
	    "ADT+IDT+OFG+OVL messages.\n");

    if(NULL == (folp = fopen(rg->Dump_File_Name,"w"))){
      fprintf(stderr,"* Can not open output file %s\n",rg->Dump_File_Name);
      exit(1);
    }

    output_mesgs (heapva->frags,
                  heapva->edges,
                  folp);
    fclose(folp);
  }
  
  return( status);
}
