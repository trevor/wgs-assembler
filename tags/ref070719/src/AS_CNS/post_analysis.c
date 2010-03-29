
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
/* $Id: post_analysis.c,v 1.17 2007-05-29 10:54:28 brianwalenz Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "assert.h"
#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "MultiAlignStore_CNS.h"
#include "MultiAlignment_CNS.h"

void usage (char *pgmname){
  fprintf(stderr,"USAGE: %s -f <frgStore> -o <outprefix> [-O] [-m <minlen>] [-l <iidlist>]\n"
	  "\t-o outprefix specifies prefix of output file names\n"
	  "\t-O instructs to write out protoIO messages of selected IDs\n"
	  "\t-m instructs to ignore contigs less than minlen columns\n"
	  "\t-l gives a list of contig IIDs to process\n"
	  ,pgmname);
  exit(-1);
}


int main(int argc, char *argv[])
{ GenericMesg *pmesg;
 IntScaffoldMesg *isf;
 IntContigPairs *pairs;
 IntConConMesg *contig;
 IntUnitigMesg *unitig;
 MultiAlignT *ma;
 int i;
 int isplaced = 1;
 GateKeeperStore *frag_store;
 FILE *pcs = NULL;
 FILE *pfs = NULL;
 FILE *out = NULL;
 FILE *sublist = NULL;
 char buffer[256];
 char *frgstore_name=NULL;
 char *sublist_file=NULL;
 char *outputfile_prefix=NULL;
 HashTable_AS  *tig_iids;
 HashTable_AS  *tig_iids_found;
 int64  this_id;
 int do_all = 1;
 int errflg=0;
 char ch;
 int min_len=0;

 while ( !errflg && 
	 ( (ch = getopt(argc, argv, "b:f:l:m:o:O")) != EOF)) {
   switch(ch){
   case 'f':
     frag_store = openGateKeeperStore(optarg, FALSE);
     break;
   case 'l':
     do_all=0;
     sublist_file = optarg;
     break;
   case 'm':
     min_len = atoi(optarg);
     assert(min_len>=0);
     break;
   case 'o':
     outputfile_prefix=optarg;
     break;
   case 'O':
     out = fopen("post_analysis.out","w");
     assert(out!=NULL);
     break;
   default:
     errflg=1;
     break;
   }
 }
 if(errflg || frag_store==NULL){
   usage(argv[0]);
 }

 assert(outputfile_prefix!=NULL);
 sprintf(buffer,"%s.pcs",outputfile_prefix);
 pcs = fopen(buffer,"w");
 sprintf(buffer,"%s.pfs",outputfile_prefix);
 pfs = fopen(buffer,"w");
 assert(pfs && pcs );

   
 if ( !do_all ) {
   char   string[1000];
   int    num_iids=0;
   sublist = fopen(sublist_file,"r");
   if( sublist == NULL )
     {
       fprintf( stderr, "Failed to open list file %s for reading.\n", argv[2] );
       exit(1);
     }
   num_iids = 0;
   while( fgets( string, 1000, sublist ) )
     {
       num_iids++;
     }
   rewind( sublist );
   tig_iids       = CreateScalarHashTable_AS( num_iids );
   tig_iids_found = CreateScalarHashTable_AS( num_iids );
   if( tig_iids == NULL || tig_iids_found == NULL ) return 1;
   for( this_id = 0; this_id < num_iids; this_id++ )
     {
       fgets( string, 1000, sublist );
       InsertInHashTable_AS(tig_iids, STR_TO_UID(string, NULL, 10), 0, 0, 0);
     }
  
   fclose( sublist );
 }

 while (ReadProtoMesg_AS(stdin,&pmesg) != EOF){
   if (pmesg->t ==MESG_ICM)  {
     contig = (IntConConMesg *) pmesg->m;
     if(contig->length < min_len) continue;
     if( do_all || ExistsInHashTable_AS(tig_iids, contig->iaccession, 0)) {
       if ( ! do_all )
         InsertInHashTable_AS(tig_iids_found, contig->iaccession, 0, 0, 0);

       ma = CreateMultiAlignTFromICM(contig, contig->iaccession,  0);
        
       if (contig->placed == AS_PLACED) {
	 CollectStats(ma, frag_store, pcs, pfs,AS_READ_CLEAR_LATEST);
       } 
         
       //      PrintMultiAlignT(out,ma,frag_store,0,0);
       if(out!=NULL)fflush(out);
       if ( ! do_all && out != NULL) {
	 WriteProtoMesg_AS(out,pmesg);
       }
       if(out!=NULL)fflush(out);
       DeleteMultiAlignT(ma);
     }
   }
   if (pmesg->t ==MESG_ISF)  {
     break;
   }
 }
 fclose(pcs);
 fclose(pfs);
 if(out!=NULL)fclose(out);
 exit (0);
}