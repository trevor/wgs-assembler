
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

static char CM_ID[] = "$Id: createFrgDeletes.c,v 1.5 2006-10-08 08:47:39 brianwalenz Exp $";


/*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>


#include "AS_global.h"

#include "AS_PER_gkpStore.h"

#include "SYS_UIDcommon.h"
#include "SYS_UIDclient.h"

#include "MultiAlignment_CNS.h"

#define MAXSEQLEN 20000

/* Output text field item with 3-code field-name "tag". */

int main( int argc, char *argv[])
{
  char *inputPath;
  char *prefix;

  int setIIDLIST = FALSE;
  int setGatekeeperStore = FALSE;
  int fragIID,mateIID;
  FILE *iidlist;
  Fragment_ID fragUID,mateUID;
  char iidlist_name[2000];
  char GKP_Store_Name[2000];
  GateKeeperStore gkpStore;
  GateKeeperFragmentRecord gkpFrag,gkpMate;
  uint64 uid, mateuid;
  ReadStructp fsread=new_ReadStruct();
  ReadStructp fsmate=new_ReadStruct();
  int realUID=1;
  int UIDstart=1230000;
  int firstUID=1;
  CDS_UID_t       interval_UID[4];
  Overlap *ovl;
  IntUnitigMesg ium;
  IntMultiPos the_imps[2];
  CDS_UID_t mergeUid;
  char seq[MAXSEQLEN], qlt[MAXSEQLEN];
  int clr_bgn,clr_end;
  int iid;
  VA_TYPE(int32) *deltas=CreateVA_int32(1);
  VA_TYPE(char) *sequence=CreateVA_char(200000);
  VA_TYPE(char) *quality=CreateVA_char(200000);

  //  setbuf(stdout,NULL);

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,
				    "i:g:")) != EOF)){
      switch(ch) {
        case 'i':
          strcpy( iidlist_name, argv[optind - 1]);
          setIIDLIST = TRUE;
          break;
        case 'g':
          strcpy( GKP_Store_Name, argv[optind - 1]);
          setGatekeeperStore = TRUE;
          break;	  
        default :
          errflg++;
      }
    }

    if( setIIDLIST==0 || (setGatekeeperStore == 0) || errflg>0)
      {
	fprintf(stderr,"* argc = %d optind = %d setIIDLIST = %d setGatekeeperStore = %d\n",
		argc, optind, setIIDLIST,setGatekeeperStore);
	fprintf (stderr, "USAGE:  %s -i <file of UIDs> -g <GatekeeperStoreName>\n",argv[0]);
	exit (EXIT_FAILURE);
      }

  }

  InitGateKeeperStore(&gkpStore,GKP_Store_Name);
  assert(TestOpenGateKeeperStore(&gkpStore) == TRUE);
  OpenReadOnlyGateKeeperStore(&gkpStore);

  iidlist = fopen(iidlist_name,"r");

  /*************************/
  // Construct a BAT message
  /*************************/
  {

    /*************************/
    // Get a UID to use
    /*************************/
    {
      int32 blockSize = 300;
      int32  uidStatus;
      CDS_UID_t interval_UID[4];
      if(firstUID){
	firstUID=0;
	set_start_uid(UIDstart); /* used if readUID == FALSE */
	get_uids(blockSize,interval_UID,realUID);
      }

      uidStatus = get_next_uid(&mergeUid,realUID);
      if( uidStatus != UID_CODE_OK )
	{
	  get_uids(blockSize,interval_UID,realUID);
	  uidStatus = get_next_uid(&mergeUid,realUID);
	}	  
      if( UID_CODE_OK != uidStatus )
	{ 
          fprintf(stderr, "Could not get UID \n");
          assert(0);
	}
    }
    /***********************/
    // Print a BAT message
    /***********************/
    printf("{BAT\n");
    printf("bna:(Batch name)\n");
    printf("crt:" F_TIME_T "\n",time(NULL));
    printf("acc:" F_UID "\n",mergeUid);
    printf("com:\nCreated by %s\n.\n",__FILE__);
    printf("}\n");
  }

  /*************************/
  // over all fragments in list
  /*************************/
  
  while(fscanf(iidlist,F_S32,&fragIID)==1){

    int rv1,rv2;

    /*************************/
    // get the fragment
    /*************************/

    //fprintf(stderr,"Working on frgIID %d\n",fragIID);
    rv1 = getGateKeeperFragmentStore(gkpStore.frgStore,fragIID,&gkpFrag);

    assert(rv1==0);
    fragUID = gkpFrag.readUID;

    /*************************/
    // check for an appropriate mate
    /*************************/

    if(gkpFrag.numLinks>0){
      GateKeeperLinkRecordIterator iterator;
      GateKeeperLinkRecord link;
      CreateGateKeeperLinkRecordIterator(gkpStore.lnkStore, gkpFrag.linkHead,fragIID, &iterator);
      while(NextGateKeeperLinkRecordIterator(&iterator, &link))
	mateIID = (link.frag1 == fragIID) ? link.frag2 : link.frag1;
      //      if(mateIID>fragIID){
      {
	rv2 = getGateKeeperFragmentStore(gkpStore.frgStore,mateIID,&gkpFrag);
	assert(rv2==0);
	mateUID = gkpFrag.readUID;
	printf("{LKG\n");
	printf("act:D\n");
	printf("typ:M\n");
	printf("fg1:" F_UID "\n",fragUID);
	printf("fg2:" F_UID "\n",mateUID);
	printf("}\n");
      }
    }

    printf("{FRG\n");
    printf("act:D\n");
    printf("acc:" F_UID "\n",fragUID);
    printf("}\n");

  }
  exit(0);
}