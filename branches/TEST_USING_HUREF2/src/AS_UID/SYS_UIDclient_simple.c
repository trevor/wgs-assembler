
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
/**********************************************************************
 Module:      SYS_UID
 Description: This file contains functions to allocate
              UIDs from the UID server.
 Assumptions:
**********************************************************************/

static char CM_ID[] = "$Id: SYS_UIDclient_simple.c,v 1.2 2005-09-30 20:15:01 eliv Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "SYS_UIDcommon.h"
#include "SYS_UIDclient.h"



static uint64 SYS_UID_uidStart = 99000000000LL;


int32 get_uids(uint64 blockSize, uint64 *interval, int32 real)
{
  /******************************************************************/
  // Allocates blockSize many UIDs from the UID server if real
  // is TRUE. Otherwise it allocates some dummy numbers.
  /*****************************************************************/
  int32  uidStatus;
  uint64 maxBlockSize;

  if( ! real )
    {
      interval[0] = 4711;
      interval[1] = blockSize;
      interval[2] = 4711+2*blockSize;
      interval[3] = blockSize;
    }
  else
    {

      /* First check whether the UID server can accomodate for our buffer */
      uidStatus = SYS_UIDgetMaxUIDSize(&maxBlockSize);
      
      if( uidStatus != UID_CODE_OK || maxBlockSize < blockSize ) {
        fprintf(stderr, "UID blocksize query failed\n");
        assert(0);
      }

      /* Now set the blocksize of th UID server appropriately */
      SYS_UIDsetUIDSize(blockSize);

      /* Finally, get the actual interval */
      uidStatus = SYS_UIDgetNewUIDInterval(interval);
      if( uidStatus != UID_CODE_OK ) {
        fprintf(stderr, "UID query failed failed\n");
        assert(1);
      }
    }  
  return UID_CODE_OK;
}


int32 get_next_uid(uint64 *uid, int32 real){
  if( real == FALSE ){
    *uid = SYS_UID_uidStart++;
    return UID_CODE_OK;
  }
  else
    return SYS_UIDgetNextUID(uid);
}







/*--------------------------------------------------------------------*/
/*  MISC routines */
/*--------------------------------------------------------------------*/

char *UidServer = "SYS_UID_SERVER_HOSTNAME=notcurrentlyset";
char *UidPort   = "SYS_UID_SERVER_PORT=5001";
char *UidFailServer = "SYS_UID_FAILSAFE_SERVER_HOSTNAME=notcurrentlyseteither";
char *UidFailPort   = "SYS_UID_FAILSAFE_SERVER_PORT=5002";

void check_environment(){

#ifdef USE_SOAP_UID
  return;
#else

  /******************************************************************/
  // Checks whether the SYS_UID* variables are set.
  // If not it exits proposing a standard value
  /*****************************************************************/

  char *currentUidServer;
  char *currentUidPort;
  char *currentUidFailServer;
  char *currentUidFailPort;

  if( (currentUidPort       = getenv("SYS_UID_SERVER_HOSTNAME"))          == NULL || 
      (currentUidServer     = getenv("SYS_UID_SERVER_PORT"))              == NULL ||
      (currentUidFailPort   = getenv("SYS_UID_FAILSAFE_SERVER_HOSTNAME")) == NULL ||
      (currentUidFailServer = getenv("SYS_UID_FAILSAFE_SERVER_PORT"))     == NULL )
    {
      char dummy[1000];
      sprintf(dummy,"UID Server environment variables not set!\nTrying \n%s\n%s\n%s\n%s\n",UidServer,UidPort,UidFailServer,UidFailPort);
      fprintf(stderr,dummy);
      putenv(UidServer);
      putenv(UidPort);
      putenv(UidFailServer);
      putenv(UidFailPort);
    }
#endif
}

void set_start_uid(uint64 s) {
  SYS_UID_uidStart = s;
}
uint64 get_start_uid(void) {
  return(SYS_UID_uidStart);
}


