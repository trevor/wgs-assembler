
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
/* SMAcount 
 *    scans all AFG messages and orints out their id and the number of screenmatches they have
 *
 * $Id: SMAcount.c,v 1.3 2005-03-22 19:08:38 jason_miller Exp $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "AS_global.h"
#include <assert.h>

int main(void)
{ 
  GenericMesg *pmesg;
  MesgReader reader = InputFileType_AS(stdin);
  int smatches = 0;

 while (reader(stdin,&pmesg) != EOF){
   if( pmesg->t == MESG_AFG )
     {
       AugFragMesg *iaf = (AugFragMesg*) pmesg->m;
       ScreenMatch *sm = iaf->screened;
       if( sm != NULL)
	 fprintf(stderr,"\nFragment " F_UID "\n",iaf->eaccession);
       
       while( sm != NULL )
	 {
	   fprintf(stderr,
                   "Screenmatch where (" F_COORD "," F_COORD ") what " F_UID " rep. id " F_UID " \n",
                   sm->where.bgn,sm->where.end,sm->what,sm->repeat_id);
	   sm = sm->next;
	   smatches++;
	 }	
 
     }

 }
 fprintf(stderr,"\n Number of Screen matches %d\n",smatches); 
 exit (0);
}
