
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
/* $Id: AS_MSG_pmesg_test.c,v 1.3 2005-03-22 19:06:17 jason_miller Exp $ */

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h> /* man 3 getopt */

#include "AS_global.h"

int main(int argc, char * argv [])
{
  int illegal = 0;
  int timer_mode = 0;
  int no_output_mode = 0;
  int binary_output_mode = 0;
  int duplicate_mode = 0;

  /**************** Process Command Line Arguments *********************/
  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv, "bdnt")) != EOF))
      switch(ch) {
      case 't':
	timer_mode = 1;
	break;
      case 'n':
	no_output_mode  = 1;
	break;
      case 'd':
	duplicate_mode  = 1;
	break;
      case 'b':
	binary_output_mode  = 1;
	break;
      case '?':
	fprintf(stderr,"Unrecognized option -%c",optopt);
      default :
	errflg++;
      }
    if((illegal == 1) || (argc - optind > 0 ))
      {
	fprintf (stderr, "USAGE: %s "
		 "[-t] "
		 "[-b] "
		 "[-n] "
		 "[-d] "
		 "< InputFileName >OutputFileName \n"
		 " -t : produce timing information\n"
		 " -b : produce binary output\n"
		 " -n : send output to /dev/null\n"
		 " -d : internally duplicate each message\n" ,
		 argv[0]);
	exit (EXIT_FAILURE);
      }
  }

  {
    /* MessageType  imesgtype; */
    GenericMesg *pmesg = NULL, *wmesg = NULL;
    MesgReader   reader;
    time_t tp1,tp2;

    if(timer_mode) {
      time(&tp1); fprintf(stderr,"Begin timing\n");
    }

    reader = InputFileType_AS(stdin);

    while (reader(stdin,&pmesg) != EOF)
      {
	if(duplicate_mode) {
	#if 0
	  wmesg = DuplicateProtoMesg_AS(pmesg);
	#else
	  fprintf(stderr,"Duplicate no longer supported\n");
	#endif
	} else {
	  wmesg = pmesg;
	}
	if(! no_output_mode) {
	  if(binary_output_mode) {
	    WriteBinaryMesg_AS(stdout,wmesg);
	  } else {
	    WriteProtoMesg_AS(stdout,wmesg);
	  }
	}
	if(duplicate_mode) {
#ifdef TODELTE
	  FreeProtoMesg_AS(wmesg);
#endif
	}
      }
    if(timer_mode) {
      time(&tp2);
      fprintf(stderr,"%10" F_TIME_TP " sec: Finished timing\n", (tp2-tp1));
    }
  }
  exit (0);
}


