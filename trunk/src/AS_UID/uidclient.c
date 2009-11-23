
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute.
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

static const char *rcsid = "$Id: uidclient.c,v 1.6 2009-11-23 00:31:38 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "SYS_UIDclient.h"

int
main(int argc, char **argv) {
  int     blockSize = 1;
  int     numUIDs   = 1;
  int     doThrash  = 0;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-p") == 0) {
      numUIDs = atoi(argv[++arg]);
      if (numUIDs < blockSize)
        blockSize = numUIDs;

    } else if (strcmp(argv[arg], "-n") == 0) {
      //SYS_UIDset_euid_namespace(argv[++arg]);

    } else if (strcmp(argv[arg], "-E") == 0) {
      //SYS_UIDset_euid_server(argv[++arg]);

    } else if (strcmp(argv[arg], "-thrash") == 0) {
      doThrash = 1;

    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if ((err) || (numUIDs == 0)) {
    fprintf(stderr, "%s: [-p n] [-n ns] [-E server] [-d ms]\n", argv[0]);
    fprintf(stderr, "  -p n       print n UIDs and exit.\n");
    fprintf(stderr, "  -n ns      use namespace ns.\n");
    fprintf(stderr, "  -E server  contact EUID server 'server'.\n");
    fprintf(stderr, "  -thrash    debug; get UIDs as fast as possible using blocksize 1.\n");
    fprintf(stderr, "             This is not what you want.  Don't use it.\n");
    exit(1);
  }

  UIDserver  *uids = UIDserverInitialize(blockSize, 0);

  while (doThrash)
    getUID(uids);  //  Forever or never.

  while (numUIDs > 0) {
    fprintf(stdout, F_U64"\n", getUID(uids));
    numUIDs--;
  }

  exit(0);
}
