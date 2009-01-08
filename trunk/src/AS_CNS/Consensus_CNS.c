
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

const char *mainid = "$Id: Consensus_CNS.c,v 1.73 2009-01-08 21:12:06 brianwalenz Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_fileIO.h"
#include "AS_MSG_pmesg.h"
#include "AS_PER_gkpStore.h"
#include "AS_SDB_SequenceDB.h"
#include "MultiAlignStore_CNS.h"
#include "MultiAlignment_CNS.h"
#include "MultiAlignment_CNS_private.h"
#include "AS_ALN_aligners.h"

#define MAX_NUM_UNITIG_FAILURES 100
#define MAX_NUM_CONTIG_FAILURES 100


extern int    numScores;
extern double lScoreAve;
extern double aScoreAve;
extern double bScoreAve;

static void
writeFailure(char *outName, GenericMesg *pmesg) {
  FILE *fout;
  char  fname[1000];

  sprintf(fname, "%s.failed", outName);

  errno = 0;
  fout = fopen(fname, "a");
  if (errno) {
    fprintf(stderr, "Failed to open '%s' for storing the failed messge: %s\n", fname, strerror(errno));
    fprintf(stderr, "------------------------------------------------------------\n");
    WriteProtoMesg_AS(stderr, pmesg);
    fprintf(stderr, "------------------------------------------------------------\n");
  } else {
    WriteProtoMesg_AS(fout,pmesg); // pass through the Unitig message and continue
    fclose(fout);
  }
}



int
main (int argc, char **argv) {

  char  *inpName = NULL;
  char   tmpName[FILENAME_MAX] = {0};
  char  *outName = NULL;
  char  *gkpName = NULL;

  FILE   *cnsinp = NULL;
  FILE   *cnsout = NULL;

  char  *sdbName = NULL;
  int    sdbVers = -1;
  int    sdbPart = -1;

  int    gkpPart     = 0;
  int    gkpInMemory = 0;

  int    extract     = -1;

  int    allow_neg_hang_retry = 0;

  int    saveUnitigMultiAlign = 0;

  CNS_Options options = { CNS_OPTIONS_SPLIT_ALLELES_DEFAULT,
                          CNS_OPTIONS_MIN_ANCHOR_DEFAULT };

  int num_unitig_failures = 0;
  int num_contig_failures = 0;
  int num_contig_skips    = 0;


  CNS_PrintKey printwhat=CNS_STATS_ONLY;

  //  Comminucate to MultiAlignment_CNS.c that we are doing consensus and not cgw.
  thisIsConsensus = 1;


  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      allow_neg_hang = 1;

    } else if (strcmp(argv[arg], "-G") == 0) {
      allow_neg_hang_retry = 1;

    } else if (strcmp(argv[arg], "-K") == 0) {
      options.split_alleles = 0;

    } else if (strcmp(argv[arg], "-v") == 0) {
      int what = atoi(argv[++arg]);
      switch (what) {
        case 0:
          fprintf(stderr, "Verbose mode disabled.\n");
          break;
        case 1:
          printwhat = CNS_CONSENSUS;
          break;
        case 2:
          printwhat = CNS_DOTS;
          break;
        case 3:
          printwhat = CNS_NODOTS;
          break;
        case 4:
          printwhat = CNS_VIEW_UNITIG;
          break;
        case 5:
          printwhat = CNS_VERBOSE;
          break;
        default:
          fprintf(stderr, "Unknown verbose mode %d\n", what);
          err++;
          break;
      }

    } else if (strcmp(argv[arg], "-o") == 0) {
      outName = argv[++arg];

    } else if (strcmp(argv[arg], "-S") == 0) {
      gkpPart = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-w") == 0) {
      options.smooth_win = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-M") == 0) {
      saveUnitigMultiAlign = 1;

    } else if (strcmp(argv[arg], "-m") == 0) {
      gkpInMemory = 1;

    } else if (strcmp(argv[arg], "-s") == 0) {
      USE_SDB = 1;
      sdbName = argv[++arg];

    } else if (strcmp(argv[arg], "-p") == 0) {
      USE_SDB = 1;
      sdbPart = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-V") == 0) {
      USE_SDB = 1;
      sdbVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-D") == 0) {
      ++arg;
      if        (strcmp(argv[arg], "dumpunitigs") == 0) {
        DUMP_UNITIGS_IN_MULTIALIGNCONTIG = 1;
      } else if (strcmp(argv[arg], "verbosemultialign") == 0) {
        VERBOSE_MULTIALIGN_OUTPUT = 1;
      } else if (strcmp(argv[arg], "forceunitigabut") == 0) {
        FORCE_UNITIG_ABUT = 1;
      } else {
        fprintf(stderr, "Unrecognized option '%s' to -D.\n", argv[arg]);
        err++;
      }

    } else if (strcmp(argv[arg], "-e") == 0) {
      arg++;
      if (argv[arg][0] == '#') {
        extract = atoi(argv[arg] + 1);
      } else {
        fprintf(stderr, "error: form: '-s #498234'\n");
      }

    } else if (strcmp(argv[arg], "-U") == 0) {
      clear_range_to_use = AS_READ_CLEAR_OBT;

    } else {
      if (argv[arg][0] == '-') {
        fprintf(stderr, "Unrecognized option %s", argv[arg]);
        err++;
      } else if (gkpName == NULL) {
        gkpName = argv[arg];
      } else if (inpName == NULL) {
        inpName = argv[arg];
      } else {
        fprintf(stderr, "Unrecognized option %s", argv[arg]);
        err++;
      }
    }
    arg++;
  }
  if (USE_SDB && ((sdbName == NULL) || ((sdbVers == -1) || (sdbPart == -1)))) {
    fprintf(stderr, "ERROR:  Unsupported!  SDB must be partitioned, version must be supplied.\n");
    err++;
  }
  if (err || (inpName == NULL) || (gkpName == NULL)) {
    fprintf(stderr, "usage: %s [opts] gkpStore input-messages\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "    -v [0-4]     Verbose:  0 = verbose off\n");
    fprintf(stderr, "                           1 = horizontal multi-alignment print in .clg\n");
    fprintf(stderr, "                           2 = 'dots'     multi-alignment print in .clg\n");
    fprintf(stderr, "                           3 = like 2, but dots are replaced with whitespace\n");
    fprintf(stderr, "                           4 = like 1, but with unitigs in  multi-alignment print in .clg\n");
    fprintf(stderr, "    -K           don't split alleles when calling consensus\n");
    fprintf(stderr, "    -w win_size  specify the size of the 'smoothing window' that will be used in consensus calling\n");
    fprintf(stderr, "                 If two SNPs are located win_size or less bases apart one from another,\n");
    fprintf(stderr, "                 then they will be treated as one block\n");
    fprintf(stderr, "    -S partition Use gkpStorePartition partition, loaded into memory\n");
    fprintf(stderr, "    -m           Load gkpStorePartition into memory (default reads from disk)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -D opt       Enable debugging option 'opt'.  One of 'dumpunitigs', 'verbosemultialign',\n");
    fprintf(stderr, "                    and 'forceunitigabut'.  (-X not needed).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -e #%%d      Extract only a single ICM/IUM by internal id\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  gkpStore       path to previously created Fragment Store\n");
    fprintf(stderr, "  input-messages previously created .cgw/.cgb file\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Output:\n");
    fprintf(stderr, "   -o <filename>     output filename\n");
    exit(1);
  }

  gkpStore = openGateKeeperStore(gkpName, FALSE);

  if (USE_SDB) {
    sequenceDB = openSequenceDB(sdbName, FALSE, sdbVers);
    openSequenceDBPartition(sequenceDB, sdbPart);
  } else {
    unitigStore = CreateMultiAlignStoreT();
  }

  if (gkpPart)
    loadGateKeeperPartition(gkpStore, gkpPart);
  else if (gkpInMemory)
    loadGateKeeperStorePartial(gkpStore, 0, 0, FRAG_S_QLT);

  //
  //  INPUT and OUTPUT
  //

  //  Remove any existing failed file.  If we don't, and we are doing
  //  a rerun, we still look like there is a failure.
  //
  sprintf(tmpName, "%s.failed", outName);
  unlink(tmpName);
  unlink(outName);

  sprintf(tmpName, "%s_tmp", outName);

  errno = 0;
  cnsinp = fopen(inpName, "r");
  if (errno) {
    fprintf(stderr, "Could not open '%s' for input: %s\n", inpName, strerror(errno));
    exit(1);
  }

  cnsout = fopen(tmpName, "w");
  if (errno) {
    fprintf(stderr, "Could not open '%s' for output: %s\n", tmpName, strerror(errno));
    exit(1);
  }

  VA_TYPE(int32) *deltas   = CreateVA_int32(1);
  VA_TYPE(char)  *sequence = CreateVA_char(200000);
  VA_TYPE(char)  *quality  = CreateVA_char(200000);

  GenericMesg   *pmesg = NULL;
  while ((ReadProtoMesg_AS(cnsinp,&pmesg) != EOF)) {

    if (pmesg->t == MESG_IUM) {
      IntUnitigMesg *iunitig = (IntUnitigMesg *)(pmesg->m);
      int            unitigfail = 0;

      if (extract > -1 && iunitig->iaccession != extract)
        break;

      if (VERBOSE_MULTIALIGN_OUTPUT)
        fprintf(stderr, "MultiAlignUnitig %d %f pieces/length\n",
                iunitig->iaccession,
                (double)iunitig->num_frags / iunitig->length);

      unitigfail = MultiAlignUnitig(iunitig, gkpStore, sequence, quality, deltas, printwhat, &options);

      if ((unitigfail == EXIT_FAILURE) &&
          (allow_neg_hang_retry) &&
          (allow_neg_hang == 0)) {
        fprintf(stderr, "MultiAlignUnitig failed for unitig %d -- try again with negative hangs allowed\n", iunitig->iaccession);
        allow_neg_hang = 1;
        unitigfail = MultiAlignUnitig(iunitig, gkpStore, sequence, quality, deltas, printwhat, &options);
        allow_neg_hang = 0;
        if (unitigfail != EXIT_FAILURE)
          NumUnitigRetrySuccess++;
      }

      if (unitigfail == EXIT_FAILURE) {
        num_unitig_failures++;
        if (num_unitig_failures <= MAX_NUM_UNITIG_FAILURES) {
          fprintf(stderr, "MultiAlignUnitig failed for unitig %d\n", iunitig->iaccession);
          writeFailure(outName, pmesg);
        } else {
          fprintf(stderr, "MultiAlignUnitig failed more than MAX_NUM_UNITIG_FAILURES times.  Fail.");
          exit(1);
        }
      }

      if (unitigfail == EXIT_SUCCESS) {
        if (saveUnitigMultiAlign)
          SetMultiAlignInStore(unitigStore,
                               iunitig->iaccession,
                               CreateMultiAlignTFromIUM(iunitig, -1, FALSE));

        pmesg->t = MESG_IUM;
        pmesg->m = iunitig;
        WriteProtoMesg_AS(cnsout, pmesg);
      }

    } else if (pmesg->t == MESG_ICM) {
      IntConConMesg *pcontig = (IntConConMesg *)(pmesg->m);
      int            contigFail = 0;

      if (extract > -1 && pcontig->iaccession != extract)
        break;

      if (VERBOSE_MULTIALIGN_OUTPUT)
        fprintf(stderr, "MultiAlignContig %d %f pieces/length\n",
                pcontig->iaccession,
                (double)pcontig->num_pieces / pcontig->length);

      pcontig->num_vars == 0;
      pcontig->v_list == NULL;

      if ((((double)pcontig->num_pieces / pcontig->length) > 1.0) &&
          (pcontig->consensus != NULL)) {
        fprintf(stderr, "WARNING:  MultiAlignContig skipping contig %d -- %f pieces/length (and it already has a consensus sequence)\n",
                pcontig->iaccession,
                (double)pcontig->num_pieces / pcontig->length);
        num_contig_skips++;
        contigFail = EXIT_SUCCESS;
      } else {
        contigFail = MultiAlignContig(pcontig, sequence, quality, deltas, printwhat, &options);
      }

      if (contigFail == EXIT_FAILURE) {
        num_contig_failures++;

        if (num_contig_failures <= MAX_NUM_CONTIG_FAILURES) {
          fprintf(stderr, "MultiAlignContig failed for contig %d\n", pcontig->iaccession);
          writeFailure(outName, pmesg);
        } else {
          fprintf(stderr, "MultiAlignContig failed more than MAX_NUM_CONTIG_FAILURES times.  Fail.");
          exit(1);
        }
      }

      if (contigFail == EXIT_SUCCESS) {
        if ( printwhat == CNS_CONSENSUS && pcontig->num_pieces > 0) {
          MultiAlignT   *ma = NULL;
          ma = CreateMultiAlignTFromICM(pcontig,-1,0);
          PrintMultiAlignT(stderr, ma, gkpStore, 1, 0, clear_range_to_use);
          DeleteMultiAlignT(ma);
        }

        pmesg->t = MESG_ICM;
        pmesg->m = pcontig;
        WriteProtoMesg_AS(cnsout, pmesg);

        if (pcontig->v_list != NULL) {
          int i;
          for (i=0; i<pcontig->num_vars; i++) {
            safe_free(pcontig->v_list[i].nr_conf_alleles);
            safe_free(pcontig->v_list[i].weights);
            safe_free(pcontig->v_list[i].var_seq);
            safe_free(pcontig->v_list[i].conf_read_iids);
          }
          safe_free(pcontig->v_list);
        }
        pcontig->num_vars = 0;
      }
    } else {
      //  Pass through messages that aren't IUM or ICM.
      WriteProtoMesg_AS(cnsout, pmesg);
    }
  }


 goodbye:

  if (unitigStore)
    DeleteMultiAlignStoreT(unitigStore);

  fclose(cnsout);

  fprintf(stderr, "\n");
  fprintf(stderr, "NumColumnsInUnitigs             = %d\n", NumColumnsInUnitigs);
  fprintf(stderr, "NumGapsInUnitigs                = %d\n", NumGapsInUnitigs);
  fprintf(stderr, "NumRunsOfGapsInUnitigReads      = %d\n", NumRunsOfGapsInUnitigReads);
  fprintf(stderr, "NumColumnsInContigs             = %d\n", NumColumnsInContigs);
  fprintf(stderr, "NumGapsInContigs                = %d\n", NumGapsInContigs);
  fprintf(stderr, "NumRunsOfGapsInContigReads      = %d\n", NumRunsOfGapsInContigReads);
  fprintf(stderr, "NumAAMismatches                 = %d\n", NumAAMismatches);
  fprintf(stderr, "NumVARRecords                   = %d\n", NumVARRecords);
  fprintf(stderr, "NumVARStringsWithFlankingGaps   = %d\n", NumVARStringsWithFlankingGaps);
  fprintf(stderr, "NumUnitigRetrySuccess           = %d\n", NumUnitigRetrySuccess);
  fprintf(stderr, "\n");

  errno = 0;
  rename(tmpName, outName);
  if (errno) {
    fprintf(stderr, "ERROR!  Failed to rename output '%s' to '%s': %s\n",
            tmpName, outName, strerror(errno));
    return(1);
  }

  if (num_contig_skips) {
    fprintf(stderr, "WARNING:  Total number of skipped contigs = %d\n", num_contig_skips);
  }

  if (num_unitig_failures || num_contig_failures) {
    fprintf(stderr, "WARNING:  Total number of unitig failures = %d\n", num_unitig_failures);
    fprintf(stderr, "WARNING:  Total number of contig failures = %d\n", num_contig_failures);
    fprintf(stderr, "\n");
    fprintf(stderr, "Consensus did NOT finish successfully.\n");
    return(1);
  }

  if (numScores > 0) {
    fprintf(stderr, "Average Scores (%d samples):\n", numScores);
    fprintf(stderr, "  lScore = %f\n", lScoreAve / numScores);
    fprintf(stderr, "  aScore = %f\n", aScoreAve / numScores);
    fprintf(stderr, "  bScore = %f\n", bScoreAve / numScores);
  }

  fprintf(stderr, "Consensus finished successfully.  Bye.\n");
  return(0);
}
