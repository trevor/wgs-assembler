
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
/* MODULE FOR READING FASTA SEQUENCES, COMPUTING OVERLAPS, AND THE COLUMN
   SETS FOR THE CORRELATED-DIFFERENCES DETECTOR.
*/

#undef INPUT

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_UTL_Var.h"
#include "AS_ALN_aligners.h"
#include "AS_CGB_all.h"


/* Get_sequence gets the next FASTA formated sequence from input.  It
   assumes that it and only it is doing input and is designed to handle
   lines and strings of arbitrary length.  It allocates memory for the
   sequence and returns a pointer to it.  */

#define LBUFLEN 512

char *get_sequence(FILE *input)
{ static char *seqbuf, linebuf[LBUFLEN];
  static int   first = 1;
  static int   top, nei;

  register char *newbuf;
  register size_t l;
  register int e, bol, beg;

  if (first)
    { first  = 0;
      top    = 2048;
      seqbuf = (char *) ckalloc(sizeof(char)*top);
      if (fgets(linebuf,LBUFLEN,input) == NULL) return (NULL);
      if (*linebuf != '>')
        { fprintf(stderr,"First line must start with an >-sign\n");
          exit (1);
        }
    }
  else
    { if (!nei) return (NULL); }

  do
    { l = strlen(linebuf);
      if (linebuf[l-1] == '\n') break;
    }
  while (fgets(linebuf,LBUFLEN,input) != NULL);

  bol = 1;
  beg = 1;
  e   = 0;
  while((nei = (fgets(linebuf,LBUFLEN,input) != NULL)) != 0)
    { if (bol && *linebuf == '>')
        if (beg)
          { do
              { l = strlen(linebuf);
                if (linebuf[l-1] == '\n') break;
              }
            while (fgets(linebuf,LBUFLEN,input) != NULL);
          }
        else
          break;
      else
        { l = strlen(linebuf);
          if (e + l >= top)
            { top = (int) (1.5*(e+l) + 200);
              newbuf = (char *) ckalloc(sizeof(char)*top);
              seqbuf[e] = '\0';
              strcpy(newbuf,seqbuf);
              free(seqbuf);
              seqbuf = newbuf;
            }
          strcpy(seqbuf+e,linebuf);
          bol = (linebuf[l-1] == '\n');
          e = (e+l) - bol;
          beg = 0;
        }
    }
  seqbuf[e] = '\0';

  newbuf = (char *) ckalloc(sizeof(char)*(e+1));
  strcpy(newbuf,seqbuf);
  
  return (newbuf);
}

/* Get_sequences gets all the FASTA formatted sequences from input, where
   it is assuming the input is a series of such sequences.  It sets *nseq
   to the number of sequences read less one, allocates space for each
   sequence, and returns an allocated array, seqa[0..k], of pointers to
   the sequences.
*/

char **get_sequences(FILE *input, int *nseq)
{ int    max, k;
  char **seqa, **seqn;

  max  = 32;
  seqa = (char **) ckalloc(max*sizeof(char *));

  k = 0;
  while (1)
    { for (; k < max; k++)
        { seqa[k] = get_sequence(input);
          if (seqa[k] == NULL) break;
        }
      if (k < max) break;
      seqn = (char **) ckalloc(2*max*sizeof(char *));
      for (k = 0; k < max; k++)
        seqn[k] = seqa[k];
      free(seqa);
      seqa = seqn;
      max *= 2;
    }

  *nseq = k-1;
  return (seqa);
}

/* Write seq on the standard output, 50 symbols to a line. */

void show_sequence(char *seq)
{ size_t len, i;

  len = strlen(seq);
  for (i = 0; i < len; i += 50)
    if (i+50 < len)
      printf("%.50s\n",seq+i);
    else
      printf("%s\n",seq+i);
}

int main(int argc, char *argv[])
{ int    K;
  char **Seqs;
  InternalFragMesg  A, B;
  OverlapMesg  *O;
  FILE *OVLFile=NULL;
  MesgWriter WriteMesg_AS = NULL;
  int ori;
  int abnd;
  int bbnd;
  int minlen=40;
  int first=1;
  double err=.06;

  Seqs = get_sequences(stdin,&K);

  if(argc>2&&strcmp(argv[1],"-P")==0){
    fprintf(stderr,"Printing OVLs to %s\n",argv[2]);
    OVLFile=fopen(argv[2],"w");
    assert(OVLFile!=NULL);
    WriteMesg_AS = OutputFileType_AS(AS_PROTO_OUTPUT);
    assert(WriteMesg_AS!=NULL);
  }

fprintf(stderr,"Read in %d sequences\n",K+1);

#ifdef INPUT
  { int i;

    printf("\nThe Sequences %d:\n\n",K+1);
    for (i = 0; i <= K; i++)
      { printf("> %d\n",i+1);
        show_sequence(Seqs[i]);
      }
  }
#endif


  { 
    int i, j, where=1;
    int olaps, tlaps;

    A.quality = NULL;
    B.quality = NULL;
    tlaps = olaps = 0;
    for (j = 0; j < K; j++){
      for (i = j+1; i <= K; i++){
	for(ori=0;ori<=1;ori++){
	  tlaps += 1;
	  A.sequence = Seqs[j];
	  B.sequence = Seqs[i];
	  A.iaccession = A.eaccession = j+1;
	  B.iaccession = B.eaccession = i+1;

	  abnd=strlen(A.sequence);
	  bbnd=-strlen(B.sequence);


	  if(argc>2){
	    if(strcmp(argv[1],"-P")!=0|| argc>=6){
	      assert(argc>=4);
	      err=atof(argv[argc-3]);
	      abnd=atoi(argv[argc-1]);
	      bbnd=atoi(argv[argc-2]);
	      if((strcmp(argv[1],"-P")!=0&&argc>4)||argc>6){
		minlen=atoi(argv[argc-4]);
	      }
	    }
	    if(first){
	      first=0;
	      fprintf(stderr,"Calling Local_Overlap_AS with band %d to %d,erate %f, minlen %d\n",
		      bbnd,abnd,err,minlen);
	    }

	  }


	  O = DP_Compare_AS(&A,&B,bbnd,abnd,
			    ori,err,1e-6,minlen,AS_FIND_ALIGN,&where);
          if (O != NULL){
            olaps += 1;
	    Print_Overlap_AS(stdout,&A,&B,O);
	    printf("Overlap quality: %f\n",O->quality);

	    { int del, sub, ins, affdel, affins, alen, blen, blockdel, blockins;
	      float errRate, errRateAffine;
	      
#define AFFINEBLOCKSIZE 4
	      Analyze_Affine_Overlap_AS(&A,&B,O,AS_ANALYZE_ALL,&alen,&blen,&del,&sub,&ins,
					&affdel,&affins,&blockdel,&blockins,AFFINEBLOCKSIZE);
	      
	      errRate = (sub+ins+del)/(double)(alen+ins);
	      
	      errRateAffine = (sub+affins+affdel)/
		(double)(alen+ins-(del-affdel+ins-affins));
	      
	      printf("Alen %d, Blen %d, del %d, sub %d, ins %d\n"
		     " affdel %d, affins %d, blockdel %d, blockins %d\n",
		     alen,blen,del,sub,ins,
		     affdel,affins,blockdel,blockins);
	      printf("Simple mismatch rate %f\n",errRate);
	      printf("Affine mismatch rate %f\n",errRateAffine);
	    }
	    O->min_offset=O->max_offset=O->ahg;
	    if(OVLFile!=NULL){
	      GenericMesg pmesg;
	      pmesg.m=O;
	      pmesg.t=MESG_OVL;
	      WriteMesg_AS(OVLFile,&pmesg);
	    }

	  }
	}
      }
    }
    fprintf(stderr,"Performed %d x 2 compares, found %d overlaps\n",
	    tlaps,olaps);
  }
  return (0);
}
