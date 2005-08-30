
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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "CA_ALN_local.h"
#include "AS_ALN_aligners.h"
#include "CA_ALN_scafcomp.h"

//int VarWindows[] = { 3, 4, 5, 6, 7, 8, 9, 10, 0 };
int VarWindows[] = { 3, 0 };
int MaxWindow = 0;

#define OVL_ERATE   .99
#define MIN_SEG_LEN  40
#define SEG_ERATE   .10

/* Debug conditional compilation flags */

#undef  DIAGNOSTICS
#undef  DEBUG_COMPARE
#undef  DEBUG_RAYSHOOT
#undef  DEBUG_ALIGN
#undef  DEBUG_ANALYSIS
#undef  DEBUG_LOCAL
#undef  DEBUG_OPTTAIL
#undef  DEBUG_SEGORDER
#undef  XFIG
#undef  STATS

#define ALLOW_NO_OVERLAP_INTERLEAVINGS

#define MAX_LINE_LEN   1000
#define MAX_SEQ_LEN  500000



/* Read a sequence of fast files for scaffolds and return a list of these. */

char *ScfCmp_Next_Line(FILE *ifile)
{ static char buffer[MAX_LINE_LEN+2];

  buffer[MAX_LINE_LEN+1] = '\n';
  if (fgets(buffer,MAX_SEQ_LEN+2,ifile) == NULL) return (NULL);
  if (buffer[MAX_LINE_LEN+1] != '\n')
    { fprintf(stderr,"Input line is longer than %d chars\n",MAX_LINE_LEN);
      exit (1);
    }
  return (buffer);
}

static Scaffold *Read_Next_Scaffold(FILE *ifile, int first)
{ static char seqbuf[MAX_SEQ_LEN+1];
  static char *line;

  Scaffold_Gap *gaplist;
  Scaffold_Tig *ctglist;
  Scaffold     *scaf;
  char         *packseq;
  int           numgaps;

  if (first)
    line = ScfCmp_Next_Line(ifile);
  if (line == NULL)
    return (NULL);

  if (line[0] != '>')
    { fprintf(stderr,"Header line should begin with '>'\n");
      exit (1);
    }

  { char *locate, *eptr;
    int   i;

    locate = strstr(line,"stddev:");
    if (locate == NULL)
      { fprintf(stderr,"Expecting gap variance in header line\n");
        exit (1);
      }
    locate += 7;

    numgaps = strtol(locate,&eptr,10);
    if (eptr == locate+7)
      { fprintf(stderr,"Expecting # of gaps after 'stddev:'\n");
        exit (1);
      }

    gaplist = (Scaffold_Gap *) malloc(sizeof(Scaffold_Gap)*numgaps);
    ctglist = (Scaffold_Tig *) malloc(sizeof(Scaffold_Tig)*(numgaps+1));

    locate = eptr;
    for (i = 0; i < numgaps; i++)
      { gaplist[i].gap_var = strtod(locate,&eptr);
        if (eptr == locate)
          { fprintf(stderr,"Expecting # of gaps after 'stddev:'\n");
            exit (1);
          }
        locate = eptr;
      }
  }
       
  { char *p;
    int   len;

    p = seqbuf;
    while ((line = ScfCmp_Next_Line(ifile)) != NULL)
      { if (line[0] == '>') break;
        len = strlen(line)-1;
        if ((p-seqbuf) + len > MAX_SEQ_LEN)
          { fprintf(stderr,"Sequence is longer than %d chars\n",
                           MAX_SEQ_LEN);
            exit (1);
          }
        strcpy(p,line);
        p += len;
      }
    *p = '\0';
  }

#ifdef DEBUG_READ
  { int len, i;

    len = strlen(seqbuf);
    for (i = 0; i < len; i += 60)
      fprintf(stderr,"%.60s\n",seqbuf+i);
  }
#endif

  { int gaps, plen;
    int i, j, p;
  
    gaps = 0;
    plen = 0;
    for (i = 0; seqbuf[i] != '\0'; i++)
      if (seqbuf[i] == 'n')
        { for (j = i+1; seqbuf[j] != '\0'; j++)
            if (seqbuf[j] != 'n')
              break;
          gaps += 1;
          i = j-1;
        }
      else
        plen += 1;

    if (gaps != numgaps)
      { fprintf(stderr,"Header gap count (%d) and N-gap count (%d) differ\n",
                       numgaps,gaps);
        exit (1);
      }

    packseq = (char *) malloc(plen+1);

    gaps = 0;
    p = 0;
    ctglist[gaps].insert_pnt = 0;
    ctglist[gaps].lft_end    = 0;
    for (i = 0; seqbuf[i] != '\0'; i++)
      if (seqbuf[i] == 'n')
        { ctglist[gaps].length = p - ctglist[gaps].insert_pnt;
          plen = 1;
          for (j = i+1; seqbuf[j] != '\0'; j++)
            if (seqbuf[j] != 'n')
              break;
            else
              plen += 1;
          gaplist[gaps].gap_length = plen;
          gaps += 1;
          ctglist[gaps].insert_pnt = p;
          ctglist[gaps].lft_end    = ctglist[gaps-1].lft_end
                                   + (ctglist[gaps-1].length + plen);
          i = j-1;
        }
      else
        packseq[p++] = seqbuf[i];
    packseq[p] = '\0';
    ctglist[gaps].length = p - ctglist[gaps].insert_pnt;
  }

  scaf = (Scaffold *) malloc(sizeof(Scaffold));
  scaf->num_gaps = numgaps;
  scaf->gaps = gaplist;
  scaf->ctgs = ctglist;
  scaf->length = ctglist[numgaps].lft_end + ctglist[numgaps].length;
  scaf->packed_seq = packseq;

  return (scaf);
}

Scaffold_List *Read_Multi_Fasta(char *fname)
{ FILE          *ifile;
  Scaffold_List *list = NULL, *fing;
  Scaffold      *S;
  int            first;

  if ((ifile = fopen(fname,"r")) == NULL)
    { fprintf(stderr,"Cannot open file '%s'\n",fname);
      exit (1);
    }

  fing = NULL;
  first = 1;
  while ((S = Read_Next_Scaffold(ifile,first)) != NULL)
    { if (fing == NULL)
        list = fing = (Scaffold_List *) malloc(sizeof(Scaffold_List));
      else
        fing = fing->next = (Scaffold_List *) malloc(sizeof(Scaffold_List));
      fing->scaffold = S;
      first = 0;
    }
  if (fing != NULL)
    fing->next = NULL;
  else
    list = NULL;

  fclose(ifile);

  return (list);
}

void Complement_Scaffold(Scaffold *S)
{ int i, j, plen, slen;
  Scaffold_Gap T;
  Scaffold_Tig U;

  plen = S->ctgs[S->num_gaps].insert_pnt + S->ctgs[S->num_gaps].length;
  slen = S->length;
  if(S->packed_seq!=NULL)  Complement_Seq(S->packed_seq);
  j = S->num_gaps-1;
  for (i = 0; i < j; i++, j--)
    { T = S->gaps[i];
      S->gaps[i] = S->gaps[j];
      S->gaps[j] = T;
    }
  j = S->num_gaps;
  for (i = 0; i < j; i++, j--)
    { U = S->ctgs[i];
      S->ctgs[i] = S->ctgs[j];
      S->ctgs[j] = U;
      S->ctgs[i].insert_pnt = plen - (S->ctgs[i].insert_pnt+S->ctgs[i].length);
      S->ctgs[j].insert_pnt = plen - (S->ctgs[j].insert_pnt+S->ctgs[j].length);
      S->ctgs[i].lft_end    = slen - (S->ctgs[i].lft_end+S->ctgs[i].length);
      S->ctgs[j].lft_end    = slen - (S->ctgs[j].lft_end+S->ctgs[j].length);
    }
  if (i == j)
    { S->ctgs[j].insert_pnt = plen - (S->ctgs[j].insert_pnt+S->ctgs[j].length);
      S->ctgs[j].lft_end    = slen - (S->ctgs[j].lft_end+S->ctgs[j].length);
    }
}

void Free_Scaffold(Scaffold *S)
{ free(S->gaps);
  free(S->ctgs);
  free(S->packed_seq);
  free(S);
}

static AggregateString *Build_Indexable_String(Scaffold_List *S)
{ int  *idx, *sfn;
  char *str;
  int   len, ctg;
  AggregateString *agg;

  { Scaffold_List *a;
    Scaffold      *s;
    int i;

    len = 0;
    ctg = 0;
    for (a = S; a != NULL; a = a->next)
      { s = a->scaffold;
        for (i = 0; i <= s->num_gaps; i++)
          { len += s->ctgs[i].length + 6;
            ctg += 1;
          }
      }
  }

  { int as;

    as = sizeof(int);
    as = ( (sizeof(AggregateString) + as - 1) / as ) * as;
    agg = (AggregateString *)
             malloc(as + 2*ctg*sizeof(int) + len + 1);
    if (agg == NULL)
      { fprintf(stderr,"Out of memory (aggregate string)\n");
        exit (1);
      }
    idx = (int *) (((char *) agg) + as);
    sfn = idx + ctg;
    str = (char *) (sfn + ctg);
  }

  { Scaffold_List *a;
    Scaffold      *s;
    int   i, c, n;
    char *p;

    c = n = 0;
    p = str;
    for (a = S; a != NULL; a = a->next)
      { s = a->scaffold;
        for (i = 0; i <= s->num_gaps; i++)
          { strcpy(p,"nnnnnn");
            p += 6;
            strncpy(p,s->packed_seq + s->ctgs[i].insert_pnt,
                      s->ctgs[i].length);
            p += s->ctgs[i].length;
            sfn[c] = n;
            idx[c] = p-str;
            c += 1;
          }
        n += 1;
      }
  }

  agg->bigstr = str;
  agg->parts  = idx;
  agg->scafs  = sfn;
  agg->slen   = len;
  agg->plen   = ctg;
  return (agg);
}

static int Locate_Local(AggregateString *cmp, int pos)
{ int l, r, m;

  l = 0;
  r = cmp->plen;
  while (l < r)
    { m = (l+r)/2;
      if (pos <= cmp->parts[m])
        r = m;
      else
        l = m+1;
    } 
  return (r);
}

Local_Address *MSORT_addr;
#define USE_ONLY_FORWARD_SEGS
#ifdef USE_ONLY_FORWARD_SEGS
Local_Segment *MSORT_segs;
#endif


int MSORT(const void *l, const void *r)
{ Local_Address *x, *y;

#ifdef USE_ONLY_FORWARD_SEGS
  Local_Segment *X, *Y; 
 int xrev,yrev;
#endif

  x = MSORT_addr + *((int *) l);
  y = MSORT_addr + *((int *) r);
  if (x->ascaf != y->ascaf)
    return (x->ascaf - y->ascaf);
  else if (x->bscaf != y->bscaf)
    return (x->bscaf - y->bscaf);
  else if (x->acntg != y->acntg)
    return (x->acntg - y->acntg);
  else if (x->bcntg != y->bcntg)
      return (x->bcntg - y->bcntg);

#ifdef USE_ONLY_FORWARD_SEGS
  X = MSORT_segs + *((int *) l);
  Y = MSORT_segs + *((int *) r);
  assert(X->abpos<X->aepos);
  assert(Y->abpos<Y->aepos);
  xrev = ( X->bbpos<X->bepos) ? 0 : 1;
  yrev = ( Y->bbpos<Y->bepos) ? 0 : 1;
  return (xrev - yrev);
#endif

  return(0);

}

static Local_Pool *Find_All_Locals(Scaffold_List *AS, Scaffold_List *BS)
{ AggregateString *acmp, *bcmp;
  int             *perm;
  int              nseg;
  Local_Segment   *segs;
  Local_Address   *addr;
  static Local_Pool pool;
 
  acmp = Build_Indexable_String(AS);
  bcmp = Build_Indexable_String(BS);
  segs = Find_Local_Segments(acmp->bigstr,acmp->slen,
                             bcmp->bigstr,bcmp->slen,LOCAL_BOTH,
                             MIN_SEG_LEN,SEG_ERATE,&nseg);

  { int i;

    addr = (Local_Address *) malloc(nseg*sizeof(Local_Address));
    perm = (int *) malloc(nseg*sizeof(int));
    for (i = 0; i < nseg; i++)
      { int ac, bc;

        perm[i] = i;
        addr[i].acntg = ac = Locate_Local(acmp,segs[i].abpos);
        addr[i].bcntg = bc = Locate_Local(bcmp,segs[i].bbpos);
        addr[i].ascaf = acmp->scafs[ac];
        addr[i].bscaf = bcmp->scafs[bc];

        if (ac == 0)
         ac = 6;
        else
          ac = acmp->parts[ac-1] + 6;
        if (bc == 0)
          bc = 6;
        else
          bc = bcmp->parts[bc-1] + 6;
        segs[i].abpos -= ac;
        segs[i].aepos -= ac;
        segs[i].bbpos -= bc;
        segs[i].bepos -= bc;
        segs[i].ldiag -= (ac-bc);
        segs[i].hdiag -= (ac-bc);
      }
  }

  MSORT_addr = addr;
#ifdef USE_ONLY_FORWARD_SEGS
  MSORT_segs = segs;
#endif
  qsort(perm,nseg,sizeof(int),MSORT);
  
#ifdef DEBUG_LOCAL
  { int i;

    fprintf(stderr,"\n  A index\n");
    for (i = 0; i < acmp->plen; i++)
      fprintf(stderr,"   %6d %3d\n",acmp->parts[i],acmp->scafs[i]);
    fprintf(stderr,"\n  B index\n");
    for (i = 0; i < bcmp->plen; i++)
      fprintf(stderr,"   %6d %3d\n",bcmp->parts[i],bcmp->scafs[i]);
    fprintf(stderr,"\nHits\n");
    for (i = 0; i < nseg; i++)
      fprintf(stderr,"%3d: %3d [%6d,%6d] vs [%6d,%6d] to (%d,%d) [%d,%d]\n",
             i,perm[i],segs[i].abpos,segs[i].aepos,
             segs[i].bbpos,segs[i].bepos,addr[i].acntg,addr[i].bcntg,
             addr[i].ascaf,addr[i].bscaf);
  }
#endif

  free(acmp);
  free(bcmp);

  { int a, w;
    Local_Segment t;
    Local_Address u;

    for (a = 0; a < nseg; a++)
      { t = segs[a];
        u = addr[a];
        while ((w = perm[a]) >= 0)
          { if (perm[w] < 0)
              { segs[a] = t;
                addr[a] = u;
              }
            else
              { segs[a] = segs[w];
                addr[a] = addr[w];
              }
            perm[a] = -1;
            a = w;
          }
      }
  }

  free(perm);

#ifdef DEBUG_LOCAL
  { int i;

    fprintf(stderr,"\nSorted Hits\n");
    for (i = 0; i < nseg; i++)
      fprintf(stderr,"%3d: %3d [%6d,%6d] vs [%6d,%6d] to (%d,%d) [%d,%d]\n",
             i,perm[i],segs[i].abpos,segs[i].aepos,
             segs[i].bbpos,segs[i].bepos,addr[i].acntg,addr[i].bcntg,
             addr[i].ascaf,addr[i].bscaf);
  }
#endif

  pool.num_locals = nseg;
  pool.locals     = segs;
  pool.address    = addr;
  return (&pool);
}

Segment *Find_All_Overlaps(Scaffold *AF, Scaffold *BF, Local_Pool *pool,
                                  int as, int bs, int ac, int bc,
                                  int *fing, int *segs, int comp)
{ int      i, alen;
  int      bf;
  Local_Address *addr;
  Local_Segment *locl;
  Segment *list;

  list  = NULL;
  *segs = 0;

  bf   = *fing;
  addr = pool->address;
  locl = pool->locals;
  for (i = 0; i <= AF->num_gaps; i++)
    { int j, blen;
      Segment *last;

      alen = AF->ctgs[i].length;
      last = list;

      for (j = 0; j <= BF->num_gaps; j++)
        { Local_Overlap *ovl;
          int n;

	  // The following looks odd, but remember, if comp==1,
	  // then the scaffold was reversed, but the hits weren't
	  //     --- ALH, 11/7/01
          if (comp)
            blen = BF->ctgs[BF->num_gaps-j].length;
          else
            blen = BF->ctgs[j].length;

#ifdef USE_ONLY_FORWARD_SEGS
	  //if reversed comparison, skip over forward-oriented matches
	  if(comp){
	    while (bf < pool->num_locals &&
		   addr[bf].ascaf == as && addr[bf].bscaf == bs &&
		   addr[bf].acntg == ac+i && addr[bf].bcntg == bc+j&&
		   locl[bf].bbpos < locl[bf].bepos)
	      bf++;
	  }
#endif

          n = bf;
          while (n < pool->num_locals &&
                 addr[n].ascaf == as && addr[n].bscaf == bs &&
                 addr[n].acntg == ac+i && addr[n].bcntg == bc+j){
#ifdef USE_ONLY_FORWARD_SEGS
	    //if forward comparison, stop on reversed match
	    if(!comp&&locl[n].bbpos > locl[bf].bepos)break;
#endif
            n += 1;
	  }

#ifdef DEBUG_LOCAL
          { int x;
            fprintf(stderr,"\nInput to Find_Local_Overlap (%d,%d) [%d,%d] fing = %d\n",
                   ac+i,bc+j,as,bs,bf);
            for (x = bf; x < n; x++)
              fprintf(stderr,"%3d: [%6d,%6d] vs [%6d,%6d] to (%d,%d) [%d,%d]\n",
                     x,locl[x].abpos,locl[x].aepos,
                     locl[x].bbpos,locl[x].bepos,addr[x].acntg,addr[x].bcntg,
                     addr[x].ascaf,addr[x].bscaf);
          }
#endif

          ovl = Find_Local_Overlap(alen,blen,comp,0,locl+bf,n-bf,20,OVL_ERATE);

#ifdef USE_ONLY_FORWARD_SEGS
	  // if forward comparison, skip over reversed matches
	  if(!comp){
	    bf=n;
	    while (bf < pool->num_locals &&
		   addr[bf].ascaf == as && addr[bf].bscaf == bs &&
		   addr[bf].acntg == ac+i && addr[bf].bcntg == bc+j)
	      bf++;
	    n=bf;
	  }
#endif


          while (ovl != NULL)
            { Segment *newseg;

              newseg = (Segment *) malloc(sizeof(Segment));
              newseg->next     = list;
              newseg->a_contig = i;
              newseg->b_contig = j;

              if (ovl->begpos < 0)
                newseg->alow = 0;
              else
                newseg->alow = ovl->begpos;
              if (ovl->begpos > 0)
                newseg->blow = 0;
              else
                newseg->blow = - ovl->begpos;

              if (ovl->endpos > 0)
                newseg->ahgh = alen;
              else
                newseg->ahgh = alen + ovl->endpos;
              if (ovl->endpos < 0)
                newseg->bhgh = blen;
              else
                newseg->bhgh = blen - ovl->endpos;

              newseg->overlap  = ovl;
              list = newseg;
              *segs += 1;

#ifdef DEBUG_LOCAL
              Print_Local_Overlap(stdout,ovl,4);
#endif

              ovl = Find_Local_Overlap(alen,blen,comp,1,
                                       locl+bf,n-bf,20,OVL_ERATE);
            }

          bf  = n;
        }

      // The following looks odd, but remember, if comp==1,
      // then the scaffold was reversed, but the hits weren't
      //     --- ALH, 11/7/01
      if (comp)
        { Segment *r, *s;
          for (r = last; list != last; list = s)
            { s = list->next;
              list->b_contig = BF->num_gaps - list->b_contig;
              list->next = r;
              r = list;
            }
          list = r;
        }
    }

  *fing = bf;
  return (list);
}

void Free_Segments_ScafComp(Segment *seglist)
{ Segment *s, *t;

  for (s = seglist; s != NULL; s = t)
    { t = s->next;
      free(s->overlap);
      free(s);
    }
}

int     MaxAlign  = -1;
int     MaxBucket = -1;
COvlps  *CtgOvls  = NULL;
COvlps **ABuckets = NULL;
COvlps **BBuckets = NULL;
#ifdef XFIG
static double scale_factor;
static FILE *figfile;

#define MAP(x)  ((int) ((x)*scale_factor + 600))
#define UNMAP(y)  (int)(((y)-600.)/scale_factor)

#endif

#ifdef XFIG
static void Draw_Matrix(Segment *seglist, int numsegs,
                 Scaffold *AF, Scaffold *BF,
                 int varwin)
{ Segment *f;
  int i, j;

  if (AF->length > BF->length)
    scale_factor = 6000. / AF->length;
  else
    scale_factor = 6000. / BF->length;
  fprintf(figfile,"#FIG 3.2\n");
  fprintf(figfile,"Landscape\n");
  fprintf(figfile,"Center\n");
  fprintf(figfile,"Inches\n");
  fprintf(figfile,"Letter\n");
  fprintf(figfile,"100.\n");
  fprintf(figfile,"Single\n");
  fprintf(figfile,"-2\n");
  fprintf(figfile,"1200 2\n");
  fprintf(figfile,"2 2 0 3 0 -1 100 0 -1 4.00 0 0 0 0 0 5\n\t");
  fprintf(figfile," %d %d",MAP(0),MAP(0));
  fprintf(figfile," %d %d",MAP(0),MAP(BF->length));
  fprintf(figfile," %d %d",MAP(AF->length),MAP(BF->length));
  fprintf(figfile," %d %d",MAP(AF->length),MAP(0));
  fprintf(figfile," %d %d\n",MAP(0),MAP(0));

  // draw boxes for contigs
  for (i = 0; i <= AF->num_gaps; i++)
    for (j = 0; j <= BF->num_gaps; j++)
      { int xl, yl, xh, yh;
        xl = MAP(AF->ctgs[i].lft_end);
        yl = MAP(BF->ctgs[j].lft_end);
        xh = MAP(AF->ctgs[i].lft_end + AF->ctgs[i].length);
        yh = MAP(BF->ctgs[j].lft_end + BF->ctgs[j].length);
        fprintf(figfile,"2 2 0 2 1 -1 101 0 -1 4.00 0 0 0 0 0 5\n\t");
        fprintf(figfile," %d %d",xl,yl);
        fprintf(figfile," %d %d",xl,yh);
        fprintf(figfile," %d %d",xh,yh);
        fprintf(figfile," %d %d",xh,yl);
        fprintf(figfile," %d %d\n",xl,yl);
      }

  // draw gap symbols?
  for (i = 0; i < AF->num_gaps; i++)
    { int xl, xh, yl;
      xl = MAP(AF->ctgs[i+1].lft_end - varwin * AF->gaps[i].gap_var);
      xh = MAP(AF->ctgs[i+1].lft_end + varwin * AF->gaps[i].gap_var);
      yl = 500;
      fprintf(figfile,"2 1 0 1 1 -1 102 0 -1 4.00 0 0 0 1 1 2\n\t");
      fprintf(figfile,"0 0 1.00 60.00 75.00\n\t");
      fprintf(figfile,"0 0 1.00 60.00 75.00\n\t");
      fprintf(figfile," %d %d",xl,yl);
      fprintf(figfile," %d %d\n",xh,yl);
    }
  // draw gap symbols?
  for (i = 0; i < BF->num_gaps; i++)
    { int yl, yh, xl;
      yl = MAP(BF->ctgs[i+1].lft_end - varwin * BF->gaps[i].gap_var);
      yh = MAP(BF->ctgs[i+1].lft_end + varwin * BF->gaps[i].gap_var);
      xl = 500;
      fprintf(figfile,"2 1 0 1 1 -1 102 0 -1 4.00 0 0 0 1 1 2\n\t");
      fprintf(figfile,"0 0 1.00 60.00 75.00\n\t");
      fprintf(figfile,"0 0 1.00 60.00 75.00\n\t");
      fprintf(figfile," %d %d",xl,yl);
      fprintf(figfile," %d %d\n",xl,yh);
    }

  // draw contig overlaps
  for (f = seglist; f != NULL; f = f->next)
    { int xl, yl, xh, yh;
      xl = MAP(AF->ctgs[f->a_contig].lft_end + f->alow);
      xh = MAP(AF->ctgs[f->a_contig].lft_end + f->ahgh);
      yl = MAP(BF->ctgs[f->b_contig].lft_end + f->blow);
      yh = MAP(BF->ctgs[f->b_contig].lft_end + f->bhgh);
      fprintf(figfile,"2 1 0 2 2 -1 102 0 -1 4.00 0 0 0 0 0 2\n\t");
      fprintf(figfile," %d %d",xl,yl);
      fprintf(figfile," %d %d\n",xh,yh);
    }
}
#endif

#ifdef STATS
static int gapcount, blockcount, linkcount, edgecount;
static int vertcount, boundcount;

void Start_Stats()
{ gapcount = blockcount = linkcount = edgecount = 0;
  vertcount = boundcount = 0;
}

void Print_Stats()
{ fprintf(stderr,"  Block Boundaries Reached  = %3d\n",blockcount);
  fprintf(stderr,"  Gap Portals Reached       = %3d\n",gapcount);
  fprintf(stderr,"  Traverse Calls Made (Bnd) = %3d(%d)\n",linkcount,boundcount);
  fprintf(stderr,"  Edges in Graphs          <= %3d\n",edgecount);
  fprintf(stderr,"  Vertices in Graph         = %3d\n",vertcount);
}
#endif

int Link_Horizontal(Scaffold *A, Scaffold *B, int varwin,
                           int i, int j, int low, int hgh, COvlps *source)
{ int k, var;
  int terminal;
  COvlps *afing;
#ifdef XFIG
  int xl, zl, zh, gapped;

  if (low != hgh)
    { xl = A->ctgs[i].lft_end + A->ctgs[i].length + B->ctgs[j].lft_end;
      zl = MAP(xl - hgh > A->ctgs[i].lft_end + A->ctgs[i].length ? 
	       xl - hgh : A->ctgs[i].lft_end + A->ctgs[i].length);
      zh = MAP(xl - low > A->ctgs[i].lft_end + A->ctgs[i].length ? 
	       xl - low : A->ctgs[i].lft_end + A->ctgs[i].length);
      xl = MAP(B->ctgs[j].lft_end);
      gapped = 1;
    }
  else
    { xl  = MAP(A->ctgs[i].lft_end + A->ctgs[i].length);
      zl = MAP(low);
      zh = MAP(hgh);
      gapped = 0;
    }
#endif

#ifdef STATS
  linkcount += 1;
#endif
  if (i == A->num_gaps)
    {
#ifdef DEBUG_ALIGN
      fprintf(stderr,"    Ray reaches A-boundary\n"); 
#endif
      return (1);
    }

  terminal = 0;
  low += A->gaps[i].gap_length;
  hgh += A->gaps[i].gap_length;
  var  = (int)(A->gaps[i].gap_var * varwin);
  afing = ABuckets[i+1];

  for (k = j; k <= B->num_gaps; k++)
    { int x, beg, end;

      x = B->ctgs[k].lft_end;
      if (low-var > x)
        beg = low-var;
      else
        beg = x; 
      x += B->ctgs[k].length;
      if (hgh+var < x)
        end = hgh+var;
      else
        end = x;
      if (beg <= end)
        { COvlps *af;
          int     pnt, score;

#ifdef STATS
          blockcount += 1;
#endif
#ifdef XFIG
          { int xh, yl, yh;

            yl = MAP(beg);
            yh = MAP(end);
            xh  = MAP(A->ctgs[i+1].lft_end);
            fprintf(figfile,"2 3 0 1 4 4  98 0 35 4.01 1 0 0 0 0 5\n\t");
            if (gapped)
              { fprintf(figfile," %d %d",zl,xl);
                fprintf(figfile," %d %d",zh,xl);
                fprintf(figfile," %d %d",xh,yl);
                fprintf(figfile," %d %d",xh,yh);
                fprintf(figfile," %d %d\n",zl,xl);
              }
            else
              { fprintf(figfile," %d %d",xl,zl);
                fprintf(figfile," %d %d",xl,zh);
                fprintf(figfile," %d %d",xh,yh);
                fprintf(figfile," %d %d",xh,yl);
                fprintf(figfile," %d %d\n",xl,zl);
              }
          }
#endif
#ifdef DEBUG_RAYSHOOT
          fprintf(stderr,"    Ray overlaps B%d[%d,%d]\n",k,beg,end);
#endif
          while (afing != NULL && afing->seg->b_contig < k)
            afing = afing->Alink; 
          for (af = afing; af != NULL && af->seg->b_contig == k; af = af->Alink)
            if (af->seg->overlap->begpos <= 0)
              { pnt = B->ctgs[k].lft_end - af->seg->overlap->begpos;
                if (beg <= pnt && pnt <= end)
                  { score = af->seg->overlap->length;
                    if (source != NULL)
                      score += source->best;
                    if (score > af->best)
                      { af->best  = score;
                        af->trace = source;
                      }
#ifdef STATS
                    edgecount += 1;
#endif
#ifdef DEBUG_ALIGN
                    fprintf(stderr,"    Finds (%d,%d) hangs(%d,%d)= %d\n",
                           af->seg->a_contig,af->seg->b_contig,af->seg->overlap->begpos,af->seg->overlap->endpos,af->best);
#endif

// If you are reading this, it is possible you are confused as to whether we should 
// recurse on segments reached from Link_Horizontal or Link_Vertical; it appears that
// we should NOT.  Instead, Align_Scaffold proceeds as follows, in three stages:
// (a) determine all the segments that can be reached directly from starting 
// in a gap on the entry border; each such segment has its "best" set to its length
// (b) determine all segments that can be reached by starting from a segment that
// begins on the entry border; again, reachable segments have "best" set to length
// (c) consider all segments with best >= 0 (i.e. those that can be reached from
// the entry border via gaps or segments previously visited); determine where these
// segments can reach--any segment they can reach itself becomes reachable.
// To avoid redundant computation for segments in path tails, once we reach
// a segment, we do not recursively follow it but instead mark it as accessible (i.e.
// set its best > 0).  For this to work, we must be sure to evaluate segments in the
// right order -- it would work to go in increasing order of A contigs, with subsorting
// on increasing B contigs, or vice versa ... but critically we can't go in reverse order
// on either.
//
// As an aside, it appears that repeated evaluation of the same gap intervals is not
// prevented by this scheme ....  Better might be to proceed one A contig or gap at a time.
//
                  }
              }
        }

      if (k < B->num_gaps)
        { int y, l, del;

          if (low-var > x)
            beg = low-var;
          else
            beg = x; 
          y = B->ctgs[k+1].lft_end;
          if (hgh+var < y)
	    if(hgh+var<x){
	      continue;
	    } else {
	      end = hgh+var;
	    }
	  else
            end = y;
          del = (int)(B->gaps[k].gap_var * varwin);
          if (beg <= end)
            { l = A->ctgs[i+1].lft_end;
#ifdef STATS
              gapcount += 1;
#endif
#ifdef XFIG
              { int xh, yl, yh;

                yl = MAP(beg);
                yh = MAP(end);
                xh  = MAP(A->ctgs[i+1].lft_end);
                fprintf(figfile,"2 3 0 1 4 4  98 0 35 4.02 1 0 0 0 0 5\n\t");
                if (gapped)
                  { fprintf(figfile," %d %d",zl,xl);
                    fprintf(figfile," %d %d",zh,xl);
                    fprintf(figfile," %d %d",xh,yl);
                    fprintf(figfile," %d %d",xh,yh);
                    fprintf(figfile," %d %d\n",zl,xl);
                  }
                else
                  { fprintf(figfile," %d %d",xl,zl);
                    fprintf(figfile," %d %d",xl,zh);
                    fprintf(figfile," %d %d",xh,yh);
                    fprintf(figfile," %d %d",xh,yl);
                    fprintf(figfile," %d %d\n",xl,zl);
                  }
              }
#endif
#ifdef DEBUG_RAYSHOOT
              fprintf(stderr,"    Ray overlaps gap B%d:%d to B%d:%d\n",k,beg,k+1,end);
#endif

              if (Link_Vertical(A,B,varwin,
                                i+1,k,l - (end-x),l - (beg-x),source)) 
                terminal = 1;
            }
          else if (beg > y && beg-del < y)
            { l = A->ctgs[i+1].lft_end;
#ifdef STATS
              gapcount += 1;
#endif
#ifdef DEBUG_RAYSHOOT
              fprintf(stderr,"    Ray indirects gap B%d to B%d @ %d\n",k,k+1,beg-del);
#endif
              if (Link_Vertical(A,B,varwin,
                  i+1,k,l - (beg-x),l - (beg-x),source)) 
                terminal = 1;
            }
          var += del;
        }

      else if (end == x)
        { terminal = 1;
#ifdef XFIG
          { int tl, th, yh;

            tl = A->ctgs[i+1].lft_end - ((hgh+var) - end);
            if (tl < A->ctgs[i].lft_end + A->ctgs[i].length)
              tl = A->ctgs[i].lft_end + A->ctgs[i].length;
            yh = MAP(end);
            tl = MAP(tl);
            // th = MAP(A->ctgs[i+1].lft_end);
	    th = A->ctgs[i+1].lft_end - ((low-var) -end);
            if (th > A->ctgs[i+1].lft_end)
              th = A->ctgs[i+1].lft_end ;  /* tl -> th; did this fix anything? */
    	    th = MAP(th);
            fprintf(figfile,"2 3 0 1 4 4  98 0 35 4.03 1 0 0 0 0 5\n\t");
            if (gapped)
              { fprintf(figfile," %d %d",zl,xl);
                fprintf(figfile," %d %d",zh,xl);
                fprintf(figfile," %d %d",th,yh);
                fprintf(figfile," %d %d",tl,yh);
                fprintf(figfile," %d %d\n",zl,xl);
              }
            else
              { fprintf(figfile," %d %d",xl,zl);
                fprintf(figfile," %d %d",xl,zh);
                fprintf(figfile," %d %d",tl,yh);
                fprintf(figfile," %d %d",th,yh);
                fprintf(figfile," %d %d\n",xl,zl);
              }
          }
#endif
#ifdef DEBUG_ALIGN
          fprintf(stderr,"    Ray reaches B-boundary\n"); 
#endif
        }
    }

  return (terminal);
}

int Link_Vertical(Scaffold *A, Scaffold *B, int varwin,
                         int i, int j, int low, int hgh, COvlps *source)
{ int k, var;
  int terminal;
  COvlps *bfing;
#ifdef XFIG
  int yl, zl, zh, gapped;

  if (low != hgh)
    { yl = B->ctgs[j].lft_end + B->ctgs[j].length + A->ctgs[i].lft_end;
      zl = MAP(yl - hgh > B->ctgs[j].lft_end + B->ctgs[j].length ? 
	       yl - hgh : B->ctgs[j].lft_end + B->ctgs[j].length);
      zh = MAP(yl - low > B->ctgs[j].lft_end + B->ctgs[j].length ? 
	       yl - low : B->ctgs[j].lft_end + B->ctgs[j].length);
      yl = MAP(A->ctgs[i].lft_end);
      gapped = 1;
    }
  else
    { yl  = MAP(B->ctgs[j].lft_end + B->ctgs[j].length);
      zl = MAP(low);
      zh = MAP(hgh);
      gapped = 0;
    }
#endif

#ifdef STATS
  linkcount += 1;
#endif
  if (j == B->num_gaps)
    {
#ifdef DEBUG_ALIGN
      fprintf(stderr,"    Ray reaches B-boundary\n"); 
#endif
      return (1);
    }

  terminal = 0;
  low += B->gaps[j].gap_length;
  hgh += B->gaps[j].gap_length;
  var  = (int)(B->gaps[j].gap_var * varwin);
  bfing = BBuckets[j+1];

  for (k = i; k <= A->num_gaps; k++)
    { int x, beg, end;

      x = A->ctgs[k].lft_end;
      if (low-var > x)
        beg = low-var;
      else
        beg = x; 
      x += A->ctgs[k].length;
      if (hgh+var < x)
        end = hgh+var;
      else
        end = x;
      if (beg <= end)
        { COvlps *bf;
          int     pnt, score;

#ifdef STATS
          blockcount += 1;
#endif
#ifdef XFIG
          { int yh, xl, xh;

            xl = MAP(beg);
            xh = MAP(end);
            yh  = MAP(B->ctgs[j+1].lft_end);
            fprintf(figfile,"2 3 0 1 5 5  98 0 35 4.04 1 0 0 0 0 5\n\t");
            if (gapped)
              { fprintf(figfile," %d %d",yl,zl);
                fprintf(figfile," %d %d",yl,zh);
                fprintf(figfile," %d %d",xl,yh);
                fprintf(figfile," %d %d",xh,yh);
                fprintf(figfile," %d %d\n",yl,zl);
              }
            else
              { fprintf(figfile," %d %d",zl,yl);
                fprintf(figfile," %d %d",zh,yl);
                fprintf(figfile," %d %d",xh,yh);
                fprintf(figfile," %d %d",xl,yh);
                fprintf(figfile," %d %d\n",zl,yl);
              }
          }
#endif
#ifdef DEBUG_RAYSHOOT
          fprintf(stderr,"    Ray overlaps A%d[%d,%d]\n",k,beg,end);
#endif
          while (bfing != NULL && bfing->seg->a_contig < k){
#ifdef DEBUG_RAYSHOOT
	    fprintf(stderr,"Skipping over segment involving A%d, ahg:%d\n",
		    bfing->seg->a_contig,
		    bfing->seg->overlap->begpos);
#endif
            bfing = bfing->Blink; 
	  }
          for (bf = bfing; bf != NULL && bf->seg->a_contig == k; bf = bf->Blink){
#ifdef DEBUG_RAYSHOOT
	    fprintf(stderr,"Examining segment involving A%d ahg:%d\n",
		    bf->seg->a_contig,
		    bf->seg->overlap->begpos);
#endif
            if (bf->seg->overlap->begpos >= 0)
              { pnt = A->ctgs[k].lft_end + bf->seg->overlap->begpos;
                if (beg <= pnt && pnt <= end)
                  { score = bf->seg->overlap->length;
                    if (source != NULL)
                      score += source->best;
                    if (score > bf->best)
                      { bf->best  = score;
                        bf->trace = source;
                      }
#ifdef STATS
                    edgecount += 1;
#endif
#ifdef DEBUG_ALIGN
                    fprintf(stderr,"    Finds (%d,%d) hangs(%d,%d)= %d\n",
                           bf->seg->a_contig,bf->seg->b_contig,bf->seg->overlap->begpos,bf->seg->overlap->endpos,bf->best);
#endif
                  }
              }
	  }
        }

      if (k < A->num_gaps)
        { int y, l, del;

          if (low-var > x)
            beg = low-var;
          else
            beg = x; 
          y = A->ctgs[k+1].lft_end;
          if (hgh+var < y)
	    if(hgh+var < x){
	      continue;
	    } else {
	      end = hgh+var;
	    }
          else
            end = y;
          del = (int)(A->gaps[k].gap_var * varwin);
          if (beg <= end)
            { l = B->ctgs[j+1].lft_end;
#ifdef STATS
              gapcount += 1;
#endif
#ifdef XFIG
              { int yh, xl, xh;

                xl = MAP(beg);
                xh = MAP(end);
                yh  = MAP(B->ctgs[j+1].lft_end);
                fprintf(figfile,"2 3 0 1 5 5  98 0 35 4.05 1 0 0 0 0 5\n\t");
                if (gapped)
                  { fprintf(figfile," %d %d",yl,zl);
                    fprintf(figfile," %d %d",yl,zh);
                    fprintf(figfile," %d %d",xl,yh);
                    fprintf(figfile," %d %d",xh,yh);
                    fprintf(figfile," %d %d\n",yl,zl);
                  }
                else
                  { fprintf(figfile," %d %d",zl,yl);
                    fprintf(figfile," %d %d",zh,yl);
                    fprintf(figfile," %d %d",xh,yh);
                    fprintf(figfile," %d %d",xl,yh);
                    fprintf(figfile," %d %d\n",zl,yl);
                  }
              }
#endif
#ifdef DEBUG_RAYSHOOT
              fprintf(stderr,"    Ray overlaps gap A%d:%d to A%d:%d\n",k,beg,k+1,end);
#endif
              if (Link_Horizontal(A,B,varwin,
                                  k,j+1,l - (end-x),l - (beg-x),source))
                terminal = 1; 
            }
          else if (beg > y && beg-del < y)
            { l = B->ctgs[j+1].lft_end;
#ifdef STATS
              gapcount += 1;
#endif
#ifdef DEBUG_RAYSHOOT
              fprintf(stderr,"    Ray indirects gap A%d to A%d @ %d\n",k,k+1,beg-del);
#endif
              if (Link_Horizontal(A,B,varwin,
                                  k,j+1,l - (beg-x),l - (beg-x),source))
                terminal = 1;
            }
          var += del;
        }

      else if (end == x)
        { terminal = 1;
#ifdef XFIG
          { int tl, th, xh;

            tl = B->ctgs[j+1].lft_end - ((hgh+var) - end);
            if (tl < B->ctgs[j].lft_end + B->ctgs[j].length)
              tl = B->ctgs[j].lft_end + B->ctgs[j].length;
            xh = MAP(end);
            tl = MAP(tl);
            // th = MAP(B->ctgs[j+1].lft_end);
	    th = B->ctgs[j+1].lft_end - ((low-var) -end);
            if (th > B->ctgs[j+1].lft_end)
              th = B->ctgs[j+1].lft_end ; /* tl -> th : did this fix anything? */
    	    th = MAP(th);
            fprintf(figfile,"2 3 0 1 5 5  98 0 35 4.06 1 0 0 0 0 5\n\t");
            if (gapped)
              { fprintf(figfile," %d %d",yl,zl);
                fprintf(figfile," %d %d",yl,zh);
                fprintf(figfile," %d %d",xh,th);
                fprintf(figfile," %d %d",xh,tl);
                fprintf(figfile," %d %d\n",yl,zl);
              }
            else
              { fprintf(figfile," %d %d",zl,yl);
                fprintf(figfile," %d %d",zh,yl);
                fprintf(figfile," %d %d",xh,tl);
                fprintf(figfile," %d %d",xh,th);
                fprintf(figfile," %d %d\n",zl,yl);
              }
          }
#endif
#ifdef DEBUG_ALIGN
          fprintf(stderr,"    Ray reaches A boundary\n");
#endif
        }
    }

  return (terminal);
}

#ifdef AHANG_BAND_TEST
Segment *Align_Scaffold(Segment *seglist, int numsegs, int varwin,
			Scaffold *AF, Scaffold *BF, int *best,
			int bandbeg, int bandend)
#else
Segment *Align_Scaffold(Segment *seglist, int numsegs, int varwin,
                               Scaffold *AF, Scaffold *BF, int *best)
#endif
{ int optc = 0, mval=-1;
 COvlps *optco = NULL; 
 int     score=0, term;
 
#ifdef XFIG
 static int Case=0; 
 char fname[50]; 
  
 sprintf(fname,"Graphic_case_%d.fig",Case++); 
 fprintf(stderr,"Align_Scaffold writing xfig to %s\n", fname); 
 figfile = fopen(fname,"w"); 
 Draw_Matrix(seglist,numsegs,AF,BF,varwin); 
#endif


#ifdef AHANG_BAND_TEST
  assert(bandbeg<=bandend);
  assert(bandend<=AF->length);
  assert(bandbeg>=-(BF->length));
#endif

  if (numsegs > MaxAlign)
    { MaxAlign = (int)(1.3*numsegs + 100);
      CtgOvls  = (COvlps *) realloc(CtgOvls,sizeof(COvlps)*MaxAlign);
      if (CtgOvls == NULL)
        { fprintf(stderr,"Out of memory allocating DP array\n");
          exit (1);
        }
    }

  if (AF->num_gaps + BF->num_gaps + 2 > MaxBucket)
    { MaxBucket = (int)(1.3*(AF->num_gaps + BF->num_gaps + 2) + 100);
      ABuckets  = (COvlps **) realloc(ABuckets,sizeof(COvlps *)*MaxBucket);
      if (ABuckets == NULL)
        { fprintf(stderr,"Out of memory allocating segment sort arrays\n");
          exit (1);
        }
    }
  BBuckets  = ABuckets + (AF->num_gaps+1);

  { int i,c;
    Segment *s;


    for (i = 0; i <= AF->num_gaps; i++)
      ABuckets[i] = NULL;
    for (i = 0; i <= BF->num_gaps; i++)
      BBuckets[i] = NULL;

    c = numsegs;
    for (s = seglist; s != NULL; s = s->next)
      { c -= 1;
        CtgOvls[c].seg = s;
        CtgOvls[c].best = -1;
        CtgOvls[c].trace = NULL;

#ifdef DEBUG_SEGORDER
	fprintf(stderr,"CtgOvls[%d] actg: %d bctg: %d\n",
		c,CtgOvls[c].seg->a_contig,
		CtgOvls[c].seg->b_contig);
#endif
	// push segment onto Alink list; this needs to result in all 
	// segments involving s->a_contig being linked together,
	// and the order of the elements should be such that
	// s->b_contig <= s->Alink->b_contig (if s->Alink != NULL)

        CtgOvls[c].Alink = ABuckets[s->a_contig];
        ABuckets[s->a_contig] = CtgOvls+c;
	if(ABuckets[s->a_contig]->Alink!=NULL)
	  assert(ABuckets[s->a_contig]->seg->b_contig <= ABuckets[s->a_contig]->Alink->seg->b_contig);

	// original code did something similar for BBuckets and Blink,
      }


    // push segment onto Blink list; this needs to result in all 
    // segments involving s->b_contig being linked together,
    // and the order of the elements should be such that
    // s->a_contig <= s->Blink->a_contig (if s->Blink != NULL)
    
    for(i=AF->num_gaps;i>=0;i--){
      COvlps *co;
      co = ABuckets[i];
      while(co!=NULL){
	co->Blink = BBuckets[co->seg->b_contig];
        BBuckets[co->seg->b_contig] = co;
	if(co->Blink!=NULL)
	  assert(co->seg->a_contig <= co->Blink->seg->a_contig);
	co=co->Alink;
      }
    }

  }

#ifdef DEBUG_ALIGN
  { Segment *s;
    COvlps  *c;
    int      i;

    fprintf(stderr,"\nAlign Scaffolds\n\n  Seglist:\n");
    for (s = seglist; s != NULL; s = s->next)
      fprintf(stderr,"    (%d,%d)\n",s->a_contig,s->b_contig);
    fprintf(stderr,"\n  A-Buckets:\n");
    for (i = 0; i <= AF->num_gaps; i++)
      { fprintf(stderr,"    %2d:",i);
        for (c = ABuckets[i]; c != NULL; c = c->Alink)
          fprintf(stderr," %d",c->seg->b_contig);
        fprintf(stderr,"\n");
      }
    fprintf(stderr,"\n  B-Buckets:\n");
    for (i = 0; i <= BF->num_gaps; i++)
      { fprintf(stderr,"    %2d:",i);
        for (c = BBuckets[i]; c != NULL; c = c->Blink)
          fprintf(stderr," %d",c->seg->a_contig);
        fprintf(stderr,"\n");
      }
    fprintf(stderr,"\n");
  }
#endif

  { int i;

    for (i = 0; i < AF->num_gaps; i++)
      { int low, hgh;
  
        low = AF->ctgs[i].lft_end + AF->ctgs[i].length;
        hgh = AF->ctgs[i+1].lft_end;

#ifdef AHANG_BAND_TEST
	if(low>bandend||hgh<bandbeg)continue;
	if(low<bandbeg)low=bandbeg;
	if(hgh>bandend)hgh=bandend;
#endif


#ifdef DEBUG_ALIGN
        fprintf(stderr,"  Start ray A[%d,%d] at %d,%d\n",low,hgh,0,i);
#endif

	term = Link_Horizontal(AF, BF, varwin, i, 0, 
			     // begin point negative by hgh - gap begin
			     -( hgh - (AF->ctgs[i].lft_end+AF->ctgs[i].length)),
			     // end point negative by low - gap begin
			     -( low - (AF->ctgs[i].lft_end+AF->ctgs[i].length)),
			     NULL);

#ifdef ALLOW_NO_OVERLAP_INTERLEAVINGS
	// COMMENT ABC: if we allow the no-overlap interleavings, we need to 
	// set term and test whether we found a solution with
	// no overlaps

        if (term && score > mval)
          { mval = score;
  #ifdef STATS
            edgecount += 1;
  #endif
            optc = -1; // special value indicating no segments used
          }

#else
	// COMMENT CBA: with overlaps required, we don't need to check for
	// whether we reached a boundary; if we did *and* it involved
	// an overlap (that is, a successful overlapping), we will discover it
	// when we test the relevant segment(s) below
#endif



      }
    for (i = 0; i < BF->num_gaps; i++)
      { int low, hgh;
   
        low = BF->ctgs[i].lft_end + BF->ctgs[i].length;
        hgh = BF->ctgs[i+1].lft_end;

	// adjust low and hgh to conform to the banding 
	if(-hgh>bandend||-low<bandbeg)continue;
	if(-hgh<bandbeg)hgh=-bandbeg;
	if(-low>bandend)low=-bandend;

#ifdef DEBUG_ALIGN
        fprintf(stderr,"  Start ray B[%d,%d] at %d,%d\n",low,hgh,0,i);
#endif


	term = Link_Vertical(AF, BF, varwin, 0, i, 
			     // begin point negative by hgh - gap begin
			     -( hgh - (BF->ctgs[i].lft_end+BF->ctgs[i].length)),
			     // end point negative by low - gap begin
			     -( low - (BF->ctgs[i].lft_end+BF->ctgs[i].length)),
			     NULL);

#ifdef ALLOW_NO_OVERLAP_INTERLEAVINGS

	// see comment ABC above

        if (term && score > mval)
          { mval = score;
  #ifdef STATS
            edgecount += 1;
  #endif
            optc = -1; // special value indicating no segments used
          }
#else
	// see comment CBA above
#endif

      }
  }


  { COvlps *c;

    // over all segments involving the first A contig
    for (c = ABuckets[0]; c != NULL; c = c->Alink)
      // if ahang is negative, the overlap starts along the B edge, so it is a starting point
      if (c->seg->overlap->begpos <= 0)
        { score = c->seg->overlap->length;
          if (score > c->best)
            { c->best  = score;
              c->trace = NULL;
            }
#ifdef DEBUG_ALIGN
          fprintf(stderr,"  Start path (%d,%d) = %d\n",
                 c->seg->a_contig,c->seg->b_contig,c->best);
#endif
        }
    // over all segments involving the first B contig
    for (c = BBuckets[0]; c != NULL; c = c->Blink)
      // if ahang is positive, the overlap starts along the A edge, so it is a starting point
      if (c->seg->overlap->begpos >= 0)
        { score = c->seg->overlap->length;
          if (score > c->best)
            { c->best  = score;
              c->trace = NULL;
            }
#ifdef DEBUG_ALIGN
           fprintf(stderr,"  Start path (%d,%d) = %d\n",
                 c->seg->a_contig,c->seg->b_contig,c->best);
#endif
        }
  }


  {
    int ca;
    for( ca=0; ca<=AF->num_gaps;ca++){
      COvlps *co = ABuckets[ca];
      while(co!=NULL){
	Segment *s = co->seg;
	int pnt;

	if(co->Alink!=NULL){
	  Segment *t = co->Alink->seg;
	  assert(s->a_contig==t->a_contig);
	  assert(s->b_contig<=t->b_contig);
	}


#ifdef DEBUG_ALIGN
	fprintf(stderr,"working on ctg[%d,%d] score %d\n",s->a_contig,s->b_contig,co->best);
#endif

        term = 0;

	// if this is either along the starting boundary or is an internal segment that has
	// already been reached from starting in a gap ...
        if ((score = co->best) >= 0){
            if (s->overlap->endpos < 0)
              { pnt = AF->ctgs[s->a_contig].lft_end +
                      (AF->ctgs[s->a_contig].length + s->overlap->endpos); 
#ifdef DEBUG_ALIGN
                fprintf(stderr,"  Point ray B%d at %d, source (%d,%d) = %d\n",
                       s->b_contig,pnt,s->a_contig,s->b_contig,score);
#endif
                term = Link_Vertical(AF, BF, varwin,
                                     s->a_contig, s->b_contig,
                                     pnt, pnt, co);
              }
            else 
              { pnt = BF->ctgs[s->b_contig].lft_end +
                      (BF->ctgs[s->b_contig].length - s->overlap->endpos); 
#ifdef DEBUG_ALIGN
                fprintf(stderr,"  Point ray A%d at %d, source (%d,%d) = %d\n",
                       s->a_contig,pnt,s->a_contig,s->b_contig,score);
#endif
                term = Link_Horizontal(AF, BF, varwin,
                                       s->a_contig, s->b_contig,
                                       pnt, pnt, co);
              }
          }
        if (term && score > mval)
          { mval = score;
#ifdef STATS
            edgecount += 1;
#endif

            optco = co;
#ifdef DEBUG_OPTTAIL
	    fprintf(stderr, "Best path so far ends at ctgs[%d,%d} score  %d\n",
		    optco->seg->a_contig,
		    optco->seg->b_contig,
		    mval);
#endif
          }
	co=co->Alink;
        }

      }
  }



  if (mval <= 0)
    { seglist = NULL;
#ifdef DEBUG_ALIGN
      fprintf(stderr,"  NO OVERLAP\n");
#endif

      // remember to NOT free seglist members here 
      // --- they will be used in Compare_Scaffold if this fcn returns NULL


      // BUT IF WE ARE NOT USING COMPARE_SCAFFOLD, OR OTHER CODE THAT KEEPS
      // THE ORIGINAL POINTER, 
      // I.E. cgw InterleavedMerging CODE,
      // IS THIS A MEMORY LEAK? -- ALH

    }
  else
    { Segment *s, *r;
      COvlps  *c;

      for(optc=0;optc<numsegs;optc++){
	if(CtgOvls+optc==optco)break;
      }
      assert(optc<numsegs);
#ifdef DEBUG_OPTTAIL
      fprintf(stderr,"Optimal path ends at ctgs[%d,%d] score %d\n",
	      CtgOvls[optc].seg->a_contig,
	      CtgOvls[optc].seg->b_contig,
	      CtgOvls[optc].best);
#endif
      assert(optc>=0);

      for (c = CtgOvls + optc; c != NULL; c = c->trace)
	c->seg->alow = - (c->seg->alow+1);

      for (s = seglist; s != NULL; s = r)
        { r = s->next;
          if (s->alow >= 0)
            { free(s->overlap);
              free(s);
            }
          else
            s->alow = - (s->alow+1);
        }

      r = NULL;
      for (c = CtgOvls + optc; c != NULL; c = c->trace)
	{ s = c->seg;
	s->next = r;
	r = s;
	}

      seglist = r;

#ifdef DEBUG_ALIGN
      fprintf(stderr,"  Best overlap ends at (%d,%d) = %d\n",
             CtgOvls[optc].seg->a_contig,CtgOvls[optc].seg->b_contig,
             CtgOvls[optc].best);
#endif
    }

#ifdef XFIG
  fclose(figfile);
#endif

#ifdef STATS
  vertcount  += numsegs+2;
  boundcount += AF->num_gaps + BF->num_gaps;
#endif

  *best = mval;
  return (seglist);
}
  

Segment *Compare_Scaffold(Segment *seglist, int numsegs,
                                 Scaffold *AF, Scaffold *BF, int *best)
{ static int *vector = NULL;
  static int  maxvec = -1;

  int      nacs, nbcs;
  Segment *s, *r, *t;
  int     *whom, *trace;
  int lastrow, lastcol, lastadd = 0;
  int segno;
  int j, max, opt;

  nacs = AF->num_gaps+1;
  nbcs = BF->num_gaps+1;

  if (numsegs + 2*nbcs > maxvec)
    { maxvec = (int)(1.2*(numsegs + 2*nbcs) + 500);
      vector = (int *) realloc(vector,sizeof(int)*maxvec);
      if (vector == NULL)
        { fprintf(stderr,"Out of memory allocating DP array\n");
          exit (1);
        }
    }
  whom = vector + nbcs;
  trace = whom + nbcs;
  
  for (j = 0; j <= BF->num_gaps; j++)
    { vector[j] = 0;
      whom[j] = -1;
    }

  lastrow = nacs+1;
  lastcol = nbcs-1;
  segno = 0;
  r = NULL;
  for (s = seglist; s != NULL; s = t)
    { if (s->a_contig != lastrow || s->b_contig != lastcol)
        { for (j = lastcol-1; j >= s->b_contig; j--)
            if (vector[j+1] > vector[j])
              { vector[j] = vector[j+1];
                whom[j] = whom[j+1];
              }
          lastcol = s->b_contig;
          lastrow = s->a_contig;
          lastadd = s->overlap->length;
          vector[lastcol] = lastadd + vector[lastcol]; 
          trace[segno] = whom[lastcol];
          whom[lastcol] = segno;
        }
      else
        { if (s->overlap->length > lastadd)
            { vector[lastcol] += (s->overlap->length - lastadd);
              lastadd = s->overlap->length;
              trace[segno] = trace[whom[lastcol]];
              whom[lastcol] = segno;
            }
          else
            trace[segno] = -1;
        }

      segno += 1;
      t = s->next;
      s->next = r;
      r = s;

#ifdef DEBUG_COMPARE
      fprintf(stderr,"  (a,b,v) = (%d,%d,%d)\n     ",lastrow,lastcol,lastadd);
      for (j = 0; j <= BF->num_gaps; j++)
        fprintf(stderr," %d(%d)",vector[j],whom[j]);
      fprintf(stderr,"\n");
#endif
    }
  seglist = r;

#ifdef DEBUG_COMPARE
  fprintf(stderr,"\nTraceback:\n");
  for (j = 0; j < numsegs; j++)
    fprintf(stderr," %3d: -> %d\n",j,trace[j]);
#endif

  max = vector[0];
  opt = whom[0];
  for (j = 1; j <= BF->num_gaps; j++)
    if (vector[j] > max)
      { max = vector[j];
        opt = whom[j];
      }

  segno = numsegs-1;
  r = NULL;
  for (s = seglist; s != NULL; s = t)
    { t = s->next;
      if (segno == opt)
        { s->next = r;
          r = s;
          opt = trace[opt];
        }
      else
        { free(s->overlap);
          free(s);
        }
      segno -= 1;
    }
  seglist = r;

  r = NULL;
  for (s = seglist; s != NULL; s = t)
    { t = s->next;
      s->next = r;
      r = s;
    }
  seglist = r;

  *best = max;
  return (seglist);
}

void Free_Scaffold_List(Scaffold_List *SL)
{ Scaffold_List *s, *t;

  for (s = SL; s != NULL; s = t)
    { t = s->next;
      Free_Scaffold(s->scaffold);
      free(s);
    }
}

static int UnAccountedBits(int alow, int ahgh, int blow, int bhgh,
                           int i, int p, int j, int q, 
                           Scaffold *AF, Scaffold *BF)
{ int abg, abe;
  int u, v;
  int delta;

  delta = 0;
  if (alow < ahgh)
    { abg = alow; abe = ahgh; }
  else
    { abg = ahgh; abe = alow; }

#ifdef DEBUG_ANALYSIS
  fprintf(stderr,"  Checking gap [%d,%d] vs [%d,%d]\n",ahgh,alow,bhgh,blow);
  fprintf(stderr,"     Contigs   [%d,%d] vs [%d,%d]\n",p,i,q,j);
#endif

  for (u = p; u <= i; u++)
    { int b, e;
      int x, y;

      b = AF->ctgs[u].lft_end;
      e = b + AF->ctgs[u].length;
      if (abg > b) b = abg;
      if (abe < e) e = abe;
#ifdef DEBUG_ANALYSIS
      fprintf(stderr,"       X [%d,%d] = [%d,%d]\n",
             AF->ctgs[u].lft_end,AF->ctgs[u].lft_end+AF->ctgs[u].length,b,e);
#endif
      if (b >= e && abg != abe) continue;
      if (abg == abe)
        { x = bhgh;
          y = blow;
        }
      else
        { x = (int)(((blow - 1.*bhgh) / (alow - ahgh)) * (b - ahgh) + bhgh);
          y = (int)(((blow - 1.*bhgh) / (alow - ahgh)) * (e - ahgh) + bhgh);
        }
#ifdef DEBUG_ANALYSIS
      fprintf(stderr,"       yields [%d,%d]\n",x,y);
#endif
      if (x > y)
        { int m; m = x; x = y; y = m; }
      for (v = q; v <= j; v++)
        { int s, t;

          s = BF->ctgs[v].lft_end;
          t = s + BF->ctgs[v].length;
          if (s < x) s = x;
          if (t > y) t = y;
#ifdef DEBUG_ANALYSIS
          fprintf(stderr,"       X [%d,%d] = [%d,%d] of size %d\n",
                 BF->ctgs[v].lft_end,
                 BF->ctgs[v].lft_end+BF->ctgs[v].length,
                 s,t,t-s);
#endif
          if (s < t)
            delta += t-s;
        }
    }

  return (delta);
}

Scaffold_Overlap *Analyze_Overlap(Scaffold *AF, Scaffold *BF,
                                         Segment *seglist, int score, int comp)
{ int i, j;
  int diag, dmax = 0, dmin = 0;
  int firstime, lastime;
  int alow, ahgh = 0;
  int blow, bhgh = 0;
  int adelta, bdelta;
  int aone, bone;
  int alst, blst;
  int totdiff, totleng, segdiff;
  Segment *seg;

  Scaffold_Overlap *sovl;

  seg = seglist;
  aone = seg->a_contig;
  bone = seg->b_contig;

  segdiff = 0;
  totdiff = 0;
  totleng = 0;
  adelta = bdelta = 0;
  firstime = 1;
  lastime  = 0;
  i = j = 0;
  while (seg != NULL || ! lastime)
    { int p, q;

      alst = p = i;
      blst = q = j;
      if (seg == NULL)
        { i = AF->num_gaps;
          j = BF->num_gaps;
          if (AF->length - ahgh < BF->length - bhgh)
            { blow = bhgh + (AF->length-ahgh);
              alow = AF->length;
            }
          else
            { alow = ahgh + (BF->length-bhgh);
              blow = BF->length;
            }
        } 
      else
        { i = seg->a_contig;
          j = seg->b_contig;
          if (seg->overlap->begpos > 0)
            { alow = AF->ctgs[i].lft_end + seg->overlap->begpos;
              blow = BF->ctgs[j].lft_end;
            }
          else
            { alow = AF->ctgs[i].lft_end;
              blow = BF->ctgs[j].lft_end - seg->overlap->begpos;
            }
          totdiff += seg->overlap->diffs;
          totleng += seg->overlap->length;

          if (comp)
            { int m;
#ifdef DPC_BASED
              seg->overlap->comp = 1;
#endif
              m = BF->ctgs[j].length - seg->bhgh;
              seg->bhgh = BF->ctgs[j].length - seg->blow;
              seg->blow = m;
            }
        }

      if (firstime)
        { if (alow < blow)
            { ahgh = 0;
              bhgh = blow - alow;
            }
          else
            { ahgh = alow - blow;
              bhgh = 0;
            }
         }

      { int a_bits, b_bits;

        a_bits = UnAccountedBits(blow,bhgh,alow,ahgh,j,q,i,p,BF,AF);
        b_bits = UnAccountedBits(alow,ahgh,blow,bhgh,i,p,j,q,AF,BF);
        adelta += a_bits;
        bdelta += b_bits;
        if (a_bits > b_bits)
          segdiff += a_bits;
        else
          segdiff += b_bits;
      }

      if (seg != NULL)
        { if (seg->overlap->endpos > 0)
            { bhgh = BF->ctgs[j].lft_end
                   + (BF->ctgs[j].length - seg->overlap->endpos);
              ahgh = AF->ctgs[i].lft_end + AF->ctgs[i].length;
            }
          else
            { bhgh = BF->ctgs[j].lft_end + BF->ctgs[j].length;
              ahgh = AF->ctgs[i].lft_end
                   + (AF->ctgs[i].length + seg->overlap->endpos);
            }

          diag = ( (alow - blow) + (ahgh - bhgh) ) / 2;
          if (firstime || diag > dmax)
            dmax = diag;
          if (firstime || diag < dmin)
            dmin = diag;

#ifdef DEBUG_ANALYSIS
          fprintf(stderr,"  diag = %d\n",diag);
          fprintf(stderr,"  [%d,%d]  [%d,%d]\n",alow,ahgh,blow,bhgh);
#endif
        }

      firstime = 0;
      if (seg == NULL)
        lastime = 1;
      else
        seg = seg->next;
    }

  { int      i, diag;
    double   slack;
    Overlap *ovl;

    slack = 0.;
    for (i = aone; i < alst; i++)
      slack += AF->gaps[i].gap_var;
    for (i = bone; i < blst; i++)
      slack += BF->gaps[i].gap_var;

    diag = (dmin + dmax) / 2;

    ovl = (Overlap *) malloc(sizeof(Overlap));
    ovl->begpos = diag;
    ovl->endpos = diag - (AF->length - BF->length);
    ovl->diffs  = totdiff + segdiff;
    ovl->length = (AF->length + BF->length -
                        (abs(diag) + abs(ovl->endpos))) / 2;
    ovl->comp   = comp;
    //    ovl->aseq   = ovl->bseq = NULL;
    //    ovl->trace  = NULL;

    sovl = (Scaffold_Overlap *) malloc(sizeof(Scaffold_Overlap));
    sovl->score   = score;
    sovl->erate   = (1.*ovl->diffs) / ovl->length;
    sovl->ascaf   = AF;
    sovl->bscaf   = BF;
    sovl->D_delta = dmax-dmin;
    sovl->D_var   = (int)slack;
    sovl->A_delta = adelta;
    sovl->B_delta = bdelta;
    sovl->overlap = ovl;
    sovl->seglist = seglist;
  }

  return (sovl);
}

Scaffold_Overlap *Compare_Scaffolds(Scaffold_List *AS, Scaffold_List *BS)
{ int i, ip;
  int j, jp;
  int tab;
  Local_Pool       *pool;
  Scaffold_Overlap *list, *list2, *fing;
  Scaffold_List    *a, *b;

  pool = Find_All_Locals(AS,BS);
  
  i  = 0;
  ip = 0;
  list = fing = NULL;
  tab  = 0;
  for (a = AS; a != NULL; a = a->next)
    { j  = 0;
      jp = 0;
      for (b = BS; b != NULL; b = b->next)
        { Segment *seglist, *scflist = NULL;
          int      vw, score, confirmed, numsegs;
	  int first=1;

          seglist = Find_All_Overlaps(a->scaffold,b->scaffold,pool,
                                      i,j,ip,jp,&tab,&numsegs,0);

          for (vw = 0; VarWindows[vw] != 0; vw++)
            { 
#ifdef AHANG_BAND_TEST
	      scflist = Align_Scaffold(seglist,numsegs,VarWindows[vw],
                                       a->scaffold,b->scaffold,&score,
				       -b->scaffold->length,
				       a->scaffold->length);
#else
	      scflist = Align_Scaffold(seglist,numsegs,VarWindows[vw],
                                       a->scaffold,b->scaffold,&score);
#endif

	      if(score==0&&scflist==NULL&&first){
		fprintf(stderr, "It's possible to interleave scaffolds %d and %d with no overlapping contigs at %d sigma\n",i,j,VarWindows[vw]);
		first=0;
	      }


              if (scflist != NULL) break;
            }
          if (scflist != NULL)
            { seglist = scflist;
              confirmed = VarWindows[vw];
            }
          else
            { seglist = Compare_Scaffold(seglist,numsegs,
                                         a->scaffold,b->scaffold,&score);
              confirmed = 0;
            }
          if (seglist != NULL)
            { Scaffold_Overlap *overlap;

              overlap = Analyze_Overlap(a->scaffold,b->scaffold,
                                        seglist,score,0);
              if (list == NULL)
                list = overlap;
              else
                fing->next = overlap;
              fing = overlap;

              fing->asnum = i;
              fing->bsnum = j;
              fing->abase = ip;
              fing->bbase = jp;
              fing->firm  = confirmed;

              { Segment *s;
                for (s = seglist; s != NULL; s = s->next)
                  { s->a_contig += ip;
                    s->b_contig += jp;
                  }
              }
            }

          j  += 1;
          jp += (b->scaffold->num_gaps+1);
        }
      i  += 1;
      ip += (a->scaffold->num_gaps+1);
    }

  if (list != NULL)
    fing->next = NULL;

  for (b = BS; b != NULL; b = b->next)
    Complement_Scaffold(b->scaffold);

  i  = 0;
  ip = 0;
  list2 = fing;
  tab   = 0;
  for (a = AS; a != NULL; a = a->next)
    { j  = 0;
      jp = 0;
      for (b = BS; b != NULL; b = b->next)
        { Segment *seglist, *scflist = NULL;
          int      vw, score, confirmed, numsegs;
	  int first=1;

          seglist = Find_All_Overlaps(a->scaffold,b->scaffold,pool,
                                      i,j,ip,jp,&tab,&numsegs,1);

          for (vw = 0; VarWindows[vw] != 0; vw++)
            { 
#ifdef AHANG_BAND_TEST
	      scflist = Align_Scaffold(seglist,numsegs,VarWindows[vw],
                                       a->scaffold,b->scaffold,&score,
				       -b->scaffold->length,
				       a->scaffold->length);
#else
	      scflist = Align_Scaffold(seglist,numsegs,VarWindows[vw],
                                       a->scaffold,b->scaffold,&score);
#endif

	      if(score==0&&scflist==NULL&&first){
		fprintf(stderr, "It's possible to interleave scaffolds %d and rev_comp(%d) with no overlapping contigs at %d sigma\n",i,j,VarWindows[vw]);
		first=0;
	      }

              if (scflist != NULL) break;
            }
          if (scflist != NULL)
            { seglist = scflist;
              confirmed = VarWindows[vw];
            }
          else
            { seglist = Compare_Scaffold(seglist,numsegs,
                                         a->scaffold,b->scaffold,&score);
              confirmed = 0;
            }
          if (seglist != NULL)
            { Scaffold_Overlap *overlap;

              overlap = Analyze_Overlap(a->scaffold,b->scaffold,
                                        seglist,score,1);
              if (list == NULL)
                list = overlap;
              else
                fing->next = overlap;
              fing = overlap;

              fing->asnum = i;
              fing->bsnum = j;
              fing->abase = ip;
              fing->bbase = jp;
              fing->firm  = confirmed;

              { Segment *s;
                for (s = seglist; s != NULL; s = s->next)
                  { s->a_contig += ip;
                    s->b_contig  = jp + b->scaffold->num_gaps - s->b_contig;
                  }
              }
            }

          j  += 1;
          jp += (b->scaffold->num_gaps+1);
        }
      i  += 1;
      ip += (a->scaffold->num_gaps+1);
    }

  if (list != NULL)
    fing->next = NULL;

  for (b = BS; b != NULL; b = b->next)
    Complement_Scaffold(b->scaffold);

  free(pool->address);
  free(pool->locals);

  if (list2 != NULL && list2->next != NULL)
    { Scaffold_Overlap *x, *y, *n;

      y = list2 = list2->next;
      for (x = list; x != list2; x = x->next)
        { while (y != NULL && (y->asnum < x->asnum || (y->asnum == x->asnum &&
                                                       y->bsnum  < x->bsnum)) )
            y = y->next;
          if (y != NULL && y->asnum == x->asnum && y->bsnum == x->bsnum)
            { int single_rev;
	      Local_Overlap *o;
              single_rev = 0;

	      /* revisions by ALH 12/14/00: 
		 single_rev must consider possibility of y
		 being a single reversed segment as well as x being so;
		 originally, a containment test was applied as well, but
		 this isn't really logically required in order for an
		 overlap consisting of a single reversed segment to be
		 rejected.

		 The logic is ... if an overlap in one direction consists
		 of a single reversed segment, then the overlap in the
		 other orientation had the potential to consist of the
		 very same segment (in forward orientation).  Either that
		 is the case (in which case the original, reversed overlap
		 should be rejected), or the second overlap contains the
		 single segment of the first overlap (in which case, the
		 second overlap was considered to be superior to the single
		 segment alone, and the first, reversed overlap should be
		 rejected) OR the forward overlap doesn't contain the reversed
		 segment BUT was chosen OVER any overlap containing that
		 segment (and thus, we should once again reject the overlap
		 consisting of the single, reversed segment).

	      */

	      o = x->seglist->overlap;
	      if (x->seglist->next == NULL&& o->num_pieces == 1 && o->chain[0].reversed)
		{ single_rev = 1; }
	      else 
		{ o = y->seglist->overlap;
		  if(y->seglist->next == NULL&&o->num_pieces == 1 && o->chain[0].reversed)
		    { single_rev = 2; }
		}

	      if(single_rev>0){
		if(single_rev==1){
		  x->abase = -1;
		} else {
		  y->abase = -1;
		}
	      } else {
		if (y->erate < x->erate ){
		  x->abase = -1;
		} else {

		  /* Revision, ALH, 1/5/01:
		     Additional testing in the case of ties in erate;
		     I came across a case where the forward and reverse
		     overlaps consisted of the same sets of segments,
		     but the statistics were quite different, since all the
		     segments were forward in one overlap and reverse
		     in the other.  Obviously, we prefer that version of
		     the overlap which agrees with the orientation of the
		     segments--and this seems to be the one with the better
		     statistics (better deltas, better confirmed level, etc).
		     In the one known case of this, any of the measures
		     other than erate would have given the right result;
		     I chose to test on D_delta relatively arbitrarily.
		  */
		  if (y->D_delta < x->D_delta) {
		    x->abase = -1;
		  } else {
		    y->abase = -1;
		  }
		}
	      }

            }


        }

      y = NULL;
      for (x = list; x != NULL; x = n)
        { n = x->next;
          if (x->abase < 0)
            { if (y == NULL)
                list = n;
              else
                y->next = n;
              Free_Segments_ScafComp(x->seglist);
              free(x->overlap);
              free(x);
            }
          else
            y = x;
        }
    }

  return (list);
}


/* These are the lines that the AnnoGraph edge_command script parses
   to appropriately highlight the celamy files for display 
   Note that the contig ids are incremented, since the celamy
   view starts numbering at 1.                                 

   REPORT a_ctg_id [low,high] ib_ctg_id [low,high]                  
*/

void Print_Anno_Info(Scaffold_Overlap *sovl)
{ printf("\nAnnograp Info:\n\n");
  for (; sovl != NULL; sovl = sovl->next)
    { Segment  *s;

      for (s = sovl->seglist; s != NULL; s = s->next)
        printf("REPORT %d [%d,%d] %d [%d,%d]\n",
               s->a_contig+1,s->alow,s->ahgh,s->b_contig+1,s->blow,s->bhgh);
    }
}

void Free_Scaffold_Overlaps(Scaffold_Overlap *SO)
{ Scaffold_Overlap *s, *t;

  for (s = SO; s != NULL; s = t)
    { t = s->next;
      Free_Segments_ScafComp(s->seglist);
      free(s->overlap);
      free(s);
    }
}



