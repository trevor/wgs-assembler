
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
static char CM_ID[]= "$Id: AS_MSG_pmesg.c,v 1.31 2007-02-08 06:48:53 brianwalenz Exp $";

//  reads old and new AFG message (with and w/o chaff field)
#define AFG_BACKWARDS_COMPATIBLE
#define IAF_BACKWARDS_COMPATIBLE

//  FreeBSD 6.1 fgets() sporadically replaces \n with \0, which
//  horribly breaks this reader.  Defined this to replace
//  fgets() with fgetc().
#undef FGETS_IS_BROKEN

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <ctype.h>
#include <assert.h>

#include "AS_global.h"

// time_t field formats
#define CRT_FORMAT   "crt:" F_TIME_T
#define CTM_FORMAT   "ctm:" F_TIME_T
#define ETM_FORMAT   "etm:" F_TIME_T

// action field formats
#define ACT_FORMAT   "act:%c"
#define ACT1_FORMAT  "act:%1"

// type field formats
#define TYP_FORMAT   "typ:%c"
#define TYP1_FORMAT  "typ:%1"

// orientation field formats
#define ORI_FORMAT   "ori:%c"
#define ORI1_FORMAT  "ori:%1"

// floating point formats
#define MEA_IN_FORMAT   "mea:%f"
#define MEA_OUT_FORMAT  "mea:%.3f"
#define STD_IN_FORMAT   "std:%f"
#define STD_OUT_FORMAT  "std:%.3f"
#define VAR_FORMAT      "var:%f"
#define QUA_FORMAT      "qua:%f"
#define RAT_IN_FORMAT   "rat:%f"
#define RAT_OUT_FORMAT  "rat:%.3f"

//position formats
#define POS1_FORMAT  "pos:" F_COORD
#define POS2_FORMAT  "pos:" F_COORD "," F_COORD
#define CLR_FORMAT   "clr:" F_COORD "," F_COORD

#define MINC_FORMAT "min:" F_COORD
#define MAXC_FORMAT "max:" F_COORD

#define EFRAG1_FORMAT "fg1:" F_UID
#define EFRAG2_FORMAT "fg2:" F_UID
#define IFRAG1_FORMAT "fg1:" F_IID
#define IFRAG2_FORMAT "fg2:" F_IID

#define EACC_FORMAT  "acc:" F_UID
#define EBTG_FORMAT  "btd:" F_UID
#define EDST_FORMAT  "dst:" F_UID
#define ELEN_FORMAT  "len:" F_UID
#define ELOC_FORMAT  "loc:" F_UID
#define EREF_FORMAT  "ref:" F_UID
#define ESEQ_FORMAT  "sid:" F_UID

#define F_UID_IID    "(" F_UID "," F_IID ")"

#define IACCS_FORMAT "acc:" F_UID_IID
#define IDST_FORMAT  "dst:" F_IID
#define ILENS_FORMAT "len:" F_UID_IID
#define ILOCS_FORMAT "loc:" F_UID_IID
#define IREFS_FORMAT "ref:" F_UID_IID
#define ISEQS_FORMAT "sid:" F_UID_IID

#define ZERO 0
#define MAX_LINE_LEN (128 * 1024)    /* Maximum input line length (checked) */

#define ROUNDUP(n,u) ((((n)-1)/(u) + 1)*(u))  /* Round n up to nearest
                                                 multiple of u */

static int LineNum = 0;   /* Current line number */

static const char *Mcode;       /* 3-code of current read/write routine */

static char *MemBuffer = NULL;      /* Memory allocation buffer for messages */
static int   MemMax = -1, MemTop;   /* Memory ceiling and current top */

int novar = 0;  /*  Output or not the variation records */





////////////////////////////////////////
//
//  For lack of a better place, this is here...
//
//
/* ---------------------------------------------------
 * Convert values of type: char
 * to the enumerated type: FragType
 *
 * This violates the abstraction provided by enumerated types.
 * But our storage requirements demand this byte-saving step.
 *
 * Param 'strict' controls behavior on bad input.
 * If strict, function dies when input is not valid.
 * Else, function returns AS_EXTR when input is not valid.
   ----------------------------------------------------*/
FragType AS_MSG_SafeConvert_charToFragType (const char input, bool strict) {
  FragType output = AS_EXTR;

  // The vast majority of data will satisfy the first test.

  if      (input=='R') output=AS_READ;
  else if (input=='X') output=AS_EXTR;
  else if (input=='T') output=AS_TRNR;
  else if (strict)     assert(FALSE);

  return output;
}
//
//
//
////////////////////////////////////////


/* Make sure there is a block of size bytes left in memory buffer starting
   at an index that is a multiple of boundary.  Return the *index* into the
   array, so that the realloc does not blow structures in the process of
   being built.  All pointers in such structures are saved as integers and
   converted to pointers after all allocation has taken place.             */

static void MakeSpace(const int size){
  size_t newsize=1;
  char *newbufr;

  if(MemMax> size)
    return;

  if(MemMax < 0)
    newsize = 2 * 2048 * 2048; // This may be excessive, but this code is BRITTLE!
  else
    newsize = 2 * size;

  newsize = ROUNDUP(newsize,8);
  newbufr = (char *)safe_realloc(MemBuffer,newsize);
  MemBuffer = newbufr;
  MemMax    = newsize;
}


static long MoreSpace(const int size, const int boundary)
{ 
  MemTop = ROUNDUP(MemTop,boundary);

  MakeSpace(MemTop + size);
  { int alloc;

    alloc   = MemTop;
    MemTop += size;
    return (alloc);
  }
}

static char  CurLine[MAX_LINE_LEN];     /* Line buffer for reading messages. */

static char *ReadLine(FILE *fin) {

  CurLine[MAX_LINE_LEN-2] = '\n';
  CurLine[MAX_LINE_LEN-1] = 0;

  errno = 0;

  LineNum++;

#ifdef FGETS_IS_BROKEN
  int p=0;
  for (p=0; p<MAX_LINE_LEN-1; p++) {
    CurLine[p] = fgetc(fin);
    if (CurLine[p] == 0)
      CurLine[p] = '\n';
    if (CurLine[p] == '\n') {
      CurLine[p+1] = 0;
      break;
    }
  }
#else
  if (fgets(CurLine, MAX_LINE_LEN-1, fin) == NULL) {
    fprintf(stderr,"ERROR: AS_MSG_pmesg.c::ReadLine()-- Premature end of input at line %d (%s)\n", LineNum, Mcode);
    fprintf(stderr,"       '%s'\n", CurLine);
    exit(1);
  }
#endif

  if (errno) {
    fprintf(stderr,"ERROR: AS_MSG_pmesg.c::ReadLine()-- Read error at line %d: '%s'\n", LineNum, strerror(errno));
    fprintf(stderr,"       '%s'\n", CurLine);
    exit(1);
  }

  if (CurLine[MAX_LINE_LEN-2] != '\n') {
    fprintf(stderr,"ERROR: Input line %d is too long (%s)\n", LineNum, Mcode);
    fprintf(stderr,"       '%s'\n", CurLine);
    exit(1);
  }

  //fprintf(stderr, "READLINE -- %d %d %s", CurLine[0], CurLine[1], CurLine);

  return(CurLine);
}


static char *GetLine(FILE *fin, int skipComment)   /* Get next input line (there must be one). */
{ 
    do {
        ReadLine(fin);
    }
    while (skipComment && CurLine[0] == '#');
    return (CurLine);
}

   /* Found an enum out-of-range error. */

static void MtypeError(const char * const name)
{   
    fprintf(stderr,"ERROR: Illegal %s type value \"%c\" (%s) at line %d\n",
            name,CurLine[4],Mcode, LineNum);
    exit (1);
}

   /* Found a bad 3-code field name. */

static void MtagError(const char * const tag)
{ 
    fprintf(stderr,"ERROR: Illegal tag \"%s\" (expected \"%s\") (%s) at line %d\n",
            CurLine, tag, Mcode, LineNum);
    exit (1);
}

   /* Field content area did not have correct syntax. */

static void MfieldError(const char * const mesg)
{ 
    int len;

    len = strlen(CurLine)-1;
    if (CurLine[len] == '\n')
        CurLine[len] = '\0';
    fprintf(stderr,"ERROR: %s \"%s\" (%s) at line %d\n",
            mesg,CurLine,Mcode,LineNum);
    exit (1);
}

   /* General error message exit. */

static void MgenError(const char * const mesg)
{ 
    fprintf(stderr,"ERROR: %s (%s) at line %d\n",
            mesg,Mcode,LineNum);
    exit (1);
}

   /* Get a text field item: syntax "tag\n(%s\n)*.\n".
      Tricky part is length is not known ahead of time.  */

static long GetText(const char * const tag, FILE *fin, const int delnewlines)
{ 
  long text, idx; 
  int len;
 
  if (strncmp(GetLine(fin, TRUE),tag,4) != 0)
      MtagError(tag);
  text = MemTop;
  while (1)
  { 
      char *line = GetLine(fin, FALSE);
      if ((line[0] == '.') && (line[1] == '\n'))
        break;
      len = strlen(CurLine);
      if (delnewlines && CurLine[len-1] == '\n') len -= 1;
      idx = MoreSpace(len,1);
      strncpy(MemBuffer+idx,CurLine,len);
  }
  idx = MoreSpace(1,1);
  MemBuffer[idx] = '\0';
  return (text);
}

   /* Get a string field item: syntax "tag%s\n".
      Tricky part is length is not known ahead of time.  */

static long GetString(const char * const tag, FILE *fin)
{ char *str;
  int   eos, len, text, idx;

  errno = 0;

  text = MemTop;
  do
    {
      ReadLine(fin);
    }
  while (CurLine[0] == '#');
  if (strncmp(CurLine,tag,4) != 0)
    MtagError(tag);
  str = CurLine + 4;

  while (1)
    { len = strlen(str);
      eos = (str[len-1] == '\n');
      if (eos) len -= 1;
      idx = MoreSpace(len,1);
      strncpy(MemBuffer+idx,str,len);
      if (eos) break;
      str = ReadLine(fin);
    }
  idx = MoreSpace(1,1);
  MemBuffer[idx] = '\0';
  return (text);
}

   /* Output text field item with 3-code field-name "tag". */

static void PutText(FILE *fout, const char * const tag, 
		     char * text, const int format)
{
    // Note that the data of "text" is modified!!!
    int i, len;
    fprintf(fout,"%s\n",tag);
    if ((text != NULL) && (text[0] != 0)) {
        len = strlen(text);
        if (format) { 
            for (i = 0; i < len; i += 70)
            { 
                fprintf(fout,"%.*s\n",70,text);
                text += 70;
            }
            fprintf(fout,".\n");
        } else{ 
            if (text[len-1] == '\n')     /* Strip trailing new line if prez. */
            text[len-1] = '\0';
            fprintf(fout,"%s\n.\n", text);
        }
    } else{
        fprintf(fout,".\n");
    }
}

/* Macros to get each type of 3-code field: obvious combo's of above */

#define GET_TYPE(lvalue,format,name)		\
{ char value[2];				\
  if (sscanf(GetLine(fin, TRUE),format,value) != 1)	\
    MtypeError(name);				\
  lvalue =  *value;				\
}

/** DO NOT USE FOR SETTING ENUM TYPES FROM CHARACTER FLAGS!!!! **/
#define GET_FIELD(lvalue,format,emesg)			\
{ if (sscanf(GetLine(fin, TRUE),format,&(lvalue)) != 1)	\
    MfieldError(emesg);					\
}

#define GET_PAIR(lvalue1,lvalue2,format,emesg)			\
{ if (sscanf(GetLine(fin,TRUE),format,&(lvalue1),&(lvalue2)) != 2)	\
    MfieldError(emesg);						\
}

// lvalue3 is _ALWAYS_ type char coupled with a %1 format
#define GET_TRIPLE(lvalue1,lvalue2,lvalue3,format,emesg)	           \
{                                                                          \
  char value[2];                                                           \
  if(sscanf(GetLine(fin,TRUE),format,&(lvalue1),&(lvalue2), value) != 3)   \
    MfieldError(emesg);						           \
  lvalue3 = *value;                                                        \
}

#define GET_EOM					\
{ if (GetLine(fin,TRUE)[0] != '}')			\
    MgenError("Expecting end of message");	\
}

/******************** INPUT ROUTINES ***************************/

/*  Routine to input each type of proto-IO message. */


static void *Read_Dist_Mesg(FILE *fin, int external)
{ static InternalDistMesg dmesg;

  if (external)
    { char ch;
      GET_TYPE(ch,ACT_FORMAT,"action");
      dmesg.action = (ActionType) ch;
      GET_FIELD(dmesg.eaccession,EACC_FORMAT,"accession field")
    }
  else
  { char ch;
      GET_TYPE(ch,ACT1_FORMAT "[ADR]","action");
      dmesg.action = (ActionType) ch;
      GET_PAIR(dmesg.eaccession,dmesg.iaccession,
               IACCS_FORMAT,"accession field pair");
    }
  if (dmesg.action == AS_ADD ||
      dmesg.action == AS_REDEFINE)
    { GET_FIELD(dmesg.mean,MEA_IN_FORMAT,"mean field");
      GET_FIELD(dmesg.stddev ,STD_IN_FORMAT,"stddev field");
    }
  GET_EOM;
  return ((void *) (&dmesg));
}

static void *Read_DST_Mesg(FILE *fin)
{ return Read_Dist_Mesg(fin,1); }

static void *Read_IDT_Mesg(FILE *fin)
{ return Read_Dist_Mesg(fin,0); }

// This is tricky, since we are trying to update the fields
// of a structure that may be realloced.  Store offsets in
// temps, and stuff them in in one go.

static int Read_ADL_Struct(long last, FILE *fin)
{ AuditLine mesg; // temporary
  long crnt;
  AuditLine *newMesg;

  mesg.name    = (char *)GetString("who:",fin);
  GET_FIELD(mesg.complete,CTM_FORMAT,"completion field");
  mesg.version = (char *) GetString("vsn:",fin);
  mesg.comment = (char *)GetText("com:",fin, FALSE);
  mesg.next = (AuditLine *)last;
  // this looks SUSPICIOUS
  GET_EOM;

  // Now allocate space for this guy, and copy him to the allocated space
  crnt = MoreSpace(sizeof(AuditLine),8);
  newMesg = (AuditLine *)(MemBuffer + crnt);
  *newMesg = mesg;

  return crnt;
}

static void *Read_ADT_Mesg(FILE *fin)
{ static AuditMesg amesg;
  AuditLine *cptr, *tail;
  long       last,  crnt;

  /* First build up list (in reverse order) using indices and not pointers */

  last = crnt = -1;
  while (strncmp(GetLine(fin,TRUE),"{ADL",4) == 0)
    { 
      last =  Read_ADL_Struct(last,fin);
    }
  if (CurLine[0] != '.')
    MgenError("Expecting end of ADL list");
  GET_EOM;

  /* Traverse again, reversing list order and converting indices to ptrs. */

  //  fprintf(stderr,"* crnt = %lx last = %lx \n", crnt, last);

  tail = NULL;

  crnt = last;

  while (crnt >= 0)
    { cptr = (AuditLine *) (MemBuffer + crnt);
      crnt = (long) (cptr->next);
      cptr->next    = tail;
      cptr->name    = MemBuffer + ((long) (cptr->name));
      cptr->version = MemBuffer + ((long) (cptr->version));
      cptr->comment = MemBuffer + ((long) (cptr->comment));
      
      tail = cptr;
    }
  amesg.list = tail;

  return (&amesg);
}

static void *Read_Frag_Mesg(FILE *fin, int frag_class)
{ // frag_class
  //   FragMesg 0
  //   InternalFragMesg 1
  //   OFGMesg 2
  static FragMesg fmesg;
  char   ch;

  assert((frag_class == 0) || (frag_class == 1) || (frag_class == 2));
  
  if (frag_class == 0)
    { 
      GET_TYPE(ch,ACT_FORMAT,"action");
      fmesg.action = (ActionType) ch;
      GET_FIELD(fmesg.eaccession,EACC_FORMAT,"accession field")
    }
  else
    {
      GET_TYPE(ch,ACT1_FORMAT "[AD]","action");
      fmesg.action = (ActionType) ch;
      GET_PAIR(fmesg.eaccession,fmesg.iaccession,
               IACCS_FORMAT,"accession field pair");
    } 

  fmesg.sequence = NULL;
  fmesg.quality  = NULL;

  if (fmesg.action == AS_ADD) { 
    if (frag_class == 0){
      GET_TYPE(ch,TYP_FORMAT,"type");
      fmesg.type = (FragType) ch;
      // We want to succeed on all reads, and let the gatekeeper do its stuff
    }else{
      GET_TYPE(ch,TYP1_FORMAT "[RXELTFUCBW]","type");
      fmesg.type = (FragType) ch;
    }

    fmesg.eseq_id = 0;

    fmesg.source   = (char *) GetText("src:",fin,FALSE);

    GET_FIELD(fmesg.entry_time,ETM_FORMAT,"time field");

    if( frag_class != 2 ) {
      fmesg.sequence = (char *) GetText("seq:",fin,TRUE);
      fmesg.quality  = (char *) GetText("qlt:",fin,TRUE);
    }
    GET_PAIR(fmesg.clear_rng.bgn,fmesg.clear_rng.end,
             CLR_FORMAT,"clear range field");
    fmesg.source   = MemBuffer + ((long) (fmesg.source));
    if( frag_class != 2 ) {
      // Convert from an index to a pointer.
      fmesg.sequence = MemBuffer + ((long) (fmesg.sequence));
      fmesg.quality  = MemBuffer + ((long) (fmesg.quality));
    }
  }  //  action is AS_ADD
  GET_EOM;
  return ((void *) (&fmesg));
}

static void *Read_FRG_Mesg(FILE *fin)
{ return Read_Frag_Mesg(fin,0); } 

static void *Read_IFG_Mesg(FILE *fin)
{ return Read_Frag_Mesg(fin,1); }

static void *Read_OFG_Mesg(FILE *fin)
{ return Read_Frag_Mesg(fin,2); }


static void *Read_OVL_Mesg(FILE *fin)
{ static OverlapMesg omesg;
  int    idx;
  char ch;
  
  GET_FIELD(omesg.aifrag,"afr:" F_IID,"a-fragment field");
  GET_FIELD(omesg.bifrag,"bfr:" F_IID,"b-fragment field");
  GET_TYPE(ch,ORI1_FORMAT "[NAIO]","orientation");
  omesg.orientation = (OrientType) ch;
  GET_TYPE(ch,"olt:%1[DCSXdc]","overlap");
  omesg.overlap_type = (OverlapType) ch;
  GET_FIELD(omesg.ahg,"ahg:" F_COORD,"a-hang field");
  GET_FIELD(omesg.bhg,"bhg:" F_COORD,"b-hang field");
  GET_FIELD(omesg.quality,QUA_FORMAT,"quality field");
  GET_FIELD(omesg.min_offset,"mno:" F_COORD,"min-offset field");
  GET_FIELD(omesg.max_offset,"mxo:" F_COORD,"max-offset field");
  GET_FIELD(omesg.polymorph_ct,"pct:" F_S32,"poly-count field");

  if (strncmp(GetLine(fin,TRUE),"del:",4) != 0)
    MgenError("delta tag label");
  idx = MoreSpace(2*AS_READ_MAX_LEN,1);
  omesg.delta = (signed char *) (MemBuffer + idx);

  { int i, n;     /* Read a delta item (only one of its kind) */
    char *t, *u;

    i = 0;
    while ((t = GetLine(fin,TRUE))[0] != '.')
      while (1)
        { n = strtol(t,&u,10);
          if (u == t) break;
          t = u;
          if (! isspace((int)*t))
            MgenError("Delta is not a sequence of digits");
          omesg.delta[i++] = n;
        }
    omesg.delta[i] = 0;
  }

  GET_EOM;
  return ((void *) (&omesg));
}

static void *Read_LKG_Mesg(FILE *fin)
{ static LinkMesg lmesg;
  char ch; 
  GET_TYPE(ch,ACT_FORMAT,"action");
  lmesg.action = (ActionType) ch;
  GET_TYPE(ch,TYP_FORMAT,"link");
  lmesg.type = (LinkType) ch;
  GET_FIELD(lmesg.frag1,EFRAG1_FORMAT,"fragment 1 field");
  GET_FIELD(lmesg.frag2,EFRAG2_FORMAT,"fragment 2 field");
  if (lmesg.action == AS_ADD)
    {
      GET_FIELD(lmesg.entry_time,ETM_FORMAT,"entry time field");
      GET_FIELD(lmesg.distance,EDST_FORMAT,"distance field");
      GET_TYPE(ch,ORI1_FORMAT "[NAIOU]","link orientation");
      lmesg.link_orient = (OrientType) ch;
    }
  GET_EOM;
  return ((void *) (&lmesg));
}

static void *Read_ILK_Mesg(FILE *fin)
{ static InternalLinkMesg lmesg;
  char ch;
  
  GET_TYPE(ch,ACT1_FORMAT "[AD]","action");
  lmesg.action = (ActionType) ch;
  GET_TYPE(ch,TYP1_FORMAT "[MBSRYTG]","link");
  lmesg.type = (LinkType) ch;
  GET_FIELD(lmesg.ifrag1,IFRAG1_FORMAT,"fragment 1 field");
  GET_FIELD(lmesg.ifrag2,IFRAG2_FORMAT,"fragment 2 field");
  if (lmesg.action == AS_ADD)
    { GET_FIELD(lmesg.entry_time,ETM_FORMAT,"entry time field");
      GET_FIELD(lmesg.idistance,IDST_FORMAT,"distance field");
      GET_TYPE(ch,ORI1_FORMAT "[NAIOU]","link orientation");
      lmesg.link_orient = (OrientType) ch;
    }
  GET_EOM;
  return ((void *) (&lmesg));
}

static void *Read_UOM_Mesg(FILE *fin)
{ static UnitigOverlapMesg	mesg;
  char ch; 
  GET_FIELD(mesg.chunk1,"ck1:" F_IID,"chunk 1 id field");
  GET_FIELD(mesg.chunk2,"ck2:" F_IID,"chunk 2 id field");
  GET_TYPE(ch,ORI1_FORMAT "[NAIO]","orientation");
  mesg.orient = (ChunkOrientationType) ch;
  GET_TYPE(ch,"ovt:%1[NOTCIMXdcYZ]","overlap type");
  mesg.overlap_type = (UnitigOverlapType) ch;

  #ifdef AS_ENABLE_SOURCE
  {
    long		sindx;
    sindx = GetText("src:",fin,FALSE);
    mesg.source = MemBuffer + sindx;
  }
  #endif
  GET_FIELD(mesg.best_overlap_length,"len:" F_COORD,"best overlap");
  GET_FIELD(mesg.min_overlap_length,MINC_FORMAT,"min overlap");
  GET_FIELD(mesg.max_overlap_length,MAXC_FORMAT,"max overlap");
  GET_FIELD(mesg.quality,QUA_FORMAT,"quality field");
  GET_EOM;
  return ((void *) (&mesg));
}

static void *Read_FOM_Mesg(FILE *fin)
{ static FragOverlapMesg	mesg;
  char ch;
  
  GET_FIELD(mesg.afrag,"afr:" F_IID,"fragment A id field");
  GET_FIELD(mesg.bfrag,"bfr:" F_IID,"fragment B id field");
  GET_TYPE(ch,ORI1_FORMAT "[NAIO]","orientation");
  mesg.orient = (ChunkOrientationType) ch;
#if 0
  GET_TYPE(ch,"ovt:%1[NOCIXTBHW]","overlap type");
#else
  GET_TYPE(ch,"ovt:%1[NOTCIMXYZ]","overlap type");
#endif
  mesg.overlap_type = (UnitigOverlapType) ch;
  #ifdef AS_ENABLE_SOURCE
  {
    long		sindx;
    sindx = GetText("src:",fin,FALSE);
    mesg.source = MemBuffer + sindx;
  }
  #endif
  GET_FIELD(mesg.best_overlap_length,"len:" F_COORD,"best overlap");
  GET_FIELD(mesg.min_overlap_length,MINC_FORMAT,"min overlap");
  GET_FIELD(mesg.max_overlap_length,MAXC_FORMAT,"max overlap");
  GET_FIELD(mesg.quality,QUA_FORMAT,"quality field");
  GET_EOM;
  return ((void *) (&mesg));
}

static void Read_IMP_Mesg(FILE *fin, long indx)
{
  IntMultiPos	*imp;
  int		 i;
  int32		 n, *delta;
  long		 tindx;
  char		*line, *u;
  char ch;
  
  imp = (IntMultiPos *) (MemBuffer + indx);
  GET_TYPE(ch,TYP1_FORMAT "[RXTELUFSUcBCG]","multipos$");
  imp->type = (FragType) ch;
  GET_FIELD(imp->ident,    "mid:" F_IID,"multipos id");
  GET_FIELD(imp->contained,"con:" F_IID,"contained id");
  #ifdef NEW_UNITIGGER_INTERFACE
  GET_FIELD(imp->ident2,   "bid:" F_IID,"multipos id");
  #endif
#ifdef AS_ENABLE_SOURCE
  imp->sourceInt = -1;
#endif
  GET_PAIR(imp->position.bgn,imp->position.end,
           POS2_FORMAT,"position field");
  #ifdef NEW_UNITIGGER_INTERFACE
  GET_FIELD(imp->ahang       ,"ahg:" F_S32,"ahang");
  GET_FIELD(imp->bhang       ,"bhg:" F_S32,"bhang");
  #endif
  GET_FIELD(imp->delta_length,"dln:" F_S32,"delta length");
  if (strncmp(GetLine(fin,TRUE),"del:",4) != 0)
    MgenError("Missing del: field");
  if (imp->delta_length > 0) {
    tindx = MoreSpace(sizeof(int32)*imp->delta_length,8);
    imp = (IntMultiPos *) (MemBuffer + indx);	// in case of realloc
    imp->delta = (int32 *) tindx;
    delta = (int32 *) (MemBuffer + (long) imp->delta);
    i = 0;
    while (i < imp->delta_length) {
      line = GetLine(fin,TRUE);
      n = (int32) strtol(line,&u,10);
      while (u != line) {
	line = u;
	delta[i++] = n;
        n = (int32) strtol(line,&u,10);
      }
    }
  } 
  else
    imp->delta = NULL;
  GET_EOM;
  return;
}

static void Read_IMV_Mesg(FILE *fin, long indx)
{
  IntMultiVar *imv = (IntMultiVar *) (MemBuffer + indx);

  GET_PAIR(imv->position.bgn,imv->position.end, POS2_FORMAT,"position field");
  GET_FIELD(imv->num_reads,"nrd:" F_S32,"number of reads"); 
  GET_FIELD(imv->num_conf_alleles,"nca:" F_S32,"number of confirmed alleles");
  GET_FIELD(imv->anchor_size,"anc:" F_S32,"anchor size");
  GET_FIELD(imv->var_length,"len:" F_S32,"length field");
  imv->nr_conf_alleles = (char *)GetText("nra:",fin,FALSE);
  imv->weights         = (char *)GetText("wgt:",fin,FALSE);
  imv->var_seq =         (char *)GetText("seq:",fin,FALSE);
  GET_EOM;
  return;
}

static void Read_VAR_Mesg(FILE *fin, long indx)
{
  IntMultiVar *smv = (IntMultiVar *) (MemBuffer + indx);

  GET_PAIR(smv->position.bgn,smv->position.end, POS2_FORMAT,"position field");
  GET_FIELD(smv->num_reads,"nrd:" F_S32,"number of reads");
  GET_FIELD(smv->num_conf_alleles,"nca:" F_S32,"number of confirmed alleles");
  GET_FIELD(smv->anchor_size,"anc:" F_S32,"anchor size");
  GET_FIELD(smv->var_length,"len:" F_S32,"length field");
  smv->nr_conf_alleles = (char *)GetText("nra:",fin,FALSE);
  smv->weights         = (char *)GetText("wgt:",fin,FALSE);
  smv->var_seq         = (char *)GetText("seq:",fin,FALSE);
  GET_EOM;
  return;
}

static void Read_IUP_Mesg(FILE *fin, long indx)
{
  IntUnitigPos		*iup;
  int			i;
  int32			n, *delta;
  long			tindx;
  char			*line, *u;
  char ch;
  
  iup = (IntUnitigPos *) (MemBuffer + indx);
  GET_TYPE(ch,TYP1_FORMAT "[URSPsX]","unitigpos type");
  iup->type = (UnitigType) ch;
  GET_FIELD(iup->ident,"lid:" F_IID,"unitigpos id");
  GET_PAIR(iup->position.bgn,iup->position.end,
           POS2_FORMAT,"position field");
  GET_FIELD(iup->delta_length,"dln:" F_S32,"delta length");
  if (strncmp(GetLine(fin,TRUE),"del:",4) != 0)
    MgenError("Missing del: field");
  if (iup->delta_length > 0) {
    tindx = MoreSpace(sizeof(int32)*iup->delta_length,8);
    iup = (IntUnitigPos *) (MemBuffer + indx);	// in case of realloc
    iup->delta = (int32 *) tindx;
    delta = (int32 *) (MemBuffer + (long) iup->delta);
    i = 0;
    while (i < iup->delta_length) {
      line = GetLine(fin,TRUE);
      n = (int32) strtol(line,&u,10);
      while (u != line) {
	line = u;
	delta[i++] = n;
        n = (int32) strtol(line,&u,10);
      }
    }
  } 
  else
    iup->delta = NULL;
  GET_EOM;
  return;
}

static void *Read_IUM_Mesg(FILE *fin)
{ static IntUnitigMesg		mesg;
  int				i;
  long				cindx, qindx, mpindx, indx;
# ifdef AS_ENABLE_SOURCE
  long				sindx;
# endif
  char ch;

  GET_FIELD(mesg.iaccession,"acc:" F_IID,"accession field");
# ifdef AS_ENABLE_SOURCE
  sindx = GetText("src:",fin,FALSE);
# endif
  GET_FIELD(mesg.coverage_stat,"cov:%f","coverage stat");
  GET_TYPE(ch,"sta:%1[UCNSX]","status");
  mesg.status = (UnitigStatus) ch;
  GET_FIELD(mesg.a_branch_point,"abp:" F_COORD,"a branch point");
  GET_FIELD(mesg.b_branch_point,"bbp:" F_COORD,"b branch point");
  GET_FIELD(mesg.length,"len:" F_COORD,"length field");
  cindx = GetText("cns:",fin,TRUE);
  qindx = GetText("qlt:",fin,TRUE);
  GET_FIELD(mesg.forced,"for:" F_S32,"forced booleon");
  GET_FIELD(mesg.num_frags,"nfr:" F_S32,"num frags field");
  if (mesg.num_frags > 0) {
    indx = mpindx = MoreSpace(mesg.num_frags*sizeof(IntMultiPos),8);
    /* Why use 8 boundary above? Damned if I know - ela */
    for (i=0; i < mesg.num_frags; ++i) {
      if (strncmp(GetLine(fin,TRUE),"{IMP",4) != 0)
	MgenError("Expecting IMP record");
      Read_IMP_Mesg(fin, indx);
      indx += sizeof(IntMultiPos);
    }
    mesg.f_list = (IntMultiPos *) (MemBuffer + mpindx);
  }
  else
    mesg.f_list = NULL;
  GET_EOM;
  mesg.consensus = MemBuffer + cindx;
  mesg.quality  = MemBuffer + qindx;

  assert(strlen(mesg.consensus) == strlen(mesg.quality));
  assert((strlen(mesg.consensus) == mesg.length) ||
	 (strlen(mesg.consensus) == 0) );

  for (i=0; i < mesg.num_frags; ++i) {
    if (mesg.f_list[i].delta_length > 0)
      mesg.f_list[i].delta = (int32 *) (MemBuffer+(long) mesg.f_list[i].delta);
  }
# ifdef AS_ENABLE_SOURCE
  mesg.source = MemBuffer + sindx;
# endif
  return ((void *) (&mesg));
}

static void *Read_IUL_Mesg(FILE *fin)
{ static IntUnitigLinkMesg	mesg;
  int				i,size;
  long				indx;
  IntMate_Pairs			*imp;
  char ch;
  
  GET_FIELD(mesg.unitig1,"ut1:" F_IID,"unitig 1 field");
  GET_FIELD(mesg.unitig2,"ut2:" F_IID,"unitig 2 field");
  GET_TYPE(ch,ORI1_FORMAT "[NAOI]","orientation");
  mesg.orientation = (ChunkOrientationType) ch;
#if 0
  GET_TYPE(ch,"ovt:%1[NOCIXTBHW]","overlap type");
#else
  GET_TYPE(ch,"ovt:%1[NOTCIMXYZ]","overlap type");
#endif
  mesg.overlap_type = (UnitigOverlapType) ch;
  GET_FIELD(mesg.is_possible_chimera,"ipc:" F_S32,"warning");
  GET_FIELD(mesg.includes_guide,"gui:" F_S32,"guide flag");
  GET_FIELD(mesg.mean_distance,MEA_IN_FORMAT,"mean distance");
  GET_FIELD(mesg.std_deviation,STD_IN_FORMAT,"standard deviation");
  GET_FIELD(mesg.num_contributing,"num:" F_S32,"number of links");
  GET_TYPE(ch,"sta:%1[APBCU]","placement status");
  mesg.status = (PlacementStatusType) ch;
  if (strncmp(GetLine(fin,TRUE),"jls:",4) != 0)
      MgenError("Expecting jls field");
  size = mesg.num_contributing;
  if (mesg.overlap_type != AS_NO_OVERLAP)
    --size;
  if (size > 0) {
    indx = MoreSpace(sizeof(IntMate_Pairs)*size,8);
    imp = mesg.jump_list = (IntMate_Pairs *) (MemBuffer + indx);
    for (i=0; i < size; ++i) {
      GET_TRIPLE(imp->in1,imp->in2,ch, F_IID "," F_IID ",%1[MSBRYT]",
                 "mate pair");
      imp->type = (LinkType)ch;
      ++imp;
    }
  }
  else
    mesg.jump_list = NULL;
  GET_EOM;
  return ((void *) (&mesg));
}

static void *Read_ICL_Mesg(FILE *fin)
{ static IntContigLinkMesg	mesg;
  int				i,size;
  long				indx;
  IntMate_Pairs			*imp;
  char ch;
  
  GET_FIELD(mesg.contig1,"co1:" F_IID,"contig 1 field");
  GET_FIELD(mesg.contig2,"co2:" F_IID,"contig 2 field");
  GET_TYPE(ch,ORI1_FORMAT "[NAOI]","orientation");
  mesg.orientation = (ChunkOrientationType) ch;
#if 0
  GET_TYPE(ch,"ovt:%1[NOCIXTBHW]","overlap type");
  
#else
  GET_TYPE(ch,"ovt:%1[NOTCIMXYZ]","overlap type");
#endif
  mesg.overlap_type = (UnitigOverlapType) ch;
  GET_FIELD(mesg.is_possible_chimera,"ipc:" F_S32,"warning");
  GET_FIELD(mesg.includes_guide,"gui:" F_S32,"guide flag");
  GET_FIELD(mesg.mean_distance,MEA_IN_FORMAT,"mean distance");
  GET_FIELD(mesg.std_deviation,STD_IN_FORMAT,"standard deviation");
  GET_FIELD(mesg.num_contributing,"num:" F_S32,"number of links");
  GET_TYPE(ch,"sta:%1[APBCU]","placement status");
  mesg.status = (PlacementStatusType) ch;
  if (strncmp(GetLine(fin,TRUE),"jls:",4) != 0)
      MgenError("Expecting jls field");
  size = mesg.num_contributing;
  if (mesg.overlap_type != AS_NO_OVERLAP)
    --size;
  if (size > 0) {
    indx = MoreSpace(sizeof(IntMate_Pairs)*size,8);
    imp = mesg.jump_list = (IntMate_Pairs *) (MemBuffer + indx);
    for (i=0; i < size; ++i) {
      GET_TRIPLE(imp->in1,imp->in2,ch, F_IID "," F_IID ",%1[MSBRYT]",
                 "mate pair");
      imp->type = (LinkType)ch;
      ++imp;
    }
  }
  else
    mesg.jump_list = NULL;
  GET_EOM;
  return ((void *) (&mesg));
}

static void *Read_ISL_Mesg(FILE *fin)
{ static InternalScaffoldLinkMesg	mesg;
  int				i,size;
  long				indx;
  IntMate_Pairs			*imp;
  char ch;
  
  GET_FIELD(mesg.iscaffold1,"sc1:" F_IID,"scaffold 1 field");
  GET_FIELD(mesg.iscaffold2,"sc2:" F_IID,"scaffold 2 field");
  GET_TYPE(ch,ORI1_FORMAT "[NAOI]","orientation");
  mesg.orientation = (ChunkOrientationType) ch;
  GET_FIELD(mesg.includes_guide,"gui:" F_S32,"guide flag");
  GET_FIELD(mesg.mean_distance,MEA_IN_FORMAT,"mean distance");
  GET_FIELD(mesg.std_deviation,STD_IN_FORMAT,"standard deviation");
  GET_FIELD(mesg.num_contributing,"num:" F_S32," number of links");
  if (strncmp(GetLine(fin,TRUE),"jls:",4) != 0)
      MgenError("Expecting jls field");
  size = mesg.num_contributing;
  assert(size > 0);
  indx = MoreSpace(sizeof(IntMate_Pairs)*size,8);
  imp = mesg.jump_list = (IntMate_Pairs *) (MemBuffer + indx);
  for (i=0; i < size; ++i) {
    GET_TRIPLE(imp->in1,imp->in2,ch, F_IID "," F_IID ",%1[MSBRYT]",
               "mate pair");
    imp->type = (LinkType) ch;
    ++imp;
  }
  GET_EOM;
  return ((void *) (&mesg));
}

static void *Read_AFG_Mesg(FILE *fin)
{ static AugFragMesg		mesg;
  char *line;
  int32 i,j;
  char ch;
  
  GET_PAIR(mesg.eaccession,mesg.iaccession,IACCS_FORMAT,"accession field");
  GET_TYPE(ch,"mst:%1[ZGCLSONHADEURF]","mate status");
  mesg.mate_status = (MateStatType) ch;
  GET_FIELD(mesg.chimeric,"chi:" F_S32,"chimeric flag");
#ifdef AFG_BACKWARDS_COMPATIBLE
  line = GetLine(fin, TRUE);
  if(sscanf(line,"cha:" F_S32,&i)==1){
    mesg.chaff=i;
    line = GetLine(fin, TRUE);
  }
  else
    mesg.chaff=0;

  if(sscanf(line,CLR_FORMAT,&i,&j)==2){
    mesg.clear_rng.bgn=i;
    mesg.clear_rng.end=j;
  }
  else
    MfieldError("chaff flag/clear range");	
#else
  GET_FIELD(mesg.chaff,"cha:" F_S32,"chaff flag");
  GET_PAIR(mesg.clear_rng.bgn,mesg.clear_rng.end,
           CLR_FORMAT,
	   "clear range");
#endif
  GET_EOM;
  return ((void *) (&mesg));
}

static void Read_ICP_Mesg(FILE *fin, IntContigPairs *icp)
{
  char ch;
  GET_FIELD(icp->contig1,"ct1:" F_IID,"contig 1 id");
  GET_FIELD(icp->contig2,"ct2:" F_IID,"contig 2 id");
  GET_FIELD(icp->mean,MEA_IN_FORMAT,"mean distance");
  GET_FIELD(icp->stddev,STD_IN_FORMAT,"standard deviation");
  GET_TYPE(ch,ORI1_FORMAT "[NAIOU]","link orientation");
  icp->orient = (ChunkOrientationType) ch;
  GET_EOM;
  return;
}

static void *Read_ISF_Mesg(FILE *fin)
{ static IntScaffoldMesg	mesg;
  int				i, num;
  long				indx;
  IntContigPairs		*icp;

  GET_FIELD(mesg.iaccession,"acc:" F_IID,"ISF accession");
  GET_FIELD(mesg.num_contig_pairs,"noc:" F_S32,"number of contigs");
  num = MAX(1,mesg.num_contig_pairs);
  if (num > 0) {
    indx = MoreSpace(num*sizeof(IntContigPairs),8);
    icp = mesg.contig_pairs = (IntContigPairs *) (MemBuffer + indx);
    for (i=0; i < num; ++i) {
      if (strncmp(GetLine(fin,TRUE),"{ICP",4) != 0)
	MgenError("Expecting ICP record");
      Read_ICP_Mesg(fin,icp);
      ++icp;
    }
  }
  else
    mesg.contig_pairs = NULL;
  GET_EOM;
  return ((void *) (&mesg));
}

static void *Read_IMD_Mesg(FILE *fin)
{ static IntMateDistMesg	mesg;
  long				indx;
  int				i;

  GET_FIELD(mesg.refines,"ref:" F_IID,"distance id");
  GET_FIELD(mesg.mean,MEA_IN_FORMAT,"mean distance");
  GET_FIELD(mesg.stddev,STD_IN_FORMAT,"standard deviation");
  GET_FIELD(mesg.min,MINC_FORMAT,"min distance");
  GET_FIELD(mesg.max,MAXC_FORMAT,"max distance");
  GET_FIELD(mesg.num_buckets,"buc:" F_S32,"number of buckets");
  if (strncmp(GetLine(fin,TRUE),"his:",4) != 0)
    MgenError("Expecting his field");
  if (mesg.num_buckets > 0) {
    indx = MoreSpace(mesg.num_buckets*sizeof(int32),8);
    mesg.histogram = (int32 *) (MemBuffer + indx);
    for (i=0; i < mesg.num_buckets; ++i)
      GET_FIELD(mesg.histogram[i],F_S32,"histogram entry");
  }
  else
    mesg.histogram = NULL;
  GET_EOM;
  return ((void *) (&mesg));
}

static void *Read_ICM_Mesg(FILE *fin)
{ static IntConConMesg		mesg;
  long	 cindx, qindx, mpindx, upindx, indx, uindx, vindx, vpindx;
  int				i;
  char ch;
  
  GET_FIELD(mesg.iaccession,"acc:" F_IID,"accession number");
  GET_TYPE(ch,"pla:%1[PU]"," placed flag");
  mesg.placed = (ContigPlacementStatusType) ch;
  GET_FIELD(mesg.length,"len:" F_COORD," contig length");
  cindx = GetText("cns:",fin,TRUE);
  qindx = GetText("qlt:",fin,TRUE);
  GET_FIELD(mesg.forced,"for:" F_S32," forced flag");
  GET_FIELD(mesg.num_pieces,"npc:" F_S32," number of pieces");
  GET_FIELD(mesg.num_unitigs,"nou:" F_S32," number of unitigs");
  if (!novar)
  {
      GET_FIELD(mesg.num_vars,"nvr:" F_S32,"num vars field");
  }
  /* Why use 8 boundary above & below? Damned if I know - ela */
  vindx = vpindx = MoreSpace(mesg.num_vars   *sizeof(IntMultiVar),8); 
  indx  = mpindx = MoreSpace(mesg.num_pieces *sizeof(IntMultiPos),8);
  uindx = upindx = MoreSpace(mesg.num_unitigs*sizeof(IntUnitigPos),8);

  if (!novar)
  {
     if (mesg.num_vars > 0)
     {
       mesg.v_list = (IntMultiVar *) (MemBuffer + vpindx);
       for (i=0; i < mesg.num_vars; ++i) {
         if (strncmp(GetLine(fin,TRUE),"{IMV",4) != 0)
           MgenError("Expecting IMV record");
         Read_IMV_Mesg(fin, vindx);
         vindx += sizeof(IntMultiVar);
       }
     }
  else
    mesg.v_list = NULL;
  }

  if (mesg.num_pieces > 0)
  {
    mesg.pieces = (IntMultiPos *) (MemBuffer + mpindx);
    for (i=0; i < mesg.num_pieces; ++i) {
      if (strncmp(GetLine(fin,TRUE),"{IMP",4) != 0)
        MgenError("Expecting IMP record");
      Read_IMP_Mesg(fin, indx);
      indx += sizeof(IntMultiPos);
    }
  }  
  else
    mesg.pieces = NULL;

  if (mesg.num_unitigs > 0) {
    mesg.unitigs = (IntUnitigPos *) (MemBuffer + upindx);
    for (i=0; i < mesg.num_unitigs; ++i) {
      if (strncmp(GetLine(fin,TRUE),"{IUP",4) != 0)
	MgenError("Expecting IUP record");
      Read_IUP_Mesg(fin, uindx);
      uindx += sizeof(IntUnitigPos);
    }
  }
  else
    mesg.unitigs = NULL;
  GET_EOM;

  mesg.consensus = MemBuffer + cindx;
  mesg.quality = MemBuffer + qindx;
  if (!novar)
  {
     if (mesg.num_vars > 0)
       mesg.v_list = (IntMultiVar *) (MemBuffer + vpindx);
     else
       mesg.v_list = NULL;
     for (i=0; i < mesg.num_vars; ++i) {
       mesg.v_list[i].nr_conf_alleles = MemBuffer + 
           (long) mesg.v_list[i].nr_conf_alleles;
       mesg.v_list[i].weights = MemBuffer + (long) mesg.v_list[i].weights;
       mesg.v_list[i].var_seq = MemBuffer + (long) mesg.v_list[i].var_seq;
     }
  }

  if (mesg.num_pieces > 0)
    mesg.pieces = (IntMultiPos *) (MemBuffer + mpindx);
  else
    mesg.pieces = NULL;
  for (i=0; i < mesg.num_pieces; ++i) {
    if (mesg.pieces[i].delta_length > 0)
      mesg.pieces[i].delta = (int32 *) (MemBuffer +
				 (long) mesg.pieces[i].delta);
  }
  if (mesg.num_unitigs > 0)
    mesg.unitigs = (IntUnitigPos *) (MemBuffer + upindx);
  else
    mesg.unitigs = NULL;
  for (i=0; i < mesg.num_unitigs; ++i) {
    if (mesg.unitigs[i].delta_length > 0)
      mesg.unitigs[i].delta = (int32 *) (MemBuffer +
				 (long) mesg.unitigs[i].delta);
  }
  return ((void *) (&mesg));
}


static void *Read_IAF_Mesg(FILE *fin)
{ static IntAugFragMesg		mesg;
  char ch;
  
  GET_FIELD(mesg.iaccession,"acc:" F_IID,"accession field");
  GET_TYPE(ch,TYP1_FORMAT "[RXELTFSUCBWG]","type");
  mesg.type = (FragType) ch;
  GET_FIELD(mesg.chimeric,"chi:" F_S32,"chimeric flag");
#ifdef IAF_BACKWARDS_COMPATIBLE
	{
		int i, j;
  		char *line;

		line = GetLine(fin, TRUE);
		if(sscanf(line,"cha:" F_S32,&i)==1)
		{
			mesg.chaff=i;
                        line = GetLine(fin, TRUE);
		}
		else
			mesg.chaff=0;
		if(sscanf(line,CLR_FORMAT,&i,&j)==2)
		{
			mesg.clear_rng.bgn=i;
			mesg.clear_rng.end=j;
		}
		else
    		MgenError("IAF: choked on cha: or chr:");
	}
#else
  GET_FIELD(mesg.chaff,"cha:" F_S32,"chaff flag");
  GET_PAIR(mesg.clear_rng.bgn,mesg.clear_rng.end,CLR_FORMAT,
	   "clear range");
#endif
  GET_TYPE(ch,"mst:%1[ZGCLSONHADEURF]","mate status");
  mesg.mate_status = (MateStatType) ch;
  GET_EOM;
  return ((void *) (&mesg));
}


static void *Read_EOF_Mesg(FILE *fin)
{
  static EndOfFileMesg mesg;
  int commentidx;

  GET_FIELD(mesg.status, "sta:" F_S32, "status field");
  GET_FIELD(mesg.created, CRT_FORMAT, "entry time field");
  commentidx = GetText("com:",fin,FALSE);
  mesg.comment = MemBuffer + commentidx;
  GET_EOM;

  return ((void *) &mesg);
}



/* Genome snapshot input routines */
/**********************************/

static void Read_MPS_Mesg(FILE *fin, long indx)
{
  SnapMultiPos		*imp;
  int			i;
  int32			n, *delta;
  long			tindx;
  char			*line, *u;
  char ch;
  
  imp = (SnapMultiPos *) (MemBuffer + indx);
  GET_TYPE(ch,TYP1_FORMAT "[RXTEFUSLuBG]","multipos type");
  imp->type = (FragType) ch;
  GET_FIELD(imp->eident,"mid:" F_UID,"multipos id");
#ifdef AS_ENABLE_SOURCE
  tindx = GetText("src:",fin,FALSE);
  imp = (SnapMultiPos *) (MemBuffer + indx);	// in case of realloc
  imp->source = (char *) tindx;
#endif
  GET_PAIR(imp->position.bgn,imp->position.end,
           POS2_FORMAT,"position field");
  GET_FIELD(imp->delta_length,"dln:" F_S32,"delta length");
  if (strncmp(GetLine(fin,TRUE),"del:",4) != 0)
    MgenError("Missing del: field");
  if (imp->delta_length > 0) {
    tindx = MoreSpace(sizeof(int32)*imp->delta_length,8);
    imp = (SnapMultiPos *) (MemBuffer + indx);	// in case of realloc
    imp->delta = (int32 *) tindx;
    delta = (int32 *) (MemBuffer + (long) imp->delta);
    i = 0;
    while (i < imp->delta_length) {
      line = GetLine(fin,TRUE);
      n = (int32) strtol(line,&u,10);
      while (u != line) {
	line = u;
	delta[i++] = n;
        n = (int32) strtol(line,&u,10);
      }
    }
  } 
  else
    imp->delta = NULL;
  GET_EOM;
  return;
}

static void Read_UPS_Mesg(FILE *fin, long indx)
{
  UnitigPos		*iup;
  int			i;
  int32			n, *delta;
  long			tindx;
  char			*line, *u;
  char ch;
  
  iup = (UnitigPos *) (MemBuffer + indx);
  GET_TYPE(ch,TYP1_FORMAT "[URSPs]","unitigpos type");
  iup->type = (UnitigType) ch;
  GET_FIELD(iup->eident,"lid:" F_UID,"unitig id");
  GET_PAIR(iup->position.bgn,iup->position.end,
           POS2_FORMAT,"position field");
  GET_FIELD(iup->delta_length,"dln:" F_S32,"delta length");
  if (strncmp(GetLine(fin,TRUE),"del:",4) != 0)
    MgenError("Missing del: field");
  if (iup->delta_length > 0) {
    tindx = MoreSpace(sizeof(int32)*iup->delta_length,8);
    iup = (UnitigPos *) (MemBuffer + indx);	// in case of realloc
    iup->delta = (int32 *) tindx;
    delta = (int32 *) (MemBuffer + (long) iup->delta);
    i = 0;
    while (i < iup->delta_length) {
      line = GetLine(fin,TRUE);
      n = (int32) strtol(line,&u,10);
      while (u != line) {
	line = u;
	delta[i++] = n;
        n = (int32) strtol(line,&u,10);
      }
    }
  } 
  else
    iup->delta = NULL;
  GET_EOM;
  return;
}


static void Read_CTP_Mesg(FILE *fin, SnapContigPairs *icp)
{
  char ch;
  GET_FIELD(icp->econtig1,"ct1:" F_UID ,"contig 1 id");
  GET_FIELD(icp->econtig2,"ct2:" F_UID ,"contig 2 id");
  GET_FIELD(icp->mean,MEA_IN_FORMAT,"mean distance");
  GET_FIELD(icp->stddev,STD_IN_FORMAT,"standard deviation");
  GET_TYPE(ch,ORI1_FORMAT "[NAIOU]","link orientation");
  icp->orient = (ChunkOrientationType) ch;
  GET_EOM;
  return;
}


static void *Read_UTG_Mesg(FILE *fin)
{ static SnapUnitigMesg		mesg;
  int				i;
  long				cindx, qindx, mpindx, indx;
  #ifdef AS_ENABLE_SOURCE
  long				sindx;
  #endif
  char ch;

  GET_PAIR(mesg.eaccession,mesg.iaccession,IACCS_FORMAT,"accession field");
  #ifdef AS_ENABLE_SOURCE
  sindx = GetText("src:",fin,FALSE);
  #endif
  GET_FIELD(mesg.coverage_stat,"cov:%f","coverage stat");
  GET_TYPE(ch,"sta:%1[UCNSX]","status");
  mesg.status = (UnitigStatus) ch;
  GET_FIELD(mesg.a_branch_point,"abp:" F_COORD,"a branch point");
  GET_FIELD(mesg.b_branch_point,"bbp:" F_COORD,"b branch point");
  GET_FIELD(mesg.length,"len:" F_COORD,"length field");
  cindx = GetText("cns:",fin,TRUE);
  qindx = GetText("qlt:",fin,TRUE);
  GET_FIELD(mesg.forced,"for:" F_S32,"forced booleon");
  GET_FIELD(mesg.num_frags,"nfr:" F_S32,"num frags field");
  if (mesg.num_frags > 0) {
    indx = mpindx = MoreSpace(mesg.num_frags*sizeof(SnapMultiPos),8);
    /* Why use 8 boundary above? Damned if I know - ela */
    for (i=0; i < mesg.num_frags; ++i) {
      if (strncmp(GetLine(fin,TRUE),"{MPS",4) != 0)
	MgenError("Expecting MPS record");
      Read_MPS_Mesg(fin, indx);
      indx += sizeof(SnapMultiPos);
    }
    mesg.f_list = (SnapMultiPos *) (MemBuffer + mpindx);
  }
  else
    mesg.f_list = NULL;
  GET_EOM;
  mesg.consensus = MemBuffer + cindx;
  mesg.quality  = MemBuffer + qindx;
  for (i=0; i < mesg.num_frags; ++i) {
    #ifdef AS_ENABLE_SOURCE
    mesg.f_list[i].source = MemBuffer + (long) mesg.f_list[i].source;
    #endif
    if (mesg.f_list[i].delta_length > 0)
      mesg.f_list[i].delta = (int32 *) (MemBuffer+(long) mesg.f_list[i].delta);
  }
  #ifdef AS_ENABLE_SOURCE
  mesg.source = MemBuffer + sindx;
  #endif
  return ((void *) (&mesg));
}


static void *Read_ULK_Mesg(FILE *fin)
{ static SnapUnitigLinkMesg	mesg;
  int				i,size;
  long				indx;
  SnapMate_Pairs	       	*imp;
  char ch;

  GET_FIELD(mesg.eunitig1,"ut1:" F_UID,"unitig 1 field");
  GET_FIELD(mesg.eunitig2,"ut2:" F_UID,"unitig 2 field");
  GET_TYPE(ch,ORI1_FORMAT "[NAOI]","orientation");
  mesg.orientation = (ChunkOrientationType) ch;
#if 0
  GET_TYPE(ch,"ovt:%1[NOCIXTBHW]","overlap type");
#else
  GET_TYPE(ch,"ovt:%1[NOTCIMXYZ]","overlap type");
#endif
  mesg.overlap_type = (UnitigOverlapType) ch;
  GET_FIELD(mesg.is_possible_chimera,"ipc:" F_S32,"warning");
  GET_FIELD(mesg.includes_guide,"gui:" F_S32,"guide flag");
  GET_FIELD(mesg.mean_distance,MEA_IN_FORMAT,"mean distance");
  GET_FIELD(mesg.std_deviation,STD_IN_FORMAT,"standard deviation");
  GET_FIELD(mesg.num_contributing,"num:" F_S32,"number of links");
  GET_TYPE(ch,"sta:%1[APBCU]","placement status");
  mesg.status = (PlacementStatusType) ch;
  if (strncmp(GetLine(fin,TRUE),"jls:",4) != 0)
      MgenError("Expecting jls field");
  size = mesg.num_contributing;
  if (mesg.overlap_type != AS_NO_OVERLAP)
    --size;
  if (size > 0) {
    indx = MoreSpace(sizeof(SnapMate_Pairs)*size,8);
    imp = mesg.jump_list = (SnapMate_Pairs *) (MemBuffer + indx);
    for (i=0; i < size; ++i) {
      GET_TRIPLE(imp->in1,imp->in2,ch,
		 F_UID ","
		 F_UID ",%1[MSBRYT]","mate pair");
      imp->type = (LinkType)ch;
      ++imp;
    }
  }
  else
    mesg.jump_list = NULL;
  GET_EOM;
  return ((void *) (&mesg));
}


static void *Read_CCO_Mesg(FILE *fin)
{ static SnapConConMesg		mesg;
  long	 cindx, qindx, mpindx, upindx, indx, uindx, vindx, vpindx;
  int  	 i;
  char   ch;
  
  GET_PAIR(mesg.eaccession,mesg.iaccession,IACCS_FORMAT,"accession number");
  GET_TYPE(ch,"pla:%1[PU]"," placed flag");
  mesg.placed = (ContigPlacementStatusType) ch;
  GET_FIELD(mesg.length,"len:" F_COORD,"contig length");
  cindx = GetText("cns:",fin,TRUE);
  qindx = GetText("qlt:",fin,TRUE);
  GET_FIELD(mesg.forced,"for:" F_S32,"forced flag");
  GET_FIELD(mesg.num_pieces,"npc:" F_S32,"number of pieces");
  GET_FIELD(mesg.num_unitigs,"nou:" F_S32,"number of unitigs");
  if (!novar)
  {
     GET_FIELD(mesg.num_vars,"nvr:" F_S32,"number of vars");
  }
  /* Why use 8 boundary above & below? Damned if I know - ela */
  vindx = vpindx = MoreSpace(mesg.num_vars  *sizeof(IntMultiVar),8);
  indx  = mpindx = MoreSpace(mesg.num_pieces*sizeof(SnapMultiPos),8);
  uindx =upindx = MoreSpace(mesg.num_unitigs*sizeof(UnitigPos),8);

  if (!novar)
  {
     if (mesg.num_vars > 0)
     {
        mesg.vars = (IntMultiVar *) (MemBuffer + vpindx);
        for (i=0; i < mesg.num_vars; ++i) {
          if (strncmp(GetLine(fin,TRUE),"{VAR",4) != 0)
            MgenError("Expecting VAR record");
          Read_VAR_Mesg(fin, vindx);
          vindx += sizeof(IntMultiVar);
        }
      }
      else
        mesg.vars = NULL;
  }
// **************************************************
  for (i=0; i < mesg.num_pieces; ++i) {
    if (strncmp(GetLine(fin,TRUE),"{MPS",4) != 0)
      MgenError("Expecting MPS record");
    Read_MPS_Mesg(fin, indx);
    indx += sizeof(SnapMultiPos);
  }
// **************************************************
  if (mesg.num_unitigs > 0) {
    mesg.unitigs  = (UnitigPos *) (MemBuffer + upindx);
    for (i=0; i < mesg.num_unitigs; ++i) {
      if (strncmp(GetLine(fin,TRUE),"{UPS",4) != 0)
	MgenError("Expecting UPS record");
      Read_UPS_Mesg(fin, uindx);
      uindx += sizeof(UnitigPos);
    }
  }
  else
    mesg.unitigs = NULL;
  GET_EOM;
  mesg.consensus = MemBuffer + cindx;
  mesg.quality = MemBuffer + qindx;
  if (!novar)
  {
    if (mesg.num_vars > 0)
      mesg.vars = (IntMultiVar *) (MemBuffer + vpindx);
    else
      mesg.vars = NULL;
    for (i=0; i < mesg.num_vars; ++i) {
      mesg.vars[i].nr_conf_alleles = MemBuffer + 
          (long) mesg.vars[i].nr_conf_alleles;
      mesg.vars[i].weights = MemBuffer + (long) mesg.vars[i].weights;
      mesg.vars[i].var_seq = MemBuffer + (long) mesg.vars[i].var_seq;
    } 
  }
  if (mesg.num_pieces > 0)
    mesg.pieces = (SnapMultiPos *) (MemBuffer + mpindx);
  else
    mesg.pieces = NULL;
  for (i=0; i < mesg.num_pieces; ++i) {
    #ifdef AS_ENABLE_SOURCE
    mesg.pieces[i].source = MemBuffer + (long) mesg.pieces[i].source;
    #endif
    if (mesg.pieces[i].delta_length > 0)
      mesg.pieces[i].delta = (int32 *) (MemBuffer +
				 (long) mesg.pieces[i].delta);
  }
  if (mesg.num_unitigs > 0)
    mesg.unitigs = (UnitigPos *) (MemBuffer + upindx);
  else
    mesg.unitigs = NULL;
  for (i=0; i < mesg.num_unitigs; ++i) {
    if (mesg.unitigs[i].delta_length > 0)
      mesg.unitigs[i].delta = (int32 *) (MemBuffer +
				 (long) mesg.unitigs[i].delta);
  }

  return ((void *) (&mesg));
}



static void *Read_CLK_Mesg(FILE *fin)
{ static SnapContigLinkMesg	mesg;
  int				i,size;
  long				indx;
  SnapMate_Pairs		*imp;
  char ch;
  // char*                         line;
  
  GET_FIELD(mesg.econtig1,"c%*1[ot]1:" F_UID,"unitig 1 field");
  GET_FIELD(mesg.econtig2,"c%*1[ot]2:" F_UID,"unitig 2 field");

  GET_TYPE(ch,ORI1_FORMAT "[NAOI]","orientation");
  mesg.orientation = (ChunkOrientationType) ch;
#if 0
  GET_TYPE(ch,"ovt:%1[NOCIXTBHW]","overlap type");
#else
  GET_TYPE(ch,"ovt:%1[NOTCIMXYZ]","overlap type");
#endif
  mesg.overlap_type = (UnitigOverlapType) ch;
  GET_FIELD(mesg.is_possible_chimera,"ipc:" F_S32,"warning");
  GET_FIELD(mesg.includes_guide,"gui:" F_S32,"guide flag");
  GET_FIELD(mesg.mean_distance,MEA_IN_FORMAT,"mean distance");
  GET_FIELD(mesg.std_deviation,STD_IN_FORMAT,"standard deviation");
  GET_FIELD(mesg.num_contributing,"num:" F_S32,"number of links");
  GET_TYPE(ch,"sta:%1[APBCU]","placement status");
  mesg.status = (PlacementStatusType) ch;
  if (strncmp(GetLine(fin,TRUE),"jls:",4) != 0)
      MgenError("Expecting jls field");
  size = mesg.num_contributing;
  if (mesg.overlap_type != AS_NO_OVERLAP)
    --size;
  if (size > 0) {
    indx = MoreSpace(sizeof(SnapMate_Pairs)*size,8);
    imp = mesg.jump_list = (SnapMate_Pairs *) (MemBuffer + indx);
    for (i=0; i < size; ++i) {
      GET_TRIPLE(imp->in1,imp->in2,ch,
		 F_UID ","
		 F_UID ",%1[MSBRYT]","mate pair");
      imp->type = (LinkType)ch;
      ++imp;
    }
  }
  else
    mesg.jump_list = NULL;
  GET_EOM;
  return ((void *) (&mesg));
}

static void *Read_SLK_Mesg(FILE *fin)
{ static SnapScaffoldLinkMesg	mesg;
  int				i,size;
  long				indx;
  SnapMate_Pairs		*imp;
  char ch;
  
  GET_FIELD(mesg.escaffold1,"sc1:" F_UID,"scaffold 1 field");
  GET_FIELD(mesg.escaffold2,"sc2:" F_UID,"scaffold 2 field");

  GET_TYPE(ch,ORI1_FORMAT "[NAOI]","orientation");
  mesg.orientation = (ChunkOrientationType) ch;
  GET_FIELD(mesg.includes_guide,"gui:" F_S32,"guide flag");
  GET_FIELD(mesg.mean_distance,MEA_IN_FORMAT,"mean distance");
  GET_FIELD(mesg.std_deviation,STD_IN_FORMAT,"standard deviation");
  GET_FIELD(mesg.num_contributing,"num:" F_S32,"number of links");
  if (strncmp(GetLine(fin,TRUE),"jls:",4) != 0)
      MgenError("Expecting jls field");
  size = mesg.num_contributing;
  assert(size > 0) ;
  indx = MoreSpace(sizeof(SnapMate_Pairs)*size,8);
  imp = mesg.jump_list = (SnapMate_Pairs *) (MemBuffer + indx);
  for (i=0; i < size; ++i) {
    GET_TRIPLE(imp->in1,imp->in2,ch,
	       F_UID ","
	       F_UID ",%1[MSBRYT]","mate pair");
    imp->type = (LinkType)ch;
    ++imp;
  }
  GET_EOM;
  return ((void *) (&mesg));
}


static void *Read_SCF_Mesg(FILE *fin)
{ static SnapScaffoldMesg	mesg;
  int				i, num;
  long				indx;
  SnapContigPairs		*icp;

  GET_PAIR(mesg.eaccession,mesg.iaccession,IACCS_FORMAT,"accession number");
  GET_FIELD(mesg.num_contig_pairs,"noc:" F_S32,"number of contigs");
  num = MAX(mesg.num_contig_pairs, 1);
  if (num > 0) {
    indx = MoreSpace(num*sizeof(SnapContigPairs),8);
    icp = mesg.contig_pairs = (SnapContigPairs *) (MemBuffer + indx);
    for (i=0; i < num; ++i) {
      if (strncmp(GetLine(fin,TRUE),"{CTP",4) != 0)
	MgenError("Expecting CTP record");
      Read_CTP_Mesg(fin,icp);
      ++icp;
    }
  }
  else
    mesg.contig_pairs = NULL;
  GET_EOM;
  return ((void *) (&mesg));
}

static void *Read_DSC_Mesg(FILE *fin)
{ static SnapDegenerateScaffoldMesg	mesg;

  GET_FIELD(mesg.eaccession,EACC_FORMAT, "scaffold ID");
  GET_FIELD(mesg.econtig,"ctg:" F_UID, "contig ID");
  GET_EOM;
  return ((void *) (&mesg));
}

static void *Read_IDS_Mesg(FILE *fin)
{ static IntDegenerateScaffoldMesg	mesg;

  GET_FIELD(mesg.icontig,"ctg:" F_IID,"contig ID");
  GET_EOM;
  return ((void *) (&mesg));
}



static void *Read_MDI_Mesg(FILE *fin)
{ static SnapMateDistMesg	mesg;
  long				indx;
  int				i;

  GET_PAIR(mesg.erefines,mesg.irefines,IREFS_FORMAT,"distance id");
  GET_FIELD(mesg.mean,MEA_IN_FORMAT,"mean distance");
  GET_FIELD(mesg.stddev,STD_IN_FORMAT,"standard deviation");
  GET_FIELD(mesg.min,MINC_FORMAT,"min distance");
  GET_FIELD(mesg.max,MAXC_FORMAT,"max distance");
  GET_FIELD(mesg.num_buckets,"buc:" F_S32,"number of buckets");
  if (strncmp(GetLine(fin,TRUE),"his:",4) != 0)
    MgenError("Expecting his field");
  if (mesg.num_buckets > 0) {
    indx = MoreSpace(mesg.num_buckets*sizeof(int32),8);
    mesg.histogram = (int32 *) (MemBuffer + indx);
    for (i=0; i < mesg.num_buckets; ++i)
      GET_FIELD(mesg.histogram[i],F_S32,"histogram entry");
  }
  else
    mesg.histogram = NULL;
  GET_EOM;
  return ((void *) (&mesg));
}


static void *Read_BAT_Mesg(FILE *fin){
  static BatchMesg mesg;
  int nameidx, commentidx;
  
  nameidx = GetString("bna:",fin);
  GET_FIELD(mesg.created,CRT_FORMAT,"entry time field");
  GET_FIELD(mesg.eaccession,EACC_FORMAT,"accession number");
  commentidx = GetText("com:",fin, FALSE);

  mesg.comment = MemBuffer + commentidx;
  mesg.name = MemBuffer + nameidx;

  GET_EOM;
  return ((void *) (&mesg));
}


static void *Read_IBA_Mesg(FILE *fin){
  static BatchMesg mesg;
  int nameidx, commentidx;
  
  nameidx = GetString("bna:",fin);
  GET_FIELD(mesg.created,CRT_FORMAT,"entry time field");
  GET_PAIR(mesg.eaccession, mesg.iaccession,IACCS_FORMAT, "accession field pair");
  commentidx = GetText("com:",fin, FALSE);

  mesg.comment = MemBuffer + commentidx;
  mesg.name = MemBuffer + nameidx;

  GET_EOM;
  return ((void *) (&mesg));
}

/******************** OUTPUT ROUTINES ***************************/

/*  Routine to output each type of proto-IO message. */


static void Write_Dist_Mesg(FILE *fout, void *vmesg, int external)
{ InternalDistMesg *mesg = (InternalDistMesg *) vmesg;

  fprintf(fout,"{%s\n",(external?"DST":"IDT"));
  fprintf(fout,ACT_FORMAT "\n",mesg->action);
  if (external)
    fprintf(fout,EACC_FORMAT "\n",mesg->eaccession);
  else
    fprintf(fout,IACCS_FORMAT "\n",mesg->eaccession,mesg->iaccession);
  if (mesg->action != AS_DELETE)
    { fprintf(fout,MEA_OUT_FORMAT "\n",mesg->mean);
      fprintf(fout,STD_OUT_FORMAT "\n",mesg->stddev);
    }
  fprintf(fout,"}\n");
}

static void Write_DST_Mesg(FILE *fout, void *vmesg)
{ Write_Dist_Mesg(fout,vmesg,1); }

static void Write_IDT_Mesg(FILE *fout, void *vmesg)
{ Write_Dist_Mesg(fout,vmesg,0); }

static void Write_ADL_Struct(FILE *fout, AuditLine *mesg)
{ fprintf(fout,"{ADL\n");
  fprintf(fout,"who:%s\n",mesg->name);
  fprintf(fout,CTM_FORMAT "\n",mesg->complete);
  fprintf(fout,"vsn:%s\n",mesg->version);
  fflush(fout);
  PutText(fout,"com:",mesg->comment,FALSE);
  //  fprintf(fout,"com:%s\n",mesg->comment);
  fflush(fout);
  fprintf(fout,"}\n");
  fflush(fout);
}

static void Write_ADT_Mesg(FILE *fout, void *vmesg)
{ AuditMesg *mesg = (AuditMesg *) vmesg;
  AuditLine *a;

  fprintf(fout,"{ADT\n");
  for (a = mesg->list; a != NULL; a = a->next)
    Write_ADL_Struct(fout,a);
  fprintf(fout,".\n");
  fprintf(fout,"}\n");
}

static void Write_LKG_Mesg(FILE *fout, void *vmesg)
{ LinkMesg *mesg = (LinkMesg *) vmesg;

  fprintf(fout,"{LKG\n");
  fprintf(fout,ACT_FORMAT "\n",mesg->action);
  fprintf(fout,TYP_FORMAT "\n",(char) mesg->type);
  fprintf(fout,EFRAG1_FORMAT "\n",mesg->frag1);
  fprintf(fout,EFRAG2_FORMAT "\n",mesg->frag2);
  if(mesg->action == AS_ADD)
    { fprintf(fout,ETM_FORMAT "\n",mesg->entry_time);
      fprintf(fout,EDST_FORMAT "\n",mesg->distance);
      fprintf(fout,ORI_FORMAT "\n",mesg->link_orient);
    }
  fprintf(fout,"}\n");
}

static void Write_ILK_Mesg(FILE *fout, void *vmesg)
{ InternalLinkMesg *mesg = (InternalLinkMesg *) vmesg;

  fprintf(fout,"{ILK\n");
  fprintf(fout,ACT_FORMAT "\n",mesg->action);
  fprintf(fout,TYP_FORMAT "\n",(char) mesg->type);
  fprintf(fout,IFRAG1_FORMAT "\n",mesg->ifrag1);
  fprintf(fout,IFRAG2_FORMAT "\n",mesg->ifrag2);
  if(mesg->action == AS_ADD)
    { fprintf(fout,ETM_FORMAT "\n",mesg->entry_time);
      fprintf(fout,IDST_FORMAT "\n",mesg->idistance);
      fprintf(fout,ORI_FORMAT "\n",mesg->link_orient);
    }
  fprintf(fout,"}\n");
}

static void Write_Frag_Mesg(FILE *fout, void *vmesg, int frag_class) {
  FragMesg *mesg = (FragMesg *) vmesg;
  static const char * const header[]  = { "FRG", "IFG", "OFG" };

  assert((frag_class == 0) || (frag_class == 1) || (frag_class == 2));

  fprintf(fout,"{%s\n",header[frag_class]);
  fprintf(fout,ACT_FORMAT "\n",mesg->action);
  if (frag_class == 0)
    fprintf(fout,EACC_FORMAT "\n",mesg->eaccession);
  else
    fprintf(fout,IACCS_FORMAT "\n",mesg->eaccession,mesg->iaccession);

  if (mesg->action == AS_ADD) {
    fprintf(fout,TYP_FORMAT "\n",(char) mesg->type);
    PutText(fout,"src:",mesg->source,FALSE);
    fprintf(fout,ETM_FORMAT "\n",mesg->entry_time);
    if( frag_class != 2 ) {
      PutText(fout,"seq:",mesg->sequence,TRUE);
      PutText(fout,"qlt:",mesg->quality,TRUE);
    }
    fprintf(fout,CLR_FORMAT "\n",
            mesg->clear_rng.bgn,mesg->clear_rng.end);
  }

  fprintf(fout,"}\n");
}

static void Write_FRG_Mesg(FILE *fout, void *vmesg)
{ Write_Frag_Mesg(fout,vmesg,0); }

static void Write_IFG_Mesg(FILE *fout, void *vmesg)
{ Write_Frag_Mesg(fout,vmesg,1); }

static void Write_OFG_Mesg(FILE *fout, void *vmesg)
{ Write_Frag_Mesg(fout,vmesg,2); }


static void Write_OVL_Mesg(FILE *fout, void *vmesg)
{ OverlapMesg *omesg = (OverlapMesg *) vmesg;
  int i;

  fprintf(fout,"{OVL\n");
  fprintf(fout,"afr:" F_IID "\n",omesg->aifrag);
  fprintf(fout,"bfr:" F_IID "\n",omesg->bifrag);
  fprintf(fout,ORI_FORMAT "\n",omesg->orientation);
  fprintf(fout,"olt:%c\n",omesg->overlap_type);
  fprintf(fout,"ahg:" F_COORD "\n",omesg->ahg);
  fprintf(fout,"bhg:" F_COORD "\n",omesg->bhg);
  fprintf(fout,"qua:%.6f\n",omesg->quality);
  fprintf(fout,"mno:" F_COORD "\n",omesg->min_offset);
  fprintf(fout,"mxo:" F_COORD "\n",omesg->max_offset);
  fprintf(fout,"pct:" F_S32 "\n",omesg->polymorph_ct);
  fprintf(fout,"del:");
  for (i = 0; omesg->delta[i] != AS_ENDOF_DELTA_CODE; i++)
    { if ((i%15) == 0) fprintf(fout,"\n");
      fprintf(fout,"%4d ",omesg->delta[i]);
    }
  fprintf(fout,"\n");
  fprintf(fout,".\n");
  fprintf(fout,"}\n");
}

static void Write_UOM_Mesg(FILE *fout, void *vmesg)
{ UnitigOverlapMesg *mesg = (UnitigOverlapMesg *) vmesg;

  fprintf(fout,"{UOM\n");
  fprintf(fout,"ck1:" F_IID "\n",mesg->chunk1);
  fprintf(fout,"ck2:" F_IID "\n",mesg->chunk2);
  fprintf(fout,ORI_FORMAT "\n",mesg->orient);
  fprintf(fout,"ovt:%c\n",mesg->overlap_type);
  #ifdef AS_ENABLE_SOURCE
  PutText(fout,"src:",mesg->source,FALSE);
  #endif
  fprintf(fout,"len:" F_COORD "\n",mesg->best_overlap_length);
  fprintf(fout,MINC_FORMAT "\n",mesg->min_overlap_length);
  fprintf(fout,MAXC_FORMAT "\n",mesg->max_overlap_length);
  fprintf(fout,"qua:%.6f\n",mesg->quality);
  fprintf(fout,"}\n");
  return;
}

static void Write_FOM_Mesg(FILE *fout, void *vmesg)
{ FragOverlapMesg *mesg = (FragOverlapMesg *) vmesg;

  fprintf(fout,"{FOM\n");
  fprintf(fout,"afr:" F_IID "\n",mesg->afrag);
  fprintf(fout,"bfr:" F_IID "\n",mesg->bfrag);
  fprintf(fout,ORI_FORMAT "\n",mesg->orient);
  fprintf(fout,"ovt:%c\n",mesg->overlap_type);
  #ifdef AS_ENABLE_SOURCE
  PutText(fout,"src:",mesg->source,FALSE);
  #endif
  fprintf(fout,"len:" F_COORD "\n",mesg->best_overlap_length);
  fprintf(fout,MINC_FORMAT "\n",mesg->min_overlap_length);
  fprintf(fout,MAXC_FORMAT "\n",mesg->max_overlap_length);
  fprintf(fout,"qua:%.6f\n",mesg->quality);
  fprintf(fout,"}\n");
  return;
}

static void Write_IMP_Mesg(FILE *fout, IntMultiPos *mlp)
{ int i;

  fprintf(fout,"{IMP\n");
  fprintf(fout,TYP_FORMAT "\n",(char) mlp->type);
  fprintf(fout,"mid:" F_IID "\n",mlp->ident);
  fprintf(fout,"con:" F_IID "\n",mlp->contained);
  #ifdef NEW_UNITIGGER_INTERFACE
  fprintf(fout,"bid:" F_IID "\n",mlp->ident2);
  #endif
  fprintf(fout,POS2_FORMAT "\n",
          mlp->position.bgn,mlp->position.end);
  #ifdef NEW_UNITIGGER_INTERFACE
  fprintf(fout,"ahg:" F_S32 "\n",mlp->ahang);
  fprintf(fout,"bhg:" F_S32 "\n",mlp->bhang); 
  #endif
  fprintf(fout,"dln:" F_S32 "\n",mlp->delta_length);
  fprintf(fout,"del:\n");
  if (mlp->delta_length > 0 ) {
    for(i=0; i < mlp->delta_length; i++) {
      fprintf(fout,F_S32,mlp->delta[i]);
      if (i%20 == 19) fprintf(fout,"\n");
      else  fprintf(fout," ");
    }
    if (mlp->delta_length%20 != 0) fprintf(fout,"\n");
  }
  fprintf(fout,"}\n");
  return;
}

static void Write_IMV_Mesg(FILE *fout, IntMultiVar *imv)
{
  fprintf(fout,"{IMV\n");
  fprintf(fout, POS2_FORMAT "\n", imv->position.bgn,imv->position.end);
  fprintf(fout,"nrd:" F_S32 "\n",imv->num_reads);
  fprintf(fout,"nca:" F_S32 "\n",imv->num_conf_alleles);
  fprintf(fout,"anc:" F_S32 "\n",imv->anchor_size);
  fprintf(fout,"len:" F_S32 "\n",imv->var_length);
  PutText(fout,"nra:",imv->nr_conf_alleles,FALSE);
  PutText(fout,"wgt:",imv->weights,FALSE);
  PutText(fout,"seq:",imv->var_seq,FALSE);
  fprintf(fout,"}\n");
  return;
}

static void Write_VAR_Mesg(FILE *fout, IntMultiVar *smv)
{
  fprintf(fout,"{VAR\n");
  fprintf(fout, POS2_FORMAT "\n", smv->position.bgn,smv->position.end);
  fprintf(fout,"nrd:" F_S32 "\n",smv->num_reads);
  fprintf(fout,"nca:" F_S32 "\n",smv->num_conf_alleles);
  fprintf(fout,"anc:" F_S32 "\n",smv->anchor_size);
  fprintf(fout,"len:" F_S32 "\n",smv->var_length);
  PutText(fout,"nra:",smv->nr_conf_alleles,FALSE);
  PutText(fout,"wgt:",smv->weights,FALSE);
  PutText(fout,"seq:",smv->var_seq,FALSE);
  fprintf(fout,"}\n");
  return;
}

static void Write_IUP_Mesg(FILE *fout, IntUnitigPos *up)
{ int i;

  fprintf(fout,"{IUP\n");
  fprintf(fout,TYP_FORMAT "\n",(char) up->type);
  fprintf(fout,"lid:" F_IID "\n",up->ident);
  fprintf(fout,POS2_FORMAT "\n",up->position.bgn,up->position.end);
  fprintf(fout,"dln:" F_S32 "\n",up->delta_length);
  fprintf(fout,"del:\n");
  if (up->delta_length > 0 ) {
    for(i=0; i < up->delta_length; i++) {
      fprintf(fout,F_S32,up->delta[i]);
      if (i%20 == 19) fprintf(fout,"\n");
      else  fprintf(fout," ");
    }
    if (up->delta_length%20 != 0) fprintf(fout,"\n");
  }
  fprintf(fout,"}\n");
  return;
}

static void Write_IUM_Mesg(FILE *fout, void *vmesg)
{ IntUnitigMesg *mesg = (IntUnitigMesg *) vmesg;
  int			i;

  assert(strlen(mesg->consensus) == strlen(mesg->quality));
  assert((strlen(mesg->consensus) == mesg->length) || 
	 (strlen(mesg->consensus) == 0) );
  fprintf(fout,"{IUM\n");
  fprintf(fout,"acc:" F_IID "\n",mesg->iaccession);
# ifdef AS_ENABLE_SOURCE
  PutText(fout,"src:",mesg->source,FALSE);
# endif
  fprintf(fout,"cov:%.3f\n",mesg->coverage_stat);
  fprintf(fout,"sta:%c\n",mesg->status);
  fprintf(fout,"abp:" F_COORD "\n",mesg->a_branch_point);
  fprintf(fout,"bbp:" F_COORD "\n",mesg->b_branch_point);
  fprintf(fout,"len:" F_COORD "\n",mesg->length);
  PutText(fout,"cns:",mesg->consensus,TRUE);
  PutText(fout,"qlt:",mesg->quality,TRUE);
  fprintf(fout,"for:" F_S32 "\n",mesg->forced);
  fprintf(fout,"nfr:" F_S32 "\n",mesg->num_frags);
  for (i=0; i < mesg->num_frags; ++i)
    Write_IMP_Mesg(fout,&(mesg->f_list[i]));
  fprintf(fout,"}\n");
  return;
}

static void Write_IUL_Mesg(FILE *fout, void *vmesg)
{ IntUnitigLinkMesg *mesg = (IntUnitigLinkMesg *) vmesg;
  int		i, npairs;

  fprintf(fout,"{IUL\n");
  fprintf(fout,"ut1:" F_IID "\n",mesg->unitig1);
  fprintf(fout,"ut2:" F_IID "\n",mesg->unitig2);
  fprintf(fout,ORI_FORMAT "\n",mesg->orientation);
  fprintf(fout,"ovt:%c\n",mesg->overlap_type);
  fprintf(fout,"ipc:" F_S32 "\n",mesg->is_possible_chimera);
  fprintf(fout,"gui:" F_S32 "\n",mesg->includes_guide);
  fprintf(fout,MEA_OUT_FORMAT "\n",mesg->mean_distance);
  fprintf(fout,STD_OUT_FORMAT "\n",mesg->std_deviation);
  fprintf(fout,"num:" F_S32 "\n",mesg->num_contributing);
  fprintf(fout,"sta:%c\n",mesg->status);
  npairs = mesg->num_contributing;
  if (mesg->overlap_type != AS_NO_OVERLAP)
    --npairs;
  fprintf(fout,"jls:\n");
  for (i=0; i < npairs; ++i)
    fprintf(fout,F_IID "," F_IID ",%c\n",
            mesg->jump_list[i].in1,mesg->jump_list[i].in2,
            (char)(mesg->jump_list[i].type));
  fprintf(fout,"}\n");
  return;
}

static void Write_ICL_Mesg(FILE *fout, void *vmesg)
{ IntContigLinkMesg *mesg = (IntContigLinkMesg *) vmesg;
  int		i, npairs;

  fprintf(fout,"{ICL\n");
  fprintf(fout,"co1:" F_IID "\n",mesg->contig1);
  fprintf(fout,"co2:" F_IID "\n",mesg->contig2);
  fprintf(fout,ORI_FORMAT "\n",mesg->orientation);
  fprintf(fout,"ovt:%c\n",mesg->overlap_type);
  fprintf(fout,"ipc:" F_S32 "\n",mesg->is_possible_chimera);
  fprintf(fout,"gui:" F_S32 "\n",mesg->includes_guide);
  fprintf(fout,MEA_OUT_FORMAT "\n",mesg->mean_distance);
  fprintf(fout,STD_OUT_FORMAT "\n",mesg->std_deviation);
  fprintf(fout,"num:" F_S32 "\n",mesg->num_contributing);
  fprintf(fout,"sta:%c\n",mesg->status);
  npairs = mesg->num_contributing;
  if (mesg->overlap_type != AS_NO_OVERLAP)
    --npairs;
  fprintf(fout,"jls:\n");
  for (i=0; i < npairs; ++i)
    fprintf(fout,F_IID "," F_IID ",%c\n",
            mesg->jump_list[i].in1,mesg->jump_list[i].in2,
            (char)(mesg->jump_list[i].type));
  fprintf(fout,"}\n");
  return;
}

static void Write_ISL_Mesg(FILE *fout, void *vmesg)
{ InternalScaffoldLinkMesg *mesg = (InternalScaffoldLinkMesg *) vmesg;
  int		i, npairs;

  fprintf(fout,"{ISL\n");
  fprintf(fout,"sc1:" F_IID "\n",mesg->iscaffold1);
  fprintf(fout,"sc2:" F_IID "\n",mesg->iscaffold2);
  fprintf(fout,ORI_FORMAT "\n",mesg->orientation);
  fprintf(fout,"gui:" F_S32 "\n",mesg->includes_guide);
  fprintf(fout,MEA_OUT_FORMAT "\n",mesg->mean_distance);
  fprintf(fout,STD_OUT_FORMAT "\n",mesg->std_deviation);
  fprintf(fout,"num:" F_S32 "\n",mesg->num_contributing);
  npairs = mesg->num_contributing;
  assert(npairs > 0);
  fprintf(fout,"jls:\n");
  for (i=0; i < npairs; ++i)
    fprintf(fout,F_IID "," F_IID ",%c\n",
            mesg->jump_list[i].in1,mesg->jump_list[i].in2,
            (char)(mesg->jump_list[i].type));
  fprintf(fout,"}\n");
  return;
}

static void Write_AFG_Mesg(FILE *fout, void *vmesg)
{ AugFragMesg *mesg = (AugFragMesg *) vmesg;
  
  fprintf(fout,"{AFG\n");
  fprintf(fout,IACCS_FORMAT "\n",mesg->eaccession,mesg->iaccession);
  fprintf(fout,"mst:%c\n",mesg->mate_status);
  fprintf(fout,"chi:" F_S32 "\n",mesg->chimeric);
  fprintf(fout,"cha:" F_S32 "\n",mesg->chaff);
  fprintf(fout,CLR_FORMAT "\n",
          mesg->clear_rng.bgn,mesg->clear_rng.end);
  fprintf(fout,"}\n");
  return;
}

static void Write_ICP_Mesg(FILE *fout, IntContigPairs *mesg)
{
  fprintf(fout,"{ICP\n");
  fprintf(fout,"ct1:" F_IID "\n",mesg->contig1);
  fprintf(fout,"ct2:" F_IID "\n",mesg->contig2);
  fprintf(fout,MEA_OUT_FORMAT "\n",mesg->mean);
  fprintf(fout,STD_OUT_FORMAT "\n",mesg->stddev);
  fprintf(fout,ORI_FORMAT "\n",mesg->orient);
  fprintf(fout,"}\n");
  return;
}

static void Write_ISF_Mesg(FILE *fout, void *vmesg)
{ IntScaffoldMesg *mesg = (IntScaffoldMesg *) vmesg;
  int		i;
  int num = MAX(1, mesg->num_contig_pairs);
  fprintf(fout,"{ISF\n");
  fprintf(fout,"acc:" F_IID "\n", mesg->iaccession);
  fprintf(fout,"noc:" F_S32 "\n",mesg->num_contig_pairs);
  for (i=0; i < num; ++i)
    Write_ICP_Mesg(fout,&mesg->contig_pairs[i]);
  fprintf(fout,"}\n");
  return;
}

static void Write_IMD_Mesg(FILE *fout, void *vmesg)
{ IntMateDistMesg *mesg = (IntMateDistMesg *) vmesg;
  int		i;

  fprintf(fout,"{IMD\n");
  fprintf(fout,"ref:" F_IID "\n",mesg->refines);
  fprintf(fout,MEA_OUT_FORMAT "\n",mesg->mean);
  fprintf(fout,STD_OUT_FORMAT "\n",mesg->stddev);
  fprintf(fout,MINC_FORMAT "\n",mesg->min);
  fprintf(fout,MAXC_FORMAT "\n",mesg->max);
  fprintf(fout,"buc:" F_S32 "\n",mesg->num_buckets);
  fprintf(fout,"his:\n");
  for (i=0; i < mesg->num_buckets; ++i)
    fprintf(fout,F_S32 "\n",mesg->histogram[i]);
  fprintf(fout,"}\n");
  return;
}

static void Write_ICM_Mesg(FILE *fout, void *vmesg)
{ IntConConMesg *mesg = (IntConConMesg *) vmesg;
  int		i;

  fprintf(fout,"{ICM\n");
  fprintf(fout,"acc:" F_IID "\n",mesg->iaccession);
  fprintf(fout,"pla:%c\n",mesg->placed);
  fprintf(fout,"len:" F_COORD "\n",mesg->length);
  PutText(fout,"cns:",mesg->consensus,TRUE);
  PutText(fout,"qlt:",mesg->quality,TRUE);
  fprintf(fout,"for:" F_S32 "\n",mesg->forced);
  fprintf(fout,"npc:" F_S32 "\n",mesg->num_pieces);
  fprintf(fout,"nou:" F_S32 "\n",mesg->num_unitigs);
  if (!novar)
  {
     fprintf(fout,"nvr:" F_S32 "\n",mesg->num_vars);
     fflush(NULL); 
     for (i=0; i < mesg->num_vars; ++i)
       Write_IMV_Mesg(fout, &mesg->v_list[i]);
  }
  fflush(NULL);
  for (i=0; i < mesg->num_pieces; ++i)
    Write_IMP_Mesg(fout, &mesg->pieces[i]);
  fflush(NULL);
  for (i=0; i < mesg->num_unitigs; ++i)
    Write_IUP_Mesg(fout, &(mesg->unitigs[i]));
  fprintf(fout,"}\n");
  fflush(NULL);
  return;
}


static void Write_IAF_Mesg(FILE *fout, void *vmesg)
{ IntAugFragMesg *mesg = (IntAugFragMesg *) vmesg;
  
  fprintf(fout,"{IAF\n");
  fprintf(fout,"acc:" F_IID "\n",mesg->iaccession);
  fprintf(fout,TYP_FORMAT "\n",(char) mesg->type);
  fprintf(fout,"chi:" F_S32 "\n",mesg->chimeric);
  fprintf(fout,"cha:" F_S32 "\n",mesg->chaff);
  fprintf(fout,CLR_FORMAT "\n",
          mesg->clear_rng.bgn,mesg->clear_rng.end);
  fprintf(fout,"mst:%c\n",mesg->mate_status);
  fprintf(fout,"}\n");
  return;
}

/* Genome Snapshot output routines */
/***********************************/


static void Write_UPS_Mesg(FILE *fout, UnitigPos *up)
{ int i;

  fprintf(fout,"{UPS\n");
  fprintf(fout,TYP_FORMAT "\n",(char) up->type);
  fprintf(fout,"lid:" F_UID "\n",up->eident);
  fprintf(fout,POS2_FORMAT "\n",up->position.bgn,up->position.end);
  fprintf(fout,"dln:" F_S32 "\n",up->delta_length);
  fprintf(fout,"del:\n");
  if (up->delta_length > 0 ) {
    for(i=0; i < up->delta_length; i++) {
      fprintf(fout,F_S32,up->delta[i]);
      if (i%20 == 19) fprintf(fout,"\n");
      else  fprintf(fout," ");
    }
    if (up->delta_length%20 != 0) fprintf(fout,"\n");
  }
  fprintf(fout,"}\n");
  return;
}
static void Write_MPS_Mesg(FILE *fout, SnapMultiPos *mlp)
{ int i;

  fprintf(fout,"{MPS\n");
  fprintf(fout,TYP_FORMAT "\n",(char) mlp->type);
  fprintf(fout,"mid:" F_UID "\n",mlp->eident);
  #ifdef AS_ENABLE_SOURCE
  PutText(fout,"src:",mlp->source,FALSE);
  #endif
  fprintf(fout,POS2_FORMAT "\n",
          mlp->position.bgn,mlp->position.end);
  fprintf(fout,"dln:" F_S32 "\n",mlp->delta_length);
  fprintf(fout,"del:\n");
  if (mlp->delta_length > 0 ) {
    for(i=0; i < mlp->delta_length; i++) {
      fprintf(fout,F_S32,mlp->delta[i]);
      if (i%20 == 19) fprintf(fout,"\n");
      else  fprintf(fout," ");
    }
    if (mlp->delta_length%20 != 0) fprintf(fout,"\n");
  }
  fprintf(fout,"}\n");
  return;
}


static void Write_UTG_Mesg(FILE *fout, void *vmesg)
{ SnapUnitigMesg *mesg = (SnapUnitigMesg *) vmesg;
  int			i;

  fprintf(fout,"{UTG\n");
  fprintf(fout,IACCS_FORMAT "\n",
          mesg->eaccession,mesg->iaccession);
  #ifdef AS_ENABLE_SOURCE
  PutText(fout,"src:",mesg->source,FALSE);
  #endif
  fprintf(fout,"cov:%.3f\n",mesg->coverage_stat);
  fprintf(fout,"sta:%c\n",mesg->status);
  fprintf(fout,"abp:" F_COORD "\n",mesg->a_branch_point);
  fprintf(fout,"bbp:" F_COORD "\n",mesg->b_branch_point);
  fprintf(fout,"len:" F_COORD "\n",mesg->length);
  PutText(fout,"cns:",mesg->consensus,TRUE);
  PutText(fout,"qlt:",mesg->quality,TRUE);
  fprintf(fout,"for:" F_S32 "\n",mesg->forced);
  fprintf(fout,"nfr:" F_S32 "\n",mesg->num_frags);
  for (i=0; i < mesg->num_frags; ++i)
    Write_MPS_Mesg(fout,&(mesg->f_list[i]));
  fprintf(fout,"}\n");
  return;
}


static void Write_ULK_Mesg(FILE *fout, void *vmesg)
{ SnapUnitigLinkMesg *mesg = (SnapUnitigLinkMesg *) vmesg;
  int		i, npairs;

  fprintf(fout,"{ULK\n");
  fprintf(fout,"ut1:" F_UID "\n",mesg->eunitig1);
  fprintf(fout,"ut2:" F_UID "\n",mesg->eunitig2);
  fprintf(fout,ORI_FORMAT "\n",mesg->orientation);
  fprintf(fout,"ovt:%c\n",mesg->overlap_type);
  fprintf(fout,"ipc:" F_S32 "\n",mesg->is_possible_chimera);
  fprintf(fout,"gui:" F_S32 "\n",mesg->includes_guide);
  fprintf(fout,MEA_OUT_FORMAT "\n",mesg->mean_distance);
  fprintf(fout,STD_OUT_FORMAT "\n",mesg->std_deviation);
  fprintf(fout,"num:" F_S32 "\n",mesg->num_contributing);
  fprintf(fout,"sta:%c\n",mesg->status);
  npairs = mesg->num_contributing;
  if (mesg->overlap_type != AS_NO_OVERLAP)
    --npairs;
  fprintf(fout,"jls:\n");
  for (i=0; i < npairs; ++i)
    fprintf(fout, F_UID "," F_UID ",%c\n",
            mesg->jump_list[i].in1,
            mesg->jump_list[i].in2,
            (char)(mesg->jump_list[i].type));
  fprintf(fout,"}\n");
  return;
}


static void Write_CCO_Mesg(FILE *fout, void *vmesg)
{ SnapConConMesg *mesg = (SnapConConMesg *) vmesg;
  int		i;

  fprintf(fout,"{CCO\n");
  fprintf(fout,IACCS_FORMAT "\n",mesg->eaccession,mesg->iaccession);
  fprintf(fout,"pla:%c\n",mesg->placed);
  fprintf(fout,"len:" F_COORD "\n",mesg->length);
  PutText(fout,"cns:",mesg->consensus,TRUE);
  PutText(fout,"qlt:",mesg->quality,TRUE);
  fprintf(fout,"for:" F_S32 "\n",mesg->forced);
  fprintf(fout,"npc:" F_S32 "\n",mesg->num_pieces);
  fprintf(fout,"nou:" F_S32 "\n",mesg->num_unitigs);
  if (!novar)
  {
     fprintf(fout,"nvr:" F_S32 "\n",mesg->num_vars);
     for (i=0; i < mesg->num_vars; ++i)
       Write_VAR_Mesg(fout, &(mesg->vars[i]));
  }
  for (i=0; i < mesg->num_pieces; ++i)
    Write_MPS_Mesg(fout, &mesg->pieces[i]);
  for (i=0; i < mesg->num_unitigs; ++i)
    Write_UPS_Mesg(fout, &(mesg->unitigs[i]));
  fprintf(fout,"}\n");
  return;
}


static void Write_CLK_Mesg(FILE *fout, void *vmesg)
{ SnapContigLinkMesg *mesg = (SnapContigLinkMesg *) vmesg;
  int		i, npairs;

  fprintf(fout,"{CLK\n");
  fprintf(fout,"co1:" F_UID "\n",mesg->econtig1);
  fprintf(fout,"co2:" F_UID "\n",mesg->econtig2);
  fprintf(fout,ORI_FORMAT "\n",mesg->orientation);
  fprintf(fout,"ovt:%c\n",mesg->overlap_type);
  fprintf(fout,"ipc:" F_S32 "\n",mesg->is_possible_chimera);
  fprintf(fout,"gui:" F_S32 "\n",mesg->includes_guide);
  fprintf(fout,MEA_OUT_FORMAT "\n",mesg->mean_distance);
  fprintf(fout,STD_OUT_FORMAT "\n",mesg->std_deviation);
  fprintf(fout,"num:" F_S32 "\n",mesg->num_contributing);
  fprintf(fout,"sta:%c\n",mesg->status);
  npairs = mesg->num_contributing;
  if (mesg->overlap_type != AS_NO_OVERLAP)
    --npairs;
  fprintf(fout,"jls:\n");
  for (i=0; i < npairs; ++i)
    fprintf(fout, F_UID "," F_UID ",%c\n",
            mesg->jump_list[i].in1,
            mesg->jump_list[i].in2,
            (char)(mesg->jump_list[i].type));
  fprintf(fout,"}\n");
  return;
}

static void Write_SLK_Mesg(FILE *fout, void *vmesg)
{ SnapScaffoldLinkMesg *mesg = (SnapScaffoldLinkMesg *) vmesg;
  int		i, npairs;

  fprintf(fout,"{SLK\n");
  fprintf(fout,"sc1:" F_UID "\n",mesg->escaffold1);
  fprintf(fout,"sc2:" F_UID "\n",mesg->escaffold2);
  fprintf(fout,ORI_FORMAT "\n",mesg->orientation);
  fprintf(fout,"gui:" F_S32 "\n",mesg->includes_guide);
  fprintf(fout,MEA_OUT_FORMAT "\n",mesg->mean_distance);
  fprintf(fout,STD_OUT_FORMAT "\n",mesg->std_deviation);
  fprintf(fout,"num:" F_S32 "\n",mesg->num_contributing);
  npairs = mesg->num_contributing;
  assert(npairs > 0);
  fprintf(fout,"jls:\n");
  for (i=0; i < npairs; ++i)
    fprintf(fout, F_UID "," F_UID ",%c\n",
            mesg->jump_list[i].in1,
            mesg->jump_list[i].in2,
            (char)(mesg->jump_list[i].type));
  fprintf(fout,"}\n");
  return;
}


static void Write_CTP_Mesg(FILE *fout, SnapContigPairs *mesg)
{
  fprintf(fout,"{CTP\n");
  fprintf(fout,"ct1:" F_UID "\n",mesg->econtig1);
  fprintf(fout,"ct2:" F_UID "\n",mesg->econtig2);
  fprintf(fout,MEA_OUT_FORMAT "\n",mesg->mean);
  fprintf(fout,STD_OUT_FORMAT "\n",mesg->stddev);
  fprintf(fout,ORI_FORMAT "\n",mesg->orient);
  fprintf(fout,"}\n");
  return;
}

static void Write_SCF_Mesg(FILE *fout, void *vmesg)
{ SnapScaffoldMesg *mesg = (SnapScaffoldMesg *) vmesg;
  int		i;
  int num = MAX(1,mesg->num_contig_pairs);
  fprintf(fout,"{SCF\n");
  fprintf(fout,IACCS_FORMAT "\n",mesg->eaccession,mesg->iaccession);
  fprintf(fout,"noc:" F_S32 "\n",mesg->num_contig_pairs);
  for (i=0; i < num; ++i)
    Write_CTP_Mesg(fout,&mesg->contig_pairs[i]);
  fprintf(fout,"}\n");
  return;
}

static void Write_DSC_Mesg(FILE *fout, void *vmesg)
{ SnapDegenerateScaffoldMesg *mesg = (SnapDegenerateScaffoldMesg *) vmesg;
  fprintf(fout,"{DSC\n");
  fprintf(fout,EACC_FORMAT "\n",mesg->eaccession);
  fprintf(fout,"ctg:" F_UID "\n",mesg->econtig);
  fprintf(fout,"}\n");
  return;
}

static void Write_IDS_Mesg(FILE *fout, void *vmesg)
{ IntDegenerateScaffoldMesg *mesg = (IntDegenerateScaffoldMesg *) vmesg;
  fprintf(fout,"{IDS\n");
  fprintf(fout,"ctg:" F_IID "\n",mesg->icontig);
  fprintf(fout,"}\n");
  return;
}


static void Write_MDI_Mesg(FILE *fout, void *vmesg)
{ SnapMateDistMesg *mesg = (SnapMateDistMesg *) vmesg;
  int		i;

  fprintf(fout,"{MDI\n");
  fprintf(fout,IREFS_FORMAT "\n",mesg->erefines,mesg->irefines);
  fprintf(fout,MEA_OUT_FORMAT "\n",mesg->mean);
  fprintf(fout,STD_OUT_FORMAT "\n",mesg->stddev);
  fprintf(fout,MINC_FORMAT "\n",mesg->min);
  fprintf(fout,MAXC_FORMAT "\n",mesg->max);
  fprintf(fout,"buc:" F_S32 "\n",mesg->num_buckets);
  fprintf(fout,"his:\n");
  for (i=0; i < mesg->num_buckets; ++i)
    fprintf(fout,F_S32 "\n",mesg->histogram[i]);
  fprintf(fout,"}\n");
  return;
}

static void Write_BAT_Mesg(FILE *fout, void *vmesg){
  BatchMesg *mesg = (BatchMesg *)vmesg;
  fprintf(fout,"{BAT\n");
  fprintf(fout,"bna:%s\n",mesg->name);
  fprintf(fout,CRT_FORMAT "\n",mesg->created);
  fprintf(fout,EACC_FORMAT "\n",mesg->eaccession);
  PutText(fout,"com:",mesg->comment, FALSE);
  fprintf(fout,"}\n");

}
static void Write_IBA_Mesg(FILE *fout, void *vmesg){
  InternalBatchMesg *mesg = (InternalBatchMesg *)vmesg;
  fprintf(fout,"{IBA\n");
  fprintf(fout,"bna:%s\n",mesg->name);
  fprintf(fout,CRT_FORMAT "\n",mesg->created);
  fprintf(fout,IACCS_FORMAT "\n",mesg->eaccession, mesg->iaccession);
  PutText(fout,"com:",mesg->comment, FALSE);
  //  fprintf(stderr,"* Write_IBA_Mesg comment = %s\n", mesg->comment);
  fprintf(fout,"}\n");
}

static void Write_EOF_Mesg(FILE *fout, void *vmesg)
{
  EndOfFileMesg *mesg = (EndOfFileMesg *) vmesg;

  fprintf(fout, "{EOF\n");
  fprintf(fout, "sta:" F_S32 "\n", mesg->status);
  fprintf(fout, CRT_FORMAT "\n", mesg->created);
  PutText(fout, "com:", mesg->comment, FALSE);
  fprintf(fout, "}\n");
}


/******************** TRANSFER ROUTINES ***************************/

/*  Routines to transfer a message type through the phases.  */

/* the next 2 transfer routines exist only for backward compatability */
/* They are not needed! */

void Transfer_FRG_to_IFG_AS(FragMesg *fmg, InternalFragMesg *img)
{ memcpy(img,fmg,sizeof(FragMesg)); }

void Transfer_IFG_to_OFG_AS(InternalFragMesg *ifg, OFGMesg *ofg)
{ memcpy(ofg,ifg,sizeof(FragMesg)); }

void Transfer_DST_to_IDT_AS(DistanceMesg *dst, InternalDistMesg *idt)
{ idt->action      = dst->action;
  idt->eaccession  = dst->eaccession;
  idt->iaccession  = 0;
  idt->mean        = dst->mean;
  idt->stddev      = dst->stddev;
}

void Transfer_LKG_to_ILK_AS(LinkMesg *lkg, 
			    InternalLinkMesg *ilk){
  ilk->action      = lkg->action;
  ilk->type        = lkg->type;
  ilk->ifrag1      = 0;
  ilk->ifrag2      = 0;
  ilk->entry_time  = lkg->entry_time;
  ilk->idistance   = 0;
  ilk->link_orient = lkg->link_orient;
}


void AppendAuditLine_AS(AuditMesg *adt, AuditLine *adl,
                        time_t t, char *name, char *version, char *comment)
{ AuditLine *list;

// Skip over the entire list, stopping at the last element
 if(adt->list) {
   for (list = adt->list; list->next != NULL; list = list->next);

   list->next = adl;
 }else{
   adt->list = adl;
 }

 adl->complete = t;
 adl->name     = name;
 adl->version  = version;
 adl->comment  = comment;
 adl->next     = NULL;
}

/******************** FREE ROUTINES ***************************/

/*  Routines to free the second-level parts of a message.  */

static void Clear_ADT_Mesg(void *mesg, int typ)
{ AuditLine *alm, *nxt;

  for (alm = ((AuditMesg *) mesg)->list; alm != NULL; alm = nxt)
    { nxt = alm->next;
      free(alm->name);
      free(alm->version);
      free(alm->comment);
      free(alm);
    }
}

static void Clear_FRG_Mesg(void *mesg, int typ)
{ FragMesg *fgm = (FragMesg *) mesg;

  if (fgm->action == AS_ADD)
    { free(fgm->source);
      if (typ != MESG_OFG)
        { free(fgm->sequence);
          free(fgm->quality);
        }
    }
}

static void Clear_OVL_Mesg(void *mesg, int typ)
{
  free(((OverlapMesg *) mesg)->delta);
}

static void Clear_IUM_Mesg(void *vmesg, int typ)
{ 
  IntUnitigMesg *mesg = (IntUnitigMesg *) vmesg;
  int i;

  free(mesg->consensus);
  free(mesg->quality);
  for (i=0; i < mesg->num_frags; i++)
    free(mesg->f_list[i].delta);
  free (mesg->f_list);
//free (mesg->v_list);
}

static void Clear_IUL_Mesg(void *vmesg, int typ)
{
  IntUnitigLinkMesg *mesg = (IntUnitigLinkMesg *) vmesg;
  free(mesg->jump_list);
}

static void Clear_ICL_Mesg(void *vmesg, int typ)
{
  IntContigLinkMesg *mesg = (IntContigLinkMesg *) vmesg;
  free(mesg->jump_list);
}

static void Clear_ISF_Mesg(void *vmesg, int typ)
{
  IntScaffoldMesg *mesg = (IntScaffoldMesg *) vmesg;
  free(mesg->contig_pairs);
}

static void Clear_IMD_Mesg(void *vmesg, int typ)
{
  IntMateDistMesg *mesg = (IntMateDistMesg *) vmesg;
  free (mesg->histogram);
}

static void Clear_ICM_Mesg(void *vmesg, int typ)
{
  IntConConMesg *mesg = (IntConConMesg *) vmesg;
  int	i;

  free(mesg->consensus);
  free(mesg->quality);
  for (i=0; i < mesg->num_pieces; ++i)
      free(mesg->pieces[i].delta);
  free(mesg->pieces);
  free(mesg->unitigs);
  for (i=0; i < mesg->num_vars; ++i)
  {
      FREE(mesg->v_list[i].nr_conf_alleles);
      FREE(mesg->v_list[i].weights); 
      FREE(mesg->v_list[i].var_seq);
  }
  free(mesg->v_list);           
}



/******************** EXTERNAL ENTRY POINTS ***************************/

/*  Routines to duplicate the second-level parts of a message.  */

typedef struct {
  char const * const header;
  void *(*reader)(FILE *);
  void (*writer)(FILE *, void *);
  void (*clearer)(void *, int);
  size_t size;
} callrecord;

static const callrecord CallTable[] = {
  {"", NULL, NULL, NULL, 0l},
  {"{ADT", Read_ADT_Mesg, Write_ADT_Mesg, Clear_ADT_Mesg,   sizeof(AuditMesg) },
  {"{FRG", Read_FRG_Mesg, Write_FRG_Mesg, Clear_FRG_Mesg,   sizeof(FragMesg)  },
  {"{IFG", Read_IFG_Mesg, Write_IFG_Mesg, Clear_FRG_Mesg,   sizeof(InternalFragMesg) },
  {"", NULL, NULL, NULL, 0l },
  {"{OFG", Read_OFG_Mesg, Write_OFG_Mesg, Clear_FRG_Mesg,   sizeof(OFGMesg) },
  {"{LKG", Read_LKG_Mesg, Write_LKG_Mesg, NULL,             sizeof(LinkMesg) },
  {"{ILK", Read_ILK_Mesg, Write_ILK_Mesg, NULL,             sizeof(InternalLinkMesg) },
  {"{DST", Read_DST_Mesg, Write_DST_Mesg, NULL,             sizeof(DistanceMesg) },
  {"{IDT", Read_IDT_Mesg, Write_IDT_Mesg, NULL,             sizeof(InternalDistMesg) },
  {"", NULL, NULL, NULL, 0l },
  {"", NULL, NULL, NULL, 0l },
  {"", NULL, NULL, NULL, 0l },
  {"{OVL", Read_OVL_Mesg, Write_OVL_Mesg, Clear_OVL_Mesg,   sizeof(OverlapMesg) },
  {"", NULL, NULL, NULL, 0l },
  {"{UOM", Read_UOM_Mesg, Write_UOM_Mesg, NULL,             sizeof(UnitigOverlapMesg) },
  {"{IUM", Read_IUM_Mesg, Write_IUM_Mesg, Clear_IUM_Mesg,   sizeof(IntUnitigMesg) },
  {"{IUL", Read_IUL_Mesg, Write_IUL_Mesg, Clear_IUL_Mesg,   sizeof(IntUnitigLinkMesg) },
  {"{ICL", Read_ICL_Mesg, Write_ICL_Mesg, Clear_ICL_Mesg,   sizeof(IntContigLinkMesg) },
  {"{AFG", Read_AFG_Mesg, Write_AFG_Mesg, NULL,             sizeof(AugFragMesg) },
  {"{ISF", Read_ISF_Mesg, Write_ISF_Mesg, Clear_ISF_Mesg,   sizeof(IntScaffoldMesg) },
  {"{IMD", Read_IMD_Mesg, Write_IMD_Mesg, Clear_IMD_Mesg,   sizeof(IntMateDistMesg) },
  {"{IAF", Read_IAF_Mesg, Write_IAF_Mesg, NULL,     	    sizeof(IntAugFragMesg) },
  {"{UTG", Read_UTG_Mesg, Write_UTG_Mesg, NULL,             sizeof(SnapUnitigMesg) },
  {"{ULK", Read_ULK_Mesg, Write_ULK_Mesg, NULL,  	    sizeof(SnapUnitigLinkMesg) },
  {"{ICM", Read_ICM_Mesg, Write_ICM_Mesg, NULL,  	    sizeof(IntConConMesg) },
  {"{CCO", Read_CCO_Mesg, Write_CCO_Mesg, NULL,  	    sizeof(SnapConConMesg) },
  {"{CLK", Read_CLK_Mesg, Write_CLK_Mesg, NULL, 	    sizeof(SnapContigLinkMesg) },
  {"{SCF", Read_SCF_Mesg, Write_SCF_Mesg, NULL,  	    sizeof(SnapScaffoldMesg) },
  {"{MDI", Read_MDI_Mesg, Write_MDI_Mesg, NULL, 	    sizeof(SnapMateDistMesg) },
  {"{BAT", Read_BAT_Mesg, Write_BAT_Mesg, NULL, 	    sizeof(BatchMesg) },
  {"{IBA", Read_IBA_Mesg, Write_IBA_Mesg, NULL, 	    sizeof(InternalBatchMesg) },
  {"", NULL, NULL, NULL, 0l },
  {"", NULL, NULL, NULL, 0l },
  {"", NULL, NULL, NULL, 0l },
  {"", NULL, NULL, NULL, 0l },
  {"", NULL, NULL, NULL, 0l },
  {"", NULL, NULL, NULL, 0l },
  {"", NULL, NULL, NULL, 0l },
  {"", NULL, NULL, NULL, 0l },
  {"", NULL, NULL, NULL, 0l },
  {"{IDS", Read_IDS_Mesg, Write_IDS_Mesg, NULL,  	    sizeof(IntDegenerateScaffoldMesg) },
  {"{DSC", Read_DSC_Mesg, Write_DSC_Mesg, NULL,  	    sizeof(SnapDegenerateScaffoldMesg) },
  {"{SLK", Read_SLK_Mesg, Write_SLK_Mesg, NULL,  	    sizeof(SnapScaffoldLinkMesg) },
  {"{ISL", Read_ISL_Mesg, Write_ISL_Mesg, NULL,  	    sizeof(InternalScaffoldLinkMesg) },
  {"{FOM", Read_FOM_Mesg, Write_FOM_Mesg, NULL,             sizeof(FragOverlapMesg) },
  {"", NULL, NULL, NULL, 0l },
  {"", NULL, NULL, NULL, 0l },
  {"", NULL, NULL, NULL, 0l },
  {"", NULL, NULL, NULL, 0l },
  {"{EOF", Read_EOF_Mesg, Write_EOF_Mesg, NULL,             sizeof(EndOfFileMesg) }
};

static GenericMesg ReadMesg;

int GetMessageType(char *string){
  int t;
  for(t = 1; t < NUM_OF_REC_TYPES;t++){
    if (strncmp(string,CallTable[t].header + 1,3) == 0)
      return t;
  }
  return 0;
}

const char *GetMessageName(int type){
  if(type >= NUM_OF_REC_TYPES || type < 1)
    return NULL;

   return CallTable[type].header + 1;
}

int ReadProtoMesg_AS(FILE *fin, GenericMesg **pmesg)
{ int t;

  *pmesg = &ReadMesg;

  errno = 0;

  CurLine[MAX_LINE_LEN-2] = '\n';
  MemTop    = 0;
  do {
    //  Can't use ReadLine() here, because we want to return EOF if we
    //  read an empty line.
    LineNum++;
    if (fgets(CurLine,MAX_LINE_LEN,fin) == NULL)
      return (EOF);
  }
  while (CurLine[0] == '#');

  if (errno)
    fprintf(stderr, "ERROR: Read Failure looking for message type: %s\n",
            strerror(errno)), exit(1);

  for(t = 1; t <= NUM_OF_REC_TYPES; t++)
    if (strncmp(CurLine,CallTable[t].header,4) == 0)
      break;
  if (t > NUM_OF_REC_TYPES)
    { int len;

      len = strlen(CurLine)-1;
      if (CurLine[len] == '\n')
        CurLine[len] = '\0';
      fprintf(stderr,"ERROR: Unrecognized message type (%d > %d) \"%s\" at line %d\n",
	      t, NUM_OF_REC_TYPES, CurLine,LineNum);
      exit (1);
    }     
  Mcode = CallTable[t].header+1;

  ReadMesg.t = (MessageType) t;
  ReadMesg.m = CallTable[t].reader(fin);
  ReadMesg.s = MemTop;

  if (errno) {
    fprintf(stderr, "ERROR: Read Failure reading message %s: %s\n",
            CallTable[t].header+1, strerror(errno));
    exit(1);
  }

  return (0);
}

int WriteProtoMesg_AS(FILE *fout, GenericMesg *pmesg)
{ errno = 0;
  CallTable[pmesg->t].writer(fout,pmesg->m);
  if (errno) {
    fprintf(stderr, "ERROR: Write Failure: %s\n", strerror(errno));
    exit(1);
  }
  return (0);
}

int GetProtoLineNum_AS(void)     /* Current Line Number */
{ return (LineNum); }
