
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007-2008, J. CRaig Venter Institute.
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

const char *mainid = "$Id: sffToCA.c,v 1.25 2009-06-08 20:36:09 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <sys/stat.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "AS_UTL_reverseComplement.h"
#include "AS_UTL/AS_UTL_fasta.h"
#include "AS_GKP_include.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_encodeSequenceQuality.h"
#include "AS_ALN_bruteforcedp.h"


#define CLEAR_ALL        0x00
#define CLEAR_454        0x01
#define CLEAR_N          0x02
#define CLEAR_PAIR_N     0x04
#define CLEAR_DISCARD_N  0x08
#define CLEAR_ERRR       0xff

#define TRIM_NONE        0
#define TRIM_SOFT        1
#define TRIM_HARD        2
#define TRIM_CHOP        3
#define TRIM_ERRR        9

#define AS_LINKER_MAX_SEQS       50
#define AS_LINKER_CUSTOM_OFFSET  24

GateKeeperStore *gkpStore    = NULL;
FILE            *logFile     = NULL;
uint32           clearAction = CLEAR_454;
uint32           clearSet    = 0;
uint32           trimAction  = TRIM_HARD;

typedef struct {
  uint32   numReadsInSFF;

  uint32   numReadsWithPartialLinker;  //  A single unmated read
  uint32   numReadsWithFullLinker;     //  A pair of reads joined by a mate
  uint32   numReadsWithNoLinker;       //  A single unmated read

  uint32   numTooLong;
  uint32   numTooShort;
  uint32   numWithUnallowedN;
  uint32   numDeDuped;
} statistics;

statistics st = {0};

static
void
writeStatistics(FILE *statOut) {
  fprintf(statOut, "numReadsInSFF             "F_U32"\n", st.numReadsInSFF);
  fprintf(statOut, "numReadsWithPartialLinker "F_U32"\n", st.numReadsWithPartialLinker);
  fprintf(statOut, "numReadsWithFullLinker    "F_U32"\n", st.numReadsWithFullLinker);
  fprintf(statOut, "numReadsWithNoLinker      "F_U32"\n", st.numReadsWithNoLinker);
  fprintf(statOut, "numTooLong                "F_U32"\n", st.numTooLong);
  fprintf(statOut, "numTooShort               "F_U32"\n", st.numTooShort);
  fprintf(statOut, "numWithUnallowedN         "F_U32"\n", st.numWithUnallowedN);
  fprintf(statOut, "numDeDuped                "F_U32"\n", st.numDeDuped);
}



static
void
addReadToStore(GateKeeperStore *gkp,
               fragRecord      *fr) {
  int        encodedlength;

  fr->gkfr.readIID = getLastElemStore(gkpStore->frg) + 1;

  fr->gkfr.seqLen = strlen(fr->seq);
  fr->gkfr.hpsLen = 0;
  fr->gkfr.srcLen = 0;

  fr->gkfr.seqOffset = getLastElemStore(gkpStore->seq) + 1;
  fr->gkfr.qltOffset = getLastElemStore(gkpStore->qlt) + 1;
  fr->gkfr.hpsOffset = getLastElemStore(gkpStore->hps) + 1;
  fr->gkfr.srcOffset = getLastElemStore(gkpStore->src) + 1;

  setGatekeeperUIDtoIID(gkpStore, fr->gkfr.readUID, fr->gkfr.readIID, AS_IID_FRG);
  appendIndexStore(gkpStore->frg, &fr->gkfr);

  encodedlength = encodeSequence(fr->enc, fr->seq);
  appendStringStore(gkpStore->seq, fr->enc, encodedlength);

  encodeSequenceQuality(fr->enc, fr->seq, fr->qlt);
  appendStringStore(gkpStore->qlt, fr->enc, fr->gkfr.seqLen);

  appendStringStore(gkpStore->hps, NULL, 0);
  appendStringStore(gkpStore->src, NULL, 0);

  gkpStore->gkp.sffLoaded++;
}



//  Ugly hack that disappears in CA 6.x
static
int
fileExists(const char *path,
           int directory,
           int readwrite) {
  struct stat  s;
  int          r;

  errno = 0;
  r = stat(path, &s);
  if (errno) {
    //fprintf(stderr, "Failed to stat() file '%s'; assumed to not exist.\n", path);
    return(0);
  }

  if (directory == 1) {
    if ((readwrite == 0) &&
        (s.st_mode & S_IFDIR) &&
        (s.st_mode & (S_IRUSR | S_IRGRP | S_IROTH)) &&
        (s.st_mode & (S_IXUSR | S_IXGRP | S_IXOTH)))
      return(1);
    if ((readwrite == 1) &&
        (s.st_mode & S_IFDIR) &&
        (s.st_mode & (S_IRUSR | S_IRGRP | S_IROTH)) &&
        (s.st_mode & (S_IWUSR | S_IWGRP | S_IWOTH)) &&
        (s.st_mode & (S_IXUSR | S_IXGRP | S_IXOTH)))
      return(1);
  }

  if (directory == 0) {
    if ((readwrite == 0) &&
        (s.st_mode & (S_IRUSR | S_IRGRP | S_IROTH)))
      return(1);
    if ((readwrite == 1) &&
        (s.st_mode & (S_IRUSR | S_IRGRP | S_IROTH)) &&
        (s.st_mode & (S_IWUSR | S_IWGRP | S_IWOTH)))
      return(1);
  }

  return(0);
}



////////////////////////////////////////////////////////////////////////////////
//
//  Reads an SFF file, inserts all the reads (that are of the proper length)
//  into a gkpStore.
//

typedef struct {
  //  The next block is read in one swoop from the sff file.  DO NOT MODIFY!
  uint32   magic_number;
  char     version[4];
  uint64   index_offset;
  uint32   index_length;
  uint32   number_of_reads;
  uint16   header_length;
  uint16   key_length;
  uint16   number_of_flows_per_read;
  uint8    flowgram_format_code;

  char    *flow_chars;        //  h->number_of_flows_per_read
  char    *key_sequence;      //  h->key_length

  void    *data_block;
  uint32   data_block_len;

  uint32   swap_endianess;
} sffHeader;


typedef struct {
  //  The next block is read in one swoop from the sff file.  DO NOT MODIFY!
  uint32    magic_number;
  char      version[4];
  uint32    manifest_length;
  uint32    nothing;

  char     *manifest;
} sffManifest;


typedef struct {
  //  The next block is read in one swoop from the sff file.  DO NOT MODIFY!
  uint16   read_header_length;
  uint16   name_length;
  uint32   number_of_bases;
  uint16   clip_quality_left;
  uint16   clip_quality_right;
  uint16   clip_adapter_left;
  uint16   clip_adapter_right;

  char    *name;                 //  r->name_length
  uint16  *flowgram_values;      //  h->number_of_flows_per_read
  uint8   *flow_index_per_base;  //  r->number_of_bases
  char    *bases;                //  r->number_of_bases
  uint8   *quality_scores;       //  r->number_of_bases
  char    *quality;              //  quality_scores converted to CA-format qv

  int      final_length;         //  trimmed, processed read, ready for
  char    *final_bases;          //  loading.  NOT zero terminated.
  char    *final_quality;        //  DO NOT ZERO TERMINATE.

  void    *data_block;
  uint32   data_block_len;
} sffRead;


inline
uint64
uint64Swap(uint64 x) {
  x = ((x >>  8) & 0x00ff00ff00ff00ffLLU) | ((x <<  8) & 0xff00ff00ff00ff00LLU);
  x = ((x >> 16) & 0x0000ffff0000ffffLLU) | ((x << 16) & 0xffff0000ffff0000LLU);
  x = ((x >> 32) & 0x00000000ffffffffLLU) | ((x << 32) & 0xffffffff00000000LLU);
  return(x);
}

inline
uint32
uint32Swap(uint32 x) {
  x = ((x >>  8) & 0x00ff00ff) | ((x <<  8) & 0xff00ff00);
  x = ((x >> 16) & 0x0000ffff) | ((x << 16) & 0xffff0000);
  return(x);
}

inline
uint16
uint16Swap(uint16 x) {
  x = ((x >>  8) & 0x000000ff) | ((x <<  8) & 0x0000ff00);
  return(x);
}



static
void
readsff_manifest(FILE *sff, sffHeader *h, sffManifest *m) {

  if (h->index_length == 0)
    //  No manifest.
    return;

  if (AS_UTL_ftell(sff) != h->index_offset)
    //  Not at the manifest.
    return;

  if (m->manifest)
    //  Already got it?!
    return;

  AS_UTL_safeRead(sff, m, "readsff_manifest", sizeof(char), 16);

  if (h->swap_endianess) {
    m->magic_number    = uint32Swap(m->magic_number);
    m->manifest_length = uint32Swap(m->manifest_length);
  }

  m->manifest = (char *)safe_malloc(sizeof(char) * m->manifest_length + 1);
  AS_UTL_safeRead(sff, m->manifest, "readsff_manifest_text", sizeof(char), m->manifest_length);

  m->manifest[m->manifest_length] = 0;

  //  We only read the manifest.  There is still an index in there.

  uint64  padding_length = h->index_length - 16 - m->manifest_length;
  if (padding_length > 0) {
    //fprintf(stderr, "manifest pad "F_U64"\n", padding_length);
    char *junk = (char *)safe_malloc(sizeof(char) * padding_length);
    AS_UTL_safeRead(sff, junk, "readsff_manifest_pad", sizeof(char), padding_length);
    safe_free(junk);
  }
}


static
void
readsff_header(FILE *sff, sffHeader *h, sffManifest *m) {

  AS_UTL_safeRead(sff, h, "readsff_header_1", 31, 1);

  if (h->magic_number != 0x2e736666) {
    h->swap_endianess           = 1;
    h->magic_number             = uint32Swap(h->magic_number);
    h->index_offset             = uint64Swap(h->index_offset);
    h->index_length             = uint32Swap(h->index_length);
    h->number_of_reads          = uint32Swap(h->number_of_reads);
    h->header_length            = uint16Swap(h->header_length);
    h->key_length               = uint16Swap(h->key_length);
    h->number_of_flows_per_read = uint16Swap(h->number_of_flows_per_read);
  }

  assert(h->magic_number == 0x2e736666);

  uint32 newlen = h->number_of_flows_per_read + h->key_length + 2;
  if (h->data_block_len < newlen) {
    h->data_block_len = newlen;
    h->data_block     = safe_realloc(h->data_block, h->data_block_len);
  }

  memset(h->data_block, 0, h->data_block_len);

  h->flow_chars   = h->data_block;
  h->key_sequence = h->data_block + (h->number_of_flows_per_read + 1) * sizeof(char);

  AS_UTL_safeRead(sff,  h->flow_chars,   "readsff_header_2", sizeof(char), h->number_of_flows_per_read);
  AS_UTL_safeRead(sff,  h->key_sequence, "readsff_header_3", sizeof(char), h->key_length);

  uint64  padding_length = h->header_length - 31 - h->number_of_flows_per_read - h->key_length;
  if (padding_length > 0) {
    //fprintf(stderr, "header pad "F_U64"\n", padding_length);
    char *junk = (char *)safe_malloc(sizeof(char) * padding_length);
    AS_UTL_safeRead(sff, junk, "readsff_header_4", sizeof(char), padding_length);
    safe_free(junk);
  }

  //  The spec says the index might be here, however, all files I've
  //  seen have the index at the end of the file.
  //
  readsff_manifest(sff, h, m);
}


static
void
readsff_read(FILE *sff, sffHeader *h, sffRead *r) {

  AS_UTL_safeRead(sff, r, "readsff_read_1", 16, 1);

  if (h->swap_endianess) {
    r->read_header_length = uint16Swap(r->read_header_length);
    r->name_length        = uint16Swap(r->name_length);
    r->number_of_bases    = uint32Swap(r->number_of_bases);
    r->clip_quality_left  = uint16Swap(r->clip_quality_left);
    r->clip_quality_right = uint16Swap(r->clip_quality_right);
    r->clip_adapter_left  = uint16Swap(r->clip_adapter_left);
    r->clip_adapter_right = uint16Swap(r->clip_adapter_right);
  }

  //  Can you say UGLY?  Hey, it's a lot better than what I originally came up with.

  uint32 ss[6];
  ss[0] = (r->name_length + 1)          * sizeof(char);
  ss[1] = (h->number_of_flows_per_read) * sizeof(uint16) + ss[0];
  ss[2] = (r->number_of_bases)          * sizeof(uint8)  + ss[1];
  ss[3] = (r->number_of_bases + 1)      * sizeof(char)   + ss[2];
  ss[4] = (r->number_of_bases + 1)      * sizeof(uint8)  + ss[3];
  ss[5] = (r->number_of_bases + 1)      * sizeof(char)   + ss[4];

  if (r->data_block_len < ss[5]) {
    r->data_block_len = ss[5];
    r->data_block     = safe_realloc(r->data_block, r->data_block_len);
  }

  memset(r->data_block, 0, r->data_block_len);

  r->name                 = r->data_block;
  r->flowgram_values      = r->data_block + ss[0];
  r->flow_index_per_base  = r->data_block + ss[1];
  r->bases                = r->data_block + ss[2];
  r->quality_scores       = r->data_block + ss[3];
  r->quality              = r->data_block + ss[4];

  AS_UTL_safeRead(sff, r->name, "readsff_read_2", sizeof(char), r->name_length);
  r->name[r->name_length] = 0;

  uint64  padding_length = r->read_header_length - 16 - r->name_length;
  if (padding_length > 0) {
    //fprintf(stderr, "read pad 1 "F_U64"\n", padding_length);
    uint64  junk;
    AS_UTL_safeRead(sff, &junk, "readsff_read_3", sizeof(char), padding_length);
  }

  AS_UTL_safeRead(sff, r->flowgram_values,     "readsff_read_4", sizeof(uint16), h->number_of_flows_per_read);
  AS_UTL_safeRead(sff, r->flow_index_per_base, "readsff_read_5", sizeof(uint8),  r->number_of_bases);
  AS_UTL_safeRead(sff, r->bases,               "readsff_read_6", sizeof(char),   r->number_of_bases);
  AS_UTL_safeRead(sff, r->quality_scores,      "readsff_read_7", sizeof(uint8),  r->number_of_bases);

  int i;
  for (i=0; i<r->number_of_bases; i++)
    r->quality[i] = r->quality_scores[i] + '0';

  r->bases[r->number_of_bases] = 0;
  r->quality[r->number_of_bases] = 0;

  //  The padding_length is the number of bytes to make the above four
  //  chunks of data be of size that is divisible by 8.  The
  //  padding_length we compute directly below is the number of bytes
  //  we read past the last multiple of 8, and if that is non-zero, we
  //  need to read 8-padding_length bytes.
  //
  padding_length = (h->number_of_flows_per_read * sizeof(uint16) +
                    r->number_of_bases * sizeof(uint8) +
                    r->number_of_bases * sizeof(char) +
                    r->number_of_bases * sizeof(uint8)) % 8;
  if (padding_length > 0) {
    //fprintf(stderr, "read pad 2 "F_U64"\n", 8-padding_length);
    char *junk = (char *)safe_malloc(sizeof(char) * padding_length);
    AS_UTL_safeRead(sff, junk, "readsff_read_8", sizeof(char), 8 - padding_length);
    safe_free(junk);
  }
}


static
void
processRead(sffHeader *h,
            sffRead   *r, fragRecord *fr) {
  int  which;

  st.numReadsInSFF++;

  clearGateKeeperFragmentRecord(&fr->gkfr);

  //  Construct a UID from the 454 read name
  fr->gkfr.readUID = AS_UID_load(r->name);

  //  Read already loaded?  Can't load again.  Set UID;s and IID's
  //  to zero to indicate this -- we'll catch it at the end.
  //
  if (getGatekeeperUIDtoIID(gkpStore, fr->gkfr.readUID, NULL)) {
    fprintf(stderr, "ERROR:  Read '%s' already exists.\n", AS_UID_toString(fr->gkfr.readUID));
    fr->gkfr.readUID = AS_UID_undefined();
    return;
  }

  fr->gkfr.libraryIID  = 1;
  fr->gkfr.orientation = AS_READ_ORIENT_UNKNOWN;

  ////////////////////////////////////////
  //
  //  Chop off any N's at the end of the read.  Titanium likes to do
  //  this to us.
  //
  while ((r->bases[r->number_of_bases-1] == 'N') ||
         (r->bases[r->number_of_bases-1] == 'n')) {
    r->number_of_bases--;
    r->bases[r->number_of_bases] = 0;
    r->quality[r->number_of_bases] = 0;
  }

  if (r->clip_adapter_right > r->number_of_bases)
    r->clip_adapter_right = r->number_of_bases;

  if (r->clip_quality_right > r->number_of_bases)
    r->clip_quality_right = r->number_of_bases;

  ////////////////////////////////////////
  //
  //  Attempt to make sense of 454 supplied clear ranges.
  //
  //  These are base-based.  If either value is 0, that means the
  //  value was not computed.  In that case, we set it to the extent
  //  (max or min).
  //
  int clq = h->key_length;
  int crq = r->number_of_bases;

  int cla = h->key_length;
  int cra = r->number_of_bases;

  if (clearAction & CLEAR_454) {
    //  Left point should be zero or after the key
    assert((r->clip_quality_left == 0) || (h->key_length <= r->clip_quality_left));
    assert((r->clip_adapter_left == 0) || (h->key_length <= r->clip_adapter_left));

    //  Right point should be zero or before the end
    assert((r->clip_quality_right == 0) || (r->clip_quality_left <= r->number_of_bases));
    assert((r->clip_adapter_right == 0) || (r->clip_adapter_left <= r->number_of_bases));

    clq = MAX(r->clip_quality_left,  h->key_length + 1) - 1;
    cla = MAX(r->clip_adapter_left,  h->key_length + 1) - 1;

    crq = (r->clip_quality_right > 0) ? r->clip_quality_right : r->number_of_bases;
    cra = (r->clip_adapter_right > 0) ? r->clip_adapter_right : r->number_of_bases;
  }

  ////////////////////////////////////////
  //
  //  Find the CLEAR_N and CLEAR_PAIR_N points.  If we're allowing the
  //  use of the 454 clear ranges, don't check for N's before that.
  //  qThis assumes that clq and cla are NOT set if the 454 clear is
  //  NOT used.
  //
  int  cln = h->key_length;
  int  crn = r->number_of_bases;
  int  frn = r->number_of_bases;  //  first-n

  if ((clearAction & CLEAR_N) ||
      (clearAction & CLEAR_DISCARD_N)) {
    int  f  = 0;
    int  b  = MAX(clq, cla);
    int  e  = r->number_of_bases;
    char *s = r->bases;

    for (f=b; f<e; f++)
      if ((s[f] == 'n') || (s[f] == 'N'))
        break;

    if (clearAction & CLEAR_N) {
      cln = b;
      crn = f;
    }

    frn = f;
  }

  int  clp = h->key_length;
  int  crp = r->number_of_bases;

  if (clearAction & CLEAR_PAIR_N) {
    int   f = 0;
    int   b = MAX(clq, cla);
    int   e = r->number_of_bases - 1;
    char *s = r->bases;

    for (f=b; f<e; f++)
      if (((s[f+0] == 'n') || (s[f+0] == 'N')) &&
          ((s[f+1] == 'n') || (s[f+1] == 'N')))
        break;

    clp = b;
    crp = f;
  }


  ////////////////////////////////////////
  //
  //  Make sense of all these by blindly intersecting them together.
  //
  int clf = MAX(MAX(clq, cla), MAX(cln, clp));
  int crf = MIN(MIN(crq, cra), MIN(crn, crp));


  ////////////////////////////////////////
  //
  //  If told to, and there is still an N in the sequence, trash the
  //  whole thing.
  //
  if ((clearAction & CLEAR_DISCARD_N) && (frn < crf)) {
    st.numWithUnallowedN++;

    if (logFile)
      fprintf(logFile, "Read '%s' contains an N at position %d.  Read deleted.\n",
              AS_UID_toString(fr->gkfr.readUID), crn);

    fr->gkfr.readUID = AS_UID_undefined();
  }


  ////////////////////////////////////////
  //
  //  Now, decide how to set the clear ranges.  Ignore, soft, hard, chop?
  //
  switch (trimAction) {
    int which;

    case TRIM_NONE:
      //  Set clear range to the whole untrimmed read.
      for (which=0; which <= AS_READ_CLEAR_LATEST; which++) {
        fr->gkfr.clearBeg[which] = h->key_length;
        fr->gkfr.clearEnd[which] = r->number_of_bases;
      }

      //  No point setting the max.
      fr->gkfr.clearBeg[AS_READ_CLEAR_MAX] = 1;
      fr->gkfr.clearEnd[AS_READ_CLEAR_MAX] = 0;
      break;

    case TRIM_SOFT:
      //  Set clear ranges to whatever we discovered above, but do not
      //  limit OBT to anything.
      for (which=0; which <= AS_READ_CLEAR_LATEST; which++) {
        fr->gkfr.clearBeg[which] = clf;
        fr->gkfr.clearEnd[which] = crf;
      }

      //  No point setting the max.
      fr->gkfr.clearBeg[AS_READ_CLEAR_MAX] = 1;
      fr->gkfr.clearEnd[AS_READ_CLEAR_MAX] = 0;
      break;

    case TRIM_HARD:
      //  Set clear ranges to whatever we discovered above, and limit
      //  OBT to those ranges.
      for (which=0; which <= AS_READ_CLEAR_LATEST; which++) {
        fr->gkfr.clearBeg[which] = clf;
        fr->gkfr.clearEnd[which] = crf;
      }

      //  The only one that sets the max.
      fr->gkfr.clearBeg[AS_READ_CLEAR_MAX] = clf;
      fr->gkfr.clearEnd[AS_READ_CLEAR_MAX] = crf;
      break;

    case TRIM_CHOP:
      //  Rewrite the read to remove the non-clear sequence.  We keep
      //  in the usually four base long key at the start (with some
      //  amount of pain since clf,clr include those four bases).

      memmove(r->bases   + h->key_length, r->bases   + clf, sizeof(char) * (crf - clf));
      memmove(r->quality + h->key_length, r->quality + clf, sizeof(char) * (crf - clf));

      r->number_of_bases = h->key_length + crf - clf;

      r->bases  [h->key_length + crf - clf] = 0;
      r->quality[h->key_length + crf - clf] = 0;

      for (which=0; which <= AS_READ_CLEAR_LATEST; which++) {
        fr->gkfr.clearBeg[which] = h->key_length;
        fr->gkfr.clearEnd[which] = h->key_length + crf - clf;
      }

      //  No point setting the max.
      fr->gkfr.clearBeg[AS_READ_CLEAR_MAX] = 1;
      fr->gkfr.clearEnd[AS_READ_CLEAR_MAX] = 0;
      break;
    default:
      break;
  }

  //  There is no vector clear defined for 454 reads.
  fr->gkfr.clearBeg[AS_READ_CLEAR_VEC] = 1;
  fr->gkfr.clearEnd[AS_READ_CLEAR_VEC] = 0;

  //  Make sure that we aren't claiming contamination in the middle.
  fr->gkfr.contaminationBeg = 0;
  fr->gkfr.contaminationEnd = 0;

  ////////////////////////////////////////////////////////////////////////////////

  //  Reads too long cannot be physically loaded into the store.  Trim
  //  off anything over the maximum size.
  //
  if (r->number_of_bases - h->key_length > AS_READ_MAX_LEN) {
    st.numTooLong++;

    if (logFile)
      fprintf(logFile, "Read '%s' of length %d is too long.  Truncating to %d bases.\n",
              r->name, r->number_of_bases - h->key_length, AS_READ_MAX_LEN);

    r->number_of_bases = AS_READ_MAX_LEN;

    r->bases  [AS_READ_MAX_LEN + h->key_length] = 0;
    r->quality[AS_READ_MAX_LEN + h->key_length] = 0;

    for (which=0; which <= AS_READ_CLEAR_LATEST; which++) {
      if (fr->gkfr.clearBeg[which] > AS_READ_MAX_LEN)
        fr->gkfr.clearBeg[which] = AS_READ_MAX_LEN;
      if (fr->gkfr.clearEnd[which] > AS_READ_MAX_LEN)
        fr->gkfr.clearEnd[which] = AS_READ_MAX_LEN;
    }
  }

  //  Likewise, reads too short will never be of any use, and they're
  //  not loaded.  We'll give the benefit of the doubt and use the
  //  read length here, not the trimmed length.
  //
  if (r->number_of_bases - h->key_length < AS_READ_MIN_LEN) {
    st.numTooShort++;

    if (logFile)
      fprintf(logFile, "Read '%s' of length %d clear %d,%d is too short.  Not loaded.\n",
              r->name,
              r->number_of_bases - h->key_length,
              fr->gkfr.clearEnd[AS_READ_CLEAR_LATEST] - h->key_length,
              fr->gkfr.clearBeg[AS_READ_CLEAR_LATEST] - h->key_length);

    fr->gkfr.readUID = AS_UID_undefined();
  }

  //  Finally, adjust everything to remove the key_length bases from the start.
  //
  for (which=0; which <= AS_READ_CLEAR_LATEST; which++) {
    if (fr->gkfr.clearBeg[which] < fr->gkfr.clearEnd[which]) {
      fr->gkfr.clearBeg[which] -= h->key_length;
      fr->gkfr.clearEnd[which] -= h->key_length;
    }
  }

  r->final_bases    = r->bases   + h->key_length;
  r->final_quality  = r->quality + h->key_length;
  r->final_length   = r->number_of_bases - h->key_length;

  //  Copy sequence to the fragRecord
  //
  memcpy(fr->seq, r->final_bases,   sizeof(char) * (r->final_length + 1));
  memcpy(fr->qlt, r->final_quality, sizeof(char) * (r->final_length + 1));
}



int
loadSFF(char *sffName) {
  FILE                      *sff  = NULL;
  int                        fic  = 0;
  sffHeader                  h    = {0};
  sffManifest                m    = {0};
  sffRead                    r    = {0};
  fragRecord                 fr   = {0};
  int                        rn   = 0;

  errno = 0;

  fprintf(stderr, "loadSFF()-- Loading '%s'.\n", sffName);

  if        (strcasecmp(sffName + strlen(sffName) - 3, ".gz") == 0) {
    char  cmd[1024];
    sprintf(cmd, "gzip -dc %s", sffName);
    sff   = popen(cmd, "r");
    fic   = 1;
  } else if (strcasecmp(sffName + strlen(sffName) - 4, ".bz2") == 0) {
    char  cmd[1024];
    sprintf(cmd, "bzip2 -dc %s", sffName);
    sff   = popen(cmd, "r");
    fic   = 1;
  } else {
    sff   = fopen(sffName, "r");
    fic   = 0;
  }
  if (errno)
    fprintf(stderr, "ERROR!  Failed to open '%s': %s\n", sffName, strerror(errno)), exit(1);


  readsff_header(sff, &h, &m);

  for (rn=0; rn < h.number_of_reads; rn++) {
    readsff_read(sff, &h, &r);
    processRead(&h, &r, &fr);
    if (AS_UID_compare(fr.gkfr.readUID, AS_UID_undefined()) != 0)
      addReadToStore(gkpStore, &fr);
  }

  //  Read the manifest if we haven't already done so.
  if (m.manifest_length == 0)
    readsff_manifest(sff, &h, &m);

  //  Make sure that the reads have been rescored.
  if (m.manifest != NULL)
    if (strstr(m.manifest, "<qualityScoreVersion>1.1.03</qualityScoreVersion>") == NULL)
      fprintf(stderr, "WARNING:  Fragments not rescored!\n");

  errno = 0;

  if (fic)
    pclose(sff);
  else
    fclose(sff);
  if (errno)
    fprintf(stderr, "WARNING!  Failed to close '%s': %s\n", sffName, strerror(errno));

  safe_free(h.data_block);
  safe_free(m.manifest);
  safe_free(r.data_block);

  return(0);
}





////////////////////////////////////////////////////////////////////////////////
//
//  Removes all reads that are a perfect prefix of some other read.
//
//  The algorithm builds a 64-bit value from the first N bases (N <=
//  AS_FRAG_MIN_LEN, currently 48 <= 64), sorts the hashes, then
//  examines any clique of hash collisions for perfect prefixes.

typedef struct {
  uint64    hash;
  uint32    iid;
} fragHash;

static
int
fragHashCompare(const void *a, const void *b) {
  fragHash const *A = (fragHash const *)a;
  fragHash const *B = (fragHash const *)b;

  if (A->hash < B->hash) return(-1);
  if (A->hash > B->hash) return( 1);
  return(0);
}

void
removeDuplicateReads(void) {

  uint32        firstElem = getFirstElemFragStore(gkpStore);
  uint32        lastElem  = getLastElemFragStore(gkpStore) + 1;

  uint32        fragsLen = 0;
  uint32        fragsMax = 16384;
  fragRecord    fr;
  fragRecord   *frags    = (fragRecord *)safe_malloc(sizeof(fragRecord) * fragsMax);

  fragHash   *fh    = safe_malloc(sizeof(fragHash) * (lastElem - firstElem + 1));
  uint32      fhLen = 0;

  uint64 map[256] = { 0 };
  uint32 s, h, n;

  uint32 thisElem;

  for (s=0; s<256; s++)
    map[s] = 0;
  map['A'] = map['a'] = 0x00;
  map['C'] = map['c'] = 0x01;
  map['G'] = map['g'] = 0x02;
  map['T'] = map['t'] = 0x03;

  fprintf(stderr, "removeDuplicateReads()-- from %d to %d\n", firstElem, lastElem);

  for (thisElem=firstElem; thisElem<lastElem; thisElem++) {
    getFrag(gkpStore, thisElem, &fr, FRAG_S_INF | FRAG_S_SEQ);

    char *seq1      = getFragRecordSequence(&fr);

    uint32 seqLen   = getFragRecordSequenceLength(&fr);
    uint64 hash     = 0;

    assert(seqLen >= 48);

#if 1
    //  Our "hash" is just the spaced seed "101" (repeating).  It
    //  covers the first 48 bases, picking out 32.
    //
    for (s=0, n=0; n<16; n++) {
      hash <<= 2;
      hash  |= map[seq1[s]];
      s++;
      s++;
      hash <<= 2;
      hash  |= map[seq1[s]];
      s++;
    }
#else
    //  Our "hash" is just the spaced seed "1010" (repeating).  It
    //  covers the first 64 bases, picking out 32.
    //
    for (s=0, n=0; n<16; n++) {
      hash <<= 2;
      hash  |= map[seq1[s]];
      s++;
      s++;
      hash <<= 2;
      hash  |= map[seq1[s]];
      s++;
      s++;
    }
#endif

    fh[fhLen].hash     = hash;
    fh[fhLen].iid      = thisElem;

    fhLen++;
  }

  qsort(fh, fhLen, sizeof(fragHash), fragHashCompare);

  uint32  beg = 0;
  uint32  end = 0;

  while (beg < fhLen) {

    //  We DO need to examine the whole clique (pairwise).  We cannot
    //  simply sort by size, because if we get three frags of the same
    //  size, it could be that #1 is a prefix of #3, and #2 is just of
    //  the same size.  Even there, we'd need to examine all pairs.

    end = beg + 1;

    //  First, find a pair of adjacent matches
    //
    while ((end < fhLen) &&
           (fh[beg].hash != fh[end].hash)) {
      beg++;
      end++;
    }

    //  Got a match?
    //
    if (end < fhLen) {
      uint32 b, e;

      //  Advance end to the end of the matches
      //
      while ((fh[beg].hash == fh[end].hash) && (end < fhLen))
        end++;

      //  Yeah, we could extend scope of this test to include the for
      //  loops, but those will stop quick enough.

      if ((beg + 1 < end) && (end-beg > 1000))
        fprintf(stderr, "Large potential duplicate set from "F_U32" to "F_U32" ("F_U32" things)\n", beg, end, end - beg);

      //  Load the fragments
      //
      if (end - beg > fragsMax) {
        fragsMax = end - beg + 1;
        frags    = (fragRecord *)safe_realloc(frags, sizeof(fragRecord) * fragsMax);
      }

      for (b=beg; b<end; b++)
        getFrag(gkpStore, fh[b].iid, &frags[b-beg], FRAG_S_SEQ);

      //  Compare all-vs-all in the range
      //
      for (b=beg; b<end; b++) {
        for (e=b+1; e<end; e++) {

          AS_IID     iid1 = fh[b].iid;
          AS_IID     iid2 = fh[e].iid;

          fragRecord  *fr1 = frags + b - beg;
          fragRecord  *fr2 = frags + e - beg;

          assert(iid1 == getFragRecordIID(fr1));
          assert(iid2 == getFragRecordIID(fr2));

          uint32 del1 = getFragRecordIsDeleted(fr1);
          uint32 del2 = getFragRecordIsDeleted(fr2);

          uint32 len1 = getFragRecordSequenceLength(fr1);
          uint32 len2 = getFragRecordSequenceLength(fr2);

          if ((del1) && (len1 < len2))
            continue;
          if ((del2) && (len2 < len1))
            continue;

          if (len1 == len2) {
            if ((del1) && (iid1 < iid2))
              continue;
            if ((del2) && (iid2 < iid1))
              continue;
          }

          if (del1 && del2)
            continue;

          char *seq1 = getFragRecordSequence(fr1);
          char *seq2 = getFragRecordSequence(fr2);

          uint32 len = MIN(len1, len2);

          if (strncmp(seq1, seq2, len) == 0) {

            //  A real collision.  Delete smaller of the two (either
            //  smaller sequence length or smaller iid).  We can skip
            //  the delete if it's already deleted.

            AS_UID     deletedUID = AS_UID_undefined();
            AS_IID     deletedIID = 0;
            uint32     deleted    = 0;

            if ((len == getFragRecordSequenceLength(fr1)) &&
                (len == getFragRecordSequenceLength(fr2))) {
              deletedIID = (iid1 < iid2) ? iid1 : iid2;
              deletedUID = (iid1 < iid2) ? getFragRecordUID(fr1) : getFragRecordUID(fr2);
              deleted    = (iid1 < iid2) ? del1 : del2;
            } else if (len == getFragRecordSequenceLength(fr1)) {
              deletedIID = iid1;
              deletedUID = getFragRecordUID(fr1);
              deleted    = del1;
            } else {
              deletedIID = iid2;
              deletedUID = getFragRecordUID(fr2);
              deleted    = del2;
            }

            //  If we need to delete something, delete it, then update
            //  our cached copy.  We still need the sequence, as an
            //  even shorter fragment can be deleted by the one we
            //  just deleted.

            if (deleted == 0) {
              st.numDeDuped++;

              if (logFile)
                fprintf(logFile, "Delete read %s,%d a prefix of %s,%d\n",
                        AS_UID_toString(deletedUID), deletedIID,
                        (deletedIID == iid1) ? AS_UID_toString(getFragRecordUID(fr2)) : AS_UID_toString(getFragRecordUID(fr1)),
                        (deletedIID == iid1) ? iid2 : iid1);

              delFrag(gkpStore, deletedIID);
              getFrag(gkpStore, deletedIID, (deletedIID == iid1) ? fr1 : fr2, FRAG_S_SEQ);
            }
          }
        }
      }
    }

    beg = end;
  }

  safe_free(fh);

  fprintf(stderr, "removeDuplicateReads()-- finished\n");
}


////////////////////////////////////////////////////////////////////////////////
//
//  For a given fragRecord, scan the sequence for a linker.  If found,
//  generate two new mated reads and delete the original read.
//
//

typedef struct {
  char     h_alignA[AS_READ_MAX_LEN + AS_READ_MAX_LEN + 2];
  char     h_alignB[AS_READ_MAX_LEN + AS_READ_MAX_LEN + 2];
  dpCell   h_matrix[AS_READ_MAX_LEN + 1][AS_READ_MAX_LEN + 1];
} dpMatrix;

dpMatrix  *globalMatrix = NULL;

static
int
processMate(fragRecord *fr,
            fragRecord *m1,
            fragRecord *m2,
            char       *linker[AS_LINKER_MAX_SEQS],
            int         search[AS_LINKER_MAX_SEQS]) {

   alignLinker_s  al = {0};

   //  Did we find enough of the linker to do something?  We just need
   //  to throw out the obviously bad stuff.  When we get shorter and
   //  shorter, it's hard to define reasonable cutoffs.
   //
   //  Things that are called good here, but are actually bad, will be
   //  examined in OBT's chimera.  If there are no overlaps spanning,
   //  they'll be trimmed out, usually by being called chimeric.
   //
   int  goodAlignment = 0;
   int  bestAlignment = 0;
   int  linkerID      = 0;

   int  lSize = 0;
   int  rSize = 0;

   while (linkerID < AS_LINKER_MAX_SEQS && goodAlignment == 0) {
      if (search[linkerID] == TRUE) {
         if (linker[linkerID] == NULL) {
            fprintf(stderr, "error requested to search using linker in position %d which is not initialized\n", linkerID);
            assert(0);
         }

         alignLinker(globalMatrix->h_alignA,
                     globalMatrix->h_alignB,
                     linker[linkerID],
                     fr->seq,
                     globalMatrix->h_matrix,
                     &al,
                     FALSE, FALSE, 0, 0);

         lSize = al.begJ;
         rSize = al.lenB - al.endJ;
   
         if ((lSize < 0) || (lSize > AS_FRAG_MAX_LEN) || (rSize < 0) || (rSize > AS_FRAG_MAX_LEN)) {
            fprintf(stderr, "fragment %d root=%s\n", getLastElemStore(gkpStore->frg) + 1, AS_UID_toString(fr->gkfr.readUID));
            fprintf(stderr, "alignLen=%d matches=%d lSize=%d rSize=%d\n", al.alignLen, al.matches, lSize, rSize);
            fprintf(stderr, "I %3d-%3d %s\n", al.begI, al.endI, globalMatrix->h_alignA);
            fprintf(stderr, "J %3d-%3d %s\n", al.begJ, al.endJ, globalMatrix->h_alignB);
         }
   
         assert(lSize >= 0);
         assert(lSize <= AS_FRAG_MAX_LEN);
         assert(rSize >= 0);
         assert(rSize <= AS_FRAG_MAX_LEN);
   
         if ((al.alignLen >=  5) && (al.matches + 1 >= al.alignLen))
            goodAlignment = 1;
         if ((al.alignLen >= 15) && (al.matches + 2 >= al.alignLen))
            goodAlignment = 1;
         if ((al.alignLen >= 30) && (al.matches + 3 >= al.alignLen))
            goodAlignment = 1;
         if ((al.alignLen >= 40) && (al.matches + 4 >= al.alignLen))
            goodAlignment = 1;
      }
      
      linkerID++;
  }

   // we went through all our possible linker sequences and couldn't find a match, give up
   if (goodAlignment == 0) {
      if ((logFile) && (m1 != NULL) && (m2 != NULL))
         //  m1, m2 are NULL if we are checking a second time for linker;
         //  do not report no linker found for those as we have already
         //  or will soon write a log mesasge.
         fprintf(logFile, "No linker detected in %s.\n", AS_UID_toString(fr->gkfr.readUID));
      return(0);
   }

   if ((al.alignLen >= 42) && (al.matches + 2 >= al.alignLen))
      bestAlignment = 1;

   int  which;


  //  Not enough on either side of the read to make a mate. Chuck the read, return that we didn't trim
  if ((bestAlignment) && (lSize < 64) && (rSize < 64)) {
    if (logFile)
      fprintf(logFile, "Linker too close to left side and right side (position %d-%d); no mate formed for %s (len %d).\n",
              al.begJ, al.endJ, AS_UID_toString(fr->gkfr.readUID), fr->gkfr.seqLen);
    
    fr->gkfr.readUID = AS_UID_undefined();
    fr->gkfr.deleted = 1;
    return(1);
  }
  
  //  Adapter found on the left, but not enough to make a read.  Trim it out.
  //
  if ((bestAlignment) && (lSize < 64)) {
    if (logFile)
      fprintf(logFile, "Linker to close to left side (position %d-%d); no mate formed for %s (len %d).\n",
              al.begJ, al.endJ, AS_UID_toString(fr->gkfr.readUID), fr->gkfr.seqLen);

    fr->gkfr.seqLen = rSize;
 
    memmove(fr->seq, fr->seq + al.endJ, rSize);
    memmove(fr->qlt, fr->qlt + al.endJ, rSize);

    fr->seq[rSize] = 0;
    fr->qlt[rSize] = 0;

    for (which=0; which <= AS_READ_CLEAR_LATEST; which++) {
      if (fr->gkfr.clearBeg[which] < fr->gkfr.clearEnd[which]) {
        if (fr->gkfr.clearBeg[which] < al.endJ)
          fr->gkfr.clearBeg[which] = 0;
        else
          fr->gkfr.clearBeg[which] -= al.endJ;

        if (fr->gkfr.clearEnd[which] < al.endJ)
          fr->gkfr.clearEnd[which] = 0;
        else
          fr->gkfr.clearEnd[which] -= al.endJ;

        assert(fr->gkfr.clearEnd[which] <= rSize);
      }
    }

    //  Recursively search for another copy of the linker.
    processMate(fr, NULL, NULL, linker, search);

    return(1);
  }


  //  Adapter found on the right, but not enough to make a read.  Trim it out.
  //
  if ((bestAlignment) && (rSize < 64)) {
    if (logFile)
      fprintf(logFile, "Linker to close to right side (position %d-%d); no mate formed for %s (len %d).\n",
              al.begJ, al.endJ, AS_UID_toString(fr->gkfr.readUID), fr->gkfr.seqLen);

    fr->gkfr.seqLen = lSize;
    fr->seq[lSize] = 0;
    fr->qlt[lSize] = 0;

    for (which=0; which <= AS_READ_CLEAR_LATEST; which++) {
      if (fr->gkfr.clearEnd[which] > lSize)
        fr->gkfr.clearEnd[which] = lSize;
    }

    //  Recursively search for another copy of the linker.
    processMate(fr, NULL, NULL, linker, search);

    return(1);
  }


  //  Adapter found in the middle, and enough to make two mated reads.
  //
  if ((bestAlignment) && (lSize >= 64) && (rSize >= 64)) {

    //  If we get here, and the two mate reads are null, we have
    //  found a second complete linker -- the original read is
    //  "seq-link-seq-link-seq" and we don't know where to break.
    //  The whole read is deleted in this case.
    //
    //  The mate is deleted later.
    //
    if ((m1 == NULL) || (m2 == NULL)) {
      fr->gkfr.readUID = AS_UID_undefined();
      fr->gkfr.deleted = 1;
      return(0);
    }

    //  0.  Copy the fragments to new mated fragments
    {
      memcpy(m1, fr, sizeof(fragRecord));
      memcpy(m2, fr, sizeof(fragRecord));
    }

    //  1.  Make new UIDs for the two mated reads.  Nuke the old
    //  read.  Make the mates.
    //
    //  WARNING!  See those getLastElemStore() below?  It forces us to
    //  load the gkm1 read before the gkm2 read.
    {
      char  uid[64];
      int   len;

      strcpy(uid, AS_UID_toString(fr->gkfr.readUID));
      len = strlen(uid);
      uid[len+1] = 0;

      uid[len] = 'a';
      m1->gkfr.readUID = AS_UID_load(uid);
      m1->gkfr.readIID = getLastElemStore(gkpStore->frg) + 1;

      uid[len] = 'b';
      m2->gkfr.readUID = AS_UID_load(uid);
      m2->gkfr.readIID = getLastElemStore(gkpStore->frg) + 2;

      fr->gkfr.readUID = AS_UID_undefined();
      fr->gkfr.deleted = 1;

      m1->gkfr.mateIID = m2->gkfr.readIID;
      m2->gkfr.mateIID = m1->gkfr.readIID;

      m1->gkfr.orientation = AS_READ_ORIENT_INNIE;
      m2->gkfr.orientation = AS_READ_ORIENT_INNIE;
    }

    //  2.  Propagate clear ranges.  Math.
    //
    //  m1 is reverse complemented, so the start of m1 is next to the
    //  linker, and the end can extend into low quality sequence at
    //  the start of the read.
    //
    //  lSize - size of the left half of the read, including X's
    //  rSize - size of the right half of the read, including X's
    //
    //       v clearBeg                   clearEnd v
    //  XXXXXX-------------------[linker]----------XXXXXXXXXXX
    //                   al.begJ ^      ^ al.endJ
    //
    //  
    {
      for (which=0; which <= AS_READ_CLEAR_LATEST; which++) {
        if (fr->gkfr.clearBeg[which] < fr->gkfr.clearEnd[which]) {
          uint  lend = al.begJ - fr->gkfr.clearBeg[which];
          uint  rend = fr->gkfr.clearEnd[which] - al.endJ;

          if (lSize < fr->gkfr.clearBeg[which])
            lend = 0;

          if (fr->gkfr.clearEnd[which] < al.endJ)
            rend = 0;

          m1->gkfr.clearBeg[which] = 0;
          m1->gkfr.clearEnd[which] = lend;

          m2->gkfr.clearBeg[which] = 0;
          m2->gkfr.clearEnd[which] = rend;
        }
      }
    }

    //  3.  Construct new rm1, rm2.  Nuke the linker.  Reverse
    //  complement -- inplace -- the left mate.
    {
      int j;

      reverseComplement(m1->seq, m1->qlt, lSize);

      m1->gkfr.seqLen = lSize;

      m1->seq[lSize] = 0;
      m1->qlt[lSize] = 0;

      m2->gkfr.seqLen    = rSize;

      memmove(m2->seq, m2->seq + al.endJ, rSize);
      memmove(m2->qlt, m2->qlt + al.endJ, rSize);

      m2->seq[rSize] = 0;
      m2->qlt[rSize] = 0;

      if (logFile)
        fprintf(logFile, "Mates '%s' (%d-%d) and '%s' (%d-%d) created.\n",
                AS_UID_toString(m1->gkfr.readUID), 0, al.begJ,
                AS_UID_toString(m2->gkfr.readUID), al.endJ, al.lenB);
    }

    //  Recursively search for another copy of the linker.
    processMate(m1, NULL, NULL, linker, search);
    processMate(m2, NULL, NULL, linker, search);

    //  If either is deleted now, then we found another linker in
    //  the split read.
    //
    if ((m1->gkfr.deleted) || (m2->gkfr.deleted)) {
      if (logFile)
        fprintf(logFile, "Multiple linker detected in mates '%s' and '%s'.  Deleted.\n", AS_UID_toString(m1->gkfr.readUID), AS_UID_toString(m2->gkfr.readUID));
      m1->gkfr.readUID = AS_UID_undefined();
      m2->gkfr.readUID = AS_UID_undefined();
      m1->gkfr.deleted = 1;
      m2->gkfr.deleted = 1;
    }

    return(1);
  }  //  if match is in middle


  //  Significant alignment, but we didn't do anything with it.  Why?
  //  Overload some of the later clear ranges to convey this
  //  information to downstream processes.
  //
  //  We have 64 bits split into 4 16 bit words.  We want to store the
  //  position of the match in both the linker and the read, as well
  //  as length and number of matches.  Six things.  If we assume the
  //  linker is 255bp or smaller, we can pack.
  //
  //  This could benefit from having 64-bits of generic unioned data
  //  in the fragment record.

  if (logFile)
    fprintf(logFile, "Linker detected but not trimmed in '%s' linker: %d-%d read: %d-%d alignLen: %d matches: %d.  Passed to OBT.\n",
            AS_UID_toString(fr->gkfr.readUID), al.begI, al.endI, al.begJ, al.endJ, al.alignLen, al.matches);

  fr->gkfr.contaminationBeg = al.begJ;
  fr->gkfr.contaminationEnd = al.endJ;

  return(0);
}



int
detectMates(char *linker[AS_LINKER_MAX_SEQS], int search[AS_LINKER_MAX_SEQS]) {
  fragRecord    fr        = {0};
  fragRecord    m1        = {0};
  fragRecord    m2        = {0};

  uint32        firstElem = getFirstElemFragStore(gkpStore);
  uint32        lastElem  = getLastElemFragStore(gkpStore) + 1;
  uint32        thisElem  = 0;

  globalMatrix = safe_malloc(sizeof(dpMatrix));

  fprintf(stderr, "detectMates()-- from %d to %d\n", firstElem, lastElem);

  for (thisElem=firstElem; thisElem<lastElem; thisElem++) {
    if ((thisElem % 1000000) == 0)
      fprintf(stderr, "removeDuplicateReads()--  at %d\n", thisElem);

    getFrag(gkpStore, thisElem, &fr, FRAG_S_ALL);

    if (fr.gkfr.deleted)
      continue;
    
    m1.gkfr.readUID = AS_UID_undefined();
    m2.gkfr.readUID = AS_UID_undefined();

    //  If processMate returns true, something changed.  Delete the
    //  original read, and add either the new single read, or the two
    //  mates.
    //
    //  WARNING!  The mates MUST be added in this order, otherwise,
    //  the UID<->IID mapping will be invalid.

    if (processMate(&fr, &m1, &m2, linker, search)) {
      delFrag(gkpStore, thisElem);

      if (AS_UID_isDefined(fr.gkfr.readUID)) {
        st.numReadsWithPartialLinker++;

        assert(AS_UID_isDefined(m1.gkfr.readUID) == 0);
        assert(AS_UID_isDefined(m2.gkfr.readUID) == 0);
        addReadToStore(gkpStore, &fr);
      }

      if (AS_UID_isDefined(m1.gkfr.readUID) &&
          AS_UID_isDefined(m2.gkfr.readUID)) {
        st.numReadsWithFullLinker++;

        assert(AS_UID_isDefined(fr.gkfr.readUID) == 0);

        addReadToStore(gkpStore, &m1);
        addReadToStore(gkpStore, &m2);
      }
    } else {
      st.numReadsWithNoLinker++;
    }
  }

  safe_free(globalMatrix);
  globalMatrix = NULL;

  return(0);
}






void
addLibrary(char *libraryName,
           int   insertSize,
           int   insertStdDev,
           int   haveLinker)  {
  GateKeeperLibraryRecord   gkl;

  clearGateKeeperLibraryRecord(&gkl);

  gkl.libraryUID = AS_UID_load(libraryName);

  //  This crud is documented in AS_PER/AS_PER_gkpStore.h
  //  Zero is the default, we set to make it explicit

  gkl.spare2                      = 0;
  gkl.spare1                      = 0;

  gkl.forceBOGunitigger           = 1;

  gkl.doNotQVTrim                 = 1;
  gkl.goodBadQVThreshold          = 1;  //  Effectively, don't QV trim, redundant

  gkl.unused1                     = 0;
  gkl.doNotTrustHomopolymerRuns   = 1;
  gkl.doNotOverlapTrim            = 0;
  gkl.isNotRandom                 = 0;

  gkl.hpsIsSomethingElse          = 0;
  gkl.hpsIsFlowGram               = 1;
  gkl.hpsIsPeakSpacing            = 0;

  if (haveLinker == FALSE) {
    gkl.orientation = AS_READ_ORIENT_UNKNOWN;
    gkl.mean        = 0;
    gkl.stddev      = 0;
  } else {
    gkl.orientation = AS_READ_ORIENT_INNIE;
    gkl.mean        = insertSize;
    gkl.stddev      = insertStdDev;
  }

  appendIndexStore(gkpStore->lib, &gkl);
  setGatekeeperUIDtoIID(gkpStore, gkl.libraryUID, getLastElemStore(gkpStore->lib), AS_IID_LIB);

  assert(getLastElemStore(gkpStore->lib) == 1);
}


//  This is an efficient version of dumpGateKeeperAsFRG() in AS_GKP_dump.c
void
dumpFragFile(char *outName, FILE *outFile) {
  fragRecord        fr;

  unsigned int      firstElem = 0;
  unsigned int      lastElem = 0;

  GenericMesg       pmesg;

  LibraryMesg       libMesg;
  FragMesg          frgMesg;
  LinkMesg          lnkMesg;

  int               i;

  AS_UID           *frgUID = safe_calloc(getLastElemStore(gkpStore->frg) + 1, sizeof(AS_UID));

  char              source[256];

  //  Dump the format message
  //
  {
    VersionMesg  vmesg;

    AS_MSG_setFormatVersion(2);

    vmesg.version = 2;

    pmesg.m = &vmesg;
    pmesg.t = MESG_VER;

    WriteProtoMesg_AS(outFile, &pmesg);
  }

  //  Exactly one library here.

  {
    GateKeeperLibraryRecord  *gkpl = getGateKeeperLibrary(gkpStore, 1);

    frgUID[0] = gkpl->libraryUID;

    pmesg.m = &libMesg;
    pmesg.t = MESG_LIB;

    libMesg.action       = AS_ADD;
    libMesg.eaccession   = gkpl->libraryUID;
    libMesg.mean         = gkpl->mean;
    libMesg.stddev       = gkpl->stddev;
    libMesg.source       = gkpl->comment;
    libMesg.link_orient  = AS_READ_ORIENT_NAMES[gkpl->orientation][0];

    AS_PER_encodeLibraryFeatures(gkpl, &libMesg);

    WriteProtoMesg_AS(outFile, &pmesg);

    AS_PER_encodeLibraryFeaturesCleanup(&libMesg);
  }

  //  Dump fragments -- as soon as both reads in a mate are defined,
  //  we dump the mate relationship.

  FragStream *fs = openFragStream(gkpStore, FRAG_S_ALL);

  while (nextFragStream(fs, &fr)) {
    if (getFragRecordIsDeleted(&fr))
      continue;

    frgUID[getFragRecordIID(&fr)] = getFragRecordUID(&fr);

    pmesg.m = &frgMesg;
    pmesg.t = MESG_FRG;

    //  This code used in AS_GKP_dump.c (dumpFRG).
    frgMesg.action            = getFragRecordIsDeleted(&fr) ? AS_DELETE : AS_ADD;
    frgMesg.eaccession        = getFragRecordUID(&fr);
    frgMesg.library_uid       = frgUID[0];
    frgMesg.library_iid       = getFragRecordLibraryIID(&fr);
    frgMesg.plate_uid         = fr.gkfr.plateUID;
    frgMesg.plate_location    = fr.gkfr.plateLocation;
    frgMesg.type              = AS_READ;
    frgMesg.is_random         = (getFragRecordIsNonRandom(&fr)) ? 0 : 1;
    frgMesg.status_code       = AS_READ_STATUS_NAMES[fr.gkfr.status][0];
    frgMesg.clear_rng.bgn     = getFragRecordClearRegionBegin(&fr, AS_READ_CLEAR_LATEST);
    frgMesg.clear_rng.end     = getFragRecordClearRegionEnd  (&fr, AS_READ_CLEAR_LATEST);
    frgMesg.clear_vec.bgn     = getFragRecordClearRegionBegin(&fr, AS_READ_CLEAR_VEC);
    frgMesg.clear_vec.end     = getFragRecordClearRegionEnd  (&fr, AS_READ_CLEAR_VEC);
    frgMesg.clear_max.bgn     = getFragRecordClearRegionBegin(&fr, AS_READ_CLEAR_MAX);
    frgMesg.clear_max.end     = getFragRecordClearRegionEnd  (&fr, AS_READ_CLEAR_MAX);
    frgMesg.contamination.bgn = fr.gkfr.contaminationBeg;
    frgMesg.contamination.end = fr.gkfr.contaminationEnd;
    frgMesg.source            = 0L;
    frgMesg.sequence          = getFragRecordSequence(&fr);
    frgMesg.quality           = getFragRecordQuality(&fr);
    frgMesg.hps               = getFragRecordHPS(&fr);
    frgMesg.iaccession        = getFragRecordIID(&fr);

    WriteProtoMesg_AS(outFile, &pmesg);

    if ((getFragRecordMateIID(&fr) > 0) &&
        (getFragRecordMateIID(&fr) < getFragRecordIID(&fr))) {
      pmesg.m = &lnkMesg;
      pmesg.t = MESG_LKG;

      lnkMesg.action      = AS_ADD;
      lnkMesg.type        = AS_MATE;
      lnkMesg.link_orient = AS_READ_ORIENT_NAMES[fr.gkfr.orientation][0];
      lnkMesg.frag1       = frgUID[getFragRecordMateIID(&fr)];
      lnkMesg.frag2       = getFragRecordUID(&fr);
      lnkMesg.distance    = frgUID[0];

      WriteProtoMesg_AS(outFile, &pmesg);
    }
  }

  safe_free(frgUID);

  closeFragStream(fs);
}


int
main(int argc, char **argv) {
  int       insertSize       = 0;
  int       insertStdDev     = 0;
  char     *libraryName      = 0L;
  FILE     *outputFile       = 0L;
  int       firstFileArg     = 0;
  char      outputName[FILENAME_MAX]   = {0};
  char      gkpStoreName[FILENAME_MAX] = {0};
  char     *logFileName      = 0L;
  char     *statsFileName    = 0L;

  bool      doDeDup          = 1;

  // initialize linker search structure
  // One array stores the character sequences of the linker
  // A boolean array stores which linkers are to be used in the search

  int       haveLinker        = FALSE;
  int       invalidLinkerSeq  = FALSE;
  char     *linker[AS_LINKER_MAX_SEQS] = { NULL };
  int       search[AS_LINKER_MAX_SEQS] = { 0    };
  
  // the first slot of the linker array is the FLX mate pair linker (which is a palindrome)
  char   *linkerFLX     = linker[0] = "GTTGGAACCGAAAGGGTTTGAATTCAAACCCTTTCGGTTCCAAC";  // palindrome
  // the next two slots are the Titanium linker. It requires two linkers because they are not palindromes
  char   *linkerFIX     = linker[1] = "TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG";    // linker for Titanium reads
  char   *linkerXIF     = linker[2] = "CGTAATAACTTCGTATAGCATACATTATACGAAGTTATACGA";    // rc of linker for Titanium reads
  // subsequent linkers will be used for future barcoding
  // final linkers are custom, provided by the user, filled in when parsing parameters

  int       bogusOptions[256] = {0};
  int       bogusOptionsLen   = 0;

  int arg = 1;
  int err = 0;

  while (arg < argc) {
    if        (strcmp(argv[arg], "-insertsize") == 0) {
      insertSize   = atoi(argv[++arg]);
      insertStdDev = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-libraryname") == 0) {
      libraryName = argv[++arg];

    } else if (strcmp(argv[arg], "-clear") == 0) {
      arg++;

      //  If this is the first time we get a -clear switch, set
      //  clearAction to exactly that value.  Later times through,
      //  we'll add in more options.
      //
      if      (strcasecmp(argv[arg], "all") == 0)
        clearAction = (clearSet == 0) ? (CLEAR_ALL) : (clearAction | CLEAR_ALL);
      else if (strcasecmp(argv[arg], "454") == 0)
        clearAction = (clearSet == 0) ? (CLEAR_454) : (clearAction | CLEAR_454);
      else if (strcasecmp(argv[arg], "n") == 0)
        clearAction = (clearSet == 0) ? (CLEAR_N) : (clearAction | CLEAR_N);
      else if (strcasecmp(argv[arg], "pair-of-n") == 0)
        clearAction = (clearSet == 0) ? (CLEAR_PAIR_N) : (clearAction | CLEAR_PAIR_N);
      else if (strcasecmp(argv[arg], "discard-n") == 0)
        clearAction = (clearSet == 0) ? (CLEAR_DISCARD_N) : (clearAction | CLEAR_DISCARD_N);
      else {
        clearAction = (clearSet == 0) ? (CLEAR_ERRR) : (clearAction | CLEAR_ERRR);
        err++;
      }

    } else if (strcmp(argv[arg], "-trim") == 0) {
      arg++;

      if      (strcasecmp(argv[arg], "none") == 0)
        trimAction = TRIM_NONE;
      else if (strcasecmp(argv[arg], "soft") == 0)
        trimAction = TRIM_SOFT;
      else if (strcasecmp(argv[arg], "hard") == 0)
        trimAction = TRIM_HARD;
      else if (strcasecmp(argv[arg], "chop") == 0)
        trimAction = TRIM_CHOP;
      else {
        trimAction = TRIM_ERRR;
        err++;
      }

    } else if (strcmp(argv[arg], "-linker") == 0) {
      arg++;

      if      (strcasecmp(argv[arg], "flx") == 0) {
        search[0]    = TRUE;
        haveLinker   = TRUE;
      }
      else if (strcasecmp(argv[arg], "titanium") == 0) {
        search[1]    = TRUE;
        search[2]    = TRUE;
        haveLinker   = TRUE;
      }
      else {
        int start = AS_LINKER_CUSTOM_OFFSET;
        if (AS_UTL_isValidSequence(argv[arg], strlen(argv[arg]))) {
          linker[start]      = argv[arg];
          search[start++]    = TRUE;
          haveLinker         = TRUE;
        } else {
          invalidLinkerSeq     = TRUE;
          err++;
        }
      }

    } else if (strcmp(argv[arg], "-nodedup") == 0) {
      doDeDup = 0;

    } else if (strcmp(argv[arg], "-output") == 0) {
      strcpy(outputName, argv[++arg]);

    } else if (strcmp(argv[arg], "-log") == 0) {
      logFileName = argv[++arg];

    } else if (strcmp(argv[arg], "-stats") == 0) {
      statsFileName = argv[++arg];

    } else {
      if (argv[arg][0] == '-') {
        bogusOptions[bogusOptionsLen++] = arg;
        err++;
      } else {
        firstFileArg = arg;
        arg          = argc;
      }
    }

    arg++;
  }

  //  Have a linker but no insert size?  Error.
  if ((haveLinker) && ((insertSize == 0) || (insertStdDev == 0)))
    err++;

  //  Have an insert size but no linker?  Error.
  if ((!haveLinker) && ((insertSize != 0) || (insertStdDev != 0)))
    err++;

  if ((err) || (libraryName == 0L) || (outputName[0] == 0) || (firstFileArg == 0)) {
    fprintf(stderr, "usage: %s [opts] -libraryname n -output f.frg in.sff ...\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -insertsize i d        Mates are on average i +- d bp apart.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -libraryname n         The UID of the library these reads are added to.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -clear all             Use the whole read.\n");
    fprintf(stderr, "  -clear 454             Use the 454 clear ranges as is (default).\n");
    fprintf(stderr, "  -clear n               Use the whole read up to the first N.\n");
    fprintf(stderr, "  -clear pair-of-n       Use the whole read up to the frist pair of Ns.\n");
    fprintf(stderr, "  -clear discard-n       Delete the read if there is an N in the clear range.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  If multiple -clear options are supplied, the intersection is used.  For\n");
    fprintf(stderr, "  'discard-n', the clear range is first computed, then if there is still an\n");
    fprintf(stderr, "  N in the clear range, the read is deleted.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Caution!  Even though the default is '454', when any -clear option is used,\n");
    fprintf(stderr, "  the list of clear ranges to intersect is reset.  To get both '454' and 'n',\n");
    fprintf(stderr, "  BOTH '-clear 454' and '-clear n' must be supplied on the command line.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -trim none             Use the whole read regardless of -clear settings.\n");
    fprintf(stderr, "  -trim soft             OBT and ECR can increase the clear range.\n");
    fprintf(stderr, "  -trim hard             OBT can only shrink the clear range, but ECR can extend (default).\n");
    fprintf(stderr, "  -trim chop             Erase sequence outside the clear range.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  'none' will emit the whole read, and reset clear ranges to cover the whole read.\n");
    fprintf(stderr, "  'soft' will emit the whole read, and leave clear ranges as set.\n");
    fprintf(stderr, "  'hard' is like soft, with the addition of a 'clm' message to stop OBT.\n");
    fprintf(stderr, "  'chop' is like none, but after the read is chopped down to just the clear bases.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -linker [name | seq]   Search for linker, create mated reads.\n");
    fprintf(stderr, "                         Name is one of:\n");
    fprintf(stderr, "                           'flx'      == %s\n",     linkerFLX);
    fprintf(stderr, "                           'titanium' == %s and\n", linkerFIX);
    fprintf(stderr, "                                         %s\n",     linkerXIF);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -nodedup               Do not remove reads that are a perfect prefix of another read.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -output f.frg          Write the CA formatted fragments to this file.\n");
    fprintf(stderr, "  -log    l.txt          Human readable log of what happened.\n");
    fprintf(stderr, "  -stats  s.txt          Human readable statistics; summarizes the log file.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "See http://apps.sourceforge.net/mediawiki/wgs-assembler/index.php?title=Formatting_Inputs\n");
    fprintf(stderr, "\n");

    for (err=0; err<bogusOptionsLen; err++)
      fprintf(stderr, "ERROR:  Unknown option '%s'\n", argv[bogusOptions[err]]);

    if (libraryName == 0L)
      fprintf(stderr, "ERROR:  Need to supply -libraryname.\n");

    if (outputName[0] == 0)
      fprintf(stderr, "ERROR:  Need to supply -output.\n");

    if (firstFileArg == 0)
      fprintf(stderr, "ERROR:  Need to supply some SFF files.\n");

    if ((haveLinker) && ((insertSize == 0) ||
                         (insertStdDev == 0)))
      fprintf(stderr, "ERROR:  Have a linker sequence, but no insert size set with -insertsize.\n");

    if ((!haveLinker) && ((insertSize != 0) ||
                          (insertStdDev != 0)))
      fprintf(stderr, "ERROR:  Have an insert size, bu no linker sequence set with -linker.\n");

    if (clearAction == CLEAR_ERRR)
      fprintf(stderr, "ERROR:  Unknown -clear value.\n");

    if (trimAction == TRIM_ERRR)
      fprintf(stderr, "ERROR:  Unknown -trim value.\n");

    if (invalidLinkerSeq == TRUE)
      fprintf(stderr, "ERROR:  Invalid -linker value. It must be one of titanium, flx, or a valid ACGT string.\n");
    
    exit(1);
  }

  logFile = NULL;

  if (logFileName) {
    errno = 0;
    logFile = fopen(logFileName, "w");
    if (errno)
      fprintf(stderr, "ERROR: Failed to open the log file '%s': %s\n", logFileName, strerror(errno)), exit(1);
  }

  if ((outputName[strlen(outputName)-4] != '.') || 
      (outputName[strlen(outputName)-3] != 'f') || 
      (outputName[strlen(outputName)-2] != 'r') || 
      (outputName[strlen(outputName)-1] != 'g'))
    strcat(outputName, ".frg");

  if (fileExists(outputName, FALSE, FALSE))
    fprintf(stderr, "ERROR: Output file '%s' exists; I will not clobber it.\n", outputName), exit(1);

  errno = 0;
  outputFile = fopen(outputName, "w");
  if (errno)
    fprintf(stderr, "ERROR: Failed to open the output file '%s': %s\n", outputName, strerror(errno)), exit(1);

  strcpy(gkpStoreName, outputName);
  strcat(gkpStoreName, ".tmpStore");

  if (fileExists(gkpStoreName, TRUE, FALSE)) {
    fprintf(stderr, "ERROR: Temporary Gatekeeper Store '%s' exists; I will not clobber it.\n", gkpStoreName);
    fprintf(stderr, "       If this is NOT from another currently running sffToCA, simply remove this directory.\n");
    exit(1);
  }

  gkpStore      = createGateKeeperStore(gkpStoreName);
  gkpStore->frg = convertStoreToMemoryStore(gkpStore->frg);

  addLibrary(libraryName, insertSize, insertStdDev, haveLinker);

  for (; firstFileArg < argc; firstFileArg++)
    loadSFF(argv[firstFileArg]);

  if (doDeDup)
    removeDuplicateReads();

  if (haveLinker)
    detectMates(linker, search);

  dumpFragFile(outputName, outputFile);

  closeGateKeeperStore(gkpStore);
  deleteGateKeeperStore(gkpStoreName);

  errno = 0;
  fclose(outputFile);
  if (errno)
    fprintf(stderr, "Failed to close '%s': %s\n", outputName, strerror(errno)), exit(1);

  errno = 0;
  if (logFile)
    fclose(logFile);
  if (errno)
    fprintf(stderr, "Failed to close '%s': %s\n", logFileName, strerror(errno)), exit(1);

  errno = 0;
  logFile = (statsFileName) ? fopen(statsFileName, "w") : NULL;
  if (errno)
    fprintf(stderr, "ERROR: Failed to open the stats file '%s': %s\n", statsFileName);
  if (logFile) {
    writeStatistics(logFile);
    fclose(logFile);
  } else {
    writeStatistics(stderr);
  }

  return(0);
}
