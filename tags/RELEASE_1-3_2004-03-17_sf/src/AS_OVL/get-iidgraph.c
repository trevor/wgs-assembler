
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
/*************************************************
* Module:  get-subgraph.c
* Description:
*   Reads the stream of messages produced by the overlapper
*   and filters those that do not refer to fragments contained
*   int the list of internal fragment ID's specified on the
*   the command line.
*
*    Programmer:  A. Delcher
*       Written:  15 Jul 99
*  Last Revised:  
* 
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: get-iidgraph.c,v 1.1.1.1 2004-04-14 13:52:39 catmandew Exp $
 * $Revision: 1.1.1.1 $
*/

static char fileID[] = "$Id: get-iidgraph.c,v 1.1.1.1 2004-04-14 13:52:39 catmandew Exp $";

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
#include "AS_MSG_pmesg.h"
#include  "AS_OVL_delcher.h"


#define  HASH_LOAD_FACTOR     0.3333
    //  Portion of hash table to be occupied by keys.
#define  MAX_LIST_LEN             1000000
    //  Length of insert list
#define  SKIP_MODULUS         11
    //  Modulus for secondary hash function to determine skip distance
    //  for linear probing.


int32  Extra_List [MAX_LIST_LEN];
int32  List_Len = 0;

static int  Cmp
    (const void * a, const void * b)
  {
   int32  a_val, b_val;

   a_val = * ((int32 *) a);
   b_val = * ((int32 *) b);

   if  (a_val < b_val)
       return  -1;
   else if  (b_val < a_val)
       return  1;
     else
       return  0;
  }
static int  Hash_Find
    (int32 key, int32 tab [], int32 size);
static void  Insert
    (int32 key);
static int32  Next_Odd_Prime
    (int32 N);
static void  Read_Frag_IDs
    (FILE * fidfile, int32 * * hash_table, int32 * hash_table_size);



int  main  (int argc, char * argv [])

  {
   FILE  * ovlfile, * outfile, * fidfile, * odd_file, * odd_olap_file;
   char  * infile_name, * outfile_name;
   GenericMesg  * gmesg = NULL;
   MesgReader  read_msg_fn;
   MesgWriter  write_msg_fn;
   GenericMesg  * pmesg;
   AuditLine  audit_line;
   AuditMesg  * new_adt_mesg;
   int32  * hash_table, hash_table_size;
   char  label_line [1000];
   int  use_binary_output = TRUE;
   int  ch, error_flag; 

   optarg = NULL;
   error_flag = FALSE;
   while  (! error_flag && ((ch = getopt (argc, argv, "P")) != EOF))
     switch  (ch)
       {
        case  'P':
          use_binary_output = FALSE;
          break;
        case  '?':
          fprintf (stderr, "Unrecognized option \"-%c\"\n", optopt);
        default :
          error_flag ++;
        }

   if  (optind >= argc)
       {
        fprintf (stderr, 
                 "USAGE:  %s [-P] <ovl-file> <frag-id-file>\n", 
                 argv [0]);
        exit (EXIT_FAILURE);
       }

   infile_name = strdup (argv [optind ++]);
   assert (infile_name != NULL);
   fprintf (stderr, "Input Overlap File = %s\n", infile_name);
   {
     int len = strlen (infile_name);
     outfile_name = (char *) Safe_malloc (len + 20);
   }

   { // Process the input file neam to create the output file name.
     int i1 = strrchr(infile_name, '/') - infile_name;
     int i2 = strrchr(infile_name, '.') - infile_name;

     strncpy (outfile_name, &(infile_name[i1+1]), i2-i1 );
     
     strcat (outfile_name, "sub.");
     strcat (outfile_name, infile_name + i2 + 1);
     fprintf (stderr, "Output Overlap File = %s\n", outfile_name);
   }
   
   ovlfile = File_Open (infile_name, "r");
   outfile = File_Open (outfile_name, "w");
   fidfile = File_Open (argv [optind], "r");

   Read_Frag_IDs (fidfile, & hash_table, & hash_table_size);
   fprintf (stderr, "Hash table size = %d\n", hash_table_size);

   read_msg_fn = InputFileType_AS (ovlfile);
   if  (use_binary_output)
       write_msg_fn = OutputFileType_AS (AS_BINARY_OUTPUT);
     else
       write_msg_fn = OutputFileType_AS (AS_PROTO_OUTPUT);

   pmesg = (GenericMesg *) Safe_malloc (sizeof (GenericMesg));
   pmesg -> t = MESG_ADT;
   pmesg -> m = (AuditMesg *) Safe_malloc (sizeof (AuditMesg));
   new_adt_mesg = pmesg -> m;
   new_adt_mesg -> list = & audit_line;
      
   odd_olap_file = fopen ("oddity.olaps", "w");
   assert (odd_olap_file != NULL);

   while  (read_msg_fn (ovlfile, & gmesg) != EOF && gmesg != NULL)
     switch  (gmesg -> t)
       {
        case  MESG_ADT :
          {
           AuditMesg  * adt_mesg = gmesg -> m;

           sprintf (label_line, "%s %s %s", argv [0], infile_name,
                    argv [optind]);
           AppendAuditLine_AS (adt_mesg, & audit_line, time (0), "get-subgraph",
                               "$Revision: 1.1.1.1 $", label_line);
           write_msg_fn (outfile, gmesg);
           break;
          }

        case  MESG_ILK :
          {
           InternalLinkMesg  * ilk_mesg = gmesg -> m;
           int  find1, find2;
          
           find1 = Hash_Find (ilk_mesg -> ifrag1, hash_table, hash_table_size);
           find2 = Hash_Find (ilk_mesg -> ifrag2, hash_table, hash_table_size);

#if  0
           if  (find1 && ! find2)
               Insert (ilk_mesg -> ifrag2);
           if  (! find1 && find2)
               Insert (ilk_mesg -> ifrag1);
#endif
           if  (find1 && find2)
               write_msg_fn (outfile, gmesg);
           break;
          }

        case  MESG_OFG :
        case  MESG_OFR :
        case  MESG_SFG :
          {
           OFGMesg  * ofg_mesg = gmesg -> m;
          
           if  (Hash_Find (ofg_mesg -> iaccession, hash_table,
                           hash_table_size))
               write_msg_fn (outfile, gmesg);
           break;
          }

        case  MESG_OVL :
          {
           OverlapMesg  * ovl_mesg = gmesg -> m;
           int  finda, findb;

           finda = Hash_Find (ovl_mesg -> aifrag, hash_table, hash_table_size);
           findb = Hash_Find (ovl_mesg -> bifrag, hash_table, hash_table_size);

           if  (finda && ! findb)
               {
                Insert (ovl_mesg -> bifrag);
                fprintf (odd_olap_file, "%7d%c  %7d%c  %c  %4d  %4d  %4.2f\n",
                         ovl_mesg -> aifrag, ' ',
                         ovl_mesg -> bifrag, '*',
                         (char) (ovl_mesg -> orientation),
                         ovl_mesg -> ahg,
                         ovl_mesg -> bhg,
                         100.0 * ovl_mesg -> quality);
               }
           if  (! finda && findb)
               {
                Insert (ovl_mesg -> aifrag);
                fprintf (odd_olap_file, "%7d%c  %7d%c  %c  %4d  %4d  %4.2f\n",
                         ovl_mesg -> aifrag, '*',
                         ovl_mesg -> bifrag, ' ',
                         (char) (ovl_mesg -> orientation),
                         ovl_mesg -> ahg,
                         ovl_mesg -> bhg,
                         100.0 * ovl_mesg -> quality);
               }
           if  (finda && findb) /* (finda || findb) */
               write_msg_fn (outfile, gmesg);
           break;
          }

        default :
          write_msg_fn (outfile, gmesg);
       }

   fclose (odd_olap_file);


   fprintf (stderr, "List_Len = %d\n", List_Len);

   qsort ((void *) Extra_List, (size_t) List_Len, sizeof (int32), Cmp);

   odd_file = fopen ("oddity.out", "w");
   assert (odd_file != NULL);
   {
     int i;
     for  (i = 0;  i < List_Len;  i ++)
       fprintf (odd_file, "%d\n", Extra_List [i]);
   }
   fclose (odd_file);


   rewind (ovlfile);

   { 
     int sub = 0;
     while  (read_msg_fn (ovlfile, & gmesg) != EOF && gmesg != NULL)
       switch  (gmesg -> t)
         {
         case  MESG_OFG :
         case  MESG_OFR :
         case  MESG_SFG :
           {
             OFGMesg  * ofg_mesg = gmesg -> m;
             
             while  (sub < List_Len && Extra_List [sub] < ofg_mesg -> iaccession)
               sub ++;
             if  (ofg_mesg -> iaccession == Extra_List [sub])
               {
                 write_msg_fn (outfile, gmesg);
               }
             
             break;
           }
           
         default :
           break;
         }
   }
   
   fclose (ovlfile);
   fclose (outfile);

   return  0;
  }



static void  Insert
    (int32 key)

//  Insert  key  into global  Extra_List  in order if
//  it's not there already

  {
   int  i;

   for  (i = 0;  i < List_Len;  i ++)
     if  (Extra_List [i] == key)
         return;

   assert (List_Len < MAX_LIST_LEN - 1);
   Extra_List [List_Len ++] = key;

   return;
  }



static int  Hash_Find
    (int32 key, int32 tab [], int32 size)

//  Return  TRUE  iff  key  occurs in hash table  tab [0 .. (size - 1)] .

  {
   int  skip, sub;

   sub = key % size;
   skip = 1 + (key % SKIP_MODULUS);

   while  (tab [sub] != -1 && tab [sub] != key)
     sub = (sub + skip) % size;

   return  (tab [sub] == key);
  }



static int32  Next_Odd_Prime
    (int32 N)

//  Return the first odd prime  >= N .  Return  0  if can't find
//  one.

  {
   int32  Div, Last;

   if  (N % 2 == 0)
       N ++;
   while  (N < INT_MAX)
     {
      Last = (int64) (sqrt ((double) N));
      for  (Div = 3;  Div <= Last;  Div += 2)
        if  (N % Div == 0)
            break;
      if  (Div > Last)
          return  N;
      N += 2;
     }

   return  0;
  }



static void  Read_Frag_IDs
    (FILE * fp, int32 * * tab, int32 * size)

//  Read integers from  fp  and insert them into a hash table  (* tab) .
//  Set  (* size)  to the size of the hash table.  The hash function is
//  just  key % size .

  {
   int32  i, ct, key, skip, sub;

   for  (ct = 0;  fscanf (fp, "%d", & key) != EOF;  ct ++)
     ;

// fprintf (stderr, "ct = %d\n", ct);
   (* size) = Next_Odd_Prime ((int32) (ct / HASH_LOAD_FACTOR));
   assert (* size > 0);
   (* tab) = (int32 *) Safe_malloc ((* size) * sizeof (int32));

// fprintf (stderr, "(* size) = %d\n", (* size));
   for  (i = 0;  i < (* size);  i ++)
     (* tab) [i] = -1;

   rewind (fp);

   while  (fscanf (fp, "%d", & key) != EOF)
     {
      sub = key % (* size);
      skip = 1 + (key % SKIP_MODULUS);

      while  ((* tab) [sub] != -1 && (* tab) [sub] != key)
{
// fprintf (stderr, "try sub = %d\n", sub);
        sub = (sub + skip) % (* size);
}

// fprintf (stderr, "key = %d  sub = %d\n", key, sub);
      if  ((* tab) [sub] == -1)
          (* tab) [sub] = key;
     }

   return;
  }
