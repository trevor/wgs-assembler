
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
/*********************************************************************
 * $Id: AS_CGB_testhisto.c,v 1.6 2006-11-14 19:58:21 eliv Exp $
 * Module:
 * Description:
 * Assumptions:
 *********************************************************************/

#include <stdio.h>
#include "AS_CGB_histo.h"

#define MIN(a,b) (a<b ? a : b)
#define MAX(a,b) (a>b ? a : b)

/*************************************************************************/
static char CM_ID[] = "$Id: AS_CGB_testhisto.c,v 1.6 2006-11-14 19:58:21 eliv Exp $";
/*************************************************************************/

typedef struct {
  int n;
  int x;
  int y;
  int z;
  int p;
} MyDataType;

void setdata(MyDataType *a,int i,MyDataType *b) {
  a[i] = *b;
  a[i].p = 0;
}
MyDataType *indexdata(MyDataType *a,int i) {
  return &(a[i]);
}
void aggregate(MyDataType *a,int i,MyDataType *b) {
  a[i].n += b->n;
  a[i].x += b->x;
  a[i].y = MIN(a[i].y,b->y);
  a[i].z = MAX(a[i].z,b->z);
}
int score(MyDataType *data) {
  int i;
  i = data->x / data-> n;
  return i;
}
void printdata(FILE *fout,
	       MyDataType *data,
	       MyDataType *scan_data,
	       MyDataType *aggr_data)
{
  fprintf(fout,"%10d %10.5f %10d %10d",
	  data->n,((float)data->x)/((float)data->n),data->y,data->z);
  fprintf(fout," %10d %10.5f %10d %10.5f ",
	  (scan_data->n),
	  ((float)scan_data->n)/((float)aggr_data->n),
	  (scan_data->x),
	  ((float)scan_data->x)/((float)aggr_data->x));
}


#define NSAMPLE 100
#define NBUCKET 100

main()
{ 
  HISTOGRAM *h;
  MyDataType d;
  MyDataType samples[NSAMPLE],buckets[NBUCKET];
  int i;
  int j;

  h = create_histogram(NSAMPLE,NBUCKET,0);
#ifndef NEVER
  extend_histogram(h,sizeof(MyDataType),
		   indexdata,setdata,aggregate,printdata);
#endif

#ifdef NEVER
  for(i=-105;i<103; i++) {
    j = bucket_from_score(h,i);
    printf("  %5d %5d \n",i,j);
  }
#endif

  for(i=0; i<NBUCKET; i++){
    j = (i * 127 + 10000 % 35);
    j = i*i;
    d.x = 0;
    d.y = 0;
    d.z = 0;
    buckets[i] = d;
  }
  for (i = -25; i < 25; i++) {
    int j;
    j = (i * 127 + 10000 % 35);
    j = i*i - 10;
    j = i;
    d.n = 1;
    d.x = j;
    d.y = i;
    d.z = i;
    printf("iteration %d\n",i);
    add_to_histogram(h,j,&d);
  }
  printf("here\n");
#ifdef NEVER
  print_histogram(stdout,h,NBUCKET,1);
#endif
  printf("howdy\n");
  print_histogram(stdout,h,0,1);
  free_histogram(h);
}
