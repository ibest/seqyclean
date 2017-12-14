/****************************************************************************
*
* Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
* Maryland, U.S.A.  All rights reserved.
*
****************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "abi.h"


#define VERYBAD 16

static struct stack_struct {
  int i, j, m, c;
} stack[100];
static int abi_size=16;
static unsigned abi_mask;
    
struct abi_struct {
  unsigned tag;
  int index;
} *abi_well;
struct hit_struct {
  int diff, count;
} *hit_well;

static int condition[]={
-1, -1, -1, -1, -1, -1, -1, -1, -1,  0, /*  0-19 */
 0,  0,  0,  0,  0,  0,  0,  0,  0,  3,
 3,  3,  3,  3,  3,  3,  3,  3,  3,  3, /* 20-39 */
 3,  3,  3,  3,  3,  3,  3,  3,  3,  6,
 6,  6,  6,  6,  6,  6,  6,  6,  6,  6, /* 40-59 */
 6,  6,  6,  6,  6,  6,  6,  6,  6,  9,
 9,  9,  9,  9,  9,  9,  9,  9,  9,  9, /* 60-79 */
 9,  9,  9,  9,  9,  9,  9,  9,  9, 12,
12, 12, 12, 12, 12, 12, 12, 12, 12, 12, /* 80-99 */
12, 12, 12, 12, 12, 12, 12, 12, 12, 15};

static int badness[16][16]= {
/*       A   C   G   T   U   R   Y   M   W   S   K   D   H   V   B   N */
/*A*/ {  0,  3,  3,  3,  3,  1,  3,  1,  1,  3,  3,  1,  1,  1,  3,  2},
/*C*/ {  3,  0,  3,  3,  3,  3,  1,  1,  3,  1,  3,  3,  1,  1,  1,  2},
/*G*/ {  3,  3,  0,  3,  3,  1,  3,  3,  3,  1,  1,  1,  3,  1,  1,  2},
/*T*/ {  3,  3,  3,  0,  0,  3,  1,  3,  1,  3,  1,  1,  1,  3,  1,  2},
/*U*/ {  3,  3,  3,  0,  0,  3,  1,  3,  1,  3,  1,  1,  1,  3,  1,  2},
/*R*/ {  1,  3,  1,  3,  3,  1,  3,  2,  2,  2,  2,  2,  3,  2,  3,  3},
/*Y*/ {  3,  1,  3,  1,  1,  3,  1,  2,  2,  2,  2,  3,  2,  3,  2,  3},
/*M*/ {  1,  1,  3,  3,  3,  2,  2,  1,  2,  2,  3,  3,  2,  2,  3,  3},
/*W*/ {  1,  3,  3,  1,  1,  2,  2,  2,  1,  3,  2,  2,  2,  3,  3,  3},
/*S*/ {  3,  1,  1,  3,  3,  2,  2,  2,  3,  1,  2,  3,  3,  2,  2,  3},
/*K*/ {  3,  3,  1,  1,  1,  2,  2,  3,  2,  2,  1,  2,  3,  3,  2,  3},
/*D*/ {  1,  3,  1,  1,  1,  2,  3,  3,  2,  3,  2,  2,  3,  3,  3,  3},
/*H*/ {  1,  1,  3,  1,  1,  3,  2,  2,  2,  3,  3,  3,  2,  3,  3,  3},
/*V*/ {  1,  1,  1,  3,  3,  2,  3,  2,  3,  2,  3,  3,  3,  2,  3,  3},
/*B*/ {  3,  1,  1,  1,  1,  3,  2,  3,  3,  2,  2,  3,  3,  3,  2,  3},
/*N*/ {  2,  2,  2,  2,  2,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3},
/*       A   C   G   T   U   R   Y   M   W   S   K   D   H   V   B   N */
};

static int Dlist[]={ 0, 2, 3}, Hlist[]={ 0, 1, 3};

void abi_sort(abi_struct *l, abi_struct *r);

int abi_code(int c)
{
  int tmp = 0;

  switch (c) {
        case 0: tmp=0; break;
        case 1: tmp=1; break;
        case 2: tmp=2; break;
        case 3: 
        case 4: tmp=3; break;
        case 5: tmp=rand()%2 ? 2 : 0; break;
        case 6: tmp=rand()%2 ? 3 : 1; break;
        case 7: tmp=rand()%2        ; break;
        case 8: tmp=rand()%2 ? 3 : 0; break;
        case 9: tmp=rand()%2 ? 2 : 1; break;
        case 10: tmp=rand()%2 ? 3 : 2; break;
        case 11: tmp=Dlist[rand()%3]; break;
        case 12: tmp=Hlist[rand()%3]; break;
        case 13: tmp=rand()%3; break;
        case 14: tmp=rand()%3+1; break;
        case 15: tmp=rand()%4; break;
  //default:
  //  giveup("how can other cases happen in tag_code?");
  }
  return tmp;
}

void abi_sort(abi_struct *l, abi_struct *r)
//struct abi_struct *l, *r;
{
  register unsigned v;
  register struct abi_struct *i, *j, tmp;
 
  if (r>l) {
    v = r->tag; i = l-1; j = r;
    while (1) {
      while ((++i)->tag < v) ;
      while (j>l && (--j)->tag > v) ;
      if (i>=j) break;
      tmp = *i; *i = *j; *j = tmp;
    }
    tmp = *i; *i = *r; *r = tmp;
    abi_sort(l, i-1);
    abi_sort(i+1, r);
  }
}

void prepare_abi_mask()
{
  register int i;
  register unsigned acc;

  /* construct abi_mask */
  for (acc=i=0; i<abi_size; i++) {
    acc<<=2;
    acc|=3;
  }
  abi_mask=acc;
}

//void abi_align(char *a, int width, char *b, int height, int *start, int *end)
void abi_align(unsigned char *a, int width, unsigned char *b, int height, int *start, int *end)
//char *a, *b;
//int width, height;
//int *start, *end;
{
  register int i, j, l, r, m, diff, count;
  register unsigned x;
  int left, span, abi_len, hit_len;

  /* construct 'a' sequence tags */
  for (x=i=0; i<abi_size-1; i++, x<<=2)
    x|=abi_code(a[i]);
  for (j=0; i<width; i++, j++, x<<=2) {
    x|=abi_code(a[i]);
    x&=abi_mask;
    abi_well[j].tag=x;
    abi_well[j].index=i;
  }
  abi_len=j;
  abi_sort(abi_well, &abi_well[abi_len-1]);
  for (i=j=0; j<abi_len; j++)
    if (abi_well[i].tag!=abi_well[j].tag)
      abi_well[++i]=abi_well[j];
  abi_len=i+1;

  /* construct 'b' sequence tags, and search against 'a' sequence tags */
  for (x=i=0; i<abi_size-1; i++, x<<=2)
    x|=abi_code(b[i]);
  for (count=hit_len=0; i<height; i++, x<<=2) {
    x|=abi_code(b[i]);
    x&=abi_mask;
    for (l=0, r=abi_len-1; l<=r; ) {
      m=(l+r)/2;
      if (x<abi_well[m].tag)
	r=m-1;
      else if (x>abi_well[m].tag)
	l=m+1;
      else {
	l=abi_well[m].index-i;
	for (j=0; j<hit_len; j++)
	  if (hit_well[j].diff==l) {
	    hit_well[j].count++;
	    if (hit_well[j].count>count) {
	      count=hit_well[j].count;
	      diff=hit_well[j].diff;
	      if (j) {
		hit_well[j]=hit_well[0];
		hit_well[0].count=count;
		hit_well[0].diff=diff;
	      }
	    }
	    break;
	  }
	if (j>=hit_len) {
	  hit_well[j].count=1;
	  hit_well[j].diff=l;
	  hit_len=j+1;
	}
	break;
      }
    }
  }

  /* calculate starting points for sequence a and b */
  if (count<=0) {
    *start=*end=0;
    return;
  } else if (diff>=0) {
    i=diff; j=0;
  } else {
    i=0; j=-diff;
  }

  /* find the primary alignment region first */
  for (span=0, left=x=i; i<width && j<height; i++, j++)
    if (badness[a[i]][b[j]]) {
      /*
      for (count=0, m=0, l=i+1, r=j+1; 
	   l<width && r<height && count<VERYBAD && count>condition[m];  
	   m++, l++, r++)
	if (badness[a[l]][b[r]]) count+=badness[a[l]][b[r]];
      if (l<width && r<height && count<VERYBAD)
	continue;
	*/
      if (i-x>span) {
	left=x;
	span=i-x;
      }
      x=i+1;
    }
  if (i-x>span) {
    left=x;
    span=i-x;
  }
  if (span<=0) {
    *start=*end=0;
    return;
  }
  
  /* extend alignment region toward the right, if possible */
  for (i=left+span, j=i-diff; i<width && j<height; i++, j++) 
    if (badness[a[i]][b[j]]) {
      x=0;
      stack[x].i=i+1; stack[x].j=j+1; stack[x].m=0; stack[x++].c=0;
      while (x>0) {
	for (count=stack[--x].c, m=stack[x].m, l=stack[x].i, r=stack[x].j; 
	     l<width && r<height && badness[a[l]][b[r]]==0
	       && count>condition[m];  
	     m++, l++, r++) ;
	if (l>=width || r>=height)
	  continue;
	if (count<=condition[m])
	  break;
	if ((count+=badness[a[l]][b[r]])<VERYBAD) {
	  stack[x].i=l; stack[x].j=r+1; stack[x].m=m; stack[x++].c=count;
	  stack[x].i=l+1; stack[x].j=r; stack[x].m=m+1; stack[x++].c=count;
	  stack[x].i=l+1; stack[x].j=r+1; stack[x].m=m+1; stack[x++].c=count;
	}
      }
      if (l<width && r<height && count<=condition[m])
	continue;
      x=0;
      stack[x].i=i+1; stack[x].j=j; stack[x].m=0; stack[x++].c=0;
      while (x>0) {
	for (count=stack[--x].c, m=stack[x].m, l=stack[x].i, r=stack[x].j; 
	     l<width && r<height && badness[a[l]][b[r]]==0
	       && count>condition[m];  
	     m++, l++, r++) ;
	if (l>=width || r>=height)
	  continue;
	if (count<=condition[m])
	  break;
	if ((count+=badness[a[l]][b[r]])<VERYBAD) {
	  stack[x].i=l; stack[x].j=r+1; stack[x].m=m; stack[x++].c=count;
	  stack[x].i=l+1; stack[x].j=r; stack[x].m=m+1; stack[x++].c=count;
	  stack[x].i=l+1; stack[x].j=r+1; stack[x].m=m+1; stack[x++].c=count;
	}
      }
      if (l<width && r<height && count<=condition[m]) {
	j--;
	continue;
      }
      x=0;
      stack[x].i=i; stack[x].j=j+1; stack[x].m=0; stack[x++].c=0;
      while (x>0) {
	for (count=stack[--x].c, m=stack[x].m, l=stack[x].i, r=stack[x].j; 
	     l<width && r<height && badness[a[l]][b[r]]==0
	       && count>condition[m];  
	     m++, l++, r++) ;
	if (l>=width || r>=height)
	  continue;
	if (count<=condition[m])
	  break;
	if ((count+=badness[a[l]][b[r]])<VERYBAD) {
	  stack[x].i=l; stack[x].j=r+1; stack[x].m=m; stack[x++].c=count;
	  stack[x].i=l+1; stack[x].j=r; stack[x].m=m+1; stack[x++].c=count;
	  stack[x].i=l+1; stack[x].j=r+1; stack[x].m=m+1; stack[x++].c=count;
	}
      }
      if (l<width && r<height && count<=condition[m]) {
	i--;
	continue;
      }
      break;
    }
  span=i-left;

  /* extend toward the left, if possible */
  for (i=left-1, j=i-diff; i>=0 && j>=0; i--, j--) 
    if (badness[a[i]][b[j]]) {
      x=0;
      stack[x].i=i-1; stack[x].j=j-1; stack[x].m=0; stack[x++].c=0;
      while (x>0) {
	for (count=stack[--x].c, m=stack[x].m, l=stack[x].i, r=stack[x].j; 
	     l>=0 && r>=0 && badness[a[l]][b[r]]==0
	       && count>condition[m];  
	     m++, l--, r--) ;
	if (l<0 || r<0)
	  continue;
	if (count<=condition[m])
	  break;
	if ((count+=badness[a[l]][b[r]])<VERYBAD) {
	  stack[x].i=l; stack[x].j=r-1; stack[x].m=m; stack[x++].c=count;
	  stack[x].i=l-1; stack[x].j=r; stack[x].m=m+1; stack[x++].c=count;
	  stack[x].i=l-1; stack[x].j=r-1; stack[x].m=m+1; stack[x++].c=count;
	}
      }
      if (l>=0 && r>=0 && count<=condition[m])
	continue;
      x=0;
      stack[x].i=i-1; stack[x].j=j; stack[x].m=0; stack[x++].c=0;
      while (x>0) {
	for (count=stack[--x].c, m=stack[x].m, l=stack[x].i, r=stack[x].j; 
	     l>=0 && r>=0 && badness[a[l]][b[r]]==0
	       && count>condition[m];  
	     m++, l--, r--) ;
	if (l<0 || r<0)
	  continue;
	if (count<=condition[m])
	  break;
	if ((count+=badness[a[l]][b[r]])<VERYBAD) {
	  stack[x].i=l; stack[x].j=r-1; stack[x].m=m; stack[x++].c=count;
	  stack[x].i=l-1; stack[x].j=r; stack[x].m=m+1; stack[x++].c=count;
	  stack[x].i=l-1; stack[x].j=r-1; stack[x].m=m+1; stack[x++].c=count;
	}
      }
      if (l>=0 && r>=0 && count<=condition[m]) {
	j++;
	continue;
      }
      x=0;
      stack[x].i=i; stack[x].j=j-1; stack[x].m=0; stack[x++].c=0;
      while (x>0) {
	for (count=stack[--x].c, m=stack[x].m, l=stack[x].i, r=stack[x].j; 
	     l>=0 && r>=0 && badness[a[l]][b[r]]==0
	       && count>condition[m];  
	     m++, l--, r--) ;
	if (l<0 || r<0)
	  continue;
	if (count<=condition[m])
	  break;
	if ((count+=badness[a[l]][b[r]])<VERYBAD) {
	  stack[x].i=l; stack[x].j=r-1; stack[x].m=m; stack[x++].c=count;
	  stack[x].i=l-1; stack[x].j=r; stack[x].m=m+1; stack[x++].c=count;
	  stack[x].i=l-1; stack[x].j=r-1; stack[x].m=m+1; stack[x++].c=count;
	}
      }
      if (l>=0 && r>=0 && count<=condition[m]) {
	i++;
	continue;
      }
      break;
    }
  span=left+span-i-1;
  left=i+1;
  
  *start=left;
  *end=left+span-1;
}
