/* 
 * File:   abi.h
 * Author: ilya
 *
 * Created on 4 Ноябрь 2012 г., 12:22
 */


#ifndef ABI_H
#define	ABI_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

//using namespace std;

//#include <ctype.h>


//#ifdef	__cplusplus
//extern "C" {
//#endif


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

int abi_code(int c);
void abi_sort(abi_struct *l, abi_struct *r);
void prepare_abi_mask();
void abi_align(char *a, int width, char *b, int height, int *start, int *end);



//#ifdef	__cplusplus
//}
//#endif





#endif	/* ABI_H */

