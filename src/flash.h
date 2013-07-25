/* 
 * File:   flash.h
 * Author: ilya
 *
 * Created on July 9, 2013, 10:33 PM
 */

#ifndef FLASH_H
#define	FLASH_H

#include <stdio.h>
#include <string.h>
#include "Read.h"

using namespace std;

extern int adapterlength;
extern double overlap_t;
extern int minoverlap;
extern bool overlap_flag;

int find_overlap_pos(string seq1, string seq2, int adapterlength, bool flag);
int strdist(string s1, string s2);
Read *make_consensus(Read *seq1, Read *seq2);
#endif	/* FLASH_H */

