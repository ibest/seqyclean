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

extern unsigned int adapterlength;
extern double overlap_t;
extern unsigned int minoverlap;
extern bool overlap_flag;
extern bool i64_flag;

int find_overlap_pos(std::string seq1, std::string seq2, int adapterlength);
int strdist(std::string s1, std::string s2);
Read *make_consensus(Read *seq1, Read *seq2);
#endif	/* FLASH_H */

