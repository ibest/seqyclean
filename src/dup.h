// -*- C++ -*-
/* 
 * File:   dup.h
 * Author: ilya
 *
 * Created on July 16, 2013, 3:10 PM
 * 
 * Header file for duplicates screening
 */

#ifndef DUP_H
#define	DUP_H

#include <stdio.h>
#include "Read.h"
#include "util.h"
#include <map>

using namespace std;

void screen_duplicates(Read *read1, Read *read2, long long &duplicates);
void screen_duplicates(Read *read1, long long &duplicates);

// Starting position and a window size for searching for duplicates:
extern int start_dw, size_dw, max_dup;

#endif	/* DUP_H */

