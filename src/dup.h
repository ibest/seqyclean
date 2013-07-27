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

extern map<string, int> DupDict;
void screen_duplicates(Read *read1, Read *read2, unsigned long &duplicates);

#endif	/* DUP_H */

