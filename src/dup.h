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

#include "Read.h"

extern map<string, int> DupDict;

void screen_duplicates(Read read1, Read read2);

#endif	/* DUP_H */

