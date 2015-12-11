#ifndef SFFREADER_H
#define	SFFREADER_H

/* 
 * File:   sffreader_lin.h
 * Author: ilya
 */

#define MATHLIB_STANDALONE
#define _BSD_SOURCE

#include <endian.h>

// System includes
#include <zlib.h>
#include <stdint.h>		// uint64_t, uint32_t, uint16_t
#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "util.h"
#include "sff_lin.h"
#include "Read.h"

using namespace std;

extern vector<Read*> reads;
extern bool keep_fastq_orig;

void process_sff_to_fastq(char *sff_file, int trim_flag);
void construct_fastq_entry(FILE *fp,
                           char *name,
                           char *bases,
                           uint8_t *quality,
                           int nbases);

void process_fastq_to_sff(char *sff_file);

#endif
