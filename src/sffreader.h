/* 
 * File:   sffreader.h
 * Author: ilya
 *
 * Created on 28 Август 2012 г., 9:59
 */

#ifndef SFFREADER_H
#define	SFFREADER_H

//#include <R.h>
//#include <Rdefines.h>
//#include <Rinternals.h> // Rprintf, SEXP
//#include <R_ext/Rdynload.h>
//#include "IRanges_interface.h"
//#include "Biostrings_interface.h"

#define MATHLIB_STANDALONE
// System includes
#include <zlib.h>
#include <stdint.h>		// uint64_t, uint32_t, uint16_t
#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "util.h"
#include "sff.h"
#include "Read.h"

//#include <RInside.h> 

using namespace std;

extern vector<Read*> reads;

sff_common_header h;
int discarded_reads = 0;
extern bool keep_fastq_orig;




void process_sff_to_fastq(char *sff_file, int trim_flag);
void construct_fastq_entry(FILE *fp,
                           char *name,
                           char *bases,
                           uint8_t *quality,
                           int nbases);

void process_fastq_to_sff(char *sff_file);

#endif