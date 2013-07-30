/* 
 * File:   Roche.h
 * Author: kwt
 *
 * Created on July 28, 2013, 11:20 AM
 */

#ifndef ROCHE_H
#define	ROCHE_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include "Read.h"
#include "sffreader.h"
#include "QualTrim.h"
#include "Report.h"
#include "MainPipeLine.h"
#include "Dictionary.h"
#include "gzstream.h"

using namespace std;

extern bool trim_adapters_flag;
extern vector<char*> roche_names;
extern string roche_output_file_name;
extern string roche_rep_file_name;
extern vector<Read*> reads;
extern vector<RL_MID> rlmids;
extern bool qual_trim_flag;
extern string output_prefix;
extern bool debug_flag;
extern bool custom_rlmids_flag;
extern char* rlmids_file;
extern long discard_counter;
extern long accept_counter;
extern long trim_counter;
extern float max_a_error;
extern float max_e_at_ends;
extern bool sff_file_flag;
extern bool fastq_file_flag;
extern fstream sum_stat, sum_stat_tsv;
extern bool lucy_only_flag;

void RocheRoutine();
void ParseFastqFile(char* fastq_file, vector<Read*> &reads);
void QualityTrimming( vector<Read*>& reads );
void RemoveContaminants454(vector<Read*>& reads454);

#endif	/* ROCHE_H */

