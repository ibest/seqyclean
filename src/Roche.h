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
#include "poly.h"
#include "dup.h"

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
extern float max_a_error;
extern float max_e_at_ends;
extern bool sff_file_flag;
extern bool fastq_file_flag;
extern fstream sum_stat, sum_stat_tsv;
extern bool lucy_only_flag;
extern unsigned short cdna, c_err, crng;
extern bool fasta_output;
extern bool detailed_report;
extern bool rem_dup;

int RocheRoutine();
void ParseFastqFile(char* fastq_file, vector<Read*> &reads);
void QualityTrimming( vector<Read*>& reads );
void RemoveContaminants454(vector<Read*>& reads454);
void MakeClipPoints();
string PrintRocheStatisticsTSV(unsigned long cnt,
                                    unsigned long long bases_anal, 
                                    unsigned long left_mid_tag,
                                    unsigned long right_mid_tag,
                                    unsigned long num_vectors, 
                                    unsigned long num_contaminants, 
                                    unsigned long left_trimmed_by_adapter,
                                    unsigned long left_trimmed_by_quality,
                                    unsigned long left_trimmed_by_vector,
                                    double avg_left_trim_len, 
                                    unsigned long right_trimmed_by_adapter, 
                                    unsigned long right_trimmed_by_quality,
                                    unsigned long right_trimmed_by_vector,
                                    double avg_right_trim_len,
                                    unsigned long discarded, 
                                    unsigned long discarded_by_contaminant, 
                                    unsigned long discarded_by_read_length,
                                    unsigned long accepted, 
                                    double avg_trim_len,
                                    unsigned long left_trimmed_by_polyat,
                                    unsigned long right_trimmed_by_polyat
                               );
void MakeFinalStatistics(fstream &sum_stat);
void WriteToFASTQ(string file_name);
void WriteToSFF(string file_name);

#endif	/* ROCHE_H */

