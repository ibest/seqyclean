/* 
 * File:   Illumina.h
 * Author: kwt
 *
 * Created on July 28, 2013, 2:16 PM
 */

#ifndef ILLUMINA_H
#define	ILLUMINA_H

#include <stdio.h>
#include <iostream>
#include "Read.h"
#include "util.h"
#include "Dictionary.h"
#include "Report.h"
#include "MainPipeLine.h"
#include "gzstream.h"
#include "QualTrim.h"
#include "flash.h"
#include "dup.h"
#include "Dictionary.h"
#include "poly.h"

using namespace std;

extern bool illumina_flag;
extern bool illumina_flag_se;
extern char* illumina_file_name_R1;
extern char* illumina_file_name_R2;
extern char* illumina_file_name_se;
extern string adapter_type_R1;
extern string adapter_type_R2;
extern string query_str1;
extern string query_str2;
extern bool polyat_flag;
/*Illumina data structures*/



 
extern  bool output_sfffile_flag;
extern bool output_fastqfile_flag;
extern char *output_file_name;
extern char* custom_output_filename;
extern bool keep_fastq_orig;
extern bool lucy_only_flag;

/*-----LUCY parameters------*/
extern float max_a_error;
extern float max_e_at_ends;

extern bool old_style_illumina_flag;
extern int phred_coeff_illumina; //by default assume new illumina (1.8)
extern bool i64_flag;

extern unsigned int adapterlength;
extern double overlap_t;
extern int minoverlap;
extern bool overlap_flag;

extern bool overwrite_flag;

extern bool rem_dup;

extern bool dynflag;

extern vector<string> file_list;

extern fstream sum_stat, sum_stat_tsv;

extern string output_prefix;

extern bool VectorOnlyFlag;
extern bool new2old_illumina;



extern bool serial_flag;

extern volatile int shared_var;

extern bool shuffle_flag;

/*Report files*/
extern string rep_file_name1, rep_file_name2, pe_output_filename1, pe_output_filename2, shuffle_filename, se_filename, se_output_filename, overlap_file_name;

extern  unsigned long long se_bases_kept, se_bases_discarded;
extern  unsigned long se_discard_cnt;
extern  unsigned long long se_bases_anal;        
extern  unsigned long avg_trim_len_se;

extern  bool wildcart_search_flag;

extern  vector<char*> pe1_names, pe2_names, roche_names, se_names;

extern string stat_str, tsv_stat_str;

extern int window0;
extern int window1;

extern unsigned short cdna, c_err, crng;


void IlluminaDynamic();
void IlluminaDynamic();
int IlluminaDynRoutine(Read* read, bool& adapter_found, string& query_str);
void WritePEFile(fstream &pe_output_file, Read *read);
void WriteShuffleFile(fstream &shuffle_output_file, Read *read1, Read *read2);
void WriteSEFile(fstream &se_output_file, Read *read);
string New2OldNbl(string header);
void IlluminaDynamicSE();
void MakeClipPointsIllumina(Read* read);
void WriteSEOverlap(fstream &overlap_file, Read *read);


string PrintIlluminaStatistics(unsigned long cnt1, unsigned long cnt2, 
                                    unsigned long long pe1_bases_anal, unsigned long long pe2_bases_anal, 
                                    unsigned long ts_adapters1, unsigned long ts_adapters2, 
                                    unsigned long num_vectors1, unsigned long num_vectors2, 
                                    unsigned long num_contaminants1, unsigned long num_contaminants2, 
                                    unsigned long left_trimmed_by_quality1, unsigned long left_trimmed_by_quality2,
                                    unsigned long left_trimmed_by_vector1, unsigned long left_trimmed_by_vector2, 
                                    double avg_left_trim_len_pe1, double avg_left_trim_len_pe2, 
                                    unsigned long right_trimmed_by_adapter1, unsigned long right_trimmed_by_adapter2, 
                                    unsigned long right_trimmed_by_quality1,unsigned long right_trimmed_by_quality2,
                                    unsigned long right_trimmed_by_vector1,unsigned long right_trimmed_by_vector2,
                                    double avg_right_trim_len_pe1, double avg_right_trim_len_pe2,
                                    unsigned long discarded1, unsigned long discarded2,
                                    unsigned long discarded_by_contaminant1, unsigned long discarded_by_contaminant2,
                                    unsigned long discarded_by_read_length1, unsigned long discarded_by_read_length2,
                                    unsigned long pe_accept_cnt, unsigned long long  pe_bases_kept, 
                                    unsigned long pe_discard_cnt,unsigned long long pe_bases_discarded, 
                                    unsigned long se_pe1_accept_cnt, unsigned long long  se_pe1_bases_kept,
                                    unsigned long se_pe2_accept_cnt, unsigned long long se_pe2_bases_kept,
                                    double avg_trim_len_pe1, double avg_trim_len_pe2,
                                    double avg_len_pe1, double avg_len_pe2,
                                    unsigned long perfect_ov_cnt, unsigned long partial_ov_cnt,
                                    unsigned long duplicates,
                                    unsigned long left_trimmed_by_polyat1, unsigned long right_trimmed_by_polyat1,
                                    unsigned long left_trimmed_by_polyat2, unsigned long right_trimmed_by_polyat2
                                    );

string PrintIlluminaStatisticsTSV(unsigned long cnt1, unsigned long cnt2, 
                                    unsigned long long pe1_bases_anal, unsigned long long pe2_bases_anal, 
                                    unsigned long ts_adapters1, unsigned long ts_adapters2, 
                                    unsigned long num_vectors1, unsigned long num_vectors2, 
                                    unsigned long num_contaminants1, unsigned long num_contaminants2, 
                                    unsigned long left_trimmed_by_quality1, unsigned long left_trimmed_by_quality2,
                                    unsigned long left_trimmed_by_vector1, unsigned long left_trimmed_by_vector2, 
                                    double avg_left_trim_len_pe1, double avg_left_trim_len_pe2, 
                                    unsigned long right_trimmed_by_adapter1, unsigned long right_trimmed_by_adapter2, 
                                    unsigned long right_trimmed_by_quality1,unsigned long right_trimmed_by_quality2,
                                    unsigned long right_trimmed_by_vector1,unsigned long right_trimmed_by_vector2,
                                    double avg_right_trim_len_pe1,double avg_right_trim_len_pe2,
                                    unsigned long discarded1, unsigned long discarded2,
                                    unsigned long discarded_by_contaminant1, unsigned long discarded_by_contaminant2,
                                    unsigned long discarded_by_read_length1, unsigned long discarded_by_read_length2,
                                    unsigned long pe_accept_cnt, unsigned long long pe_bases_kept, 
                                    unsigned long pe_discard_cnt,unsigned long long pe_bases_discarded, 
                                    unsigned long se_pe1_accept_cnt, unsigned long long se_pe1_bases_kept,
                                    unsigned long se_pe2_accept_cnt, unsigned long long se_pe2_bases_kept,
                                    double avg_trim_len_pe1, double avg_trim_len_pe2,
                                    double avg_len_pe1, double avg_len_pe2,
                                    unsigned long perfect_ov_cnt, unsigned long partial_ov_cnt,
                                    unsigned long left_trimmed_by_polyat1, unsigned long right_trimmed_by_polyat1,
                                    unsigned long left_trimmed_by_polyat2, unsigned long right_trimmed_by_polyat2,
                                    unsigned long duplicates);

string PrintIlluminaStatisticsSE(unsigned long cnt, unsigned long long se_bases_anal, 
                                    unsigned long ts_adapters,
                                    unsigned long num_vectors,
                                    unsigned long num_contaminants, 
                                    unsigned long left_trimmed_by_quality,
                                    unsigned long left_trimmed_by_vector, 
                                    double avg_left_trim_len_se,
                                    unsigned long right_trimmed_by_adapter,
                                    unsigned long right_trimmed_by_quality,
                                    unsigned long right_trimmed_by_vector,
                                    double avg_right_trim_len_se,
                                    unsigned long discarded, 
                                    unsigned long discarded_by_contaminant,
                                    unsigned long discarded_by_read_length,
                                    unsigned long se_accept_cnt, unsigned long long se_bases_kept, 
                                    unsigned long se_discard_cnt,unsigned long long  se_bases_discarded, 
                                    double avg_trim_len_se,
                                    double avg_len_se,
                                    unsigned long left_trimmed_by_polyat, unsigned long right_trimmed_by_polyat,
                                    unsigned long discarded_by_polyAT
                                    );

string PrintIlluminaStatisticsTSVSE(unsigned long cnt,
                                    unsigned long long se_bases_anal, 
                                    unsigned long ts_adapters, 
                                    unsigned long num_vectors,  
                                    unsigned long num_contaminants, 
                                    unsigned long left_trimmed_by_quality, 
                                    unsigned long left_trimmed_by_vector, 
                                    double avg_left_trim_len_se, 
                                    unsigned long right_trimmed_by_adapter, 
                                    unsigned long right_trimmed_by_quality,
                                    unsigned long right_trimmed_by_vector,
                                    double avg_right_trim_len_se,
                                    unsigned long discarded, 
                                    unsigned long discarded_by_contaminant, 
                                    unsigned long discarded_by_read_length,
                                    unsigned long se_accept_cnt, 
                                    double avg_trim_len_se,
                                    unsigned long left_trimmed_by_polyat, unsigned long right_trimmed_by_polyat,
                                    unsigned long discarded_by_polyAT
                                   );






#endif	/* ILLUMINA_H */

