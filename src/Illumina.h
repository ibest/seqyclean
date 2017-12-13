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
#include "Dictionary.h"
#include "Report.h"
#include "MainPipeLine.h"
#include "gzstream.h"
#include "QualTrim.h"
#include "flash.h"
#include "dup.h"
#include "Dictionary.h"
#include "poly.h"
#include <algorithm>
#include "util.h"

using namespace std;

extern bool illumina_flag;
extern bool illumina_flag_se;
extern char* illumina_file_name_R1;
extern char* illumina_file_name_R2;
extern char* illumina_file_name_se;
extern std::string adapter_type_R1;
extern std::string adapter_type_R2;
extern std::string query_str1;
extern std::string query_str2;
extern bool polyat_flag;
extern bool verbose;
extern bool detailed_report;
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
extern short phred_coeff_illumina; //by default assume new illumina (1.8)
extern bool i64_flag;


extern double overlap_t;
extern unsigned int minoverlap;
extern bool overlap_flag;
extern unsigned int adapterlen;

extern bool overwrite_flag;

extern bool rem_dup;

extern bool dynflag;

extern std::vector<std::string> file_list;

extern std::fstream sum_stat, sum_stat_tsv;

extern std::string output_prefix;

extern bool VectorOnlyFlag;
extern bool new2old_illumina;



extern bool serial_flag;

extern volatile int shared_var;

extern bool shuffle_flag;

/*Report files*/
extern std::string rep_file_name1, rep_file_name2, pe_output_filename1, pe_output_filename2, shuffle_filename, se_filename, se_output_filename, overlap_file_name;

extern  bool wildcart_search_flag;

extern  std::vector<char*> pe1_names, pe2_names, roche_names, se_names;

extern int window0;
extern int window1;

extern unsigned short cdna, c_err, crng;

extern bool trim_adapters_flag;

void IlluminaDynamic();
int TrimIlluminaSE(Read* read, bool trim_adapter);
void WritePEFile(std::fstream &pe_output_file, Read *read);
void WritePEFileGZ(ogzstream &pe_output_file, Read *read);
void WriteShuffleFile(std::fstream &shuffle_output_file, Read *read1, Read *read2);
void WriteShuffleFileGZ(ogzstream &shuffle_output_file, Read *read1, Read *read2);
void WriteSEFile(std::fstream &se_output_file, Read *read);
void WriteSEFileGZ(ogzstream &se_output_file, Read *read);
void WriteSEOverlap(std::fstream &overlap_file, Read *read);
void WriteSEOverlapGZ(ogzstream &overlap_file, Read *read);
string New2OldNbl(std::string header);
void IlluminaDynamicSE();
void trim_read(Read *read);
void update_statistics(Read *read1, Read *read2);
void update_counters_and_print_statistics(Read *read1, Read *read2);
bool align_ssaha(Read* read, std::string &query_str );

int IlluminaDynRoutine_post(Read* read);


std::string PrintIlluminaStatistics(long long cnt1, long long cnt2, 
                                    unsigned long long pe1_bases_anal, unsigned long long pe2_bases_anal, 
                                    long long ts_adapters1, long long ts_adapters2, 
                                    long long num_vectors1, long long num_vectors2, 
                                    long long num_contaminants1, long long num_contaminants2, 
                                    long long left_trimmed_by_quality1, long long left_trimmed_by_quality2,
                                    long long left_trimmed_by_vector1, long long left_trimmed_by_vector2, 
                                    double avg_left_trim_len_pe1, double avg_left_trim_len_pe2, 
                                    long long right_trimmed_by_adapter1, long long right_trimmed_by_adapter2, 
                                    long long right_trimmed_by_quality1,long long right_trimmed_by_quality2,
                                    long long right_trimmed_by_vector1,long long right_trimmed_by_vector2,
                                    double avg_right_trim_len_pe1, double avg_right_trim_len_pe2,
                                    long long discarded1, long long discarded2,
                                    long long discarded_by_contaminant1, long long discarded_by_contaminant2,
                                    long long discarded_by_read_length1, long long discarded_by_read_length2,
                                    long long pe_accept_cnt, unsigned long long  pe_bases_kept, 
                                    long long pe_discard_cnt,unsigned long long pe_bases_discarded, 
                                    long long se_pe1_accept_cnt, unsigned long long  se_pe1_bases_kept,
                                    long long se_pe2_accept_cnt, unsigned long long se_pe2_bases_kept,
                                    double avg_trim_len_pe1, double avg_trim_len_pe2,
                                    double avg_len_pe1, double avg_len_pe2,
                                    long long perfect_ov_cnt, long long partial_ov_cnt,
                                    long long duplicates,
                                    long long left_trimmed_by_polyat1, long long right_trimmed_by_polyat1,
                                    long long left_trimmed_by_polyat2, long long right_trimmed_by_polyat2,
                                    long long left_trimmed_by_adapter1, long long left_trimmed_by_adapter2
                                    );

std::string PrintIlluminaStatisticsTSV(long long cnt1, long long cnt2, 
                                    unsigned long long pe1_bases_anal, unsigned long long pe2_bases_anal, 
                                    long long ts_adapters1, long long ts_adapters2, 
                                    long long num_vectors1, long long num_vectors2, 
                                    long long num_contaminants1, long long num_contaminants2, 
                                    long long left_trimmed_by_quality1, long long left_trimmed_by_quality2,
                                    long long left_trimmed_by_vector1, long long left_trimmed_by_vector2, 
                                    double avg_left_trim_len_pe1, double avg_left_trim_len_pe2, 
                                    long long right_trimmed_by_adapter1, long long right_trimmed_by_adapter2, 
                                    long long right_trimmed_by_quality1,long long right_trimmed_by_quality2,
                                    long long right_trimmed_by_vector1,long long right_trimmed_by_vector2,
                                    double avg_right_trim_len_pe1,double avg_right_trim_len_pe2,
                                    long long discarded1, long long discarded2,
                                    long long discarded_by_contaminant1, long long discarded_by_contaminant2,
                                    long long discarded_by_read_length1, long long discarded_by_read_length2,
                                    long long pe_accept_cnt, unsigned long long pe_bases_kept, 
                                    long long pe_discard_cnt,unsigned long long pe_bases_discarded, 
                                    long long se_pe1_accept_cnt, unsigned long long se_pe1_bases_kept,
                                    long long se_pe2_accept_cnt, unsigned long long se_pe2_bases_kept,
                                    double avg_trim_len_pe1, double avg_trim_len_pe2,
                                    double avg_len_pe1, double avg_len_pe2,
                                    long long perfect_ov_cnt, long long partial_ov_cnt,
                                    long long left_trimmed_by_polyat1, long long right_trimmed_by_polyat1,
                                    long long left_trimmed_by_polyat2, long long right_trimmed_by_polyat2,
                                    long long duplicates);

std::string PrintIlluminaStatisticsSE(long long cnt, unsigned long long se_bases_anal, 
                                    long long ts_adapters,
                                    long long num_vectors,
                                    long long num_contaminants, 
                                    long long left_trimmed_by_quality,
                                    long long left_trimmed_by_vector, 
                                    double avg_left_trim_len_se,
                                    long long right_trimmed_by_adapter,
                                    long long right_trimmed_by_quality,
                                    long long right_trimmed_by_vector,
                                    double avg_right_trim_len_se,
                                    long long discarded, 
                                    long long discarded_by_contaminant,
                                    long long discarded_by_read_length,
                                    long long se_accept_cnt, unsigned long long se_bases_kept, 
                                    long long se_discard_cnt,unsigned long long  se_bases_discarded, 
                                    double avg_trim_len_se,
                                    double avg_len_se,
                                    long long left_trimmed_by_polyat, long long right_trimmed_by_polyat,
                                    long long discarded_by_polyAT,
                                    long long duplicates
                                    );

std::string PrintIlluminaStatisticsTSVSE(long long cnt,
                                    unsigned long long se_bases_anal, 
                                    long long ts_adapters, 
                                    long long num_vectors,  
                                    long long num_contaminants, 
                                    long long left_trimmed_by_quality, 
                                    long long left_trimmed_by_vector, 
                                    double avg_left_trim_len_se, 
                                    long long right_trimmed_by_adapter, 
                                    long long right_trimmed_by_quality,
                                    long long right_trimmed_by_vector,
                                    double avg_right_trim_len_se,
                                    long long discarded, 
                                    long long discarded_by_contaminant, 
                                    long long discarded_by_read_length,
                                    long long se_accept_cnt, 
                                    double avg_trim_len_se,
                                    long long left_trimmed_by_polyat, long long right_trimmed_by_polyat,
                                    long long discarded_by_polyAT,
                                    long long duplicates
                                   );




bool TrimAdapterPE(Read *read1, Read *read2);
int TrimIllumina(Read* read1, Read* read2);
int TrimAdapterSE(Read* read);
void LoadAdapters(std::string filename, bool custom_adapters);

#endif	/* ILLUMINA_H */

