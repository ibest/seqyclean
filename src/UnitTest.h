/* 
 * File:   UnitTest.h
 * Author: kwt
 *
 * Created on July 28, 2013, 11:18 AM
 */

#ifndef UNITTEST_H
#define	UNITTEST_H

#include <stdio.h>
#include <iostream>
#include "Roche.h"
#include "Illumina.h"

using namespace std;

/*For Roche*/
extern bool trim_adapters_flag;
extern string roche_output_file_name;
extern string roche_rep_file_name;
extern bool qual_trim_flag;
extern string output_prefix;
extern long discard_counter;
extern long accept_counter;
extern long trim_counter;
extern float max_a_error;
extern float max_e_at_ends;
extern bool sff_file_flag;
extern bool fastq_file_flag;
extern fstream sum_stat, sum_stat_tsv;
extern bool lucy_only_flag;


/*Illumina*/
extern char* illumina_file_name_R1;
extern char* illumina_file_name_R2;
extern char* illumina_file_name_se;
/*Illumina data structures*/
extern char *output_file_name;
extern bool old_style_illumina_flag;
extern int phred_coeff_illumina; //by default assume new illumina (1.8)
extern bool i64_flag;

extern int adapterlength;
extern double overlap_t;
extern int minoverlap;
extern bool overlap_flag;

extern bool overwrite_flag;

extern bool rem_dup;

extern vector<string> file_list;

extern bool VectorOnlyFlag;
extern bool new2old_illumina;
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




/*----------End of output data definition------------------*/



class UnitTest {
public:
    UnitTest();
    UnitTest(const UnitTest& orig);
    virtual ~UnitTest();
    
    void Test454();
    void TestIlluminaPE();
    void TestIlluminaSE();
    
private:

};

#endif	/* UNITTEST_H */

