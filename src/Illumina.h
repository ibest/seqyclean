/* 
 * File:   Illumina.h
 * Author: ilya
 *
 * Created on 20 Сентябрь 2012 г., 16:11
 */

#ifndef ILLUMINA_H
#define	ILLUMINA_H


#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <streambuf>
#include <vector>
#include <pthread.h>
#include <exception>
#include "iz_SSAHA.h"
#include "pairwise.h"
#include "Dictionary.h"
#include "util.h"
#include "Read.h"


//#include <sparsehash/dense_hash_map>

using namespace std;
//using google::dense_hash_map; 
/*Illumina data structures*/
//extern dense_hash_map<string, Read> reads_1; extern dense_hash_map<string, Read>::iterator it_reads_1;
//extern dense_hash_map<string, Read> reads_2; extern dense_hash_map<string, Read>::iterator it_reads_2;
//extern map<string, Read> reads_1; extern map<string, Read>::iterator it_reads_1;
//extern map<string, Read> reads_2; extern map<string, Read>::iterator it_reads_2;
extern vector<Read*> reads_1;
extern vector<Read*> reads_2;
extern int max_al_mism;
extern bool qual_trim_flag;
extern long discard_counter;
extern bool contaminants_flag;
extern bool vector_flag;
string adapter_type_R1;
string adapter_type_R2;
string query_str1;
string query_str2;

typedef struct {
    map<string, Read>::iterator it_reads_;
    int pos;
    string readID;
} Targs_ill;


class Illumina {
public:
    Illumina();
    Illumina(const Illumina& orig);
    virtual ~Illumina();
    
    void CleanSeq_PE();
    void CleanSeq_SE();
    
    char* illumina_file_name_R1;
    char* illumina_file_name_R2;
    
private:
    
    void t_CleanSeq_PE();
    void t_CleanSeq_SE();
    
    
    void* tt_CleanSeq_1(void* targ);
    static void* tt_StaticThreadProc_1(void *arg) {
        return reinterpret_cast<Illumina*>(arg)->tt_CleanSeq_1(arg);
    }
    
    void* tt_CleanSeq_2(void* targ);
    static void* tt_StaticThreadProc_2(void *arg) {
        return reinterpret_cast<Illumina*>(arg)->tt_CleanSeq_2(arg);
    }
    
    int CheckContaminants(string seq);

    void t_CleanSeq_PE_Dynamic(Read &read);
};

#endif	/* ILLUMINA_H */

