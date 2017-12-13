#pragma once

#ifndef DICTIONARY_H
#define	DICTIONARY_H

/* 
 * File:   Dictionary.h
 * Author: ilya
 *
 * Created on 7 Август 2012 г., 0:48
 */



#include <stdio.h>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <pthread.h>
#include "util.h"
#include <streambuf>
#include <exception>
#include <pthread.h>
#include <math.h>
//#include <sparsehash/sparse_hash_map>
//#include <sparsehash/sparse_hash_set>
//#include <sparsehash/dense_hash_set>
#include "Read.h"
#include "KMerRoutine.h"
//#include "rlmid.h"

using namespace std;

/*Extern variables*/
extern map<string, vector<k_mer_struct> > VectorDict;
extern map<long /*seq_id*/, string /*sequence*/ > VectorSeqs;
extern map<string, vector<k_mer_struct> > ContDict;//[24];
extern map<string, vector<k_mer_struct> >::iterator it_ContDict;
extern vector<RL_MID> rlmids;
extern short KMER_SIZE;
extern unsigned short NUM_THREADS;
extern map<int, string > vector_names;
extern unsigned short KMER_SIZE_CONT;
extern fstream sum_stat;


typedef struct {
   string line;
   long tid;
   long rec_id;
   //string rec_id;
   string seq;
   int w;
   long line_cnt;
   bool *clip_found;
   bool reverse_complement;
   bool reverse;
   bool complement;
} t_args;

 
/*Builds a new Vector Dictionary. Here it assumes that frequency is 1*/
int BuildVectorDictionary(string filename);
void ParseString(string str, int rec_id);

/*Builds a new Contaminants Dictionary. Here it assumes that frequency is 1*/
int BuildContDictionary(string filename);

void PutContKmer(string str, int rec_id);

void Build_RLMIDS_Dictionary(char* rlmids_file);
void Build_RLMIDS_Dictionary();
/*For Illumina TrueSeq adapters*/
string GetTruSeqAdapter(string type, short index);

#endif	/* DICTIONARY_H */

