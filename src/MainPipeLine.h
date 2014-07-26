#ifndef MAINPIPELINE_H
#define	MAINPIPELINE_H

/* 
 * File:   MainPipeLine.h
 * Author: ilya
 *
 * Created on 10 Август 2012 г., 20:48
 */



#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include "timer.h"
#include <pthread.h>
#include "util.h"
#include <streambuf>
#include <exception>
#include "pairwise.h"
#include <pthread.h>
#include <dirent.h>
#include "iz_SSAHA.h"
#include "Read.h"
#include "KMerRoutine.h"


/*Extern variables*/
extern vector<Read*> reads;
extern long line_counter;
extern short KMER_SIZE;
extern vector<RL_MID> rlmids;
extern unsigned short NUM_THREADS;
extern bool contaminants_flag;
extern bool vector_flag;
extern bool amplicon_flag;
extern long counter;
extern long discard_counter;
extern long accept_counter;
extern long trim_counter;
extern bool pcr_flag;
extern char *pcr_file_name;
extern unsigned short minimum_read_length;
extern int max_al_mism;
/*Map that holds a whole vector genomes :*/
extern map<long /*seq_id*/, string /*sequence*/ > VectorSeqs;


void MainPipeLine();
void TrimRightEnds();
void TrimLeftEnds();
//static void *t_FindRLClip(void *targs);
//static void *t_TrimRightEnds(void *targs);
//static void *t_TrimLeftEnds(void *targs);
void GetLClip2(Read* read, bool pflag);

void SSAHA(string query_str, int read_len, string ref_str, bool reverse, bool complement, bool reverse_complement,  bool clip_found, long rec_id, int tid );

extern bool qual_trim_flag;

#endif	/* MAINPIPELINE_H */

