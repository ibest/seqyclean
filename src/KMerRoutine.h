#ifndef KMERROUTINE_H
#define	KMERROUTINE_H

/* 
 * File:   KMerRoutine.h
 * Author: ilya
 *
 * Created on 7 Август 2012 г., 23:10
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
#include "Read.h"
#include "pairwise.h"

using namespace std;



typedef struct  {
    //string seq_id; //seq id related to k_mer in the screening file
    long seq_id;
    int pos; //position of the k_mer in the screening file
} k_mer_struct;


/*Structure that hold seq_id and position of k_mer*/
typedef struct {
    vector <k_mer_struct> kmers;
    long pos;
    string k_mer_string;
    string string_to_align;
    
} HitData;


typedef struct {
    int pos;
    int flag;
    int i;
} DistanceStruct;

/*Extern variables*/
extern map<string, vector<k_mer_struct> > VectorDict;
extern map<long /*seq_id*/, string /*sequence*/ > VectorSeqs;
extern map<string, vector<k_mer_struct> > ContDict;
extern map<string, vector<k_mer_struct> >::iterator it_ContDict;
extern short KMER_SIZE;
extern unsigned short KMER_SIZE_CONT;
extern unsigned short NUM_THREADS;
extern short DISTANCE;
extern unsigned short vmr;
extern unsigned short vml;
extern unsigned short L_limit;
extern unsigned short R_limit;
extern unsigned short allowable_distance;
extern unsigned short pmax;

int CheckVectorRight(Read &read);
int CheckVectorLeft(Read &read);
DistanceStruct CheckDistance(vector <HitData> matches/*, int dir*/);
bool CheckContig(vector <HitData> matches);
void GetRClip(HitData &hit_data, string sss, Read &read);
void GetLClip(HitData &hit_data, string sss, Read &read);
int CheckContaminants(string seq);
int CheckVectorLeftIllumina(Read &read) ;
int CheckVector(Read* read) ;

#endif	/* KMERROUTINE_H */

