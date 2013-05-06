#ifndef PAIRWISE_H
#define	PAIRWISE_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string.h>
#include <assert.h>
#include <algorithm>
#include <stdexcept>

#define EG 7;//End gap penalty
#define IG 10;//internal gap

using namespace std;

struct AlignResult
{
  short scores;
  string seq_1_al;
  string seq_2_al;
  string seq_1;
  string seq_2;
  long tid;
  int w;
  short read_len;
  short n_mismatches;
  int pos;
  bool found_flag;
  int pos_left;
  int pos_right;
};

struct AlignScores
{
  short scores;
  short mismatches;
  int trim_pos;
};

void  print_matrix( int ** S, string seq_1, string seq_2 );
bool InsideBand(int i, int j, int k);		   
int match(char &a, char &b);
int get_max2(int a, int b);

//Banded alignment method:
/*
 * seq_1 : read
 * seq_2 : sequence to search for
 * seq_1_al : aligned sequence 1 (seq_1)
 * seq_2_al : aligned sequence 2 (seq_2)
 * k : width of the band
 * d : gap penalty
 */

//Score calculations:
AlignScores CalcScores(string &seq_1, string &seq_2,int lim,int dir);
AlignScores CalcScores2(string &seq_1, int lim, int dir);
AlignResult CalcPos(string &seq);

//AlignResult banded(string seq_1, string seq_2, int k, int d, int tid,int w);

AlignResult banded(string &seq_1, string &seq_2, int k, int d);

#endif