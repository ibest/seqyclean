/* 
 * File:   util.h
 * Author: ilya
 *
 * Created on 31 Май 2012 г., 10:27
 */

#ifndef UTIL_H
#define	UTIL_H

#include <stdio.h>
#include <stdint.h>
#include <arpa/inet.h> // htons(), htonl()
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <sys/stat.h>
#include <errno.h>
#include <sys/types.h>
#include <dirent.h>
#include <map>
#include <math.h>

#define SFF_MAGIC_NUM   0x2e736666 /* ".sff" */


#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })


using namespace std;

typedef struct {
    uint32_t  magic; /*Magic number: identical for all sff files (0x2E736666, the newbler manual explains that this is ‘the uint32_t encoding of the string „.sff“ ‘)*/
    char      version[4]; /*Version: also identical for all sff files (0001)*/
    uint64_t  index_offset; /*Index offset and length: has to do with the index of the binary sff file (points to the location of the index in the file)*/
    uint32_t  index_len;
    uint32_t  nreads; /* # of reads: stored in the sff file*/
    uint16_t  header_len;/*Header length: looks like it is 440 for GS FLX reads, 840 for GS FLX Titanium reads*/
    uint16_t  key_len;/*Key length: the length (in bases) of the key sequence that each read starts with, so far always 4*/
    uint16_t  flow_len;/*# of Flows: each flow consists of a base that is flowed over the plate; for GS20, there were 168 flows (42 cycles of all four nucleotides), 400 for GS FLX (100 cycles) and 800 for Titanium (200 cycles)*/
    uint8_t   flowgram_format;/*Flowgram code: kind of the version of coding the flowgrams (signal strengths); so far, ’1′ for all sff files*/
    char     *flow;/*Flow Chars: a string consisting of’ # of flow’ characters (168, 400 or 800) of the bases in flow order (‘TACG’ up to now)*/
    char     *key; /* Key Sequence: the first four bases of reads are either added during library preparation (they are the last bases of the ‘A’ adaptor) or they are a part of the control beads. For example, Titanium sample beads have key sequence TACG (default library protocol) or GACT (rapid library protocol), control beads have CATG or ATGC. Control reads never make it into sff files*/
} sff_c_header;



void stoupper(std::string& s);
char* itoa(int value, char* result, int base);
string MakeSeqComplement(string init_str);
string MakeRevComplement(string init_str);
//int GetRandomInt (int from, int to);
//double GetRandomDouble (double from, double to);
void split_str(const string& str, vector<string>& tokens, const string& delimiters);
void TrimNs(string &read);
void TrimNs2(string &read);
string GenNs(int num, char* letter);
int min4( int x, int y, int z, int k );
int max4( int x, int y, int z, int k );
int max3( int x, int y, int z );
int min3( int x, int y, int z );
bool exists(char* filePath);
double GetAvg( double past_avg, long n, int cur_diff, int first_avg );
int MakeDirectory(string path_to_create);
void GetDirectories(std::vector<string> &out, char *directory);
vector< vector<string> > GetPEfilenames(string prefix1, string prefix2, char *directory);
string i2str(unsigned long long value, char* result, int base);
string double2str(double num);
string int2str(int num);
short GetFormat(char* filename);
void  get_sff_common_header(FILE *fp, sff_c_header *h);
short  check_sff_common_header(sff_c_header *h);
void free_sff_c_header(sff_c_header *h);
uint64_t BE64toNE(uint64_t bigEndian);


#endif	/* UTIL_H */

