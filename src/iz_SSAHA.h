/* 
 * File:   iz_SSAHA.h
 * Author: ilya
 *
 * Created on 1 Август 2012 г., 14:50
 */

#ifndef IZ_SSAHA_H
#define	IZ_SSAHA_H

#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <math.h>
#include <exception>
#include <math.h>
//#include <list>
#include <algorithm>
#include "pairwise.h"
//#include "util.h"

using namespace std;


typedef struct {
        int index;
        int shift;
        int offset;
} Master_struct;

typedef struct {
        int seq_id;
        int pos;
        int t;
        int cnt;
} Hit_struct;

extern void stoupper(std::string& s);

class iz_SSAHA {
public:
    iz_SSAHA();
    iz_SSAHA(const iz_SSAHA& orig);
    virtual ~iz_SSAHA();
    
    void AddElementToMasterStruct( Master_struct a);
    
    AlignResult Find(string &read, string &query_str);
    
    
private:
    void RunSort();
    
    vector<Master_struct> m_structs;
    //vector<Master_struct> h_structs;
    
    map<string, vector<Hit_struct> > Table1;
    map<string, vector<Hit_struct> >::iterator it_Table1;
    
    void MakeTable1(string read);

};

#endif	/* IZ_SSAHA_H */

