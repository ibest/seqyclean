#include "dup.h"

map<string, int > DupDict;

//Version with strings & substrings. Not computationally efficient.
void screen_duplicates(Read *seq1, Read *seq2, unsigned long long &duplicates) {
    int c = 0;
    
    if( (seq1->initial_length < start_dw+size_dw) && (seq2->initial_length < start_dw+size_dw) ) return;
    
    string comb = seq1->read.substr(start_dw,size_dw) + seq2->read.substr(start_dw,size_dw);
    
    string rcomb = MakeRevComplement(comb);
    map<string, int>::iterator it_DupDict;
    it_DupDict = DupDict.find(comb);
    if(it_DupDict != DupDict.end()) {
     it_DupDict->second += 1;
     c = it_DupDict->second;
    }else{
     DupDict.insert(std::pair<string, int>(comb, 1));
     c = 1;
    }
    it_DupDict = DupDict.find(rcomb);
    if(it_DupDict != DupDict.end()) {
     it_DupDict->second += 1;
     c = it_DupDict->second;
    }else{
     DupDict.insert(std::pair<string, int>(rcomb, 1));
     c = 1;
    }
    if( c == 1 ) {
     //cout << "seq1 fastq" << endl << "seq2 fastq" << endl;
    }else {
     duplicates += 1;
     seq1->discarded = 1;
     seq2->discarded = 1;
    }
}
