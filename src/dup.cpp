#include "dup.h"

//Version with strings & substrings. Not computationally efficient.
void screen_duplicates(Read *seq1, Read *seq2, unsigned long &duplicates) {
    int c = 0;
    
    if( (seq1->initial_length < 35) && (seq2->initial_length < 35) ) return;
    
    string comb = seq1->read.substr(10,25) + seq2->read.substr(10,25);
    
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
