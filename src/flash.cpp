#include "flash.h"

//modify this to down-weight difference in low-quality bases
int strdist(string s1, string s2) {
    int sum = 0;
    if((s1.length() == s2.length()) && (s1.length() > 0)) {
        for(unsigned int i=0; i < s1.length(); i++) {
            if(s1[i] != s2[i]) sum++;
        }
    } else {
        cout << "ERROR lengths of barcodes and index read do not match!" << endl;
        cout << "Target: " << s1 << endl;
        cout << "Index read: " << s2 << endl;
        sum = -1;
    }
    return sum;
}


int find_overlap_pos(std::string seq1, std::string seq2, int adapterlength) {
    //compare sequences starting at a dovetailed overlap defined by adapterlength
    //Note: this assumes untrimmed sequences of equal length and reports an overlap
    //when 90% identity otherwise it checks until overlap is < minoverlap
    std::string s1 = seq1;
    std::string s2 = seq2;
    unsigned int rlen = s1.length();
    //first check for dovetail:
    if( s1.length() != s2.length() ) {
        std::cout << "!!! " << s1.length() << " " << s2.length() << std::endl; 
        std::cout << s1 << "\n" << s2 << "\n";
        return -10000;
    }
    
    if(((int)s1.length() < adapterlength) || ((int)s2.length() < adapterlength) ) return -10000;
    
    for(unsigned int i = rlen; i >= minoverlap; i--) {
        if ((double)(i - strdist(s1.substr(0,i), s2.substr(rlen-i,i)))/(double)(i) >= overlap_t ) {
            //std::cout << "i=" << i << " " << (double)( (i) - strdist(s1.substr(0,i), s2.substr(rlen-i,i)))/(double)(i)  <<  '\n';
            return i;
        }
        //std::cout << "i=" << i << " " << (double)( (i) - strdist(s1.substr(0,i), s2.substr(rlen-i,i)))/(double)(i)  <<  '\n';
    }
    
    
    return -10000;
}

Read *make_consensus(Read *seq1, Read *seq2) {
    //given two sequences of equal length (these are overlaps only) call a consensus
    //set qualities, and return a new sequence
    string new_seq = "";
    string new_qual = "";
    for(unsigned int i = 0; i < seq1->read.length(); ++i ) {
        if( seq1->read[i] == seq2->read[i] ) {
            new_seq += seq1->read[i];
            new_qual += max(seq1->illumina_quality_string[i], seq2->illumina_quality_string[i]);
        } else {
            if(seq1->illumina_quality_string[i] == seq2->illumina_quality_string[i] ) {
                new_seq += seq1->read[i];
                new_qual += seq1->illumina_quality_string[i];
            } else if(seq1->illumina_quality_string[i] > seq2->illumina_quality_string[i] ) {
                new_seq += seq1->read[i];
                new_qual += (char)( i64_flag ? 64 : 33 + (unsigned int)seq1->illumina_quality_string[i] - (unsigned int)seq2->illumina_quality_string[i]);
            } else {
                new_seq += seq2->read[i];
                new_qual += (char)( i64_flag ? 64 : 33 + (unsigned int)seq2->illumina_quality_string[i] - (unsigned int)seq1->illumina_quality_string[i]);
            }
        }
    }
    Read *consensus = new Read();
    consensus->read = new_seq;
    consensus->illumina_quality_string = new_qual;
    return consensus;
}

