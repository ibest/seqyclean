#include "flash.h"

//modify this to down-weight difference in low-quality bases
int strdist(string s1, string s2) {
    int sum = 0;
    if((s1.length() == s2.length()) && (s1.length() > 0)) {
        for(unsigned int i=0; i < s1.length(); ++i) {
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


int find_overlap_pos(string seq1, string seq2, int adapterlength, bool flag) {
    //compare sequences starting at a dovetailed overlap defined by adapterlength
    //Note: this assumes untrimmed sequences of equal length and reports an overlap
    //when 90% identity otherwise it checks until overlap is < minoverlap
    string s1 = seq1;
    string s2 = seq2;
    unsigned int rlen = s1.length();
    //first check for dovetail:
    if( s1.length() != s2.length() ) return -10000;
    if((int)s1.length() < adapterlength) return -10000;
    if((int)s2.length() < adapterlength) return -10000;
    
    if(!flag) {
        //print "checking dovetail"
        for(int i = adapterlength; i>=0; i--) {
            if((double)((rlen - i) - strdist(s1.substr(0,rlen-i),s2.substr(i,rlen-i)))/(double)(rlen-i) >= overlap_t ) {
               //found dovetail overlap
               return -i;
            }
        }
    }
    if(flag) {
        //Next check for perfect overlap
        if((double)(rlen - strdist(s1, s2))/(double)rlen >= overlap_t ) {
                //found perfect overlap
                return 0;
        }
        //Finally check for partial overlap
        for(unsigned int i = 1; i <= rlen-minoverlap; i++) {
            if((double)( (rlen-i) - strdist(s1.substr(i,rlen-i), s2.substr(0,rlen-i)))/(double)(rlen-i) >= overlap_t ) {
                return i;
            }
            
        }
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
                new_qual += (char)( 33 + (unsigned int)seq1->illumina_quality_string[i] - (unsigned int)seq2->illumina_quality_string[i]);
            } else {
                new_seq += seq2->read[i];
                new_qual += (char)( 33 + (unsigned int)seq2->illumina_quality_string[i] - (unsigned int)seq1->illumina_quality_string[i]);
            }
        }
    }
    Read *consensus = new Read();
    consensus->read = new_seq;
    consensus->illumina_quality_string = new_qual;
    return consensus;
}

