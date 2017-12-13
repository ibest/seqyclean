#include "flash.h"

//modify this to down-weight difference in low-quality bases
int strdist(std::string s1, std::string s2) {
    int sum = 0;
    if((s1.length() == s2.length()) && (s1.length() > 0)) {
        for(unsigned int i=0; i < s1.length(); i++) {
            if(s1[i] != s2[i]) sum++;
        }
    } else {
        std::cout << "ERROR lengths of barcodes and index read do not match!\n";
        std::cout << "Target: " << s1 << "\n";
        std::cout << "Index read: " << s2 << "\n";
        sum = -1;
    }
    return sum;
}


int find_overlap_pos(std::string seq1, std::string seq2, int minoverlap) {
    //compare sequences starting at a dovetailed overlap defined by adapterlength
    //Note: this assumes untrimmed sequences of equal length and reports an overlap
    //when 75% identity otherwise it checks until overlap is < minoverlap
    std::string s1 = seq1;
    std::string s2 = seq2;
    unsigned int rlen = s1.length();
    
    //first check for dovetail:
    /*if( s1.length() != s2.length() ) {
        std::cout << "!!! " << s1.length() << " " << s2.length() << std::endl; 
        std::cout << s1 << "\n" << s2 << "\n";
        return -10000;
    }*/
    
    //if(((int)s1.length() < adapterlength) || ((int)s2.length() < adapterlength) ) return -10000;
    
    for(unsigned int i = 0; i < rlen-minoverlap; i++) 
    {
        if ((1-(double)(strdist(s1.substr(i,rlen-i), s2.substr(0,rlen-i)))/(double)(rlen-i)) >= overlap_t ) 
        {
            return i;
        }
    }
    
    return -10000;
}

int find_overlap_pos_adapter(std::string seq1, std::string seq2, int adaplen) {
    // adaplen - maximum allowable adapter length
    //compare sequences starting at a dovetailed overlap defined by adapterlength
    //Note: this assumes untrimmed sequences of equal length and reports an overlap
    //when 75% identity otherwise it checks until overlap is < minoverlap
    std::string s1 = seq1;
    std::string s2 = seq2;
    unsigned int rlen = s1.length();
    
    //first check for dovetail:
    /*if( s1.length() != s2.length() ) {
        std::cout << "!!! " << s1.length() << " " << s2.length() << std::endl; 
        std::cout << s1 << "\n" << s2 << "\n";
        return -10000;
    }*/
    
    //if(((int)s1.length() < adapterlength) || ((int)s2.length() < adapterlength) ) return -10000;
    
    for(unsigned int i = 0; i < adaplen+1; i++) 
    {
        if ((1-(double)(strdist(s1.substr(i,rlen-i), s2.substr(0,rlen-i)))/(double)(rlen-i)) >= overlap_t ) 
        {
            return i;
        }
    }
    
    return -10000;
}

Read *make_consensus(Read *seq1, Read *seq2) {
    //given two sequences of equal length (these are overlaps only) call a consensus
    //set qualities, and return a new sequence
    std::string new_seq;
    std::string new_qual;
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
                char q = seq1->illumina_quality_string[i] - seq2->illumina_quality_string[i];
                new_qual += static_cast<char>( i64_flag ? phred_coeff_illumina64 + q : phred_coeff_illumina + q);
            } else {
                new_seq += seq2->read[i];
                char q = seq2->illumina_quality_string[i] - seq1->illumina_quality_string[i];
                new_qual += static_cast<char>( i64_flag ? phred_coeff_illumina64 + q : phred_coeff_illumina + q);
            }
        }
    }
    Read *consensus = new Read();
    consensus->read = new_seq;
    consensus->illumina_quality_string = new_qual;
    return consensus;
}

