#include "flash.h"

//modify this to down-weight difference in low-quality bases
int strdist(string s1, string s2) {
    int sum = 0;
    if((s1.length() == s2.length()) && (s1.length() > 0)) {
        for(unsigned int i=0; i < s1.length(); ++i) {
            if(s1[i] != s2[i] /*&& (s1[i] != 'N' || s1[i] != 'n' || s2[i] != 'N' || s2[i] != 'n')*/) sum++;
        }
    } else {
        cout << "ERROR lengths of barcodes and index read do not match!" << endl;
        cout << "Target " << s1 << endl;
        cout << "Index read: " << s2 << endl;
        sum = -1;
    }
    return sum;
}


int find_overlap_pos(string seq1, string seq2, int adapterlength) {
    //compare sequences starting at a dovetailed overlap defined by adapterlength
    //Note: this assumes untrimmed sequences of equal length and reports an overlap
    //when 90% identity otherwise it checks until overlap is < minoverlap
    string s1 = seq1;
    string s2 = seq2;
    unsigned int rlen = s1.length();
    //first check for dovetail:
    if(s1.length() != s2.length()) return 0;
    //print "checking dovetail"
    for(int i = adapterlength; i>=0; i--) {
        //cout << (double)(rlen - strdist(s1.substr(0,s1.length()-i), s2.substr(i,s2.length())))/(double)rlen << endl;
        if((double)(rlen - strdist(s1.substr(0,s1.length()-i), s2.substr(i,s2.length())))/(double)rlen >= overlap_t ) {
            //found dovetail overlap
            //cout << (double)(rlen - strdist(s1.substr(0,s1.length()-i), s2.substr(i,s2.length())))/(double)rlen << endl;
            return -i;
        }
    }
    return 0;
}

