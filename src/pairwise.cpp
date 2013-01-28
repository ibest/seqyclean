#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string.h>
#include <assert.h>
#include <algorithm>
#include "pairwise.h"
#include <stdio.h>	
#include <malloc.h>
#include <new>

extern long counter; 
			 
void init( int **S, int a, int b, int d)//function to initialize scoring matrix and alignment tracing matrix
{
    //d - gap penalty
	S[ 0 ][ 0 ] =  0;//Scoring matrix
	
	int i=0, j=0;

	for( j = 1; j <= a; j++ )	{
		S[ 0 ][ j ] =  -j * d;
                
        }
        
        for( i = 1; i <= b; i++ ) {
		S[ i ][ 0 ] = -i * d;
	}
        
        i=0, j=0;
        for(i=1;i<=b;i++) {
            for(j=1;j<=a;j++) {
                S[ i ][ j ] = 0;
            }
        }
}
			  
//print score and traceback
void  print_matrix( int ** F, string seq_1, string seq_2 )
{
	int  L1 = seq_1.length();
	int  L2 = seq_2.length();

	cout << "       ";
	for( int j = 0; j < L1; j++ )
	{
		cout << seq_1[ j ] << "   ";
	}
	cout << "\n  ";

	for( int i = 0; i <= L2; i++ )
	{
		if( i > 0 )
		{
			cout << seq_2[ i-1 ] << " ";
		}
		for( int j = 0; j <= L1; j++ )
		{
			cout.width( 3 );
			cout << F[ i ][ j ] << " ";
		}
		cout << endl;
	}
}




void  print_al( string& seq_1_al, string& seq_2_al )
{
	cout << seq_1_al << endl;
	cout << seq_2_al << endl;
}
				  


//Banded alignment method:
/*
 * seq_1 : read
 * seq_2 : sequence to search for
 * seq_1_al : aligned sequence 1 (seq_1)
 * seq_2_al : aligned sequence 2 (seq_2)
 * k : width of the band
 * d : gap penalty
 */
AlignResult banded(string &seq_1, string &seq_2,int k, int d) {
    int  L1 = seq_1.length();
    int  L2 = seq_2.length();
    string seq_1_al = "";
    string seq_2_al = "";
    AlignResult output;
   
    // Making Scoring Matrix
    int size = L2+1;
    
    int **S = NULL;
    try {
        S = new int*[size];
    } catch (bad_alloc& ba) {
        cerr << "bad_alloc caught: " << ba.what() << endl;
    }
    
    try {
        for(int ii = 0; ii <= size-1; ii++) {
                S[ii] = new int[L2+1];//L1];
        }
    } catch(bad_alloc& ba) {
        cerr << "bad_alloc caught: " << ba.what() << endl;
    }
    
    init( S, L1, L2,d);
    
    int i,j;
    
    
    //Actual algorithm:
//    int i, j;
    for(i=1; i<=L2; i++) {
        for(int h=-k;h<=k;h++) {
            j = i+h;
            if((j>=1) && (j<=L2)) {
                int mm = match(seq_1[ i-1 ],seq_2[ j-1 ]);
                S[i][j] = S[i-1][j-1]+ mm ;
                //if((mm == -1) && (i==j))
                //    output.n_mismatches++;
                
                if(InsideBand(i-1,j,k)) {
                    //int tmp1 = S[i-1][j]-1;
                    S[i][j]=get_max2(S[i][j],S[i-1][j]-1);
                }
                if(InsideBand(i,j-1,k)) {
                    //int tmp2 = S[i][j-1] - 1;
                    S[i][j] = get_max2(S[i][j],S[i][j-1] - 1);
                }
            }
        }
    }
    
    
    
    //print_matrix(S, seq_1,seq_2);
    j=L2; 
    i=L1;
    int a=0;
    
    //Main routine
    do {
       a=0;
        
       if(InsideBand(i-1,j,k)) {
         if (S[i][j] <= S[i-1][j]) {
           a=1;
         }
       }
       if(InsideBand(i,j+1-1,k)) {
         if(S[i][j] <= S[i][j-1]) {
           a=2;
         }
       } 
        
       switch(a) {
            case 1 : 
                seq_1_al += seq_1[ i-1]; 
		seq_2_al += '-';
                i--; 
                //output.n_mismatches++;
                break;
            case 2 : 
                seq_1_al += '-'; 
		seq_2_al += seq_2[ j-1 ];
                //output.n_mismatches++;
                j--;
                break;
            default:
                seq_1_al += seq_1[ i-1 ]; 
		seq_2_al += seq_2[ j-1 ];
                i--; 
                j--;
       }
        
    } while ((i>0) && (j>0)); 
    
      
    if(i>0) {
        while(i>0) {
                seq_1_al += seq_1[ i-1 ];
                seq_2_al += '-';
                i--;
        }
    }
    if(j>0) {
        while(j>0) {
                seq_2_al += seq_2[ j-1 ];
                seq_1_al += '-';
                j--;
        }
    }
    
 //   print_matrix(S, seq_1,seq_2);
    
    for( int jj = 0; jj <= L2; jj++ )  
    {
        //for( int kk = 0; kk <= L2; kk++ )  
        //{
        	free(S[jj]);// = NULL;
                //S[jj] = NULL;
        //}
    }
 
    free(S);
    //S = NULL;
   
    
    reverse( seq_1_al.begin(), seq_1_al.end() );
    reverse( seq_2_al.begin(), seq_2_al.end() );
   
    
    
    output.seq_1_al = seq_1_al;
    output.seq_2_al = seq_2_al;
    
    //AlignScores a_scor = CalcScores(seq_1_al,seq_2_al,seq_1_al.length(),0);
    //output.scores = a_scor.scores;
    //output.n_mismatches = a_scor.mismatches;
    
//    output.tid = tid;
//    output.w = w;
    
    
   // delete seq_1_al;//  = NULL;
 //   seq_1_al.erase();
  //  delete seq_2_al; 
 //   seq_2_al.erase();
    //seq_1.erase();
    //seq_1.erase();
    //output.seq_1 = seq_1;
    
    
    return output;
}

int get_max2(int a, int b) {
    if(a>b) {
        return a;
    } else {
        return b;
    }
}

int match(char &a, char &b) {
    const int  ma =  4;   // Match
    const int  mis = -1;//4;   // Mismatch
    int x,y;
    //Match-mismatch matrix
    const int  mm[ 5 ][ 5 ] = {{ ma, mis, mis, mis },
				{ mis, ma, mis, mis },
				{ mis, mis, ma, mis },
				{ mis, mis, mis, ma },
                                { mis, mis, mis, mis } 
                                };
    
    switch( a ) {
	case 'A':  x = 0;  break;
	case 'C':  x = 1;  break;
	case 'G':  x = 2;  break;
	case 'T':  x = 3; break;
        default : x = 4;
    }

    switch( b )	{
	case 'A':  y = 0;  break;
	case 'C':  y = 1;  break;
	case 'G':  y = 2;  break;
	case 'T':  y = 3; break;
        default : y = 4; break;
    }
    
    return mm[x][y];
    
}

//Checking if the matrix value is inside of the band
bool InsideBand(int i, int j, int k) {
    if(((i-j)>=-k) && ((i-j) <=k)) {
        return true;
    } else {
        return false;
    }
}



/*Score calculation*/
/*
 Here is the algorithm:
 1. If mismatch: -1 penalty;
 2. If match: +2;
 dir : direction. 0 - from left to right , 1 - from right to left
 */
AlignScores CalcScores(string &seq_1, string &seq_2, int lim, int dir) {
    
    AlignScores al_scores;
    
    short scores = 0;
    short mismatches = 0;
    
    //cout << seq_1 << endl;
    //cout << seq_2 << endl;
    if(dir == 0) {
        for(int i = 0; i<lim; i++) {
                /*Mismatch*/
                if((seq_1[i] == '-') && (seq_2[i] == '-')) {
                        scores--; mismatches++; continue;
            
                }
                /*Mismatch*/
                if((seq_1[i] == '-')&& (seq_2[i] != '-')) {
                        scores--; mismatches++; continue;
            
                }
                /*Mismatch*/
                if((seq_1[i] != '-') && (seq_2[i] == '-')) {
                        scores--; mismatches++; continue;
                }
                /*Match*/
                if(seq_1[i] == seq_2[i]) {
                        scores = scores + 2; continue;
                }
                /*Mismatch*/
                if((seq_1[i] != '-') && (seq_2[i] != '-') && (seq_1[i] != seq_2[i])) {
                        scores--; mismatches++; continue;
                }
        
        }
    }
    
    if(dir == 1) {
        for(int i = (int)seq_1.length(); i>(int)seq_1.length() - lim; i--) {
                /*Mismatch*/
                if((seq_1[i] == '-') && (seq_2[i] == '-')) {
                        scores--; mismatches++; continue;
            
                }
                /*Mismatch*/
                if((seq_1[i] == '-')&& (seq_2[i] != '-')) {
                        scores--; mismatches++; continue;
            
                }
                /*Mismatch*/
                if((seq_1[i] != '-') && (seq_2[i] == '-')) {
                        scores--; mismatches++; continue;
                }
                /*Match*/
                if(seq_1[i] == seq_2[i]) {
                        scores = scores + 2; continue;
                }
                /*Mismatch*/
                if((seq_1[i] != '-') && (seq_2[i] != '-') && (seq_1[i] != seq_2[i])) {
                        scores--; mismatches++; continue;
                }
        
        }
    }
    
   
    al_scores.scores = scores;
    al_scores.mismatches = mismatches;
    
  //  seq_1.erase();
  //  seq_2.erase();
    
    return al_scores;
}

/*Calculating mismatches from gaps: */
AlignScores CalcScores2(string &seq_1, int lim, int dir) {
    
    AlignScores al_scores;
    
    short scores = 0;
    short mismatches = 0;
//    char prev_symb = '-';
    
    //cout << seq_1 << endl;
    //cout << seq_2 << endl;
    if(dir == 0) {
        for(int i = 0; i<lim; i++) {
                /*Mismatch*/
                if(seq_1[i] == '-') {
                    scores-=2; 
                    mismatches++; 
                    
                    continue;
                } else {
                    scores++;
                }
        }
    }
    
    if(dir == 1) {
        for(int i = (int)seq_1.length(); i>(int)seq_1.length() - lim; i--) {
                /*Mismatch*/
                if(seq_1[i] == '-') {
                    scores-=2; mismatches++; continue;
                } else {
                    scores++;
                }
        }
    }
    
   
    al_scores.scores = scores;
    al_scores.mismatches = mismatches;
    
    
    return al_scores;
}

/*Calculating mismatches from gaps: */
AlignResult CalcPos(string &seq) {
    
    AlignResult pos;
    
    short pos_left = 0;
    short pos_right = 0;
    
    
    for(int i = 0; i<(int)seq.length(); i++) {
       if(seq[i] == '-') {
           pos_left++;
           continue;
       } else {
          break;
       }
    }
    
    
    
        for(int i = seq.length()-1; i>=0; i--) {
                /*Mismatch*/
                if(seq[i] == '-') {
                    pos_left++;; continue;
                } else {
                    
                }
        }
    
    
   
    pos.pos_left = pos_left;
    pos.pos_right = pos_right;
    
    
    return pos;
}