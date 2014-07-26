/* 
 * File:   Read.h
 * Author: ilya
 *
 * Created on 3 Декабрь 2012 г., 10:12
 */

#ifndef READ_H
#define	READ_H

#include <stdio.h>
#include <iostream>
#include <string>
#include <stdint.h>


using namespace std;

typedef struct {
    /*Left MID*/
    string lmid_name;
    string lmid_value;
    short lmid_id;
     int lmid_start;
     int lmid_end;
    short lmid_err;
    /*Right MID*/
    string rmid_name;
    string rmid_value;
    short rmid_id;
     int rmid_start;
     int rmid_end;
    short rmid_err;
} RL_MID;

class Read {
public:
   Read();
   Read(const Read& orig);
   virtual ~Read();
    
   char* readID;
   string read;
   string illumina_readID;
   
   /*----Left clip point----*/
   unsigned short lclip;
   short lclip_errors;
   int L_lclip;
   /*----Right clip point---*/
   short rclip_errors;
   short rclip;
   int R_rclip;
   /*----------------------*/
   short contaminants;
   int initial_length;
   int clear_length;
   uint8_t* quality;
   bool clip_found;
   RL_MID rlmid;
   /*----Lucy clips---*/
   int lucy_lclip;
   int lucy_rclip;
   /*SFF parameters*/
   uint16_t *flowgram;  
   uint8_t  *flow_index; /* relative to last */
   /*Primer*/
   string fwd_bc;
   string lbtype;
   int fwdP_start;
   int fwdP_end;
   short fwdP_errors;
   string rev_bc;
   string rbtype;
   int revP_start;
   int revP_end;
   short revP_errors;
   /*Adapter B*/
   string b_adapter;
    int b_adapter_pos;
   short b_adapter_err;
   /*Vector*/
   int r_vec_start;
   int r_vec_end;
   int l_vec_start;
   int l_vec_end;
   int vec_err;
   int vec_len;
   int vec_start_pos;
   char *vec_id;
   int v_start;
   int v_end;
   int vector_found;
   /*Contaminants*/
   string cont_id;
   int contam_found;
   
   
   
   /*Illumina parameters*/
   int tru_sec_pos; /*Coordinates of TruSec adapter*/
   int index_i5; /*index of TrueSec adapter i5*/
   int index_i7; /*index of TrueSex adapter i7*/
   bool fwd_adapter; /*does the read contain a forward or reverse complement of the adapter*/
   bool partial_adapter; /*does read contain partial or full TrueSec adaptor*/
   
   
   int tru_sec_found;
    
   
   /*Statistical parameters*/
   short discarded;
   short discarded_by_quality;
   short discarded_by_contaminant;
   short discarded_by_vector;
   short discarded_by_read_length;
   short discarded_by_polyAT;
   /*Left trims*/
   short left_trimmed_by_quality;
   short left_trimmed_by_adapter;
   short left_trimmed_by_vector;
   short left_trimmed_by_polyat;
   /*Right trims/discards*/
   short right_trimmed_by_quality;
   short right_trimmed_by_adapter;
   short right_trimmed_by_vector;
   short right_trimmed_by_polyat;
   
   int poly_T_clip;
   int poly_A_clip;
   bool poly_A_found;
   bool poly_T_found;
   
   //Roche clip points
   int roche_left_clip;
   int roche_right_clip;
   
   string illumina_quality_string;
    
    
private:

};

#endif	/* READ_H */

