/* 
 * File:   Read.cpp
 * Author: ilya
 * 
 * Created on 3 Декабрь 2012 г., 10:12
 */

#include "Read.h"

Read::Read() 
{
    this->illumina_readID = "NA";
    this->contaminants = 0;
                        
    this->initial_length = 0;
    this->lclip = 0;
    this->rclip = 0;
    this->read = "NA";
    this->illumina_quality_string = "NA";
    this->clip_found = false;
                        
                        
    this->lucy_lclip = 0;
    this->lucy_rclip = 0;
    this->r_vec_start = 0;
    this->r_vec_end = 0;
    this->l_vec_start = 0;
    this->l_vec_end = 0;
    this->vec_err = 0;
    this->vec_len = 0;
    this->vec_start_pos = 0;
    this->vec_id = (char*)"NA";
    
    this->rlmid.lmid_start = 0;
    this->rlmid.lmid_end = 0;
    this->rlmid.lmid_id = 0;
    this->rlmid.lmid_err = 0;
    
    this->rlmid.rmid_start = 100000;
    this->rlmid.rmid_end = 100000;
    this->rlmid.rmid_id = 0;
    this->rlmid.rmid_err = 0;
    
    this->b_adapter_pos = 100000;
            
    this->discarded = 0;
    this->discarded_by_quality = 0;
    this->discarded_by_contaminant = 0;
    this->discarded_by_vector = 0;
    this->discarded_by_read_length = 0;
    this->discarded_by_polyAT = 0;
   
    //Left trims
    this->left_trimmed_by_quality = 0;
    this->left_trimmed_by_adapter = 0;
    this->left_trimmed_by_vector = 0;
    this->left_trimmed_by_polyat = 0;
    //Right trims/discards
    this->right_trimmed_by_quality = 0;
    this->right_trimmed_by_adapter = 0;
    this->right_trimmed_by_vector = 0;
    this->right_trimmed_by_polyat = 0;
    
    this->poly_T_clip = 0;
    this->poly_A_clip = 0;
    
    this->roche_left_clip = 0;
    this->roche_right_clip = 0;
    
    this->tru_sec_pos = 0;
    this->tru_sec_found = 0;
    
    this->vector_found = 0;
    this->contam_found = 0;
    
    this->merge = false;
}

Read::Read(const Read& orig) {
}

Read::~Read() {
    
}

