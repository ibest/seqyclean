#include "Report.h"

void MakeClipPoints()
{
    for(int ff = 0; ff < (int)reads.size(); ff++)
    {
        //Read ID
       if( reads[ff]->rlmid.rmid_name.length() == 0 ) reads[ff]->rlmid.rmid_start = 0;
       if( reads[ff]->b_adapter.length() == 0 ) reads[ff]->b_adapter_pos = 0;
       
       //Clip points
       if( (vector_flag == true) && (qual_trim_flag == true) )
       {
           if(reads[ff]->discarded == 0) 
           {
                int a = reads[ff]->v_start;
                int b = reads[ff]->read.length() - reads[ff]->v_end;
                if( a >= b ) //Vector is on the right side
                {
                        reads[ff]->rclip = min(reads[ff]->rlmid.rmid_start == 0 ? reads[ff]->read.length() : reads[ff]->rlmid.rmid_start, min(reads[ff]->b_adapter_pos == 0 ? reads[ff]->read.length() : reads[ff]->b_adapter_pos, min(reads[ff]->lucy_rclip, reads[ff]->v_start) ) );

                        if( reads[ff]->rclip == reads[ff]->lucy_rclip )
                        {
                                reads[ff]->right_trimmed_by_quality = 1;
                        }
                        else if( reads[ff]->rclip == reads[ff]->rlmid.rmid_start )
                        {
                                reads[ff]->right_trimmed_by_adapter = 1;
                        }
                        else if(reads[ff]->rclip == reads[ff]->v_start)
                        {
                                reads[ff]->right_trimmed_by_vector = 1;
                        }
                        else if(reads[ff]->rclip == reads[ff]->b_adapter_pos)
                        {
                                reads[ff]->right_trimmed_by_adapter = 1;
                        }
                   
                        reads[ff]->lclip = max(reads[ff]->lucy_lclip,reads[ff]->rlmid.lmid_end );
                        if(reads[ff]->lclip == reads[ff]->lucy_lclip)
                        {
                                reads[ff]->left_trimmed_by_quality = 1;
                        }
                        if(reads[ff]->lclip == reads[ff]->rlmid.lmid_end)
                        {
                                reads[ff]->left_trimmed_by_adapter = 1;
                        }
                }
                else //Vector is on the left side
                {
                        reads[ff]->rclip = min(reads[ff]->rlmid.rmid_start == 0 ? reads[ff]->read.length() : reads[ff]->rlmid.rmid_start, min(reads[ff]->b_adapter_pos == 0 ? reads[ff]->read.length() : reads[ff]->b_adapter_pos, reads[ff]->lucy_rclip) );

                        if( reads[ff]->rclip == reads[ff]->lucy_rclip )
                        {
                                reads[ff]->right_trimmed_by_quality = 1;
                        }
                        else if( reads[ff]->rclip == reads[ff]->rlmid.rmid_start )
                        {
                                reads[ff]->right_trimmed_by_adapter = 1;
                        }
                        else if(reads[ff]->rclip == reads[ff]->b_adapter_pos)
                        {
                                reads[ff]->right_trimmed_by_adapter = 1;
                        }  
                
                        reads[ff]->lclip = max3(reads[ff]->lucy_lclip,reads[ff]->rlmid.lmid_end, reads[ff]->v_end );

                        if(reads[ff]->lclip == reads[ff]->lucy_lclip)
                        {
                                reads[ff]->left_trimmed_by_quality = 1;
                        }
                        if(reads[ff]->lclip == reads[ff]->v_end)
                        {
                                reads[ff]->left_trimmed_by_vector = 1;
                        }
                        if(reads[ff]->lclip == reads[ff]->rlmid.lmid_end)
                        {
                                reads[ff]->left_trimmed_by_adapter = 1;
                        }
                }
                        
                                
                if(reads[ff]->rclip <= 1) reads[ff]->discarded = 1; 
                                
            
           } 
           else
           {
               reads[ff]->lclip = reads[ff]->rclip = 1;
           }
       }
       else if( (vector_flag == true) && (qual_trim_flag == false) ) 
       {
           if ( reads[ff]->discarded == 0 )
           {
                int a = reads[ff]->v_start;
                int b = reads[ff]->read.length() - reads[ff]->v_end;
                    
                if( a >= b ) //Vector is on the right side
                {
                    reads[ff]->rclip = min(reads[ff]->rlmid.rmid_start == 0 ? reads[ff]->read.length() : reads[ff]->rlmid.rmid_start, min(reads[ff]->b_adapter_pos == 0 ? reads[ff]->read.length() : reads[ff]->b_adapter_pos, reads[ff]->v_start) );

                    if( reads[ff]->rclip == reads[ff]->rlmid.rmid_start )
                    {
                        reads[ff]->right_trimmed_by_adapter = 1;
                    }
                    else if(reads[ff]->rclip == reads[ff]->v_start)
                    {
                        reads[ff]->right_trimmed_by_vector = 1;
                    }
                    else if(reads[ff]->rclip == reads[ff]->b_adapter_pos)
                    {
                        reads[ff]->right_trimmed_by_adapter = 1;
                    }
                    
                    reads[ff]->lclip = reads[ff]->rlmid.lmid_end;
                    reads[ff]->lclip > 0 ? reads[ff]->left_trimmed_by_adapter = 1 : reads[ff]->left_trimmed_by_adapter = 0;
                    
                }
                else //Vector is on the left side
                {
                    reads[ff]->rclip = min( reads[ff]->read.length(), min( reads[ff]->rlmid.rmid_start == 0 ? reads[ff]->read.length() : reads[ff]->rlmid.rmid_start ,reads[ff]->b_adapter_pos == 0 ? reads[ff]->rlmid.rmid_start : reads[ff]->b_adapter_pos ) );

                    if( reads[ff]->rclip == reads[ff]->rlmid.rmid_start )
                    {
                        reads[ff]->right_trimmed_by_adapter = 1;
                    }
                    else if(reads[ff]->rclip == reads[ff]->b_adapter_pos)
                    {
                        reads[ff]->right_trimmed_by_adapter = 1;
                    }        
                    
                    
                    reads[ff]->lclip = max(reads[ff]->rlmid.lmid_end, reads[ff]->v_end );

                    if(reads[ff]->lclip == reads[ff]->v_end)
                    {
                       reads[ff]->left_trimmed_by_vector = 1;
                    }
                    if(reads[ff]->lclip == reads[ff]->rlmid.lmid_end)
                    {
                       reads[ff]->left_trimmed_by_adapter = 1;
                    }
                }
                        
                if(reads[ff]->rclip <= 1) reads[ff]->discarded = 1; 
                            
           }
           else
           {
               reads[ff]->lclip = reads[ff]->rclip = 1;
           }
       }
       else if( (vector_flag == false) && (qual_trim_flag == true) ) 
       {
           if (reads[ff]->discarded == 0)
           {
                reads[ff]->rclip = min( reads[ff]->rlmid.rmid_start == 0 ? reads[ff]->read.length() : reads[ff]->rlmid.rmid_start, min(reads[ff]->b_adapter_pos == 0 ? reads[ff]->read.length()  : reads[ff]->b_adapter_pos, reads[ff]->lucy_rclip) );
                
                if( reads[ff]->rclip == reads[ff]->rlmid.rmid_start )
                {
                   reads[ff]->right_trimmed_by_adapter = 1;
                }
                else if(reads[ff]->rclip == reads[ff]->b_adapter_pos)
                {
                   reads[ff]->right_trimmed_by_adapter = 1;
                }
                else
                {
                    reads[ff]->right_trimmed_by_quality = 1;
                }
                    
                
                reads[ff]->lclip = max( reads[ff]->lucy_lclip,reads[ff]->rlmid.lmid_end );
                if(reads[ff]->lclip == reads[ff]->lucy_lclip)
                {
                   reads[ff]->left_trimmed_by_quality = 1;
                }
                else
                {
                   reads[ff]->left_trimmed_by_adapter = 1;
                }
                
                if(reads[ff]->rclip <= 1) reads[ff]->discarded = 1; 
                if(reads[ff]->rclip == (int)reads[ff]->read.length()) reads[ff]->right_trimmed_by_adapter = 0;
                
           }
           else
           {
               reads[ff]->lclip = reads[ff]->rclip = 1;
           }
       }
       else if( (vector_flag == false) && (qual_trim_flag == false) ) 
       {
           if(reads[ff]->discarded == 0)
           {    
               
                reads[ff]->rclip = min( reads[ff]->read.length(), min( reads[ff]->rlmid.rmid_start == 0 ? reads[ff]->read.length() : reads[ff]->rlmid.rmid_start ,reads[ff]->b_adapter_pos == 0 ? reads[ff]->read.length() : reads[ff]->b_adapter_pos ) );
                
                if( reads[ff]->rclip == reads[ff]->rlmid.rmid_start )
                {
                   reads[ff]->right_trimmed_by_adapter = 1;
                }
                else if(reads[ff]->rclip == reads[ff]->b_adapter_pos)
                {
                   reads[ff]->right_trimmed_by_adapter = 1;
                }
                
                
                reads[ff]->lclip = reads[ff]->rlmid.lmid_end;
                reads[ff]->lclip > 0 ? reads[ff]->left_trimmed_by_adapter = 1 : reads[ff]->left_trimmed_by_adapter = 0;
                
                if(reads[ff]->rclip <= 1) reads[ff]->discarded = 1; 
                if(reads[ff]->rclip == (int)reads[ff]->read.length()) reads[ff]->right_trimmed_by_adapter = 0;
           } 
           else
           {
               reads[ff]->lclip = reads[ff]->rclip = 1;
           }
           
       }
    }
}

void MakeReport(string rep_file_name) {
    FILE* rep_file;
    
    rep_file = fopen( rep_file_name.c_str(), "w");
    
    //Making a Table :
    fputs("ReadID\tlclip\trclip\tRaw_read_length\tLMHamming\tLMStart\tLMEnd\tLMErr\tLMId\tVecStart\tVecEnd\tVecErr\tVecLen\tVecStartPos\tVecID\tCont\tRMHamming\tRMStart\tRMEnd\tRMErr\tRMId\tRAHamming\tRAStart\tRAEnd\tRAErr\tDiscarded\tLUCY_lclip\tLUCY_rclip\tRoche_left_clip\tRoche_right_clip\n", rep_file );
       
    for(int ff = 0; ff < (int)reads.size(); ff++)
    {
                
       //reads[ff]->readID[0] == '@' ?  fputs((reads[ff]->readID.substr(1, reads[ff]->readID.length() - 1 )).c_str(), rep_file ) : fputs(reads[ff]->readID.c_str(), rep_file );
       fputs(reads[ff]->readID, rep_file );
       fputs("\t", rep_file );
                
       fputs( reads[ff]->lclip == 0 ? "1" : itoa(reads[ff]->lclip,new char[5],10),  rep_file ); //Left clip point
       fputs("\t", rep_file );
       fputs( reads[ff]->rclip == 0 ? "1" : itoa(reads[ff]->rclip,new char[5],10),  rep_file ); //Right clip point
       fputs("\t", rep_file );        
            
       fputs(itoa(reads[ff]->initial_length,new char[5],10),  rep_file ); /*Raw read length*/
       fputs("\t", rep_file );   
            
            //Left midtag
            fputs( reads[ff]->rlmid.lmid_err != 0 ? itoa(reads[ff]->rlmid.lmid_err,new char[5],10) : "NA",  rep_file ); //LMHamming
            fputs("\t", rep_file );
            fputs( reads[ff]->rlmid.lmid_start != 0 ? itoa(reads[ff]->rlmid.lmid_start,new char[5],10) : "NA",  rep_file ); //LMStart
            fputs("\t", rep_file );
            fputs( reads[ff]->rlmid.lmid_end != 0 ? itoa(reads[ff]->rlmid.lmid_end,new char[5],10) : "NA",  rep_file ); //LMEnd
            fputs("\t", rep_file );
            fputs( reads[ff]->rlmid.lmid_err != 0 ? itoa(reads[ff]->rlmid.lmid_err,new char[5],10) : "NA",  rep_file ); //LMErr
            fputs("\t", rep_file );
            fputs(reads[ff]->rlmid.lmid_name == "" ? "NA" : reads[ff]->rlmid.lmid_name.c_str(),  rep_file );
            fputs("\t", rep_file );
            
            //Vector
            if( vector_flag == true )
            {
               fputs(itoa(reads[ff]->r_vec_start,new char[5],10),  rep_file ); //VecStart
               fputs("\t", rep_file );
               fputs(itoa(reads[ff]->r_vec_end,new char[5],10),  rep_file ); //VecEnd
               fputs("\t", rep_file );
               fputs( itoa(reads[ff]->vec_err,new char[5],10) ,  rep_file ); //VecErr
               fputs("\t", rep_file );
               if(reads[ff]->r_vec_start != 0 && reads[ff]->r_vec_end != 0) 
               {
                  fputs( itoa( abs( reads[ff]->r_vec_end - reads[ff]->r_vec_start ),new char[5],10),  rep_file  ); //VecEnd
               } else if( reads[ff]->r_vec_start != 0 ) 
               {
                  fputs( itoa( abs( reads[ff]->clear_length - reads[ff]->r_vec_start ),new char[5],10),  rep_file ); //VecLen
               }
               fputs("\t", rep_file );
               fputs( itoa( reads[ff]->vec_start_pos,new char[5],10),  rep_file  ); //VecStartPos
               fputs("\t", rep_file );
               fputs( reads[ff]->vec_id ,  rep_file  );
               fputs("\t", rep_file );
            } 
            else 
            {
               fputs("NA", rep_file ); //VecStart
               fputs("\t", rep_file ); 
               fputs("NA", rep_file ); //VecEnd
               fputs("\t", rep_file ); 
               fputs("NA", rep_file );//VecErr
               fputs("\t", rep_file );
               fputs("NA", rep_file );//VecLen
               fputs("\t", rep_file );
               fputs("NA", rep_file );//VecStartPos
               fputs("\t", rep_file );
               fputs("NA", rep_file );//VecID
               fputs("\t", rep_file );
            }
            
            //Contaminants
            fputs(itoa(reads[ff]->discarded_by_contaminant, new char[5],10),  rep_file ); //ContID
            fputs("\t", rep_file );
            
            //Right mid tag
            fputs( reads[ff]->rlmid.rmid_err != 0 ? itoa(reads[ff]->rlmid.rmid_err,new char[5],10) : "NA",  rep_file ); //RMHamming
            fputs("\t", rep_file );
            fputs( reads[ff]->rlmid.rmid_start != 0 ? itoa(reads[ff]->rlmid.rmid_start,new char[5],10) : "NA",  rep_file ); //RMStart
            fputs("\t", rep_file );
            fputs( reads[ff]->rlmid.rmid_end != 0 ? itoa(reads[ff]->rlmid.rmid_end,new char[5],10) : "NA",  rep_file ); //RMEnd
            fputs("\t", rep_file );
            fputs( reads[ff]->rlmid.rmid_err != 0 ? itoa(reads[ff]->rlmid.rmid_err,new char[5],10) : "NA",  rep_file ); //RMErr
            fputs("\t", rep_file );  
            fputs( reads[ff]->rlmid.rmid_name == "" ? "NA" : reads[ff]->rlmid.rmid_name.c_str(),  rep_file ); //RMId
            fputs("\t", rep_file );
            
            //Adapter B (also called a barcode)
            if( reads[ff]->b_adapter.length() == 0 ) 
            {  
                fputs("NA",  rep_file ); //AdpBHamming
                fputs("\t", rep_file );
                fputs("NA",  rep_file ); //AdpBStart
                fputs("\t",  rep_file ); 
                fputs("NA", rep_file ); //AdpBEnd
                fputs("\t", rep_file );
                fputs("NA",  rep_file ); //AdpBErr
                fputs("\t", rep_file );
            } 
            else 
            {
                fputs("NA",  rep_file ); //AdpBHamming
                fputs("\t", rep_file );
                fputs(itoa(reads[ff]->b_adapter_pos,new char[5],10),  rep_file ); //AdpBStart
                fputs("\t",  rep_file ); 
                fputs(itoa(reads[ff]->clear_length,new char[5],10), rep_file ); //AdpBEnd
                fputs("\t", rep_file );
                fputs(itoa(reads[ff]->b_adapter_err,new char[5],10),  rep_file ); //AdpBErr
                fputs("\t", rep_file );
            }
                
                
            fputs(itoa(reads[ff]->discarded,new char[5],10),  rep_file ); //RevPHamming
            fputs("\t", rep_file );   
            
            if( qual_trim_flag == true ) 
            {
                fputs(itoa(reads[ff]->lucy_lclip,new char[5],10),  rep_file ); //lucy_lclip
                fputs("\t", rep_file );
                fputs(itoa(reads[ff]->lucy_rclip,new char[5],10),  rep_file ); //lucy_rclip
                fputs("\t", rep_file );
            } 
            else
            {
                fputs("NA",  rep_file );
                fputs("\t", rep_file ); 
                fputs("NA",  rep_file );
                fputs("\t", rep_file ); 
            }
            
            //Roche clip points
            fputs(itoa(reads[ff]->roche_left_clip,new char[5],10),  rep_file ); //roche_left_clip
            fputs("\t", rep_file );
            fputs(itoa(reads[ff]->roche_right_clip,new char[5],10),  rep_file ); //roche_right_clip
                
            
            fputc( '\n', rep_file );
                
        }
    
    
    fclose(rep_file);
    
}


void WriteToFASTQ(string file_name) {
    FILE* output_file;
    
    output_file = fopen( file_name.c_str(), "w");
    
    for(int i=0; i<(int)reads.size(); i++) 
    {
        if( (reads[i]->discarded == 0) && (reads[i]->rclip > 1) && (reads[i]->lclip > 1) ) 
        {
            if(reads[i]->lclip >= (int)reads[i]->read.length()  ) 
            {
                reads[i]->discarded = 1;
                reads[i]->discarded_by_read_length = 1;
                reads[i]->lclip = reads[i]->rclip = 1;
                continue;
            }
            
            string read_id_to_write = string(reads[i]->readID);
            if(read_id_to_write[0] != '@')
                read_id_to_write = "@" + read_id_to_write;
            
            if( reads[i]->lclip >= reads[i]->rclip ) {reads[i]->discarded = 1; reads[i]->discarded_by_read_length = 1; reads[i]->lclip = reads[i]->rclip = 1; continue;}
            if( reads[i]->lclip >= (int)reads[i]->read.length() ) {reads[i]->discarded = 1; reads[i]->discarded_by_read_length = 1; reads[i]->lclip = reads[i]->rclip = 1; continue;}
            if( reads[i]->rclip > (int)reads[i]->read.length() ) {reads[i]->rclip = reads[i]->discarded_by_read_length = 1; reads[i]->lclip = reads[i]->rclip = 1; continue;}
            
            reads[i]->read = reads[i]->read.substr(0 , reads[i]->rclip );
            string quality = string((char*)reads[i]->quality); 
            quality = quality.substr(0,reads[i]->rclip) ; 
            reads[i]->read = reads[i]->read.substr( reads[i]->lclip-1, reads[i]->read.length() - reads[i]->lclip-1 );
            quality = quality.substr( reads[i]->lclip-1, quality.length() - reads[i]->lclip-1 );
            
            if( (int)reads[i]->read.length() < minimum_read_length ) {reads[i]->discarded = 1; reads[i]->discarded_by_read_length = 1; reads[i]->lclip = reads[i]->rclip = 1; continue;}
            
            fputs( read_id_to_write.c_str(), output_file );
            fputc( '\n', output_file );
            fputs( reads[i]->read.c_str(), output_file );
            fputc( '\n', output_file );
            fputc( '+', output_file );
            fputc( '\n', output_file );
            fputs( quality.c_str(), output_file );
            fputc( '\n', output_file );
        } 
    }
    
    fclose(output_file);
}

void WriteToSFF(string file_name) {
    
    for(int i=0; i<(int)reads.size(); i++) 
    {
        if(reads[i]->discarded == 0) 
        {
            if(reads[i]->lclip >= (int)reads[i]->read.length()  ) {
                reads[i]->discarded = 1;
                reads[i]->discarded_by_read_length = 1;
                reads[i]->lclip = reads[i]->rclip = 1;
                continue;
            }
            
        /*    if(reads[i]->rclip > (int)reads[i]->read.length()  ) {
                reads[i]->discarded = 1;
                reads[i]->discarded_by_read_length = 1;
                reads[i]->lclip = reads[i]->rclip = 1;
                continue;
            }
        **/    
            if( reads[i]->lclip >= reads[i]->rclip  ) {
                reads[i]->discarded = 1;
                reads[i]->discarded_by_read_length = 1;
                reads[i]->lclip = reads[i]->rclip = 1;
                continue;
            }
            
            if( (reads[i]->rclip - reads[i]->lclip) < minimum_read_length ) {
                reads[i]->discarded = 1;
                reads[i]->discarded_by_read_length = 1;
                reads[i]->lclip = reads[i]->rclip = 1;
                continue;
            }
        }
        else
        {
            reads[i]->lclip = reads[i]->rclip = 1;
        }
    }
    
    process_fastq_to_sff( (char*)file_name.c_str() );
}


void MakeLucyReport2(char *filename, vector<Read*>& reads) {
    
    
    /*Making a report*/
    FILE* rep_file;
    rep_file = fopen( filename, "w");
    fputs("ReadID\tlucy_lclip\tlucy_rclip\tpoly_T_clip\tpoly_A_clip\tdiscarded\n", rep_file );
    
    for(int i=0; i<(int)reads.size(); i++) 
    {
        /*if(reads[i]->readID[0] == '@') {
                fputs((reads[i]->readID.substr(1, reads[i]->readID.length() - 1 )).c_str(), rep_file );
        } else {
                fputs(reads[i]->readID.c_str(), rep_file );
        }**/
        fputs(reads[i]->readID, rep_file );
        
        fputs("\t", rep_file );
        fputs(itoa(reads[i]->lucy_lclip == 0 ? 1 : reads[i]->lucy_lclip, new char[5],10),  rep_file );
        fputs("\t", rep_file );
        fputs(itoa(reads[i]->lucy_rclip == 0 ? 1 : reads[i]->lucy_rclip, new char[5],10),  rep_file );
        fputs("\t", rep_file );
        fputs(itoa(reads[i]->poly_T_clip == 0 ? 1 : reads[i]->poly_T_clip, new char[5],10),  rep_file );
        fputs("\t", rep_file );
        fputs(itoa(reads[i]->poly_A_clip == 0 ? 1 : reads[i]->poly_A_clip, new char[5],10),  rep_file );
        fputs("\t", rep_file );
        
        fputs(itoa(reads[i]->discarded,new char[5],10),  rep_file );
        
        fputc( '\n', rep_file );
        
    }
    
    fclose(rep_file);
    
}



void MakeLucyFastq(string custom_file_name) {
    FILE* output_file;
    
    output_file =  fopen( (char*)custom_file_name.c_str(), "w" );
    
    for(int i=0; i<(int)reads.size(); i++) {
        if(reads[i]->discarded == 0) {
            if(reads[i]->lclip >= (int)reads[i]->read.length()  ) {
                reads[i]->discarded = 1;
                reads[i]->discarded_by_read_length = 1;
                continue;
            }
            
            string read_id_to_write = reads[i]->readID;
            if(read_id_to_write[0] != '@')
                read_id_to_write = "@" + read_id_to_write;
            
            reads[i]->lclip = reads[i]->lucy_lclip;
            reads[i]->rclip = reads[i]->lucy_rclip;
            
            if( reads[i]->lclip >= reads[i]->rclip ) 
            {
                reads[i]->discarded = 1; 
                reads[i]->discarded_by_read_length = 1;
                continue;
            }
            if( reads[i]->lclip >= (int)reads[i]->read.length() ) 
            {
                reads[i]->discarded = 1;
                reads[i]->discarded_by_read_length = 1;
                continue;
            }
            if( reads[i]->rclip > (int)reads[i]->read.length() ) 
            {
                reads[i]->rclip = (int)reads[i]->read.length();
            }
            
            
            reads[i]->read = reads[i]->read.substr(0 , reads[i]->rclip );
            string quality = string((char*)reads[i]->quality);
            quality = quality.substr(0,reads[i]->rclip) ; 
            reads[i]->read = reads[i]->read.substr( reads[i]->lclip-1, reads[i]->read.length() - reads[i]->lclip-1 );
            quality = quality.substr( reads[i]->lclip-1, quality.length() - reads[i]->lclip-1 );
            
            if( reads[i]->read.length() < 50 ) 
            {
                reads[i]->discarded = 1; 
                reads[i]->discarded_by_read_length = 1;
                continue;
            }
            
            if( reads[i]->rclip < reads[i]->initial_length ) 
            {
                reads[i]->right_trimmed_by_quality = 1;
            }
            if( reads[i]->lclip > 1 ) 
            {
                reads[i]->left_trimmed_by_quality = 1;
            }
            
            fputs( read_id_to_write.c_str(), output_file );
            fputc( '\n', output_file );
            fputs( reads[i]->read.c_str(), output_file );
            fputc( '\n', output_file );
            fputc( '+', output_file );
            fputc( '\n', output_file );
            fputs( quality.c_str(), output_file );
            fputc( '\n', output_file );
        } 
    }
    
    fclose(output_file);
}

void MakeFinalStatistics( fstream &sum_stat )
{
    long bases_anal = 0;
    long num_vectors = 0;
    long num_contaminants = 0;
    long discarded = 0;
    long accepted = 0;
    long discarded_by_contaminant = 0;
    long discarded_by_read_length = 0;
    /*Left trims*/
    long left_trimmed_by_quality = 0;
    long left_trimmed_by_adapter = 0;
    long left_trimmed_by_vector = 0;
    /*Right trims/discards*/
    long right_trimmed_by_quality = 0;
    long right_trimmed_by_adapter = 0;
    long right_trimmed_by_vector = 0;
    
    long left_mid_tag, right_mid_tag;
    left_mid_tag = right_mid_tag = 0;
    
    double avg_trim_len, avg_read_len, cnt_avg, cnt_avg_len, avg_left_trim_len, avg_right_trim_len, cnt_avg_left_trim_len, cnt_avg_right_trim_len ;
    avg_trim_len = avg_read_len = cnt_avg = cnt_avg_len = avg_left_trim_len = avg_right_trim_len = cnt_avg_left_trim_len = cnt_avg_right_trim_len = 0;
    
    for(int i=0; i<(int)reads.size(); i++)
    {
        bases_anal+=reads[i]->initial_length;
        
       if(reads[i]->rlmid.lmid_start != 0 ) left_mid_tag+=1;
       if( (reads[i]->rlmid.rmid_start != 0) && (reads[i]->rlmid.rmid_start < (int)reads[i]->read.length()) ) right_mid_tag+=1;
       
       if(reads[i]->vector_found == 1) num_vectors+=1;
       if(reads[i]->contaminants == 1) num_contaminants+=1;
        
       if(reads[i]->discarded == 1) discarded+=1;
       if(reads[i]->discarded == 0) accepted+=1;
       
       if(reads[i]->discarded_by_contaminant == 1) discarded_by_contaminant+=1;
       if(reads[i]->discarded_by_read_length == 1) discarded_by_read_length+=1;
                
       if(reads[i]->left_trimmed_by_quality == 1) left_trimmed_by_quality+=1;
       if(reads[i]->left_trimmed_by_adapter == 1) left_trimmed_by_adapter+=1;
       if(reads[i]->left_trimmed_by_vector == 1) left_trimmed_by_vector+=1;
                
       if(reads[i]->right_trimmed_by_quality == 1) right_trimmed_by_quality+=1;
       if(reads[i]->right_trimmed_by_adapter == 1) right_trimmed_by_adapter+=1;
       if(reads[i]->right_trimmed_by_vector == 1) right_trimmed_by_vector+=1;
       
       if( reads[i]->discarded == 0 )
       {
           //cnt_avg_len+=1;
           //avg_read_len = GetAvg( avg_read_len, cnt_avg_len, reads[i]->read.length() );
           
           //Average right and left trimmed lengths
           if( reads[i]->initial_length > reads[i]->rclip )
           {
               cnt_avg+=1;
               avg_trim_len = GetAvg( avg_trim_len, cnt_avg, reads[i]->rclip - reads[i]->lclip );
               
               cnt_avg_right_trim_len+=1;
               avg_right_trim_len = GetAvg( avg_right_trim_len, cnt_avg_right_trim_len, reads[i]->initial_length - reads[i]->rclip );
           }
           if( reads[i]->lclip > 0 )
           {
               cnt_avg_left_trim_len+=1;
               avg_left_trim_len = GetAvg( avg_left_trim_len, cnt_avg_left_trim_len, reads[i]->lclip );
           }
           
           
       }
       
       
    }
    
    /*
    stat_str = "====================Summary Statistics====================\n" +
                ("Reads analyzed: " +  i2str(reads.size(),new char[15],10)  + ", Bases:" +  i2str(bases_anal, new char[25],10)  + "\n") +
                "Found ->\n" +
                "Left mid tag: : " + i2str(left_mid_tag, new char[15], 10) + ", " + str2double( (double)left_mid_tag/(double)reads.size()*100.0, 10) + "%\n" + 
                        "Right mid tag: : " + i2str(right_mid_tag, new char[15], 10) + ", " + str2double( (double)right_mid_tag/(double)reads.size()*100.0, 10) + "%\n" + 
                        ( vector_flag ? "# of reads with vector: " + i2str(num_vectors,new char[15],10) + ", " + str2double( (double)num_vectors/(double)reads.size()*100.0, 10) + "%\n" : "") +
                        ( contaminants_flag ? "# of reads with contaminants: " + i2str(num_contaminants,new char[15],10) + ", " + str2double( (double)num_contaminants/(double)reads.size()*100.0, 10) + "%\n" : "") +
                        "Reads left trimmed ->\n" +
                        ( qual_trim_flag ? "By quality: " +  i2str(left_trimmed_by_quality,new char[15], 10) + "\n" : "" ) +
                        ( vector_flag ? "By vector: " +  i2str(left_trimmed_by_vector,new char[15], 10) + "\n" : "" ) +
                        "Average left trimmed length: " + str2double(avg_left_trim_len, 10) + " bp\n" +
                        "Reads right trimmed ->\n" +
                        "By adapter: " +  i2str(right_trimmed_by_adapter,new char[15],10) + "\n" +
                        ( qual_trim_flag ? "By quality: " +  i2str(right_trimmed_by_quality,new char[15],10) + "\n" : "") +
                        ( vector_flag ? "By vector: " +  i2str(right_trimmed_by_vector,new char[15],10) + "\n" : "" ) +
                        "Average right trimmed length: " + str2double(avg_right_trim_len,10) + " bp\n" +
                        "Reads discarded: " + i2str(discarded,new char[15],10) + "\n" +
                        ( contaminants_flag ? "By contaminants: " +  i2str(discarded_by_contaminant,new char[15],10) + "\n" : "" ) +
                        "By read length: " +  i2str(discarded_by_read_length,new char[15],10) + "\n" +
                        "-----------------------------------------------------------\n" +
                        "----------------------Summary for Roche 454----------------------\n" +
                        ("Reads accepted: " + i2str(accepted,new char[15],10) + ", " + str2double( (double)accepted/(double)reads.size()*100) + "%\n") +
                        ("Average trimmed length: " + str2double(avg_trim_len,10) + " bp\n") +
                        ("Average read length: " + str2double(avg_read_len,10) + " bp\n");
    
    cout << stat_str;
    sum_stat << stat_str;
    **/
    
    
    cout << "====================Summary Statistics====================\n";
    sum_stat << "====================Summary Statistics====================\n";
    cout << "Reads analyzed: " << reads.size() << ", Bases:" << bases_anal << "\n";
    sum_stat << "Reads analyzed: " << reads.size() << ", Bases:" << bases_anal << "\n";
    
    cout << "Found ->\n";
    sum_stat << "Found ->\n";
    
    cout << "Left mid tag: " << left_mid_tag << ", " << ( (double)left_mid_tag/(double)reads.size()*100.0) << "%\n";
    sum_stat << "Left mid tag: " << left_mid_tag << ", " << ( (double)left_mid_tag/(double)reads.size()*100.0) << "%\n";
    
    cout << "Right mid tag: " << right_mid_tag << ", " << ( (double)right_mid_tag/(double)reads.size()*100.0) << "%\n";
    sum_stat << "Right mid tag: " << right_mid_tag << ", " << ( (double)right_mid_tag/(double)reads.size()*100.0) << "%\n";
    
    if(vector_flag)
    {
        cout << "# of reads with vector: " << num_vectors << ", " << ( (double)num_vectors/(double)reads.size()*100.0) << "%\n";
        sum_stat << "# of reads with vector: " << num_vectors << ", " << ( (double)num_vectors/(double)reads.size()*100.0) << "%\n";
    }
    
    if(contaminants_flag)
    {
        cout << "# of reads with contaminants: " << num_contaminants << ", " << ( (double)num_contaminants/(double)reads.size()*100.0) << "%\n";
        sum_stat << "# of reads with contaminants: " << num_contaminants << ", " << ( (double)num_contaminants/(double)reads.size()*100.0) << "%\n";
    }
    
    cout << "Reads left trimmed ->" << "\n";
    sum_stat << "Reads left trimmed ->" << "\n";
    
    cout << "By adapter: " <<  left_trimmed_by_adapter << "\n";
    sum_stat << "By adapter: " <<  left_trimmed_by_adapter << "\n";
    
    
    if(qual_trim_flag)
    {
        cout << "By quality: " <<  left_trimmed_by_quality << "\n";
        sum_stat << "By quality: " <<  left_trimmed_by_quality << "\n";
    }
    
    if(vector_flag)
    {
        cout << "By vector: " <<  left_trimmed_by_vector << "\n";
        sum_stat << "By vector: " <<  left_trimmed_by_vector << "\n";
    }
    
    cout << "Average left trim length: " << avg_left_trim_len << " bp\n";
    sum_stat << "Average left trim length: " << avg_left_trim_len << " bp\n";
    
    cout << "Reads right trimmed ->" << "\n";
    sum_stat << "Reads right trimmed ->" << "\n";
    
    cout << "By adapter: " <<  right_trimmed_by_adapter << "\n";
    sum_stat << "By adapter: " <<  right_trimmed_by_adapter << "\n";
    
    if(qual_trim_flag)
    {
        cout << "By quality: " <<  right_trimmed_by_quality << "\n";
        sum_stat << "By quality: " <<  right_trimmed_by_quality << "\n";
    }
    
    if(vector_flag)
    {
        cout << "By vector: " <<  right_trimmed_by_vector << "\n";
        sum_stat << "By vector: " <<  right_trimmed_by_vector << "\n";
    }
    
    cout << "Average right trim length: " << avg_right_trim_len << " bp\n";
    sum_stat << "Average right trim length: " << avg_right_trim_len << " bp\n";
    
    cout << "Reads discarded: " << discarded << " ->\n";
    sum_stat << "Reads discarded: " << discarded << " ->\n";
    
    if(contaminants_flag)
    {
        cout << "By contaminants: " <<  discarded_by_contaminant << "\n";
        sum_stat << "By contaminant: " <<  discarded_by_contaminant << "\n";
    }
    
    cout << "By read length: " <<  discarded_by_read_length << "\n";
    sum_stat << "By read length: " <<  discarded_by_read_length << "\n";
    
    cout << "--------------------------------------------------------\n";
    sum_stat << "--------------------------------------------------------\n";
    
    cout << "Reads accepted: " << accepted << ", %" << (double)accepted/(double)reads.size()*100 << "\n";
    sum_stat << "Reads accepted: " << accepted << ", %" << (double)accepted/(double)reads.size()*100 << "\n";
    
    cout << "Average trimmed length: " << avg_trim_len << " bp\n";
    sum_stat << "Average trimmed length: " << avg_trim_len << " bp\n";
    
    cout << "==========================================================\n";
    sum_stat << "==========================================================\n";
}


void MakeRocheReport(char* file_name, char* roche_fname) {
    //Creates a report from a given XML file obtained from SFFFile
    /*Report file*/
    FILE* rep_file;
    
    //Roche clip points
    cout << "Roche clip points:\n";
    
    int roche_lclip=0;
    int roche_rclip=0;
        
    printf("Parsing XML file...\n");
    int ii=0;
    std::string str;
    string rec_id = "";
    
    ifstream infile_Roche;
    
    //Making a Table :
    infile_Roche.open(roche_fname,ifstream::in);
    rep_file = fopen( file_name, "w");
    /*Header*/
    fputs("ReadID\tROCHE_lclip\tROCHE_rclip\n", rep_file );
    while (getline(infile_Roche, str) ) 
    {
       std::remove(str.begin(), str.end(), ' ');
       if(str.substr(0,12)=="<trace_name>" )
       {
           rec_id = str.substr(12,14); 
           ii++;
           continue;        
       }
                
       if(str.substr(0,18) == "<clip_vector_left>") 
       {
            int f1 = str.find("</clip_vector_left>") - 18 ;
                roche_lclip = atoi(str.substr(18,f1).c_str());
                ii++;
                continue;
        }
         
       if(str.substr(0,19) == "<clip_vector_right>") 
       {
            roche_rclip = atoi(str.substr(19,str.find("</clip_vector_right>")-19).c_str());
            
            ii = 0;
            
            fputs( rec_id.c_str(), rep_file );
            fputs("\t",rep_file);
            /*Left clip point*/
            fputs( itoa(roche_lclip, new char[5], 10), rep_file );
            fputs("\t",rep_file);
            /*Right clip point*/
            fputs( itoa(roche_rclip, new char[5], 10), rep_file );
            fputs("\n",rep_file);
            
            rec_id.clear();
            roche_lclip=0;
            roche_rclip=0;
            
            continue;
        }
    }
    fclose(rep_file);
    infile_Roche.close();
    
}

