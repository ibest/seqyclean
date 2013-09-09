#include "Report.h"

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
            
       fputs(itoa(reads[ff]->initial_length,new char[5],10),  rep_file ); 
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
        if( reads[i]->discarded == 0 ) 
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
            
            
            reads[i]->read = reads[i]->read.substr(0 , reads[i]->rclip );
            string quality = string((char*)reads[i]->quality); 
            quality = quality.substr(0,reads[i]->rclip) ; 
            reads[i]->read = reads[i]->read.substr( reads[i]->lclip, reads[i]->read.length() - reads[i]->lclip );
            quality = quality.substr( reads[i]->lclip, quality.length() - reads[i]->lclip );
            
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

