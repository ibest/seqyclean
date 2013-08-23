#include "Roche.h"

void RocheRoutine()
{
    if( !trim_adapters_flag ) 
    {
           
        for(int i=0; i<(int)roche_names.size(); ++i)
        {
        
           if( string(roche_names[i]).substr( strlen(roche_names[i])-5, 5 ) == "fastq" ) 
           { /*FASTQ file given. Process it.*/
              ParseFastqFile(roche_names[i], reads);
           } 
           else if( string(roche_names[i]).substr( strlen(roche_names[i])-3, 3 ) == "sff" ) 
           {
              process_sff_to_fastq( roche_names[i], 0 );
           }
           else
           {
               cout << "Unknown file format\n";
               return;
           }
        }
            
        if( qual_trim_flag  )
        {
           cout << "Only LUCY clipping.\n";
           QualityTrimming(reads);
        }
           
        cout << "Making cleaned output data file...\n";
        MakeLucyFastq( output_prefix + "_qual.fastq" );
           
        string report_filename =  output_prefix + "_qual_Report.tsv"; 
        MakeLucyReport2( (char*)report_filename.c_str(), reads );
        cout << "Done Lucy ONLY quality trimming.\n";
           
        return;
        
    } 
    else 
    {
       //DEBUG - для тестирования SFF конвертации 
       /*for(int i=0; i<(int)roche_names.size(); ++i)
       {
        if(debug_flag ) 
        {
          process_sff_to_fastq(roche_names[i], 0);
          process_fastq_to_sff((char*)"out.sff"); 
          process_sff_to_fastq((char*)"out.sff",  0);
          return;
        }
       }*/
        
       //Building dictionary for RL MIDS : 
       if(custom_rlmids_flag ) 
       {
          Build_RLMIDS_Dictionary(rlmids_file);
       } 
       else 
       {
          Build_RLMIDS_Dictionary();
       }
        
       for(int i=0; i<(int)roche_names.size(); ++i)
       {
            cout << "Parsing file: " << roche_names[i] << "..." << endl;
        
            //If SFF format is given -> process it
            if( string(roche_names[i]).substr( strlen(roche_names[i])-3, 3 ) == "sff" ) 
            {
               cout << "File is in SFF format, starting conversion...\n" ;
               process_sff_to_fastq( roche_names[i], 0 );
               if(output_fastqfile_flag == false)
               {
                 output_sfffile_flag = true;
               }
            } 
            else if(string(roche_names[i]).substr( strlen(roche_names[i])-5, 5 ) == "fastq") 
            {
               //FASTQ file given. Process it.
               cout << "File is in FASTQ format, starting conversion...\n" ;
               ParseFastqFile(roche_names[i], reads);
               output_fastqfile_flag = true;
               output_sfffile_flag = false;
            }
       }
        
       //reads_total += reads.size();
       cout << "Conversion finished. Total number of reads read from given file(s): " << reads.size() << endl;
       
       if(polyat_flag) {
           //If poly A/T flag is set:
           for(unsigned long i=0; i< reads.size(); i++) {
              PolyAT_Trim(reads[i]);
           }
       }
       
       /*If quality trimming flag is set up -> perform the quality trimming before vector/contaminants/adaptors clipping.*/
       if( qual_trim_flag  ) 
       {
           for(unsigned long i=0; i< reads.size(); i++) 
           {
              QualTrim( reads[i], max_a_error, max_e_at_ends );/*This function generates LUCY clips of the read. Later they should be compared and read should be trimmed based on the results of comparison.*/
           }
       }
    
       if(contaminants_flag )
       {
           RemoveContaminants454(reads);
       }
           
       /*Run the main routine: Adaptor + Vector/Contaminants trimming or only Adaptors*/
       cout << "Running the main pipeline..." << endl;
       MainPipeLine();
       cout << "Making clip points..." << endl;
       MakeClipPoints();
       
       cout << "Making output files..." << endl;
       
       WriteToSFF( roche_output_file_name );
        
       if( output_fastqfile_flag ) 
       {
           WriteToFASTQ( output_prefix + ".fastq" );
       }
       
       
       cout << "Making a report..." << endl;
       MakeReport( roche_rep_file_name );
       
       //---End of making statistics----
       MakeFinalStatistics(sum_stat);
        
       reads.clear();
        
       /*Reset counters*/
       accept_counter = 0;
       discard_counter = 0;
       trim_counter = counter = 0;
    }        
       
    reads.clear();
    
}

void ParseFastqFile(char* fastq_file, vector<Read*> &reads) 
{
   std::string line;
    
   printf("Parsing .fastq file...\n");
    
   vector<string> record_block;
   int ii = 0;
        
   igzstream in(fastq_file);
        
   while ( getline(in,line) ) 
   {
       if(ii==0) record_block.push_back(line); /*Read ID*/
                
       if(ii==1) record_block.push_back(line); /*actual data*/
                
       if(ii==2) record_block.push_back(line);/*a symbol "+"*/
                
       if(ii==3) 
       {
           record_block.push_back(line); /*Quality scores*/
           
           Read* read = new Read();
           read->readID = (char*)record_block[0].c_str();
           read->initial_length = record_block[1].length();
           read->read = record_block[1];
           read->quality = (uint8_t*)malloc(sizeof(uint8_t)*record_block[3].length());
           for(unsigned int jj=0; jj<record_block[3].length();++jj)
           {
               read->quality[jj] = record_block[3][jj];
           }
                        
           counter++;
           
           reads.push_back(read);
           
       }
        
       ii++;
       if(ii == 4) 
       {
          ii = 0;
          record_block.clear();
       }
   }
    
   in.close();
     
}

void QualityTrimming( vector<Read*>& reads ) 
{
    for(unsigned int ii=0; ii< reads.size(); ii++) 
    {
        if(reads[ii]->discarded == 0)
        {
           QualTrim( reads[ii], max_a_error, max_e_at_ends );/*This function generates LUCY clips of the read. Later they should be compared and read should be trimmed based on the results of comparison.*/
        }
        
        if( lucy_only_flag  ) 
        {
           reads[ii]->lclip = reads[ii]->lucy_lclip;
           reads[ii]->rclip = reads[ii]->lucy_rclip;
           
        }
    }
          
}

void RemoveContaminants454(vector<Read*>& reads454)
{
    for(unsigned int index = 0; index < (unsigned int)reads454.size(); index++)
    {
        
        if(CheckContaminants(reads454[index]->read) == 0) 
        {
           reads454[index]->discarded = 1;
           reads454[index]->discarded_by_contaminant = 1;
           reads454[index]->contaminants = 1;
           reads454[index]->contam_found = 1;
        }
    }
}


void MakeClipPoints()
{
    for(unsigned int ff = 0; ff < reads.size(); ff++)
    {
       //cout << "badp: " << reads[ff]->b_adapter_pos << endl;
  
       //Read ID
       if( reads[ff]->rlmid.rmid_name.length() == 0 ) reads[ff]->rlmid.rmid_start = 0;
       //if( reads[ff]->b_adapter.length() == 0 ) reads[ff]->b_adapter_pos = 0;
       
       //Clip points
       if( (vector_flag == true) && (qual_trim_flag == true) )
       {
           if(reads[ff]->discarded == 0) 
           {
                int a = reads[ff]->v_start;
                int b = reads[ff]->read.length() - reads[ff]->v_end;
                if( a >= b ) //Vector is on the right side
                {
                    reads[ff]->rclip = min(reads[ff]->rlmid.rmid_start == 0 ? (int)reads[ff]->read.length() : reads[ff]->rlmid.rmid_start, min(reads[ff]->b_adapter_pos == 0 ? (int)reads[ff]->read.length() : reads[ff]->b_adapter_pos, min(reads[ff]->lucy_rclip, reads[ff]->v_start) ) );

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
                        reads[ff]->rclip = min(reads[ff]->rlmid.rmid_start == 0 ? (int)reads[ff]->read.length() : reads[ff]->rlmid.rmid_start, min(reads[ff]->b_adapter_pos == 0 ? (int)reads[ff]->read.length() : reads[ff]->b_adapter_pos, reads[ff]->lucy_rclip) );

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
                
                        reads[ff]->lclip = max(reads[ff]->lucy_lclip,max(reads[ff]->rlmid.lmid_end, reads[ff]->v_end ) );

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
                    reads[ff]->rclip = min(reads[ff]->rlmid.rmid_start == 0 ? (int)reads[ff]->read.length() : reads[ff]->rlmid.rmid_start, min(reads[ff]->b_adapter_pos == 0 ? (int)reads[ff]->read.length() : reads[ff]->b_adapter_pos, reads[ff]->v_start) );

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
                    reads[ff]->rclip = min( (int)reads[ff]->read.length(), min( reads[ff]->rlmid.rmid_start == 0 ? (int)reads[ff]->read.length() : reads[ff]->rlmid.rmid_start ,reads[ff]->b_adapter_pos == 0 ? reads[ff]->rlmid.rmid_start : reads[ff]->b_adapter_pos ) );

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
                reads[ff]->rclip = min( reads[ff]->rlmid.rmid_start == 0 ? (int)reads[ff]->read.length() : reads[ff]->rlmid.rmid_start, min(reads[ff]->b_adapter_pos == 0 ? (int)reads[ff]->read.length()  : reads[ff]->b_adapter_pos, reads[ff]->lucy_rclip) );                
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
               
                reads[ff]->rclip = min( (int)reads[ff]->read.length(), min( reads[ff]->rlmid.rmid_start == 0 ? (int)reads[ff]->read.length() : reads[ff]->rlmid.rmid_start ,reads[ff]->b_adapter_pos == 0 ? (int)reads[ff]->read.length() : reads[ff]->b_adapter_pos ) );
                //cout << "B adapter:\n" << reads[ff]->b_adapter_pos << endl;
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
                
                if(reads[ff]->rclip <= 1) 
                {
                    reads[ff]->discarded = 1;
                }
                if(reads[ff]->rclip == (int)reads[ff]->read.length()) reads[ff]->right_trimmed_by_adapter = 0;
                
                //cout << reads[ff]->lclip << endl << reads[ff]->rclip << endl;

           } 
           else
           {
               reads[ff]->lclip = reads[ff]->rclip = 1;
               reads[ff]->discarded = 1;
           }
           
       }
       
       if(polyat_flag) {
         if(reads[ff]->discarded == 0) {
            if( (reads[ff]->rclip > reads[ff]->poly_A_clip) && (reads[ff]->poly_A_clip > 0)) {
              if(reads[ff]->rclip == reads[ff]->b_adapter_pos) {
                  reads[ff] ->right_trimmed_by_adapter = 0;
              } else if(reads[ff]->rclip == reads[ff]->rlmid.rmid_start) {
                  reads[ff] ->right_trimmed_by_adapter = 0;
              } else if(reads[ff]->rclip == reads[ff]->lucy_rclip) {
                  reads[ff]->right_trimmed_by_quality = 0;
              } else if(reads[ff]->rclip == reads[ff]->v_start) {
                  reads[ff]->right_trimmed_by_vector = 0;
              }
                    
              reads[ff]->rclip = reads[ff]->poly_A_clip;
              reads[ff]->right_trimmed_by_polyat = 1;
            }
         if( (reads[ff]->lclip < reads[ff]->poly_T_clip) && (reads[ff]->poly_T_clip > 0)) {
             if(reads[ff]->lclip == reads[ff]->rlmid.lmid_end) {
                reads[ff] ->left_trimmed_by_adapter = 0;
             } else if(reads[ff]->lclip == reads[ff]->lucy_lclip) {
                reads[ff]->left_trimmed_by_quality = 0;
             } else if(reads[ff]->lclip == reads[ff]->v_end) {
                reads[ff]->left_trimmed_by_vector = 0;
              } 
             
             reads[ff]->lclip = reads[ff]->poly_T_clip;
             reads[ff]->left_trimmed_by_polyat = 0;
           }
         }
       }
 }
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
    
    unsigned long right_trimmed_by_polyat, left_trimmed_by_polyat; right_trimmed_by_polyat = left_trimmed_by_polyat = 0;
    
    long left_mid_tag, right_mid_tag;
    left_mid_tag = right_mid_tag = 0;
    unsigned long long avg_bases, avg_right_clip_1, avg_left_clip_1;
    avg_bases = avg_right_clip_1 = avg_left_clip_1 = 0;
    
    double avg_trim_len, avg_read_len, cnt_avg, cnt_avg_len, avg_left_trim_len, avg_right_trim_len, cnt_avg_left_trim_len, cnt_avg_right_trim_len ;
    avg_trim_len = avg_read_len = cnt_avg = cnt_avg_len = avg_left_trim_len = avg_right_trim_len = cnt_avg_left_trim_len = cnt_avg_right_trim_len = 0;
    
    for(int i=0; i<(int)reads.size(); i++)
    {
        bases_anal+=reads[i]->initial_length;
        
       if(reads[i]->rlmid.lmid_start != 0 ) left_mid_tag+=1;
       if( (reads[i]->rlmid.rmid_start != 0) && (reads[i]->rlmid.rmid_start < (int)reads[i]->read.length()) ) right_mid_tag+=1;
       
       if(reads[i]->vector_found == 1) num_vectors+=1;
       if(reads[i]->contaminants == 1) num_contaminants+=1;
        
       if(reads[i]->discarded == 1) {
           discarded+=1;
           cout << reads[i]->rclip << " " << reads[i]->lclip << endl;
       }
       if(reads[i]->discarded == 0) accepted+=1;
       
       if(reads[i]->discarded_by_contaminant == 1) discarded_by_contaminant+=1;
       if(reads[i]->discarded_by_read_length == 1) discarded_by_read_length+=1;
                
       if(reads[i]->left_trimmed_by_quality == 1) left_trimmed_by_quality+=1;
       if(reads[i]->left_trimmed_by_adapter == 1) left_trimmed_by_adapter+=1;
       if(reads[i]->left_trimmed_by_vector == 1) left_trimmed_by_vector+=1;
                
       if(reads[i]->right_trimmed_by_quality == 1) right_trimmed_by_quality+=1;
       if(reads[i]->right_trimmed_by_adapter == 1) right_trimmed_by_adapter+=1;
       if(reads[i]->right_trimmed_by_vector == 1) right_trimmed_by_vector+=1;
        
        if(reads[i]->right_trimmed_by_polyat == 1) right_trimmed_by_polyat+=1;
        if(reads[i]->left_trimmed_by_polyat == 1) left_trimmed_by_polyat += 1;
        
       
       if( reads[i]->discarded == 0 )
       {
           cnt_avg_len+=1;
           avg_bases += reads[i]->rclip - reads[i]->lclip;
           avg_trim_len = avg_bases/cnt_avg_len;
           
           //Average right and left trimmed lengths
           if( reads[i]->initial_length > reads[i]->rclip )
           {
               //cnt_avg+=1;
               avg_right_clip_1 += reads[i]->initial_length - reads[i]->rclip;
               cnt_avg_right_trim_len+=1;
               avg_right_trim_len = avg_right_clip_1/cnt_avg_right_trim_len;
           }
           if( reads[i]->lclip > 0 )
           {
               cnt_avg_left_trim_len+=1;
               avg_left_clip_1 += reads[i]->lclip;
               avg_left_trim_len = avg_left_clip_1/cnt_avg_left_trim_len;
           }
       }
       
       
    }
    
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
    
    if(polyat_flag) {
        cout << "By poly A/T: " <<  left_trimmed_by_polyat << "\n";
        sum_stat << "By poly A/T: " <<  left_trimmed_by_polyat << "\n";
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
    
    if(polyat_flag) {
        cout << "By poly A/T: " <<  right_trimmed_by_polyat << "\n";
        sum_stat << "By poly A/T: " <<  right_trimmed_by_polyat << "\n";
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
    
    sum_stat_tsv << PrintRocheStatisticsTSV(reads.size(),
                                    bases_anal, 
                                    left_mid_tag,
                                    right_mid_tag,
                                    num_vectors, 
                                    num_contaminants, 
                                    left_trimmed_by_adapter,
                                    left_trimmed_by_quality,
                                    left_trimmed_by_vector,
                                    avg_left_trim_len, 
                                    right_trimmed_by_adapter, 
                                    right_trimmed_by_quality,
                                    right_trimmed_by_vector,
                                    avg_right_trim_len,
                                    discarded, 
                                    discarded_by_contaminant, 
                                    discarded_by_read_length,
                                    accepted, 
                                    avg_trim_len,
                                    left_trimmed_by_polyat,
                                    right_trimmed_by_polyat);
}



string PrintRocheStatisticsTSV(unsigned long cnt,
                                    unsigned long long bases_anal, 
                                    unsigned long left_mid_tag,
                                    unsigned long right_mid_tag,
                                    unsigned long num_vectors, 
                                    unsigned long num_contaminants, 
                                    unsigned long left_trimmed_by_adapter,
                                    unsigned long left_trimmed_by_quality,
                                    unsigned long left_trimmed_by_vector,
                                    double avg_left_trim_len, 
                                    unsigned long right_trimmed_by_adapter, 
                                    unsigned long right_trimmed_by_quality,
                                    unsigned long right_trimmed_by_vector,
                                    double avg_right_trim_len,
                                    unsigned long discarded, 
                                    unsigned long discarded_by_contaminant, 
                                    unsigned long discarded_by_read_length,
                                    unsigned long accepted, 
                                    double avg_trim_len,
                                    unsigned long left_trimmed_by_polyat,
                                    unsigned long right_trimmed_by_polyat
                               )
{
        string filename_str;
    
        for(int i=0; i<(int)roche_names.size(); ++i)
        {
            filename_str += string(roche_names[i]);
        }
        string stat_str_tsv =   version + "\t" + 
                                filename_str + "\t" +
                                int2str(NUM_THREADS) + "\t" +
                                (trim_adapters_flag  ? "YES" : "NO") + "\t" +
                                ( vector_flag ? "YES" : "NO") + "\t" +
                                ( vector_flag ? int2str(KMER_SIZE) : "NA") + "\t" +
                                ( vector_flag ? int2str(DISTANCE) : "NA") +"\t" +
                                ( contaminants_flag ? "YES" : "NO") +"\t" +
                                ( contaminants_flag ? int2str(KMER_SIZE_CONT) : "NA") +"\t" +
                                (qual_trim_flag ? "YES" : "NO") +"\t" +
                                (qual_trim_flag ? double2str(-10*log10(max_a_error)) : "NA")+ "\t" +
                                (qual_trim_flag ? double2str(-10*log10(max_e_at_ends)) : "NA")+ "\t" +
                                output_prefix +"\t" +
                                roche_rep_file_name + "\t"+
                                roche_output_file_name + (output_fastqfile_flag ? ", " + output_prefix + ".fastq" : "" ) + "\t" +
                                int2str(max_al_mism) +"\t" +
                                int2str(minimum_read_length)+ "\t"; 
                   
    
    
        stat_str_tsv += int2str(cnt)   + "\t" + //reads analyzed
                        int2str(bases_anal) + "\t"  + //bases
                        int2str(left_mid_tag) + "\t" 
                        + double2str( (double)left_mid_tag/(double)cnt*100.0) + "\t" + //perc 
                        int2str(right_mid_tag) + "\t" +
                        double2str( (double)right_mid_tag/(double)cnt*100.0) + "\t" + //perc 
                       ( vector_flag ? i2str(num_vectors,new char[15],10) + "\t" + double2str( (double)num_vectors/(double)cnt*100.0) + "\t" : "NA\tNA\t" ) + //perc vectors
                       ( contaminants_flag ? i2str(num_contaminants,new char[15],10) + "\t" //cont
                        + double2str( (double)num_contaminants/(double)cnt*100.0) + "\t" : "NA\tNA\t" ) + //perc cont
                        int2str(left_trimmed_by_adapter) + "\t" +
                       ( qual_trim_flag ? i2str(left_trimmed_by_quality,new char[15],10) + "\t" : "NA\t" ) +  //left trimmed qual
                       ( vector_flag ? i2str(left_trimmed_by_vector,new char[15],10) + "\t" : "NA\t" ) + //left trimmed vect
                       double2str(avg_left_trim_len) + "\t" + //avg left trim len
                       i2str(right_trimmed_by_adapter,new char[15],10) + "\t" + 
                       ( qual_trim_flag ? i2str(right_trimmed_by_quality,new char[15],10) + "\t" : "NA\t") +
                       ( vector_flag ? i2str(right_trimmed_by_vector,new char[15],10) + "\t" : "NA\t" ) +
                       double2str(avg_right_trim_len) + "\t" +
                       i2str(discarded,new char[15],10) + "\t" + //discard
                       ( contaminants_flag ? i2str(discarded_by_contaminant,new char[15],10) + "\t" : "NA\t" ) +
                       i2str(discarded_by_read_length,new char[15],10) + "\t" +
                       i2str(accepted,new char[15],10) + "\t" + //se reads kept
                        double2str( (double)accepted/(double)cnt*100.0) + "\t" +//perc kept
                       double2str(avg_trim_len) + 
                        (polyat_flag ? "\tYes\t" + int2str(cdna) + "\t" + int2str(c_err) + "\t" + int2str(crng) + "\t" + int2str(left_trimmed_by_polyat) + "\t" + int2str(right_trimmed_by_polyat) : "");
            
    
     return stat_str_tsv;
    
    
}
