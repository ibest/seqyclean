#include "MainPipeLine.h"

void MainPipeLine() {
    
    /*Here we are making threads*/
    /*========================================================================*/
    pthread_t threads[NUM_THREADS];
    
    /*Remove Keys/Adaptors/Primers/Barcodes*/
    int i=0;
    for( i=0; i<(int)reads.size()-NUM_THREADS; i+=NUM_THREADS ) {
        
        for(int t=0; t<NUM_THREADS; t++) {
           pthread_create( &threads[t], NULL, &t_FindRLClip, (void *)(t+i) );
        }
    
        for(int t=0; t<NUM_THREADS; t++) {
           pthread_join(threads[t], NULL);
        }
        
    }
    
    /*Processing the rest of the set of reads: */
    for(int t=0; t< (int)reads.size() - i; t++) {
       pthread_create( &threads[t], NULL, &t_FindRLClip, (void *)(t+i) );
       
    }
   
    for(int t=0; t<(int)reads.size() - i; t++) {
       pthread_join(threads[t], NULL);
    } 
    
    /*If vector file is set up, perform vector tails cleaning*/
    if(vector_flag == true) {
        cout << "Trimming small pieces of vector..." << endl;
        cout << "Right ends..." << endl;
        TrimRightEnds();
        cout << "Left ends..." << endl;
        TrimLeftEnds();
    }
    
    /*If PCR/Amplicon flags are set -> perform trimming the primers*/
    if(pcr_flag == true) {
        cout << "Trimming primers..." << endl;
      
        //Building a new dictionary for RL MIDS : 
        Build_RLMIDS_Dictionary(pcr_file_name);
        i=0;
        for( i=0; i<(int)reads.size()-NUM_THREADS; i+=NUM_THREADS ) {
        
                for(int t=0; t<NUM_THREADS; t++) {
                        pthread_create( &threads[t], NULL, &t_FindClipAmplicon, (void *)(t+i) );
                }
    
                for(int t=0; t<NUM_THREADS; t++) {
                        pthread_join(threads[t], NULL);
                }
        }
    
        //Processing the rest of the set of reads:
        for(int t=0; t< (int)reads.size() - i; t++) {
                pthread_create( &threads[t], NULL, &t_FindClipAmplicon, (void *)(t+i) );
       
        }
   
        for(int t=0; t<(int)reads.size() - i; t++) {
                pthread_join(threads[t], NULL);
        }
    }
   
}

bool GetRochBAdapter(Read* read) {
    /*
        >Common_Adapter_A
        CCATCTCATCCCTGCGTGTCTCCGACGACT
        >Common_Adapter_B
        CCTATCCCCTGTGTGGACC
    */
    
    /*First, search for the whole b-adaptor*/
    iz_SSAHA *izssaha = new iz_SSAHA();
    string b_adapter = "CCTATCCCCTGTGTGGACC";
    int read_len = read->read.length();
    int from; 
    if(read_len <= 100) {
      from = 0;
    } else if( (read_len <=300) && (read_len > 100) ) {
      from = read_len - read_len/2; //divide by 2
    } else if( (read_len <= 600) && (read_len > 300) ) {
      from = read_len - read_len/4; //divide by 4
    } else  {
      from = read_len - read_len/8; //divide by 8
    } 
  
    string ref_str = read->read.substr( from, read_len - from );
  
    /*Normal*/
    AlignResult al_res = izssaha->Find( ref_str , b_adapter );
    AlignScores scores;
    if(al_res.found_flag == true) {
        scores = CalcScores(al_res.seq_1_al, al_res.seq_2_al, al_res.seq_1_al.length(), 0);
        
        if(scores.mismatches <= max_al_mism  ) {
           read->b_adapter_pos = al_res.pos + from;
           read->b_adapter = b_adapter;
           read->b_adapter_err = scores.mismatches;
           
           delete izssaha;
           al_res.seq_1_al.clear();
           al_res.seq_2_al.clear();
           b_adapter.clear();
           
           return true;
        } 
        
    }
    
    /*Reverse complement*/
    reverse( b_adapter.begin(), b_adapter.end() );
    b_adapter = MakeSeqComplement(b_adapter);
    al_res = izssaha->Find( ref_str , b_adapter );
    if(al_res.found_flag == true) {
        scores = CalcScores(al_res.seq_1_al, al_res.seq_2_al, al_res.seq_1_al.length(), 0);
        
  
        if(scores.mismatches <= max_al_mism  ) {
           read->b_adapter_pos = al_res.pos + from;
           read->b_adapter = b_adapter;
           read->b_adapter_err = scores.mismatches;
           
           delete izssaha;
           al_res.seq_1_al.clear();
           al_res.seq_2_al.clear();
           //b_adapter.clear();
           
           return true;
        } 
        
    }
    
    
    /*Then for a half of it*/
    /*Normal*/
    b_adapter = "CCTATCCCC";
    al_res = izssaha->Find( ref_str , b_adapter );
    if(al_res.found_flag == true) {
        scores = CalcScores(al_res.seq_1_al, al_res.seq_2_al, al_res.seq_1_al.length(), 0);
        
        if(scores.mismatches <= max_al_mism  ) {
           read->b_adapter_pos = al_res.pos + from;
           read->b_adapter = b_adapter;
           read->b_adapter_err = scores.mismatches;
           
           delete izssaha;
           al_res.seq_1_al.clear();
           al_res.seq_2_al.clear();
           //b_adapter.clear();
           
           return true;
        } 
        
    }
    
    /*Reverse complement*/
    reverse( b_adapter.begin(), b_adapter.end() );
    b_adapter = MakeSeqComplement(b_adapter);
    al_res = izssaha->Find( ref_str , b_adapter );
    if(al_res.found_flag == true) {
        scores = CalcScores(al_res.seq_1_al, al_res.seq_2_al, al_res.seq_1_al.length(), 0);
        
        if(scores.mismatches <= max_al_mism  ) {
           read->b_adapter_pos = al_res.pos + from;
           read->b_adapter = b_adapter;
           read->b_adapter_err = scores.mismatches;
           
           delete izssaha;
           al_res.seq_1_al.clear();
           al_res.seq_2_al.clear();
           b_adapter.clear();
           
           return true;
        } 
        
    }
    
    
    delete izssaha;
    al_res.seq_1_al.clear();
    al_res.seq_2_al.clear();
    b_adapter.clear();
    
    return false;
}

void TrimRightEnds() {
    pthread_t threads[NUM_THREADS];
    
    /*First stage: preprocess the read by using SSAHA: */
    long lim = reads.size() - NUM_THREADS;
    for(int i=0; i<(int)reads.size(); i+=NUM_THREADS) {
        
        if(i > lim ) break;
        
        for(int t=0; t<NUM_THREADS; t++) {
           pthread_create( &threads[t], NULL, &t_TrimRightEnds, (void *)(t+i) );//(void *)&t_args);
        }
    
        for(int t=0; t<NUM_THREADS; t++) {
           pthread_join(threads[t], NULL);
        }
    }
    
}

/*Looks for the left position of the vector: */
void TrimLeftEnds() {
    pthread_t threads[NUM_THREADS];
    
    long lim = reads.size() - NUM_THREADS;
    for(int i=0; i<(int)reads.size(); i+=NUM_THREADS) {
        
        if(i > lim ) break;
        
       for(int t=0; t<NUM_THREADS; t++) {
           pthread_create( &threads[t], NULL, &t_TrimLeftEnds, (void *)(t+i) );//(void *)&t_args);
       }
    
        for(int t=0; t<NUM_THREADS; t++) {
           pthread_join(threads[t], NULL);
        
        }
    }
  
}

/*Threaded method for Contaminants Dictionary*/
static void *t_FindRLClip(void *targs) {
  /*Getting thread ID tid:*/
  long i = 0;
  i = (long)targs;
  
  /*if( reads[i].discarded == 1 )
  {
      pthread_exit(NULL);
      return 0;
  }**/
  
  /* Remove not-needed Ns: */
  for(int s=0; s<10; s++) {
        TrimNs(reads[i]->read);
        reads[i]->quality = reads[i]->quality.substr(0, reads[i]->read.length());
        reads[i]->clear_length = reads[i]->read.length();
  }
  
  
  /*Obtaing a left clip point*/
  GetLClip2(reads[i], false);
  
  reads[i]->clip_found = false;
  //Right clip point:
  for(int kk=0; kk<(int)rlmids.size(); kk++ ) {
      size_t found;
      
      found = reads[i]->read.rfind(rlmids[kk].rmid_value);
      if( found != string::npos )
      {
          reads[i]->rlmid.rmid_name = rlmids[kk].rmid_name;
          reads[i]->rlmid.rmid_id = rlmids[kk].rmid_id;
          reads[i]->rlmid.rmid_value = rlmids[kk].rmid_value;
                
          reads[i]->clip_found = true;
               
          reads[i]->rlmid.rmid_start = found;
          reads[i]->rlmid.rmid_end = found + rlmids[kk].rmid_value.length();
          reads[i]->rlmid.rmid_err = 0;
      }
      else
      {
          //try reverse complement:
          string mid = rlmids[kk].rmid_value;
          reverse( mid.begin(), mid.end() );
          found = reads[i]->read.rfind(mid);
          if( found != string::npos )
          {
                reads[i]->rlmid.rmid_name = rlmids[kk].rmid_name;
                reads[i]->rlmid.rmid_id = rlmids[kk].rmid_id;
                reads[i]->rlmid.rmid_value = rlmids[kk].rmid_value;
                
                reads[i]->clip_found = true;
               
                reads[i]->rlmid.rmid_start = found;
                reads[i]->rlmid.rmid_end = found + rlmids[kk].rmid_value.length();
                reads[i]->rlmid.rmid_err = 0;
          }
      }
  }
  
  if( reads[i]->clip_found == false )
  {
  
        if( GetRochBAdapter( reads[i] ) == true ) {
                reads[i]->rlmid.rmid_name = rlmids[reads[i]->rlmid.lmid_id].rmid_name;
                reads[i]->rlmid.rmid_id = rlmids[reads[i]->rlmid.lmid_id].rmid_id;
                reads[i]->rlmid.rmid_value = rlmids[reads[i]->rlmid.lmid_id].rmid_value;
                
                reads[i]->clip_found = true;
                reads[i]->R_rclip = reads[i]->rclip + rlmids[0].rmid_value.length();
                
                reads[i]->rlmid.rmid_start = reads[i]->b_adapter_pos - rlmids[0].rmid_value.length();
                reads[i]->rlmid.rmid_end = reads[i]->b_adapter_pos;
                reads[i]->rlmid.rmid_err = reads[i]->b_adapter_err;
        } 
        else 
        {
  
                pthread_t rlthreads[rlmids.size()];
                /*Look for RCLIP*/
                //Serial realization  
                iz_SSAHA *izssaha = new iz_SSAHA();
                string rlstr = rlmids[reads[i]->rlmid.lmid_id].rmid_value;
                //string BADAPTER = "GGTCGGCGTCTCTCAAGGCACACAGGGGATAGG";
                int read_len = reads[i]->read.length();
                int from; 
                if(read_len <= 100) {
                        from = 0;
                } else if( (read_len <=300) && (read_len > 100) ) {
                        from = read_len - read_len/2; //divide by 2
                } else if( (read_len <= 600) && (read_len > 300) ) {
                        from = read_len - read_len/4; //divide by 4
                } else  {
                        from = read_len - read_len/8; //divide by 8
                } 
  
                string ref_str = reads[i]->read.substr( from, read_len - from );
  
                /*Normal*/
                AlignResult al_res = izssaha->Find( ref_str , rlstr );
                AlignScores scores;
                if(al_res.found_flag == true) 
                {
                        scores = CalcScores(al_res.seq_1_al, al_res.seq_2_al, al_res.seq_1_al.length(), 0);
                        delete izssaha;
  
                        if(scores.mismatches <= max_al_mism  ) 
                        {
                                reads[i]->rclip = al_res.pos + from;
                                if(reads[i]->rclip <= 16) reads[i]->rclip = reads[i]->read.length();
                                
                                /*Right mid*/
                                reads[i]->rlmid.rmid_name = rlmids[reads[i]->rlmid.lmid_id].rmid_name;
                                reads[i]->rlmid.rmid_id = rlmids[reads[i]->rlmid.lmid_id].rmid_id;
                                reads[i]->rlmid.rmid_value = rlmids[reads[i]->rlmid.lmid_id].rmid_value;
                
                                reads[i]->clip_found = true;
                                reads[i]->R_rclip = reads[i]->rclip + rlmids[0].rmid_value.length();
                
                                reads[i]->rlmid.rmid_start = al_res.pos + from;
                                reads[i]->rlmid.rmid_end = al_res.pos + from + rlmids[0].rmid_value.length();
                                reads[i]->rlmid.rmid_err = scores.mismatches;
                        } else 
                        {
                                //Preprocess the read by using SSAHA:
                                int jj =0;
                                for(jj=0; jj<(int)rlmids.size() - 4; jj+=4 ) 
                                {
                                        if(reads[i]->clip_found == true) break;
      
                                        t_args rlargs0;
                                        rlargs0.rec_id = i;
                                        rlargs0.tid = jj;
                                        rlargs0.reverse = false;
                                        rlargs0.complement = false;
                                        rlargs0.reverse_complement = false;
                                        pthread_create( &rlthreads[jj], NULL, &tt_SSAHA, (void*)&rlargs0 );
      
                                        t_args rlargs1;
                                        rlargs1.rec_id = i;
                                        rlargs1.tid = jj+1;
                                        rlargs1.reverse = false;
                                        rlargs1.complement = false;
                                        rlargs1.reverse_complement = false;
                                        pthread_create( &rlthreads[jj+1], NULL, &tt_SSAHA, (void*)&rlargs1 );
      
                                        t_args rlargs2;
                                        rlargs2.rec_id = i;
                                        rlargs2.tid = jj+2;
                                        rlargs2.reverse = false;
                                        rlargs2.complement = false;
                                        rlargs2.reverse_complement = false;
                                        pthread_create( &rlthreads[jj+2], NULL, &tt_SSAHA, (void*)&rlargs2 );
      
                                        t_args rlargs3;
                                        rlargs3.rec_id = i;
                                        rlargs3.tid = jj+3;
                                        rlargs3.reverse = false;
                                        rlargs3.complement = false;
                                        rlargs3.reverse_complement = false;
                                        pthread_create( &rlthreads[jj+3], NULL, &tt_SSAHA, (void*)&rlargs3 );
      
                                        pthread_join(rlthreads[jj], NULL);
                                        pthread_join(rlthreads[jj+1], NULL);
                                        pthread_join(rlthreads[jj+2], NULL);
                                        pthread_join(rlthreads[jj+3], NULL);
                                }
                                if(reads[i]->clip_found == false) 
                                { 
                                        for(int t=0; t< (int)rlmids.size() - jj; t++) 
                                        {
                                                t_args rlargs0;
                                                rlargs0.rec_id = i;
                                                rlargs0.tid = jj+1;
                                                rlargs0.reverse = false;
                                                rlargs0.complement = false;
                                                rlargs0.reverse_complement = false;
                                                pthread_create( &rlthreads[t], NULL, &tt_SSAHA, (void*)&rlargs0 );
               
                                                pthread_join(rlthreads[t], NULL);
                                        }
           
                                }
                        }
                }
  
  
                /*Try Reverse Complement*/
                if( reads[i]->clip_found == false ) 
                {
      
                        reverse( rlstr.begin(), rlstr.end() );
                        rlstr = MakeSeqComplement(rlstr);
                        iz_SSAHA *izssaha = new iz_SSAHA();
                        al_res = izssaha->Find( ref_str , rlstr );
                        if(al_res.found_flag == true) 
                        {
                                scores = CalcScores(al_res.seq_1_al, al_res.seq_2_al, al_res.seq_1_al.length(), 0);
        
                                if(scores.mismatches <= max_al_mism  ) 
                                {
                                        reads[i]->rclip = al_res.pos + from;
                                        if(reads[i]->rclip <= 16) reads[i]->rclip = reads[i]->read.length();
                                }
    
                                //reads[i]->rev_bc = rlmids[reads[i]->rlmid.lmid_id]->rname;
                                reads[i]->rlmid.rmid_name = rlmids[reads[i]->rlmid.lmid_id].rmid_name;
                                reads[i]->rlmid.rmid_id = rlmids[reads[i]->rlmid.lmid_id].rmid_id;
                                reads[i]->rlmid.rmid_value = rlmids[reads[i]->rlmid.lmid_id].rmid_value;
                        
                                reads[i]->clip_found = true;
                                //reads[i]->rclip_errors = scores.mismatches;
                
                                reads[i]->R_rclip = reads[i]->rclip + rlmids[0].rmid_value.length();
                
                                reads[i]->rlmid.rmid_start = al_res.pos + from;
                                reads[i]->rlmid.rmid_end = al_res.pos + from + rlmids[0].rmid_value.length();
                                reads[i]->rlmid.rmid_err = scores.mismatches;
                        } else 
                        {
                                //Preprocess the read by using SSAHA:
        
                                int jj =0;
                                for(jj=0; jj<(int)rlmids.size() - 4; jj+=4 ) 
                                {
                                        //cout << jj << endl;
                                        if(reads[i]->clip_found == true) break;
      
                                        t_args rlargs0;
                                        rlargs0.rec_id = i;
                                        rlargs0.tid = jj;
                                        rlargs0.reverse = false;
                                        rlargs0.complement = false;
                                        rlargs0.reverse_complement = true;
                                        pthread_create( &rlthreads[jj], NULL, &tt_SSAHA, (void*)&rlargs0 );
      
                                        t_args rlargs1;
                                        rlargs1.rec_id = i;
                                        rlargs1.tid = jj+1;
                                        rlargs1.reverse = false;
                                        rlargs1.complement = false;
                                        rlargs1.reverse_complement = true;
                                        pthread_create( &rlthreads[jj+1], NULL, &tt_SSAHA, (void*)&rlargs1 );
      
                                        t_args rlargs2;
                                        rlargs2.rec_id = i;
                                        rlargs2.tid = jj+2;
                                        rlargs2.reverse = false;
                                        rlargs2.complement = false;
                                        rlargs2.reverse_complement = true;
                                        pthread_create( &rlthreads[jj+2], NULL, &tt_SSAHA, (void*)&rlargs2 );
      
                                        t_args rlargs3;
                                        rlargs3.rec_id = i;
                                        rlargs3.tid = jj+3;
                                        rlargs3.reverse = false;
                                        rlargs3.complement = false;
                                        rlargs3.reverse_complement = true;
                                        pthread_create( &rlthreads[jj+3], NULL, &tt_SSAHA, (void*)&rlargs3 );
      
                                        pthread_join(rlthreads[jj], NULL);
                                        pthread_join(rlthreads[jj+1], NULL);
                                        pthread_join(rlthreads[jj+2], NULL);
                                        pthread_join(rlthreads[jj+3], NULL);
                
                                }
                                if(reads[i]->clip_found == false) 
                                { 
                                        for(int t=0; t< (int)rlmids.size() - jj; t++) 
                                        {
                                                t_args rlargs0;
                                                rlargs0.rec_id = i;
                                                rlargs0.tid = jj+1;
                                                rlargs0.reverse = false;
                                                rlargs0.complement = false;
                                                rlargs0.reverse_complement = true;
                                                pthread_create( &rlthreads[t], NULL, &tt_SSAHA, (void*)&rlargs0 );
               
                                                pthread_join(rlthreads[t], NULL);
                                        }
                                }
                        }
        
                        delete izssaha;
        
                }
                
                al_res.seq_1_al.clear();
                al_res.seq_2_al.clear();
                rlstr.clear();
        }
  }
  
  if(vector_flag == true) 
  {
     CheckVector(reads[i]); 
       
     if( ( reads[i]->read.length() - (reads[i]->v_end - reads[i]->v_start) ) < minimum_read_length ) 
     {
        reads[i]->discarded = 1;
        reads[i]->discarded_by_read_length = 1;
        //pthread_exit(NULL); return 0;
     }
      
  } 
  
  
  line_counter++;
  if(line_counter%50000 == 0)
     cout << "Line No: " << line_counter << endl;
         
  pthread_exit(NULL);
    
}

static void *tt_SSAHA(void *targs) {
    t_args foo = *((t_args*)targs);
  
    iz_SSAHA *izssaha = new iz_SSAHA();
    string query_str = rlmids[foo.tid].rmid_value;// + "GGTCGGCGTCT";
    int read_len = reads[foo.rec_id]->read.length();
    int from; 
    if(read_len <= 100) {
        from = 0;
    } else if( (read_len <=300) && (read_len > 100) ) {
        from = read_len - read_len/2; //divide by 2
    } else if( (read_len <= 600) && (read_len > 300) ) {
        from = read_len - read_len/4; //divide by 4
    } else  {
        from = read_len - read_len/8; //divide by 8
    } 
    
    string ref_str = reads[foo.rec_id]->read.substr( from, read_len - from );
    
    if( foo.reverse == true ) {
        reverse( query_str.begin(), query_str.end() );
    }
    if( foo.complement == true ) {
        query_str =  MakeSeqComplement(query_str);
    }
    if( foo.reverse_complement == true ) {
        query_str =  MakeSeqComplement(query_str);
        reverse( query_str.begin(), query_str.end() );
    }
    
    //cout << from << " " << read_len << endl << reads[foo.rec_id]->readID << endl << ref_str << endl;
    if (reads[foo.rec_id]->clip_found == false) {
        AlignResult al_res = izssaha->Find( ref_str , query_str );
        if (al_res.found_flag == true) {
                AlignScores scores;
                //cout << al_res.seq_1_al  << endl << al_res.seq_2_al << endl;
                scores = CalcScores(al_res.seq_1_al, al_res.seq_2_al, al_res.seq_1_al.length(), 0);
                
                if( scores.mismatches <= max_al_mism ) {
                        //if(reads[foo.rec_id]->readID=="@HLIQV0O02IFCPA") 
                        //        cout << reads[foo.rec_id]->rclip << endl << al_res.seq_1_al << endl << al_res.seq_2_al << endl;
          
                        if (reads[foo.rec_id]->clip_found == false) 
                        {
                                reads[foo.rec_id]->rclip = al_res.pos + from;
                                if(reads[foo.rec_id]->rclip <= 16) {
                                        reads[foo.rec_id]->rclip = reads[foo.rec_id]->read.length();;
                                }
                                
                                /*Right MID*/
                                reads[foo.rec_id]->rlmid.rmid_name = rlmids[foo.tid].rmid_name;
                                reads[foo.rec_id]->rlmid.rmid_id = rlmids[foo.tid].rmid_id;
                                reads[foo.rec_id]->rlmid.rmid_value = rlmids[foo.tid].rmid_value;
                                
                                reads[foo.rec_id]->clip_found = true;
                                reads[foo.rec_id]->rclip_errors = scores.mismatches;
                                reads[foo.rec_id]->R_rclip = reads[foo.rec_id]->rclip + rlmids[0].rmid_value.length();
                            
                                reads[foo.rec_id]->rlmid.rmid_start = al_res.pos + from;
                                reads[foo.rec_id]->rlmid.rmid_end = al_res.pos + from + rlmids[0].rmid_value.length();
                                reads[foo.rec_id]->rlmid.rmid_err = scores.mismatches;
                                
                                
                        }
                } 
                if(reads[foo.rec_id]->clip_found == false) reads[foo.rec_id]->rclip_errors = scores.mismatches;
        
        }
        
        al_res.seq_1_al.clear();
        al_res.seq_2_al.clear();
        query_str.clear();
        ref_str.clear();
        delete izssaha;
    }
      
    pthread_exit(NULL);
}

/*Threaded method for Contaminants Dictionary*/
static void *t_TrimRightEnds(void *targs) { 
    bool rclip_found = false;
    long i = 0;
    i = (long)targs;
    
     if( (reads[i]->rlmid.rmid_end <= 16) && (reads[i]->rlmid.rmid_start == reads[i]->clear_length) && (reads[i]->discarded == 0) ) {
            rclip_found = false;
            //cout << reads[i]->readID << endl;
            //Good candidate to trim ends
            //Trimming the right end:
            //Doing pairwise alignment, starting from the (KMER_SIZE-1) and the right end of the read:
            for(int j=1; j<KMER_SIZE-2; j++) {
                
                string right_end = reads[i]->read.substr(reads[i]->read.length() - (KMER_SIZE-j), KMER_SIZE-j );
                
                
                //Looking at heads (beginnings) of vector genome string: 
                for(int k=0; k<(int)VectorSeqs.size(); k++) {
                    string vector_head = VectorSeqs[k].substr(0,KMER_SIZE-j);
                    //Banded pairwise alignment: 
                    //AlignResult al_res = banded(right_end, vector_head, 2, 1, 0, 0);
                    //Calculating number of mismatches : 
                    //AlignScores scores = CalcScores(al_res.seq_1_al,al_res.seq_2_al, KMER_SIZE-j, 0);
                    
                    //!!!WARNING!!! 
                     // Number of mismatches can be tunable!
                    
                    //if(scores.mismatches == 0) {
                    if(right_end == vector_head) {
                        //reads[i]->rclip -= KMER_SIZE-j;
                        //reads[i]->read = reads[i]->read.substr(0, reads[i]->rclip );
                        rclip_found = true;
                        //cout << reads[i]->rclip << " " << reads[i]->lclip << endl;
                        //cout << reads[i]->readID << endl;
                        reads[i]->r_vec_start -= KMER_SIZE-j;
                        break;
                    } 
                }
                
                if(rclip_found == true) 
                        break;
                
                //Looking at tail (left end) of the vector genome string: 
                for(int k=0; k<(int)VectorSeqs.size(); k++) {
                    string vector_tail = VectorSeqs[k].substr(VectorSeqs[k].length() - (KMER_SIZE-j) ,KMER_SIZE-j );
                    //Banded pairwise alignment: 
                    //AlignResult al_res = banded(right_end, vector_tail, 2, 1, 0, 0);
                    //Calculating number of mismatches : 
                    //AlignScores scores = CalcScores(al_res.seq_1_al,al_res.seq_2_al, KMER_SIZE-j, 0);
                    
                    //!!!WARNING!!! 
                     // Number of mismatches can be tunable!
                    
                    //if(scores.mismatches == 0 ) {
                    if(right_end == vector_tail) {
                        //reads[i].rclip -= KMER_SIZE-j;
                        //reads[i].read = reads[i].read.substr(0, reads[i].rclip );
                        rclip_found = true;
                        reads[i]->r_vec_start -= KMER_SIZE-j;
                        //cout << reads[i].rclip << " " << reads[i].lclip << endl;
                        //cout << reads[i].readID << endl;
                        break;
                    } 
                }
                
             if(rclip_found == true) {
                 break;
             } else {
                 rclip_found = false;
                 continue;
             }
                
            }
        }
    
    pthread_exit(NULL);
    
}


/*Threaded method for Contaminants Dictionary*/
static void *t_TrimLeftEnds(void *targs) { 
    bool lclip_found = false;
    
    long i = 0;
    i = (long)targs;
    
     if( reads[i]->discarded == 0 ) 
     {
        
        lclip_found = false;
            
         /*Good candidate to trim ends*/
         /*Trimming the right end:*/
         /*Doing pairwise alignment, starting from the (KMER_SIZE-1) and the right end of the read:*/
         for(int j=1; j<KMER_SIZE-2; j++) 
         {
            string left_end = reads[i]->read.substr(reads[i]->lclip, KMER_SIZE-j );
                
            /*Looking at heads (beginnings) of vector genome string: */
            for(int k=0; k<(int)VectorSeqs.size(); k++) 
            {
               string vector_head = VectorSeqs[k].substr(0,KMER_SIZE-j);
               /*Banded pairwise alignment: */
               /*Calculating number of mismatches : */
              /*!!!WARNING!!! 
               * Number of mismatches can be tunable!
               */
              if(left_end == vector_head) 
              {
                 //reads[i].lclip += KMER_SIZE-j;
                 lclip_found = true;
                 reads[i]->l_vec_end += KMER_SIZE-j;
                 reads[i]->l_vec_start = reads[i]->rlmid.lmid_end;
                 //if (i == 1) cout << "! " << reads[i].vec_end << endl;
                 break;
              } 
            }
                
            if(lclip_found == true) 
               break;
                
            /*Looking at tail (left end) of the vector genome string: */
            for(int k=0; k<(int)VectorSeqs.size(); k++) 
            {
                string vector_tail = VectorSeqs[k].substr(VectorSeqs[k].length() - (KMER_SIZE-j) ,KMER_SIZE-j );
                /*Banded pairwise alignment: */
                //AlignResult al_res = banded(left_end, vector_tail, 2, 1, 0, 0);
                /*Calculating number of mismatches : */
                /*!!!WARNING!!! 
                 * Number of mismatches can be tunable!
                 */
                if(left_end == vector_tail) 
                {
                   //reads[i].lclip += KMER_SIZE-j;
                   lclip_found = true;
                   reads[i]->l_vec_end += KMER_SIZE-j;
                   reads[i]->l_vec_start = reads[i]->rlmid.lmid_end;
                   //if (i == 1) cout << "!! " << reads[i].vec_end << endl;
                   break;
                } 
            }
                
            if(lclip_found == true) 
            {
              break;
            } else 
            {
              lclip_found = false;
              continue;
            }
                
        }
     }
     
    
    return 0;
}

void GetLClip2(Read* read,bool pflag) {
    
    string sss;
    bool lclip_found = false;
  
    if(lclip_found == false) {
        /*SSAHA job begins:*/
        int best_pos = 0;
        int best_k = 0;
        int best_mism = 20;
        
        for(int k=0; k<(int)rlmids.size(); k++) {
            sss = rlmids[k].lmid_value;
            iz_SSAHA *izssaha = new iz_SSAHA();
            
            string ref_str = read->read.substr( 0, 40 ) ;
            stoupper(ref_str);
            AlignResult al_res = izssaha->Find( ref_str , sss );
            if(al_res.found_flag == true) {
                AlignScores scores;
                scores = CalcScores2( al_res.seq_1_al, al_res.seq_1_al.length(), 0 );
            
                delete izssaha;
                //cout << scores.mismatches << endl;
                if(scores.mismatches <= 2  ) {
                        if( best_mism > scores.mismatches ) {
                                best_pos = al_res.pos;
                                best_k = k;
                                best_mism = scores.mismatches;
                                
                                if(pflag == false) {
                                        read->lclip = best_pos + rlmids[best_k].lmid_value.length() + 1;
                                        read->rlmid.lmid_name = rlmids[best_k].lmid_name;
                                        //read.L_lclip = best_pos  + 1 + pos.pos_left;
                                        read->rlmid.lmid_start = best_pos + 1;
                                        read->rlmid.lmid_end = best_pos + rlmids[best_k].lmid_value.length() + 1;
                                        read->rlmid.lmid_err = best_mism;
                                } else {
                                        read->fwd_bc = rlmids[best_k].lmid_name;
                                        read->fwdP_start = best_pos + 1;
                                        read->fwdP_end =  best_pos + rlmids[best_k].lmid_value.length() + 1;     
                                        read->fwdP_errors = best_mism;        
                                }
                                
                                lclip_found = true;
                                
                                AlignResult pos = CalcPos(al_res.seq_2_al);
                   
                                
                        }
                
                        al_res.seq_1_al.clear();
                        al_res.seq_2_al.clear();
                        sss.clear();
                
                
                } //else {
                read->lclip_errors = best_mism;//scores.mismatches;
                
            //}
            }
            al_res.seq_1_al.clear();
            al_res.seq_2_al.clear();
            sss.clear();
            
            //if(lclip_found == true) break;
        }
    }
    
    /*Reverse complement*/
    if(lclip_found == false) {
        /*SSAHA job begins:*/
        int best_pos = 0;
        int best_k = 0;
        int best_mism = 20;
        
        for(int k=0; k<(int)rlmids.size(); k++) {
            sss = rlmids[k].lmid_value;
            reverse( sss.begin(), sss.end() );
            sss = MakeSeqComplement(sss);
            
            iz_SSAHA *izssaha = new iz_SSAHA();
            
            string ref_str = read->read.substr( 0, 40 ) ;
            stoupper(ref_str);
            
            AlignResult al_res = izssaha->Find( ref_str , sss );
            if(al_res.found_flag == true) {
                AlignScores scores;
                scores = CalcScores2( al_res.seq_1_al, al_res.seq_1_al.length(), 0 );
            
                delete izssaha;
                //cout << scores.mismatches << endl;
                if(scores.mismatches <= 2  ) {
                        if( best_mism > scores.mismatches ) {
                                best_pos = al_res.pos;
                                best_k = k;
                                best_mism = scores.mismatches;
                   
                                if(pflag == false) {
                                        read->lclip = best_pos + rlmids[best_k].lmid_value.length() + 1;
                                        read->rlmid.lmid_name = rlmids[best_k].lmid_name;
                                        //read.L_lclip = best_pos  + 1 + pos.pos_left;
                                        read->rlmid.lmid_start = best_pos + 1;
                                        read->rlmid.lmid_end = best_pos + rlmids[best_k].lmid_value.length() + 1;
                                        read->rlmid.lmid_err = best_mism;
                                } else {
                                        read->fwd_bc = rlmids[best_k].lmid_name;
                                        read->fwdP_start = best_pos + 1;
                                        read->fwdP_end =  best_pos + rlmids[best_k].lmid_value.length() + 1;     
                                        read->fwdP_errors = best_mism;        
                                }
                                
                                lclip_found = true;
                                read->lclip_errors = best_mism;
                   
                                AlignResult pos = CalcPos(al_res.seq_2_al);
                   
                                
                        }
                
                        al_res.seq_1_al.clear();
                        al_res.seq_2_al.clear();
                        sss.clear();
                
                
                } //else {
                read->lclip_errors = best_mism;//scores.mismatches;
                
            //}
            }
            al_res.seq_1_al.clear();
            al_res.seq_2_al.clear();
            sss.clear();
            
            if(lclip_found == true) break;
        }
    }
}

/*Threaded method for Contaminants Dictionary*/
static void *t_FindClipAmplicon(void *targs) {
  /*Getting thread ID tid:*/
  long i = 0;
  i = (long)targs;
  
  
  reads[i]->clip_found = false;
  reads[i]->fwd_bc = "";
  reads[i]->rev_bc = "";
  if(reads[i]->read.length() > 50) {
        
        GetLClip2(reads[i],true);
  /*if(reads[i].readID=="@HRCRL2004IH0W8") 
           cout << reads[i].read.length() << endl << reads[i].quality.length() << endl;**/
  
        pthread_t rlthreads[rlmids.size()];
        //Look for RCLIP
        //Serial realization  
        iz_SSAHA *izssaha = new iz_SSAHA();
        string rlstr = rlmids[reads[i]->rlmid.lmid_id].rmid_value;
  
        int read_len = reads[i]->read.length();
        int from; 
        if(read_len <= 100) {
                from = 0;
        } else if( (read_len <=300) && (read_len > 100) ) {
                from = read_len - read_len/2; //divide by 2
        } else if( (read_len <= 600) && (read_len > 300) ) {
                from = read_len - read_len/4; //divide by 4
        } else  {
                from = read_len - read_len/8; //divide by 8
        } 
  
        string ref_str = reads[i]->read.substr( from, read_len - from );
  
        //Normal
        AlignResult al_res = izssaha->Find( ref_str , rlstr );
        AlignScores scores;
        if( al_res.found_flag == true ) {
                
  
                //scores = CalcScores(al_res.seq_1_al, al_res.seq_1_al.length(), 0);
                scores = CalcScores(al_res.seq_1_al, al_res.seq_2_al, al_res.seq_1_al.length(), 0);
  
                if(scores.mismatches <= max_al_mism  ) {
                        reads[i]->rclip = al_res.pos + from;
                        if(reads[i]->rclip <= 16) {
                                reads[i]->rclip = reads[i]->read.length();
                        }
                        reads[i]->rev_bc = rlmids[reads[i]->rlmid.lmid_id].rmid_name;//string(itoa(reads[i].rlmid.lmid_id+1, new char[5], 10));
                        reads[i]->clip_found = true;
                        reads[i]->revP_errors = scores.mismatches;
                        reads[i]->R_rclip = reads[i]->rclip + rlmids[0].rmid_value.length();
                        
                        reads[i]->revP_start = al_res.pos + from; //cout << al_res.pos << " " << from << endl;
                        reads[i]->revP_end = reads[i]->revP_start + rlmids[0].rmid_value.length();
          
                } else {
                        //Preprocess the read by using SSAHA:
                        int jj =0;
                        for(jj=0; jj<(int)rlmids.size() - 4; jj+=4 ) {
                                if(reads[i]->clip_found == true) break;
      
                                t_args rlargs0;
                                rlargs0.rec_id = i;
                                rlargs0.tid = jj;
                                rlargs0.reverse = false;
                                rlargs0.complement = false;
                                rlargs0.reverse_complement = false;
                                pthread_create( &rlthreads[jj], NULL, &tt_SSAHA_PCR, (void*)&rlargs0 );
      
                                t_args rlargs1;
                                rlargs1.rec_id = i;
                                rlargs1.tid = jj+1;
                                rlargs1.reverse = false;
                                rlargs1.complement = false;
                                rlargs1.reverse_complement = false;
                                pthread_create( &rlthreads[jj+1], NULL, &tt_SSAHA_PCR, (void*)&rlargs1 );
      
                                t_args rlargs2;
                                rlargs2.rec_id = i;
                                rlargs2.tid = jj+2;
                                rlargs2.reverse = false;
                                rlargs2.complement = false;
                                rlargs2.reverse_complement = false;
                                pthread_create( &rlthreads[jj+2], NULL, &tt_SSAHA_PCR, (void*)&rlargs2 );
      
                                t_args rlargs3;
                                rlargs3.rec_id = i;
                                rlargs3.tid = jj+3;
                                rlargs3.reverse = false;
                                rlargs3.complement = false;
                                rlargs3.reverse_complement = false;
                                pthread_create( &rlthreads[jj+3], NULL, &tt_SSAHA_PCR, (void*)&rlargs3 );
      
                                pthread_join(rlthreads[jj], NULL);
                                pthread_join(rlthreads[jj+1], NULL);
                                pthread_join(rlthreads[jj+2], NULL);
                                pthread_join(rlthreads[jj+3], NULL);
                        }
                
                        if(reads[i]->clip_found == false) { 
                                for(int t=0; t< (int)rlmids.size() - jj; t++) {
                                        t_args rlargs0;
                                        rlargs0.rec_id = i;
                                        rlargs0.tid = jj+1;
                                        rlargs0.reverse = false;
                                        rlargs0.complement = false;
                                        rlargs0.reverse_complement = false;
                                        pthread_create( &rlthreads[t], NULL, &tt_SSAHA_PCR, (void*)&rlargs0 );
               
                                        pthread_join(rlthreads[t], NULL);
                        }
           
                }
            }
        }
  
         delete izssaha;
        
          //Try Reverse Complement
  if( reads[i]->clip_found == false ) {
        //Reverse again (it is already a complement)
        reverse( rlstr.begin(), rlstr.end() );
        iz_SSAHA *izssaha = new iz_SSAHA();
        al_res = izssaha->Find( ref_str , rlstr );
        if( al_res.found_flag == true ) {
                scores = CalcScores(al_res.seq_1_al, al_res.seq_2_al, al_res.seq_1_al.length(), 0);
        
                if(scores.mismatches <= max_al_mism  ) {
            
                        //if(reads[i]->readID == "@HLIQV0O02IIS4V")
                        //        cout << reads[i]->rclip << endl;
                        reads[i]->rclip = al_res.pos + from;
                
                        if(reads[i]->rclip <= 16) {
                                reads[i]->rclip = reads[i]->read.length();
                        }
    
                        reads[i]->rev_bc = rlmids[reads[i]->rlmid.lmid_id].rmid_name;
      
                        reads[i]->clip_found = true;
                        reads[i]->revP_errors = scores.mismatches;
                
                        reads[i]->R_rclip = reads[i]->rclip + rlmids[0].rmid_value.length();
                
                        reads[i]->revP_start = al_res.pos + from; //cout << al_res.pos << " " << from << endl;
                        reads[i]->revP_end = reads[i]->revP_start + rlmids[0].rmid_value.length();
          
                } else {
            //Preprocess the read by using SSAHA:
        
                        int jj =0;
                        for(jj=0; jj<(int)rlmids.size() - 4; jj+=4 ) {
                        //cout << jj << endl;
                                if(reads[i]->clip_found == true) break;
      
                                t_args rlargs0;
                                rlargs0.rec_id = i;
                                rlargs0.tid = jj;
                                rlargs0.reverse = false;
                                rlargs0.complement = false;
                                rlargs0.reverse_complement = true;
                                pthread_create( &rlthreads[jj], NULL, &tt_SSAHA_PCR, (void*)&rlargs0 );
      
                                t_args rlargs1;
                                rlargs1.rec_id = i;
                                rlargs1.tid = jj+1;
                                rlargs1.reverse = false;
                                rlargs1.complement = false;
                                rlargs1.reverse_complement = true;
                                pthread_create( &rlthreads[jj+1], NULL,&tt_SSAHA_PCR, (void*)&rlargs1 );
      
                                t_args rlargs2;
                                rlargs2.rec_id = i;
                                rlargs2.tid = jj+2;
                                rlargs2.reverse = false;
                                rlargs2.complement = false;
                                rlargs2.reverse_complement = true;
                                pthread_create( &rlthreads[jj+2], NULL, &tt_SSAHA_PCR, (void*)&rlargs2 );
      
                                t_args rlargs3;
                                rlargs3.rec_id = i;
                                rlargs3.tid = jj+3;
                                rlargs3.reverse = false;
                                rlargs3.complement = false;
                                rlargs3.reverse_complement = true;
                                pthread_create( &rlthreads[jj+3], NULL, &tt_SSAHA_PCR, (void*)&rlargs3 );
      
                                pthread_join(rlthreads[jj], NULL);
                                pthread_join(rlthreads[jj+1], NULL);
                                pthread_join(rlthreads[jj+2], NULL);
                                pthread_join(rlthreads[jj+3], NULL);
                
                
    
                        }
                        if(reads[i]->clip_found == false) { 
                                for(int t=0; t< (int)rlmids.size() - jj; t++) {
                                        t_args rlargs0;
                                        rlargs0.rec_id = i;
                                        rlargs0.tid = jj+1;
                                        rlargs0.reverse = false;
                                        rlargs0.complement = false;
                                        rlargs0.reverse_complement = true;
                                        pthread_create( &rlthreads[t], NULL, &tt_SSAHA_PCR, (void*)&rlargs0 );
               
                                        pthread_join(rlthreads[t], NULL);
                                }
                        }
                }
        }
        delete izssaha;
        
   
  }
  
  
  
  
        al_res.seq_1_al.clear();
        al_res.seq_2_al.clear();
        rlstr.clear();
     
  
      
  
        
        if(reads[i]->discarded == 0) accept_counter++;
  
  } else {
      reads[i]->discarded = 1;
      discard_counter++;
  }
  
  
  line_counter++;
  if(line_counter%50000 == 0)
     cout << "Line No: " << line_counter << endl;
         
  pthread_exit(NULL);
    
}

static void *tt_SSAHA_PCR(void *targs) {
    t_args foo = *((t_args*)targs);
  
    iz_SSAHA *izssaha = new iz_SSAHA();
    string query_str = rlmids[foo.tid].rmid_value;// + "GGTCGGCGTCT";
    int read_len = reads[foo.rec_id]->read.length();
    int from; 
    if(read_len <= 100) {
        from = 0;
    } else if( (read_len <=300) && (read_len > 100) ) {
        from = read_len - read_len/2; //divide by 2
    } else if( (read_len <= 600) && (read_len > 300) ) {
        from = read_len - read_len/4; //divide by 4
    } else  {
        from = read_len - read_len/8; //divide by 8
    } 
    
    string ref_str = reads[foo.rec_id]->read.substr( from, read_len - from );
    
    if( foo.reverse == true ) {
        reverse( query_str.begin(), query_str.end() );
    }
    if( foo.complement == true ) {
        query_str =  MakeSeqComplement(query_str);
    }
    if( foo.reverse_complement == true ) {
        query_str =  MakeSeqComplement(query_str);
        reverse( query_str.begin(), query_str.end() );
    }
    
    //cout << from << " " << read_len << endl << reads[foo.rec_id]->readID << endl << ref_str << endl;
    if (reads[foo.rec_id]->clip_found == false) {
        AlignResult al_res = izssaha->Find( ref_str , query_str );
        if (al_res.found_flag == true) {
                AlignScores scores;
                scores = CalcScores(al_res.seq_1_al, al_res.seq_2_al, al_res.seq_1_al.length(), 0);
                if( scores.mismatches <= max_al_mism  ) {
            
                        if (reads[foo.rec_id]->clip_found == false) {
                        //pthread_spin_lock(&spinlock1);
                
                                reads[foo.rec_id]->rclip = al_res.pos + from;
                                if(reads[foo.rec_id]->rclip <= 16) {
                                        reads[foo.rec_id]->rclip = reads[foo.rec_id]->read.length();;
                                }
                                reads[foo.rec_id]->rev_bc = rlmids[foo.tid].rmid_name;//string(itoa(foo.tid+1, new char[5], 10));
                                //cout << rlmids[foo.tid].name << endl;
                                //cout << al_res.seq_1 << endl << al_res.seq_2 << endl;
                                reads[foo.rec_id]->clip_found = true;
                                reads[foo.rec_id]->revP_errors = scores.mismatches;
                                reads[foo.rec_id]->R_rclip = reads[foo.rec_id]->rclip + rlmids[0].rmid_value.length();
                        //pthread_spin_unlock(&spinlock1);
                                reads[foo.rec_id]->revP_start = al_res.pos + from; //cout << al_res.pos << " " << from << endl;
                                reads[foo.rec_id]->revP_end = reads[foo.rec_id]->revP_start + rlmids[0].rmid_value.length();
                        }
                } 
                if(reads[foo.rec_id]->clip_found == false) reads[foo.rec_id]->rclip_errors = scores.mismatches;
        
        }
        
        al_res.seq_1_al.clear();
        al_res.seq_2_al.clear();
        query_str.clear();
        ref_str.clear();
        delete izssaha;
    }
      
    pthread_exit(NULL);
}
