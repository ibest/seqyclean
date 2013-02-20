/* 
 * File:   Illumina.cpp
 * Author: ilya
 * 
 * Created on 20 Сентябрь 2012 г., 16:11
 */

#include "Illumina.h"

Illumina::Illumina() {
}

Illumina::Illumina(const Illumina& orig) {
}

Illumina::~Illumina() {
}


void Illumina::CleanSeq_PE() {
    cout << "# of reads in R1 file: " << reads_1.size() << endl << "# of reads in R2 file: " << reads_2.size() << endl;
    /*TruSeq adapters detecting*/
    /*First file, detecting the type of adapters it contains*/
    bool adapter1_found = false;
    bool adapter2_found = false;
    
    //for(it_reads_1 = reads_1.begin(); it_reads_1 != reads_1.end(); it_reads_1++) {
    for(int i=0; i<(int)reads_1.size(); i++) {
        /*First 20 bases of i5 adapter forward*/
        size_t found;
        string ts_adapter = tmpl_i5_1.substr(0,15);
        //string ts_adapter = tmpl_i5_1.substr(0,10);
        found = reads_1[i]->read.find( ts_adapter );
        if( found != string::npos ) {
            cout << found << " i5 adapter in forward" << endl; 
            adapter_type_R1 = "i5f";
            adapter1_found = true;
            query_str1 = ts_adapter;
            break;
        }
       
        /*First 20 bases of i5 adapter in reverse complement*/
        ts_adapter = MakeRevComplement(tmpl_i5_2).substr(0,15);
        //ts_adapter = MakeRevComplement(tmpl_i5_2).substr(0,10);
        found = reads_1[i]->read.find( ts_adapter );
        if( found != string::npos ) 
        {
            cout << found << " i5 adapter in reverse complement" << endl; 
            adapter_type_R1 = "i5r";
            adapter1_found = true;
            query_str1 = ts_adapter;
            break; 
        }
        
        /*First 20 bases of i7 adapter forward*/
        ts_adapter = tmpl_i7_1.substr(0,15);
        //ts_adapter = tmpl_i7_1.substr(0,10);
        found = reads_1[i]->read.find( ts_adapter );
        if( found != string::npos ) 
        {
            cout << found << " i7 adapter in forward" << endl; 
            adapter_type_R1 = "i7f";
            adapter1_found = true;
            query_str1 = ts_adapter;
            break;
        }
       
        /*First 20 bases of i5 adapter in reverse complement*/
        ts_adapter = MakeRevComplement(tmpl_i7_2).substr(0,15);
        //ts_adapter = MakeRevComplement(tmpl_i7_2).substr(0,10);
        found = reads_1[i]->read.find( ts_adapter );
        if( found != string::npos ) {
            cout << found << " i7 adapter in reverse complement" << endl; 
            adapter_type_R1 = "i7r";
            adapter1_found = true;
            query_str1 = ts_adapter;
            break; 
        }
        
    }
    
    /*For the second file*/
    //for(it_reads_2 = reads_2.begin(); it_reads_2 != reads_2.end(); it_reads_2++) {
    for(int i=0; i<(int)reads_2.size(); i++) {
        /*First 20 bases of i5 adapter forward*/
        size_t found;
        string ts_adapter = tmpl_i5_1.substr(0,15);
        //string ts_adapter = tmpl_i5_1.substr(0,10);
        found = reads_2[i]->read.find( ts_adapter );
        if( found != string::npos ) {
            cout << found << " i5 adapter in forward" << endl; 
            adapter_type_R2 = "i5f";
            adapter2_found = true;
            query_str2 = ts_adapter;
            break;
        }
       
        /*First 20 bases of i5 adapter in reverse complement*/
        ts_adapter = MakeRevComplement(tmpl_i5_2).substr(0,15);
        //ts_adapter = MakeRevComplement(tmpl_i5_2).substr(0,10);
        found = reads_2[i]->read.find( ts_adapter );
        if( found != string::npos ) {
            cout << found << " i5 adapter in reverse complement" << endl; 
            adapter_type_R2 = "i5r";
            adapter2_found = true;
            query_str2 = ts_adapter;
            break; 
        }
        
        /*First 20 bases of i7 adapter forward*/
        ts_adapter = tmpl_i7_1.substr(0,15);
        //ts_adapter = tmpl_i7_1.substr(0,10);
        found = reads_2[i]->read.find( ts_adapter );
        if( found != string::npos ) {
            cout << found << " i7 adapter in forward" << endl; 
            adapter_type_R2 = "i7f";
            adapter2_found = true;
            query_str2 = ts_adapter;
            break;
        }
       
        /*First 20 bases of i5 adapter in reverse complement*/
        ts_adapter = MakeRevComplement(tmpl_i7_2).substr(0,15);
        //ts_adapter = MakeRevComplement(tmpl_i7_2).substr(0,10);
        found = reads_2[i]->read.find( ts_adapter );
        if( found != string::npos ) {
            cout << found << " i7 adapter in reverse complement" << endl; 
            adapter_type_R2 = "i7r";
            adapter2_found = true;
            query_str2 = ts_adapter;
            break; 
        }
        
    }
    
    
    
    if( (adapter1_found == true) && (adapter2_found == true) ) 
    {
        query_str1 = query_str1.substr(0, 15);
        query_str2 = query_str2.substr(0, 15);
        
        /*Looking for TrueSec adapters in both of files with help of SSAHA algorithm*/
        //pthread_t threads[2];
        //pthread_create( &threads[0], NULL, &Illumina::StaticThreadProc, (void *)(1) );
        //pthread_create( &threads[1], NULL, &Illumina::StaticThreadProc, (void *)(2) );
        //pthread_join(threads[0], NULL);
        //pthread_join(threads[1], NULL);
        this->t_CleanSeq_PE();
    
    }
    else 
    {
        cout << "TruSeq adapters not found\n";
        //return 0;
        query_str1 = "";
        query_str2 = "" ; 
        this->t_CleanSeq_PE();
    }
    /*Making a report and output files*/
    
    
}


void Illumina::t_CleanSeq_PE() {
    
    pthread_t threads1[NUM_THREADS];
    pthread_t threads2[NUM_THREADS]; 
    int lim;
    lim = reads_1.size();
     for(int i = 0; i < lim - NUM_THREADS; i+=NUM_THREADS) {
         for(int t=0; t<NUM_THREADS; t++) {
             pthread_create( &threads1[t], NULL, &Illumina::tt_StaticThreadProc_1, (void *)(i+t) );
             pthread_create( &threads2[t], NULL, &Illumina::tt_StaticThreadProc_2, (void *)(i+t) );
         }
    
         for(int t=0; t<NUM_THREADS; t++) {
           pthread_join(threads1[t], NULL);
           pthread_join(threads2[t], NULL);
         }
      }
    
    
}


void* Illumina::tt_CleanSeq_1(void* targ) {
    
    long index;
    index = (long)targ;
    bool adapter_found = false;
    
    if(reads_1[index]->discarded == 1) 
    {
        pthread_exit(NULL); return 0;
    }
    
    if(vector_flag == true) 
    {
       CheckVector(reads_1[index]); 
       
       if( ( reads_1[index]->read.length() - (reads_1[index]->v_end - reads_1[index]->v_start) ) < 50 ) 
       {
          reads_1[index]->discarded = 1;
          reads_1[index]->discarded_by_vector = 1;
          pthread_exit(NULL); return 0;
       }
    }
    
    
    adapter_found = false;
    size_t found;
    found = reads_1[index]->read.rfind( query_str1 );
    if( found != string::npos ) 
    {
        adapter_found = true;
        reads_1[index]->tru_sec_pos = found;
        reads_1[index]->b_adapter = adapter_type_R1;
             
        reads_1[index]->lclip = 1;
        reads_1[index]->rclip = reads_1[index]->tru_sec_pos;
             
        if(reads_1[index]->rclip <= 50) 
        { 
            reads_1[index]->discarded = 1;  
            reads_1[index]->discarded_by_read_length = 1;
            pthread_exit(NULL);
            return 0;
        }
    } else 
    {
        //SSAHA job starts here
        iz_SSAHA *izssaha = new iz_SSAHA();
        AlignResult al_res = izssaha->Find( reads_1[index]->read , query_str1 );
        AlignScores scores;
        if( al_res.found_flag == true ) 
        {
            scores = CalcScores(al_res.seq_1_al, al_res.seq_2_al, al_res.seq_1_al.length(), 0);
            if(scores.mismatches <= max_al_mism  ) 
            {
               adapter_found = true;
               reads_1[index]->tru_sec_pos = al_res.pos;
               reads_1[index]->b_adapter = adapter_type_R1;// cout << "SSAHA " << reads_1[i]->tru_sec_pos << " " << reads_1[i]->b_adapter << endl;
               reads_1[index]->lclip = 1;
               reads_1[index]->rclip = reads_1[index]->tru_sec_pos;
                  
               if(reads_1[index]->rclip <= 50) 
               { 
                   reads_1[index]->discarded = 1;
                   reads_1[index]->discarded_by_read_length = 1;
                   delete izssaha;
                   pthread_exit(NULL); return 0;
               }
            }
         }
             
         if(adapter_found == false) reads_1[index]->tru_sec_pos = reads_1[index]->clear_length;
             
         reads_1[index]->lclip = 1;
         reads_1[index]->rclip = reads_1[index]->tru_sec_pos;
             
             
         if(reads_1[index]->rclip <= 50) 
         { 
            reads_1[index]->discarded = 1;
            reads_1[index]->discarded_by_read_length = 1;
            delete izssaha;
            pthread_exit(NULL); return 0;
         }
             
       delete izssaha;
    }
    
    pthread_exit(NULL); 
    return 0;
    
}

void* Illumina::tt_CleanSeq_2(void* targ) {
    long index;
    index = (long)targ;
    bool adapter_found = false;
   
    if(reads_2[index]->discarded == 1) {
        pthread_exit(NULL); return 0;
    }
    
    if(vector_flag == true) 
    {
       CheckVector(reads_2[index]); 
       
       if( ( reads_2[index]->read.length() - (reads_2[index]->v_end - reads_2[index]->v_start) ) < 50 ) 
       {
          reads_2[index]->discarded = 1;
          reads_2[index]->discarded_by_vector = 1;
          pthread_exit(NULL); return 0;
       }
    }
    
    adapter_found = false;
    size_t found = reads_2[index]->read.rfind( query_str2 );
    if( found != string::npos ) 
    {
        adapter_found = true;
        reads_2[index]->tru_sec_pos = found;
        reads_2[index]->b_adapter = adapter_type_R2;
        reads_2[index]->lclip = 1;
        reads_2[index]->rclip = reads_2[index]->tru_sec_pos;
               
        if(reads_2[index]->rclip < 50) 
        { 
           reads_2[index]->discarded = 1;
           reads_2[index]->discarded_by_read_length = 1;
           pthread_exit(NULL); return 0;
        }
             
        if(qual_trim_flag == true) 
        {
           if(reads_2[index]->lucy_rclip < reads_2[index]->rclip) 
           {
               reads_2[index]->read = reads_2[index]->read.substr(0 , reads_2[index]->lucy_rclip );
               reads_2[index]->illumina_quality_string = reads_2[index]->illumina_quality_string.substr(0,reads_2[index]->lucy_rclip) ; 
           } else 
           {
               reads_2[index]->read = reads_2[index]->read.substr(0 , reads_2[index]->rclip );
               reads_2[index]->illumina_quality_string = reads_2[index]->illumina_quality_string.substr(0,reads_2[index]->rclip) ; 
           }
        } 
        else 
        {
           reads_2[index]->read = reads_2[index]->read.substr(0 , reads_2[index]->rclip );
           reads_2[index]->illumina_quality_string = reads_2[index]->illumina_quality_string.substr(0,reads_2[index]->rclip) ; 
        }
    } else 
    {
        //SSAHA job starts here
        iz_SSAHA *izssaha = new iz_SSAHA();
        AlignResult al_res = izssaha->Find( reads_2[index]->read , query_str2 );
        AlignScores scores;
        if( al_res.found_flag == true ) 
        {
           
            scores = CalcScores(al_res.seq_1_al, al_res.seq_2_al, al_res.seq_1_al.length(), 0);
           if(scores.mismatches <= max_al_mism  ) 
           {
               adapter_found = true;
               reads_2[index]->tru_sec_pos = al_res.pos;
               reads_2[index]->b_adapter = adapter_type_R2;
                     
               reads_2[index]->lclip = 1;
               reads_2[index]->rclip = reads_2[index]->tru_sec_pos;
                     
                     
               //
             
               if(reads_2[index]->rclip < 50) 
               { 
                   delete izssaha;
                   reads_2[index]->discarded = 1;
                   reads_2[index]->discarded_by_read_length = 1;
                   pthread_exit(NULL); return 0;
               }
           }
        } 
        
        if(adapter_found == false) reads_2[index]->tru_sec_pos = reads_2[index]->clear_length;
        
        reads_2[index]->lclip = 1;
        reads_2[index]->rclip = reads_2[index]->tru_sec_pos;
               
        if(reads_2[index]->rclip < 50) 
        { 
          reads_2[index]->discarded = 1;
          reads_2[index]->discarded_by_read_length = 1;
          delete izssaha;
          pthread_exit(NULL); return 0;
        }
             
        delete izssaha;
               
    }
    
             
    pthread_exit(NULL); 
    return 0;
    
               
}


void Illumina::CleanSeq_SE() {
    cout << reads_1.size() << " " << reads_2.size() << endl;
    /*TruSeq adapters detecting*/
    /*First file, detecting the type of adapters it contains*/
    bool adapter1_found = false;
    
    //for(it_reads_1 = reads_1.begin(); it_reads_1 != reads_1.end(); it_reads_1++) {
    for(int i=0; i<(int)reads_1.size(); i++) {
        /*First 20 bases of i5 adapter forward*/
        size_t found;
        string ts_adapter = tmpl_i5_1.substr(0,20);
        found = reads_1[i]->read.find( ts_adapter );
        if( found != string::npos ) {
            cout << found << " i5 adapter in forward" << endl; 
            adapter_type_R1 = "i5f";
            adapter1_found = true;
            query_str1 = ts_adapter;
            break;
        }
       
        /*First 20 bases of i5 adapter in reverse complement*/
        ts_adapter = MakeRevComplement(tmpl_i5_2).substr(0,20);
        found = reads_1[i]->read.find( ts_adapter );
        if( found != string::npos ) {
            cout << found << " i5 adapter in reverse complement" << endl; 
            adapter_type_R1 = "i5r";
            adapter1_found = true;
            query_str1 = ts_adapter;
            break; 
        }
        
        /*First 20 bases of i7 adapter forward*/
        ts_adapter = tmpl_i7_1.substr(0,20);
        found = reads_1[i]->read.find( ts_adapter );
        if( found != string::npos ) {
            cout << found << " i7 adapter in forward" << endl; 
            adapter_type_R1 = "i7f";
            adapter1_found = true;
            query_str1 = ts_adapter;
            break;
        }
       
        /*First 20 bases of i5 adapter in reverse complement*/
        ts_adapter = MakeRevComplement(tmpl_i7_2).substr(0,20);
        found = reads_1[i]->read.find( ts_adapter );
        if( found != string::npos ) {
            cout << found << " i7 adapter in reverse complement" << endl; 
            adapter_type_R1 = "i7r";
            adapter1_found = true;
            query_str1 = ts_adapter;
            break; 
        }
        
    }
    
    
    
    if( adapter1_found == true ) 
    {
        /*Looking for TrueSec adapters in both of files with help of SSAHA algorithm*/
        this->t_CleanSeq_SE();
    }
    /*Making a report and output files*/
    
    
}

void Illumina::t_CleanSeq_SE() {
    
    pthread_t threads1[NUM_THREADS];
    int lim;
    lim = reads_1.size();
     for(int i = 0; i < lim - NUM_THREADS; i+=NUM_THREADS) {
         for(int t=0; t<NUM_THREADS; t++) {
             pthread_create( &threads1[t], NULL, &Illumina::tt_StaticThreadProc_1, (void *)(i+t) );
         }
    
         for(int t=0; t<NUM_THREADS; t++) {
           pthread_join(threads1[t], NULL);
         }
      }
    
    
}

int Illumina::CheckContaminants(string seq)
{
    //Sampling with some special frequency:
    string read = string(seq);
    int line_len = seq.length();
    stoupper(read);
    //cout <<  line_len << endl;
    bool good_flag = false; //false means this k_mer belongs to some known sequence, true - the opposite way
    bool stop_flag = false;
    
    good_flag = true;
    vector <HitData> matches;
    int i = 0;
    
    while( i< line_len ) 
    {
       
       string kmer_str;
       
       kmer_str = read.substr(i,KMER_SIZE);
       
       it_ContDict = ContDict.find(kmer_str);
       if(it_ContDict != ContDict.end()) 
       {
               /*Got hit:*/
           HitData t_match_struct;
           t_match_struct.pos = i;
           t_match_struct.kmers = (*it_ContDict).second;
           t_match_struct.k_mer_string = kmer_str;
           matches.push_back(t_match_struct);
                    
           if(matches.size() == 2) 
           {
                         
               stop_flag = true;
               good_flag = false;
               break;
                       
           }
                        
                        
       }
               
       if(stop_flag == true) break;   
              
       i = i+ KMER_SIZE+DISTANCE;   
   }
    
   //cout <<  i << endl;  
   matches.clear();
    
   return good_flag;
}



//Dynamic methods
void Illumina::t_CleanSeq_PE_Dynamic(Read &read) {
    
    pthread_t threads1[NUM_THREADS];
    pthread_t threads2[NUM_THREADS]; 
    int lim;
    lim = reads_1.size();
     for(int i = 0; i < lim - NUM_THREADS; i+=NUM_THREADS) {
         for(int t=0; t<NUM_THREADS; t++) {
             pthread_create( &threads1[t], NULL, &Illumina::tt_StaticThreadProc_1, (void *)(i+t) );
             pthread_create( &threads2[t], NULL, &Illumina::tt_StaticThreadProc_2, (void *)(i+t) );
         }
    
         for(int t=0; t<NUM_THREADS; t++) {
           pthread_join(threads1[t], NULL);
           pthread_join(threads2[t], NULL);
         }
      }
    
    
}