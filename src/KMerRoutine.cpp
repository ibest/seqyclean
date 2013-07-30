#include "KMerRoutine.h"

/*Holds distance differences*/
vector<int> dist_diff;


/*
 KMerRoutine:
 * This function is responsible for making decision whether it is a "good" of a "bad" sequence.
 * Parameters:
 *  string _seq -> input sequence that we want to check,
 *  int k_mer -> k_mer size (length of piece in bases),
 *  int freq -> sampling frequency, i.e. how far we place k_mers from each other,
 *  int num_attempts -> actually this is the required number of k_mers we want to put in the sequence to check,
 *  int num_of_matches -> number of k_mer matches in sequence "seq",
 *  int start_pos -> starting point (in bases) in checking sequence "seq" (default start_pos = 0).
 * int flag -> 1 : parsing contaminants, 0 : parsing vector
 */
int check_cont_call = 0;

int CheckContaminants(string seq)
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
       
       kmer_str = read.substr(i,KMER_SIZE_CONT);
       
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
              
       i = i+ KMER_SIZE_CONT+DISTANCE;   
   }
    
   matches.clear();
    
   return good_flag;
}

/*Vectors : parsing from the right end*/
int CheckVectorRight(Read &read) {
    
    /*
     * Short notes how it works:
     * 1. Obtain the sequence. Let's call it "seq"
     * 2. Starting from the position 0 (pos0) from that seq, take the piece with length equal to k_mer. Let's call this pease a "kmer_str".
     * 3. Look for the kmer_str in the dictionary Dict. 
     * 4. If it was successful (if we found that kmer_str in the Dict) -> go to 2.  but start with position pos1 = pos0+freq.
     * 5. If not (this can be considered as good sequence) -> put this seq to the file "good.txt".
     * 6. Otherwise, mar the seq as bad, put it into the "bad.txt" file.
     * 7. Obtain the next string.
     */
    string seq = read.read;
    int line_len = seq.length();
    bool good_flag = false; //false means this k_mer belongs to some known sequence, true - the opposite way
//    bool stop_flag = false;
    
    //Sampling with some special frequency:
    stoupper(seq);
    
    
    map<string, vector<k_mer_struct> >::iterator it_VectorDict2;
    //Serial implementation
    vector<HitData> matches;
    for( int i=line_len; i>KMER_SIZE; i-= DISTANCE ) 
    {
       /*Obtaining a new piece of KMER_SIZE:*/
       string kmer_str = seq.substr(i-KMER_SIZE, KMER_SIZE);
       it_VectorDict2 = VectorDict.find(kmer_str);
                
       if(it_VectorDict2 != VectorDict.end()) 
       { 
          /*Got hit: */
          /*Taking a piece from data read :
           * From i-(KMER_SIZE)
           * To (tail) : end of the read 
           */
           /*All of these variables are to be saved into HitData : */
          HitData hit_data;
          hit_data.kmers = (*it_VectorDict2).second;
          hit_data.k_mer_string = kmer_str;
          
          if (i-KMER_SIZE-10 < 5) continue;
          
          hit_data.string_to_align = seq.substr( i-KMER_SIZE-10 /*From*/, line_len - (i-KMER_SIZE-10) );
          hit_data.pos = i; /*Position where match occured*/
          matches.push_back(hit_data);
       } 
                
       if (i <= (KMER_SIZE+DISTANCE) ) 
       {
          int msize = matches.size();
                   
          if ( msize > 3 ) 
          {
             /*Doing a pairwise alignment : */
             /*Check if match positions of the first and last hits are in the range between [2*(KMER_SIZE+FREQ), (FREQ-KMER_SIZE)], i.e. approximately [50, line_len -25] */
             if( (matches[msize-1].pos <= 2*(KMER_SIZE+DISTANCE) /*~50*/ ) && (matches[0].pos >= line_len-25) /*~ (line_len - 25)*/ ) 
             {
                 /*Check if no long gaps between matches :*/
                 DistanceStruct dist_struct = CheckDistance(matches);
                 if ( dist_struct.flag == false ) 
                 {
                    //Completely throw away this read :  
                    good_flag = false;
                    break;
                 } 
                 else 
                 {
                    /*Sometimes, the part of vector can be found in one of the ends*/
                    if(dist_struct.i > 3) 
                    {
                       GetRClip(matches[dist_struct.i], read.readID, read);
                    } 
                    else 
                    {
                        read.r_vec_start = 0;
                        read.r_vec_end = 0;
                        read.vec_len = 0;
                        read.vec_err = 0;
                        read.vec_start_pos = 0;
                        read.vec_id = (char*)"NA";
                    }
                        
                    good_flag = true;
                    break;
                 }
             } else 
             {
                 /*If first match : */
                 if(matches[0].pos >= line_len-50) 
                 {
                    //Doing a pairwise alignment :
                    good_flag = true;
                                
                    DistanceStruct dist_struct = CheckDistance(matches);
                    if( dist_struct.flag == true ) 
                    {
                      /*Looking for the right clip point : */
                      GetRClip(matches[dist_struct.i], read.readID, read);
                                                //cout << sss << endl << matches[dist_struct.i].pos << endl;
                    } 
                    else 
                    {
                      GetRClip(matches[msize-1], read.readID, read);
                    }
                        
                    break;
                 } else 
                 {
                   good_flag = true;
                   read.r_vec_start = 0;
                   read.r_vec_end = 0;
                   read.vec_len = 0;
                   read.vec_err = 0;
                   read.vec_start_pos = 0;
                   read.vec_id = (char*)"NA";
                   break;
                 }    
             }
          } 
          else if( (msize > 0) && (msize<=3) ) 
          {
             if(matches[0].pos >= line_len-10) 
             {
                 GetRClip(matches[0], read.readID, read);
             } 
             else 
             {
                 read.r_vec_start = 0;
                 read.r_vec_end = 0;
                 read.vec_len = 0;
                 read.vec_err = 0;
                 read.vec_start_pos = 0;
                 read.vec_id = (char*)"NA";
             }
             
             good_flag = true;
             break;
          } else 
          {
             read.r_vec_start = 0;
             read.r_vec_end = 0;
             read.vec_len = 0;
             read.vec_err = 0;
             read.vec_start_pos = 0;
             read.vec_id = (char*)"NA";
             good_flag = true;
             break;
          }
       }
    }
                    
    return (good_flag == true) ? 1 : 0 ;
}

/*Vectors : parsing from the left end*/
int CheckVectorLeft(Read &read) 
{
    
    string seq = read.read;
    int line_len = seq.length();
    bool good_flag = false; //false means this k_mer belongs to some known sequence, true - the opposite way
    //bool stop_flag = false;
    
    map<string, vector<k_mer_struct> >::iterator it_VectorDict2;
    //Serial implementation
    vector<HitData> matches;

    for( int i=0; i<line_len-KMER_SIZE-10; i+=DISTANCE ) 
    {
        /*Obtaining a new piece of KMER_SIZE:*/
        string kmer_str = seq.substr(i,KMER_SIZE); //cout << "!!!\n";
        it_VectorDict2 = VectorDict.find(kmer_str);
        if(it_VectorDict2 != VectorDict.end()) 
        {
          /*Got hit: */
          /*Taking a piece from data read : */
          /*All of these variables are to be saved into HitData : */
          string tmpstr;
          /*To : */
          int to_pos =  i+KMER_SIZE+10;
          /*The whole string : */ //cout << seq << endl << to_pos << endl;
          tmpstr = seq.substr( 0, to_pos ); 
          HitData hit_data;
          hit_data.kmers = (*it_VectorDict2).second;
          hit_data.k_mer_string = kmer_str;
          hit_data.string_to_align = tmpstr;
          hit_data.pos = i;
          matches.push_back(hit_data);
        }
                
        if (i >= (line_len-KMER_SIZE-30)) 
        {
            int msize = matches.size();
            if(msize > 1) 
            {
               /*Doing a pairwise alignment : */
               if(matches[msize-1].pos >= (line_len-KMER_SIZE - 30) &&  (matches[0].pos <= (15 + KMER_SIZE)) ) 
               {
                    DistanceStruct dist_struct = CheckDistance(matches);
                    if ( dist_struct.flag == false ) 
                    {
                       //Completely throw away this read :  
                       good_flag = false;
                       break;
                    } else 
                    {
                       GetLClip(matches[dist_struct.i],read.readID, read);
                       good_flag = true;
                       break;
                    }
               } else 
               {
                  /*If first match : */
                  if(matches[0].pos <= (15 + KMER_SIZE) ) 
                  {
                     //Doing a pairwise alignment :
                     good_flag = true;
                     DistanceStruct dist_struct = CheckDistance(matches);
                     if( dist_struct.flag == true ) 
                     {
                        /*Looking for the right clip point : */
                        GetLClip(matches[dist_struct.i],read.readID, read);
                     } else 
                     {
                        GetLClip(matches[msize-1],read.readID, read);
                     }
                  } else 
                  {
                        good_flag = true;
                        if (read.l_vec_start > 10000)
                        {
                                read.l_vec_start = 0;
                                read.l_vec_end = 0;
                                read.vec_len = 0;
                                read.vec_err = 0;
                                read.vec_start_pos = 0;
                                read.vec_id = (char*)"NA";
                        }
                  }
               }
            } else 
            {
                if (read.l_vec_start > 10000)
                {
                        read.l_vec_start = 0;
                        read.l_vec_end = 0;
                        read.vec_len = 0;
                        read.vec_err = 0;
                        read.vec_start_pos = 0;
                        read.vec_id = (char*)"NA";
                }
                good_flag = true;
            }
            break;
        }
                
    }
    
    return (good_flag == true) ? 1 : 0 ;
}

DistanceStruct CheckDistance(vector <HitData> matches/*, int dir*/) {
    DistanceStruct dist_struct;
    dist_struct.i = 0;
    dist_struct.flag = false;
    dist_struct.pos = 0;
    for(int i=1; i<(int)matches.size(); i++) {
        if( abs(matches[i-1].pos - matches[i].pos) > 2*DISTANCE ) {
           dist_struct.flag = true;
           dist_struct.pos = matches[i-1].pos;
           dist_struct.i = i-1;
           break;
        } else {
            dist_struct.flag = false;
        }
    }
    
    return dist_struct;
}

bool CheckContig(vector <HitData> matches) 
{
    bool contig_flag = false;
    int tmp_err = 0;
    /*Calculating the difference in size between two hits:*/
    int tmp_dif = matches[0].kmers.size() - matches[1].kmers.size();
    /*---------------------------------------------------------------------------*/
    /*Now let us compare distances:*/
    if(tmp_dif < 0) {
         /*If  matches[0].kmers.size() < matches[1].kmers.size() : */  
         for(int ff=0; ff<(int)matches[0].kmers.size() ; ff++) {
              /*Look for similar sequence ids:*/
              if(matches[0].kmers[ff].seq_id == matches[1].kmers[ff].seq_id) {
                   /*Calculating difference in distance:*/
                   tmp_err = abs(abs(matches[0].kmers[ff].pos - matches[1].kmers[ff].pos) - abs(matches[0].pos - matches[1].pos));
                   dist_diff.push_back(tmp_err);
                   if(tmp_err <= 5) {
                        //cout << "diff pos in dictionary: "<< abs(matches[0].kmers[ff].pos - matches[1].kmers[ff].pos) << ", diff pos in read: " << abs(matches[0].pos - matches[1].pos) << endl;
                      // cout << matches[0].kmers[ff].pos - matches[1].kmers[ff].pos << endl;
                      // cout << tmp_err << endl; 
                       contig_flag = true;
                        break;
                   } 
                                             
                   contig_flag = false;
                                             
               }
         }
   } else {
         /*If  matches[0].kmers.size() > matches[1].kmers.size() : */  
         for(int ff=0; ff<(int)matches[1].kmers.size() ; ff++) {
             /*Look for similar sequence ids:*/
             if(matches[0].kmers[ff].seq_id == matches[1].kmers[ff].seq_id) {
                  /*Calculating difference in distance:*/
                  tmp_err = abs(abs(matches[0].kmers[ff].pos - matches[1].kmers[ff].pos) - abs(matches[0].pos - matches[1].pos));
                  dist_diff.push_back(tmp_err);
                  if(tmp_err <= 5) {
                      
                  //    cout << matches[0].kmers[ff].pos - matches[1].kmers[ff].pos << endl;
                  //    cout << tmp_err << endl;
                    // cout << "diff pos in dictionary: "<< abs(matches[0].kmers[ff].pos - matches[1].kmers[ff].pos) << ", diff pos in read: " << abs(matches[0].pos - matches[1].pos) << endl;  
                     contig_flag = true;
                     break;
                  } 
                                            
                  contig_flag = false;
                                       
             }
         }  
  }
    
    return contig_flag;
}

/*Input variables : */
/*
 structure that holds:
 * tail - the length of that piece
 * tmostr - a piece from read
 * k_mer_struct - hold a position and sequence id
*/
void GetRClip(HitData &hit_data, string sss, Read &read) {
    string vector_string_to_align;
    string string_to_align = hit_data.string_to_align;
    
    map<long /*seq_id*/, string /*sequence*/ >::iterator it_VectorSeqs2;
    
    if(hit_data.kmers.size() > 1) 
    {
        
        
        for(int s=0; s< (int)hit_data.kmers.size(); s++) 
        {
           /*In the vector : */
           int vpos = hit_data.kmers[s].pos; /*Position where match occured in the Vector string*/
           if(vpos < 10) 
           {
              string_to_align = hit_data.string_to_align.substr(vpos, hit_data.string_to_align.length() - hit_data.kmers[s].pos ) ;
           } else 
           {
              vpos = vpos - 10;
           }
              
              
           it_VectorSeqs2 = VectorSeqs.find(hit_data.kmers[s].seq_id);
           if(it_VectorSeqs2 != VectorSeqs.end()) 
           {
            
                if(string_to_align.length() > ((*it_VectorSeqs2).second.length() - vpos)) {
                //if(string_to_align.length() > ((*it_VectorSeqs).second.length() - vpos)) {
                     //vector_string_to_align = (*it_VectorSeqs).second.substr(vpos, (*it_VectorSeqs).second.length() - vpos );
                     vector_string_to_align = (*it_VectorSeqs2).second.substr(vpos, (*it_VectorSeqs2).second.length() - vpos );
                     string_to_align = string_to_align.substr(0, vector_string_to_align.length() );
                } else {
                     //vector_string_to_align = (*it_VectorSeqs).second.substr(vpos, string_to_align.length() );
                     vector_string_to_align = (*it_VectorSeqs2).second.substr(vpos, string_to_align.length() );
                }
              
                AlignResult al_res;
                al_res = banded(string_to_align, vector_string_to_align, 2, 1);
                          
                          //cout << al_res.n_mismatches << endl;
                          //cout << al_res.seq_1_al << endl;
                          //cout << al_res.seq_2_al << endl;
                          
              /*Calculating number of mismatches : */
              AlignScores scores;
              scores = CalcScores(al_res.seq_1_al,al_res.seq_2_al, 20, 0);
                          
              if(scores.mismatches <= 3) {
                  
                //read.rclip = hit_data.pos - (KMER_SIZE + 10);
                read.r_vec_start = hit_data.pos - (KMER_SIZE + 10);
                read.r_vec_end = read.rlmid.rmid_start;
                read.vec_start_pos = vpos;
                
                read.vec_err = scores.mismatches;
                
                read.vec_id = itoa(hit_data.kmers[s].seq_id, new char[8], 10) ;
                //cout << hit_data.kmers[s].seq_id << endl;
              } else {
                  
                //read.rclip = hit_data.pos - KMER_SIZE;
                read.r_vec_start = hit_data.pos - KMER_SIZE;
                
                if( read.rlmid.rmid_start == 0 ) {
                    read.r_vec_end = read.clear_length;
                } else {
                    read.r_vec_end = read.rlmid.rmid_start;
                }
                
                read.vec_start_pos = vpos;
                
                read.vec_err = scores.mismatches;
                
                read.vec_id = itoa(hit_data.kmers[s].seq_id, new char[8], 10) ;
                //cout << hit_data.kmers[s].seq_id << endl;
                //cout << read.readID << " " << read.vec_end << endl;
              }
              
              al_res.seq_1_al.erase();
              al_res.seq_2_al.erase();
           } 
      }  
    } else {
                 
        /*The size is 1 : */
        /*In the vector : */
        int vpos = hit_data.kmers[0].pos; /*Position where match occured in the Vector string*/
                        
        if(vpos < 10) {
            string_to_align = (string_to_align).substr(vpos, (string_to_align).length() - hit_data.kmers[0].pos );
        } else {
            vpos = vpos - 10;
        }
                        
                               
        //it_VectorSeqs = VectorSeqs.find(hit_data.kmers[0].seq_id);
        it_VectorSeqs2 = VectorSeqs.find(hit_data.kmers[0].seq_id);
        //if(it_VectorSeqs != VectorSeqs.end()) {
        if(it_VectorSeqs2 != VectorSeqs.end()) 
        {
            //if(string_to_align.length() > ((*it_VectorSeqs).second.length() - vpos)) {
            if((string_to_align).length() > ((*it_VectorSeqs2).second.length() - vpos)) {
                 //vector_string_to_align = (*it_VectorSeqs).second.substr(vpos, (*it_VectorSeqs).second.length() - vpos );
                 vector_string_to_align = (*it_VectorSeqs2).second.substr(vpos, (*it_VectorSeqs2).second.length() - vpos );
                 string_to_align = string_to_align.substr(0, vector_string_to_align.length() );
            } else {
                 //vector_string_to_align = (*it_VectorSeqs).second.substr(vpos, string_to_align.length() );
                 vector_string_to_align = (*it_VectorSeqs2).second.substr(vpos, string_to_align.length() );
            }
            
            AlignResult al_res;
            al_res = banded(string_to_align, vector_string_to_align, 2, 1);
                          
            /*Calculating number of mismatches : */
            AlignScores scores;
            scores = CalcScores(al_res.seq_1_al,al_res.seq_2_al, 20, 0);
         //   cout << sss << endl;
         //     cout << al_res.seq_1_al << endl;
         //     cout << al_res.seq_2_al << endl;
                          
            if(scores.mismatches <= 3) {
                  
                //rclip = hit_data.pos - 10;
                //read.rclip = hit_data.pos - (KMER_SIZE + 10);
                read.r_vec_start = hit_data.pos - (KMER_SIZE + 10);
                
                
                if( read.rlmid.rmid_start == 0 ) {
                    read.r_vec_end = read.clear_length;
                } else {
                    read.r_vec_end = read.rlmid.rmid_start;
                }
                
                read.vec_start_pos = vpos;
                
                read.vec_err = scores.mismatches;
                
                read.vec_id = itoa(hit_data.kmers[0].seq_id, new char[8], 10) ;
                //cout << hit_data.kmers[0].seq_id << endl;
            //    if(sss == "@ARTIF97") {
            //      cout << sss << endl;
           //       cout << al_res.seq_1_al << endl;
           //       cout << al_res.seq_2_al << endl;
          //        cout << hit_data.pos << endl;
          //      }
            } else {
                  
                //read.rclip = hit_data.pos - KMER_SIZE;
                read.r_vec_start = hit_data.pos - KMER_SIZE;
                
                if( read.rlmid.rmid_start == 0 ) {
                    read.r_vec_end = read.clear_length;
                } else {
                    read.r_vec_end = read.rlmid.rmid_start;
                }
                
                read.vec_err = scores.mismatches;
                
                read.vec_start_pos = vpos;
                
                read.vec_id = itoa(hit_data.kmers[0].seq_id, new char[8], 10) ;
                //cout << hit_data.kmers[0].seq_id << endl;
            //    if(sss == "@ARTIF97") {
            //      cout << sss << endl;
             //     cout << al_res.seq_1_al << endl;
             //     cout << al_res.seq_2_al << endl;
             //     cout << hit_data.pos << endl;
             //   }
                
            }
            al_res.seq_1_al.clear();
            al_res.seq_2_al.clear();
       } 
                        
    }
    
    string_to_align.clear();
    vector_string_to_align.clear();
    //delete string_to_align; delete vector_string_to_align;
    
     
}

void GetLClip(HitData &hit_data, string sss, Read &read) {
    
    map<long /*seq_id*/, string /*sequence*/ >::iterator it_VectorSeqs2;
    string tmpstr = hit_data.string_to_align; /*The piece of read*/
    int tail = hit_data.string_to_align.length(); /*The whole length of piece of read*/
    int i = hit_data.pos; /*Position where match occured*/
    
    if(hit_data.kmers.size() > 1) {
           for(int s=0; s< (int)hit_data.kmers.size(); s++) 
           {
              //string seqid = hit_data.kmers[s].seq_id;
              long seqid = hit_data.kmers[s].seq_id;
              int kpos = hit_data.kmers[s].pos; /*Position of that k-mer in the dictionary*/
              kpos+= KMER_SIZE + 10;
               //cout << seqid << endl;
              string vector_string;    
              //it_VectorSeqs = VectorSeqs.find(seqid);
              
              
              it_VectorSeqs2 = VectorSeqs.find(seqid);
              //if(it_VectorSeqs != VectorSeqs.end()) {
              if(it_VectorSeqs2 != VectorSeqs.end()) {
                                
                   if(tail < kpos) {
                                 //   cout << "b1" << endl;
                       //vector_string = (*it_VectorSeqs).second.substr(kpos-tail, tail );
                       vector_string = (*it_VectorSeqs2).second.substr(kpos-tail, tail );
                   } else {
                       //vector_string = (*it_VectorSeqs).second.substr(0, tail );
                       vector_string = (*it_VectorSeqs2).second.substr(0, tail );
                   }
                                //cout << tmpstr << endl;
                                //cout << vector_string << endl;
                                /*Performing banded alignment in order to find the best rclip : */
                   AlignResult al_res;
                                /*Checking for lengths of the sequences : */
                   int vector_string_len = vector_string.length();
                   int tmpstr_len = tmpstr.length();
                          
                   if(vector_string_len == tmpstr_len) {
                                        /*Lengths are equal : */
                        al_res = banded(tmpstr,vector_string,2,1);
                                 
                   } else {
                        if(vector_string_len > tmpstr_len) {
                            string s1 = vector_string.substr(0,tmpstr_len);
                            al_res = banded(tmpstr,s1,2,1);
                            s1.erase();
                        } else {
                            string s1 = tmpstr.substr(0,vector_string_len);
                            al_res = banded(s1,vector_string,2,1);
                            s1.erase();
                        }
                   }
                                
                                /*Calculating number of mismatches : */
                   AlignScores scores;
                   scores = CalcScores(al_res.seq_1_al,al_res.seq_2_al, 10, 1);
                   
                //   cout << sss << endl;
                //   cout << al_res.seq_1_al << endl;
                //   cout << al_res.seq_2_al << endl;
                   
                   if(scores.mismatches <= 1) {
              //                      cout << sss << endl;
              //            cout << hit_data.k_mer_string << endl;
              //     cout << al_res.seq_1_al << endl;
             //      cout << al_res.seq_2_al << endl;
             //      cout << tmpstr << endl << vector_string << endl;
                        //read.lclip = i+KMER_SIZE+10;
                        read.l_vec_end = i+KMER_SIZE+10;
                        read.l_vec_start = read.rlmid.lmid_end;
                        read.vec_start_pos = kpos;
                        read.vec_id = itoa(hit_data.kmers[s].seq_id, new char[8], 10) ;
                        //read.lclip+=1;
            //            cout << lclip << endl;   
                   } else {
                       
                      read.vec_start_pos = kpos;
                                    //good_flag = true;
                        //read.lclip = i+KMER_SIZE ;
                       read.l_vec_end = i+KMER_SIZE;
                       read.l_vec_start = read.rlmid.lmid_end;
                       read.vec_id = itoa(hit_data.kmers[s].seq_id, new char[8], 10) ; 
                   }
                   al_res.seq_1_al.clear();
                   al_res.seq_2_al.clear();
                                
              }
              
              vector_string.clear();
          }  
    } else {
                        
          /*The size is 1 : */
           //string seqid = hit_data.kmers[0].seq_id;
           long seqid = hit_data.kmers[0].seq_id;
           int kpos = hit_data.kmers[0].pos; /*Position of that k-mer in the dictionary*/
           kpos+= KMER_SIZE + 10;
                        
           
                                
           string vector_string;    
           //it_VectorSeqs = VectorSeqs.find(seqid);
           it_VectorSeqs2 = VectorSeqs.find(seqid);
           //if(it_VectorSeqs != VectorSeqs.end()) {
           if(it_VectorSeqs2 != VectorSeqs.end()) {
              
               if(tail < kpos) {
                    //vector_string = (*it_VectorSeqs).second.substr(kpos-tail, tail );
                    vector_string = (*it_VectorSeqs2).second.substr(kpos-tail, tail );
               } else {
                    //vector_string = (*it_VectorSeqs).second.substr(0, tail );
                    vector_string = (*it_VectorSeqs2).second.substr(0, tail );
               }
               
             //  cout << "kmer : " << endl;
             //  cout << hit_data.k_mer_string << endl;
             //             cout << "tmpstr : " << endl;
             //             cout << tmpstr << endl;
             //             cout << "Vector_string :" << endl;
             //             cout << vector_string << endl;
                          /*Performing banded alignment in order to find the best rclip : */
               AlignResult al_res;
                          /*Checking for lengths of the sequences : */
               int vector_string_len = vector_string.length();
               int tmpstr_len = tmpstr.length();
                          
               if(vector_string_len == tmpstr_len) {
                    /*Lengths are equal : */
                    al_res = banded(tmpstr,vector_string,2,1);
                                 
               } else {
                    if(vector_string_len > tmpstr_len) {
                        string s1 = vector_string.substr(0,tmpstr_len);
                        al_res = banded(tmpstr,s1,2,1);
                        s1.clear();
                    } else {
                        string s1 = tmpstr.substr(0,vector_string_len);
                        al_res = banded(s1,vector_string,2,1);
                        s1.clear();
                    }
               }
                          
                          //cout << al_res.n_mismatches << endl;
                 //         cout << al_res.seq_1_al << endl;
                //          cout << al_res.seq_2_al << endl;
                          
               /*Calculating number of mismatches : */
               AlignScores scores;
               scores = CalcScores(al_res.seq_1_al,al_res.seq_2_al, 10, 1);
                          
              // cout << sss << endl;
              //     cout << al_res.seq_1_al << endl;
              //     cout << al_res.seq_2_al << endl;
               
               
                if(scores.mismatches <= 1) {
                                   
                   read.vec_start_pos = kpos;
                   
                        //read.lclip = i+KMER_SIZE+10;
                        read.l_vec_end = i+KMER_SIZE+10;
                        read.l_vec_start = read.rlmid.lmid_end;
                       read.vec_id = itoa(hit_data.kmers[0].seq_id, new char[8], 10) ;
                } else {
               //           cout << sss << endl;
               //           cout << hit_data.k_mer_string << endl;
               //    cout << al_res.seq_1_al << endl;
               //    cout << al_res.seq_2_al << endl;
               //    cout << tmpstr << endl << vector_string << endl;
                   read.vec_start_pos = kpos;
                        //read.lclip = i+KMER_SIZE ;
                        read.l_vec_end = i+KMER_SIZE;
                        read.l_vec_start = read.rlmid.lmid_end;
                  read.vec_id = itoa(hit_data.kmers[0].seq_id, new char[8], 10) ;
                }
               
                al_res.seq_1_al.clear();
                al_res.seq_2_al.clear();
           }
           
           vector_string.clear();
                        
    }
    //cout << read.readID << endl;
    
    
    tmpstr.clear();
}


/*Vectors : parsing from the right end*/
int CheckVector(Read* read) {
    
    /*
     * Short notes how it works:
     * 1. Obtain the sequence. Let's call it "seq"
     * 2. Starting from the position 0 (pos0) from that seq, take the piece with length equal to k_mer. Let's call this pease a "kmer_str".
     * 3. Look for the kmer_str in the dictionary Dict. 
     * 4. If it was successful (if we found that kmer_str in the Dict) -> go to 2.  but start with position pos1 = pos0+freq.
     * 5. If not (this can be considered as good sequence) -> put this seq to the file "good.txt".
     * 6. Otherwise, mar the seq as bad, put it into the "bad.txt" file.
     * 7. Obtain the next string.
     */
    map<string, vector<k_mer_struct> >::iterator it_VectorDict2;
    string seq = read->read;
    stoupper(seq);
    vector<int> matches;
    int p = 0;
    for( int i=0; i < (int)seq.length() - R_limit*KMER_SIZE; i+=DISTANCE ) 
    {
       //Obtaining a new piece of KMER_SIZE:
       it_VectorDict2 = VectorDict.find( seq.substr(i, KMER_SIZE) );
                
       if(it_VectorDict2 != VectorDict.end()) 
       { 
          if (matches.size() > 1)
          {
              if ( abs(matches[matches.size()-1] - matches[matches.size()-2]) > allowable_distance ) 
              {
                  p+= 1;
                  if( p > pmax ) break;
                  
              }
          }
          
          matches.push_back(i);
          
       } 
    }            
       
    if ( matches.size() > 1 ) 
    {
        read->v_start = matches[0]; //Vector start
        if( read->v_start > vml )
            read->v_start -= vml;
        
        read->v_end = matches[matches.size()-1] + KMER_SIZE; //Vector end
        if( (read->v_end + vmr) < (int)read->read.length() )
            read->v_end += vmr;
            
        read->vector_found = 1;
    } else 
    {
        read->v_start = -1;
        read->v_end = -1;
        
    }
                    
    return 0;
}


//Vectors : parsing from the left end
int CheckVectorLeftIllumina(Read &read) 
{
    
    string seq = read.read;
    int line_len = seq.length();
    bool good_flag = false; //false means this k_mer belongs to some known sequence, true - the opposite way
    //bool stop_flag = false;
    
    map<string, vector<k_mer_struct> >::iterator it_VectorDict2;
    //Serial implementation
    vector<HitData> matches;

    for( int i=0; i<line_len-R_limit*KMER_SIZE; i+=DISTANCE ) 
    {
        /*Obtaining a new piece of KMER_SIZE:*/
        string kmer_str = seq.substr(i,KMER_SIZE); //cout << "!!!\n";
        it_VectorDict2 = VectorDict.find(kmer_str);
        if(it_VectorDict2 != VectorDict.end()) 
        {
          /*Got hit: */
          /*Taking a piece from data read : */
          /*All of these variables are to be saved into HitData : */
          string tmpstr;
          /*To : */
          int to_pos =  i+KMER_SIZE+10;
          /*The whole string : */ //cout << seq << endl << to_pos << endl;
          tmpstr = seq.substr( 0, to_pos ); 
          HitData hit_data;
          hit_data.kmers = (*it_VectorDict2).second;
          hit_data.k_mer_string = kmer_str;
          hit_data.string_to_align = tmpstr;
          hit_data.pos = i;
          
          if (matches.size() > 1)
          {
              if ( abs(matches[matches.size()-1].pos - matches[matches.size()-2].pos) > allowable_distance ) 
              {
                  break;
              }
          }
          
          matches.push_back(hit_data);
          
        }
    }           
        
    if(matches.size() > 1) 
    {
       //Doing a pairwise alignment : 
       if( ( matches[matches.size()-1].pos >= (line_len-(R_limit+1)*KMER_SIZE) ) &&  (matches[0].pos <= KMER_SIZE) ) 
       {
          good_flag = false; return good_flag;
       } else 
       {
          if(matches[0].pos <= KMER_SIZE ) 
          {
              read.l_vec_start = matches[0].pos;
              read.l_vec_end = matches[matches.size()-1].pos + KMER_SIZE+vml;
              read.vec_len = 0;
              read.vec_err = 0;
              read.vec_start_pos = 0;
              read.vec_id = (char*)"NA";
              good_flag = true; return good_flag;
          } 
          else
          {
              read.l_vec_start = 0;
              read.l_vec_end = 0;
              read.vec_len = 0;
              read.vec_err = 0;
              read.vec_start_pos = 0;
              read.vec_id = (char*)"NA";
              good_flag = true; return good_flag;
              
          }
        }
    }
    else if(matches.size() == 1)
    {
        if(matches[0].pos <= KMER_SIZE ) 
        {
                read.l_vec_start = matches[0].pos;
                read.l_vec_end = matches[0].pos + KMER_SIZE+vml;
                read.vec_len = 0;
                read.vec_err = 0;
                read.vec_start_pos = 0;
                read.vec_id = (char*)"NA";
                good_flag = true; return good_flag;
        }
        else
        {
                read.l_vec_start = 0;
                read.l_vec_end = 0;
                read.vec_len = 0;
                read.vec_err = 0;
                read.vec_start_pos = 0;
                read.vec_id = (char*)"NA";
                good_flag = true; return good_flag;
        }
    }
    else
    {
        read.l_vec_start = 0;
        read.l_vec_end = 0;
        read.vec_len = 0;
        read.vec_err = 0;
        read.vec_start_pos = 0;
        read.vec_id = (char*)"NA";
        good_flag = true; return good_flag;
    }
    
                
    return (good_flag == true) ? 1 : 0 ;
}
