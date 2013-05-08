#include "Dictionary.h"

/*Builds a new Vector Dictionary. Here it assumes that frequency is 1*/
int BuildVectorDictionary(string filename) {
    string str;
    
    std::fstream infile;
    /*Open given file:*/
    
    infile.open(filename.c_str(), std::fstream::in);
    
    cout << "Parsing screening file " << filename << endl;
    sum_stat << "Parsing screening file " << filename << endl;
    
    //Length of the read:
    int line_len;
    
    //Counts the number of lines (reads) in parsed file:
    long line_cnt = 0;
    
    //Record (read) ID
    //string rec_id;
    long vector_id = 0;
    //This holds a set of t_args - a structures to pass via threaded function:
    //vector<t_args> args_t;
    
    string vector_read;
    string vector_actual_id;
    //Loop thru all lines in input file:
    while ( getline(infile, str) ) {
        
        line_len = str.length();
        if (line_len <= 2) {
            continue;
        }
        
        //Determining a record (read) ID:
        remove( str.begin(), str.end(), ' ' );
        if(str[0]== '>') 
        {
            if( vector_read.length() > 15) 
            {
                VectorSeqs.insert(std::pair<long, string >(vector_id, vector_read));
                
                ParseString(vector_read, vector_id);
                
                vector_names.insert( std::pair<long, string >(vector_id, vector_actual_id) );
                
                vector_actual_id = str.substr(1,str.length()-1);
                //cout << vector_actual_id << endl;
                
                vector_id++;
                
                line_cnt = 0;
                vector_read.clear();
                
                
            }
            
            continue;
        }
        
        //To uppercase:
        stoupper(str);
        
        //Getting the whole read : 
        vector_read.append(string(str));
        
        line_cnt+=line_len;
        str.clear();
    }
    

    
    
    //The last one : 
    VectorSeqs.insert(std::pair<long, string >(vector_id, vector_read));
    ParseString(vector_read, vector_id);
    
    vector_names.insert( std::pair<long, string >(vector_id, vector_actual_id) );
    
    vector_read.clear();
    infile.close();
    
    //free(vector_read);
    //cout << VectorDict.size() << endl;
    return 0;
}

void ParseString(string str, int rec_id)
{
    
    int line_len = str.length();
    /*k_mer_struct - the structure wich holds information about k-mers. This structure is saved in the dictionary.*/
    k_mer_struct k_struct;
    vector<k_mer_struct> rec_id_set;
        
    for(int w=0; w< line_len-KMER_SIZE; w++) {
    
        string seq0;
   
        if((w+KMER_SIZE) < line_len) {
            seq0 = str.substr(w,KMER_SIZE);
        } else {
            seq0 = str.substr(w, line_len-w);
        }
    
         /*Making complement : */
         string seq_complement = MakeSeqComplement(seq0);
                
                
         /*Sequence (read) ID:*/
         k_struct.seq_id = rec_id;
         /*Calculating a position in the read:*/
         k_struct.pos = w;
      
      
         map<string, vector<k_mer_struct> >::iterator it_VectorDict = VectorDict.find(seq0);
    
         if(it_VectorDict == VectorDict.end()) {
              rec_id_set.push_back(k_struct);
              VectorDict.insert(std::pair<string, vector<k_mer_struct> >(seq0, rec_id_set));
              rec_id_set.clear();
         } else {
              (*it_VectorDict).second.push_back(k_struct);
         }
         
         
         it_VectorDict = VectorDict.find(seq_complement);
    
         if(it_VectorDict == VectorDict.end()) {
             rec_id_set.push_back(k_struct);
             VectorDict.insert(std::pair<string, vector<k_mer_struct> >(seq_complement, rec_id_set));
             rec_id_set.clear();
         } else {
             (*it_VectorDict).second.push_back(k_struct);
         }
         
         /*Making reverse complement*/
         reverse( seq_complement.begin(), seq_complement.end() );
         it_VectorDict = VectorDict.find(seq_complement);
    
         if(it_VectorDict == VectorDict.end()) {
             rec_id_set.push_back(k_struct);
             VectorDict.insert(std::pair<string, vector<k_mer_struct> >(seq_complement, rec_id_set));
             rec_id_set.clear();
         } else {
             (*it_VectorDict).second.push_back(k_struct);
         }
         
    }
    
    
}

/*Builds a new Contaminants Dictionary. Here it assumes that frequency is 1*/
int BuildContDictionary(string filename) {
    
    std::fstream infile;
    /*Open given file:*/
    infile.open(filename.c_str(),std::fstream::in);
    std::string str;
    cout << "Parsing screening file " << filename << endl;
    sum_stat << "Parsing screening file " << filename << endl;
    
    /*Length of the read:*/
    int line_len;
    
    /*Counts the number of lines (reads) in parsed file:*/
    long line_cnt = 0;
    
    /*Record (read) ID*/
    long contaminant_id = 0;
    
    string contaminant_read;
    /*Loop thru all lines in input file:*/
    while ( getline(infile, str) ) {
        
        line_len = str.length();
        if (line_len == 0) {
            continue;
        }
        
        /*Determining a record (read) ID:*/
        remove( str.begin(), str.end(), ' ' );
        if(str.substr(0,1)== ">") {
            
            if( contaminant_read.length() > 0) 
            {
               PutContKmer(contaminant_read, contaminant_id);
               contaminant_id++;
               line_cnt = 0;
               contaminant_read.clear();
                
            }
            continue;
        }
        
        /*To uppercase:*/
        stoupper(str);
        
        /*Getting the whole read : */
        contaminant_read+=str;
        
        line_cnt+=line_len;
        
    }
    
    infile.close();
    
    /*The last one : */
    PutContKmer(contaminant_read, contaminant_id);
    
    contaminant_read.clear();
    
    return 0;
}


void PutContKmer(string str, int rec_id)
{
  
  int line_len = str.length();
  /*k_mer_struct - the structure wich holds information about k-mers. This structure is saved in the dictionary.*/
  k_mer_struct k_struct;
  vector<k_mer_struct> rec_id_set;
        
  for(int w=0; w< line_len-KMER_SIZE_CONT; w++) {
    
     string seq0;
   
     if((w+KMER_SIZE_CONT) < line_len) 
     {
       seq0 = str.substr(w,KMER_SIZE_CONT);
     } else 
     {
       seq0 = str.substr(w, line_len-w);
     }
    
     /*Making complement : */
     string seq_complement = MakeSeqComplement(seq0);
                
                
     /*Sequence (read) ID:*/
     k_struct.seq_id = rec_id;
     /*Calculating a position in the read:*/
     k_struct.pos = w;
      
      
     map<string, vector<k_mer_struct> >::iterator it_ContDict = ContDict.find(seq0);
    
     if(it_ContDict == ContDict.end()) 
     {
         rec_id_set.push_back(k_struct);
         ContDict.insert(std::pair<string, vector<k_mer_struct> >(seq0, rec_id_set));
         rec_id_set.clear();
     } else 
     {
         (*it_ContDict).second.push_back(k_struct);
     }
         
         
     it_ContDict = ContDict.find(seq_complement);
    
     if(it_ContDict == ContDict.end()) 
     {
        rec_id_set.push_back(k_struct);
        ContDict.insert(std::pair<string, vector<k_mer_struct> >(seq_complement, rec_id_set));
        rec_id_set.clear();
     } else 
     {
        (*it_ContDict).second.push_back(k_struct);
     }
         
     /*Making reverse complement*/
     reverse( seq_complement.begin(), seq_complement.end() );
     it_ContDict = ContDict.find(seq_complement);
    
     if(it_ContDict == ContDict.end()) 
     {
        rec_id_set.push_back(k_struct);
        ContDict.insert(std::pair<string, vector<k_mer_struct> >(seq_complement, rec_id_set));
        rec_id_set.clear();
     } else 
     {
        (*it_ContDict).second.push_back(k_struct);
     }
         
    }
    
  
}

void Build_RLMIDS_Dictionary(char* rlmids_file) {
    rlmids.erase(rlmids.begin(), rlmids.end());
    std::ifstream infile;
    
    cout << "Parsing custom RL MIDS file ... " << rlmids_file;
    
    vector<string> line_holder;
    string str;
    /*Loop thru all lines in input file:*/
    /*Open given file:*/
    
    while(1) {
        if(infile.is_open() == false) break;
    }
    
    infile.open(rlmids_file,std::fstream::in);
    while ( getline(infile, str) ) {
        
        split_str(str, line_holder, ",");
        
        RL_MID rlmid;
        rlmid.lmid_name = line_holder[0];
        rlmid.lmid_value = line_holder[1];
        rlmid.rmid_name = line_holder[2];
        rlmid.rmid_value = line_holder[3];
        rlmids.push_back( rlmid );
        //cout << rlmid.lmid_name << endl << rlmid.lmid_value << endl;
        line_holder.clear();
    }
    
    infile.close();
    //delete rlmids_file;
    
    cout << "Done!" << endl;
}

void Build_RLMIDS_Dictionary() 
{
    string mids[36] = {   "RL1,ACACGACGACT,RL1,AGTCGTGGTGT",
                        "RL2,ACACGTAGTAT,RL2,ATACTAGGTGT",
                        "RL3,ACACTACTCGT,RL3,ACGAGTGGTGT",
                        "RL4,ACGACACGTAT,RL4,ATACGTGGCGT",
                        "RL5,ACGAGTAGACT,RL5,AGTCTACGCGT",
                        "RL6,ACGCGTCTAGT,RL6,ACTAGAGGCGT",
                        "RL7,ACGTACACACT,RL7,AGTGTGTGCGT",
                        "RL8,ACGTACTGTGT,RL8,ACACAGTGCGT",
                        "RL9,ACGTAGATCGT,RL9,ACGATCTGCGT",
                        "RL10,ACTACGTCTCT,RL10,AGAGACGGAGT",
                        "RL11,ACTATACGAGT,RL11,ACTCGTAGAGT",
                        "RL12,ACTCGCGTCGT,RL12,ACGACGGGAGT",
                        "RL13,AGACTCGACGT,RL13,ACGTCGGGTCT",
                        "RL14,AGTACGAGAGT,RL14,ACTCTCGGACT",
                        "RL15,AGTACTACTAT,RL15,ATAGTAGGACT",
                        "RL16,AGTAGACGTCT,RL16,AGACGTCGACT",
                        "RL17,AGTCGTACACT,RL17,AGTGTAGGACT",
                        "RL18,AGTGTAGTAGT,RL18,ACTACTAGACT",
                        "RL19,ATAGTATACGT,RL19,ACGTATAGTAT",
                        "RL20,CAGTACGTACT,RL20,AGTACGTGCTG",
                        "RL21,CGACGACGCGT,RL21,ACGCGTGGTCG",
                        "RL22,CGACGAGTACT,RL22,AGTACTGGTCG",
                        "RL23,CGATACTACGT,RL23,ACGTAGTGTCG",
                        "RL24,CGTACGTCGAT,RL24,ATCGACGGACG",
                        "RL25,CTACTCGTAGT,RL25,ACTACGGGTAG",
                        "RL26,GTACAGTACGT,RL26,ACGTACGGTAC",
                        "RL27,GTCGTACGTAT,RL27,ATACGTAGGAC",
                        "RL28,GTGTACGACGT,RL28,ACGTCGTGCAC",
                        "RL29,ACACAGTGAGT,RL29,ACTCACGGTGT",
                        "RL30,ACACTCATACT,RL30,AGTATGGGTGT",
                        "RL31,ACAGACAGCGT,RL31,ACGCTGTGTGT",
                        "RL32,ACAGACTATAT,RL32,ATATAGTGTGT",
                        "RL33,ACAGAGACTCT,RL33,AGAGTCTGTGT",
                        "RL34,ACAGCTCGTGT,RL34,ACACGAGGTGT",
                        "RL35,ACAGTGTCGAT,RL35,ATCGACAGTGT",
                        "RL36,ACGAGCGCGCT,RL36,AGCGCGCGCGT"   
                };
    
    vector<string> line_holder;
    for ( int i = 0; i<36; ++i ) {
        
        split_str(mids[i], line_holder, ",");
        
        RL_MID rlmid;
        rlmid.lmid_name = line_holder[0];
        rlmid.lmid_value = line_holder[1];
        rlmid.rmid_name = line_holder[2];
        rlmid.rmid_value = line_holder[3];
        rlmids.push_back( rlmid );
        //cout << rlmid.lmid_name << endl << rlmid.lmid_value << endl;
        line_holder.clear();
    }
    
}

string GetTruSeqAdapter(string type, short index) {
    /*i5  adapter*/
    string tru_seq_adapter = "";
   
    string D500[8] = {"TATAGCCT",
                       "ATAGAGGC",
                       "CCTATCCT",
                       "GGCTCTGA",
                       "AGGCGAAG",
                       "TAATCTTA",
                       "CAGGACGT",
                       "GTACTGAC"};
    
    string D700[12] = {"ATTACTCG",
                       "TCCGGAGA",
                       "CGCTCATT",
                       "GAGATTCC",
                       "ATTCAGAA",
                       "GAATTCGT",
                       "CTGAAGCT",
                       "TAATGCGC",
                       "CGGCTATG",
                       "TCCGCGAA",
                       "TCTCGCGC",
                       "AGCGATAG"};
    
    if( (type == "i5") && (index < 8) ) {
        tru_seq_adapter = tmpl_i5_1 + D500[index] + tmpl_i5_2;
    }
    if( (type == "i7") && (index) < 12 ) {
        tru_seq_adapter = tmpl_i7_1 + D700[index] + tmpl_i7_2;
    }
    
    return tru_seq_adapter;
}
