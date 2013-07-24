/* 
 * File:   iz_SSAHA.cpp
 * Author: ilya
 * 
 * Created on 1 Август 2012 г., 14:50
 */

#include "iz_SSAHA.h"

bool CompFunctIndex(Master_struct a, Master_struct b) { 
    return (a.index < b.index); 
}
bool CompFunctShift(Master_struct a, Master_struct b) { 
    return (a.shift < b.shift) ; 
}
bool CompFunctOffset(Master_struct a, Master_struct b) { 
    
    return (a.offset < b.offset) ; 
}

bool CompFunctLargestArea(Hit_struct a, Hit_struct b) { 
    
    return (a.cnt > b.cnt)  ; 
}

iz_SSAHA::iz_SSAHA() {
}

iz_SSAHA::iz_SSAHA(const iz_SSAHA& orig) {
}

iz_SSAHA::~iz_SSAHA() {
}

void iz_SSAHA::AddElementToMasterStruct( Master_struct a) {
    this->m_structs.push_back(a);
}

void iz_SSAHA::RunSort() {
    // using function as comp
    
  std::sort ( this->m_structs.begin(), this->m_structs.end() ,CompFunctIndex);
  
  int index = 0;
  int i=1;
  for(i=1; i<(int)m_structs.size(); i++) {
      if(m_structs[i-1].index != m_structs[i].index) {
        std::sort ( this->m_structs.begin()+index, this->m_structs.end() - (m_structs.size() - i) ,CompFunctShift);
        index = i;
      }
      
  }
  
  std::sort ( this->m_structs.begin()+index, this->m_structs.end()  ,CompFunctShift);
  
  index = 0;
  for(i=1; i<(int)m_structs.size(); i++) {
      if((m_structs[i-1].index == m_structs[i].index) && (m_structs[i-1].shift != m_structs[i].shift) ) {
        std::sort ( this->m_structs.begin()+index, this->m_structs.end() - (m_structs.size() - i) ,CompFunctOffset);
        index = i;
      }
  }
  
  std::sort ( this->m_structs.begin()+index, this->m_structs.end() ,CompFunctOffset);
  
  
}

void iz_SSAHA::MakeTable1(string read) {
    /*KMER_SIZE = 2*/
    int seqid = 0;
    for(int i=0; i<(int)read.length()-1; i++) {
    
        Hit_struct h_struct;
        h_struct.seq_id = seqid;
        if(i == 60) seqid++;
        h_struct.pos = i;
        
        vector<Hit_struct> tmp_array;
        
        string k_str = read.substr(i,2);
        
        this->it_Table1 = this->Table1.find(k_str);
    
        if(this->it_Table1 == this->Table1.end()) {
              tmp_array.push_back(h_struct);
              this->Table1.insert(std::pair<string, vector<Hit_struct> >(k_str, tmp_array));
              tmp_array.clear();
        } else {
              (*(this->it_Table1)).second.push_back(h_struct);
        }
        
        k_str.clear();
        tmp_array.clear();
    }
    
    read.clear();
}


AlignResult iz_SSAHA::Find(string &read, string &query_str) {
    
    /*To uppercase:*/
    stoupper(read);
    stoupper(query_str);
        
    MakeTable1(read);
    
    //Building a H-array
    int t = 0; 
    int i;
    for(i=0; i<(int)query_str.length()-1; i++) {
        string k_str = query_str.substr(i,2);
        
        it_Table1 = Table1.find(k_str);
    
        if(this->it_Table1 != this->Table1.end()) {
            for(int j=0; j<(int)(*it_Table1).second.size(); j++) {
              Master_struct m_struct;
              m_struct.index = (*it_Table1).second[j].seq_id;
              m_struct.shift = (*it_Table1).second[j].pos - t;
              m_struct.offset = (*it_Table1).second[j].pos;
              
              m_structs.push_back(m_struct);
            } 
            
            
            
        }
        
        k_str.clear();
        
        t++;
        
    }
    
    RunSort();
    
    vector<Hit_struct> tmp_array;
    //Now we have to find the longest match:
    int cnt=0;
    
    AlignResult al_res;
    
    if(m_structs.size() == 0) {
        al_res.found_flag = false;
        //cout << al_res.found_flag << "ddddd" << endl;
        return al_res; /*No mathes!*/
    } else {
        al_res.found_flag = true;
    }
    
    for(i=0; i<(int)m_structs.size()-1; i++) {
        cnt++;
        if(m_structs[i].shift != m_structs[i+1].shift) {
            Hit_struct hs;
            hs.cnt = cnt;
            hs.pos = m_structs[i].shift;
            tmp_array.push_back(hs);
            cnt = 0;
        } 
    }
    Hit_struct hs;
    hs.cnt = ++cnt;
    hs.pos = m_structs[i].shift;
    tmp_array.push_back(hs);
    
    //for(i=0; i<tmp_array.size(); i++) {
    //    cout << tmp_array[i].pos << " " << tmp_array[i].cnt << endl;
    //}
    
    //Largest
    std::sort ( tmp_array.begin(), tmp_array.end() ,CompFunctLargestArea);
    
    tmp_array[0].pos = tmp_array[0].pos < 0 ? 0 : tmp_array[0].pos;
    
    string *s1, *s2, seq1, seq2;
    
    if (tmp_array[0].pos + query_str.length() > read.length() ) {
        //s1 = new string( read.substr(tmp_array[0].pos, read.length() - tmp_array[0].pos) );
        al_res.found_flag = false;
        return al_res;
        
    } else {
        s1 = new string(read.substr(tmp_array[0].pos, query_str.length()));
        
    }
    
    s2 = new string(tmp_array[0].pos + query_str.length() > read.length() ? query_str.substr(0, read.length() - tmp_array[0].pos ) : query_str );
    seq2 = string(tmp_array[0].pos + query_str.length() > read.length() ? query_str.substr(0, read.length() - tmp_array[0].pos ) : query_str) ;
    //cout << al_res.seq_1 << endl << al_res.seq_2 << endl;
    al_res = banded( *s1, *s2 ,2,1 );
    al_res.pos = tmp_array[0].pos ;
    al_res.found_flag = true;
    al_res.seq_1 = seq1;
    al_res.seq_2 = seq2;
    //cout << read.substr(tmp_array[0].pos, query_str.length()) << endl;
    //cout << query_str << endl;
    
    //cout << al_res.seq_1_al << endl;
    //cout << al_res.seq_2_al << endl;
    s1->clear();
    s2->clear();
    delete s1;
    delete s2;
    
    tmp_array.clear();
    m_structs.clear();
    this->Table1.clear();
    //read.clear();
    //query_str.erase(query_str.begin(),query_str.end());
    
    return al_res;
}