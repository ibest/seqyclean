

#include "TrimLucy.h"

using namespace std;
void Trim(char* filename) {
    std::ifstream lucyfile;
    
    std::string line;
    
   
    
    printf("Parsing .fastq file...\n");
    
    vector<string> record_block;
    int ii = 0;
    while(1) {
        if(lucyfile.is_open() == false) break;
    }
    lucyfile.open( filename );
    
    while ( getline(lucyfile,line) ) {
       
       if(ii==0) {
         record_block.push_back(line); /*Read ID*/
       }
       if(ii==1) {
         record_block.push_back(line); /*actual data*/
       }
       if(ii==2) {
         record_block.push_back(line);/*a symbol "+"*/
       }
       if(ii==3) {
           
           record_block.push_back(line); /*Quality scores*/
           
           Read *read = new Read();
           read->readID = record_block[0];
           read->contaminants = 0;
           read->discarded = 0;
           read->initial_length = record_block[1].length();
           read->lclip = 0;
           read->rclip = 0;
           read->read = record_block[1];
           read->quality = record_block[3];
           read->clip_found = false;
           read->rlmid.lmid_id = 0;
           read->rlmid.rmid_id = 0;
           //counter++;
           
           reads.push_back(read);
           
       }
        
       ii++;
       if(ii == 4) {
           ii = 0;
           record_block.clear();
       }
    }
    
    lucyfile.close();
    //sleep(5);
    
    /*Actual trimiing*/
    for(int i=0; i<(int)reads.size(); i++) {
        vector<string> tmp;        
        split_str(reads[i]->readID, tmp, " ");
        /*Right clipping*/
        reads[i]->read = reads[i]->read.substr( 0, atoi( tmp[tmp.size()-1].c_str() ) - 1 );
        reads[i]->quality = reads[i]->quality.substr( 0, atoi( tmp[tmp.size()-1].c_str() ) - 1 );
        /*Left clipping*/
        reads[i]->read = reads[i]->read.substr( atoi(tmp[tmp.size()-2].c_str()) - 1, reads[i]->read.length() - (atoi(tmp[tmp.size()-2].c_str()) - 1) );
        reads[i]->quality = reads[i]->quality.substr( atoi(tmp[tmp.size()-2].c_str()) - 1, reads[i]->quality.length() - (atoi(tmp[tmp.size()-2].c_str()) - 1) );
        
    }
    
    /*Write results*/
    FILE* output_file;
    
    vector<string> tmp_rep;
    split_str(string(filename), tmp_rep, "//");
    output_file =  fopen( (tmp_rep[tmp_rep.size() - 1] + "_trimmed.fastq").c_str(), "w" );
    cout << (tmp_rep[tmp_rep.size() - 1] + "_trimmed.fastq").c_str() << endl;
    
    for(int i=0; i<(int)reads.size(); i++) {
        if(reads[i]->discarded == 0) {
            if(reads[i]->lclip >= (int)reads[i]->read.length()  ) {
                reads[i]->discarded = 1;
//                discard_counter = discard_counter+1;
                continue;
            }
            fputs( reads[i]->readID.c_str(), output_file );
            fputc( '\n', output_file );
            fputs( reads[i]->read.c_str(), output_file );
            fputc( '\n', output_file );
            fputc( '+', output_file );
            fputc( '\n', output_file );
            fputs( reads[i]->quality.c_str(), output_file );
            fputc( '\n', output_file );
        } 
    }
    
    fclose(output_file);
    
}
