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
       for(int i=0; i<(int)roche_names.size(); ++i)
       {
        if(debug_flag ) 
        {
          process_sff_to_fastq(roche_names[i], 0);
          process_fastq_to_sff((char*)"out.sff"); 
          process_sff_to_fastq((char*)"out.sff",  0);
          return;
        }
       }
        
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
