#include "Illumina.h"

/*i5 adapter*/
string tmpl_i5_1 = "AATGATACGGCGACCACCGAGATCTACAC";
string tmpl_i5_2 = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT";
 /*i7 adapter*/
string tmpl_i7_1 = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
string tmpl_i7_2 = "ATCTCGTATGCCGTCTTCTGCTTG";

//Dynamic Illumina: does not need space to store reads:
void IlluminaDynamic()
{
    unsigned long long pe1_bases_anal, pe2_bases_anal, pe_bases_kept, pe_bases_discarded, se_pe1_bases_kept, se_pe2_bases_kept;
    unsigned long pe_discard_cnt;
    double avg_trim_len_pe1, avg_trim_len_pe2;
    
    /*Raw implementation of average. Later I will come with a better algorithm*/
    unsigned long long avg_bases_pe1 = 0;
    unsigned long long avg_bases_pe2 = 0;
    unsigned long long avg_left_clip_1 = 0;
    unsigned long long avg_left_clip_2 = 0;
    unsigned long long avg_right_clip_1 = 0;
    unsigned long long avg_right_clip_2 = 0;
    
    pe_bases_kept = pe_bases_discarded = se_pe1_bases_kept = se_pe2_bases_kept = 0;
    pe_discard_cnt = 0;
    pe1_bases_anal = pe2_bases_anal = 0;        
    avg_trim_len_pe1 = avg_trim_len_pe2 = 0;
    
    unsigned long cnt1_avg, cnt2_avg; cnt1_avg = cnt2_avg = 0; //Counters needed for calculating the average trimming length
    unsigned long cnt_avg_len1, cnt_avg_len2; cnt_avg_len1 = cnt_avg_len2 = 0;
                 
    double avg_len_pe1, avg_len_pe2; avg_len_pe1 = avg_len_pe2 = 0.0;
    double cnt_right_trim_pe1, avg_right_trim_len_pe1, cnt_right_trim_pe2, avg_right_trim_len_pe2; 
    double cnt_left_trim_pe1, avg_left_trim_len_pe1, cnt_left_trim_pe2, avg_left_trim_len_pe2;
    
    cnt_right_trim_pe1 = avg_right_trim_len_pe1 = cnt_right_trim_pe2 = avg_right_trim_len_pe2 = 0;
    cnt_left_trim_pe1 = avg_left_trim_len_pe1 = cnt_left_trim_pe2 = avg_left_trim_len_pe2 = 0;
    
    unsigned long cnt1, cnt2; cnt1 = cnt2 = 0;
    unsigned long pe_accept_cnt, se_pe1_accept_cnt, se_pe2_accept_cnt; pe_accept_cnt = se_pe1_accept_cnt = se_pe2_accept_cnt = 0;
    unsigned long ts_adapters1, ts_adapters2; ts_adapters1 = ts_adapters2 = 0;
    unsigned long num_vectors1, num_vectors2; num_vectors1 = num_vectors2 = 0;
    unsigned long num_contaminants1, num_contaminants2; num_contaminants1 = num_contaminants2 = 0;
    unsigned long accepted1,accepted2; accepted1 = accepted2 = 0;
    unsigned long discarded1,discarded2; discarded1 = discarded2 = 0;
    unsigned long perfect_ov_cnt, partial_ov_cnt; perfect_ov_cnt = partial_ov_cnt = 0;
//    unsigned long discarded_by_quality1, discarded_by_quality2; discarded_by_quality1 = discarded_by_quality2 = 0;
    unsigned long discarded_by_contaminant1, discarded_by_contaminant2; discarded_by_contaminant1 = discarded_by_contaminant2 = 0;
    unsigned long discarded_by_read_length1, discarded_by_read_length2; discarded_by_read_length1 = discarded_by_read_length2 = 0;
//    unsigned long discarded_by_vector1 , discarded_by_vector2; discarded_by_vector1 = discarded_by_vector2 = 0;
    /*Left trims*/
    unsigned long left_trimmed_by_quality1 , left_trimmed_by_quality2; left_trimmed_by_quality1 = left_trimmed_by_quality2 = 0;
    unsigned long left_trimmed_by_vector1 , left_trimmed_by_vector2; left_trimmed_by_vector1 = left_trimmed_by_vector2 = 0;
    /*Right trims/discards*/
    unsigned long right_trimmed_by_quality1 , right_trimmed_by_quality2; right_trimmed_by_quality1 = right_trimmed_by_quality2 = 0;
    unsigned long right_trimmed_by_adapter1 , right_trimmed_by_adapter2; right_trimmed_by_adapter1 = right_trimmed_by_adapter2 = 0;
    unsigned long right_trimmed_by_vector1 , right_trimmed_by_vector2;  right_trimmed_by_vector1 = right_trimmed_by_vector2 = 0;
    
    unsigned long right_trimmed_by_polyat1, right_trimmed_by_polyat2; right_trimmed_by_polyat1 = right_trimmed_by_polyat2 = 0;
    unsigned long left_trimmed_by_polyat1, left_trimmed_by_polyat2; left_trimmed_by_polyat1 = left_trimmed_by_polyat2 = 0;
    
            
    unsigned long duplicates = 0;
    
    fstream rep_file1, rep_file2, pe_output_file1, pe_output_file2, shuffle_file, se_file, overlap_file;
    rep_file1.open(rep_file_name1.c_str(),ios::out);
    rep_file2.open(rep_file_name2.c_str(),ios::out);
    
    if(overlap_flag)
        overlap_file.open(overlap_file_name.c_str(),ios::out);
    
    rep_file1 << "ReadID\tlclip\trclip\tTruSeq_pos\tRaw_read_length\tLlucy\tRlucy\tDiscarded\tContaminants\tVecStart\tVecEnd\tVecLen\n";
    rep_file2 << "ReadID\tlclip\trclip\tTruSeq_pos\tRaw_read_length\tLlucy\tRlucy\tDiscarded\tContaminants\tVecStart\tVecEnd\tVecLen\n";
    
    cout << "Running the Illumina cleaning process..." << endl;
    sum_stat << "Running the Illumina cleaning process..." << endl;
    
    
    
    vector<string> record_block1, record_block2;
    
    if (!shuffle_flag)
    {
        pe_output_file1.open( pe_output_filename1.c_str(), ios::out );
        pe_output_file2.open( pe_output_filename2.c_str(), ios::out );
    } 
    else
    {
        shuffle_file.open( shuffle_filename.c_str(), ios::out );
    }    
    se_file.open( se_filename.c_str(), ios::out );
    
    string st_str; //output statistics
    
    for(int jj=0; jj<(int)pe1_names.size(); ++jj)
    {
    
        bool adapter_found1 = false;
        bool adapter_found2 = false;
        
        string query_string1 = "NA";
        string query_string2 = "NA";
        
        int ii = 0;
        
        std::string line1, line2;
        igzstream in1(/*fastq_file1*/pe1_names[jj]); //for R1
        igzstream in2(/*fastq_file2*/pe2_names[jj]); //for R2
        
        cout << "Processing files: " << pe1_names[jj] << ", " << pe2_names[jj] << "\n";
        sum_stat << "Processing files: " << pe1_names[jj] << ", " << pe2_names[jj] << "\n";
        
        while ( getline(in1,line1) && getline(in2,line2) )
        {
                /*Read ID*/
                if(ii==0) 
                {
                        
                        
                        //Check for order
                        vector <string> fields1, fields2;
                        split_str( line1, fields1, " " );
                        split_str( line2, fields2, " " );
                        //cout << line1 << endl;
                        if( (fields1[0] != fields2[0] ) && !old_style_illumina_flag)
                        {
                            cout << "Warning: read IDs do not match in input files: PE1-> " << pe1_names[jj] << ", PE2-> " << pe2_names[jj] << endl;
                            sum_stat << "Warning: read IDs do not match in input files: PE1-> " << pe1_names[jj] << ", PE2-> " << pe2_names[jj] << endl;
                        }
                        
                        fields1.clear();
                        fields2.clear();
                        
                        if ( new2old_illumina && !old_style_illumina_flag ) //if convert to old-style illumina headers is true and not old illumina files.
                        {
                            split_str( line1, fields1, " " );
                            split_str( fields1[0], fields2, ":" );
                            line1 = string(fields2[0] + "_" + fields2[2] + ":" + fields2[3] + ":" + fields2[4] + ":" + fields2[5] + ":" + fields2[6] + "#0/1"); //+ fields1[1].substr(0,1) );//+ " (" + line1 + ")");
                            
                            fields1.clear();
                            fields2.clear();
                            
                            split_str( line2, fields1, " " );
                            split_str( fields1[0], fields2, ":" );
                            
                            line2 = string(fields2[0] + "_" + fields2[2] + ":" + fields2[3] + ":" + fields2[4] + ":" + fields2[5] + ":" + fields2[6] + "#0/2");// + fields1[1].substr(0,1) ); // + " (" + line2 + ")");
                            
                            fields1.clear();
                            fields2.clear();
                        }
                        
                        record_block1.push_back(line1); 
                        record_block2.push_back(line2);
                        
                        ii++;
                        continue;
                }
                /*DNA string*/
                if(ii==1) 
                {
                        record_block1.push_back(line1); /*DNA string*/
                        record_block2.push_back(line2);
                        ii++;
                        continue;
                }
                /*a symbol "+"*/
                if(ii==2) 
                {
                        record_block1.push_back(line1);
                        record_block2.push_back(line2);
                        ii++;
                        continue;
                }
                if(ii==3) 
                {
                        
                    cnt1+=1; cnt2+=1;
                    ii=0;
           
                    Read *read1 = new Read();
                    read1->illumina_readID = record_block1[0];
                    read1->initial_length = record_block1[1].length();
                    read1->read = record_block1[1];
                    read1->illumina_quality_string = line1;
                    pe1_bases_anal += (unsigned long long)read1->read.length();
                        
                        if(read1->initial_length <= minimum_read_length)
                        {
                            cout << "Warming: in PE1 file raw read length is less or equal then minimum_read_length\n" ; 
                            sum_stat << "Warming: in PE1 file raw read length is less or equal then minimum_read_length\n" ; 
                        }
                        
                        
          
                        Read *read2 = new Read();
                        read2->illumina_readID = record_block2[0];
                        read2->initial_length = record_block2[1].length();
                        read2->read = record_block2[1];
                        read2->illumina_quality_string = line2;
                        pe2_bases_anal += (unsigned long long)read2->read.length();
                        
                        if(read2->initial_length <= minimum_read_length)
                        {
                            cout << "Warming: in PE2 file, the raw read length is less or equal than minimum_read_length\n" ; 
                            sum_stat << "Warming: in PE2 file, the raw read length is less or equal than minimum_read_length\n" ; 
                        }
                        
                        //Serial realization - useful for debugging if something does not work as expected
                        //Duplicates removal
                        if(rem_dup)
                            screen_duplicates(read1, read2, duplicates);
                        
                        
                        IlluminaDynRoutine(read1, adapter_found1, query_string1);
                        
                        if(read1->discarded_by_contaminant == 0) {
                            IlluminaDynRoutine(read2, adapter_found2, query_string2);
                        } else {
                            read2->discarded_by_contaminant = 1;
                            read2->discarded = 1;
                        }
                        
                        if(read2->discarded_by_contaminant == 1) {
                            read1->discarded_by_contaminant = 1;
                            read1->discarded = 1;
                        }
                        
                        //Establishing the most conservative adapters
                        if( (read1->tru_sec_found == 1) && (read2->tru_sec_found == 1))
                        {
                            //Take the most conservative position:
                            if( read1->tru_sec_pos < read2->tru_sec_pos )
                            {
                                read2->tru_sec_pos = read1->tru_sec_pos;
                                read2->tru_sec_found = 1;
                            } 
                            else
                            {
                                read1->tru_sec_pos = read2->tru_sec_pos;
                                read1->tru_sec_found = 1;
                            }
                        }
                        else if( (read1->tru_sec_found == 1) && (read2->tru_sec_found == 0) )
                        {
                            read2->tru_sec_pos = read1->tru_sec_pos;
                            read2->tru_sec_found = 1;
                        }
                        else if( (read1->tru_sec_found == 0) && (read2->tru_sec_found == 1) )
                        {
                            read1->tru_sec_pos = read2->tru_sec_pos;
                            read1->tru_sec_found = 1;
                        }
                        else if( (read1->tru_sec_found == 0) && (read2->tru_sec_found == 0) && trim_adapters_flag)
                        {
                            int o = find_overlap_pos(read1->read, MakeRevComplement(read2->read), adapterlength, false);
                            if( (o < 0) && (o != -10000)) {
                               //dovetails, remove adapters:
                               read1->tru_sec_found = 1; read2->tru_sec_found = 1;
                               read1->tru_sec_pos = read1->initial_length + o - 1; read2->tru_sec_pos = read2->initial_length + o - 1;
                            }
                        }
                        
          
                        //Checking for overlap:
                        bool overlap_found = false;
                        Read *c = new Read();
                        if( overlap_flag && (read1->discarded == 0) && (read2->discarded == 0) ) {
                              Read *s1 = new Read(); s1->read = read1->read; s1->illumina_quality_string = read1->illumina_quality_string;
                              Read *s2 = new Read(); s2->read = read2->read; s2->illumina_quality_string = read2->illumina_quality_string;
                                        
                              int ov;
                              string tread = MakeRevComplement(s2->read);
                              ov = find_overlap_pos(s1->read, tread, adapterlength, true);
                              
                              if(ov == 0) {
                                  perfect_ov_cnt += 1;
                                  //Perfect overlap, make consensus sequence:
                                  s2->read = MakeRevComplement(s2->read);
                                  c = make_consensus(s1, s2);  //don't cut, just overlap
                                  overlap_found = true;
                              }
                              if( ov > 0 ) {
                                  partial_ov_cnt += 1;
                                  //Overlap, make consensus sequence:
                                  s1->read = s1->read.substr(ov, s1->read.length() - ov);  // cut off left end
                                  s2->read = MakeRevComplement(s2->read.substr(0,s2->read.length() - ov));  //cut off right end
                                    
                                  c = make_consensus(s1, s2);
                                    
                                  string c_se = read1->read.substr(0,ov) + c->read + MakeRevComplement(read2->read.substr(0,read2->read.length() - ov));
                                                
                                  string c_qual = read2->illumina_quality_string.substr(0,read2->illumina_quality_string.length() - ov);
                                  c_qual = read1->illumina_quality_string.substr(0,ov) + c->illumina_quality_string + string ( c_qual.rbegin(), c_qual.rend() );
                                    
                                  //cout << read1->read.substr(0,ov) << endl << read1->illumina_quality_string.substr(0,ov) << endl;
                                    
                                  c->read = c_se;
                                  c->illumina_quality_string = c_qual;
                                  overlap_found = true;
                                    //cout << c->read << endl << c->illumina_quality_string << endl;
                              }
                        }
                        
                        if(overlap_found) {
                            //Establis new clip points for SE overlap:
                            MakeClipPointsIllumina(read1);
                            MakeClipPointsIllumina(read2);
                            
                            c->lclip = read1->lclip;
                            c->rclip = c->read.length() - (read2->initial_length - read2->rclip);
                            
                            c->read = c->read.substr(0 , c->rclip );
                            c->illumina_quality_string = c->illumina_quality_string.substr(0,c->rclip) ; 
                            c->read = c->read.substr( c->lclip, c->rclip - c->lclip );
                            c->illumina_quality_string = c->illumina_quality_string.substr( c->lclip, c->rclip - c->lclip );
                            
                            vector<string> temp_id;
                            split_str( read1->illumina_readID, temp_id, " " );
                            vector<string> temp_id1;
                            split_str( temp_id[0], temp_id1, " " );
                            c->illumina_readID = temp_id1[0];
                            WriteSEOverlap(overlap_file, c);
                            temp_id.clear();
                            temp_id1.clear();
                        } else {
                                        
                                if(read1->discarded == 0)
                                {
                                        MakeClipPointsIllumina(read1);
                                }
                                
                                if(read2->discarded == 0)
                                {
                                        MakeClipPointsIllumina(read2);
                                }
                                
                                if( read1->discarded_by_contaminant == 0)
                                {    
                                        if( read1->lclip >= read1->rclip ) { read1->discarded = 1; read1->discarded_by_read_length = 1; } 
                                        if( read1->lclip >= (int)read1->read.length() ) { read1->discarded = 1; read1->discarded_by_read_length = 1; }
                                        if( read1->rclip > (int)read1->read.length() ) { read1->rclip = read1->read.length(); }
                                        if( (int)read1->read.length() < minimum_read_length ) { read1->discarded = 1; read1->discarded_by_read_length = 1; }
                                        if( (read1->rclip - read1->lclip) < minimum_read_length ) { read1->discarded = 1; read1->discarded_by_read_length = 1; }
                                }
                            
                                if( read2->discarded_by_contaminant == 0)
                                {
                                        if( read2->lclip >= read2->rclip ) {read2->discarded = 1; read2->discarded_by_read_length = 1; }
                                        if( read2->lclip >= (int)read2->read.length() ) {read2->discarded = 1; read2->discarded_by_read_length = 1; }
                                        if( read2->rclip > (int)read2->read.length() ) {read2->rclip = read2->read.length(); }
                                        if( (int)read2->read.length() < minimum_read_length ) {read2->discarded = 1; read2->discarded_by_read_length = 1; }
                                        if( (read2->rclip - read2->lclip) < minimum_read_length ) { read2->discarded = 1; read2->discarded_by_read_length = 1; }
                                }

                                if( (read1->discarded == 0) && (read2->discarded == 0) )
                                {
                                        
                                    if(  read1->rclip < read1->initial_length  )
                                    {
                                        cnt_right_trim_pe1 += 1;
                                        //////////avg_right_trim_len_pe1 = GetAvg( avg_right_trim_len_pe1, cnt_right_trim_pe1, read1->initial_length - read1->rclip );
                                        avg_right_clip_1 += read1->initial_length - read1->rclip;
                                        avg_right_trim_len_pe1 = avg_right_clip_1/cnt_right_trim_pe1;
                                    }
                                    if(read1->lclip > 0)
                                    {
                                        cnt_left_trim_pe1 += 1;
                                        /////////////avg_left_trim_len_pe1 = GetAvg( avg_left_trim_len_pe1, cnt_left_trim_pe1, read1->lclip );
                                        avg_left_clip_1 += read1->lclip;
                                        avg_left_trim_len_pe1 = avg_left_clip_1/cnt_left_trim_pe1;
                                    }
                                    
                                    read1->read = read1->read.substr(0 , read1->rclip );
                                    read1->illumina_quality_string = read1->illumina_quality_string.substr(0,read1->rclip) ; 
                                    read1->read = read1->read.substr( read1->lclip, read1->rclip - read1->lclip );
                                    read1->illumina_quality_string = read1->illumina_quality_string.substr( read1->lclip, read1->rclip - read1->lclip );
                 
                                    if(  read2->rclip < read2->initial_length  )
                                    {
                                        cnt_right_trim_pe2 += 1;
                                        //////////avg_right_trim_len_pe2 = GetAvg( avg_right_trim_len_pe2, cnt_right_trim_pe2, read2->initial_length - read2->rclip );
                                        avg_right_clip_2 += read2->initial_length - read2->rclip;
                                        avg_right_trim_len_pe2 = avg_right_clip_2/cnt_right_trim_pe2;
                                    }
                                    if(read2->lclip > 0)
                                    {
                                        cnt_left_trim_pe2 += 1;
                                        /////////////avg_left_trim_len_pe2 = GetAvg( avg_left_trim_len_pe2, cnt_left_trim_pe2, read2->lclip );
                                        avg_left_clip_2 += read2->lclip;
                                        avg_left_trim_len_pe2 = avg_left_clip_2/cnt_left_trim_pe2;
                                    }
                                        
                                    read2->read = read2->read.substr(0 , read2->rclip );
                                    read2->illumina_quality_string = read2->illumina_quality_string.substr(0,read2->rclip) ; 
                                    read2->read = read2->read.substr( read2->lclip, read2->read.length() - read2->lclip );
                                    read2->illumina_quality_string = read2->illumina_quality_string.substr( read2->lclip, read2->illumina_quality_string.length() - read2->lclip );
        	
                                    if (!shuffle_flag)
                                        {
                                                WritePEFile(pe_output_file1, read1);
                                                WritePEFile(pe_output_file2, read2);
                                                pe_accept_cnt+=1;
                                                pe_bases_kept += read1->read.length();
                                                pe_bases_kept += read2->read.length();
                                        } else 
                                        {
                                                WriteShuffleFile( shuffle_file, read1, read2 );
                                                pe_accept_cnt+=1;
                                                pe_bases_kept += read1->read.length();
                                                pe_bases_kept += read2->read.length();
                                        }
                                    
                 
                                    if( read1->initial_length > (int)read1->read.length() )
                                    {
                                        cnt1_avg+=1;
                                        ///////////avg_trim_len_pe1 = GetAvg( avg_trim_len_pe1, cnt1_avg,  read1->rclip - read1->lclip );//read1->initial_length - read1->read.length()
                                        avg_bases_pe1 += read1->rclip - read1->lclip;
                                        avg_trim_len_pe1 = avg_bases_pe1/cnt1_avg;
                                    }
                 
                                    if( read2->initial_length > (int)read2->read.length() )
                                    {
                                        cnt2_avg+=1;
                                        /////////////avg_trim_len_pe2 = GetAvg( avg_trim_len_pe2, cnt2_avg, read2->rclip - read2->lclip );//read2->initial_length - read2->read.length()*
                                        avg_bases_pe2 += read2->rclip - read2->lclip;
                                        avg_trim_len_pe2 = avg_bases_pe2/cnt2_avg;
                                    }
                 
                                    cnt_avg_len1 +=1; cnt_avg_len2 +=1;
                 
                                    ///////////avg_len_pe1 = GetAvg( avg_len_pe1, cnt_avg_len1, read1->read.length() );
                                    /////////////avg_len_pe2 = GetAvg( avg_len_pe2, cnt_avg_len2, read2->read.length() );
                
                                        
                                } else if ((read1->discarded == 0) && (read2->discarded == 1)) 
                                {
                                    
                                    read1->read = read1->read.substr(0 , read1->rclip );
                                    read1->illumina_quality_string = read1->illumina_quality_string.substr(0,read1->rclip) ; 
                                    read1->read = read1->read.substr( read1->lclip, read1->rclip - read1->lclip );
                                    read1->illumina_quality_string = read1->illumina_quality_string.substr( read1->lclip, read1->rclip - read1->lclip );
                 
                                    
                                    if( new2old_illumina && !old_style_illumina_flag ) //if convert to old-style illumina headers is true and not old illumina files.
                                                read1->illumina_readID = read1->illumina_readID.substr(0,read1->illumina_readID.length()-2);                                  
                                    
                                    WriteSEFile( se_file, read1 );
                                    se_pe1_accept_cnt+=1;
                                    se_pe1_bases_kept += read1->read.length();
                 
                                    
                                } else if( (read1->discarded == 1) && (read2->discarded == 0) )
                                {
                                      
                                    read2->read = read2->read.substr(0 , read2->rclip );
                                    read2->illumina_quality_string = read2->illumina_quality_string.substr(0,read2->rclip) ; 
                                    read2->read = read2->read.substr( read2->lclip, read2->read.length() - read2->lclip );
                                    read2->illumina_quality_string = read2->illumina_quality_string.substr( read2->lclip, read2->illumina_quality_string.length() - read2->lclip );
        	 
                                    if( new2old_illumina && !old_style_illumina_flag ) //if convert to old-style illumina headers is true and not old illumina files.
                                                read2->illumina_readID = read2->illumina_readID.substr(0,read2->illumina_readID.length()-3);                                  
                                    
                                    WriteSEFile( se_file, read2 );
                                    se_pe2_accept_cnt +=1;
                                    se_pe2_bases_kept += read2->read.length();
                                } else 
                                {
                                        pe_discard_cnt+=1;
                                        pe_bases_discarded += read1->read.length();
                                        pe_bases_discarded += read2->read.length();
                                }
                        }
                        
                        rep_file1 << read1->illumina_readID.substr(1,read1->illumina_readID.length()-1) << "\t" << read1->lclip << "\t" << read1->rclip << "\t" << (read1->tru_sec_pos == -1 ? "NA" : int2str(read1->tru_sec_pos))  << "\t" << read1->initial_length << "\t" << (read1->lucy_lclip <= 1 ? 1 : read1->lucy_lclip) << "\t" << (read1->lucy_rclip <= 1 ? 1 : read1->lucy_rclip) << "\t" << read1->discarded << "\t" << read1->contaminants << "\t" << (vector_flag == true ? int2str(read1->v_start) : "NA") << "\t" << (vector_flag == true ? int2str(read1->v_end) : "NA") << "\t" << (vector_flag == true ? int2str(read1->vec_len) : "NA") << "\n";
                        rep_file2 << read2->illumina_readID.substr(1,read2->illumina_readID.length()-1) << "\t" << read2->lclip << "\t" << read2->rclip << "\t" << (read2->tru_sec_pos == -1 ? "NA" : int2str(read2->tru_sec_pos)) << "\t"  << read2->initial_length << "\t" << (read2->lucy_lclip <= 1 ? 1 : read2->lucy_lclip) << "\t" << (read2->lucy_rclip <= 1 ? 1 : read2->lucy_rclip) << "\t" << read2->discarded << "\t" << read2->contaminants << "\t" << (vector_flag == true ? int2str(read2->v_start) : "NA") << "\t" << (vector_flag == true ? int2str(read2->v_end) : "NA") << "\t" << (vector_flag == true ? int2str(read2->vec_len) : "NA") << "\n";
          
          
                        if (read1->tru_sec_found == 1) ts_adapters1++;
                        if (read1->vector_found == 1) num_vectors1++;
                        if (read1->contam_found == 1) num_contaminants1++;
                        if (read1->discarded == 0) accepted1++;
                        if (read1->discarded == 1) discarded1++;
                        if (read1->discarded_by_contaminant == 1) discarded_by_contaminant1++;
                        if (read1->discarded_by_read_length == 1) discarded_by_read_length1++;
                        if (read1->left_trimmed_by_quality == 1) left_trimmed_by_quality1++;
                        if (read1->left_trimmed_by_vector == 1) left_trimmed_by_vector1++;
                        if (read1->right_trimmed_by_quality == 1) right_trimmed_by_quality1++;
                        if (read1->right_trimmed_by_adapter == 1) right_trimmed_by_adapter1++;
                        if (read1->right_trimmed_by_vector == 1) right_trimmed_by_vector1++;
                        if (read1->right_trimmed_by_polyat == 1) right_trimmed_by_polyat1++;
                        if (read1->left_trimmed_by_polyat == 1) left_trimmed_by_polyat1++;
          
                        if (read2->tru_sec_found == 1) ts_adapters2++;
                        if (read2->vector_found == 1) num_vectors2++;
                        if (read2->contam_found == 1) num_contaminants2++;
                        if (read2->discarded == 0) accepted2++;
                        if (read2->discarded == 1) discarded2++;
                        if (read2->discarded_by_contaminant == 1) discarded_by_contaminant2++;
                        if (read2->discarded_by_read_length == 1) discarded_by_read_length2++;
                        if (read2->left_trimmed_by_quality == 1) left_trimmed_by_quality2++;
                        if (read2->left_trimmed_by_vector == 1) left_trimmed_by_vector2++;
                        if (read2->right_trimmed_by_quality == 1) right_trimmed_by_quality2++;
                        if (read2->right_trimmed_by_adapter == 1) right_trimmed_by_adapter2++;
                        if (read2->right_trimmed_by_vector == 1) right_trimmed_by_vector2++;
                        if (read2->right_trimmed_by_polyat == 1) right_trimmed_by_polyat2++;
                        if (read2->left_trimmed_by_polyat == 1) left_trimmed_by_polyat2++;
                        
                        record_block1.clear();
                        read1->illumina_readID.clear(); 
                        read1->illumina_quality_string.clear();
                        read1->read.clear();
          
                        record_block2.clear();
                        read2->illumina_readID.clear(); 
                        read2->illumina_quality_string.clear();
                        read2->read.clear();
          
                        delete read1;
                        delete read2;
                        
                        if( (cnt1 % 1000 ) == 0)
                        {
                            st_str = PrintIlluminaStatistics(cnt1, cnt2, 
                                    pe1_bases_anal, pe2_bases_anal, 
                                    ts_adapters1, ts_adapters2, 
                                    num_vectors1, num_vectors2, 
                                    num_contaminants1, num_contaminants2, 
                                    left_trimmed_by_quality1, left_trimmed_by_quality2,
                                    left_trimmed_by_vector1, left_trimmed_by_vector2, 
                                    avg_left_trim_len_pe1, avg_left_trim_len_pe2, 
                                    right_trimmed_by_adapter1, right_trimmed_by_adapter2, 
                                    right_trimmed_by_quality1,right_trimmed_by_quality2,
                                    right_trimmed_by_vector1,right_trimmed_by_vector2,
                                    avg_right_trim_len_pe1,avg_right_trim_len_pe2,
                                    discarded1, discarded2,
                                    discarded_by_contaminant1, discarded_by_contaminant2,
                                    discarded_by_read_length1, discarded_by_read_length2,
                                    pe_accept_cnt, pe_bases_kept, 
                                    pe_discard_cnt,pe_bases_discarded, 
                                    se_pe1_accept_cnt, se_pe1_bases_kept,
                                    se_pe2_accept_cnt, se_pe2_bases_kept,
                                    avg_trim_len_pe1, avg_trim_len_pe2,
                                    avg_len_pe1, avg_len_pe2,
                                    perfect_ov_cnt, partial_ov_cnt,
                                    duplicates,
                                    left_trimmed_by_polyat1, right_trimmed_by_polyat1,
                                    right_trimmed_by_polyat2, right_trimmed_by_polyat2
                                   );
                            
                            if (cnt1 > 1000)
                            {
                                vector<string> t;
                                split_str(st_str, t, "\n");
                                for(int kk=0; kk<(int)t.size(); ++kk)
                                {
                                   cout << "\033[A\033[2K";
                                }
                                t.clear();
                            }
                            
                            cout << st_str;
                            
                        }
          
                }
        }
        in1.close();
        in2.close();
    }
    
    
    st_str = PrintIlluminaStatistics(cnt1, cnt2, 
                            pe1_bases_anal, pe2_bases_anal, 
                            ts_adapters1, ts_adapters2, 
                            num_vectors1, num_vectors2, 
                            num_contaminants1, num_contaminants2, 
                            left_trimmed_by_quality1, left_trimmed_by_quality2,
                            left_trimmed_by_vector1, left_trimmed_by_vector2, 
                            avg_left_trim_len_pe1, avg_left_trim_len_pe2, 
                            right_trimmed_by_adapter1, right_trimmed_by_adapter2, 
                            right_trimmed_by_quality1,right_trimmed_by_quality2,
                            right_trimmed_by_vector1,right_trimmed_by_vector2,
                            avg_right_trim_len_pe1,avg_right_trim_len_pe2,
                            discarded1, discarded2,
                            discarded_by_contaminant1, discarded_by_contaminant2,
                            discarded_by_read_length1, discarded_by_read_length2,
                            pe_accept_cnt, pe_bases_kept, 
                            pe_discard_cnt,pe_bases_discarded, 
                            se_pe1_accept_cnt, se_pe1_bases_kept,
                            se_pe2_accept_cnt, se_pe2_bases_kept,
                            avg_trim_len_pe1, avg_trim_len_pe2,
                            avg_len_pe1, avg_len_pe2,
                            perfect_ov_cnt, partial_ov_cnt,
                            duplicates,
                            left_trimmed_by_polyat1, right_trimmed_by_polyat1,
                            right_trimmed_by_polyat2, right_trimmed_by_polyat2
                            );
    
    
    
    vector<string> t;
    split_str(st_str, t, "\n");
    for(int kk=0; kk<(int)t.size()+1; ++kk)
    {
       cout << "\033[A\033[2K";
       
    }
    t.clear();
    
    cout << st_str;
    sum_stat << st_str;
    
    
    sum_stat_tsv << PrintIlluminaStatisticsTSV(cnt1, cnt2, 
                            pe1_bases_anal, pe2_bases_anal, 
                            ts_adapters1, ts_adapters2, 
                            num_vectors1, num_vectors2, 
                            num_contaminants1, num_contaminants2, 
                            left_trimmed_by_quality1, left_trimmed_by_quality2,
                            left_trimmed_by_vector1, left_trimmed_by_vector2, 
                            avg_left_trim_len_pe1, avg_left_trim_len_pe2, 
                            right_trimmed_by_adapter1, right_trimmed_by_adapter2, 
                            right_trimmed_by_quality1,right_trimmed_by_quality2,
                            right_trimmed_by_vector1,right_trimmed_by_vector2,
                            avg_right_trim_len_pe1,avg_right_trim_len_pe2,
                            discarded1, discarded2,
                            discarded_by_contaminant1, discarded_by_contaminant2,
                            discarded_by_read_length1, discarded_by_read_length2,
                            pe_accept_cnt, pe_bases_kept, 
                            pe_discard_cnt,pe_bases_discarded, 
                            se_pe1_accept_cnt, se_pe1_bases_kept,
                            se_pe2_accept_cnt, se_pe2_bases_kept,
                            avg_trim_len_pe1, avg_trim_len_pe2,
                            avg_len_pe1, avg_len_pe2,
                            perfect_ov_cnt, partial_ov_cnt,
                            left_trimmed_by_polyat1, right_trimmed_by_polyat1,
                            right_trimmed_by_polyat2, right_trimmed_by_polyat2,
                            duplicates
                            ) 
                 << endl;
                
    cout << "====================Done cleaning====================\n";  
    sum_stat << "====================Done cleaning====================\n";  
    stat_str.clear();
    
    pe_output_file1.close();
    pe_output_file2.close();
    se_file.close();
    shuffle_file.close();
    
    rep_file1.close();
    rep_file2.close();
}

int IlluminaDynRoutine(Read* read, bool& adapter_found, string &query_str)
{
    //Remove not-needed Ns:
    //TrimNs( read->read );
    
    if((int)read->read.length() > minimum_read_length)
    {
        read->illumina_quality_string = read->illumina_quality_string.substr(0, read->read.length());
        read->clear_length = read->read.length();
    }
    else
    {
        read->discarded = 1;
        read->discarded_by_read_length = 1;
        return -1;
    }
    
    if(polyat_flag) {
       //If poly A/T flag is set:
       PolyAT_Trim(read);
    }
    
    if(contaminants_flag )
    {
       if(CheckContaminants(read->read) == 0) 
       {
           read->contam_found = 1;
           read->discarded_by_contaminant = 1;
           read->contaminants = 1;
           read->discarded = 1;
           return -1;
       }
    }
            
        
    //If quality trimming flag is set up -> perform the quality trimming before vector/contaminants/adaptors clipping.
    if( qual_trim_flag  ) 
    {
       QualTrimIllumina( read, max_a_error, max_e_at_ends );//This function generates LUCY clips of the read. Later they should be compared and read should be trimmed based on the results of comparison.
       if (read->discarded_by_quality == 1)
       {
          read->discarded = 1;
          return -1;
       }
       if (lucy_only_flag)
       {
           read->rclip = read->lucy_rclip;
           read->lclip = read->lucy_lclip;
           return 0;
       }
    }
    
    if( vector_flag ) 
            CheckVector(read); 
       
    //Run the main routine: Adapter + Vector/Contaminants trimming or only Adaptors
    //First 15 bases of i5 adapter forward
    if(trim_adapters_flag) {
        size_t found;
        if (!adapter_found)
        {
                string ts_adapter = tmpl_i5_1.substr(0,15);
                found = read->read.find( ts_adapter );
                if( found != string::npos ) 
                {
                        cout << "i5 adapter in forward first found in the read " << read->illumina_readID << ", in the position: " << found << endl;
                        sum_stat << "i5 adapter in forward first found in the read " << read->illumina_readID << ", in the position: " << found << endl;
                        adapter_found = true;
                        query_str = ts_adapter;
                        read->tru_sec_pos = found;
                        read->tru_sec_found = 1;
                }
                else
                {
                        //First 20 bases of i5 adapter in reverse complement
                        ts_adapter = MakeRevComplement(tmpl_i5_2).substr(0,15);
                        found = read->read.find( ts_adapter );
                        if( found != string::npos ) 
                        {
                                cout << "i5 adapter in forward first found in the read " << read->illumina_readID << ", in the position: " << found << endl;
                                sum_stat << "i5 adapter in forward first found in the read " << read->illumina_readID << ", in the position: " << found << endl;
                                adapter_found = true;
                                query_str = ts_adapter;
                                read->tru_sec_pos = found;
                                read->tru_sec_found = 1;
                        }
                        else
                        {
                                //First 20 bases of i7 adapter forward
                                ts_adapter = tmpl_i7_1.substr(0,15);
                                found = read->read.find( ts_adapter );
                                if( found != string::npos ) 
                                {
                                        cout << "i7 adapter in forward first found in the read " << read->illumina_readID << ", in the position: " << found << endl;
                                        sum_stat << "i7 adapter in forward first found in the read " << read->illumina_readID << ", in the position: " << found << endl;
                                        adapter_found = true;
                                        query_str = ts_adapter;
                                        read->tru_sec_pos = found;
                                        read->tru_sec_found = 1;
                                } 
                                else
                                {
                                        //First 20 bases of i5 adapter in reverse complement
                                        ts_adapter = MakeRevComplement(tmpl_i7_2).substr(0,15);
                                        found = read->read.find( ts_adapter );
                                        if( found != string::npos ) 
                                        {
                                                cout << "i7 adapter in forward first found in the read " << read->illumina_readID << ", in the position: " << found << endl;
                                                sum_stat << "i7 adapter in forward first found in the read " << read->illumina_readID << ", in the position: " << found << endl;
                                                adapter_found = true;
                                                query_str = ts_adapter;
                                                read->tru_sec_pos = found;
                                                read->tru_sec_found = 1;
                                        } 
                                        else
                                        {
                                                read->tru_sec_pos = -1;
                                                read->tru_sec_found = 0;
                                        }
                                }
                        }
                }
        }
        else
        {
                bool adp_found = false;
                found = read->read.rfind( query_str );
                if( found != string::npos ) 
                {
                        adp_found = true;
                        read->tru_sec_pos = found;
                        read->tru_sec_found = 1;
                } 
                else 
                {
                        //SSAHA job starts here
                        iz_SSAHA *izssaha = new iz_SSAHA();
                        AlignResult al_res = izssaha->Find( read->read , query_str );
                        AlignScores scores;
                        if( al_res.found_flag  ) 
                        {
                                scores = CalcScores(al_res.seq_1_al, al_res.seq_2_al, al_res.seq_1_al.length(), 0);
                                if(scores.mismatches <= max_al_mism  ) 
                                {
                                        adp_found = true;
                                        read->tru_sec_pos = al_res.pos;
                                        read->tru_sec_found = 1;
                                }
                        }
                        delete izssaha;
                }
        
                if(!adp_found) 
                {
                        read->tru_sec_pos = -1;
                        read->tru_sec_found = 0;
                }
        }
    }
    return 0;
}

void MakeClipPointsIllumina(Read* read) 
{
    //Clip points
   if( (qual_trim_flag ) && (vector_flag ) )
   {
        if(read->vector_found == 1)
        {
           if( read->v_start >= (int)(read->read.length() - read->v_end) ) //Vector is on the right side
           {
               read->lclip = read->lucy_lclip;
               
               if(read->lclip > 0)
                  read->left_trimmed_by_quality = 1;
               
               
               read->rclip = min(trim_adapters_flag ? (read->tru_sec_pos == -1 ? (int)read->read.length() : read->tru_sec_pos) : (int)read->read.length(), min(read->lucy_rclip, read->v_start) );
                    
               if( (read->rclip == read->lucy_rclip) && (read->rclip < read->initial_length ) )
               {
                 read->right_trimmed_by_quality = 1;
               }
               else if((read->rclip == read->tru_sec_pos) && trim_adapters_flag)
               {
                 if( (read->rclip < (int)read->read.length()) && (read->tru_sec_found == 1) && (read->rclip >= minimum_read_length))
                        read->right_trimmed_by_adapter = 1;
               }
               else if(read->rclip == read->v_start)
               {
                 read->right_trimmed_by_vector = 1;
               }
           }
           else //Vector is on the left side or the whole read is vector
           {
               read->lclip = max(read->lucy_lclip,read->v_end);//max(read->lucy_lclip,max(1, read->v_end ) );
           
               if( (read->lclip == read->lucy_lclip) && (read->lclip > 0) )//&& (read->lucy_lclip > 1)) 
               {
                 read->left_trimmed_by_quality = 1;
               }
               if(read->lclip == read->v_end)
               {
                 read->left_trimmed_by_vector = 1;
               }
           
               read->rclip = min(trim_adapters_flag ? (read->tru_sec_pos == -1 ? (int)read->read.length() : read->tru_sec_pos) : (int)read->read.length(), read->lucy_rclip );
               if( (read->rclip == read->lucy_rclip) && (read->rclip < read->initial_length ) )
               {
                 read->right_trimmed_by_quality = 1;
               }
               else
               {
                   if( (read->rclip < (int)read->read.length()) && (read->tru_sec_found == 1) && (read->rclip >= minimum_read_length) && trim_adapters_flag)
                        read->right_trimmed_by_adapter = 1;
               }
           }
           
        } 
        else
        {
            read->lclip = read->lucy_lclip;//max(read->lucy_lclip, 1);
            if(read->lclip > 0)
                read->left_trimmed_by_quality = 1;
            
            read->rclip = min(trim_adapters_flag ? (read->tru_sec_pos == -1 ? (int)read->read.length() : read->tru_sec_pos) : (int)read->read.length(), read->lucy_rclip );
            if( (read->rclip == read->lucy_rclip) && (read->rclip < read->initial_length ) )
            {
              read->right_trimmed_by_quality = 1;
            }
            else
            {
              if( (read->rclip < (int)read->read.length()) && (read->tru_sec_found == 1) && (read->rclip >= minimum_read_length) && trim_adapters_flag)
                        read->right_trimmed_by_adapter = 1;
            }
        }
    }
    else if( (qual_trim_flag ) && (!vector_flag) )
    {
        read->lclip = read->lucy_lclip;//max(read->lucy_lclip,1);
        if(read->lclip > 0)
           read->left_trimmed_by_quality = 1;
        
        read->rclip = min(trim_adapters_flag ? (read->tru_sec_pos == -1 ? (int)read->read.length() : read->tru_sec_pos) : (int)read->read.length(),read->lucy_rclip );
        if( (read->rclip == read->lucy_rclip) && (read->rclip < read->initial_length ) )
        {
           read->right_trimmed_by_quality = 1;
        }
        else if( (read->rclip == read->tru_sec_pos) && trim_adapters_flag)
        {
           
            if( (read->rclip < (int)read->read.length()) && (read->tru_sec_found == 1) && (read->rclip >= minimum_read_length))
                        read->right_trimmed_by_adapter = 1;
        }
        
        
                
        if(read->rclip >= read->clear_length)
        {
            read->rclip = read->clear_length; read->right_trimmed_by_adapter = 0;
        } 
        
        return;
        
        
    }
    else if( (!qual_trim_flag) && (vector_flag ) )
    {
       if( read->v_start >= ((int)read->read.length() - read->v_end) ) //Vector is on the right side
       {
          read->lclip = 0;
            
          read->rclip = min(trim_adapters_flag ? (read->tru_sec_pos == -1 ? (int)read->read.length() : read->tru_sec_pos) : (int)read->read.length(), read->v_start == -1 ? read->tru_sec_pos : read->v_start );
                    
          if( (read->rclip == (unsigned short)read->tru_sec_pos) && trim_adapters_flag)
          {
             if( (read->rclip < (int)read->read.length()) && (read->tru_sec_found == 1) && (read->rclip >= minimum_read_length))
                        read->right_trimmed_by_adapter = 1;
          }
          else if(read->rclip == read->v_start)
          {
             read->right_trimmed_by_vector = 1;
          }
           
       }
       else 
       {
          read->lclip = max(0, read->v_end == -1 ? 0 : read->v_end);
           
          if(read->lclip == read->v_end)
          {
             read->left_trimmed_by_vector = 1;
          }
           
          read->rclip = trim_adapters_flag ? (read->tru_sec_pos == -1 ? (int)read->read.length() : read->tru_sec_pos) : (int)read->read.length();
          if( (read->rclip < (int)read->read.length()) && (read->tru_sec_found == 1) && (read->rclip >= minimum_read_length) && trim_adapters_flag)
          {
                        read->right_trimmed_by_adapter = 1;
          }
       }
                
       if(read->rclip >= read->clear_length)
       {
          read->rclip = read->clear_length; read->right_trimmed_by_adapter = 0;
       }    
    }
    else if((!qual_trim_flag) && (!vector_flag))
    {
       read->rclip = trim_adapters_flag ? (read->tru_sec_pos == -1 ? (int)read->read.length() : read->tru_sec_pos) : (int)read->read.length();
       if( (read->rclip < (int)read->read.length()) && (read->tru_sec_found == 1) && (read->rclip >= minimum_read_length) && trim_adapters_flag)
          read->right_trimmed_by_adapter = 1;
            
       read->lclip = 0;
    }
    
    
    if(polyat_flag) {
       if( (read->rclip > read->poly_A_clip) && (read->poly_A_clip > 0)) {
          if(read->rclip == read->tru_sec_pos) {
             read ->right_trimmed_by_adapter = 0;
          } else if(read->rclip == read->lucy_rclip) {
             read->right_trimmed_by_quality = 0;
          } else if(read->rclip == read->v_start) {
             read->right_trimmed_by_vector = 0;
          }
                    
          read->rclip = read->poly_A_clip;
          read->right_trimmed_by_polyat = 1;
        }
        if( (read->lclip < read->poly_T_clip) && (read->poly_T_clip > 0)) {
          if(read->lclip == read->lucy_lclip) {
                read->left_trimmed_by_quality = 0;
          } else if(read->lclip == read->v_end) {
                read->left_trimmed_by_vector = 0;
          } 
             
             read->lclip = read->poly_T_clip;
             read->left_trimmed_by_polyat = 1;
          }
   }
   
   return;
    
}

string New2OldNbl(string header)
{
    vector <string> fields1, fields2;
    
    split_str( header, fields1, " " );
    split_str( fields1[0], fields2, ":" );
    
    string tmp = string(fields2[0] + "_" + fields2[2] + ":" + fields2[3] + ":" + fields2[4] + ":" + fields2[5] + ":" + fields2[6] + "#0/" + fields1[1].substr(0,1) + " (" + header + ")");
    
    fields2.clear();
    fields1.clear(); 
    
    return tmp;
    
}

void WriteSEOverlap(fstream &overlap_file, Read *read)
{
    overlap_file << read->illumina_readID << endl;
    overlap_file << read->read << endl;
    overlap_file << '+' << endl;
    overlap_file << read->illumina_quality_string << endl;
}


void WritePEFile(fstream &pe_output_file, Read *read)
{
    pe_output_file << read->illumina_readID << endl;
    pe_output_file << read->read << endl;
    pe_output_file << '+' << endl;
    pe_output_file << read->illumina_quality_string << endl;
}

void WriteShuffleFile(fstream &shuffle_output_file, Read *read1, Read *read2)
{
    shuffle_output_file << read1->illumina_readID << endl;
    shuffle_output_file << read1->read << endl;
    shuffle_output_file << '+' << endl;
    shuffle_output_file << read1->illumina_quality_string << endl;
    
    shuffle_output_file << read2->illumina_readID << endl;
    shuffle_output_file << read2->read << endl;
    shuffle_output_file << '+' << endl;
    shuffle_output_file << read2->illumina_quality_string << endl;
}

void WriteSEFile(fstream &se_output_file, Read *read)
{
    se_output_file << read->illumina_readID << endl;
    se_output_file << read->read << endl;
    se_output_file << '+' << endl;
    se_output_file << read->illumina_quality_string << endl;
    
    //fields1.clear(); fields2.clear();
}

//Dynamic Illumina: does not need space to store reads:
void IlluminaDynamicSE()
{
    se_bases_kept = se_bases_discarded = 0;
    se_discard_cnt = 0;
    se_bases_anal = 0;        
    avg_trim_len_se = 0;
    /*Raw implementation of average. Later I will come with a better algorithm*/
    unsigned long long avg_bases_se = 0;
    unsigned long long avg_left_clip = 0;
    unsigned long long avg_right_clip = 0;
    
    unsigned long cnt_avg; cnt_avg = 0; //Counters needed for calculating the average trimming length
    unsigned long cnt_avg_len; cnt_avg_len = 0;
                 
    double avg_len_se; avg_len_se = 0.0;
    double cnt_right_trim_se, avg_right_trim_len_se; 
    double cnt_left_trim_se, avg_left_trim_len_se;
    
    cnt_right_trim_se = avg_right_trim_len_se = 0;
    cnt_left_trim_se = avg_left_trim_len_se = 0;
    
    unsigned long cnt; cnt = 0;
    unsigned long se_accept_cnt; se_accept_cnt = 0;
    unsigned long ts_adapters; ts_adapters = 0;
    unsigned long num_vectors; num_vectors = 0;
    unsigned long num_contaminants; num_contaminants = 0;
    unsigned long accepted; accepted = 0;
    unsigned long discarded; discarded = 0;
//    unsigned long discarded_by_quality1, discarded_by_quality2; discarded_by_quality1 = discarded_by_quality2 = 0;
    unsigned long discarded_by_contaminant; discarded_by_contaminant = 0;
    unsigned long discarded_by_read_length; discarded_by_read_length = 0;
//    unsigned long discarded_by_vector1 , discarded_by_vector2; discarded_by_vector1 = discarded_by_vector2 = 0;
    /*Left trims*/
    unsigned long left_trimmed_by_quality; left_trimmed_by_quality = 0;
    unsigned long left_trimmed_by_vector; left_trimmed_by_vector = 0;
    /*Right trims/discards*/
    unsigned long right_trimmed_by_quality; right_trimmed_by_quality = 0;
    unsigned long right_trimmed_by_adapter; right_trimmed_by_adapter = 0;
    unsigned long right_trimmed_by_vector;  right_trimmed_by_vector = 0;
    
    unsigned long right_trimmed_by_polyat, left_trimmed_by_polyat; right_trimmed_by_polyat = left_trimmed_by_polyat = 0;
    unsigned long discarded_by_polyAT = 0;
    
    fstream rep_file, se_output_file;
    rep_file.open(rep_file_name1.c_str(),ios::out);
    rep_file << "ReadID\tlclip\trclip\tTruSeq_pos\tTruSeq_type\tRaw_read_length\tLlucy\tRlucy\tDiscarded\tContaminants\tVectorID\tVecStart\tVecEnd\tVecLen\n";
    
    cout << "Running the Illumina cleaning process..." << endl;
    sum_stat << "Running the Illumina cleaning process..." << endl;
    
    
    
    vector<string> record_block;
    
    
    
    se_output_file.open( se_output_filename.c_str(), ios::out );
    
    string st_str;
    //int first_avg = 0;
    for(int jj=0; jj<(int)se_names.size(); ++jj)
    {
        
        bool adapter_found = false;
        
        string query_string = "NA";
        
        int ii = 0;
        
        std::string line;
        igzstream in(/*fastq_file1*/se_names[jj]); //for R1
        
        cout << "Processing files: " << se_names[jj] << "\n";
        sum_stat << "Processing files: " << se_names[jj] << "\n";
        
        while ( getline(in,line) )
        {
                /*Read ID*/
                if(ii==0) 
                {
                    //Check for order
                    
                    if ( new2old_illumina && !old_style_illumina_flag )
                    {
                        vector <string> fields1, fields2;     
                        split_str( line, fields1, " " );
                        split_str( fields1[0], fields2, ":" );
                        line = string(fields2[0] + "_" + fields2[2] + ":" + fields2[3] + ":" + fields2[4] + ":" + fields2[5] + ":" + fields2[6] + "#0" ) ;
			
                        fields1.clear();
                        fields2.clear();
                    }
                    
                    record_block.push_back(line); 
                        
                    ii++;
                    continue;
                }
                /*DNA string*/
                if(ii==1) 
                {
                        record_block.push_back(line); /*DNA string*/
                        ii++;
                        continue;
                }
                /*a symbol "+"*/
                if(ii==2) 
                {
                        record_block.push_back(line);
                        ii++;
                        continue;
                }
                if(ii==3) 
                {
                        ii=0;
           
                        Read *read = new Read();
                        read->illumina_readID = record_block[0];
                        read->initial_length = record_block[1].length();
                        read->read = record_block[1];
                        read->illumina_quality_string = line;
                        se_bases_anal += read->read.length();
          
                        //Serial realization - useful for debugging if something does not work as expected
          
                        IlluminaDynRoutine(read, adapter_found, query_string);
                        
                        if(read->discarded == 0)
                        {
                           MakeClipPointsIllumina(read);
                        }
                        
                        cnt+=1;
          
                        //Report
                        //Read ID
                        rep_file << read->illumina_readID.substr(1,read->illumina_readID.length()-1) << "\t" << read->lclip << "\t" << read->rclip << "\t" << read->tru_sec_pos << "\t" << read->b_adapter << "\t" << read->initial_length << "\t" << (read->lucy_lclip <= 1 ? 1 : read->lucy_lclip) << "\t" << (read->lucy_rclip <= 1 ? 1 : read->lucy_rclip) << "\t" << read->discarded << "\t" << read->contaminants << "\t" << "NA" << "\n";
                        
                        
                                if( read->lclip >= read->rclip ) { read->discarded = 1; read->discarded_by_read_length = 1; } 
                        //cout <<  read->lclip << " " << read->rclip << endl;
                                if( read->lclip >= (int)read->read.length() ) { read->discarded = 1; read->discarded_by_read_length = 1; }
                                if( read->rclip > (int)read->read.length() ) { read->rclip = read->read.length(); }
                                if( (int)read->read.length() < minimum_read_length ) { read->discarded = 1; read->discarded_by_read_length = 1; }
                                if( (read->rclip - read->lclip) < minimum_read_length ) { read->discarded = 1; read->discarded_by_read_length = 1; }
              
                                if( read->discarded == 0 )
                                {
                                        
                                    if(  read->rclip < read->initial_length  )
                                    {
                                        cnt_right_trim_se += 1;
                                        avg_right_clip += read->initial_length - read->rclip;
                                        avg_right_trim_len_se = avg_right_clip/cnt_right_trim_se;
                                        //avg_right_trim_len_se = GetAvg( avg_right_trim_len_se, cnt_right_trim_se, read->initial_length - read->rclip );
                                    }
                                    if(read->lclip > 0)
                                    {
                                        cnt_left_trim_se += 1;
                                        //avg_left_trim_len_se = GetAvg( avg_left_trim_len_se, cnt_left_trim_se, read->lclip );
                                        avg_left_clip += read->lclip; 
                                        avg_left_trim_len_se = avg_left_clip/cnt_left_trim_se;
                                    }
                                    
                                    read->read = read->read.substr(0 , read->rclip );
                                    read->illumina_quality_string = read->illumina_quality_string.substr(0,read->rclip) ; 
                                    read->read = read->read.substr( read->lclip, read->rclip - read->lclip );
                                    read->illumina_quality_string = read->illumina_quality_string.substr( read->lclip, read->rclip - read->lclip );
                 
                                    WriteSEFile(se_output_file, read);
                                    se_accept_cnt+=1;
                                    se_bases_kept += read->read.length();
                                    
                                    //if( read->initial_length > (read->rclip - read->lclip) )
                                    if( read->discarded == 0 )
                                    {
                                        cnt_avg+=1;
                                        /*if (cnt_avg == 1) {
                                            first_avg = read->rclip - read->lclip;
                                            avg_trim_len_se = first_avg;
                                        } else {
                                            avg_trim_len_se = GetAvg( avg_trim_len_se, cnt_avg, read->rclip - read->lclip, first_avg );
                                            cout << avg_trim_len_se << endl;
                                        }*/
                                        avg_bases_se += read->rclip - read->lclip;
                                        avg_trim_len_se = avg_bases_se/cnt_avg;
                                        
                                    }
                 
                                    cnt_avg_len+=1; 
                                    //avg_len_se = GetAvg( avg_len_se, cnt_avg_len, read->read.length() );
                                    
                 
                                } 
                                
                         
                        
                        if (read->tru_sec_found == 1) ts_adapters++;
                        if (read->vector_found == 1) num_vectors++;
                        if (read->contam_found == 1) num_contaminants++;
                        if (read->discarded == 0) accepted++;
                        if (read->discarded == 1) discarded++;
                        //if (read1->discarded_by_quality == 1) discarded_by_quality1++;    
                        if (read->discarded_by_contaminant == 1) discarded_by_contaminant++;
                        if (read->discarded_by_read_length == 1) discarded_by_read_length++;
                        //if (read1->discarded_by_vector == 1) discarded_by_vector1++;
                        if (read->left_trimmed_by_quality == 1) left_trimmed_by_quality++;
                        if (read->left_trimmed_by_vector == 1) left_trimmed_by_vector++;
                        if (read->right_trimmed_by_quality == 1) right_trimmed_by_quality++;
                        if (read->right_trimmed_by_adapter == 1) right_trimmed_by_adapter++;
                        if (read->right_trimmed_by_vector == 1) right_trimmed_by_vector++;
                        if (read->right_trimmed_by_polyat == 1) right_trimmed_by_polyat++;
                        if (read->left_trimmed_by_polyat == 1)  left_trimmed_by_polyat++;
                        if(read->discarded_by_polyAT == 1) discarded_by_polyAT++;
                        
                        record_block.clear();
                        read->illumina_readID.clear(); 
                        read->illumina_quality_string.clear();
                        read->read.clear();
          
                        delete read;
                        
                        
                        if( (cnt % 1000 ) == 0)
                        {
                            st_str = PrintIlluminaStatisticsSE(cnt, 
                                    se_bases_anal, 
                                    ts_adapters,
                                    num_vectors,
                                    num_contaminants, 
                                    left_trimmed_by_quality,
                                    left_trimmed_by_vector, 
                                    avg_left_trim_len_se,
                                    right_trimmed_by_adapter,
                                    right_trimmed_by_quality,
                                    right_trimmed_by_vector,
                                    avg_right_trim_len_se,
                                    discarded, 
                                    discarded_by_contaminant,
                                    discarded_by_read_length,
                                    se_accept_cnt, se_bases_kept, 
                                    se_discard_cnt,se_bases_discarded, 
                                    avg_trim_len_se,
                                    avg_len_se,
                                    right_trimmed_by_polyat,
                                    left_trimmed_by_polyat,
                                    discarded_by_polyAT
                                   );
                            
                            if (cnt > 1000)
                            {
                                vector<string> t;
                                split_str(st_str, t, "\n");
                                for(int kk=0; kk<(int)t.size(); ++kk)
                                {
                                   cout << "\033[A\033[2K";
                                   //sum_stat << "\033[A\033[2K";
                                }
                                t.clear();
                            }
                            
                            cout << st_str;
                            
                        }
                        
                }
        }
        in.close();
        
    }
    
    
    st_str = PrintIlluminaStatisticsSE(cnt, 
                                    se_bases_anal, 
                                    ts_adapters,
                                    num_vectors,
                                    num_contaminants, 
                                    left_trimmed_by_quality,
                                    left_trimmed_by_vector, 
                                    avg_left_trim_len_se,
                                    right_trimmed_by_adapter,
                                    right_trimmed_by_quality,
                                    right_trimmed_by_vector,
                                    avg_right_trim_len_se,
                                    discarded, 
                                    discarded_by_contaminant,
                                    discarded_by_read_length,
                                    se_accept_cnt, se_bases_kept, 
                                    se_discard_cnt,se_bases_discarded, 
                                    avg_trim_len_se,
                                    avg_len_se,
                                    left_trimmed_by_polyat, 
                                    right_trimmed_by_polyat,
                                    discarded_by_polyAT
                                   );
    
    vector<string> t;
    split_str(st_str, t, "\n");
    for(int kk=0; kk<(int)t.size(); ++kk)
    {
        cout << "\033[A\033[2K";
        //sum_stat << "\033[A\033[2K";
    }
    t.clear();
    
    cout << st_str;
    sum_stat << st_str;
    
    sum_stat_tsv << PrintIlluminaStatisticsTSVSE(cnt,
                                    se_bases_anal, 
                                    ts_adapters, 
                                    num_vectors,  
                                    num_contaminants, 
                                    left_trimmed_by_quality, 
                                    left_trimmed_by_vector, 
                                    avg_left_trim_len_se, 
                                    right_trimmed_by_adapter, 
                                    right_trimmed_by_quality,
                                    right_trimmed_by_vector,
                                    avg_right_trim_len_se,
                                    discarded, 
                                    discarded_by_contaminant, 
                                    discarded_by_read_length,
                                    se_accept_cnt, 
                                   avg_trim_len_se,
                                   left_trimmed_by_polyat, 
                                   right_trimmed_by_polyat,
                                   discarded_by_polyAT
                            ) << endl;
                 
    
    cout << "====================Done cleaning====================\n";  
    sum_stat << "====================Done cleaning====================\n";  
    
    se_output_file.close();
    
    rep_file.close();
    
    
    
   
    
}


string PrintIlluminaStatistics(unsigned long cnt1, unsigned long cnt2, 
                                    unsigned long long  pe1_bases_anal, unsigned long long  pe2_bases_anal, 
                                    unsigned long ts_adapters1, unsigned long ts_adapters2, 
                                    unsigned long num_vectors1, unsigned long num_vectors2, 
                                    unsigned long num_contaminants1, unsigned long num_contaminants2, 
                                    unsigned long left_trimmed_by_quality1, unsigned long left_trimmed_by_quality2,
                                    unsigned long left_trimmed_by_vector1, unsigned long left_trimmed_by_vector2, 
                                    double avg_left_trim_len_pe1, double avg_left_trim_len_pe2, 
                                    unsigned long right_trimmed_by_adapter1, unsigned long right_trimmed_by_adapter2, 
                                    unsigned long right_trimmed_by_quality1,unsigned long right_trimmed_by_quality2,
                                    unsigned long right_trimmed_by_vector1,unsigned long right_trimmed_by_vector2,
                                    double avg_right_trim_len_pe1,double avg_right_trim_len_pe2,
                                    unsigned long discarded1, unsigned long discarded2,
                                    unsigned long discarded_by_contaminant1, unsigned long discarded_by_contaminant2,
                                    unsigned long discarded_by_read_length1, unsigned long discarded_by_read_length2,
                                    unsigned long pe_accept_cnt, unsigned long long  pe_bases_kept, 
                                    unsigned long pe_discard_cnt,unsigned long long  pe_bases_discarded, 
                                    unsigned long se_pe1_accept_cnt, unsigned long long se_pe1_bases_kept,
                                    unsigned long se_pe2_accept_cnt, unsigned long long se_pe2_bases_kept,
                                    double avg_trim_len_pe1, double avg_trim_len_pe2,
                                    double avg_len_pe1, double avg_len_pe2,
                                    unsigned long perfect_ov_cnt, unsigned long partial_ov_cnt,
                                    unsigned long duplicates,
                                    unsigned long left_trimmed_by_polyat1, unsigned long right_trimmed_by_polyat1,
                                    unsigned long left_trimmed_by_polyat2, unsigned long right_trimmed_by_polyat2
                                    )
{
    
     string stat_str = string("====================Summary Statistics====================\n");
     stat_str += string("PE1 reads analyzed: ") +  int2str(cnt1)   + string(", Bases:") +  int2str(pe1_bases_anal) + string("\n")  +
                        "Found ->\n" +
                        "Adapters: " + i2str(ts_adapters1,new char[15],10) + ", " + double2str( (double)ts_adapters1/(double)cnt1*100.0) + "%\n" + 
                        ( vector_flag ? "# of reads with vector: " + i2str(num_vectors1,new char[15],10) + ", " + double2str( (double)num_vectors1/(double)cnt1*100.0) + "%\n" : "") +
                        ( contaminants_flag ? "# of reads with contaminants: " + i2str(num_contaminants1,new char[15],10) + ", " + double2str( (double)num_contaminants1/(double)cnt1*100.0) + "%\n" : "") +
                        ( (qual_trim_flag || vector_flag) ? "Reads left trimmed ->\n" : "" ) +
                        ( qual_trim_flag ? "By quality: " +  i2str(left_trimmed_by_quality1,new char[15],10) + "\n" : "" ) +
                        ( vector_flag ? "By vector: " +  i2str(left_trimmed_by_vector1,new char[15],10) + "\n" : "" ) +
                        ( polyat_flag ? "By poly A/T: " + int2str(left_trimmed_by_polyat1) + "\n" : "") +
                        "Average left trim length: " + double2str(avg_left_trim_len_pe1) + " bp\n" +
                        "Reads right trimmed ->\n" +
                        "By adapter: " +  i2str(right_trimmed_by_adapter1,new char[15],10) + "\n" +
                        ( qual_trim_flag ? "By quality: " +  i2str(right_trimmed_by_quality1,new char[15],10) + "\n" : "") +
                        ( vector_flag ? "By vector: " +  i2str(right_trimmed_by_vector1,new char[15],10) + "\n" : "" ) +
                        ( polyat_flag ? "By poly A/T: " + int2str(right_trimmed_by_polyat1) + "\n" : "") +
                        "Average right trim length: " + double2str(avg_right_trim_len_pe1) + " bp\n" +
                        "PE1 reads discarded: " + i2str(discarded1,new char[15],10) + "\n" +
                        ( contaminants_flag ? "By contaminants: " +  i2str(discarded_by_contaminant1,new char[15],10) + "\n" : "" ) +
                        "By read length: " +  i2str(discarded_by_read_length1,new char[15],10) + "\n" +
                        "-----------------------------------------------------------\n" +
                        "PE2 reads analyzed: " + i2str(cnt2,new char[15],10) + ", Bases:" + i2str(pe2_bases_anal,new char[15],10) + "\n" +
                        "Found ->\n" + 
                        ("Adapters: " + i2str(ts_adapters2,new char[15],10) + ", " + double2str( (double)ts_adapters2/(double)cnt2*100.0) + "%\n") +
                        ( vector_flag ? ("# of reads with vector: " + i2str(num_vectors2,new char[15],10) + ", " + double2str( (double)num_vectors2/(double)cnt2*100.0) + "%\n") : "") +
                        ( contaminants_flag? "# of reads with contaminants: " + i2str(num_contaminants2,new char[15],10) + ", " + double2str( (double)num_contaminants2/(double)cnt2*100.0) + "%\n" : "") +
                        ( (qual_trim_flag || vector_flag ) ? "Reads left trimmed ->\n" : "" ) +
                        (qual_trim_flag ? "By quality: " +  i2str(left_trimmed_by_quality2,new char[15],10) + "\n" : "" ) +
                        ( vector_flag ? "By vector: " +  i2str(left_trimmed_by_vector2,new char[15],10) + "\n" : "" ) +
                        ( polyat_flag ? "By poly A/T: " + int2str(left_trimmed_by_polyat2) + "\n" : "") +
                        ("Average left trim length: " + double2str(avg_left_trim_len_pe2) + " bp\n" ) +
                        "Reads right trimmed ->\n" +
                        ( qual_trim_flag ? "By quality: " +  i2str(right_trimmed_by_quality2,new char[15],10) + "\n" : "") +
                        (vector_flag ? "By vector: " +  i2str(right_trimmed_by_vector2,new char[15],10) + "\n" : "" ) +
                        ( polyat_flag ? "By poly A/T: " + int2str(right_trimmed_by_polyat2) + "\n" : "") +                        
                        "By adapter: " +  i2str(right_trimmed_by_adapter2,new char[15],10) + "\n" +
                        ("Average right trim length: " + double2str(avg_right_trim_len_pe2) + " bp\n") +
                        "PE2 reads discarded:" + i2str(discarded2,new char[15],10) + "\n" +
                        (contaminants_flag ? "By contaminants: " +  i2str(discarded_by_contaminant2,new char[15],10) + "\n" : "" ) +
                        "By read length: " +  i2str(discarded_by_read_length2,new char[15],10) + "\n" + 
                        "----------------------Summary for PE & SE----------------------\n" +
                        ("Pairs kept: " + i2str(pe_accept_cnt,new char[15],10) + ", " + double2str( (double)pe_accept_cnt/(double)cnt1*100.0) + "%, Bases: " + i2str(pe_bases_kept,new char[15],10) + ", " + double2str( (double)pe_bases_kept/(double)(pe1_bases_anal+pe2_bases_anal)*100) +  "%\n") +
                        ("Pairs discarded: " + i2str(pe_discard_cnt,new char[15],10) + ", " + double2str( (double)pe_discard_cnt/(double)cnt1*100.0) + "%, Bases: " + i2str(pe_bases_discarded,new char[15],10) + ", " + double2str( (double)pe_bases_discarded/(double)(pe1_bases_anal+pe2_bases_anal)*100) +  "%\n") +
                        ("Single Reads PE1 kept: " + i2str(se_pe1_accept_cnt,new char[15],10) + ", Bases: " + i2str(se_pe1_bases_kept,new char[15],10) +"\n") +
                        ("Single Reads PE2 kept: " + i2str(se_pe2_accept_cnt,new char[15],10) + ", Bases: " + i2str(se_pe2_bases_kept,new char[15],10) +"\n") +
                        ("Average trimmed length PE1: " + double2str(avg_trim_len_pe1) + " bp\n") +
                        ("Average trimmed length PE2: " + double2str(avg_trim_len_pe2) + " bp\n") +
                        (overlap_flag ? "Perfect overlaps: " + int2str(perfect_ov_cnt) + "\n" : "") +
                        (overlap_flag ? "Partial overlaps: " + int2str(partial_ov_cnt) + "\n" : "") + 
                        (rem_dup ? "Duplicates: " + int2str(duplicates) + "\n" : "");
            
    
     return stat_str;
    
}

string PrintIlluminaStatisticsTSV(unsigned long cnt1, unsigned long cnt2, 
                                    unsigned long long  pe1_bases_anal, unsigned long long  pe2_bases_anal, 
                                    unsigned long ts_adapters1, unsigned long ts_adapters2, 
                                    unsigned long num_vectors1, unsigned long num_vectors2, 
                                    unsigned long num_contaminants1, unsigned long num_contaminants2, 
                                    unsigned long left_trimmed_by_quality1, unsigned long left_trimmed_by_quality2,
                                    unsigned long left_trimmed_by_vector1, unsigned long left_trimmed_by_vector2, 
                                    double avg_left_trim_len_pe1, double avg_left_trim_len_pe2, 
                                    unsigned long right_trimmed_by_adapter1, unsigned long right_trimmed_by_adapter2, 
                                    unsigned long right_trimmed_by_quality1,unsigned long right_trimmed_by_quality2,
                                    unsigned long right_trimmed_by_vector1,unsigned long right_trimmed_by_vector2,
                                    double avg_right_trim_len_pe1,double avg_right_trim_len_pe2,
                                    unsigned long discarded1, unsigned long discarded2,
                                    unsigned long discarded_by_contaminant1, unsigned long discarded_by_contaminant2,
                                    unsigned long discarded_by_read_length1, unsigned long discarded_by_read_length2,
                                    unsigned long pe_accept_cnt, unsigned long long pe_bases_kept, 
                                    unsigned long pe_discard_cnt,unsigned long long pe_bases_discarded, 
                                    unsigned long se_pe1_accept_cnt, unsigned long long se_pe1_bases_kept,
                                    unsigned long se_pe2_accept_cnt, unsigned long long se_pe2_bases_kept,
                                    double avg_trim_len_pe1, double avg_trim_len_pe2,
                                    double avg_len_pe1, double avg_len_pe2,
                                    unsigned long perfect_ov_cnt, unsigned long partial_ov_cnt,
                                    unsigned long left_trimmed_by_polyat1, unsigned long right_trimmed_by_polyat1,
                                    unsigned long left_trimmed_by_polyat2, unsigned long right_trimmed_by_polyat2,
                                    unsigned long duplicates
                                    )
{
    
        string filename_str;
    
        for(int i=0; i<(int)pe1_names.size(); ++i)
        {
            filename_str += string(pe1_names[i]) + ", " + string(pe2_names[i]);
        }
        
        string stat_str_tsv =   version + "\t" + 
                                filename_str + "\t" +
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
                                rep_file_name1+ "\t" +
                                rep_file_name2 +"\t" +
                                ( !shuffle_flag ?  pe_output_filename1 : "NA" ) +"\t"+
                                ( !shuffle_flag ?  pe_output_filename2 : "NA") +"\t"+
                                ( shuffle_flag ? shuffle_filename : "NA" ) +"\t" +
                                se_filename+ "\t" +
                                int2str(max_al_mism) +"\t" +
                                int2str(minimum_read_length)+ "\t" +
                                ( new2old_illumina ? "YES" : "NO") + "\t"; 
                   
    
    
        stat_str_tsv += int2str(cnt1)   + "\t" +  int2str(pe1_bases_anal) + "\t"  +
                       i2str(ts_adapters1,new char[15],10) + "\t" + double2str( (double)ts_adapters1/(double)cnt1*100.0) + "\t" + 
                       ( vector_flag ? i2str(num_vectors1,new char[15],10) + "\t" + double2str( (double)num_vectors1/(double)cnt1*100.0) + "\t" : "NA\tNA\t" ) +
                       ( contaminants_flag ? i2str(num_contaminants1,new char[15],10) + "\t" + double2str( (double)num_contaminants1/(double)cnt1*100.0) + "\t" : "NA\tNA\t" ) +
                       ( qual_trim_flag ? i2str(left_trimmed_by_quality1,new char[15],10) + "\t" : "NA\t" ) +
                       ( vector_flag ? i2str(left_trimmed_by_vector1,new char[15],10) + "\t" : "NA\t" ) +
                       double2str(avg_left_trim_len_pe1) + "\t" +
                       i2str(right_trimmed_by_adapter1,new char[15],10) + "\t" +
                       ( qual_trim_flag ? i2str(right_trimmed_by_quality1,new char[15],10) + "\t" : "NA\t") +
                       ( vector_flag ? i2str(right_trimmed_by_vector1,new char[15],10) + "\t" : "NA\t" ) +
                       double2str(avg_right_trim_len_pe1) + "\t" +
                       i2str(discarded1,new char[15],10) + "\t" +
                       ( contaminants_flag ? i2str(discarded_by_contaminant1,new char[15],10) + "\t" : "NA\t" ) +
                       i2str(discarded_by_read_length1,new char[15],10) + "\t" +
                       i2str(cnt2,new char[15],10) + "\t" + i2str(pe2_bases_anal,new char[15],10) + "\t" +
                       i2str(ts_adapters2,new char[15],10) + "\t" + double2str( (double)ts_adapters2/(double)cnt2*100.0) + "\t" +
                       ( vector_flag ? (i2str(num_vectors2,new char[15],10) + "\t" + double2str( (double)num_vectors2/(double)cnt2*100.0) + "\t") : "NA\tNA\t") +
                       ( contaminants_flag? i2str(num_contaminants2,new char[15],10) + "\t" + double2str( (double)num_contaminants2/(double)cnt2*100.0) + "\t" : "NA\tNA\t" ) +
                       (qual_trim_flag ? i2str(left_trimmed_by_quality2,new char[15],10) + "\t" : "NA\t" ) +
                       ( vector_flag ? i2str(left_trimmed_by_vector2,new char[15],10) + "\t" : "NA\t" ) +
                       double2str(avg_left_trim_len_pe2) + "\t"  +
                       ( qual_trim_flag ? i2str(right_trimmed_by_quality2,new char[15],10) + "\t" : "NA\t") +
                       (vector_flag ? i2str(right_trimmed_by_vector2,new char[15],10) + "\t" : "NA\t" ) +
                       i2str(right_trimmed_by_adapter2,new char[15],10) + "\t" +
                       double2str(avg_right_trim_len_pe2) + "\t" +
                       i2str(discarded2,new char[15],10) + "\t" +
                       (contaminants_flag ? i2str(discarded_by_contaminant2,new char[15],10) + "\t" : "NA\t" ) +
                       i2str(discarded_by_read_length2,new char[15],10) + "\t" + 
                       (i2str(pe_accept_cnt,new char[15],10) + "\t" + double2str( (double)pe_accept_cnt/(double)cnt1*100.0) + "\t" + i2str(pe_bases_kept,new char[15],10) + "\t" + double2str( (double)pe_bases_kept/(double)(pe1_bases_anal+pe2_bases_anal)*100) +  "\t") +
                       ( i2str(pe_discard_cnt,new char[15],10) + "\t" + double2str( (double)pe_discard_cnt/(double)cnt1*100.0) + "\t" + i2str(pe_bases_discarded,new char[15],10) + "\t" + double2str( (double)pe_bases_discarded/(double)(pe1_bases_anal+pe2_bases_anal)*100) +  "\t") +
                       (i2str(se_pe1_accept_cnt,new char[15],10) + "\t" + i2str(se_pe1_bases_kept,new char[15],10) +"\t") +
                       (i2str(se_pe2_accept_cnt,new char[15],10) + "\t" + i2str(se_pe2_bases_kept,new char[15],10) +"\t") +
                       (double2str(avg_trim_len_pe1) + "\t") +
                       double2str(avg_trim_len_pe2) +
                       (rem_dup ? "\t" + int2str(duplicates) : "\tNA") + 
                       ( overlap_flag ? "\t" + int2str(perfect_ov_cnt) + "\t" + int2str(partial_ov_cnt) : "\tNA\tNA") +
                       (polyat_flag ? "\tYES\t" + int2str(cdna) + "\t" + int2str(c_err) + "\t" + int2str(crng) + "\t" + int2str(left_trimmed_by_polyat1) + "\t" + int2str(right_trimmed_by_polyat1) + "\t" + int2str(left_trimmed_by_polyat2) + "\t" + int2str(right_trimmed_by_polyat2) : "\tNA\tNA\tNA\tNA\tNA\tNA\tNA");
                      
    return stat_str_tsv;
}

string PrintIlluminaStatisticsSE(unsigned long cnt, unsigned long long se_bases_anal, 
                                    unsigned long ts_adapters,
                                    unsigned long num_vectors,
                                    unsigned long num_contaminants, 
                                    unsigned long left_trimmed_by_quality,
                                    unsigned long left_trimmed_by_vector, 
                                    double avg_left_trim_len_se,
                                    unsigned long right_trimmed_by_adapter,
                                    unsigned long right_trimmed_by_quality,
                                    unsigned long right_trimmed_by_vector,
                                    double avg_right_trim_len_se,
                                    unsigned long discarded, 
                                    unsigned long discarded_by_contaminant,
                                    unsigned long discarded_by_read_length,
                                    unsigned long se_accept_cnt, unsigned long long se_bases_kept, 
                                    unsigned long se_discard_cnt,unsigned long long se_bases_discarded, 
                                    double avg_trim_len_se,
                                    double avg_len_se,
                                    unsigned long left_trimmed_by_polyat, unsigned long right_trimmed_by_polyat,
                                    unsigned long discarded_by_polyAT
                                    )
{
    
   
    
    string ans = "====================Summary Statistics====================\n" +
                        ("SE reads analyzed: " +  i2str(cnt,new char[15],10)  + ", Bases:" +  i2str(se_bases_anal, new char[25],10)  + "\n") +
                        "Found ->\n" +
                        "Adapters: " + i2str(ts_adapters,new char[15],10) + ", " + double2str( (double)ts_adapters/(double)cnt*100.0) + "%\n" + 
                        ( vector_flag ? "# of reads with vector: " + i2str(num_vectors,new char[15],10) + ", " + double2str( (double)num_vectors/(double)cnt*100.0) + "%\n" : "") +
                        ( contaminants_flag ? "# of reads with contaminants: " + i2str(num_contaminants,new char[15],10) + ", " + double2str( (double)num_contaminants/(double)cnt*100.0) + "%\n" : "") +
                        ( (qual_trim_flag || vector_flag) ? "Reads left trimmed ->\n" : "" ) +
                        ( qual_trim_flag ? "By quality: " +  i2str(left_trimmed_by_quality,new char[15],10) + "\n" : "" ) +
                        ( vector_flag ? "By vector: " +  i2str(left_trimmed_by_vector,new char[15],10) + "\n" : "" ) +
                        ( polyat_flag ? "By poly A/T: " + int2str(left_trimmed_by_polyat) + "\n": "") +
                        "Average left trim length: " + double2str(avg_left_trim_len_se) + " bp\n" +
                        "Reads right trimmed ->\n" +
                        "By adapter: " +  i2str(right_trimmed_by_adapter,new char[15],10) + "\n" +
                        ( qual_trim_flag ? "By quality: " +  i2str(right_trimmed_by_quality,new char[15],10) + "\n" : "") +
                        ( vector_flag ? "By vector: " +  i2str(right_trimmed_by_vector,new char[15],10) + "\n" : "" ) +
                        ( polyat_flag ? "By poly A/T: " + int2str(right_trimmed_by_polyat) + "\n" : "") +
                        "Average right trim length: " + double2str(avg_right_trim_len_se) + " bp\n" +
                        "SE reads discarded: " + i2str(discarded,new char[15],10) + "\n" +
                        ( contaminants_flag ? "By contaminants: " +  i2str(discarded_by_contaminant,new char[15],10) + "\n" : "" ) +
                        "By read length: " +  i2str(discarded_by_read_length,new char[15],10) + "\n" +
                        "----------------------Summary for SE----------------------\n" +
                        ("Reads kept: " + i2str(se_accept_cnt,new char[15],10) + ", " + double2str( (double)se_accept_cnt/(double)cnt*100.0) + "%, Bases: " + i2str(se_bases_kept,new char[15],10) + ", " + double2str( (double)se_bases_kept/(double)(se_bases_anal)*100) +  "%\n") +
                        ("Average trimmed length: " + double2str(avg_trim_len_se) + " bp\n");// +
                       // ("Average read length: " + double2str(avg_len_se) + " bp\n");
    
    return ans;
   
}


string PrintIlluminaStatisticsTSVSE(unsigned long cnt,
                                    unsigned long long se_bases_anal, 
                                    unsigned long ts_adapters, 
                                    unsigned long num_vectors,  
                                    unsigned long num_contaminants, 
                                    unsigned long left_trimmed_by_quality, 

                                    unsigned long left_trimmed_by_vector, 
                                    double avg_left_trim_len_se, 
                                    unsigned long right_trimmed_by_adapter, 
                                    unsigned long right_trimmed_by_quality,
                                    unsigned long right_trimmed_by_vector,
                                    double avg_right_trim_len_se,
                                    unsigned long discarded, 
                                    unsigned long discarded_by_contaminant, 
                                    unsigned long discarded_by_read_length,
                                    unsigned long se_accept_cnt, 
                                    double avg_trim_len_se,
                                    unsigned long left_trimmed_by_polyat, unsigned long right_trimmed_by_polyat,
                                    unsigned long discarded_by_polyAT
                                    )
{
    
        string filename_str;
    
        for(int i=0; i<(int)pe1_names.size(); ++i)
        {
            filename_str += string(se_names[i]);
        }
        
        string stat_str_tsv =   version + "\t" + 
                                filename_str + "\t" +
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
                                rep_file_name1+ "\t" +
                                ( !shuffle_flag ?  se_output_filename : "NA" ) +"\t"+
                                int2str(max_al_mism) +"\t" +
                                int2str(minimum_read_length)+ "\t" +
                                ( new2old_illumina ? "YES" : "NO") + "\t"; 
                   
    
    
        stat_str_tsv += int2str(cnt)   + "\t" + //reads analyzed
                        int2str(se_bases_anal) + "\t"  + //bases
                       i2str(ts_adapters,new char[15],10) + "\t" //adapters
                        + double2str( (double)ts_adapters/(double)cnt*100.0) + "\t" + //perc adapters
                       ( vector_flag ? i2str(num_vectors,new char[15],10) + "\t" + double2str( (double)num_vectors/(double)cnt*100.0) + "\t" : "NA\tNA\t" ) + //perc vectors
                       ( contaminants_flag ? i2str(num_contaminants,new char[15],10) + "\t" //cont
                        + double2str( (double)num_contaminants/(double)cnt*100.0) + "\t" : "NA\tNA\t" ) + //perc cont
                       ( qual_trim_flag ? i2str(left_trimmed_by_quality,new char[15],10) + "\t" : "NA\t" ) +  //left trimmed qual
                       ( vector_flag ? i2str(left_trimmed_by_vector,new char[15],10) + "\t" : "NA\t" ) + //left trimmed vect
                       double2str(avg_left_trim_len_se) + "\t" + //avg left trim len
                       i2str(right_trimmed_by_adapter,new char[15],10) + "\t" + 
                       ( qual_trim_flag ? i2str(right_trimmed_by_quality,new char[15],10) + "\t" : "NA\t") +
                       ( vector_flag ? i2str(right_trimmed_by_vector,new char[15],10) + "\t" : "NA\t" ) +
                       double2str(avg_right_trim_len_se) + "\t" +
                       i2str(discarded,new char[15],10) + "\t" + //discard
                       ( contaminants_flag ? i2str(discarded_by_contaminant,new char[15],10) + "\t" : "NA\t" ) +
                       i2str(discarded_by_read_length,new char[15],10) + "\t" +
                       i2str(se_accept_cnt,new char[15],10) + "\t" + //se reads kept
                        double2str( (double)se_accept_cnt/(double)cnt*100.0) + "\t" //perc kept
                        + i2str(se_bases_kept,new char[15],10) + "\t" + //bases kept
                        double2str( (double)se_bases_kept/(double)se_bases_anal*100.0) + "\t" + //%
                       double2str(avg_trim_len_se) +
                       (polyat_flag ? "\tYES\t" + int2str(cdna) + "\t" + int2str(c_err) + "\t" + int2str(crng) + "\t" + int2str(left_trimmed_by_polyat) + "\t" + int2str(right_trimmed_by_polyat)  + "\n" : "");
            
    
     return stat_str_tsv;
    
    
}

void RemoveContaminants(vector<Read*>& illumina_reads)
{
    
    for(unsigned int index = 0; index < illumina_reads.size(); index++)
    {
        if(contaminants_flag ) 
        {
                try 
                {//Check against contaminants :
                        if(CheckContaminants(illumina_reads[index]->read) == 0) 
                        {
                                illumina_reads[index]->discarded_by_contaminant = 1;
                                illumina_reads[index]->contaminants = 1;
                                illumina_reads[index]->discarded = 1;
                                
                        }
                }
                catch(exception& e)
                {
                        cout << e.what() << endl;
                }
        }
    }
}

void ClearNNs( vector<Read*>& reads ) 
{
    /*Clear NNs*/
    for(unsigned int ii = 0; ii < reads.size(); ii++) 
    {
        //cout << reads[ii].read.length() << " " << reads[ii].readID << endl;
        /*If Ns are are somewhere in the middle - completely discard the read:*/
         size_t found_Ns = ((reads[ii]->read).substr(0,reads[ii]->initial_length-100)).find("NNN");
         if( found_Ns != string::npos ) {
            /*If the read is discarded, its pair must be discarded also:*/
            //middle_nns_counter++;
            reads[ii]->discarded = 1;
            
            continue;
         } else {
                for(int s=0; s<10; s++) 
                {
                        TrimNs(reads[ii]->read);
                        reads[ii]->illumina_quality_string = reads[ii]->illumina_quality_string.substr(0, reads[ii]->read.length());
                        reads[ii]->clear_length = reads[ii]->read.length();
                }
         }
         
         
    }
}
