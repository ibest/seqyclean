/* 
 * File:   UnitTest.cpp
 * Author: kwt
 * 
 * Created on July 28, 2013, 11:18 AM
 */

#include "UnitTest.h"


UnitTest::UnitTest() {
     
    //Test Illumina PE
    
    //Test Illumina SE
    
    //Test overlap Illumins
}

UnitTest::UnitTest(const UnitTest& orig) {
}

UnitTest::~UnitTest() {
}

void UnitTest::Test454() {
    //Test 454
    roche_names.push_back("test_data/in.sff");
    trim_adapters_flag = true;
    qual_trim_flag = true;
    output_prefix = "UnitTest454";
            
    sum_stat_tsv << "Version\tFiles\tNUM_THREADS\tAdaptersTrimming\tVectorTrimming\tkmer_size\tDistance\tContamScr\tkmer_contam_size\tQualTrim\tQual_max_avg_error\tQual_max_err_at_ends\tOutputPrefix\tRepFilename\tOutputFilename\tMax_align_mismatches\tMinReadLen\tReadsAnalyzed\tBases\tleft_mid_tags_found\tpercentage_left_mid_tags_found\tright_mid_tags_found\tpercentage_right_mid_tags_found\tReadsWithVector_found\tpercentage_ReadsWithVector_found\tReadsWithContam_found\tpercentage_ReadsWithContam_found\tLeftTrimmedByAdapter\tLeftTrimmedByQual\tLeftTrimmedByVector\tAvgLeftTrimLen\tRightTrimmedByAdapter\tRightTrimmedByQual\tRightTrimmedByVector\tAvgRightTrimLen\tDiscardedTotal\tDiscByContam\tDiscByLength\tReadsKept\tPercentageKept\tAvgTrimmedLen\n";
    cout << "Unit test 454..." << endl;
    cout << "--------------------Basic parameters--------------------\n";
    sum_stat << "--------------------Basic parameters--------------------\n";
        
    cout << "Provided data file(s) : \n" ;
    sum_stat << "Provided data file(s) : \n" ;
       
    for(int i=0; i<(int)roche_names.size(); ++i)
    {
       cout << roche_names[i] << "\n";
       sum_stat << roche_names[i] << "\n";
    }
    
    cout << string("Adapters trimming: ") + ( trim_adapters_flag  ? "YES. Custom RLMIDS: " + ( custom_rlmids_flag  ? "YES, file_of_RL_MIDS: " + string(rlmids_file) : "NO. Default RL MIDs will be used." ) : " ")  << endl;
    sum_stat << string("Adapters trimming: ") + ( trim_adapters_flag  ? "YES. Custom RLMIDS: " + ( custom_rlmids_flag  ? "YES, file_of_RL_MIDS: " + string(rlmids_file) : "NO. Default RL MIDs will be used." ) : " " )  << endl;
        
    if(vector_flag)
    {
       cout << "Vector screening: YES. Vector_file provided: " << vector_file << endl;
       sum_stat << "Vector screening: YES. Vector_file provided: " << vector_file << endl;
       cout << "K-mer_size for for vector trimming: " <<  KMER_SIZE << " bp\n";
       sum_stat << "K-mer_size for vector trimming: " <<  KMER_SIZE << " bp\n";
       cout << "Distance between the first bases of two consequitve kmers: " << DISTANCE << " bp\n";
       sum_stat << "Distance between the first bases of two consequitve kmers: " << DISTANCE << " bp\n";
    
    } 
    else
    {
       cout << "Vector screening: NO" << endl;
       sum_stat << "Vector screening: NO" << endl;
    }
        
    if(contaminants_flag)
    {
       cout << "Contaminants screening: YES. File_of_contaminants: " << cont_file << endl;
       sum_stat << "Contaminants screening: YES. File_of_contaminants: " << cont_file << endl;
       cout << "K-mer size for contaminants: " << KMER_SIZE_CONT << " bp\n";
       sum_stat << "K-mer size for contaminants: " << KMER_SIZE_CONT << " bp\n";
    } 
    else
    {
       cout << "Contaminants screening: NO" << endl;
       sum_stat << "Contaminants screening: NO" << endl;
    }
        
    if(qual_trim_flag)
    {
           cout << "Quality trimming: YES" << endl;
           sum_stat << "Quality trimming: YES" << endl;
           cout << "Maximim error: " << -10*log10(max_a_error) << endl;
           sum_stat << "Maximim error: " << -10*log10(max_a_error) << endl;
           cout << "Maximim error at ends: " << -10*log10(max_e_at_ends) << endl;
           sum_stat << "Maximim error at ends: " << -10*log10(max_e_at_ends) << endl;
     }
     else
     {
           cout << "Quality trimming: NO" << endl;
           sum_stat << "Quality trimming: NO" << endl;
     }
        
     cout << "--------------------Output files--------------------\n";
     sum_stat << "--------------------Output files--------------------\n";
        
     cout << "Output prefix: " << output_prefix << endl;
     sum_stat << "Output prefix: " << output_prefix << endl;
        
     roche_output_file_name = output_prefix + ".sff";// + (output_fastqfile_flag ? ", " + output_prefix + ".fastq" : "" );
     roche_rep_file_name = output_prefix + "_Report.tsv" ;
        
     cout << "Report file: " << roche_rep_file_name << "\n";
     sum_stat << "Report file: " << roche_rep_file_name << "\n";
        
     cout << "Roche output file(s): " << ( roche_output_file_name + (output_fastqfile_flag ? ", " + output_prefix + ".fastq" : "" ) ) << "\n";
     sum_stat << "Roche output file(s): " << ( roche_output_file_name + (output_fastqfile_flag ? ", " + output_prefix + ".fastq" : "" ) ) << "\n";
        
        
     cout << "--------------------Other parameters--------------------\n";
     sum_stat << "--------------------Other parameters--------------------\n";
        
     cout << "k-mer_size: " <<  KMER_SIZE << " bp\n";
     sum_stat << "k-mer_size: " <<  KMER_SIZE << " bp\n";
    
     cout << "Maximum number of mismatches allowed in alignment: " <<  max_al_mism << endl;
     sum_stat << "Maximum number of mismatches allowed in alignment: " <<  max_al_mism << endl;
    
     cout << "Minimum read length to accept: " << minimum_read_length << " bp\n";
     sum_stat << "Minimum read length to accept: " << minimum_read_length << " bp\n";
        
     cout << "Distance between the first bases of two consequitve kmers: " << DISTANCE << " bp\n";
     sum_stat << "Distance between the first bases of two consequitve kmers: " << DISTANCE << " bp\n";
    
     cout << "number_of_threads: " << NUM_THREADS << endl;
     sum_stat << "number_of_threads: " << NUM_THREADS << endl;
     
    if(vector_flag ) 
    {
        BuildVectorDictionary(vector_file);
    }
    
    if(contaminants_flag ) 
    {
        /*Building dictionary for contaminants :*/
        BuildContDictionary(cont_file);
    }
     
    RocheRoutine();
       
}

void UnitTest::TestIlluminaPE() {
    cout << "--------------------Basic parameters--------------------\n";
    sum_stat << "--------------------Basic parameters--------------------\n";
        
    //Sum stat TSV header
    sum_stat_tsv << "Version\tPE1PE2\tAdapters_trim\tVectorTrim\tK_mer_size\tDistance\tContamScr\tkc_size\tQualTrim\tQual_max_avg_error\tQual_max_err_at_ends\tOutputPrefix\tRepFilename1\tRepFilename2\tPE1OutputFilename\tPE2OutputFilename\tShuffledFilename\tSEFilename\tMax_align_mismatches\tMinReadLen\tnew2old_illumina\tPE1ReadsAn\tPE1Bases\tPE1TruSeqAdap_found\tPerc_PE1TruSeq\tPE1ReadsWVector_found\tPerc_PE1ReadsWVector\tPE1ReadsWContam_found\tPerc_PE1ReadsWContam\tPE1LeftTrimmedQual\tPE1LeftTrimmedVector\tPE1AvgLeftTrimLen\tPE1RightTrimmedAdap\tPE1RightTrimmedQual\tPE1RightTrimmedVector\tPE1AvgRightTrimLen\tPE1DiscardedTotal\tPE1DiscByContam\tPE1DiscByLength\tPE2ReadsAn\tPE2Bases\tPE2TruSeqAdap_found\tPerc_PE2TruSeq\tPE2ReadsWVector_found\tPerc_PE2ReadsWVector\tPE2ReadsWContam_found\tPerc_PE2ReadsWContam\tPE2LeftTrimmedQual\tPE2LeftTrimmedVector\tPE2AvgLeftTrimLen\tPE2RightTrimmedAdap\tPE2RightTrimmedQual\tPE2RightTrimmedVector\tPE2AvgRightTrimLen\tPE2DiscardedTotal\tPE2DiscByContam\tPE2DiscByLength\tPairsKept\tPerc_Kept\tBases\tPerc_Bases\tPairsDiscarded\tPerc_Discarded\tBases\tPerc_Bases\tSE_PE1_Kept\tSE_PE1_Bases\tSE_PE2_Kept\tSE_PE2_Bases\tAvgTrimmedLenPE1\tAvgTrimmedLenPE2\tAvgLenPE1\tAvgLenPE2\tperfect_ov\tpartial_ov\n";
    
    cout << "Provided data files : " << endl;
    sum_stat << "Provided data files : " << endl;
    string filename_str = "";
    for(int i=0; i<(int)pe1_names.size(); ++i)
    {
         cout << "PE1: " << pe1_names[i] << ", PE2: " << pe2_names[i] << endl;
         sum_stat << "PE1: " << pe1_names[i] << ", PE2: " << pe2_names[i] << endl;
         filename_str += "\t" + string(pe1_names[i]) + "\t" + string(pe2_names[i]);
    }
        
    cout << "Adapters trimming: " << (trim_adapters_flag  ? "YES. " : "NO")  << endl;
    sum_stat << "Adapters trimming: " << (trim_adapters_flag  ? "YES. " : "NO")  << endl;
    
    if(vector_flag)
    {
           cout << "Vector screening: YES. Vector_file provided: " << vector_file << endl;
           sum_stat << "Vector screening: YES. Vector_file provided: " << vector_file << endl;
           cout << "K-mer_size for for vector trimming: " <<  KMER_SIZE << endl;
           sum_stat << "K-mer_size for vector trimming: " <<  KMER_SIZE << endl;
           cout << "Distance between the first bases of two consecutive kmers: " << DISTANCE << endl;
           sum_stat << "Distance between the first bases of two consecutive kmers: " << DISTANCE << endl;
    } 
    else
    {
           cout << "Vector screening: NO" << endl;
           sum_stat << "Vector screening: NO" << endl;
    }
        
        
    if(contaminants_flag)
    {
           cout << "Contaminants screening: YES. File_of_contaminants: " << cont_file << endl;
           sum_stat << "Contaminants screening: YES. File_of_contaminants: " << cont_file << endl;
           cout << "K-mer size for contaminants: " << KMER_SIZE_CONT << endl;
           sum_stat << "K-mer size for contaminants: " << KMER_SIZE_CONT << endl;
    } 
    else
    {
           cout << "Contaminants screening: NO" << endl;
           sum_stat << "Contaminants screening: NO" << endl;
    }
    if(qual_trim_flag)
    {
           cout << "Quality trimming: YES" << endl;
           sum_stat << "Quality trimming: YES" << endl;
           cout << "Maximum error: " << -10*log10(max_a_error) << endl;
           sum_stat << "Maximum error: " << -10*log10(max_a_error) << endl;
           cout << "Maximum error at ends: " << -10*log10(max_e_at_ends) << endl;
           sum_stat << "Maximum error at ends: " << -10*log10(max_e_at_ends) << endl;
    }
    else
    {
           cout << "Quality trimming: NO" << endl;
           sum_stat << "Quality trimming: NO" << endl;
    }
    cout << "--------------------Output files--------------------\n";
    sum_stat << "--------------------Output files--------------------\n";
        
    cout << "Output prefix: " << output_prefix << endl;
    sum_stat << "Output prefix: " << output_prefix << endl;
        
    rep_file_name1 = output_prefix + "_PE1_Report.tsv";
    rep_file_name2 = output_prefix + "_PE2_Report.tsv";
    pe_output_filename1 =  output_prefix + "_PE1.fastq" ;
    pe_output_filename2 =  output_prefix + "_PE2.fastq" ;
    shuffle_filename = output_prefix + "_shuffled.fastq";
    se_filename = output_prefix + "_SE.fastq";
                
    cout << "Report files: " << rep_file_name1 << ", " << rep_file_name2 << endl;
    sum_stat << "Report files: " << rep_file_name1 << ", " << rep_file_name2 << endl;
        
    if (!shuffle_flag)
    {
        cout << "PE1 file: " << pe_output_filename1 << endl;
        sum_stat << "PE1 file: " << pe_output_filename1 << endl;
        cout << "PE2 file: " << pe_output_filename2 << endl;
        sum_stat << "PE2 file: " << pe_output_filename2 << endl;
    } 
    else
    {
        cout << "Shuffled file: " << shuffle_filename << endl;
        sum_stat << "Shuffled file: " << shuffle_filename << endl;
    }    
    cout << "Single-end reads: "<< se_filename << endl;
    sum_stat << "Single-end reads: "<< se_filename << endl;
                
    if(overlap_flag) {
       overlap_file_name = output_prefix + "_SEOLP.fastq";
       sum_stat << "Single-end overlapped reads: "<< overlap_file_name << endl;
    }
        
    cout << "--------------------Other parameters--------------------\n";
    sum_stat << "--------------------Other parameters--------------------\n";
    cout << "Maximum number of mismatches allowed in alignment: " <<  max_al_mism << endl;
    sum_stat << "Maximum number of mismatches allowed in alignment: " <<  max_al_mism << endl;
    
    cout << "Minimum read length to accept: " << minimum_read_length << endl;
    sum_stat << "Minimum read length to accept: " << minimum_read_length << endl;
        
    cout << "New to old-style Illumina headers: " << (new2old_illumina == false ? "NO" : "YES") << endl;
    sum_stat << "New to old-style Illumina headers: " << (new2old_illumina == false ? "NO" : "YES") << endl;
                
    cout << "Old-style Illumina: " << (old_style_illumina_flag == false ? "NO" : "YES") << endl;
    sum_stat << "Old-style Illumina: " << (old_style_illumina_flag == false ? "NO" : "YES") << endl;
                
    cout << "Q-value: " << phred_coeff_illumina << endl;
    sum_stat << "Q-value: " << phred_coeff_illumina << endl;
    
}

void UnitTest::TestIlluminaSE() {
    sum_stat_tsv << "Version\tSE\tAdapters_trim\tVectorTrim\tK_mer_size\tDistance\tContamScr\tkc_size\tQualTrim\tQual_max_avg_error\tQual_max_err_at_ends\tOutputPrefix\tRepFilename\ttSEOutputFilename\tMax_align_mismatches\tMinReadLen\tnew2old_illumina\tSEReadsAn\tSEBases\tSETruSeqAdap_found\tPerc_SETruSeq\tSEReadsWVector_found\tPerc_SEReadsWVector\tSEReadsWContam_found\tPerc_SEReadsWContam\tSELeftTrimmedQual\tSELeftTrimmedVector\tSEAvgLeftTrimLen\tSERightTrimmedAdap\tSERightTrimmedQual\tSERightTrimmedVector\tSEAvgRightTrimLen\tSEDiscardedTotal\tSEDiscByContam\tSEDiscByLength\tSEReadsKept\tPerc_Kept\tBases\tPerc_Bases\tAvgTrimmedLenSE\n";
    
    cout << "Provided data file(s) : " << endl;
    sum_stat << "Provided data file(s) : " << endl;
    for(int i=0; i<(int)se_names.size(); ++i)
    {
       cout << "SE: " << se_names[i] << endl;
       sum_stat << "SE: " << se_names[i] << endl;
    }
        
    cout << "Adapters trimming: " << (trim_adapters_flag  ? "YES. " : "NO")  << endl;
    sum_stat << "Adapters trimming: " << (trim_adapters_flag  ? "YES. " : "NO")  << endl;
    
    if(vector_flag)
    {
         cout << "Vector screening: YES. Vector_file provided: " << vector_file << endl;
         sum_stat << "Vector screening: YES. Vector_file provided: " << vector_file << endl;
         cout << "K-mer_size for for vector trimming: " <<  KMER_SIZE << endl;
         sum_stat << "K-mer_size for vector trimming: " <<  KMER_SIZE << endl;
         cout << "Distance between the first bases of two consequitve kmers: " << DISTANCE << endl;
         sum_stat << "Distance between the first bases of two consequitve kmers: " << DISTANCE << endl;
    
    } 
    else
    {
         cout << "Vector screening: NO" << endl;
         sum_stat << "Vector screening: NO" << endl;
    }
        
        
    if(contaminants_flag)
    {
        cout << "Contaminants screening: YES. File_of_contaminants: " << cont_file << endl;
        sum_stat << "Contaminants screening: YES. File_of_contaminants: " << cont_file << endl;
        cout << "K-mer size for contaminants: " << KMER_SIZE_CONT << endl;
        sum_stat << "K-mer size for contaminants: " << KMER_SIZE_CONT << endl;
    } 
    else
    {
        cout << "Contaminants screening: NO" << endl;
        sum_stat << "Contaminants screening: NO" << endl;
    }
        
    if(qual_trim_flag)
    {
        cout << "Quality trimming: YES" << endl;
        sum_stat << "Quality trimming: YES" << endl;
        cout << "Maximim error: " << -10*log10(max_a_error) << endl;
        sum_stat << "Maximim error: " << -10*log10(max_a_error) << endl;
        cout << "Maximim error at ends: " << -10*log10(max_e_at_ends) << endl;
        sum_stat << "Maximim error at ends: " << -10*log10(max_e_at_ends) << endl;
    }
    else
    {
        cout << "Quality trimming: NO" << endl;
        sum_stat << "Quality trimming: NO" << endl;
    }
        
    cout << "--------------------Output files--------------------\n";
    sum_stat << "--------------------Output files--------------------\n";
        
    cout << "Output prefix: " << output_prefix << endl;
    sum_stat << "Output prefix: " << output_prefix << endl;
        
    rep_file_name1 = output_prefix + "_SE_Report.tsv";
    se_output_filename =  output_prefix + "_SE.fastq" ;
                
    cout << "Report file: " << rep_file_name1 << endl;
    sum_stat << "Report file: " << rep_file_name1<< endl;
        
    cout << "SE file: " << se_output_filename << endl;
    sum_stat << "SE file: " << se_output_filename << endl;
                    
    cout << "Single-end reads: "<< se_filename << endl;
    sum_stat << "Single-end reads: "<< se_filename << endl;
        
    cout << "--------------------Other parameters--------------------\n";
    sum_stat << "--------------------Other parameters--------------------\n";
    cout <//Test is the files provided are old-style illumina< "Maximum number of mismatches allowed in alignment: " <<  max_al_mism << endl;
    sum_stat << "Maximum number of mismatches allowed in alignment: " <<  max_al_mism << endl;
    
    cout << "Minimum read length to accept: " << minimum_read_length << endl;
    sum_stat << "Minimum read length to accept: " << minimum_read_length << endl;
        
    cout << "New to old-style Illumina headers: " << (new2old_illumina == false ? "NO" : "YES") << endl;
    sum_stat << "New to old-style Illumina headers: " << (new2old_illumina == false ? "NO" : "YES") << endl;
                
    cout << "Old-style Illumina: " << (old_style_illumina_flag == false ? "NO" : "YES") << endl;
    sum_stat << "Old-style Illumina: " << (old_style_illumina_flag == false ? "NO" : "YES") << endl;
                
    cout << "Q-value: " << phred_coeff_illumina << endl;
    sum_stat << "Q-value: " << phred_coeff_illumina << endl;
        
}