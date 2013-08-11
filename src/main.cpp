#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <streambuf>
#include <exception>
#include <pthread.h>
#include "timer.h"
#include "util.h"
#include "poly.h"
#include "Read.h"
#include "Dictionary.h"
#include "Report.h"
//#include "MainPipeLine.h"
#include "gzstream.h"
//#include "QualTrim.h"
#//include "flash.h"
#//include "dup.h"
#include "Roche.h"
#include "rlmid.h"
#include "Illumina.h"


using namespace std;

/*Computational parameters (default)*/
short KMER_SIZE = 15;
short DISTANCE = 1;
unsigned short NUM_THREADS = 4;

string version = "1.8.4 (2013-08-10)";

/*Data structures*/
vector<Read*> reads;
/*Left and right midtags*/
vector<RL_MID> rlmids;
/*Dictionaries*/
/*Vectors*/
map<string, vector<k_mer_struct> > VectorDict;
map<int, string > vector_names;
/*Contaminants*/
map<string, vector<k_mer_struct> > ContDict;
map<string, vector<k_mer_struct> >::iterator it_ContDict;
/*Map that holds a whole vector genomes :*/
map<long /*seq_id*/, string /*sequence*/ > VectorSeqs;


vector<Read*> reads_1;
vector<Read*> reads_2;

//For removal of duplicates
map<string, int > DupDict;

/*----------End of data structure definition------------------*/


/*----------Input data------------------------------------*/
string input_file_list;
/*Vectors : */
char *vector_file;
/*Contaminations : */
char *cont_file;
char* rlmids_file;// = "RL_MIDS.csv";
/*PCR flags and file with prinmers*/
bool pcr_flag = false;
char *pcr_file_name;

/*ROCHE*/
//char* roche_file_name;// = "";

/*LUCY flags and filename*/
//bool trim_lucy_flag = false;
//char *lucy_file_name;
/*Other parameters*/
bool contaminants_flag = false;
bool vector_flag = false;
bool amplicon_flag = false;
bool qual_trim_flag = false;
bool sff_file_flag = false;
bool fastq_file_flag = false;
bool debug_flag = false;
bool flag_454 = true;
bool roche_report_flag = false;
bool custom_rlmids_flag = false;
bool roche_flag = false;
bool output_flag = false;
bool custom_filename_flag = false;
bool gz_flag = false;
bool illumina_pe_flag = false;
bool illumina_se_flag = false;
bool polyat_flag = false;
bool trim_adapters_flag = true;

/*Illumina*/
bool illumina_flag = false;
bool illumina_flag_se = false;
char* illumina_file_name_R1;
char* illumina_file_name_R2;
char* illumina_file_name_se;
string adapter_type_R1;
string adapter_type_R2;
string query_str1;
string query_str2;


/*Maximim number of mismatches allowed in alignment operation*/
int max_al_mism = 5; /*5 mismatches by default*/

char *sffile_name;
/*----------End of input data definition------------------*/

/*Output data and parameters to observe the computational process*/
long line_counter = 0;
long counter = 0;
long discard_counter = 0;
long accept_counter = 0;
long trim_counter = 0;
long middle_nns_counter = 0;
long illumina_pe_counter = 0;

bool output_sfffile_flag = false;
bool output_fastqfile_flag = false;
char *output_file_name = (char*)"NA";
char* custom_output_filename;// = "";
bool keep_fastq_orig = false;
bool lucy_only_flag = false;
/*----------End of output data definition------------------*/

/*-----LUCY parameters------*/
float max_a_error = 0.01;
float max_e_at_ends = 0.01;

unsigned short minimum_read_length = 50;

/*Poly A/T trimming default parameters*/
unsigned short cdna = 10;
unsigned short c_err = 3;
unsigned short crng=50;
unsigned short keep;

/*Vector trimming*/
unsigned short L_limit = 1;
unsigned short R_limit = 1;
unsigned short vmr = 0;//15;
unsigned short vml = 0;//15;
unsigned short allowable_distance = 3;
unsigned short KMER_SIZE_CONT = 15;
unsigned short pmax = 2;

/*Other variables and parameters*/
std::ifstream read_file;

void PrintHelp();
//void ParseFastqFile(char *fastq_file, vector<Read*> &reads);
void ParseFastqFileIllumina(char* fastq_file, vector<Read*> &reads );
void ClearNNs( vector<Read*>& reads );
//void QualityTrimming( vector<Read*>& reads );
void IlluminaRoutine();
//void RocheRoutine();
void RemoveContaminants(vector<Read*>& illumina_reads);
//void RemoveContaminants454(vector<Read*>& reads454);
void PolyAT_Trim(Read* read);
void PolyATRoutine();
//void RocheRoutineDynamic();
void PolyATIlluminaRoutine();
void PolyATIlluminaRoutineSE();


string PrintIlluminaStatisticsPolyAT(unsigned long cnt1, unsigned long cnt2, 
                                    unsigned long long  pe1_bases_anal, unsigned long long  pe2_bases_anal, 
                                    unsigned long left_trimmed_by_quality1, unsigned long left_trimmed_by_quality2,
                                    double avg_left_trim_len_pe1, double avg_left_trim_len_pe2, 
                                    unsigned long right_trimmed_by_quality1,unsigned long right_trimmed_by_quality2,
                                    double avg_right_trim_len_pe1,double avg_right_trim_len_pe2,
                                    unsigned long discarded_by_read_length1, unsigned long discarded_by_read_length2,
                                    unsigned long pe_accept_cnt, unsigned long long  pe_bases_kept, 
                                    unsigned long pe_discard_cnt,unsigned long long  pe_bases_discarded, 
                                    unsigned long se_pe1_accept_cnt, unsigned long long se_pe1_bases_kept,
                                    unsigned long se_pe2_accept_cnt, unsigned long long se_pe2_bases_kept,
                                    double avg_trim_len_pe1, double avg_trim_len_pe2,
                                    unsigned long left_trimmed_by_polyat1, unsigned long left_trimmed_by_polyat2,
                                    unsigned long right_trimmed_by_polyat1, unsigned long right_trimmed_by_polyat2
                                    );

string PrintIlluminaStatisticsTSVPolyAT(unsigned long cnt1, unsigned long cnt2, 
                                    unsigned long long  pe1_bases_anal, unsigned long long  pe2_bases_anal, 
                                    unsigned long left_trimmed_by_quality1, unsigned long left_trimmed_by_quality2,
                                    double avg_left_trim_len_pe1, double avg_left_trim_len_pe2, 
                                    unsigned long right_trimmed_by_quality1,unsigned long right_trimmed_by_quality2,
                                    double avg_right_trim_len_pe1,double avg_right_trim_len_pe2,
                                    unsigned long discarded_by_read_length1, unsigned long discarded_by_read_length2,
                                    unsigned long pe_accept_cnt, unsigned long long pe_bases_kept, 
                                    unsigned long pe_discard_cnt,unsigned long long pe_bases_discarded, 
                                    unsigned long se_pe1_accept_cnt, unsigned long long se_pe1_bases_kept,
                                    unsigned long se_pe2_accept_cnt, unsigned long long se_pe2_bases_kept,
                                    double avg_trim_len_pe1, double avg_trim_len_pe2,
                                    unsigned long left_trimmed_by_polyat1, unsigned long left_trimmed_by_polyat2,
                                    unsigned long right_trimmed_by_polyat1, unsigned long right_trimmed_by_polyat2
                                    );

string PrintIlluminaStatisticsPolyATSE( unsigned long cnt, 
                                        unsigned long long se_bases_anal, 
                                        unsigned long left_trimmed_by_quality,
                                        double avg_left_trim_len_se,
                                        unsigned long right_trimmed_by_quality,
                                        double avg_right_trim_len_se,
                                        unsigned long discarded_by_read_length,
                                        unsigned long se_accept_cnt, unsigned long long se_bases_kept, 
                                        double avg_trim_len_se,
                                        unsigned long left_trimmed_by_polyat,
                                        unsigned long right_trimmed_by_polyat
                                    );

string PrintIlluminaStatisticsTSVPolyATSE(unsigned long cnt,
                                    unsigned long long se_bases_anal, 
                                    unsigned long left_trimmed_by_quality, 
                                    double avg_left_trim_len_se, 
                                    unsigned long right_trimmed_by_quality,
                                    double avg_right_trim_len_se,
                                    unsigned long discarded_by_read_length,
                                    unsigned long se_accept_cnt, 
                                    double avg_trim_len_se,
                                    unsigned long left_trimmed_by_polyat,
                                    unsigned long right_trimmed_by_polyat
                                   );

/*-------------------------------------*/

bool dynflag = false;

vector<string> file_list;

fstream sum_stat, sum_stat_tsv;

string output_prefix;

bool VectorOnlyFlag = false;
bool new2old_illumina = false;



bool serial_flag = false;

volatile int shared_var = 0;

bool shuffle_flag = false;

/*Report files*/
string rep_file_name1, rep_file_name2, pe_output_filename1, pe_output_filename2, shuffle_filename, se_filename, se_output_filename, overlap_file_name;



string roche_output_file_name = "";
string roche_rep_file_name = "";
char* polyat_file_name; 
string polyat_output_file_name;

unsigned long long se_bases_kept, se_bases_discarded;
unsigned long se_discard_cnt = 0;
unsigned long long se_bases_anal = 0;        
unsigned long avg_trim_len_se;

bool wildcart_search_flag = false;

vector<char*> pe1_names, pe2_names, roche_names, se_names;

string stat_str, tsv_stat_str;

int window0 = 50;
int window1 = 10;

bool old_style_illumina_flag = false;
int phred_coeff_illumina = 33; //by default assume new illumina (1.8)
bool i64_flag = false;

unsigned int adapterlength = 40;
double overlap_t = 0.9;
int minoverlap = 10;
bool overlap_flag = false;

bool overwrite_flag = false;

bool rem_dup = false;

int main(int argc, char *argv[]) 
{
    double start, finish, elapsed;
    GET_TIME(start);
    
    
    /*******************************************/
    /* Parse command line arguments */
    /*******************************************/
    if(argv[1] == NULL) {
        PrintHelp();
        return 0;
    }
    
    if( (string(argv[1]) == "-help") || (string(argv[1]) == "--help") || (string(argv[1]) == "-?") ) {
        PrintHelp();
        return 0;
    }
        
    for (int i=1; i<argc; i++) 
    {
        if( string(argv[i]) == "--version" ) 
        {
           cout << "Version: " << version << endl;
           exit(1);
        }
        if( string(argv[i]) == "-qual" ) 
        {
           qual_trim_flag = true;
           if ((i+1)<argc && isdigit(argv[i+1][0])) 
           {
               max_a_error = pow( 10 ,-1*((double)(atof(argv[++i])/10.0)) );//atof(argv[++i]);
               if ((i+1)<argc && isdigit(argv[i+1][0])) 
               {
                  max_e_at_ends = pow( 10 ,-1*((double)(atof(argv[++i])/10.0)) );//atof(argv[++i]);
                  if((i+1)<argc && (string(argv[i+1]) == "-w0") )
                  {
                        ++i;
                        if ((i+1)<argc && isdigit(argv[i+1][0])) 
                        {
                                window0 = atoi(argv[++i]);
                                if((i+1)<argc && (string(argv[i+1]) == "-w1") )
                                {
                                        ++i;
                                        window1 = atoi(argv[++i]);
                                }
                        } else {
                            cout << "Error: parameter w0 has empty value.\n";
                            PrintHelp();
                            return 0;
                        }
                  }
               } else if((i+1)<argc && (string(argv[i+1]) == "-w0") )
               {
                        ++i;
                        if ((i+1)<argc && isdigit(argv[i+1][0])) 
                        {
                                window0 = atoi(argv[++i]);
                                if((i+1)<argc && (string(argv[i+1]) == "-w1") )
                                {
                                        ++i;
                                        window1 = atoi(argv[++i]);
                                }
                        }
               } 
           } else if((i+1)<argc && (string(argv[i+1]) == "-w0") )
           {
               ++i;
               if ((i+1)<argc && isdigit(argv[i+1][0])) 
               {
                   window0 = atoi(argv[++i]);
                   if((i+1)<argc && (string(argv[i+1]) == "-w1") )
                   {
                       ++i;
                       window1 = atoi(argv[++i]);
                   } else {
                       cout << "Error: parameter w1 has empty value.\n";
                       PrintHelp();
                       return 0;
                   }
               } else {
                   cout << "Error: parameter w0 has empty value.\n";
                   PrintHelp();
                   return 0;
               }
           }
           
           continue;
        }
        if( string(argv[i]) == "--qual_only" ) 
        {
           lucy_only_flag = true;
           trim_adapters_flag = false;
           continue;
        } 
        if( string(argv[i]) == "--shuffle" ) 
        {
           shuffle_flag = true;
           continue;
        }
        if( string(argv[i]) == "--dup" ) 
        {
           rem_dup = true;
           continue;
        }
        if( string(argv[i]) == "--ow" ) 
        {
           overwrite_flag = true;
           continue;
        }
        if( string(argv[i]) == "--overlap" ) 
        {
           overlap_flag = true;
           if ( ( (i+1)<argc ) && (argv[i+1][0] != '-') ) 
           {
              minoverlap = atoi(argv[++i]);
           }
           continue;
        }
        if( string(argv[i]) == "-max_al_mism" ) 
        {
           max_al_mism = atoi(argv[i]);
           continue;
        }
        if( string(argv[i]) == "--keep_fastq_orig" ) 
        {
           keep_fastq_orig = true;
           continue;
        }
        if( string(argv[i]) == "-polyat" ) 
        {
           polyat_flag = true;
           trim_adapters_flag = false; 
           
           if ( ( (i+1)<argc ) && (argv[i+1][0] != '-') ) 
           {
                polyat_file_name = argv[++i];
           }
           
           if ((i+1)<argc && isdigit(argv[i+1][0])) 
           {
              cdna=atoi(argv[++i]);
              if ((i+1)<argc && isdigit(argv[i+1][0])) 
              {
                c_err=atoi(argv[++i]);
                if ((i+1)<argc && isdigit(argv[i+1][0])) 
                {
                    crng=atoi(argv[++i]);
                }
              }
           }
           continue;
        }
        if(string(argv[i]) == "-v" )
        {
           if ((i+1)<argc && !(isdigit(argv[i+1][0])) ) 
           {
             vector_file = argv[++i]; /*Vector file given*/
             vector_flag = true;
           }
           continue;
        }
        if(string(argv[i]) == "-c" )
        { 
           if ((i+1)<argc && !(isdigit(argv[i+1][0])) ) 
           {
             cont_file = argv[++i]; /*File with contaminants given*/
             contaminants_flag = true;
           }
           continue;
        }
        if(string(argv[i]) == "-m" )
        {
           if ( (i+1)<argc && !(isdigit(argv[i+1][0])) ) 
           {
             custom_rlmids_flag = true;
             rlmids_file = argv[++i]; /*Custom file with RL MIDS given*/
           }
           continue;
        }
        if(string(argv[i]) == "-adapter_length" )
        {
           if ( (i+1)<argc && !(isdigit(argv[i+1][0])) ) 
           {
             adapterlength = atoi(argv[++i]); 
           }
           continue;
        }
        if(string(argv[i]) == "-ot" )
        {
           if ( (i+1)<argc && !(isdigit(argv[i+1][0])) ) 
           {
             overlap_t = atof(argv[++i]);
           }
           continue;
        }
        if(string(argv[i]) == "-k" )
        {
           if ( (i+1)<argc && isdigit(argv[i+1][0]) ) 
           {
              KMER_SIZE = atoi(argv[++i]);
           } else {
               cout << "Error: parameter k has empty value.\n";
               PrintHelp();
               return 0;
           }
           continue;
        }
        if(string(argv[i]) == "-f" )
        {
           if ( (i+1)<argc && isdigit(argv[i+1][0]) ) 
           {
              DISTANCE = atoi(argv[++i]);
           }
           continue;
        }
        if(string(argv[i]) == "--test" )
        {
           if ( (i+1)<argc && isdigit(argv[i+1][0]) ) 
           {
              DISTANCE = atoi(argv[++i]);
           }
           continue;
        }
        if(string(argv[i]) == "-t" )
        {
           if ( (i+1)<argc && isdigit(argv[i+1][0]) ) 
           {
              NUM_THREADS = atoi(argv[++i]);
           } else {
               cout << "Error: parameter \'t\' has empty value.\n";
               PrintHelp();
               return 0;
           }
           continue;
        }
        if(string(argv[i]) == "-p" )
        {
           if ( (i+1)<argc && !(isdigit(argv[i+1][0])) ) 
           {
              pcr_flag = true;
              pcr_file_name = argv[++i]; 
           }
           continue;
        }
        if(string(argv[i]) == "-i64" )
        {
           i64_flag = true;
           
           continue;
        }
        if(string(argv[i]) == "--serial" )
        {
           serial_flag = true;
           continue;
        }
        if(string(argv[i]) == "-d" )
        {
           debug_flag = true;
           continue;
        }
        if(string(argv[i]) == "--new2old_illumina" )
        {
           new2old_illumina = true;
           continue;
        }
        if(string(argv[i]) == "-L_limit" )
        {
           if ( isdigit(argv[i+1][0]) ) 
                L_limit = atoi(argv[++i]);
           
           continue;
        }
        if(string(argv[i]) == "-R_limit" )
        {
           if ( isdigit(argv[i+1][0]) ) 
                R_limit = atoi(argv[++i]);
           
           continue;
        }
        if(string(argv[i]) == "-vmr" )
        {
           if ( isdigit(argv[i+1][0]) ) 
                vmr = atoi(argv[++i]);
           
           continue;
        }
        if(string(argv[i]) == "-vml" )
        {
           if ( isdigit(argv[i+1][0]) ) 
                vml = atoi(argv[++i]);
           
           continue;
        }
        if(string(argv[i]) == "-allowable_distance" )
        {
           if ( isdigit(argv[i+1][0]) ) {
                allowable_distance = atoi(argv[++i]);
           } else {
               cout << "Error: parameter \'allowable_distance\' has empty value.\n";
               PrintHelp();
               return 0;
           }
           
           if(allowable_distance > 50) allowable_distance = 15;
           
           continue; 
        }
        if(string(argv[i]) == "-minimum_read_length" )
        {
           if ( isdigit(argv[i+1][0]) ) {
                minimum_read_length = atoi(argv[++i]);
           } else {
               cout << "Error: parameter \'minimum_read_length\' has empty value.\n";
               PrintHelp();
               return 0;
           }
           
           if(minimum_read_length > 10000) {
               minimum_read_length = 50;
               cout << "warning: parameter minimum_read_length exceeded the maximum value of 10,000 bases and was set to 50 bases.\n";
           }
           
           continue;
        }
	if(string(argv[i]) == "-kc" )
        {
           if ( isdigit(argv[i+1][0]) ) {
                KMER_SIZE_CONT = atoi(argv[++i]);
           } else {
               cout << "Error: parameter \'kc\' has empty value.\n";
               PrintHelp();
               return 0;
           }
           
           continue;
        }
        if(string(argv[i]) == "--RLMIDS" )
        {
           cout << "Supported RL MIDS:\n";
           for(int i=0; i<36; i++) 
           {
              cout << mids[i] << endl;
           }
           
           return 0;
        }
        if(string(argv[i]) == "-o" ) 
        {
           if(argv[i+1] == NULL) 
           {
              cout << "Error: you have to specify the output prefix of output files." << endl;
              PrintHelp();
              return 0;
           } else
           {
               output_prefix = string(argv[++i]);
           }
        }
        if( string(argv[i]) == "--fastq" ) 
        {
           output_fastqfile_flag = true; //output file with cleaned reads in FASTQ format, for 454 mode only
           continue;
        }
        if(string(argv[i]) == "-?" )
        {
           PrintHelp();
           exit(1);
        }   
        if(string(argv[i]) == "-h" )
        {
           PrintHelp();
           exit(1);
        }
        if( string(argv[i]) == "-1" )
        {
           if ( ( (i+1)<argc ) && (argv[i+1][0] != '-') ) 
           {
              illumina_flag = true;
              illumina_file_name_R1 = argv[++i];
              pe1_names.push_back(illumina_file_name_R1);
              
              int jj=0;
              while( ( (i+1+jj)<argc ) && (argv[i+1+jj][0] != '-') )
              {
                  pe1_names.push_back(argv[i+jj+1]);
                  //printf("%s\n", argv[i+jj+1]);
                  jj+=1;
              }
              
           }
           
           continue;
        }
        if( string(argv[i]) == "-2" )
        {
           if ( ( (i+1)<argc ) && (argv[i+1][0] != '-') ) 
           {
              illumina_flag = true;
              illumina_file_name_R2 = argv[++i];
              pe2_names.push_back(illumina_file_name_R2);
              
              int jj=0;
              while( ( (i+1+jj)<argc ) && (argv[i+1+jj][0] != '-') )
              {
                  pe2_names.push_back(argv[i+1+jj]);
                  //printf("%s\n", argv[i+1+jj]);
                  jj+=1;
              }
              
           }
           
           continue;
        } 
        if( string(argv[i]) == "-U" ) //single-end file mode
        {
           if ( ( (i+1)<argc ) && (argv[i+1][0] != '-') ) 
           {
              illumina_se_flag = true;
              illumina_flag = true;
              illumina_file_name_se = argv[++i];
              se_names.push_back(illumina_file_name_se);
              
              int jj=0;
              while( ( (i+1+jj)<argc ) && (argv[i+1+jj][0] != '-') )
              {
                  se_names.push_back(argv[i+1+jj]);
                  //printf("%s\n", argv[i+1+jj]);
                  jj+=1;
              }
           }
           
           continue;
        }
        if( string(argv[i]) == "-454" ) 
        {
            if ( ( (i+1)<argc ) && (argv[i+1][0] != '-') ) 
            {
              roche_flag = true; 
              
              roche_names.push_back(argv[++i]);
              
              int jj=0;
              while( ( (i+1+jj)<argc ) && (argv[i+1+jj][0] != '-') )
              {
                  roche_names.push_back(argv[i+1+jj]);
                  jj+=1;
              }
              
            }
            
            continue;
        }
        if( string(argv[i]) == "--MakeRocheReport" ) 
        {
            MakeRocheReport((char*)"RocheClipPoints_report.tsv", argv[++i]);
            return 0;
        }
        if( string(argv[i]) == "-pmax" ) 
        {
            pmax = atoi(argv[++i]);
            continue;
        }
        if( string(argv[i]) == "--vonly" ) 
        {
            VectorOnlyFlag = true;
            continue;    
        }
    }
    
    /*******************************************/
    /* End of parsing command line arguments */
    /*******************************************/
    
    
    /*******************************************/
    /* Make desisions based on command line arguments */
    /*===============================================*/
    
    if(output_prefix == "")
    {
        cout << "No output prefix found.\n";
        PrintHelp();
        return 0;
    }
    
    
    //Check if input files exist
    if (illumina_flag)
    {
        if(illumina_se_flag)
        {
            for(int i=0; i<(int)se_names.size(); ++i)
            {
                if ( !exists( se_names[i] ) )
                {
                  cout<< "Error: file " <<  se_names[i] << " does not exist\n";
                  return 0;
                }
                
                //Test is the files provided are old-style illumina
                std::string line;
                igzstream in(se_names[i]); 
                getline(in,line);
                vector <string> fields;
                split_str( line, fields, ":" );
                if (fields.size() == 5)
                {
                    //Old headers == True
                    old_style_illumina_flag = true;
                    
                } else //if (fields.size() == 7)
                {
                    //Old headers == False
                    old_style_illumina_flag = false;
                }
                
            }
        } else 
        {
                if(pe1_names.size() != pe2_names.size())
                {
                        cout<< "Error: numbers of PE1 files and PE2 files do not match!\n";
                        return 0;
                } else
                {
                        for(int i=0; i<(int)pe1_names.size(); ++i)
                        {
        
                                if ( !exists( pe1_names[i] ) )
                                {
                                        cout<< "Error: file " <<  pe1_names[i] << " does not exist\n";
                                        return 0;
                                }
                                if (!exists( pe2_names[i] ) )
                                {
                                        cout<< "Error: file " <<  pe2_names[i] << " does not exist\n";
                                        return 0;
                                }
                                
                                
                                //Test is the files provided are old-style illumina
                                std::string line1, line2;
                                igzstream in1(pe1_names[i]); igzstream in2(pe2_names[i]); 
                                getline(in1,line1); getline(in2,line2);
                                vector <string> fields1, fields2;
                                split_str( line1, fields1, ":" ); split_str( line2, fields2, ":" );
                                if ( (fields1.size() == 5) && (fields2.size() == 5) )
                                {
                                        //Old headers == True
                                        old_style_illumina_flag = true;
                                        
                                } else if ( ( (fields1.size() == 9) && (fields2.size() == 9) ) || ((fields1.size() == 10) && (fields2.size() == 10)) )
                                {
                                        //Old headers == False
                                        old_style_illumina_flag = false;
                                } else if ( fields1.size() != fields2.size() )
                                {
                                    cout << "Error: impossible to have both old & new illumina as paired-end reads" << endl;
                                    return -1;
                                } else 
                                {
                                    cout << "Warning: unknown header format in file: " << pe1_names[i] << ", " << pe2_names[i] << endl;
                                    cout << "Header is:\n" << line1 << "(R1),\n" << line2 << "(R2)" << endl;
                                    old_style_illumina_flag = false;
                                }
                        }
                        
                }
        }
        
        
        if (i64_flag == true)
           phred_coeff_illumina = 64;
        
    }
    if(roche_flag)
    {
        for(int i=0; i<(int)roche_names.size(); ++i)
        {
            if (!exists( roche_names[i] ) )
            {
               cout<< "Error: file " <<  roche_names[i] << " does not exist\n";
               return 0;
            }
        }
    }
    
    /*if(polyat_flag && !illumina_flag)
    {
        if (!exists( polyat_file_name ) )
        {
                cout<< "Error: file " <<  polyat_file_name << " does not exist\n";
                return 0;
        }
    }**/
    
    if(vector_flag)
    {
        if (!exists( vector_file ) )
        {
                cout<< "Error: vector file provided " <<  vector_file << " does not exist\n";
                return 0;
        }
    }
    
    if(contaminants_flag)
    {
        if (!exists( cont_file ) )
        {
                cout<< "Error: contaminants file provided " <<  cont_file << " does not exist\n";
                return 0;
        }
    }
    
    //Check if output path exist
    vector<string> t;
    split_str(output_prefix, t, "/");
    string t_prefix = "";
    if(output_prefix[0] == '/') t_prefix += '/';
    for(int ii=0; ii < (int)t.size()-1; ii++)
    {
        t_prefix += t[ii] + "/";
    }
    if ( (t_prefix != "") && (!exists( (char*)t_prefix.c_str() ) ) )
    {
         cout<< "Warning: path " <<  t_prefix << " in output prefix does not exist" << endl;
         cout << "Trying to create...\n";
         if ( MakeDirectory(t_prefix) == 0 )
         {
             cout << "Sucess!\n";
         } 
         else
         {
             cout << "Could not created following path: " << t_prefix << "\n";
             return 0;
         }
    } else {
        if(illumina_flag || illumina_se_flag)
        {
                if(exists( (char*)(output_prefix + "_PE1.fastq").c_str() ) && !overwrite_flag)
                {
                    cout << "The output files you have specified already exist. Please delete these files or change your output file name and re-run SeqyClean." << endl;
                    return 0;
                }
                if(exists( (char*)(output_prefix + "_PE2.fastq").c_str() ) && !overwrite_flag)
                {
                    cout << "The output files you have specified already exist. Please delete these files or change your output file name and re-run SeqyClean." << endl;
                    return 0;
                }
        }
        if(roche_flag) 
        {
            if(exists( (char*)(output_prefix + ".sff").c_str() ) && !overwrite_flag)
            {
                cout << "The output files you have specified already exist. Please delete these files or change your output file name and re-run SeqyClean." << endl;
                return 0;
            }
            if(exists( (char*)(t_prefix + ".fastq").c_str() ) && output_fastqfile_flag && !overwrite_flag)
            {
                cout << "The output files you have specified already exist. Please delete these files or change your output file name and re-run SeqyClean." << endl;
                return 0;
            }
        }
    }
    
    t.clear();
    t_prefix.clear();        
   
    string sum_stat_filename = output_prefix + "_SummaryStatistics.txt";
    sum_stat.open(sum_stat_filename.c_str(),ios::out);// | ios::app);
    sum_stat_filename = output_prefix + "_SummaryStatistics.tsv";
    sum_stat_tsv.open(sum_stat_filename.c_str(),ios::out);
    
    if ( (!illumina_flag) && (!roche_flag) && (!polyat_flag) )
    {
        cout << "Not enough input parameters. Expected -454 <filename> or -1 <R1 filename>, -2 <R2 filename> arguments or -polyat <filename> arguments\n";
        sum_stat << "Not enough input parameters. Expected -454 <filename> or -1 <R1 filename>, -2 <R2 filename> arguments or -polyat <filename> arguments\n";
        sum_stat.close();
        
        sum_stat_tsv << "Not enough input parameters. Expected -454 <filename> or -1 <R1 filename>, -2 <R2 filename> arguments or -polyat <filename> arguments\n";
        sum_stat_tsv.close();
        
        PrintHelp();
        return 0;
    }
    
    if( (illumina_flag ) && (roche_flag ) ) 
    {
        cout << "Error! You can not use both 454 and Illumina modes in the same run. Use 454 or Illumina only, but not both!\n";
        sum_stat << "Error! You can not use both 454 and Illumina modes in the same run. Use 454 or Illumina only, but not both!\n";
        sum_stat.close();
        
        sum_stat_tsv << "Error! You can not use both 454 and Illumina modes in the same run. Use 454 or Illumina only, but not both!\n";
        sum_stat_tsv.close();
        
        PrintHelp();
        return 0;
    }
    
    //Printing parameters
    cout << "====================Parameters========================\n";
    sum_stat << "====================Parameters========================\n";
    
    
    cout << "Version: " << version << endl;
    sum_stat << "Version: " << version << endl;
    if(illumina_flag && !polyat_flag)
    {
        cout << "--------------------Basic parameters--------------------\n";
        sum_stat << "--------------------Basic parameters--------------------\n";
        
        //Sum stat TSV header
        
        
        if(!illumina_se_flag)
        {
        
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
        else
        {
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
                //cout <//Test is the files provided are old-style illumina< "Maximum number of mismatches allowed in alignment: " <<  max_al_mism << endl;
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
        
    }
    if(roche_flag & !polyat_flag)
    {
        sum_stat_tsv << "Version\tFiles\tNUM_THREADS\tAdaptersTrimming\tVectorTrimming\tkmer_size\tDistance\tContamScr\tkmer_contam_size\tQualTrim\tQual_max_avg_error\tQual_max_err_at_ends\tOutputPrefix\tRepFilename\tOutputFilename\tMax_align_mismatches\tMinReadLen\tReadsAnalyzed\tBases\tleft_mid_tags_found\tpercentage_left_mid_tags_found\tright_mid_tags_found\tpercentage_right_mid_tags_found\tReadsWithVector_found\tpercentage_ReadsWithVector_found\tReadsWithContam_found\tpercentage_ReadsWithContam_found\tLeftTrimmedByAdapter\tLeftTrimmedByQual\tLeftTrimmedByVector\tAvgLeftTrimLen\tRightTrimmedByAdapter\tRightTrimmedByQual\tRightTrimmedByVector\tAvgRightTrimLen\tDiscardedTotal\tDiscByContam\tDiscByLength\tReadsKept\tPercentageKept\tAvgTrimmedLen\n";
    
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
        
    }
    if(polyat_flag && roche_flag)
    {
        cout << "--------------------Basic parameters--------------------\n";
        sum_stat << "--------------------Basic parameters--------------------\n";
        
        cout << "Provided data file : \n";
        sum_stat << "Provided data file : \n";
        
        for(int i=0; i<(int)roche_names.size(); ++i)
        {
            cout << roche_names[i] << "\n";
            sum_stat << roche_names[i] << "\n";
        }
        
        cout << string("Poly A/T trimming: ") + ( "YES. CDNA = " + string(itoa(cdna, new char[5], 10)) + ", CERR = " + string(itoa(c_err, new char[5], 10)) + ", CRNG = " + string(itoa(crng, new char[5], 10)) ) << endl;
        sum_stat << string("Poly A/T trimming: ") + ( "YES. CDNA = " + string(itoa(cdna, new char[5], 10)) + ", CERR = " + string(itoa(c_err, new char[5], 10)) + ", CRNG = " + string(itoa(crng, new char[5], 10)) ) << endl;
        
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
        
        polyat_output_file_name = output_prefix;// + (output_fastqfile_flag ? ".fastq" : ".sff" );
    
        cout << "Poly A/T output file: " << polyat_output_file_name << "\n";
        sum_stat << "Poly A/T output file: " << polyat_output_file_name << "\n";
        
        
    }
    
    if(polyat_flag && illumina_flag)
    {
        if(!illumina_se_flag) 
        {
                sum_stat_tsv << "Version\tPE1PE2\tCDNA\tCERR\tCRNG\tQualTrim\tQual_max_avg_error\tQual_max_err_at_ends\tOutputPrefix\tRepFilename1\tRepFilename2\tPE1OutputFilename\tPE2OutputFilename\tShuffledFilename\tSEFilename\tMinReadLen\tnew2old_illumina\tPE1ReadsAn\tPE1Bases\tPE1LeftTrimmedPolyAT\tPE1LeftTrimmedQual\tPE1AvgLeftTrimLen\tPE1RightTrimmedPolyAT\tPE1RightTrimmedQual\tPE1AvgRightTrimLen\tPE1DiscByLength\tPE2ReadsAn\tPE2Bases\tPE2LeftTrimmedPolyAT\tPE2LeftTrimmedQual\tPE2AvgLeftTrimLen\tPE2RightTrimmedPolyAT\tPE2RightTrimmedQual\tPE2AvgRightTrimLen\tPE2DiscByLength\tPairsKept\tPerc_Kept\tBases\tPerc_Bases\tPairsDiscarded\tPerc_Discarded\tBases\tPerc_Bases\tSE_PE1_Kept\tSE_PE1_Bases\tSE_PE2_Kept\tSE_PE2_Bases\tAvgTrimmedLenPE1\tAvgTrimmedLenPE2\n";
            
                cout << "--------------------Basic parameters--------------------\n";
                sum_stat << "--------------------Basic parameters--------------------\n";
        
                cout << "Provided data files : \n";
                sum_stat << "Provided data files : \n";
                string filename_str = "";
        
                cout << pe1_names.size() << endl;
                for(int i=0; i<(int)pe1_names.size(); ++i)
                {
                        cout << "PE1: " << string(pe1_names[i]) << ", PE2: " << pe2_names[i] << endl;
                        sum_stat << "PE1: " << pe1_names[i] << ", PE2: " << pe2_names[i] << endl;
                        filename_str += "\t" + string(pe1_names[i]) + "\t" + string(pe2_names[i]);
                }
        
                cout << string("Poly A/T trimming: ") + ( "YES. CDNA = " + string(itoa(cdna, new char[5], 10)) + ", CERR = " + string(itoa(c_err, new char[5], 10)) + ", CRNG = " + string(itoa(crng, new char[5], 10)) ) << endl;
                sum_stat << string("Poly A/T trimming: ") + ( "YES. CDNA = " + string(itoa(cdna, new char[5], 10)) + ", CERR = " + string(itoa(c_err, new char[5], 10)) + ", CRNG = " + string(itoa(crng, new char[5], 10)) ) << endl;
        
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
        }
        else
        {
            sum_stat_tsv << "Version\tSE\tCDNA\tCERR\tCRNG\tQualTrim\tQual_max_avg_error\tQual_max_err_at_ends\tOutputPrefix\tRepFilename\tSEOutputFilename\tMinReadLen\tnew2old_illumina\tSEReadsAn\tSEBases\tSELeftTrimmedPolyAT\tSELeftTrimmedQual\tSEAvgLeftTrimLen\tSERightTrimmedPolyAT\tSERightTrimmedQual\tSEAvgRightTrimLen\tSEDiscByLength\tSEReadsKept\tPerc_Kept\tBases\tPerc_Bases\tAvgTrimmedLenSE\n";
    
            cout << "Provided data file(s) : " << endl;
            sum_stat << "Provided data file(s) : " << endl;
            for(int i=0; i<(int)se_names.size(); ++i)
            {
                 cout << "SE: " << se_names[i] << endl;
                 sum_stat << "SE: " << se_names[i] << endl;
            }
            
            cout << string("Poly A/T trimming: ") + ( "YES. CDNA = " + string(itoa(cdna, new char[5], 10)) + ", CERR = " + string(itoa(c_err, new char[5], 10)) + ", CRNG = " + string(itoa(crng, new char[5], 10)) ) << endl;
            sum_stat << string("Poly A/T trimming: ") + ( "YES. CDNA = " + string(itoa(cdna, new char[5], 10)) + ", CERR = " + string(itoa(c_err, new char[5], 10)) + ", CRNG = " + string(itoa(crng, new char[5], 10)) ) << endl;
        
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
        
            cout << "Minimum read length to accept: " << minimum_read_length << endl;
            sum_stat << "Minimum read length to accept: " << minimum_read_length << endl;
        
            cout << "New to old-style Illumina headers: " << (new2old_illumina == false ? "NO" : "YES") << endl;
            sum_stat << "New to old-style Illumina headers: " << (new2old_illumina == false ? "NO" : "YES") << endl;
                
            cout << "Old-style Illumina: " << (old_style_illumina_flag == false ? "NO" : "YES") << endl;
            sum_stat << "Old-style Illumina: " << (old_style_illumina_flag == false ? "NO" : "YES") << endl;
                
            cout << "Q-value: " << phred_coeff_illumina << endl;
            sum_stat << "Q-value: " << phred_coeff_illumina << endl;
        }
    }
    
    
    cout << "====================Starting the process===================="<< endl;
    sum_stat << "====================Starting the process===================="<< endl;
    
    // End of making desisions based on command line arguments
    
    /*---------------------Building dictionaries-----------------------------*/
    /*Vector dictionary*/
    if(vector_flag ) 
    {
        BuildVectorDictionary(vector_file);
    }
    
    if(pcr_flag ) 
    {
        /*Building dictionary for PCR :*/
        BuildVectorDictionary(pcr_file_name);
    }
    
    if(contaminants_flag ) 
    {
        /*Building dictionary for contaminants :*/
        BuildContDictionary(cont_file);
    }
    /*----------End of building the dictionaries------------------------*/
    
    if(illumina_flag && !polyat_flag) 
    {
        if(!illumina_se_flag)
        {
            IlluminaDynamic();
        }
        else
        {
            IlluminaDynamicSE();
        }
    }
    
    if( roche_flag && !polyat_flag )
    {
        RocheRoutine();
    }
    
    if( polyat_flag && roche_flag ) 
    {
        PolyATRoutine();
    }
    
    if( polyat_flag && illumina_flag) 
    {
        if(!illumina_se_flag)
        {
            PolyATIlluminaRoutine();
        } 
        else
        {
            PolyATIlluminaRoutineSE();
        }
    }
    
    cout << "Program finished.\n";
    sum_stat << "Program finished.\n";
    
    
    
    VectorSeqs.clear();
    VectorDict.clear();
    ContDict.clear();
    
    GET_TIME(finish);
    elapsed = finish - start;
    printf("Elapsed time = %e seconds\n", elapsed);
    sum_stat << "Elapsed time = " << elapsed << " seconds." << endl;
    sum_stat.close();
    sum_stat_tsv.close();
    output_prefix.clear();
    
}

void PolyATRoutine()
{
    long left_trimmed_by_polyat, right_trimmed_by_polyat, discarded_by_polyAT, bases_anal, accepted, discarded, left_trimmed_by_quality, right_trimmed_by_quality;
    left_trimmed_by_polyat = right_trimmed_by_polyat = discarded_by_polyAT = bases_anal = accepted = discarded = left_trimmed_by_quality = right_trimmed_by_quality = 0;
    
    if( string(roche_names[0]).substr( strlen(roche_names[0])-5, 5 ) == "fastq" ) 
    { /*FASTQ file given. Process it.*/
        ParseFastqFile(roche_names[0], reads);
        
    } 
    else if( string(roche_names[0]).substr( strlen(roche_names[0])-3, 3 ) == "sff" ) 
    {
        process_sff_to_fastq( roche_names[0], 0 );
    }
    else if( string(roche_names[0]).substr( strlen(roche_names[0])-2, 2 ) == "gz" ) 
    {
        ParseFastqFile(roche_names[0], reads);
    }
    else
    {
        cout << "Unknown file format\n";
        return;
    }
    
    if( qual_trim_flag  )
    {
         QualityTrimming(reads);
    }
    
    for(unsigned long i=0; i<reads.size(); i++)
    {
        bases_anal += reads[i]->initial_length;
        
        if(reads[i]->discarded == 0)
        {
            PolyAT_Trim(reads[i]);
            
            if(reads[i]->poly_T_clip > 0) 
            {
                left_trimmed_by_polyat += 1;
                reads[i]->lclip = reads[i]->lucy_lclip + reads[i]->poly_T_clip;
            }
            else
            {
                if(reads[i]->lucy_lclip > 1 ) 
                {
                    left_trimmed_by_quality += 1;
                    reads[i]->lclip = reads[i]->lucy_lclip;
                }
            }
            
            if(reads[i]->poly_A_clip > 0) 
            {
                right_trimmed_by_polyat += 1;
                reads[i]->rclip = reads[i]->lucy_rclip - reads[i]->poly_A_clip;
            }
            else
            {
                if(reads[i]->lucy_rclip > 1 ) 
                {
                    right_trimmed_by_quality += 1;
                    reads[i]->rclip = reads[i]->lucy_rclip;
                }
            }
            
            if ( (reads[i]->rclip - reads[i]->lclip < minimum_read_length) && (reads[i]->rclip > 1) && (reads[i]->lclip > 1) ) 
            {
                reads[i]->lucy_rclip = reads[i]->lucy_lclip = 0;
                reads[i]->discarded = 1;
                reads[i]->discarded_by_polyAT = 1;
                reads[i]->discarded_by_read_length = 1;
                
                discarded += 1;
            } else
            {
                accepted+=1;
            }
            
            
        } 
        else
        {
            discarded += 1;
        }
        
    }
    
    
    
    if( (string(roche_names[0]).substr( strlen(roche_names[0])-5, 5 ) == "fastq") || (string(roche_names[0]).substr( strlen(roche_names[0])-2, 2 ) == "gz")  )
    {
        WriteToFASTQ( polyat_output_file_name + ".fastq");
    } else {
        WriteToSFF( polyat_output_file_name + ".sff" );
        if( output_fastqfile_flag)
            WriteToFASTQ( polyat_output_file_name + ".fastq" );
    }
    
    
    cout << "====================Summary Statistics for Poly A/T====================\n";
    sum_stat << "====================Summary Statistics for Poly A/T====================\n";
    cout << "Reads analyzed: " << reads.size() << ", Bases:" << bases_anal << "\n";
    sum_stat << "Reads analyzed: " << reads.size() << ", Bases:" << bases_anal << "\n";
    
    cout << "Reads left trimmed ->\n ";
    sum_stat << "Reads left trimmed ->\n " ;
    
    cout << "By Poly A/T: " << left_trimmed_by_polyat << "\n";
    sum_stat << "By Poly A/T: " << left_trimmed_by_polyat << "\n";
    
    if(qual_trim_flag)
    {
        cout << "By quality: " <<  left_trimmed_by_quality << "\n";
        sum_stat << "By quality: " <<  left_trimmed_by_quality << "\n";
    }
    
    cout << "Reads right trimmed ->\n";
    sum_stat << "Reads right trimmed ->\n";
    
    cout << "By Poly A/T: " << right_trimmed_by_polyat << "\n";
    sum_stat << "By Poly A/T: " << right_trimmed_by_polyat << "\n";
    
    if(qual_trim_flag)
    {
        cout << "By quality: " <<  right_trimmed_by_quality << "\n";
        sum_stat << "By quality: " <<  right_trimmed_by_quality << "\n";
    }
    
    cout << "Reads discarded: " << discarded << "\n";
    sum_stat << "Reads discarded: " << discarded << "\n";
    
    cout << "--------------------------------------------------------\n";
    sum_stat << "--------------------------------------------------------\n";
    
    cout << "Reads accepted: " << accepted << ", %" << (double)accepted/(double)reads.size()*100 << "\n";
    sum_stat << "Reads accepted: " << accepted << ", %" << (double)accepted/(double)reads.size()*100 << "\n";
    
    
    cout << "==========================================================\n";
    sum_stat << "==========================================================\n";
}

void PolyATIlluminaRoutine()
{
    unsigned long long pe1_bases_anal, pe2_bases_anal, pe_bases_kept, pe_bases_discarded, se_pe1_bases_kept, se_pe2_bases_kept;
    unsigned long pe_discard_cnt;
    double avg_trim_len_pe1, avg_trim_len_pe2;
    
    
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
    unsigned long accepted1,accepted2; accepted1 = accepted2 = 0;
    unsigned long discarded1,discarded2; discarded1 = discarded2 = 0;
    unsigned long discarded_by_read_length1, discarded_by_read_length2; discarded_by_read_length1 = discarded_by_read_length2 = 0;
//    unsigned long discarded_by_vector1 , discarded_by_vector2; discarded_by_vector1 = discarded_by_vector2 = 0;
    /*Left trims*/
    unsigned long left_trimmed_by_quality1 , left_trimmed_by_quality2; left_trimmed_by_quality1 = left_trimmed_by_quality2 = 0;
    /*Right trims/discards*/
    unsigned long right_trimmed_by_quality1 , right_trimmed_by_quality2; right_trimmed_by_quality1 = right_trimmed_by_quality2 = 0;
    
    unsigned long left_trimmed_by_polyat1, left_trimmed_by_polyat2, right_trimmed_by_polyat1, right_trimmed_by_polyat2; 
    left_trimmed_by_polyat1 = left_trimmed_by_polyat2 = right_trimmed_by_polyat1 = right_trimmed_by_polyat2 = 0;
            
    fstream rep_file1, rep_file2, pe_output_file1, pe_output_file2, shuffle_file, se_file;
    rep_file1.open(rep_file_name1.c_str(),ios::out);
    rep_file2.open(rep_file_name2.c_str(),ios::out);
    rep_file1 << "ReadID\tlclip\trclip\tTruSeq_pos\tRaw_read_length\tLlucy\tRlucy\tDiscarded\tContaminants\tVecStart\tVecEnd\tVecLen\n";
    rep_file2 << "ReadID\tlclip\trclip\tTruSeq_pos\tRaw_read_length\tLlucy\tRlucy\tDiscarded\tContaminants\tVecStart\tVecEnd\tVecLen\n";
    
    cout << "Running Illumina poly A-T tails trimming process..." << endl;
    sum_stat << "Running Illumina poly A-T tails trimming process..." << endl;
    
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
                            line1 = string(fields2[0] + "_" + fields2[2] + ":" + fields2[3] + ":" + fields2[4] + ":" + fields2[5] + ":" + fields2[6] + "#0/" + fields1[1].substr(0,1) );//+ " (" + line1 + ")");
                            
                            fields1.clear();
                            fields2.clear();
                            
                            split_str( line2, fields1, " " );
                            split_str( fields1[0], fields2, ":" );
                            
                            line2 = string(fields2[0] + "_" + fields2[2] + ":" + fields2[3] + ":" + fields2[4] + ":" + fields2[5] + ":" + fields2[6] + "#0/" + fields1[1].substr(0,1) ); // + " (" + line2 + ")");
                            
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
                        
                        //If quality trimming flag is set up -> perform the quality trimming before vector/contaminants/adaptors clipping.
                        if( qual_trim_flag  ) 
                        {
                                QualTrimIllumina( read1, max_a_error, max_e_at_ends );//This function generates LUCY clips of the read. Later they should be compared and read should be trimmed based on the results of comparison.
                                if (read1->discarded_by_quality == 1) read1->discarded = 1;
                                 
                                QualTrimIllumina( read2, max_a_error, max_e_at_ends );
                                if (read2->discarded_by_quality == 1) read2->discarded = 1;
                        }
          
                        //Serial realization - useful for debugging if something does not work as expected
                        if(read1->discarded == 0)
                        {
                                PolyAT_Trim(read1);
            
                                if(read1->poly_T_clip > 0) 
                                {
                                        left_trimmed_by_polyat1 += 1;
                                        read1->lclip = read1->lucy_lclip + read1->poly_T_clip;
                                }
                                else
                                {
                                        if(read1->lucy_lclip > 1 ) 
                                        {
                                                left_trimmed_by_quality1 += 1;
                                                read1->lclip = read1->lucy_lclip;
                                        }
                                }
            
                                if(read1->poly_A_clip > 0) 
                                {
                                        right_trimmed_by_polyat1 += 1;
                                        read1->rclip = read1->lucy_rclip - read1->poly_A_clip;
                                }
                                else
                                {
                                        if(read1->lucy_rclip > 1 ) 
                                        {
                                                right_trimmed_by_quality1 += 1;
                                                read1->rclip = read1->lucy_rclip;
                                        }
                                }
            
                                if ( (read1->rclip - read1->lclip < minimum_read_length) && (read1->rclip > 1) && (read1->lclip > 1) ) 
                                {
                                        read1->lucy_rclip = read1->lucy_lclip = 0;
                                        read1->discarded = 1;
                                        read1->discarded_by_polyAT = 1;
                                        read1->discarded_by_read_length = 1;
                                }
                        } 
                        
                        if(read2->discarded == 0)
                        {
                                PolyAT_Trim(read2);
            
                                if(read2->poly_T_clip > 0) 
                                {
                                        left_trimmed_by_polyat2 += 1;
                                        read2->lclip = read2->lucy_lclip + read2->poly_T_clip;
                                }
                                else
                                {
                                        if(read2->lucy_lclip > 1 ) 
                                        {
                                                left_trimmed_by_quality2 += 1;
                                                read2->lclip = read2->lucy_lclip;
                                        }
                                }
            
                                if(read2->poly_A_clip > 0) 
                                {
                                        right_trimmed_by_polyat2 += 1;
                                        read2->rclip = read2->lucy_rclip - read2->poly_A_clip;
                                }
                                else
                                {
                                        if(read2->lucy_rclip > 1 ) 
                                        {
                                                right_trimmed_by_quality2 += 1;
                                                read2->rclip = read2->lucy_rclip;
                                        }
                                }
            
                                if ( (read2->rclip - read2->lclip < minimum_read_length) && (read2->rclip > 1) && (read2->lclip > 1) ) 
                                {
                                        read2->lucy_rclip = read2->lucy_lclip = 0;
                                        read2->discarded = 1;
                                        read2->discarded_by_polyAT = 1;
                                        read2->discarded_by_read_length = 1;
                                }
                        } 
                        
                        cnt1+=1; cnt2+=1;
          
                        if( (read1->discarded == 0) && (read2->discarded == 0) )
                        {
                                        
                                    if(  read1->rclip < read1->initial_length  )
                                    {
                                        cnt_right_trim_pe1 += 1;
                                        ////avg_right_trim_len_pe1 = GetAvg( avg_right_trim_len_pe1, cnt_right_trim_pe1, read1->initial_length - read1->rclip );
                                    }
                                    if(read1->lclip > 0)
                                    {
                                        cnt_left_trim_pe1 += 1;
                                        ////avg_left_trim_len_pe1 = GetAvg( avg_left_trim_len_pe1, cnt_left_trim_pe1, read1->lclip );
                                    }
                                    
                                    read1->read = read1->read.substr(0 , read1->rclip );
                                    read1->illumina_quality_string = read1->illumina_quality_string.substr(0,read1->rclip) ; 
                                    read1->read = read1->read.substr( read1->lclip, read1->rclip - read1->lclip );
                                    read1->illumina_quality_string = read1->illumina_quality_string.substr( read1->lclip, read1->rclip - read1->lclip );
                 
                                    if(  read2->rclip < read2->initial_length  )
                                    {
                                        cnt_right_trim_pe2 += 1;
                                        ///////////////avg_right_trim_len_pe2 = GetAvg( avg_right_trim_len_pe2, cnt_right_trim_pe2, read2->initial_length - read2->rclip );
                                    }
                                    if(read2->lclip > 0)
                                    {
                                        cnt_left_trim_pe2 += 1;
                                        ///////////////avg_left_trim_len_pe2 = GetAvg( avg_left_trim_len_pe2, cnt_left_trim_pe2, read2->lclip );
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
                                        //////////////avg_trim_len_pe1 = GetAvg( avg_trim_len_pe1, cnt1_avg,  read1->rclip - read1->lclip );//read1->initial_length - read1->read.length()
                                    }
                 
                                    if( read2->initial_length > (int)read2->read.length() )
                                    {
                                        cnt2_avg+=1;
                                        ////////////////avg_trim_len_pe2 = GetAvg( avg_trim_len_pe2, cnt2_avg, read2->rclip - read2->lclip );//read2->initial_length - read2->read.length()*
                                    }
                 
                                    cnt_avg_len1 +=1; cnt_avg_len2 +=1;
                 
                                    ////////////avg_len_pe1 = GetAvg( avg_len_pe1, cnt_avg_len1, read1->read.length() );
                                    //////////////avg_len_pe2 = GetAvg( avg_len_pe2, cnt_avg_len2, read2->read.length() );
                
                                        
                        } else if ((read1->discarded == 0) && (read2->discarded == 1)) 
                        {
                                    read1->read = read1->read.substr(0 , read1->rclip );
                                    read1->illumina_quality_string = read1->illumina_quality_string.substr(0,read1->rclip) ; 
                                    read1->read = read1->read.substr( read1->lclip, read1->rclip - read1->lclip );
                                    read1->illumina_quality_string = read1->illumina_quality_string.substr( read1->lclip, read1->rclip - read1->lclip );
                 
                                    vector <string> fields1, fields2;     
                                    split_str( read1->illumina_readID, fields1, " " );
                                    split_str( fields1[0], fields2, ":" );
                                    read1->illumina_readID = string(fields2[0] + "_" + fields2[2] + ":" + fields2[3] + ":" + fields2[4] + ":" + fields2[5] + ":" + fields2[6] + "#0" ) ; //+ " (" + line + ")");
                            
                                    fields1.clear();
                                    fields2.clear();                                 
                                    
                                    WriteSEFile( se_file, read1 );
                                    se_pe1_accept_cnt+=1;
                                    se_pe1_bases_kept += read1->read.length();
                 
                                    
                        } else if( (read1->discarded == 1) && (read2->discarded == 0) )
                        {
                                      
                                    read2->read = read2->read.substr(0 , read2->rclip );
                                    read2->illumina_quality_string = read2->illumina_quality_string.substr(0,read2->rclip) ; 
                                    read2->read = read2->read.substr( read2->lclip, read2->read.length() - read2->lclip );
                                    read2->illumina_quality_string = read2->illumina_quality_string.substr( read2->lclip, read2->illumina_quality_string.length() - read2->lclip );
        	 
                                    vector <string> fields1, fields2;     
                                    split_str( read2->illumina_readID, fields1, " " );
                                    split_str( fields1[0], fields2, ":" );
                                    read2->illumina_readID = string(fields2[0] + "_" + fields2[2] + ":" + fields2[3] + ":" + fields2[4] + ":" + fields2[5] + ":" + fields2[6] + "#0" ) ; //+ " (" + line + ")");
                            
                                    fields1.clear();
                                    fields2.clear();
                                    
                                    WriteSEFile( se_file, read2 );
                                    se_pe2_accept_cnt +=1;
                                    se_pe2_bases_kept += read2->read.length();
                 
                                    
                        }
                                else 
                                {
                                        pe_discard_cnt+=1;
                                        pe_bases_discarded += read1->read.length();
                                        pe_bases_discarded += read2->read.length();
                                }
                        
                        rep_file1 << read1->illumina_readID << "\t" << read1->lclip << "\t" << read1->rclip << "\t" << (read1->tru_sec_pos == -1 ? "NA" : int2str(read1->tru_sec_pos))  << "\t" << read1->initial_length << "\t" << (read1->lucy_lclip <= 1 ? 1 : read1->lucy_lclip) << "\t" << (read1->lucy_rclip <= 1 ? 1 : read1->lucy_rclip) << "\t" << read1->discarded << "\t" << read1->contaminants << "\t" << (vector_flag == true ? int2str(read1->v_start) : "NA") << "\t" << (vector_flag == true ? int2str(read1->v_end) : "NA") << "\t" << (vector_flag == true ? int2str(read1->vec_len) : "NA") << "\n";
                        rep_file2 << read2->illumina_readID << "\t" << read2->lclip << "\t" << read2->rclip << "\t" << (read2->tru_sec_pos == -1 ? "NA" : int2str(read2->tru_sec_pos)) << "\t"  << read2->initial_length << "\t" << (read2->lucy_lclip <= 1 ? 1 : read2->lucy_lclip) << "\t" << (read2->lucy_rclip <= 1 ? 1 : read2->lucy_rclip) << "\t" << read2->discarded << "\t" << read2->contaminants << "\t" << (vector_flag == true ? int2str(read2->v_start) : "NA") << "\t" << (vector_flag == true ? int2str(read2->v_end) : "NA") << "\t" << (vector_flag == true ? int2str(read2->vec_len) : "NA") << "\n";
          
          
                        if (read1->discarded == 0) accepted1++;
                        if (read1->discarded == 1) discarded1++;
                        if (read1->discarded_by_read_length == 1) discarded_by_read_length1++;
                        if (read1->left_trimmed_by_quality == 1) left_trimmed_by_quality1++;
                        if (read1->right_trimmed_by_quality == 1) right_trimmed_by_quality1++;
                        
                        if (read2->discarded == 0) accepted2++;
                        if (read2->discarded == 1) discarded2++;
                        if (read2->discarded_by_read_length == 1) discarded_by_read_length2++;
                        if (read2->left_trimmed_by_quality == 1) left_trimmed_by_quality2++;
                        if (read2->right_trimmed_by_quality == 1) right_trimmed_by_quality2++;
                        
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
                            st_str = PrintIlluminaStatisticsPolyAT(cnt1, cnt2, 
                                    pe1_bases_anal, pe2_bases_anal, 
                                    left_trimmed_by_quality1, left_trimmed_by_quality2,
                                    avg_left_trim_len_pe1, avg_left_trim_len_pe2, 
                                    right_trimmed_by_quality1,right_trimmed_by_quality2,
                                    avg_right_trim_len_pe1,avg_right_trim_len_pe2,
                                    discarded_by_read_length1, discarded_by_read_length2,
                                    pe_accept_cnt, pe_bases_kept, 
                                    pe_discard_cnt,pe_bases_discarded, 
                                    se_pe1_accept_cnt, se_pe1_bases_kept,
                                    se_pe2_accept_cnt, se_pe2_bases_kept,
                                    avg_trim_len_pe1, avg_trim_len_pe2,
                                    left_trimmed_by_polyat1, left_trimmed_by_polyat2,
                                    right_trimmed_by_polyat1, right_trimmed_by_polyat2
                                    
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
    
    
    st_str = PrintIlluminaStatisticsPolyAT(cnt1, cnt2, 
                            pe1_bases_anal, pe2_bases_anal, 
                            left_trimmed_by_quality1, left_trimmed_by_quality2,
                            avg_left_trim_len_pe1, avg_left_trim_len_pe2, 
                            right_trimmed_by_quality1,right_trimmed_by_quality2,
                            avg_right_trim_len_pe1,avg_right_trim_len_pe2,
                            discarded_by_read_length1, discarded_by_read_length2,
                            pe_accept_cnt, pe_bases_kept, 
                            pe_discard_cnt,pe_bases_discarded, 
                            se_pe1_accept_cnt, se_pe1_bases_kept,
                            se_pe2_accept_cnt, se_pe2_bases_kept,
                            avg_trim_len_pe1, avg_trim_len_pe2,
                            left_trimmed_by_polyat1, left_trimmed_by_polyat2,
                            right_trimmed_by_polyat1, right_trimmed_by_polyat2
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
    
    
    sum_stat_tsv << PrintIlluminaStatisticsTSVPolyAT(cnt1, cnt2, 
                            pe1_bases_anal, pe2_bases_anal, 
                            left_trimmed_by_quality1, left_trimmed_by_quality2,
                            avg_left_trim_len_pe1, avg_left_trim_len_pe2, 
                            right_trimmed_by_quality1,right_trimmed_by_quality2,
                            avg_right_trim_len_pe1,avg_right_trim_len_pe2,
                            discarded_by_read_length1, discarded_by_read_length2,
                            pe_accept_cnt, pe_bases_kept, 
                            pe_discard_cnt,pe_bases_discarded, 
                            se_pe1_accept_cnt, se_pe1_bases_kept,
                            se_pe2_accept_cnt, se_pe2_bases_kept,
                            avg_trim_len_pe1, avg_trim_len_pe2,
                            left_trimmed_by_polyat1, left_trimmed_by_polyat2,
                            right_trimmed_by_polyat1, right_trimmed_by_polyat2
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

void PolyATIlluminaRoutineSE()
{
    se_bases_kept = se_bases_discarded = 0;
    se_discard_cnt = 0;
    se_bases_anal = 0;        
    avg_trim_len_se = 0;
    
    unsigned long cnt_avg; cnt_avg = 0; //Counters needed for calculating the average trimming length
    unsigned long cnt_avg_len; cnt_avg_len = 0;
                 
    double avg_len_se; avg_len_se = 0.0;
    double cnt_right_trim_se, avg_right_trim_len_se; 
    double cnt_left_trim_se, avg_left_trim_len_se;
    
    cnt_right_trim_se = avg_right_trim_len_se = 0;
    cnt_left_trim_se = avg_left_trim_len_se = 0;
    
    unsigned long cnt; cnt = 0;
    unsigned long se_accept_cnt; se_accept_cnt = 0;
    unsigned long accepted; accepted = 0;
    unsigned long discarded; discarded = 0;
   
    unsigned long discarded_by_read_length; discarded_by_read_length = 0;
    //*Left trims*/
    unsigned long left_trimmed_by_quality; left_trimmed_by_quality = 0;
    /*Right trims/discards*/
    unsigned long right_trimmed_by_quality; right_trimmed_by_quality = 0;

    unsigned long left_trimmed_by_polyat, right_trimmed_by_polyat; left_trimmed_by_polyat = right_trimmed_by_polyat = 0;
    
    fstream rep_file, se_output_file;
    rep_file.open(rep_file_name1.c_str(),ios::out);
    rep_file << "ReadID\tlclip\trclip\tTruSeq_pos\tTruSeq_type\tRaw_read_length\tLlucy\tRlucy\tDiscarded\tContaminants\tVectorID\tVecStart\tVecEnd\tVecLen\n";
    
    cout << "Running the Illumina cleaning process..." << endl;
    sum_stat << "Running the Illumina cleaning process..." << endl;
    
    
    
    vector<string> record_block;
    
    
    
    se_output_file.open( se_output_filename.c_str(), ios::out );
    
    string st_str;
    
    for(int jj=0; jj<(int)se_names.size(); ++jj)
    {
    
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
                        line = string(fields2[0] + "_" + fields2[2] + ":" + fields2[3] + ":" + fields2[4] + ":" + fields2[5] + ":" + fields2[6] + "#0" ) ; //+ " (" + line + ")");
                            
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
          
                        if(read->initial_length <= minimum_read_length)
                        {
                            cout << "Warming: in the given single-end file raw read length is less or equal then minimum_read_length\n" ; 
                            sum_stat << "Warming: in the given single-end file raw read length is less or equal then minimum_read_length\n" ; 
                        }
                        
                        //If quality trimming flag is set up -> perform the quality trimming before vector/contaminants/adaptors clipping.
                        if( qual_trim_flag  ) 
                        {
                                QualTrimIllumina( read, max_a_error, max_e_at_ends );//This function generates LUCY clips of the read. Later they should be compared and read should be trimmed based on the results of comparison.
                                if (read->discarded_by_quality == 1) read->discarded = 1;
                        }
          
                        //Serial realization - useful for debugging if something does not work as expected
                        if(read->discarded == 0)
                        {
                                PolyAT_Trim(read);
            
                                if(read->poly_T_clip > 0) 
                                {
                                        left_trimmed_by_polyat += 1;
                                        read->lclip = read->lucy_lclip + read->poly_T_clip;
                                }
                                else
                                {
                                        if(read->lucy_lclip > 1 ) 
                                        {
                                                left_trimmed_by_quality += 1;
                                                read->lclip = read->lucy_lclip;
                                        }
                                }
            
                                if(read->poly_A_clip > 0) 
                                {
                                        right_trimmed_by_polyat += 1;
                                        read->rclip = read->lucy_rclip - read->poly_A_clip;
                                }
                                else
                                {
                                        if(read->lucy_rclip > 1 ) 
                                        {
                                                right_trimmed_by_quality += 1;
                                                read->rclip = read->lucy_rclip;
                                        }
                                }
            
                                if ( (read->rclip - read->lclip < minimum_read_length) && (read->rclip > 1) && (read->lclip > 1) ) 
                                {
                                        read->lucy_rclip = read->lucy_lclip = 0;
                                        read->discarded = 1;
                                        read->discarded_by_polyAT = 1;
                                        read->discarded_by_read_length = 1;
                                }
                        } 
                        
                        cnt+=1;
                        
                        rep_file << read->illumina_readID.substr(1,read->illumina_readID.length()-1) << "\t" << read->lclip << "\t" << read->rclip << "\t" << read->tru_sec_pos << "\t" << read->b_adapter << "\t" << read->initial_length << "\t" << (read->lucy_lclip <= 1 ? 1 : read->lucy_lclip) << "\t" << (read->lucy_rclip <= 1 ? 1 : read->lucy_rclip) << "\t" << read->discarded << "\t" << read->contaminants << "\t" << "NA" << "\n";
                      
                        if( read->lclip >= read->rclip ) { read->discarded = 1; read->discarded_by_read_length = 1; } 
                        if( read->lclip >= (int)read->read.length() ) { read->discarded = 1; read->discarded_by_read_length = 1; }
                        if( read->rclip > (int)read->read.length() ) { read->rclip = read->read.length(); }
                        if( (int)read->read.length() < minimum_read_length ) { read->discarded = 1; read->discarded_by_read_length = 1; }
                        if( (read->rclip - read->lclip) < minimum_read_length ) { read->discarded = 1; read->discarded_by_read_length = 1; }
              
                        if( read->discarded == 0 )
                        {
                            if(  read->rclip < read->initial_length  )
                            {
                                 cnt_right_trim_se += 1;
                                 //////////////avg_right_trim_len_se = GetAvg( avg_right_trim_len_se, cnt_right_trim_se, read->initial_length - read->rclip );
                            }
                            if(read->lclip > 0)
                            {
                                 cnt_left_trim_se += 1;
                                 /////////////avg_left_trim_len_se = GetAvg( avg_left_trim_len_se, cnt_left_trim_se, read->lclip );
                            }
                                    
                            read->read = read->read.substr(0 , read->rclip );
                            read->illumina_quality_string = read->illumina_quality_string.substr(0,read->rclip) ; 
                            read->read = read->read.substr( read->lclip, read->rclip - read->lclip );
                            read->illumina_quality_string = read->illumina_quality_string.substr( read->lclip, read->rclip - read->lclip );
                 
                            WriteSEFile(se_output_file, read);
                            se_accept_cnt+=1;
                            se_bases_kept += read->read.length();
                                    
                            if( read->initial_length > (int)read->read.length() )
                            {
                               cnt_avg+=1;
                               //////////////avg_trim_len_se = GetAvg( avg_trim_len_se, cnt_avg, /*read->initial_length - read->read.length()*/read->rclip - read->lclip );
                            }
                 
                            cnt_avg_len+=1; 
                            //////////////avg_len_se = GetAvg( avg_len_se, cnt_avg_len, read->read.length() );
                                    
                        } 
                                
                        if (read->discarded == 0) accepted++;
                        if (read->discarded == 1) discarded++;
                        if (read->left_trimmed_by_quality == 1) left_trimmed_by_quality++;
                        if (read->right_trimmed_by_quality == 1) right_trimmed_by_quality++;
                        
                        record_block.clear();
                        read->illumina_readID.clear(); 
                        read->illumina_quality_string.clear();
                        read->read.clear();
          
                        delete read;
                        
                        if( (cnt % 1000 ) == 0)
                        {
                            st_str = PrintIlluminaStatisticsPolyATSE(cnt, 
                                    se_bases_anal, 
                                    left_trimmed_by_quality,
                                    avg_left_trim_len_se,
                                    right_trimmed_by_quality,
                                    avg_right_trim_len_se,
                                    discarded_by_read_length,
                                    se_accept_cnt, se_bases_kept, 
                                    avg_len_se,
                                    left_trimmed_by_polyat,
                                    right_trimmed_by_polyat
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
    
    
    st_str = PrintIlluminaStatisticsPolyATSE(cnt, 
                                    se_bases_anal, 
                                    left_trimmed_by_quality,
                                    avg_left_trim_len_se,
                                    right_trimmed_by_quality,
                                    avg_right_trim_len_se,
                                    discarded_by_read_length,
                                    se_accept_cnt, se_bases_kept, 
                                    avg_len_se,
                                    left_trimmed_by_polyat,
                                    right_trimmed_by_polyat
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
    
    
    sum_stat_tsv << PrintIlluminaStatisticsTSVPolyATSE(cnt,
                                    se_bases_anal, 
                                    left_trimmed_by_quality, 
                                    avg_left_trim_len_se, 
                                    right_trimmed_by_quality,
                                    avg_right_trim_len_se,
                                    discarded_by_read_length,
                                    se_accept_cnt, 
                                    avg_len_se,
                                    left_trimmed_by_polyat,
                                    right_trimmed_by_polyat
                            ) << endl;
                 
    
    cout << "====================Done cleaning====================\n";  
    sum_stat << "====================Done cleaning====================\n";  
    
    se_output_file.close();
    
    rep_file.close();
    
    
}

/*
void RocheRoutineDynamic()
{
    if( !trim_adapters_flag ) 
    {
    }
    else
    {
       //Building dictionary for RL MIDS : 
       if(custom_rlmids_flag ) 
       {
          Build_RLMIDS_Dictionary(rlmids_file);
       } 
       else 
       {
          Build_RLMIDS_Dictionary();
       }
       
       //For each file in the given data set (right now for only one)
       for(int i=0; i<(int)roche_names.size(); ++i)
       {
            cout << "Parsing file: " << roche_names[i] << "..." << endl;
        
            //If SFF format is given -> process it
            if( string(roche_names[i]).substr( strlen(roche_names[i])-3, 3 ) == "sff" ) 
            {
               cout << "File is in SFF format, starting conversion...\n" ;
               if(output_fastqfile_flag == false)
               {
                 output_sfffile_flag = true;
               }
               
               sff_common_header h;
               sff_read_header rh;
               sff_read_data rd;
               FILE *sff_fp;

               if ( (sff_fp = fopen(roche_names[i], "r")) == NULL ) 
               {
                    fprintf(stderr,
                                "[err] Could not open file '%s' for reading.\n", roche_names[i]);
                    exit(1);
               }
               
               read_sff_common_header(sff_fp, &h);
               //verify_sff_common_header((char*)PRG_NAME, (char*)SFF_FILE_VERSION, &h);

               int left_clip = 0, right_clip = 0, nbases = 0;
               char *bases;
               uint8_t *qualily;
                
               unsigned int numreads = (unsigned int) h.nreads;
               cout << "Conversion finished. Total number of reads read from given file(s): " << numreads << endl;
       
               for (unsigned int i = 0; i < numreads; i++) 
               { 
                    read_sff_read_header(sff_fp, &rh);
                    read_sff_read_data(sff_fp, &rd, h.flow_len, rh.nbases);

                    get_clip_values(rh, 0, &left_clip, &right_clip);
                    nbases = right_clip - left_clip;

                    // create bases string 
                    bases = get_read_bases(rd, left_clip, right_clip);

                    // create quality array 
                    qualily = get_read_quality_values(rd, left_clip, right_clip);

                    for (int j = 0; j < nbases; j++) 
                    {
                        qualily[j] = (qualily[j] <= 93 ? qualily[j] : 93) + 33;
                    }
                    
                    
                    
                    free(bases);
                    free(qualily);
                    free_sff_read_header(&rh);
                    free_sff_read_data(&rd);


                }

                fclose(sff_fp);
        
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
    }
}
*/

void PrintHelp() {
    cout << "Version: " << version << endl;
    cout << "**********************************************************************************************************************\n";        
    cout << "Program SeqyClean\n" 
            "Main purpose of this software is to clean NGS reads. It provide adapter/key/primers searching and quality trimming (LUCY).\n" 
            "Usage:\n"  
            "Roche 454:\n"
            "./seqyclean -454 input_file_name -o output_prefix\n"
							"[-v vector_file]\n"
							"[-c file_of_contaminants]\n"
							"[-m file_of_RL_MIDS]\n" 
							"[-k k_mer_size]\n"
							"[-kc k_mer_size]\n"
							"[-f overlap ]\n"
							"[-t number_of_threads]\n" 
							"[-qual max_avg_error max_error_at_ends]\n"
							"[--qual_only]\n"
							"[--fastq]\n"
                                                        "[--ow]\n"
							"[--keep_fastq_orig]\n"
							"[-minimum_read_length <value>]\n"
							"[-polyat [cdna] [cerr] [crng] ]\n"
            "For Illumina paired-end reads:\n"
            "./seqyclean -1 input_file_name_1 -2 input_file_name_2 -o output_prefix\n"
							"[-v vector_file]\n"
							"[-c file_of_contaminants]\n"
							"[-k k_mer_size]\n"
							"[-kc k_mer_size]\n" 
							"[-qual max_avg_error max_error_at_ends]\n"
							"[--qual_only]\n"
							"[-minimum_read_length <value>]\n"
                                                        "[--shuffle]\n"
                                                        "[-i64]\n"
                                                        "[-adapter_length <value>]\n"
                                                        "[-ot <value>]\n"
                                                        "[--overlap <minoverlap=value>]\n"
                                                        "[--ow]\n"
                                                        "[--dup]\n" 
                                                        "[-polyat [cdna] [cerr] [crng] ]\n"
                                                        "[--new2old_illumina] - switch to fix read IDs ( As is detailed in: http://contig.wordpress.com/2011/09/01/newbler-input-iii-a-quick-fix-for-the-new-illumina-fastq-header/#more-342 )\n"
            "For Illumina single-end reads:\n"
            "./seqyclean -U input_file_name -o output_prefix\n"
							"[-v vector_file]\n"
							"[-c file_of_contaminants]\n"
							"[-k k_mer_size]\n"
							"[-kc k_mer_size]\n" 
							"[-qual max_avg_error max_error_at_ends]\n"
							"[--qual_only]\n"
							"[-minimum_read_length <value>]\n"
                                                        "[--shuffle]\n"
                                                        "[--ow]\n"
                                                        "[-polyat [cdna] [cerr] [crng] ]\n"
                                                        "[--new2old_illumina] - switch to fix read IDs ( As is detailed in: http://contig.wordpress.com/2011/09/01/newbler-input-iii-a-quick-fix-for-the-new-illumina-fastq-header/#more-342 )\n";
cout <<"Example:\n"
"Roche:\n"
"./seqyclean -454 test_data/in.sff -o test/Test454 -v test_data/vectors.fasta\n"
"Illumina:\n"
"./seqyclean -1 test_data/R1.fastq.gz -2 test_data/R2.fastq.gz -o test/Test_Illumina\n";
    
    
    cout << "Please ask Ilya by email: zhba3458@vandals.uidaho.edu in case of any questions.\n" ;
 
}

void PolyAT_Trim(Read* read)
{
    int left, right;
    left = right = 0;

    left = poly_at_left( (char*)read->read.substr( read->lucy_lclip, read->read.length() - read->lucy_lclip ).c_str(), read->lucy_rclip - read->lucy_lclip + 1); 
    
    if (left) 
    {
        read->poly_T_clip = left;
    }
	
    right = poly_at_right((char*)read->read.substr( 0, read->lucy_rclip).c_str(), read->lucy_rclip - read->lucy_lclip + 1);
    
    if (right) 
    {
        read->poly_A_clip = right;
    }
    
    
}






string PrintIlluminaStatisticsPolyAT(unsigned long cnt1, unsigned long cnt2, 
                                    unsigned long long  pe1_bases_anal, unsigned long long  pe2_bases_anal, 
                                    unsigned long left_trimmed_by_quality1, unsigned long left_trimmed_by_quality2,
                                    double avg_left_trim_len_pe1, double avg_left_trim_len_pe2, 
                                    unsigned long right_trimmed_by_quality1,unsigned long right_trimmed_by_quality2,
                                    double avg_right_trim_len_pe1,double avg_right_trim_len_pe2,
                                    unsigned long discarded_by_read_length1, unsigned long discarded_by_read_length2,
                                    unsigned long pe_accept_cnt, unsigned long long  pe_bases_kept, 
                                    unsigned long pe_discard_cnt,unsigned long long  pe_bases_discarded, 
                                    unsigned long se_pe1_accept_cnt, unsigned long long se_pe1_bases_kept,
                                    unsigned long se_pe2_accept_cnt, unsigned long long se_pe2_bases_kept,
                                    double avg_trim_len_pe1, double avg_trim_len_pe2,
                                    unsigned long left_trimmed_by_polyat1, unsigned long left_trimmed_by_polyat2,
                                    unsigned long right_trimmed_by_polyat1, unsigned long right_trimmed_by_polyat2
                                    )
{
    
     string stat_str = string("====================Summary Statistics====================\n");
     stat_str += string("PE1 reads analyzed: ") +  int2str(cnt1)   + string(", Bases:") +  int2str(pe1_bases_anal) + string("\n")  +
                        ( (qual_trim_flag || polyat_flag || vector_flag) ? "Reads left trimmed ->\n" : "" ) +
                        ( qual_trim_flag ? "By quality: " +  i2str(left_trimmed_by_quality1,new char[15],10) + "\n" : "" ) +
                        "By poly A/T: " + i2str(left_trimmed_by_polyat1,new char[15],10) + "\n" +
                        "Average left trim length: " + double2str(avg_left_trim_len_pe1) + " bp\n" +
                        "Reads right trimmed ->\n" +
                        ( qual_trim_flag ? "By quality: " +  i2str(right_trimmed_by_quality1,new char[15],10) + "\n" : "") +
                        "By poly A/T: " + i2str(right_trimmed_by_polyat1,new char[15],10) + "\n" +
                        "Average right trim length: " + double2str(avg_right_trim_len_pe1) + " bp\n" +
                        "PE1 reads discarded by read length: " +  i2str(discarded_by_read_length1,new char[15],10) + "\n" +
                        "-----------------------------------------------------------\n" +
                        "PE2 reads analyzed: " + i2str(cnt2,new char[15],10) + ", Bases:" + i2str(pe2_bases_anal,new char[15],10) + "\n" +
                        ( (qual_trim_flag || polyat_flag || vector_flag ) ? "Reads left trimmed ->\n" : "" ) +
                        (qual_trim_flag ? "By quality: " +  i2str(left_trimmed_by_quality2,new char[15],10) + "\n" : "" ) +
                        "By poly A/T: " + i2str(left_trimmed_by_polyat2,new char[15],10) + "\n" +
                        ("Average left trim length: " + double2str(avg_left_trim_len_pe2) + " bp\n" ) +
                        "Reads right trimmed ->\n" +
                        ( qual_trim_flag ? "By quality: " +  i2str(right_trimmed_by_quality2,new char[15],10) + "\n" : "") +
                        "By poly A/T: " + i2str(right_trimmed_by_polyat2,new char[15],10) + "\n" +
                        ("Average right trim length: " + double2str(avg_right_trim_len_pe2) + " bp\n") +
                        "PE2 reads discarded by read length: " +  i2str(discarded_by_read_length2,new char[15],10) + "\n" + 
                        "----------------------Summary for PE & SE----------------------\n" +
                        ("Pairs kept: " + i2str(pe_accept_cnt,new char[15],10) + ", " + double2str( (double)pe_accept_cnt/(double)cnt1*100.0) + "%, Bases: " + i2str(pe_bases_kept,new char[15],10) + ", " + double2str( (double)pe_bases_kept/(double)(pe1_bases_anal+pe2_bases_anal)*100) +  "%\n") +
                        ("Pairs discarded: " + i2str(pe_discard_cnt,new char[15],10) + ", " + double2str( (double)pe_discard_cnt/(double)cnt1*100.0) + "%, Bases: " + i2str(pe_bases_discarded,new char[15],10) + ", " + double2str( (double)pe_bases_discarded/(double)(pe1_bases_anal+pe2_bases_anal)*100) +  "%\n") +
                        ("Single Reads PE1 kept: " + i2str(se_pe1_accept_cnt,new char[15],10) + ", Bases: " + i2str(se_pe1_bases_kept,new char[15],10) +"\n") +
                        ("Single Reads PE2 kept: " + i2str(se_pe2_accept_cnt,new char[15],10) + ", Bases: " + i2str(se_pe2_bases_kept,new char[15],10) +"\n") +
                        ("Average trimmed length PE1: " + double2str(avg_trim_len_pe1) + " bp\n") +
                        ("Average trimmed length PE2: " + double2str(avg_trim_len_pe2) + " bp\n");
                        
     return stat_str;
    
}

string PrintIlluminaStatisticsTSVPolyAT(unsigned long cnt1, unsigned long cnt2, 
                                    unsigned long long  pe1_bases_anal, unsigned long long  pe2_bases_anal, 
                                    unsigned long left_trimmed_by_quality1, unsigned long left_trimmed_by_quality2,
                                    double avg_left_trim_len_pe1, double avg_left_trim_len_pe2, 
                                    unsigned long right_trimmed_by_quality1,unsigned long right_trimmed_by_quality2,
                                    double avg_right_trim_len_pe1,double avg_right_trim_len_pe2,
                                    unsigned long discarded_by_read_length1, unsigned long discarded_by_read_length2,
                                    unsigned long pe_accept_cnt, unsigned long long pe_bases_kept, 
                                    unsigned long pe_discard_cnt,unsigned long long pe_bases_discarded, 
                                    unsigned long se_pe1_accept_cnt, unsigned long long se_pe1_bases_kept,
                                    unsigned long se_pe2_accept_cnt, unsigned long long se_pe2_bases_kept,
                                    double avg_trim_len_pe1, double avg_trim_len_pe2,
                                    unsigned long left_trimmed_by_polyat1, unsigned long left_trimmed_by_polyat2,
                                    unsigned long right_trimmed_by_polyat1, unsigned long right_trimmed_by_polyat2
                                    )
{
    
        string filename_str;
    
        for(int i=0; i<(int)pe1_names.size(); ++i)
        {
            filename_str += string(pe1_names[i]) + ", " + string(pe2_names[i]);
        }
        
        string stat_str_tsv =   version + "\t" + 
                                filename_str + "\t" +
                                i2str(cdna,new char[15],10) + "\t" +
                                i2str(c_err,new char[15],10) + "\t" +
                                i2str(crng,new char[15],10) + "\t" +
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
                                int2str(minimum_read_length)+ "\t" +
                                ( new2old_illumina ? "YES" : "NO") + "\t"; 
                   
    
    
        stat_str_tsv += int2str(cnt1)   + "\t" +  int2str(pe1_bases_anal) + "\t"  +
                       i2str(left_trimmed_by_polyat1,new char[15],10) + "\t"  +
                       ( qual_trim_flag ? i2str(left_trimmed_by_quality1,new char[15],10) + "\t" : "NA\t" ) +
                       double2str(avg_left_trim_len_pe1) + "\t" +
                       i2str(right_trimmed_by_polyat1,new char[15],10) + "\t"  +
                       ( qual_trim_flag ? i2str(right_trimmed_by_quality1,new char[15],10) + "\t" : "NA\t") +
                       double2str(avg_right_trim_len_pe1) + "\t" +
                       i2str(discarded_by_read_length1,new char[15],10) + "\t" +
                       i2str(cnt2,new char[15],10) + "\t" + i2str(pe2_bases_anal,new char[15],10) + "\t" +
                       i2str(left_trimmed_by_polyat2,new char[15],10) + "\t"  +
                       (qual_trim_flag ? i2str(left_trimmed_by_quality2,new char[15],10) + "\t" : "NA\t" ) +
                       double2str(avg_left_trim_len_pe2) + "\t"  +
                       i2str(right_trimmed_by_polyat2,new char[15],10) + "\t"  +
                       ( qual_trim_flag ? i2str(right_trimmed_by_quality2,new char[15],10) + "\t" : "NA\t") +
                       double2str(avg_right_trim_len_pe2) + "\t" +
                       i2str(discarded_by_read_length2,new char[15],10) + "\t" + 
                       (i2str(pe_accept_cnt,new char[15],10) + "\t" + double2str( (double)pe_accept_cnt/(double)cnt1*100.0) + "\t" + i2str(pe_bases_kept,new char[15],10) + "\t" + double2str( (double)pe_bases_kept/(double)(pe1_bases_anal+pe2_bases_anal)*100) +  "\t") +
                       ( i2str(pe_discard_cnt,new char[15],10) + "\t" + double2str( (double)pe_discard_cnt/(double)cnt1*100.0) + "\t" + i2str(pe_bases_discarded,new char[15],10) + "\t" + double2str( (double)pe_bases_discarded/(double)(pe1_bases_anal+pe2_bases_anal)*100) +  "\t") +
                       (i2str(se_pe1_accept_cnt,new char[15],10) + "\t" + i2str(se_pe1_bases_kept,new char[15],10) +"\t") +
                       (i2str(se_pe2_accept_cnt,new char[15],10) + "\t" + i2str(se_pe2_bases_kept,new char[15],10) +"\t") +
                       (double2str(avg_trim_len_pe1) + "\t") +
                       (double2str(avg_trim_len_pe2) + "\t");
                      
     return stat_str_tsv;
    
    
}

string PrintIlluminaStatisticsPolyATSE( unsigned long cnt, 
                                        unsigned long long se_bases_anal, 
                                        unsigned long left_trimmed_by_quality,
                                        double avg_left_trim_len_se,
                                        unsigned long right_trimmed_by_quality,
                                        double avg_right_trim_len_se,
                                        unsigned long discarded_by_read_length,
                                        unsigned long se_accept_cnt, 
                                        unsigned long long se_bases_kept, 
                                        double avg_trim_len_se,
                                        unsigned long left_trimmed_by_polyat,
                                        unsigned long right_trimmed_by_polyat
                                    )
{
    
   
    
    string ans = "====================Summary Statistics====================\n" +
                        ("SE reads analyzed: " +  i2str(cnt,new char[15],10)  + ", Bases:" +  i2str(se_bases_anal, new char[25],10)  + "\n") +
                        "Reads left trimmed ->\n" +
                        ( qual_trim_flag ? "By quality: " +  i2str(left_trimmed_by_quality,new char[15],10) + "\n" : "" ) +
                        "By poly A/T: " +  i2str(left_trimmed_by_polyat,new char[15],10) + "\n" +
                        "Average left trim length: " + double2str(avg_left_trim_len_se) + " bp\n" +
                        "Reads right trimmed ->\n" +
                        ( qual_trim_flag ? "By quality: " +  i2str(right_trimmed_by_quality,new char[15],10) + "\n" : "") +
                        "By poly A/T: " +  i2str(right_trimmed_by_polyat,new char[15],10) + "\n" +
                        "Average right trim length: " + double2str(avg_right_trim_len_se) + " bp\n" +
                        "Discarded by read length: " +  i2str(discarded_by_read_length,new char[15],10) + "\n" +
                        "----------------------Summary for SE----------------------\n" +
                        ("Reads kept: " + i2str(se_accept_cnt,new char[15],10) + ", " + double2str( (double)se_accept_cnt/(double)cnt*100.0) + "%, Bases: " + i2str(se_bases_kept,new char[15],10) + ", " + double2str( (double)se_bases_kept/(double)(se_bases_anal)*100) +  "%\n") +
                        ("Average trimmed length: " + double2str(avg_trim_len_se) + " bp\n");// +
                       
    
    return ans;
   
}

string PrintIlluminaStatisticsTSVPolyATSE(unsigned long cnt,
                                    unsigned long long se_bases_anal, 
                                    unsigned long left_trimmed_by_quality, 
                                    double avg_left_trim_len_se, 
                                    unsigned long right_trimmed_by_quality,
                                    double avg_right_trim_len_se,
                                    unsigned long discarded_by_read_length,
                                    unsigned long se_accept_cnt, 
                                    double avg_trim_len_se,
                                    unsigned long left_trimmed_by_polyat,
                                    unsigned long right_trimmed_by_polyat
                                   )
{
    
        string filename_str;
    
        for(int i=0; i<(int)pe1_names.size(); ++i)
        {
            filename_str += string(se_names[i]);
        }
        
        string stat_str_tsv =   version + "\t" + 
                                filename_str + "\t" +
                                i2str(cdna,new char[15],10) + "\t" +
                                i2str(c_err,new char[15],10) + "\t" +
                                i2str(crng,new char[15],10) + "\t" +
                                (qual_trim_flag ? "YES" : "NO") +"\t" +
                                (qual_trim_flag ? double2str(-10*log10(max_a_error)) : "NA")+ "\t" +
                                (qual_trim_flag ? double2str(-10*log10(max_e_at_ends)) : "NA")+ "\t" +
                                output_prefix +"\t" +
                                rep_file_name1 + "\t" +
                                se_output_filename +"\t"+
                                int2str(minimum_read_length)+ "\t" +
                                ( new2old_illumina ? "YES" : "NO") + "\t"; 
                   
    
    
        stat_str_tsv += int2str(cnt)   + "\t" + //reads analyzed
                        int2str(se_bases_anal) + "\t"  + //bases
                        i2str(left_trimmed_by_polyat,new char[15],10) + "\t" +
                        ( qual_trim_flag ? i2str(left_trimmed_by_quality,new char[15],10) + "\t" : "NA\t" ) +  //left trimmed qual
                        double2str(avg_left_trim_len_se) + "\t" + //avg left trim len
                        i2str(right_trimmed_by_polyat,new char[15],10) + "\t" +
                        ( qual_trim_flag ? i2str(right_trimmed_by_quality,new char[15],10) + "\t" : "NA\t") +
                        double2str(avg_right_trim_len_se) + "\t" +
                        i2str(discarded_by_read_length,new char[15],10) + "\t" +
                        i2str(se_accept_cnt,new char[15],10) + "\t" + //se reads kept
                        double2str( (double)se_accept_cnt/(double)cnt*100.0) + "\t" //perc kept
                        + i2str(se_bases_kept,new char[15],10) + "\t" + //bases kept
                        double2str( (double)se_bases_kept/(double)se_bases_anal*100.0) + "\t" + //%
                        double2str(avg_trim_len_se);
            
    
     return stat_str_tsv;
    
    
}
