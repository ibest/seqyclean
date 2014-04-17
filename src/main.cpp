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
#include "gzstream.h"
#include "Roche.h"
#include "rlmid.h"
#include "Illumina.h"


using namespace std;

/*Computational parameters (default)*/
short KMER_SIZE = 15;
short DISTANCE = 1;
unsigned short NUM_THREADS = 4;

string version = "1.9.1 (2014-04-13)";

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
// Starting position and a window size for searching for duplicates:
extern int start_dw = 10, size_dw = 15;

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


/*Maximum number of mismatches allowed in alignment operation*/
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
void ParseFastqFileIllumina(char* fastq_file, vector<Read*> &reads );
void ClearNNs( vector<Read*>& reads );
void IlluminaRoutine();
void RemoveContaminants(vector<Read*>& illumina_reads);

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

unsigned int adapterlength = 60;
double overlap_t = 0.75;
int minoverlap = 16;
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
        else if( string(argv[i]) == "-qual" ) 
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
        else if( string(argv[i]) == "--qual_only" ) 
        {
           lucy_only_flag = true;
           trim_adapters_flag = false;
           continue;
        } 
        else if( string(argv[i]) == "--shuffle" ) 
        {
           shuffle_flag = true;
           continue;
        }
        else if( string(argv[i]) == "--dup" ) 
        {
           rem_dup = true;
           if( string(argv[i]) == "-start_dw" ) 
           {
               start_dw = atoi(argv[++i]);
           }
           if( string(argv[i]) == "-size_dw" ) 
           {
               size_dw = atoi(argv[++i]);
           }
           
           continue;
        }
        else if( string(argv[i]) == "--ow" ) 
        {
           overwrite_flag = true;
           continue;
        }
        else if( string(argv[i]) == "--overlap" ) 
        {
           overlap_flag = true;
           if ( ( (i+1)<argc ) && (argv[i+1][0] != '-') ) 
           {
              minoverlap = atoi(argv[++i]);
           }
           continue;
        }
        else if( string(argv[i]) == "-max_al_mism" ) 
        {
           max_al_mism = atoi(argv[i]);
           continue;
        }
        else if( string(argv[i]) == "--keep_fastq_orig" ) 
        {
           keep_fastq_orig = true;
           continue;
        }
        else if( string(argv[i]) == "--no_ts_adapter_trim" ) 
        {
           trim_adapters_flag = false;
           continue;
        }
        else if( string(argv[i]) == "-polyat" ) 
        {
           polyat_flag = true;
           
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
        else if(string(argv[i]) == "-v" )
        {
           if ((i+1)<argc && !(isdigit(argv[i+1][0])) ) 
           {
             vector_file = argv[++i]; /*Vector file given*/
             vector_flag = true;
           }
           continue;
        }
        else if(string(argv[i]) == "-c" )
        { 
           if ((i+1)<argc && !(isdigit(argv[i+1][0])) ) 
           {
             cont_file = argv[++i]; /*File with contaminants given*/
             contaminants_flag = true;
           }
           continue;
        }
        else if(string(argv[i]) == "-m" )
        {
           if ( (i+1)<argc && (isdigit(argv[i+1][0])) ) 
           {
             custom_rlmids_flag = true;
             rlmids_file = argv[++i]; /*Custom file with RL MIDS given*/
           }
           continue;
        }
        else if(string(argv[i]) == "-adapter_length" )
        {
           if ( (i+1)<argc && (isdigit(argv[i+1][0])) ) 
           {
             adapterlength = atoi(argv[++i]); 
           }
           continue;
        }
        else if(string(argv[i]) == "-ot" )
        {
           if ( (i+1)<argc && (isdigit(argv[i+1][0])) ) 
           {
             overlap_t = atof(argv[++i]);
           }
           continue;
        }
        else if(string(argv[i]) == "-k" )
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
        else if(string(argv[i]) == "-f" )
        {
           if ( (i+1)<argc && isdigit(argv[i+1][0]) ) 
           {
              DISTANCE = atoi(argv[++i]);
           }
           continue;
        }
        else if(string(argv[i]) == "--test" )
        {
           if ( (i+1)<argc && isdigit(argv[i+1][0]) ) 
           {
              DISTANCE = atoi(argv[++i]);
           }
           continue;
        }
        else if(string(argv[i]) == "-t" )
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
        else if(string(argv[i]) == "-p" )
        {
           if ( (i+1)<argc && !(isdigit(argv[i+1][0])) ) 
           {
              pcr_flag = true;
              pcr_file_name = argv[++i]; 
           }
           continue;
        }
        else if(string(argv[i]) == "-i64" )
        {
           i64_flag = true;
           
           continue;
        }
        else if(string(argv[i]) == "--serial" )
        {
           serial_flag = true;
           continue;
        }
        else if(string(argv[i]) == "-d" )
        {
           debug_flag = true;
           continue;
        }
        else if(string(argv[i]) == "--new2old_illumina" )
        {
           new2old_illumina = true;
           continue;
        }
        else if(string(argv[i]) == "-L_limit" )
        {
           if ( isdigit(argv[i+1][0]) ) 
                L_limit = atoi(argv[++i]);
           
           continue;
        }
        else if(string(argv[i]) == "-R_limit" )
        {
           if ( isdigit(argv[i+1][0]) ) 
                R_limit = atoi(argv[++i]);
           
           continue;
        }
        else if(string(argv[i]) == "-vmr" )
        {
           if ( isdigit(argv[i+1][0]) ) 
                vmr = atoi(argv[++i]);
           
           continue;
        }
        else if(string(argv[i]) == "-vml" )
        {
           if ( isdigit(argv[i+1][0]) ) 
                vml = atoi(argv[++i]);
           
           continue;
        }
        else if(string(argv[i]) == "-allowable_distance" )
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
        else if(string(argv[i]) == "-minimum_read_length" )
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
	else if(string(argv[i]) == "-kc" )
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
        else if(string(argv[i]) == "--RLMIDS" )
        {
           cout << "Supported RL MIDS:\n";
           for(int i=0; i<36; i++) 
           {
              cout << mids[i] << endl;
           }
           
           return 0;
        }
        else if(string(argv[i]) == "-o" ) 
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
        else if( string(argv[i]) == "--fastq" ) 
        {
           output_fastqfile_flag = true; //output file with cleaned reads in FASTQ format, for 454 mode only
           continue;
        }
        else if(string(argv[i]) == "-?" )
        {
           PrintHelp();
           exit(1);
        }   
        else if(string(argv[i]) == "-h" )
        {
           PrintHelp();
           exit(1);
        }
        else if( string(argv[i]) == "-1" )
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
        else if( string(argv[i]) == "-2" )
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
        else if( string(argv[i]) == "-U" ) //single-end file mode
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
        else if( string(argv[i]) == "-454" ) 
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
        else if( string(argv[i]) == "--MakeRocheReport" ) 
        {
            MakeRocheReport((char*)"RocheClipPoints_report.tsv", argv[++i]);
            return 0;
        }
        else if( string(argv[i]) == "-pmax" ) 
        {
            pmax = atoi(argv[++i]);
            continue;
        }
        else if( string(argv[i]) == "--vonly" ) 
        {
            VectorOnlyFlag = true;
            continue;    
        } else {
            cout << "Unknown parameter: " << argv[i] << endl;
            PrintHelp();
            return 0;
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
                /*if (fields.size() == 5)
                {
                    //Old headers == True
                    old_style_illumina_flag = true;
                    
                } else //if (fields.size() == 7)
                {
                    //Old headers == False
                    old_style_illumina_flag = false;
                }*/
                
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
    
    if ( (!illumina_flag) && (!roche_flag) )
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
    if(illumina_flag)
    {
        cout << "--------------------Basic parameters--------------------\n";
        sum_stat << "--------------------Basic parameters--------------------\n";
        
        //Sum stat TSV header
        
        
        if(!illumina_se_flag)
        {
        
            sum_stat_tsv << "Version\tPE1PE2\tAdapters_trim\tVectorTrim\tK_mer_size\tDistance\tContamScr\tkc_size\tQualTrim\tQual_max_avg_error\tQual_max_err_at_ends\tOutputPrefix\tRepFilename1\tRepFilename2\tPE1OutputFilename\tPE2OutputFilename\tShuffledFilename\tSEFilename\tMax_align_mismatches\tMinReadLen\tnew2old_illumina\tPE1ReadsAn\tPE1Bases\tPE1TruSeqAdap_found\tPerc_PE1TruSeq\tPE1ReadsWVector_found\tPerc_PE1ReadsWVector\tPE1ReadsWContam_found\tPerc_PE1ReadsWContam\tPE1LeftTrimmedQual\tPE1LeftTrimmedVector\tPE1AvgLeftTrimLen\tPE1RightTrimmedAdap\tPE1RightTrimmedQual\tPE1RightTrimmedVector\tPE1AvgRightTrimLen\tPE1DiscardedTotal\tPE1DiscByContam\tPE1DiscByLength\tPE2ReadsAn\tPE2Bases\tPE2TruSeqAdap_found\tPerc_PE2TruSeq\tPE2ReadsWVector_found\tPerc_PE2ReadsWVector\tPE2ReadsWContam_found\tPerc_PE2ReadsWContam\tPE2LeftTrimmedQual\tPE2LeftTrimmedVector\tPE2AvgLeftTrimLen\tPE2RightTrimmedAdap\tPE2RightTrimmedQual\tPE2RightTrimmedVector\tPE2AvgRightTrimLen\tPE2DiscardedTotal\tPE2DiscByContam\tPE2DiscByLength\tPairsKept\tPerc_Kept\tBases\tPerc_Bases\tPairsDiscarded\tPerc_Discarded\tBases\tPerc_Bases\tSE_PE1_Kept\tSE_PE1_Bases\tSE_PE2_Kept\tSE_PE2_Bases\tAvgTrimmedLenPE1\tAvgTrimmedLenPE2\tDuplicates\tperfect_ov\tpartial_ov\tPolyAT\tcdna\tc_err\tcrng\tleft_trimmed_by_polyat1\tright_trimmed_by_polyat1\tleft_trimmed_by_polyat2\tright_trimmed_by_polyat2\n";
    
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
                if(polyat_flag) {
                        cout << "Poly A/T trimming: YES cdna =" << cdna << " c_err = " << c_err << " crnq = " << crng << endl;
                        sum_stat << "Poly A/T trimming: YES cdna =" << cdna << " c_err = " << c_err << " crnq = " << crng << endl;
                } else {
                        cout << "Poly A/T trimming: NO" << endl;
                        sum_stat << "Poly A/T trimming: NO" << endl;
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
            sum_stat_tsv << "Version\tSE\tAdapters_trim\tVectorTrim\tK_mer_size\tDistance\tContamScr\tkc_size\tQualTrim\tQual_max_avg_error\tQual_max_err_at_ends\tOutputPrefix\tRepFilename\ttSEOutputFilename\tMax_align_mismatches\tMinReadLen\tnew2old_illumina\tSEReadsAn\tSEBases\tSETruSeqAdap_found\tPerc_SETruSeq\tSEReadsWVector_found\tPerc_SEReadsWVector\tSEReadsWContam_found\tPerc_SEReadsWContam\tSELeftTrimmedQual\tSELeftTrimmedVector\tSEAvgLeftTrimLen\tSERightTrimmedAdap\tSERightTrimmedQual\tSERightTrimmedVector\tSEAvgRightTrimLen\tSEDiscardedTotal\tSEDiscByContam\tSEDiscByLength\tSEReadsKept\tPerc_Kept\tBases\tPerc_Bases\tAvgTrimmedLenSE\tPolyAT\tcdna\tc_err\tcrng\tleft_trimmed_by_polyat\tright_trimmed_by_polyat\n";
    
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
                
                if(polyat_flag) {
                        cout << "Poly A/T trimming: YES cdna =" << cdna << " c_err = " << c_err << " crnq = " << crng << endl;
                        sum_stat << "Poly A/T trimming: YES cdna =" << cdna << " c_err = " << c_err << " crnq = " << crng << endl;
                } else {
                        cout << "Poly A/T trimming: NO" << endl;
                        sum_stat << "Poly A/T trimming: NO" << endl;
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
    //if(roche_flag & !polyat_flag)
    if(roche_flag)
    {
        sum_stat_tsv << "Version\tFiles\tNUM_THREADS\tAdaptersTrimming\tVectorTrimming\tkmer_size\tDistance\tContamScr\tkmer_contam_size\tQualTrim\tQual_max_avg_error\tQual_max_err_at_ends\tOutputPrefix\tRepFilename\tOutputFilename\tMax_align_mismatches\tMinReadLen\tReadsAnalyzed\tBases\tleft_mid_tags_found\tpercentage_left_mid_tags_found\tright_mid_tags_found\tpercentage_right_mid_tags_found\tReadsWithVector_found\tpercentage_ReadsWithVector_found\tReadsWithContam_found\tpercentage_ReadsWithContam_found\tLeftTrimmedByAdapter\tLeftTrimmedByQual\tLeftTrimmedByVector\tAvgLeftTrimLen\tRightTrimmedByAdapter\tRightTrimmedByQual\tRightTrimmedByVector\tAvgRightTrimLen\tDiscardedTotal\tDiscByContam\tDiscByLength\tReadsKept\tPercentageKept\tAvgTrimmedLen\tPolyAT\tcdna\tc_err\tcrng\tleft_trimmed_by_polyat\tright_trimmed_by_polyat\n";
    
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
        
        if(polyat_flag) {
           cout << "Poly A/T trimming: YES cdna =" << cdna << " c_err = " << c_err << " crnq = " << crng << endl;
           sum_stat << "Poly A/T trimming: YES cdna =" << cdna << " c_err = " << c_err << " crnq = " << crng << endl;
 
        } else {
           cout << "Poly A/T trimming: NO" << endl;
           sum_stat << "Poly A/T trimming: NO" << endl;
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
    
    if(illumina_flag) 
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
    
    //if( roche_flag && !polyat_flag )
    if( roche_flag )
    {
        RocheRoutine();
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
                                                        "[--dup][-start_dw][-size_dw]\n" 
                                                        "[-polyat [cdna] [cerr] [crng] ]\n"
                                                        "[--no_ts_adapter_trim]\n"
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
                                                        "[--no_ts_adapter_trim]\n"
                                                        "[--new2old_illumina] - switch to fix read IDs ( As is detailed in: http://contig.wordpress.com/2011/09/01/newbler-input-iii-a-quick-fix-for-the-new-illumina-fastq-header/#more-342 )\n";
cout <<"Example:\n"
"Roche:\n"
"./seqyclean -454 test_data/in.sff -o test/Test454 -v test_data/vectors.fasta\n"
"Illumina:\n"
"./seqyclean -1 test_data/R1.fastq.gz -2 test_data/R2.fastq.gz -o test/Test_Illumina\n";
    
    
    cout << "Please ask Ilya by email: zhba3458@vandals.uidaho.edu in case of any questions.\n" ;
 
}
