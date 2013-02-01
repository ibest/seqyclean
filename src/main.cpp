#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <streambuf>
#include <exception>
#include <pthread.h>
#include <bits/basic_string.h>
#include "timer.h"
#include "util.h"
#include "poly.h"
#include "Read.h"
#include "Dictionary.h"
#include "KMerRoutine.h"
#include "Report.h"
#include "MainPipeLine.h"
#include "TrimLucy.h"
#include "sffreader.h"
#include "Illumina.h"
#include "rlmid.h"
#include "gzstream.h"
#include <sys/types.h>
#include <sys/stat.h>


using namespace std;

/*Computational parameters (default)*/
short KMER_SIZE = 15;
short DISTANCE = 1;
short NUM_THREADS = 4;

string version = "1.2.6"; 

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

/*Illumina data structures*/
vector<Read*> reads_1;
vector<Read*> reads_2;
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
char* illumina_file_name_R1;// = "";
char* illumina_file_name_R2;// = "";
char* illumina_file_name_se;

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

/*Poly A/T trimming default parameters*/
int cdna = 10;
int c_err = 3;
int crng=50;
int keep;

/*Vector trimming*/
int L_limit = 1;
int R_limit = 1;
int vmr = 0;//15;
int vml = 0;//15;
int allowable_distance = 3;
int minimum_read_length = 50;
int KMER_SIZE_CONT = 15;
int pmax = 2;

/*Other variables and parameters*/
std::ifstream read_file;

void PrintHelp();
void ParseFastqFile(char *fastq_file, vector<Read*> &reads);
//void ParseFastqFileIllumina(char* fastq_file, map<string, Read> &reads );
//void ParseFastqFileIllumina(char* fastq_file, dense_hash_map<string, Read> &reads );
void ParseFastqFileIllumina(char* fastq_file, vector<Read*> &reads );
void ClearNNs( vector<Read*>& reads );
void QualityTrimming( vector<Read*>& reads );
void IlluminaRoutine();
void RocheRoutine();
void RemoveContaminants(vector<Read*>& illumina_reads);
void RemoveContaminants454(vector<Read*>& reads454);
void PolyAT_Trim(Read* read);
void IlluminaDynamic(char* fastq_file1, char* fastq_file2);
int IlluminaDynRoutine(Read* read, bool& adapter_found, string& query_str);
void WritePEFile(fstream &pe_output_file, Read *read);
void WriteShuffleFile(fstream &shuffle_output_file, Read *read1, Read *read2);
void WriteSEFile(fstream &se_output_file, Read *read);
void *t_IlluminaDynRoutine( void *targs );
string New2OldNbl(string header);
void PolyATRoutine();
void IlluminaDynamicSE();
/*-------------------------------------*/

bool dynflag = false;

vector<string> file_list;

fstream sum_stat;

string output_prefix;

bool VectorOnlyFlag = false;
bool new2old_illumina = false;

struct thread_data {
   Read* read;
   bool *adapter_found;
   string *query_str;
};

bool serial_flag = false;

volatile int shared_var = 0;

bool shuffle_flag = false;

/*Report files*/
string rep_file_name1, rep_file_name2, pe_output_filename1, pe_output_filename2, shuffle_filename, se_filename, se_output_filename;

long pe1_bases_anal, pe2_bases_anal, pe_bases_kept, pe_bases_discarded, se_pe1_bases_kept, se_pe2_bases_kept;
long pe_discard_cnt;
double avg_trim_len_pe1, avg_trim_len_pe2;

string roche_output_file_name = "";
string roche_rep_file_name = "";
char* polyat_file_name; 
string polyat_output_file_name;

long se_bases_kept, se_bases_discarded;
long se_discard_cnt = 0;
long se_bases_anal = 0;        
long avg_trim_len_se;

bool wildcart_search_flag = false;

vector<char*> pe1_names, pe2_names, roche_names, se_names;

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
                max_a_error = atof(argv[++i]);
                if ((i+1)<argc && isdigit(argv[i+1][0])) 
                {
                    max_e_at_ends = atof(argv[++i]);
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
        if(string(argv[i]) == "-k" )
        {
           if ( (i+1)<argc && isdigit(argv[i+1][0]) ) 
           {
              KMER_SIZE = atoi(argv[++i]);
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
        if(string(argv[i]) == "-t" )
        {
           if ( (i+1)<argc && isdigit(argv[i+1][0]) ) 
           {
              NUM_THREADS = atoi(argv[++i]);
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
           if ( isdigit(argv[i+1][0]) ) 
                allowable_distance = atoi(argv[++i]);
           
           if(allowable_distance > 50) allowable_distance = 15;
           
           continue; 
        }
        if(string(argv[i]) == "-minimum_read_length" )
        {
           if ( isdigit(argv[i+1][0]) ) 
                minimum_read_length = atoi(argv[++i]);
           
           if(minimum_read_length > 100000) minimum_read_length = 50;
           
           continue;
        }
	if(string(argv[i]) == "-kc" )
        {
           if ( isdigit(argv[i+1][0]) ) 
                KMER_SIZE_CONT = atoi(argv[++i]);
           
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
        if( string(argv[i]) == "-12" ) //single-end file mode
        {
           if ( ( (i+1)<argc ) && (argv[i+1][0] != '-') ) 
           {
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
                        }
                }
        }
        
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
    if(polyat_flag)
    {
        if (!exists( polyat_file_name ) )
        {
                cout<< "Error: file " <<  polyat_file_name << " does not exist\n";
                return 0;
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
    }
    t.clear();
    t_prefix.clear();        
   
    string sum_stat_filename = output_prefix + "_SummaryStatistics.txt";
    sum_stat.open(sum_stat_filename.c_str(),ios::out | ios::app);
    
    
    if ( (!illumina_flag) && (!roche_flag) && (!polyat_flag) )
    {
        cout << "Not enough input parameters. Expected -454 <filename> or -1 <R1 filename>, -2 <R2 filename> arguments or -polyat <filename> arguments\n";
        sum_stat << "Not enough input parameters. Expected -454 <filename> or -1 <R1 filename>, -2 <R2 filename> arguments or -polyat <filename> arguments\n";
        sum_stat.close();
        PrintHelp();
        return 0;
    }
    
    if( (illumina_flag ) && (roche_flag ) ) 
    {
        cout << "Error! You can not use both 454 and Illumina modes in the same run. Use 454 or Illumina only, but not both!\n";
        sum_stat << "Error! You can not use both 454 and Illumina modes in the same run. Use 454 or Illumina only, but not both!\n";
        sum_stat.close();
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
        
        if(!illumina_se_flag)
        {
        
                cout << "Provided data files : " << endl;
                sum_stat << "Provided data files : " << endl;
                for(int i=0; i<(int)pe1_names.size(); ++i)
                {
                        cout << "PE1: " << pe1_names[i] << ", PE2: " << pe2_names[i] << endl;
                        sum_stat << "PE1: " << pe1_names[i] << ", PE2: " << pe2_names[i] << endl;
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
                        cout << "Maximim error: " << max_a_error << endl;
                        sum_stat << "Maximim error: " << max_a_error << endl;
                        cout << "Maximim error at ends: " << max_e_at_ends << endl;
                        sum_stat << "Maximim error at ends: " << max_e_at_ends << endl;
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
        
                cout << "--------------------Other parameters--------------------\n";
                sum_stat << "--------------------Other parameters--------------------\n";
                cout << "Maximum number of mismatches allowed in alignment: " <<  max_al_mism << endl;
                sum_stat << "Maximum number of mismatches allowed in alignment: " <<  max_al_mism << endl;
    
                cout << "Minimum read length to accept: " << minimum_read_length << endl;
                sum_stat << "Minimum read length to accept: " << minimum_read_length << endl;
        
                cout << "New to old-style Illumina headers: " << (new2old_illumina == false ? "NO" : "YES") << endl;
                sum_stat << "New to old-style Illumina headers: " << (new2old_illumina == false ? "NO" : "YES") << endl;
        }
        else
        {
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
                        cout << "Maximim error: " << max_a_error << endl;
                        sum_stat << "Maximim error: " << max_a_error << endl;
                        cout << "Maximim error at ends: " << max_e_at_ends << endl;
                        sum_stat << "Maximim error at ends: " << max_e_at_ends << endl;
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
                cout << "Maximum number of mismatches allowed in alignment: " <<  max_al_mism << endl;
                sum_stat << "Maximum number of mismatches allowed in alignment: " <<  max_al_mism << endl;
    
                cout << "Minimum read length to accept: " << minimum_read_length << endl;
                sum_stat << "Minimum read length to accept: " << minimum_read_length << endl;
        
                cout << "New to old-style Illumina headers: " << (new2old_illumina == false ? "NO" : "YES") << endl;
                sum_stat << "New to old-style Illumina headers: " << (new2old_illumina == false ? "NO" : "YES") << endl;
        }
        
    }
    if(roche_flag)
    {
        cout << "--------------------Basic parameters--------------------\n";
        sum_stat << "--------------------Basic parameters--------------------\n";
        
        cout << "Provided data files : \n" ;
        sum_stat << "Provided data files : \n" ;
        
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
           cout << "Maximim error: " << max_a_error << endl;
           sum_stat << "Maximim error: " << max_a_error << endl;
           cout << "Maximim error at ends: " << max_e_at_ends << endl;
           sum_stat << "Maximim error at ends: " << max_e_at_ends << endl;
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
        
        roche_output_file_name = output_prefix + (output_fastqfile_flag ? "_.fastq" : "_.sff" );
        roche_rep_file_name = output_prefix + "_Report.tsv" ;
        
        cout << "Report file: " << roche_rep_file_name << "\n";
        sum_stat << "Report file: " << roche_rep_file_name << "\n";
        
        cout << "Roche output file: " << roche_output_file_name << "\n";
        sum_stat << "Roche output file: " << roche_output_file_name << "\n";
        
        
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
    if(polyat_flag)
    {
        cout << "--------------------Basic parameters--------------------\n";
        sum_stat << "--------------------Basic parameters--------------------\n";
        
        cout << "Provided data file : " << polyat_file_name << "\n";
        sum_stat << "Provided data file : " << polyat_file_name << "\n";
        
        cout << string("Poly A/T trimming: ") + ( "YES. CDNA = " + string(itoa(cdna, new char[5], 10)) + ", CERR = " + string(itoa(c_err, new char[5], 10)) + ", CRNG = " + string(itoa(crng, new char[5], 10)) ) << endl;
        sum_stat << string("Poly A/T trimming: ") + ( "YES. CDNA = " + string(itoa(cdna, new char[5], 10)) + ", CERR = " + string(itoa(c_err, new char[5], 10)) + ", CRNG = " + string(itoa(crng, new char[5], 10)) ) << endl;
        
        if(qual_trim_flag)
        {
           cout << "Quality trimming: YES" << endl;
           sum_stat << "Quality trimming: YES" << endl;
           cout << "Maximim error: " << max_a_error << endl;
           sum_stat << "Maximim error: " << max_a_error << endl;
           cout << "Maximim error at ends: " << max_e_at_ends << endl;
           sum_stat << "Maximim error at ends: " << max_e_at_ends << endl;
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
        
        if(qual_trim_flag)
        {
            polyat_output_file_name = output_prefix + ( (output_fastqfile_flag || (string(polyat_file_name).substr( strlen(polyat_file_name)-5, 5 ) == "fastq") || (string(polyat_file_name).substr( strlen(polyat_file_name)-2, 2 ) == "gz") ) ? "_polyAT_qual.fastq" : "_polyAT_qual.sff" );
        }
        else
        {
            polyat_output_file_name = output_prefix + ( ( output_fastqfile_flag || (string(polyat_file_name).substr( strlen(polyat_file_name)-5, 5 ) == "fastq") || (string(polyat_file_name).substr( strlen(polyat_file_name)-2, 2 ) == "gz") ) ? "_polyAT.fastq" : "_polyAT.sff" );
        }
        
        cout << "Poly A/T output file: " << polyat_output_file_name << "\n";
        sum_stat << "Poly A/T output file: " << polyat_output_file_name << "\n";
        
        
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
    
    if(illumina_flag ) 
    {
        if(!illumina_se_flag)
        {
            IlluminaDynamic(illumina_file_name_R1, illumina_file_name_R2);
        }
        else
        {
            IlluminaDynamicSE();
        }
    }
    
    if( roche_flag  )
    {
        RocheRoutine();
    }
    
    if( polyat_flag ) 
    {
        PolyATRoutine();
    }
    
    cout << "Program finished.\n";
    sum_stat << "Program finished.\n";
    
    
    
    VectorSeqs.clear();
    VectorDict.clear();
    
    GET_TIME(finish);
    elapsed = finish - start;
    printf("Elapsed time = %e seconds\n", elapsed);
    sum_stat << "Elapsed time = " << elapsed << " seconds." << endl;
    sum_stat.close();
    
    output_prefix.clear();
    
}

void PolyATRoutine()
{
    long left_trimmed_by_polyat, right_trimmed_by_polyat, discarded_by_polyAT, bases_anal, accepted, discarded, left_trimmed_by_quality, right_trimmed_by_quality;
    left_trimmed_by_polyat = right_trimmed_by_polyat = discarded_by_polyAT = bases_anal = accepted = discarded = left_trimmed_by_quality = right_trimmed_by_quality = 0;
    
    if( string(polyat_file_name).substr( strlen(polyat_file_name)-5, 5 ) == "fastq" ) 
    { /*FASTQ file given. Process it.*/
        ParseFastqFile(polyat_file_name, reads);
        
    } 
    else if( string(polyat_file_name).substr( strlen(polyat_file_name)-3, 3 ) == "sff" ) 
    {
        process_sff_to_fastq( polyat_file_name, 0 );
    }
    else if( string(polyat_file_name).substr( strlen(polyat_file_name)-2, 2 ) == "gz" ) 
    {
        ParseFastqFile(polyat_file_name, reads);
    }
    else
    {
        cout << "Unknown file format\n";
        return;
    }
    
    if( qual_trim_flag  )
    {
         cout << "Only LUCY clipping...\n";
         QualityTrimming(reads);
    }
    
    for(int i=0; i<(int)reads.size(); i++)
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
    
    if( output_fastqfile_flag || (string(polyat_file_name).substr( strlen(polyat_file_name)-5, 5 ) == "fastq") || (string(polyat_file_name).substr( strlen(polyat_file_name)-2, 2 ) == "gz") ) 
    {
        WriteToFASTQ( polyat_output_file_name );
    }
    else
    {
        WriteToSFF( polyat_output_file_name );
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
               cout << "File is in SFF format, starting convertation...\n" ;
               process_sff_to_fastq( roche_names[i], 0 );
               if(output_fastqfile_flag == false)
               {
                 output_sfffile_flag = true;
               }
            } 
            else if(string(roche_names[i]).substr( strlen(roche_names[i])-5, 5 ) == "fastq") 
            {
               //FASTQ file given. Process it.
               ParseFastqFile(roche_names[i], reads);
               output_fastqfile_flag = true;
               output_sfffile_flag = false;
            }
       }
        
       //reads_total += reads.size();
       cout << "Convertation finished. Total number of reads read from given file(s): " << reads.size() << endl;
        
       /*If quality trimming flag is set up -> perform the quality trimming before vector/contaminants/adaptors clipping.*/
       if( qual_trim_flag  ) 
       {
           for(int i=0; i< (int)reads.size(); i++) 
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
        
       if( output_fastqfile_flag ) 
       {
           WriteToFASTQ( roche_output_file_name );
       }
       else
       {
           WriteToSFF( roche_output_file_name );
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


void RemoveContaminants454(vector<Read*>& reads454)
{
    for(int index = 0; index < (int)reads454.size(); index++)
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

void RemoveContaminants(vector<Read*>& illumina_reads)
{
    
    for(int index = 0; index < (int)illumina_reads.size(); index++)
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
    for(int ii = 0; ii < (int)reads.size(); ii++) 
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
                        reads[ii]->quality = reads[ii]->quality.substr(0, reads[ii]->read.length());
                        reads[ii]->clear_length = reads[ii]->read.length();
                }
         }
         
         
    }
}

void QualityTrimming( vector<Read*>& reads ) 
{
    for(int ii=0; ii<(int)reads.size(); ii++) 
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
           read->readID = record_block[0];
           read->initial_length = record_block[1].length();
           read->read = record_block[1];
           read->quality = record_block[3];
                        
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


void PrintHelp() {
    cout << "Version: " << version << endl;
    cout << "**********************************************************************************************************************\n";        
    cout << "Program SeqyClean\n" 
            "Main purpose of this software is to clean reads. It provide adapter/key/primers searching and quality trimming (LUCY).\n" 
            "Usage:\n"  
            "Roche 454:\n"
            "./seqyclean -454 input_file_name -o output_prefix\n"
							"[-v vector_file]\n"
							"[-c file_of_contaminants]\n"
							"[-p pcr_file_name]\n"
							"[-m file_of_RL_MIDS]\n" 
							"[-k k_mer_size]\n"
							"[-kc k_mer_size]\n"
							"[-f overlap ]\n"
							"[-t number_of_threads]\n" 
							"[-qual max_avg_error max_error_at_ends]\n"
							"[--qual_only]\n"
							"[--fastq]\n"
							"[--keep_fastq_orig]\n"
							"[-minimum_read_length <value>]\n"
							"[-polyat [cdna] [cerr] [crng] ]\n"
            "For Illumina:\n"
            "./seqyclean -1 input_file_name_1 -2 input_file_name_2 -o output_prefix\n"
							"[-v vector_file]\n"
							"[-c file_of_contaminants]\n"
							"[-k k_mer_size]\n"
							"[-kc k_mer_size]\n" 
							"[-qual max_avg_error max_error_at_ends]\n"
							"[--qual_only]\n"
							"[-minimum_read_length <value>]\n"
                                                        "[--shuffle]\n"
                                                        "[--new2old_illumina] - switch to fix read IDs ( As is detailed in: http://contig.wordpress.com/2011/09/01/newbler-input-iii-a-quick-fix-for-the-new-illumina-fastq-header/#more-342 )\n";
cout <<"Example:\n"
"Roche:\n"
"./seqyclean -454 ../artif_libs/artif454_1000_0_100_0.sff -o Test -v ../vectors.fasta -m ../RL_MIDS.csv -k 15 -f 10 -t 4\n"
"Illumina:\n"
"./seqyclean -1 P01_index16_CCGTCC_L007_R1_001.fastq.gz -2 P01_index16_CCGTCC_L007_R2_001.fastq.gz -o Test1/Test1\n";
    
    
    cout << "Please ask Ilya by email: zhba3458@vandals.uidaho.edu in case of any questions.\n" ;
 
}

void PolyAT_Trim(Read* read)
{
    int left, right;
    left = right = 0;
//    int rlength = read->read.length();
    left = poly_at_left( (char*)read->read.substr( read->lucy_lclip, read->read.length() - read->lucy_lclip ).c_str(), read->lucy_rclip - read->lucy_lclip + 1); 
    
    if (left) 
    {
        //read->lucy_lclip += left;
        //read->left_trimmed_by_polyat = 1;
        read->poly_T_clip = left;
    }
	
    right = poly_at_right((char*)read->read.substr( 0, read->lucy_rclip).c_str(), read->lucy_rclip - read->lucy_lclip + 1);
    
    if (right) 
    {
        //read->lucy_rclip -= right;
        //read->right_trimmed_by_polyat = 1;
        read->poly_A_clip = right;
    }
    
    
}


//Dynamic Illumina: does not need space to store reads:
void IlluminaDynamic(char* fastq_file1, char* fastq_file2)
{
    pe_bases_kept = pe_bases_discarded = se_pe1_bases_kept = se_pe2_bases_kept = 0;
    pe_discard_cnt = 0;
    pe1_bases_anal = pe2_bases_anal = 0;        
    avg_trim_len_pe1 = avg_trim_len_pe2 = 0;
    
    long cnt1_avg, cnt2_avg; cnt1_avg = cnt2_avg = 0; //Counters needed for calculating the average trimming length
    long cnt_avg_len1, cnt_avg_len2; cnt_avg_len1 = cnt_avg_len2 = 0;
                 
    double avg_len_pe1, avg_len_pe2; avg_len_pe1 = avg_len_pe2 = 0.0;
    double cnt_right_trim_pe1, avg_right_trim_len_pe1, cnt_right_trim_pe2, avg_right_trim_len_pe2; 
    double cnt_left_trim_pe1, avg_left_trim_len_pe1, cnt_left_trim_pe2, avg_left_trim_len_pe2;
    
    cnt_right_trim_pe1 = avg_right_trim_len_pe1 = cnt_right_trim_pe2 = avg_right_trim_len_pe2 = 0;
    cnt_left_trim_pe1 = avg_left_trim_len_pe1 = cnt_left_trim_pe2 = avg_left_trim_len_pe2 = 0;
    
    long cnt1, cnt2; cnt1 = cnt2 = 0;
    long pe_accept_cnt, se_pe1_accept_cnt, se_pe2_accept_cnt; pe_accept_cnt = se_pe1_accept_cnt = se_pe2_accept_cnt = 0;
    int ts_adapters1, ts_adapters2; ts_adapters1 = ts_adapters2 = 0;
    int num_vectors1, num_vectors2; num_vectors1 = num_vectors2 = 0;
    int num_contaminants1, num_contaminants2; num_contaminants1 = num_contaminants2 = 0;
    int accepted1,accepted2; accepted1 = accepted2 = 0;
    int discarded1,discarded2; discarded1 = discarded2 = 0;
//    int discarded_by_quality1, discarded_by_quality2; discarded_by_quality1 = discarded_by_quality2 = 0;
    int discarded_by_contaminant1, discarded_by_contaminant2; discarded_by_contaminant1 = discarded_by_contaminant2 = 0;
    int discarded_by_read_length1, discarded_by_read_length2; discarded_by_read_length1 = discarded_by_read_length2 = 0;
//    int discarded_by_vector1 , discarded_by_vector2; discarded_by_vector1 = discarded_by_vector2 = 0;
    /*Left trims*/
    int left_trimmed_by_quality1 , left_trimmed_by_quality2; left_trimmed_by_quality1 = left_trimmed_by_quality2 = 0;
    int left_trimmed_by_vector1 , left_trimmed_by_vector2; left_trimmed_by_vector1 = left_trimmed_by_vector2 = 0;
    /*Right trims/discards*/
    int right_trimmed_by_quality1 , right_trimmed_by_quality2; right_trimmed_by_quality1 = right_trimmed_by_quality2 = 0;
    int right_trimmed_by_adapter1 , right_trimmed_by_adapter2; right_trimmed_by_adapter1 = right_trimmed_by_adapter2 = 0;
    int right_trimmed_by_vector1 , right_trimmed_by_vector2;  right_trimmed_by_vector1 = right_trimmed_by_vector2 = 0;
    
    fstream rep_file1, rep_file2, pe_output_file1, pe_output_file2, shuffle_file, se_file;
    rep_file1.open(rep_file_name1.c_str(),ios::out);
    rep_file2.open(rep_file_name2.c_str(),ios::out);
    rep_file1 << "ReadID\tlclip\trclip\tTruSeq_pos\tTruSeq_type\tRaw_read_length\tLlucy\tRlucy\tDiscarded\tContaminants\tVectorID\tVecStart\tVecEnd\tVecLen\n";
    rep_file2 << "ReadID\tlclip\trclip\tTruSeq_pos\tTruSeq_type\tRaw_read_length\tLlucy\tRlucy\tDiscarded\tContaminants\tVectorID\tVecStart\tVecEnd\tVecLen\n";
    
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
    
    for(int jj=0; jj<pe1_names.size(); ++jj)
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
                        
                        if(fields1[0] != fields2[0])
                        {
                            cout << "Warning: read IDs do not match in input files: PE1-> " << pe1_names[jj] << ", PE2-> " << pe2_names[jj] << endl;
                            sum_stat << "Warning: read IDs do not match in input files: PE1-> " << pe1_names[jj] << ", PE2-> " << pe2_names[jj] << endl;
                        }
                        
                        fields1.clear();
                        fields2.clear();
                        
                        if (new2old_illumina == true)
                        {
                            split_str( line1, fields1, " " );
                            split_str( fields1[0], fields2, ":" );
                            line1 = string(fields2[0] + "_" + fields2[2] + ":" + fields2[3] + ":" + fields2[4] + ":" + fields2[5] + ":" + fields2[6] + "#0/" + fields1[1].substr(0,1) + " (" + line1 + ")");
                            
                            fields1.clear();
                            fields2.clear();
                            
                            split_str( line2, fields1, " " );
                            split_str( fields1[0], fields2, ":" );
                            
                            line2 = string(fields2[0] + "_" + fields2[2] + ":" + fields2[3] + ":" + fields2[4] + ":" + fields2[5] + ":" + fields2[6] + "#0/" + fields1[1].substr(0,1) + " (" + line2 + ")");
                            
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
                        read1->readID = record_block1[0];
                        read1->initial_length = record_block1[1].length();
                        read1->read = record_block1[1];
                        read1->quality = line1;
                        pe1_bases_anal += read1->read.length();
          
                        Read *read2 = new Read();
                        read2->readID = record_block2[0];
                        read2->initial_length = record_block2[1].length();
                        read2->read = record_block2[1];
                        read2->quality = line2;
                        pe2_bases_anal += read2->read.length();
          
                        //Serial realization - useful for debugging if something does not work as expected
          
                        IlluminaDynRoutine(read1, adapter_found1, query_string1);
                        IlluminaDynRoutine(read2, adapter_found2, query_string2);
                        
                        cnt1+=1; cnt2+=1;
          
                        //Report
                        //Read ID
                        rep_file1 << read1->readID << "\t" << read1->lclip << "\t" << read1->rclip << "\t" << read1->tru_sec_pos << "\t" << read1->b_adapter << "\t" << read1->initial_length << "\t" << (read1->lucy_lclip <= 1 ? 1 : read1->lucy_lclip) << "\t" << (read1->lucy_rclip <= 1 ? 1 : read1->lucy_rclip) << "\t" << read1->discarded << "\t" << read1->contaminants << "\t" << "NA" << "\n";
                        rep_file2 << read2->readID << "\t" << read2->lclip << "\t" << read2->rclip << "\t" << read2->tru_sec_pos << "\t" << read2->b_adapter << "\t" << read2->initial_length << "\t" << (read2->lucy_lclip <= 1 ? 1 : read2->lucy_lclip) << "\t" << (read2->lucy_rclip <= 1 ? 1 : read2->lucy_rclip) << "\t" << read2->discarded << "\t" << read2->contaminants << "\t" << "NA" << "\n";
          
                        if( !VectorOnlyFlag )  
                        {
                                if( read1->lclip >= read1->rclip ) { read1->discarded = 1; read1->discarded_by_read_length = 1; } 
                                if( read1->lclip >= (int)read1->read.length() ) { read1->discarded = 1; read1->discarded_by_read_length = 1; }
                                if( read1->rclip > (int)read1->read.length() ) { read1->rclip = read1->read.length(); }
                                if( (int)read1->read.length() < minimum_read_length ) { read1->discarded = 1; read1->discarded_by_read_length = 1; }
                                if( (read1->rclip - read1->lclip) < minimum_read_length ) { read1->discarded = 1; read1->discarded_by_read_length = 1; }
              
                                if( read2->lclip >= read2->rclip ) {read2->discarded = 1; read2->discarded_by_read_length = 1; }
                                if( read2->lclip >= (int)read2->read.length() ) {read2->discarded = 1; read2->discarded_by_read_length = 1; }
                                if( read2->rclip > (int)read2->read.length() ) {read2->rclip = read2->read.length(); }
                                if( (int)read2->read.length() < minimum_read_length ) {read2->discarded = 1; read2->discarded_by_read_length = 1; }
                                if( (read2->rclip - read2->lclip) < minimum_read_length ) { read2->discarded = 1; read2->discarded_by_read_length = 1; }

                                if( (read1->discarded == 0) && (read2->discarded == 0) )
                                {
                                        
                                    if(  read1->rclip < read1->initial_length  )
                                    {
                                        cnt_right_trim_pe1 += 1;
                                        avg_right_trim_len_pe1 = GetAvg( avg_right_trim_len_pe1, cnt_right_trim_pe1, read1->initial_length - read1->rclip );
                                    }
                                    if(read1->lclip > 0)
                                    {
                                        cnt_left_trim_pe1 += 1;
                                        avg_left_trim_len_pe1 = GetAvg( avg_left_trim_len_pe1, cnt_left_trim_pe1, read1->lclip );
                                    }
                                    
                                    read1->read = read1->read.substr(0 , read1->rclip );
                                    read1->quality = read1->quality.substr(0,read1->rclip) ; 
                                    read1->read = read1->read.substr( read1->lclip, read1->rclip - read1->lclip );
                                    read1->quality = read1->quality.substr( read1->lclip, read1->rclip - read1->lclip );
                 
                                    if(  read2->rclip < read2->initial_length  )
                                    {
                                        cnt_right_trim_pe2 += 1;
                                        avg_right_trim_len_pe2 = GetAvg( avg_right_trim_len_pe2, cnt_right_trim_pe2, read2->initial_length - read2->rclip );
                                    }
                                    if(read2->lclip > 0)
                                    {
                                        cnt_left_trim_pe2 += 1;
                                        avg_left_trim_len_pe2 = GetAvg( avg_left_trim_len_pe2, cnt_left_trim_pe2, read2->lclip );
                                    }
                                        
                                    read2->read = read2->read.substr(0 , read2->rclip );
                                    read2->quality = read2->quality.substr(0,read2->rclip) ; 
                                    read2->read = read2->read.substr( read2->lclip, read2->read.length() - read2->lclip );
                                    read2->quality = read2->quality.substr( read2->lclip, read2->quality.length() - read2->lclip );
        	
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
                 
                                    if( read1->initial_length > read1->read.length() )
                                    {
                                        cnt1_avg+=1;
                                        avg_trim_len_pe1 = GetAvg( avg_trim_len_pe1, cnt1_avg, read1->initial_length - read1->read.length() );
                                    }
                 
                                    if( read2->initial_length > read2->read.length() )
                                    {
                                        cnt2_avg+=1;
                                        avg_trim_len_pe2 = GetAvg( avg_trim_len_pe2, cnt2_avg, read2->initial_length - read2->read.length() );
                                    }
                 
                                    cnt_avg_len1 +=1; cnt_avg_len2 +=1;
                 
                                    avg_len_pe1 = GetAvg( avg_len_pe1, cnt_avg_len1, read1->read.length() );
                                    avg_len_pe2 = GetAvg( avg_len_pe2, cnt_avg_len2, read2->read.length() );
                
                 
                                } else if ((read1->discarded == 0) && (read2->discarded == 1)) 
                                {
                                    
                                    read1->read = read1->read.substr(0 , read1->rclip );
                                    read1->quality = read1->quality.substr(0,read1->rclip) ; 
                                    read1->read = read1->read.substr( read1->lclip, read1->rclip - read1->lclip );
                                    read1->quality = read1->quality.substr( read1->lclip, read1->rclip - read1->lclip );
                 
                                    WriteSEFile( se_file, read1 );
                                    se_pe1_accept_cnt+=1;
                                    se_pe1_bases_kept += read1->read.length();
                 
                                    
                                } else if( (read1->discarded == 1) && (read2->discarded == 0) )
                                {
                                      
                                    read2->read = read2->read.substr(0 , read2->rclip );
                                    read2->quality = read2->quality.substr(0,read2->rclip) ; 
                                    read2->read = read2->read.substr( read2->lclip, read2->read.length() - read2->lclip );
                                    read2->quality = read2->quality.substr( read2->lclip, read2->quality.length() - read2->lclip );
        	 
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
                        } 
                        else if( VectorOnlyFlag  )  
                        {
                                read1->rclip = read1->v_end;
                                read1->lclip = read1->v_start+1;
                                read2->rclip = read2->v_end;
                                read2->lclip = read2->v_start+1;
              
                                if ( (read1->discarded_by_vector == 1) && (read1->discarded_by_quality == 0) && (read1->discarded_by_contaminant == 0) )
                                {
                                        read1->discarded = 0;
                                }
                                if ( (read2->discarded_by_vector == 1) && (read2->discarded_by_quality == 0) && (read2->discarded_by_contaminant == 0) )
                                {
                                        read2->discarded = 0;
                                }
              
                                if( read1->lclip >= (int)read1->rclip ) {read1->discarded = 1; read1->discarded_by_read_length = 1; }
                                if( read1->lclip >= (int)read1->read.length() ) {read1->discarded = 1; read1->discarded_by_read_length = 1; }
                                if( read1->rclip > (int)read1->read.length() ) {read1->rclip = read1->read.length(); }
                                //if( read1->read.length() < minimum_read_length ) {read1->discarded = 1; read1->discarded_by_read_length = 1; }
                
                                if( read2->lclip >= (int)read2->rclip ) {read2->discarded = 1; read2->discarded_by_read_length = 1; }
                                if( read2->lclip >= (int)read2->read.length() ) {read2->discarded = 1; read2->discarded_by_read_length = 1; }
                                if( read2->rclip > (int)read2->read.length() ) {read2->rclip = read2->read.length(); }
                                //if( read2->read.length() < minimum_read_length ) {read2->discarded = 1; read2->discarded_by_read_length = 1; }
              
                                if( (read1->discarded == 0) && (read2->discarded == 0) )
                                {
                                        read1->read = read1->read.substr(0 , read1->rclip );
                                        read1->quality = read1->quality.substr(0,read1->rclip) ; 
                                        read1->read = read1->read.substr( read1->lclip-1, read1->read.length() - read1->lclip-1 );
                                        read1->quality = read1->quality.substr( read1->lclip-1, read1->quality.length() - read1->lclip-1 );
                 
                                        read2->read = read2->read.substr(0 , read2->rclip );
                                        read2->quality = read2->quality.substr(0,read2->rclip) ; 
                                        read2->read = read2->read.substr( read2->lclip-1, read2->read.length() - read2->lclip-1 );
                                        read2->quality = read2->quality.substr( read2->lclip-1, read2->quality.length() - read2->lclip-1 );
        
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
                                 } else if ((read1->discarded == 0) && (read2->discarded == 1)) 
                                 {
                                        WriteSEFile( se_file, read1 ); 
                                        se_pe1_accept_cnt += 1;
                                        se_pe1_bases_kept += read1->read.length();
                                 } else if( (read1->discarded == 1) && (read2->discarded == 0) )
                                 {
                                        WriteSEFile( se_file, read2 ); 
                                        se_pe2_accept_cnt += 1;
                                        se_pe2_bases_kept += read2->read.length();
                                 }
                        }
          
                        if (read1->tru_sec_found == 1) ts_adapters1++;
                        if (read1->vector_found == 1) num_vectors1++;
                        if (read1->contam_found == 1) num_contaminants1++;
                        if (read1->discarded == 0) accepted1++;
                        if (read1->discarded == 1) discarded1++;
                        //if (read1->discarded_by_quality == 1) discarded_by_quality1++;    
                        if (read1->discarded_by_contaminant == 1) discarded_by_contaminant1++;
                        if (read1->discarded_by_read_length == 1) discarded_by_read_length1++;
                        //if (read1->discarded_by_vector == 1) discarded_by_vector1++;
                        if (read1->left_trimmed_by_quality == 1) left_trimmed_by_quality1++;
                        if (read1->left_trimmed_by_vector == 1) left_trimmed_by_vector1++;
                        if (read1->right_trimmed_by_quality == 1) right_trimmed_by_quality1++;
                        if (read1->right_trimmed_by_adapter == 1) right_trimmed_by_adapter1++;
                        if (read1->right_trimmed_by_vector == 1) right_trimmed_by_vector1++;
          
                        if (read2->tru_sec_found == 1) ts_adapters2++;
                        if (read2->vector_found == 1) num_vectors2++;
                        if (read2->contam_found == 1) num_contaminants2++;
                        if (read2->discarded == 0) accepted2++;
                        if (read2->discarded == 1) discarded2++;
                        //if (read2->discarded_by_quality == 1) discarded_by_quality2++;    
                        if (read2->discarded_by_contaminant == 1) discarded_by_contaminant2++;
                        if (read2->discarded_by_read_length == 1) discarded_by_read_length2++;
                        //if (read2->discarded_by_vector == 1) discarded_by_vector2++;
                        if (read2->left_trimmed_by_quality == 1) left_trimmed_by_quality2++;
                        if (read2->left_trimmed_by_vector == 1) left_trimmed_by_vector2++;
                        if (read2->right_trimmed_by_quality == 1) right_trimmed_by_quality2++;
                        if (read2->right_trimmed_by_adapter == 1) right_trimmed_by_adapter2++;
                        if (read2->right_trimmed_by_vector == 1) right_trimmed_by_vector2++;
          
                        record_block1.clear();
                        read1->readID.clear(); 
                        read1->quality.clear();
                        read1->read.clear();
          
                        record_block2.clear();
                        read2->readID.clear(); 
                        read2->quality.clear();
                        read2->read.clear();
          
                        delete read1;
                        delete read2;
          
          
          
                }
        }
        in1.close();
        in2.close();
    }
    
    pe_output_file1.close();
    pe_output_file2.close();
    se_file.close();
    shuffle_file.close();
    
    rep_file1.close();
    rep_file2.close();
    
    
    
    cout << "====================Summary Statistics====================\n";
    sum_stat << "====================Summary Statistics====================\n";
    cout << "PE1 reads analyzed: " << cnt1 << ", Bases:" << pe1_bases_anal << "\n";
    sum_stat << "PE1 reads analyzed: " << cnt1 << ", Bases:" << pe1_bases_anal << "\n";
    
    cout << "Found ->\n";
    sum_stat << "Found ->\n";
    
    cout << "Adapters: " << ts_adapters1 << ", " << ( (double)ts_adapters1/(double)cnt1*100.0) << "%\n";
    sum_stat << "Adapters: " << ts_adapters1 << ", " << ( (double)ts_adapters1/(double)cnt1*100.0) << "%\n";
    
    if(vector_flag)
    {
        cout << "# of reads with vector: " << num_vectors1 << ", " << ( (double)num_vectors1/(double)cnt1*100.0) << "%\n";
        sum_stat << "# of reads with vector: " << num_vectors1 << ", " << ( (double)num_vectors1/(double)cnt1*100.0) << "%\n";
    }
    
    if(contaminants_flag)
    {
        cout << "# of reads with contaminants: " << num_contaminants1 << ", " << ( (double)num_contaminants1/(double)cnt1*100.0) << "%\n";
        sum_stat << "# of reads with contaminants: " << num_contaminants1 << ", " << ( (double)num_contaminants1/(double)cnt1*100.0) << "%\n";
    }
    
    if(qual_trim_flag || vector_flag )
    {
        cout << "Reads left trimmed ->" << "\n";
        sum_stat << "Reads left trimmed ->" << "\n";
    }
    
    if(qual_trim_flag)
    {
        cout << "By quality: " <<  left_trimmed_by_quality1 << "\n";
        sum_stat << "By quality: " <<  left_trimmed_by_quality1 << "\n";
    }
    
    if(vector_flag)
    {
        cout << "By vector: " <<  left_trimmed_by_vector1 << "\n";
        sum_stat << "By vector: " <<  left_trimmed_by_vector1 << "\n";
    }
    
    cout << "Average left trimmed length: " << avg_left_trim_len_pe1 << " bp\n";
    sum_stat << "Average left trimmed length: " << avg_left_trim_len_pe1 << " bp\n";
    
    cout << "Reads right trimmed ->" << "\n";
    sum_stat << "Reads right trimmed ->" << "\n";
    
    cout << "By adapter: " <<  right_trimmed_by_adapter1 << "\n";
    sum_stat << "By adapter: " <<  right_trimmed_by_adapter1 << "\n";
    
    if(qual_trim_flag)
    {
        cout << "By quality: " <<  right_trimmed_by_quality1 << "\n";
        sum_stat << "By quality: " <<  right_trimmed_by_quality1 << "\n";
    }
    
    if(vector_flag)
    {
        cout << "By vector: " <<  right_trimmed_by_vector1 << "\n";
        sum_stat << "By vector: " <<  right_trimmed_by_vector1 << "\n";
    }
    
    cout << "Average right trimmed length: " << avg_right_trim_len_pe1 << " bp\n";
    sum_stat << "Average right trimmed length: " << avg_right_trim_len_pe1 << " bp\n";
    
    cout << "PE1 reads discarded: " << discarded1 << "\n";
    sum_stat << "PE1 reads discarded: " << discarded1 << "\n";
    
    if(contaminants_flag)
    {
        cout << "By contaminants: " <<  discarded_by_contaminant1 << "\n";
        sum_stat << "By contaminant: " <<  discarded_by_contaminant1 << "\n";
    }
    
    cout << "By read length: " <<  discarded_by_read_length1 << "\n";
    sum_stat << "By read length: " <<  discarded_by_read_length1 << "\n";
    
    cout << "-----------------------------------------------------------\n";
    sum_stat << "-----------------------------------------------------------\n";
    
    cout << "PE2 reads analyzed: " << cnt2 << ", Bases:" << pe2_bases_anal << "\n";
    sum_stat << "PE2 reads analyzed: " << cnt2 << ", Bases:" << pe2_bases_anal << "\n";
    
    cout << "Found ->\n";
    sum_stat << "Found ->\n";
    
    cout << "Adapters: " << ts_adapters2 << ", " << ( (double)ts_adapters2/(double)cnt2*100.0) << "%\n";
    sum_stat << "Adapters: " << ts_adapters2 << ", " << ( (double)ts_adapters2/(double)cnt2*100.0) << "%\n";
    
    if(vector_flag)
    {
        cout << "# of reads with vector: " << num_vectors2 << ", " << ( (double)num_vectors2/(double)cnt2*100.0) << "%\n";
        sum_stat << "# of reads with vector: " << num_vectors2 << ", " << ( (double)num_vectors2/(double)cnt2*100.0) << "%\n";
    }
    
    if(contaminants_flag)
    {
       cout << "# of reads with contaminants: " << num_contaminants2 << ", " << ( (double)num_contaminants2/(double)cnt2*100.0) << "%\n";
        sum_stat << "# of reads with contaminants: " << num_contaminants2 << ", " << ( (double)num_contaminants2/(double)cnt2*100.0) << "%\n";
    }
    
    if(qual_trim_flag || vector_flag )
    {
        cout << "Reads left trimmed ->" << "\n";
        sum_stat << "Reads left trimmed ->" << "\n";
    }
    
    if(qual_trim_flag)
    {
        cout << "By quality: " <<  left_trimmed_by_quality2 << "\n";
        sum_stat << "By quality: " <<  left_trimmed_by_quality2 << "\n";
    }
    
    if(vector_flag)
    {
        cout << "By vector: " <<  left_trimmed_by_vector2 << "\n";
        sum_stat << "By vector: " <<  left_trimmed_by_vector2 << "\n";
    }
    
    cout << "Average left trimmed length: " << avg_left_trim_len_pe2 << " bp\n";
    sum_stat << "Average left trimmed length: " << avg_left_trim_len_pe2 << " bp\n";
    
    cout << "Reads right trimmed ->" << "\n";
    sum_stat << "Reads right trimmed ->" << "\n";
    
    if(qual_trim_flag)
    {
        cout << "By quality: " <<  right_trimmed_by_quality2 << "\n";
        sum_stat << "By quality: " <<  right_trimmed_by_quality2 << "\n";
    }
    
    if(vector_flag)
    {
        cout << "By vector: " <<  right_trimmed_by_vector2 << "\n";
        sum_stat << "By vector: " <<  right_trimmed_by_vector2 << "\n";
    }
    
    cout << "By adapter: " <<  right_trimmed_by_adapter2 << "\n";
    sum_stat << "By adapter: " <<  right_trimmed_by_adapter2 << "\n";
    
    cout << "Average right trimmed length: " << avg_right_trim_len_pe2 << " bp\n";
    sum_stat << "Average right trimmed length: " << avg_right_trim_len_pe2 << " bp\n";
    
    cout << "PE2 reads discarded:" << discarded2 << "\n";
    sum_stat << "PE2 reads discarded: " << discarded2 << "\n";
    
    if(contaminants_flag)
    {
        cout << "By contaminants: " <<  discarded_by_contaminant2 << "\n";
        sum_stat << "By contaminants: " <<  discarded_by_contaminant2 << "\n";
    }
    
    cout << "By read length: " <<  discarded_by_read_length2 << "\n";
    sum_stat << "By read length: " <<  discarded_by_read_length2 << "\n";
    
    
    cout << "----------------------Summary for PE & SE----------------------\n";
    sum_stat << "----------------------Summary for PE & SE----------------------\n";
    
    cout << "Pairs kept: " << pe_accept_cnt << ", " << ( (double)pe_accept_cnt/(double)cnt1*100.0) << "%, Bases: " << pe_bases_kept << ", " << ( (double)pe_bases_kept/(double)(pe1_bases_anal+pe2_bases_anal)*100) <<  "%\n";
    sum_stat << "Pairs kept: " << pe_accept_cnt << ", " << ( (double)pe_accept_cnt/(double)cnt1*100.0) << "%, Bases: " << pe_bases_kept << ", " << ( (double)pe_bases_kept/(double)(pe1_bases_anal+pe2_bases_anal)*100) <<  "%\n";
    
    cout << "Pairs discarded: " << pe_discard_cnt << ", " << ( (double)pe_discard_cnt/(double)cnt1*100.0) << "%, Bases: " << pe_bases_discarded << ", " << ( (double)pe_bases_discarded/(double)(pe1_bases_anal+pe2_bases_anal)*100) <<  "%\n";
    sum_stat << "Pairs discarded: " << pe_discard_cnt << ", " << ( (double)pe_discard_cnt/(double)cnt1*100.0) << "%, Bases: " << pe_bases_discarded << ", " << ( (double)pe_bases_discarded/(double)(pe1_bases_anal+pe2_bases_anal)*100) <<  "%\n";
    
    cout << "Single Reads PE1 kept: " << se_pe1_accept_cnt << ", Bases: " << se_pe1_bases_kept << "\n";
    sum_stat << "Single Reads PE1 kept: " << se_pe1_accept_cnt << ", Bases: " << se_pe1_bases_kept << "\n";
    
    cout << "Single Reads PE2 kept: " << se_pe2_accept_cnt << ", Bases: " << se_pe2_bases_kept << "\n";
    sum_stat << "Single Reads PE2 kept: " << se_pe2_accept_cnt << ", Bases: " << se_pe2_bases_kept << "\n";
    
    cout << "Average trimmed length PE1: " << avg_trim_len_pe1 << " bp\n";
    sum_stat << "Average trimmed length PE1: " << avg_trim_len_pe1 << " bp\n";
    
    cout << "Average trimmed length PE2: " << avg_trim_len_pe2 << " bp\n";
    sum_stat << "Average trimmed length PE2: " << avg_trim_len_pe2 << " bp\n";
    
    cout << "Average read length PE1: " << avg_len_pe1 << " bp\n";
    sum_stat << "Average read length PE1: " << avg_len_pe1 << " bp\n";
    
    cout << "Average read length PE2: " << avg_len_pe2 << " bp\n";
    sum_stat << "Average read length PE2: " << avg_len_pe2 << " bp\n";
    
    cout << "====================Done cleaning====================\n";  
    sum_stat << "====================Done cleaning====================\n";  
   
    
}

int IlluminaDynRoutine(Read* read, bool& adapter_found, string &query_str)
{
    //Remove not-needed Ns:
    TrimNs( read->read );
    
    if(read->read.length() > minimum_read_length)
    {
        read->quality = read->quality.substr(0, read->read.length());
        read->clear_length = read->read.length();
    }
    else
    {
        read->discarded = 1;
        read->discarded_by_read_length = 1;
        return -1;
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
       QualTrim( read, max_a_error, max_e_at_ends );//This function generates LUCY clips of the read. Later they should be compared and read should be trimmed based on the results of comparison.
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
    //Run the main routine: Adapter + Vector/Contaminants trimming or only Adaptors
    //First 20 bases of i5 adapter forward
    size_t found;
    if (!adapter_found)
    {
        string ts_adapter = tmpl_i5_1.substr(0,15);
        found = read->read.find( ts_adapter );
        if( found != string::npos ) 
        {
            cout << "i5 adapter in forward first found in the read #" << found << endl; 
            sum_stat << "i5 adapter in forward first found in the read #" << found << endl; 
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
                cout << "i5 adapter in forward first found in the read #" << found << endl; 
                sum_stat << "i5 adapter in forward first found in the read #" << found << endl; 
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
                    cout << "i7 adapter in forward first found in the read #" << found << endl;
                    sum_stat << "i7 adapter in forward first found in the read #" << found << endl;
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
                        cout << "i7 adapter in forward first found in the read #" << found << endl; 
                        sum_stat << "i7 adapter in forward first found in the read #" << found << endl;
                        adapter_found = true;
                        query_str = ts_adapter;
                        read->tru_sec_pos = found;
                        read->tru_sec_found = 1;
                    } 
                    else
                    {
                        if( vector_flag ) 
                           CheckVector(read); 
                        
                        read->tru_sec_pos = read->clear_length;
                        read->tru_sec_found = 0;
         
                    }
                }
             }
         }
    }
    else
    {
        if( vector_flag ) 
            CheckVector(read); 
       
        bool adp_found = false;
        found = read->read.rfind( query_str );
        if( found != string::npos ) 
        {
            adp_found = true;
            read->tru_sec_pos = found;
            read->b_adapter = adapter_type_R1;
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
                        read->b_adapter = adapter_type_R1;
                        read->tru_sec_found = 1;
                    }
            }
         
            if(!adp_found) 
            {
                read->tru_sec_pos = read->clear_length;
                read->tru_sec_found = 0;
            }
         
            delete izssaha;
        }
    }
            
    //Clip points
    if( (qual_trim_flag ) && (vector_flag ) )
    {
        if(read->vector_found == 1)
        {
           if( read->v_start >= (read->read.length() - read->v_end) ) //Vector is on the right side
           {
               read->lclip = max(read->lucy_lclip, 1);
               if( (read->lclip == read->lucy_lclip) && (read->lucy_lclip > 1)) 
               {
                   read->left_trimmed_by_quality = 1;
               }
            
               read->rclip = min(read->tru_sec_pos, min(read->lucy_rclip, read->v_start) );
                    
               if( read->rclip == read->lucy_rclip )
               {
                 read->right_trimmed_by_quality = 1;
               }
               else if(read->rclip == read->tru_sec_pos)
               {
                 read->right_trimmed_by_adapter = 1;
               }
               else if(read->rclip == read->v_start)
               {
                 read->right_trimmed_by_vector = 1;
               }
           }
           else //Vector is on the left side or the whole read is vector
           {
               read->lclip = max(read->lucy_lclip,max(1, read->v_end ) );
           
               if( (read->lclip == read->lucy_lclip) && (read->lucy_lclip > 1)) 
               {
                 read->left_trimmed_by_quality = 1;
               }
               if(read->lclip == read->v_end)
               {
                 read->left_trimmed_by_vector = 1;
               }
           
               read->rclip = min(read->tru_sec_pos, read->lucy_rclip );
               if(read->rclip == read->lucy_rclip)
               {
                 read->right_trimmed_by_quality = 1;
               }
               else
               {
                 read->right_trimmed_by_adapter = 1;
               }
           }
                
           if(read->v_end - read->v_start > 50) 
           {
               read->discarded = 1;
               //read->discarded_by_vector = 1;
               read->discarded_by_read_length = 1;
           }
        } 
        else
        {
            read->lclip = max(read->lucy_lclip, 1);
            if( (read->lclip == read->lucy_lclip) && (read->lucy_lclip > 1))  read->left_trimmed_by_quality = 1;
            
            read->rclip = min(read->tru_sec_pos, read->lucy_rclip );
            if(read->rclip == read->lucy_rclip)
            {
              read->right_trimmed_by_quality = 1;
            }
            else
            {
              read->right_trimmed_by_adapter = 1;
            }
        }
    }
    else if( (qual_trim_flag ) && (!vector_flag) )
    {
        read->lclip = max(read->lucy_lclip,1);
        if( (read->lclip == read->lucy_lclip) && (read->lucy_lclip > 1)) 
           read->left_trimmed_by_quality = 1;
                
        read->rclip = min(read->tru_sec_pos,read->lucy_rclip );
        if(read->rclip == read->lucy_rclip)
        {
           read->right_trimmed_by_quality = 1;
        }
        else if(read->rclip == read->tru_sec_pos)
        {
           read->right_trimmed_by_adapter = 1;
        }
                
        if(read->rclip >= read->clear_length)
        {
            read->rclip = read->clear_length; read->right_trimmed_by_adapter = 0;
        }   
    }
    else if( (!qual_trim_flag) && (vector_flag ) )
    {
       if( read->v_start >= (read->read.length() - read->v_end) ) //Vector is on the right side
       {
          read->lclip = 1;
            
          read->rclip = min(read->tru_sec_pos, read->v_start == -1 ? read->tru_sec_pos : read->v_start );
                    
          if(read->rclip == read->tru_sec_pos)
          {
             read->right_trimmed_by_adapter = 1;
          }
          else if(read->rclip == read->v_start)
          {
             read->right_trimmed_by_vector = 1;
          }
           
       }
       else //Vector is on the left side
       {
          read->lclip = max(1, read->v_end == -1 ? 1 : read->v_end);
           
          if(read->lclip == read->v_end)
          {
             read->left_trimmed_by_vector = 1;
          }
           
          read->rclip = read->tru_sec_pos;
          read->right_trimmed_by_adapter = 1;
       }
                
       if(read->rclip >= read->clear_length)
       {
          read->rclip = read->clear_length; read->right_trimmed_by_adapter = 0;
       }    
    }
    else if((!qual_trim_flag) && (!vector_flag))
    {
       read->rclip = read->tru_sec_pos;
       if(read->rclip < read->clear_length)
          read->right_trimmed_by_adapter = 1;
            
       read->lclip = 1;
    }
    
    return 0;
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

void WritePEFile(fstream &pe_output_file, Read *read)
{
    pe_output_file << read->readID << endl;
    pe_output_file << read->read << endl;
    pe_output_file << '+' << endl;
    pe_output_file << read->quality << endl;
}

void WriteShuffleFile(fstream &shuffle_output_file, Read *read1, Read *read2)
{
    shuffle_output_file << read1->readID << endl;
    shuffle_output_file << read1->read << endl;
    shuffle_output_file << '+' << endl;
    shuffle_output_file << read1->quality << endl;
    
    shuffle_output_file << read2->readID << endl;
    shuffle_output_file << read2->read << endl;
    shuffle_output_file << '+' << endl;
    shuffle_output_file << read2->quality << endl;
}

void WriteSEFile(fstream &se_output_file, Read *read)
{
    /*vector <string> fields1, fields2;
    split_str( read->readID, fields1, " " );
    split_str( fields1[0], fields2, ":" );
    read->readID = fields2[0] + "_" + fields2[2] + ":" + fields2[3] + ":" + fields2[4] + ":" + fields2[5] + ":" + fields2[6] + ( fields1[1].substr(0,1) == "1" ? "_R1" : "_R2" );
    **/
    se_output_file << read->readID << endl;
    se_output_file << read->read << endl;
    se_output_file << '+' << endl;
    se_output_file << read->quality << endl;
    
    //fields1.clear(); fields2.clear();
}




void *t_IlluminaDynRoutine( void *targs )
{
    struct thread_data *t_data;
    
    t_data = (thread_data *) targs;
    
    //Remove not-needed Ns:
    TrimNs( t_data->read->read );
    t_data->read->quality = t_data->read->quality.substr(0, t_data->read->read.length());
    t_data->read->clear_length = t_data->read->read.length();
    
          //Check validity
    if (t_data->read->clear_length < minimum_read_length) 
    {
       t_data->read->discarded = 1;
       t_data->read->discarded_by_read_length = 1;
       
       pthread_exit(NULL);return 0;
    }
    
    if(contaminants_flag )
    {
       if(CheckContaminants(t_data->read->read) == 0) 
       {
           t_data->read->contam_found = 1;
           t_data->read->discarded_by_contaminant = 1;
           t_data->read->contaminants = 1;
           t_data->read->discarded = 1;
           
           pthread_exit(NULL);return 0;
       }
    }
            
        
    //If quality trimming flag is set up -> perform the quality trimming before vector/contaminants/adaptors clipping.
    if( qual_trim_flag  ) 
    {
       QualTrim( t_data->read, max_a_error, max_e_at_ends );//This function generates LUCY clips of the read. Later they should be compared and read should be trimmed based on the results of comparison.
            
       if( lucy_only_flag  ) 
       {
          t_data->read->lclip = t_data->read->lucy_lclip;
          t_data->read->rclip = t_data->read->lucy_rclip;
       }
    }
         
    //Run the main routine: Adapter + Vector/Contaminants trimming or only Adaptors
    //First 20 bases of i5 adapter forward
    size_t found;
    if (!(*t_data->adapter_found))
    {
        string ts_adapter = tmpl_i5_1.substr(0,15);
        found = t_data->read->read.find( ts_adapter );
        if( found != string::npos ) 
        {
            //cout << found << " i5 adapter in forward" << endl; 
            *(t_data->adapter_found) = true;
            *(t_data->query_str) = ts_adapter;
            t_data->read->lclip = 1;
            t_data->read->rclip = found;
            t_data->read->tru_sec_pos = found;
            t_data->read->tru_sec_found = 1;
        }
        else
        {
            //First 20 bases of i5 adapter in reverse complement
            ts_adapter = MakeRevComplement(tmpl_i5_2).substr(0,15);
            found = t_data->read->read.find( ts_adapter );
            if( found != string::npos ) 
            {
                cout << found << " i5 adapter in reverse complement" << endl; 
                *(t_data->adapter_found) = true;
                *(t_data->query_str) = ts_adapter;
                t_data->read->lclip = 1;
                t_data->read->rclip = found;
                t_data->read->tru_sec_pos = found;
                t_data->read->tru_sec_found = 1;
            }
            else
            {
                //First 20 bases of i7 adapter forward
                ts_adapter = tmpl_i7_1.substr(0,15);
                found = t_data->read->read.find( ts_adapter );
                if( found != string::npos ) 
                {
                    if (found > 0)
                    {
                        cout << found << " i7 adapter in forward" << endl;
                        *(t_data->adapter_found) = true;
                        *(t_data->query_str) = ts_adapter;
                        t_data->read->lclip = 1;
                        t_data->read->rclip = found;
                        t_data->read->tru_sec_pos = found;
                        t_data->read->tru_sec_found = 1;
                    }
                } 
                else
                {
                    //First 20 bases of i5 adapter in reverse complement
                    ts_adapter = MakeRevComplement(tmpl_i7_2).substr(0,15);
                    found = t_data->read->read.find( ts_adapter );
                    if( found != string::npos ) 
                    {
                        cout << found << " i7 adapter in reverse complement" << endl; 
                        *(t_data->adapter_found) = true;
                        *(t_data->query_str) = ts_adapter;
                    } 
                    else
                    {
                        if(vector_flag ) 
                        {
                            CheckVector(t_data->read); 
                            if( ( t_data->read->read.length() - (t_data->read->v_end - t_data->read->v_start) ) < minimum_read_length ) 
                            {
                                t_data->read->discarded = 1;
                                t_data->read->discarded_by_vector = 1;
                            }
                        }
                    
                        t_data->read->tru_sec_pos = t_data->read->clear_length;
         
                        t_data->read->lclip = 1;
                        t_data->read->rclip = t_data->read->tru_sec_pos;
         
                        if(t_data->read->rclip < minimum_read_length) 
                        {              
                            t_data->read->discarded = 1;
                            t_data->read->discarded_by_read_length = 1;
                        }
                    }
                }
             }
      }
    }
    else
    {
        if(vector_flag ) 
        {
            CheckVector(t_data->read); 
       
            if( ( t_data->read->read.length() - (t_data->read->v_end - t_data->read->v_start) ) < minimum_read_length ) 
            {
                t_data->read->discarded = 1;
                t_data->read->discarded_by_vector = 1;
            }
        }
    
        bool adp_found = false;
        found = t_data->read->read.rfind( *t_data->query_str );
        if( found != string::npos ) 
        {
            adp_found = true;
            t_data->read->tru_sec_pos = found;
            t_data->read->b_adapter = adapter_type_R1;
         
            t_data->read->lclip = 1;
            t_data->read->rclip = t_data->read->tru_sec_pos;
            t_data->read->tru_sec_found = 1;
         
            if(t_data->read->rclip < minimum_read_length) 
            { 
               t_data->read->discarded = 1;  
               t_data->read->discarded_by_read_length = 1;
                    
            }
        } else 
        {
            //SSAHA job starts here
            iz_SSAHA *izssaha = new iz_SSAHA();
            AlignResult al_res = izssaha->Find( t_data->read->read , *t_data->query_str );
            AlignScores scores;
            if( al_res.found_flag  ) 
            {
                    scores = CalcScores(al_res.seq_1_al, al_res.seq_2_al, al_res.seq_1_al.length(), 0);
                    if(scores.mismatches <= max_al_mism  ) 
                    {
                        adp_found = true;
                        t_data->read->tru_sec_pos = al_res.pos;
                        t_data->read->b_adapter = adapter_type_R1;
                        t_data->read->lclip = 1;
                        t_data->read->rclip = t_data->read->tru_sec_pos;
                        t_data->read->tru_sec_found = 1;
                        
                        if(t_data->read->rclip < minimum_read_length) 
                        { 
                           t_data->read->discarded = 1;
                           t_data->read->discarded_by_read_length = 1;
                        }
                    }
                }
         
                if(!adp_found) t_data->read->tru_sec_pos = t_data->read->clear_length;
         
                t_data->read->lclip = 1;
                t_data->read->rclip = t_data->read->tru_sec_pos;
         
         
                if(t_data->read->rclip < minimum_read_length) 
                { 
                    t_data->read->discarded = 1;
                    t_data->read->discarded_by_read_length = 1;
                }
         
                delete izssaha;
            }
      }
            
      //Clip points
      if( (qual_trim_flag ) && (vector_flag ) )
      {
          if( (t_data->read->lucy_rclip <= 1) && (t_data->read->lucy_rclip <= 1) ) 
          {
             t_data->read->discarded_by_quality = 1; 
             t_data->read->discarded = 1;
             t_data->read->lclip = 1;
             t_data->read->rclip = 1;
          } 
          else
          {
             int a = t_data->read->v_start;
             int b = t_data->read->read.length() - t_data->read->v_end;
             
             if( a >= b ) //Vector is on the right side
             {
                t_data->read->rclip = min3(t_data->read->tru_sec_pos, t_data->read->lucy_rclip, t_data->read->v_start );
                    
                if(t_data->read->rclip == t_data->read->lucy_rclip)
                {
                   t_data->read->right_trimmed_by_quality = 1;
                }
                else if(t_data->read->rclip == t_data->read->tru_sec_pos)
                {
                   t_data->read->right_trimmed_by_adapter = 1;
                }
                else if(t_data->read->rclip == t_data->read->v_start)
                {
                   t_data->read->right_trimmed_by_vector = 1;
                }
                    
                t_data->read->lclip = 1;
             }
             else //Vector is on the left side
             {
                t_data->read->rclip = min(t_data->read->tru_sec_pos, t_data->read->lucy_rclip );
                if(t_data->read->rclip == t_data->read->lucy_rclip)
                {
                   t_data->read->right_trimmed_by_quality = 1;
                }
                else
                {
                   t_data->read->right_trimmed_by_adapter = 1;
                }
                    
                t_data->read->lclip = max3(t_data->read->lucy_lclip,t_data->read->lclip, t_data->read->v_end );
                if(t_data->read->lclip == t_data->read->lucy_lclip)
                {
                  t_data->read->left_trimmed_by_quality = 1;
                }
                if(t_data->read->lclip == t_data->read->v_end)
                {
                  t_data->read->left_trimmed_by_vector = 1;
                }
             }
                
             if(t_data->read->rclip >= t_data->read->clear_length)
             {
                t_data->read->right_trimmed_by_adapter = 0;
             }
                
             if(t_data->read->rclip <= 1) 
                t_data->read->discarded = 1;
          }
      }
      else if( (qual_trim_flag ) && (!vector_flag) )
      {
            if( (t_data->read->lucy_rclip <= 1) && (t_data->read->lucy_rclip <= 1) ) 
            {
                t_data->read->discarded_by_quality = 1; 
                t_data->read->discarded = 1;// cout << t_data->read->readID;
                t_data->read->lclip = 1;//t_data->read->lucy_lclip;
                t_data->read->rclip = 1;//t_data->read->lucy_rclip;
            } 
            else
            {
                t_data->read->lclip = max(t_data->read->lucy_lclip,t_data->read->lclip);
                if(t_data->read->lclip == t_data->read->lucy_lclip)
                    t_data->read->left_trimmed_by_quality = 1;
                
                t_data->read->rclip = min(t_data->read->tru_sec_pos,t_data->read->lucy_rclip );
                if(t_data->read->rclip == t_data->read->lucy_rclip)
                {
                    t_data->read->right_trimmed_by_quality = 1;
                }
                else if(t_data->read->rclip == t_data->read->tru_sec_pos)
                {
                    t_data->read->right_trimmed_by_adapter = 1;
                }
                
                if(t_data->read->rclip >= t_data->read->clear_length)
                {
                    t_data->read->right_trimmed_by_adapter = 0;
                }
                
                if(t_data->read->rclip <= 1) 
                    t_data->read->discarded = 1;
                
                
            }
      }
      else if( (!qual_trim_flag) && (vector_flag ) )
      {
            
            int a = t_data->read->v_start;
            int b = t_data->read->read.length() - t_data->read->v_end;
                
            if( a >= b ) //Vector is on the right side
            {
                t_data->read->rclip = min(t_data->read->tru_sec_pos, t_data->read->v_start );
                    
                if(t_data->read->rclip == t_data->read->tru_sec_pos)
                {
                   t_data->read->right_trimmed_by_adapter = 1;
                }
                else if(t_data->read->rclip == t_data->read->v_start)
                {
                   t_data->read->right_trimmed_by_vector = 1;
                }
                    
                t_data->read->lclip = 1;
            }
            else //Vector is on the left side
            {
                t_data->read->rclip = t_data->read->tru_sec_pos;
                t_data->read->right_trimmed_by_adapter = 1;
                    
                t_data->read->lclip = max(t_data->read->lclip, t_data->read->v_end );
                if(t_data->read->lclip == t_data->read->v_end)
                {
                    t_data->read->left_trimmed_by_vector = 1;
                }
            }
                
            if(t_data->read->rclip >= t_data->read->clear_length)
            {
                t_data->read->right_trimmed_by_adapter = 0;
            }
                
            if(t_data->read->rclip <= 1) 
            {
                t_data->read->discarded = 1;
            }
                
      }
      else if((!qual_trim_flag) && (!vector_flag))
      {
                
            t_data->read->rclip = t_data->read->tru_sec_pos;
            if(t_data->read->rclip < t_data->read->clear_length)
                t_data->read->right_trimmed_by_adapter = 1;
            if(t_data->read->rclip <= 1) 
                t_data->read->discarded = 1;
            
            t_data->read->lclip = 1;
      }
    
    pthread_exit(NULL);
    return 0;
    
}

//Dynamic Illumina: does not need space to store reads:
void IlluminaDynamicSE()
{
    se_bases_kept = se_bases_discarded = 0;
    se_discard_cnt = 0;
    se_bases_anal = 0;        
    avg_trim_len_se = 0;
    
    long cnt_avg; cnt_avg = 0; //Counters needed for calculating the average trimming length
    long cnt_avg_len; cnt_avg_len = 0;
                 
    double avg_len_se; avg_len_se = 0.0;
    double cnt_right_trim_se, avg_right_trim_len_se; 
    double cnt_left_trim_se, avg_left_trim_len_se;
    
    cnt_right_trim_se = avg_right_trim_len_se = 0;
    cnt_left_trim_se = avg_left_trim_len_se = 0;
    
    long cnt; cnt = 0;
    long se_accept_cnt; se_accept_cnt = 0;
    int ts_adapters; ts_adapters = 0;
    int num_vectors; num_vectors = 0;
    int num_contaminants; num_contaminants = 0;
    int accepted; accepted = 0;
    int discarded; discarded = 0;
//    int discarded_by_quality1, discarded_by_quality2; discarded_by_quality1 = discarded_by_quality2 = 0;
    int discarded_by_contaminant; discarded_by_contaminant = 0;
    int discarded_by_read_length; discarded_by_read_length = 0;
//    int discarded_by_vector1 , discarded_by_vector2; discarded_by_vector1 = discarded_by_vector2 = 0;
    /*Left trims*/
    int left_trimmed_by_quality; left_trimmed_by_quality = 0;
    int left_trimmed_by_vector; left_trimmed_by_vector = 0;
    /*Right trims/discards*/
    int right_trimmed_by_quality; right_trimmed_by_quality = 0;
    int right_trimmed_by_adapter; right_trimmed_by_adapter = 0;
    int right_trimmed_by_vector;  right_trimmed_by_vector = 0;
    
    fstream rep_file, se_output_file;
    rep_file.open(rep_file_name1.c_str(),ios::out);
    rep_file << "ReadID\tlclip\trclip\tTruSeq_pos\tTruSeq_type\tRaw_read_length\tLlucy\tRlucy\tDiscarded\tContaminants\tVectorID\tVecStart\tVecEnd\tVecLen\n";
    
    cout << "Running the Illumina cleaning process..." << endl;
    sum_stat << "Running the Illumina cleaning process..." << endl;
    
    
    
    vector<string> record_block;
    
    
    
    se_output_file.open( se_output_filename.c_str(), ios::out );
    
    for(int jj=0; jj<se_names.size(); ++jj)
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
                        read->readID = record_block[0];
                        read->initial_length = record_block[1].length();
                        read->read = record_block[1];
                        read->quality = line;
                        se_bases_anal += read->read.length();
          
                        //Serial realization - useful for debugging if something does not work as expected
          
                        IlluminaDynRoutine(read, adapter_found, query_string);
                        
                        cnt+=1;
          
                        //Report
                        //Read ID
                        rep_file << read->readID << "\t" << read->lclip << "\t" << read->rclip << "\t" << read->tru_sec_pos << "\t" << read->b_adapter << "\t" << read->initial_length << "\t" << (read->lucy_lclip <= 1 ? 1 : read->lucy_lclip) << "\t" << (read->lucy_rclip <= 1 ? 1 : read->lucy_rclip) << "\t" << read->discarded << "\t" << read->contaminants << "\t" << "NA" << "\n";
                        
                        if( !VectorOnlyFlag )  
                        {
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
                                        avg_right_trim_len_se = GetAvg( avg_right_trim_len_se, cnt_right_trim_se, read->initial_length - read->rclip );
                                    }
                                    if(read->lclip > 0)
                                    {
                                        cnt_left_trim_se += 1;
                                        avg_left_trim_len_se = GetAvg( avg_left_trim_len_se, cnt_left_trim_se, read->lclip );
                                    }
                                    
                                    read->read = read->read.substr(0 , read->rclip );
                                    read->quality = read->quality.substr(0,read->rclip) ; 
                                    read->read = read->read.substr( read->lclip, read->rclip - read->lclip );
                                    read->quality = read->quality.substr( read->lclip, read->rclip - read->lclip );
                 
                                    WriteSEFile(se_output_file, read);
                                    se_accept_cnt+=1;
                                    se_bases_kept += read->read.length();
                                    
                                    if( read->initial_length > read->read.length() )
                                    {
                                        cnt_avg+=1;
                                        avg_trim_len_se = GetAvg( avg_trim_len_se, cnt_avg, read->initial_length - read->read.length() );
                                    }
                 
                                    cnt_avg_len+=1; 
                                    avg_len_se = GetAvg( avg_len_se, cnt_avg_len, read->read.length() );
                                    
                 
                                } 
                                
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
          
                        record_block.clear();
                        read->readID.clear(); 
                        read->quality.clear();
                        read->read.clear();
          
                        delete read;
                        
                }
        }
        in.close();
        
    }
    
    se_output_file.close();
    
    rep_file.close();
    
    
    cout << "====================Summary Statistics====================\n";
    sum_stat << "====================Summary Statistics====================\n";
    cout << "SE reads analyzed: " << cnt << ", Bases:" << se_bases_anal << "\n";
    sum_stat << "SE reads analyzed: " << cnt << ", Bases:" << se_bases_anal << "\n";
    
    cout << "Found ->\n";
    sum_stat << "Found ->\n";
    
    cout << "Adapters: " << ts_adapters << ", " << ( (double)ts_adapters/(double)cnt*100.0) << "%\n";
    sum_stat << "Adapters: " << ts_adapters << ", " << ( (double)ts_adapters/(double)cnt*100.0) << "%\n";
    
    if(vector_flag)
    {
        cout << "# of reads with vector: " << num_vectors<< ", " << ( (double)num_vectors/(double)cnt*100.0) << "%\n";
        sum_stat << "# of reads with vector: " << num_vectors << ", " << ( (double)num_vectors/(double)cnt*100.0) << "%\n";
    }
    
    if(contaminants_flag)
    {
        cout << "# of reads with contaminants: " << num_contaminants << ", " << ( (double)num_contaminants/(double)cnt*100.0) << "%\n";
        sum_stat << "# of reads with contaminants: " << num_contaminants << ", " << ( (double)num_contaminants/(double)cnt*100.0) << "%\n";
    }
    
    if(qual_trim_flag || vector_flag )
    {
        cout << "Reads left trimmed ->" << "\n";
        sum_stat << "Reads left trimmed ->" << "\n";
    }
    
    if(qual_trim_flag)
    {
        cout << "By quality: " <<  left_trimmed_by_quality << "\n";
        sum_stat << "By quality: " <<  left_trimmed_by_quality << "\n";
    }
    
    if(vector_flag)
    {
        cout << "By vector: " <<  left_trimmed_by_vector << "\n";
        sum_stat << "By vector: " <<  left_trimmed_by_vector << "\n";
    }
    
    cout << "Average left trimmed length: " << avg_left_trim_len_se << " bp\n";
    sum_stat << "Average left trimmed length: " << avg_left_trim_len_se << " bp\n";
    
    cout << "Reads right trimmed ->" << "\n";
    sum_stat << "Reads right trimmed ->" << "\n";
    
    cout << "By adapter: " <<  right_trimmed_by_adapter << "\n";
    sum_stat << "By adapter: " <<  right_trimmed_by_adapter << "\n";
    
    if(qual_trim_flag)
    {
        cout << "By quality: " <<  right_trimmed_by_quality << "\n";
        sum_stat << "By quality: " <<  right_trimmed_by_quality << "\n";
    }
    
    if(vector_flag)
    {
        cout << "By vector: " <<  right_trimmed_by_vector << "\n";
        sum_stat << "By vector: " <<  right_trimmed_by_vector << "\n";
    }
    
    cout << "Average right trimmed length: " << avg_right_trim_len_se << " bp\n";
    sum_stat << "Average right trimmed length: " << avg_right_trim_len_se << " bp\n";
    
    cout << "PE1 reads discarded: " << discarded << "\n";
    sum_stat << "PE1 reads discarded: " << discarded << "\n";
    
    if(contaminants_flag)
    {
        cout << "By contaminants: " <<  discarded_by_contaminant << "\n";
        sum_stat << "By contaminant: " <<  discarded_by_contaminant << "\n";
    }
    
    cout << "By read length: " <<  discarded_by_read_length << "\n";
    sum_stat << "By read length: " <<  discarded_by_read_length << "\n";
    
    cout << "-----------------------------------------------------------\n";
    sum_stat << "-----------------------------------------------------------\n";
    
    
    cout << "----------------------Summary for SE----------------------\n";
    sum_stat << "----------------------Summary for SE----------------------\n";
    
    cout << "Pairs kept: " << se_accept_cnt << ", " << ( (double)se_accept_cnt/(double)cnt*100.0) << "%, Bases: " << se_bases_kept << ", " << ( (double)se_bases_kept/(double)(se_bases_anal)*100) <<  "%\n";
    sum_stat << "Pairs kept: " << se_accept_cnt << ", " << ( (double)se_accept_cnt/(double)cnt*100.0) << "%, Bases: " << se_bases_kept << ", " << ( (double)se_bases_kept/(double)(se_bases_anal)*100) <<  "%\n";
    
    cout << "Pairs discarded: " << se_discard_cnt << ", " << ( (double)se_discard_cnt/(double)cnt*100.0) << "%, Bases: " << se_bases_discarded << ", " << ( (double)se_bases_discarded/(double)(se_bases_anal)*100) <<  "%\n";
    sum_stat << "Pairs discarded: " << se_discard_cnt << ", " << ( (double)se_discard_cnt/(double)cnt*100.0) << "%, Bases: " << se_bases_discarded << ", " << ( (double)pe_bases_discarded/(double)(se_bases_anal)*100) <<  "%\n";
    
    cout << "Single Reads PE1 kept: " << se_accept_cnt << ", Bases: " << se_bases_kept << "\n";
    sum_stat << "Single Reads PE1 kept: " << se_accept_cnt << ", Bases: " << se_bases_kept << "\n";
    
    cout << "Average trimmed length PE1: " << avg_trim_len_se << " bp\n";
    sum_stat << "Average trimmed length PE1: " << avg_trim_len_se << " bp\n";
    
    cout << "Average trimmed length PE2: " << avg_trim_len_se << " bp\n";
    sum_stat << "Average trimmed length PE2: " << avg_trim_len_se << " bp\n";
    
    cout << "Average read length PE1: " << avg_len_se << " bp\n";
    sum_stat << "Average read length PE1: " << avg_len_se << " bp\n";
    
    
    cout << "====================Done cleaning====================\n";  
    sum_stat << "====================Done cleaning====================\n";  
   
    
}
