#include "util.h"


int max4( int x, int y, int z, int k ) {
    int max = x;
    
    if(y>max) max = y;
    if(z>max) max = z;
    if(k>max) max = k;
    
    return max;
}

int max3( int x, int y, int z ) {
    int max = x;
    
    if(y>max) max = y;
    if(z>max) max = z;
    
    
    return max;
}


int min3( int x, int y, int z ) {
    int min = x;
    
    if(y<min)  min = y;
    if(z<min)  min = z;
    
    
    return min;
}

//create function to calculate matrix maximum score
int get_max(int lucy_clip,  int xml_clip,  int my_clip) {
	int max=0;
	if(lucy_clip > xml_clip) {
            if(lucy_clip >my_clip) {
                max = lucy_clip;
            } else {
                max = my_clip;
            }
	} else {
            if(xml_clip > my_clip) {
                max = xml_clip;
            } else {
                max = my_clip;
            }
        }
	
	return max;
}

int get_min(int lucy_clip,  int xml_clip,  int my_clip) {
	int min=0;
	if(lucy_clip < xml_clip) {
            if(lucy_clip < my_clip) {
                min = lucy_clip;
            } else {
                min = my_clip;
            }
	} else {
            if(xml_clip < my_clip) {
                min = xml_clip;
            } else {
                min = my_clip;
            }
        }
	
	return min;
}

//Convert string to upper case
void stoupper(std::string& s)	{
        std::string::iterator i = s.begin();
        std::string::iterator end = s.end();

        while (i != end) {
                *i = std::toupper((unsigned char)*i);
                ++i;
        }
}


char* itoa(int value, char* result, int base) {
        // check that the base if valid
	if (base < 2 || base > 36) { *result = '\0'; return result; }
	
        char* ptr = result, *ptr1 = result, tmp_char;
	int tmp_value;
	
	do {
                tmp_value = value;
		value /= base;
		*ptr++ = "zyxwvutsrqponmlkjihgfedcba9876543210123456789abcdefghijklmnopqrstuvwxyz" [35 + (tmp_value - value * base)];
	} while ( value );
	
	// Apply negative sign
	if (tmp_value < 0) *ptr++ = '-';
	*ptr-- = '\0';
	
        while(ptr1 < ptr) {
                tmp_char = *ptr;
		*ptr--= *ptr1;
		*ptr1++ = tmp_char;
	}
	
        return result;
}

std::string i2str(long long value, char* result, int base) {
        // check that the base if valid
    //std::string result;
    
    if (base < 2 || base > 36) { *result = '\0'; return string(result); }
	
        char* ptr = result, *ptr1 = result, tmp_char;
	long long tmp_value;
	
	do {
                tmp_value = value;
		value /= base;
		*ptr++ = "zyxwvutsrqponmlkjihgfedcba9876543210123456789abcdefghijklmnopqrstuvwxyz" [35 + (tmp_value - value * base)];
	} while ( value );
	
	// Apply negative sign
	if (tmp_value < 0) *ptr++ = '-';
	*ptr-- = '\0';
	
        while(ptr1 < ptr) {
                tmp_char = *ptr;
		*ptr--= *ptr1;
		*ptr1++ = tmp_char;
	}
	
        return string(result);
}

string MakeSeqComplement(string init_str) {
    
    for(int i = 0; i<(int)init_str.length(); i++) {
        if(init_str[i] == 'A') {
            init_str[i] = 'T';
            continue;
        }
        if(init_str[i] == 'T') {
            init_str[i] = 'A';
            continue;
        }
        if(init_str[i] == 'G') {
            init_str[i] = 'C';
            continue;
        }
        if(init_str[i] == 'C') {
            init_str[i] = 'G';
            continue;
        }
        
    }
    
    return init_str;
}

string MakeRevComplement(string init_str) {
    
    reverse(init_str.begin(), init_str.end());
    
    for(unsigned int i = 0; i<init_str.length(); i++) {
        if(init_str[i] == 'A') {
            init_str[i] = 'T';
            continue;
        } else if(init_str[i] == 'T') {
            init_str[i] = 'A';
            continue;
        } else if(init_str[i] == 'G') {
            init_str[i] = 'C';
            continue;
        } else if(init_str[i] == 'C') {
            init_str[i] = 'G';
            continue;
        }
    }
    
    return init_str;
}

//Split the string based on tokens:
void split_str(const string& str, vector<string>& tokens, const string& delimiters = ",")
{
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    string::size_type pos = str.find_first_of(delimiters, lastPos);
    
    while (string::npos != pos || string::npos != lastPos)
    {
       // cout << str.substr(lastPos, pos - lastPos) << endl;
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastPos);
    }
}


void TrimNs(string &read) {
    //removing un-needed n-s from line:
    for( int i = read.length()-1; i >=0; i-- ) {
        if( (read[i] =='n') || (read[i] =='N')) 
        {
            
            if ( (int)read.length() < 10 ) break;
            
            read = read.substr( i,read.length() - i );
            
        } else 
        {
            break;
        }
    }
}

void TrimNs2(string &read) {
    //removing un-needed n-s from line:
    int l_length = read.length();
    int tt=l_length - 1;
    for(int rr=tt; rr>=0; rr--) {
        if((read[rr]=='n') || (read[rr]=='N')) {
            read = read.substr(0,l_length-2);
            
            l_length--;
        } else {
            break;
        }
    }
}

string GenNs(int num, char* letter) {
    string nstr;
    
    for(int i=0; i<num; i++) {
        nstr += letter;
    }
    
    return nstr;
}

bool exists(char* filePath)
{
	//The variable that holds the file information
	struct stat fileAtt; //the type stat and function stat have exactly the same names, so to refer the type, we put struct before it to indicate it is an structure.

	//Use the stat function to get the information
	if (stat(filePath, &fileAtt) != 0) //start will be 0 when it succeeds
                return false; //So on non-zero, throw an exception
	
	//S_ISREG is a macro to check if the filepath referers to a file. If you don't know what a macro is, it's ok, you can use S_ISREG as any other function, it 'returns' a bool.
	return true;//S_ISREG(fileAtt.st_mode);
}

double GetAvg( double past_avg, long n, int cur_diff, int first_avg )
{
    //double avg = (double)( past_avg*(n-1) + cur_diff )/(double)n ;
    double avg = (double)( past_avg - first_avg/(double)n + cur_diff/(double)n ) ;
    
    return avg;
}

int MakeDirectory(string path_to_create)
{
    vector<string> path;
    split_str(path_to_create, path, "//");
    string pstr = "";
    if(path_to_create[0] == '/') pstr += '/';
    int flag = 0;
    for(int j=0; j<(int)path.size(); j++)
    {
        pstr += path[j]+"/";
        if (mkdir( (char*)(pstr).c_str(), S_IRWXU|S_IRGRP|S_IXGRP) != 0)
        {
           //perror("mkdir() error");
           flag = -1;
        } else
        {
           flag = 0;
        }
        
    }
    pstr.clear();
    path.clear();
    
    return flag;
}

/* Returns a list of directories (except the ones that begin with a dot) */

void GetDirectories(std::vector<string> &out, char *directory)
{
    DIR *dpdf;
    struct dirent *epdf;

    dpdf = opendir(directory);
    if (dpdf != NULL)
    {
        while ( (epdf = readdir(dpdf) ))
        {
            out.push_back(string(epdf->d_name));    
            printf("Filename: %s\n",epdf->d_name);
        }
    }
    
} // GetFilesInDirectory

vector< vector<string> > GetPEfilenames(string prefix1, string prefix2, char *directory)
{
    vector<string> test;
    GetDirectories(test, directory);
    
    map<string, int > pe1_names;
    map<string, int > pe2_names;
    
    size_t found;
    for(int i=0; i<(int)test.size(); ++i)
    {
        found = test[i].rfind(prefix1);
        if( found != string::npos )
        {
            pe1_names.insert( std::pair<string, int >(test[i], 1) );
        }
        
        found = test[i].rfind(prefix2);
        if( found != string::npos )
        {
            pe2_names.insert( std::pair<string, int >(test[i], 1) );
        }
    }
    
    test.clear();
    vector< vector<string> > groups;
    map<string, int>::iterator iter;
    for(iter = pe1_names.begin(); iter != pe1_names.end(); ++iter)
    {
        test.push_back( (*iter).first );
        groups.push_back( test );
        //cout << (*iter).first << endl;
        test.clear();
    }
    cout << endl;
    int kk = 0;
    for(iter = pe2_names.begin(); iter != pe2_names.end(); ++iter)
    {
        groups[kk].push_back( (*iter).first );
        //cout << (*iter).first << endl;
        kk++;
    }
    
    for(int i=0; i<(int)groups.size(); ++i)
    {
        cout << groups[i][0] << " " << groups[i][1] << endl;
    }
    
    return groups;
}

string double2str(double num)
{
    char buffer[512];  // make sure this is big enough!!!
    snprintf(buffer, sizeof(buffer), "%g", num);
    
    return string(buffer);
}

std::string longlong2str(long long num)
{
    char buffer[256];  // make sure this is big enough!!!
    snprintf(buffer, sizeof(buffer), "%llu", num);
    std::string ans = buffer;
    
    return ans;
}

std::string ulonglong2str(unsigned long long num)
{
    char buffer[256];  // make sure this is big enough!!!
    snprintf(buffer, sizeof(buffer), "%llu", num);
    std::string ans = buffer;
    
    return ans;
}


std::string int2str(int num)
{
    char buffer[256];  // make sure this is big enough!!!
    snprintf(buffer, sizeof(buffer), "%d", num);
    std::string ans = buffer;
    
    return ans;
}

short GetFormat(char* filename) {
    /*This function determines format of input file*/

    int size = 1024, pos;
    int c;
    char *buffer = (char *)malloc(size);
    short format = 3; // 0 - Fasta, 1 - Fastq, 2 - Sff, 3 - Unknown
    
    // At first, try to open a file in text mode:
    FILE *f = fopen(filename, "r");
    if (f) {
        do { // Read one lines from file
            pos = 0;
            do { // Read one line
                c = fgetc(f);
                if (c != EOF) buffer[pos++] = (char)c;
                if(pos >= size-1) { // Increase buffer length - leave room for 0
                    size *=2;
                    buffer = (char*) realloc(buffer, size);
                }
            } while(c != EOF && c != '\n');
            buffer[pos] = 0;
            // Line is now in buffer!
            // Now let's hande line:
            if (buffer[0] == '>') { // Fasta format
                format = 0;
            } else if (buffer[0] == '@') { // Fastq format
                format = 1;
            } else {
                format = 3;// Unknown file format, trying to open it in a binary mode:
            }
            break;
        } while(c != EOF);
        fclose(f);
    }
    free(buffer);
    
    if (format == 3) { // Unknown file format, try to open it in a binary mode:
        sff_c_header ch;
        FILE *sff_fp;

        if ( (sff_fp = fopen(filename, "r")) == NULL ) {
            fprintf(stderr,
                "[err] Could not open file '%s' for reading.\n", filename);
          }
    
        get_sff_common_header(sff_fp, &ch);
        short res = check_sff_common_header(&ch);
        
        if (res == 0) {
            format = 2; // Recognized as Sff file!
        }

    }
    
    return format;
}

void  get_sff_common_header(FILE *fp, sff_c_header *h) {
    char *flow;
    char *key;
    int header_size;

    size_t bytes = fread(&(h->magic)          , sizeof(uint32_t), 1, fp);
    bytes = fread(&(h->version)        , sizeof(char)    , 4, fp);
    bytes = fread(&(h->index_offset)   , sizeof(uint64_t), 1, fp);
    bytes = fread(&(h->index_len)      , sizeof(uint32_t), 1, fp);
    bytes = fread(&(h->nreads)         , sizeof(uint32_t), 1, fp);
    bytes = fread(&(h->header_len)     , sizeof(uint16_t), 1, fp);
    bytes = fread(&(h->key_len)        , sizeof(uint16_t), 1, fp);
    bytes = fread(&(h->flow_len)       , sizeof(uint16_t), 1, fp);
    bytes = fread(&(h->flowgram_format), sizeof(uint8_t) , 1, fp);

    /* sff files are in big endian notation so adjust appropriately */
    /* Linux: not in use any more
    h->magic        = htobe32(h->magic);
    h->index_offset = htobe64(h->index_offset);
    h->index_len    = htobe32(h->index_len);
    h->nreads       = htobe32(h->nreads);
    h->header_len   = htobe16(h->header_len);
    h->key_len      = htobe16(h->key_len);
    h->flow_len     = htobe16(h->flow_len);
     */
    
    h->magic        = htonl(h->magic);
    h->index_offset = BE64toNE(h->index_offset);
    h->index_len    = htonl(h->index_len);
    h->nreads       = htonl(h->nreads);
    h->header_len   = htons(h->header_len);
    h->key_len      = htons(h->key_len);
    h->flow_len     = htons(h->flow_len);

    /* finally appropriately allocate and read the flow and key strings */
    flow = (char *) malloc( h->flow_len * sizeof(char) );
    if (!flow) {
        fprintf(stderr,
                "Out of memory! Could not allocate header flow string!\n");
        exit(1);
    }

    key  = (char *) malloc( h->key_len * sizeof(char) );
    if (!key) {
        fprintf(stderr,
                "Out of memory! Could not allocate header key string!\n");
        exit(1);
    }

    bytes = fread(flow, sizeof(char), h->flow_len, fp);
    h->flow = flow;

    bytes = fread(key , sizeof(char), h->key_len , fp);
    
    h->key = key;
    
    /* the common header section should be a multiple of 8-bytes 
       if the header is not, it is zero-byte padded to make it so */

    header_size = sizeof(h->magic)
                  + sizeof( h->version )
                  + sizeof(h->index_offset)
                  + sizeof(h->index_len)
                  + sizeof(h->nreads)
                  + sizeof(h->header_len)
                  + sizeof(h->key_len)
                  + sizeof(h->flow_len)
                  + sizeof(h->flowgram_format)
                  + (sizeof(char) * h->flow_len)
                  + (sizeof(char) * h->key_len);

     
}


short  check_sff_common_header(sff_c_header *h) {
    /* ensure that the magic file type is valid */
    if (h->magic != SFF_MAGIC_NUM) {
        free_sff_c_header(h);
        return -1;
    }
    return 0;
}

void free_sff_c_header(sff_c_header *h) {
    free(h->flow);
    free(h->key);
}

/* convert 64 bit Big Endian integer to Native Endian(means current machine) */
// As far as I know, there is no standard function for 64 bit conversion on OSX
uint64_t BE64toNE(uint64_t bigEndian)
{
    uint64_t littleEndian = ((bigEndian & (0x00000000000000FF)) << 56) |
                            ((bigEndian & (0x000000000000FF00)) << 40) |
                            ((bigEndian & (0x0000000000FF0000)) << 24) |
                            ((bigEndian & (0x00000000FF000000)) << 8) |
                            ((bigEndian & (0x000000FF00000000)) >> 8) |
                            ((bigEndian & (0x0000FF0000000000)) >> 24) |
                            ((bigEndian & (0x00FF000000000000)) >> 40) |
                            ((bigEndian & (0xFF00000000000000)) >> 56);
    return littleEndian;
}