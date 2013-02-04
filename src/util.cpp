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

int min4( int x, int y, int z, int k ) {
    int min = x;
    
    if( (y<min) && (y>0) ) min = y;
    if( (z<min) && (z>0) ) min = z;
    if( (k<min) && (k>0) ) min = k;
    
    return min;
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

string i2str(int value, char* result, int base) {
        // check that the base if valid
	if (base < 2 || base > 36) { *result = '\0'; return string(result); }
	
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

double GetAvg( double past_avg, long n, int cur_diff )
{
    double avg = (double)( past_avg*(n-1) + cur_diff )/(double)n ;
    
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
        while (epdf = readdir(dpdf))
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

string int2str(int num)
{
    char buffer[256];  // make sure this is big enough!!!
    snprintf(buffer, sizeof(buffer), "%d", num);
    string ans = buffer;
    
    return ans;
}