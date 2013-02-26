#include "sffreader.h"

/* L I C E N S E *************************************************************/

/*
    Copyright (C) 2009, 2010 Indraniel Das <indraniel@gmail.com>
                             and Washington University in St. Louis

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, see <http://www.gnu.org/licenses/>
*/

/* I N C L U D E S ***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "sff.h"

/* D E F I N E S *************************************************************/
#define VERSION "0.8.0"
#define PRG_NAME "sff2fastq"
#define FASTQ_FILENAME_MAX_LENGTH 1024
#define SFF_FILENAME_MAX_LENGTH 1024

/* P R O T O T Y P E S *******************************************************/
void help_message(void);
void version_info(void);
void process_options(int argc, char *argv[]);


/* G L O B A L S *************************************************************/
char fastq_file[FASTQ_FILENAME_MAX_LENGTH] = { '\0' };
char sff_file[SFF_FILENAME_MAX_LENGTH] = { '\0' };

// trimming is enabled by default -- like 454 tools
int  opt_trim  = 1;

/* M A I N *******************************************************************/
/*int main(int argc, char *argv[]) {

    process_options(argc, argv);
    process_sff_to_fastq(sff_file, fastq_file, opt_trim);

    return 0;
}**/

vector<sff_read_header> rheaders;


/* F U N C T I O N S *********************************************************/

void process_options(int argc, char *argv[]) {
    int c;
    int index;
    char *opt_o_value = NULL;

    while( (c = getopt(argc, argv, "hvno:")) != -1 ) {
        switch(c) {
            case 'o':
                opt_o_value = optarg;
                break;
            case 'n':
                opt_trim = 0;    /* disable trimming -- like 454 tools */
                break;
            }
    }

    if ( opt_o_value != NULL ) {
        strncpy(fastq_file, opt_o_value, FASTQ_FILENAME_MAX_LENGTH);
    }

    /* process the remaining command line arguments */
    for (index = optind; index < argc; index++) {
        strncpy(sff_file, argv[index], SFF_FILENAME_MAX_LENGTH);
    }

//    /* just take the first passed in non-getopt argument as the sff file */
//    strncpy(sff_file, argv[optind], SFF_FILENAME_MAX_LENGTH);

    /* ensure that a sff file was at least passed in! */
    if ( !strlen(sff_file) ) {
        fprintf(stderr, "%s %s '%s %s' %s\n",
                "[err] Need to specify a sff file!",
                "See", PRG_NAME, "-h", "for usage!");
        exit(1);
    }
}

void process_sff_to_fastq(char *sff_file, int trim_flag) {
    //reads.clear();
    //sff_common_header h;
    sff_read_header rh;
    sff_read_data rd;
    FILE *sff_fp, *fastq_fp;

    if ( (sff_fp = fopen(sff_file, "r")) == NULL ) {
        fprintf(stderr,
                "[err] Could not open file '%s' for reading.\n", sff_file);
        exit(1);
    }
    
    //cout << stderr << endl;
    sff_file_size = get_sff_file_size(sff_fp);
    
    read_sff_common_header(sff_fp, &h);
    verify_sff_common_header((char*)PRG_NAME, (char*)VERSION, &h);


    if ( keep_fastq_orig == true ) {
        vector<string> tmp_rep;
        split_str(string(sff_file), tmp_rep, "//");
    
        if ( ( fastq_fp = fopen( (tmp_rep[tmp_rep.size() - 1].substr(0,tmp_rep[tmp_rep.size() - 1].length()-4) + ".fastq").c_str(), "w") ) == NULL ) {
            fprintf(stderr,
                    "[err] Could not open file '%s' for writing.\n",
                    "");
            exit(1);
        }
    }

    int left_clip = 0, right_clip = 0, nbases = 0;
    char *name;
    char *bases;
    uint8_t *quality;
    //register int i;
    
    int numreads = (int) h.nreads;
    
    for (int i = 0; i < numreads; i++) { //cout << i << " " << numreads << endl;
        read_sff_read_header(sff_fp, &rh);
        read_sff_read_data(sff_fp, &rd, h.flow_len, rh.nbases);
        
        //rheaders.push_back(rh);
        // get clipping points 
        get_clip_values(rh, trim_flag, &left_clip, &right_clip);
        nbases = right_clip - left_clip;

        // create bases string 
        bases = get_read_bases(rd, left_clip, right_clip);

        // create quality array 
        quality = get_read_quality_values(rd, left_clip, right_clip);

        //Create new read
        Read *read = new Read();
        
         
       
        
        //if (rh.name == "GJESF0N01AWIJ4N01BWO5P") 
        //{
        //    printf("!!!");
        //}
        
        //printf("%s\n",read->readID);
        read->initial_length = nbases;
        read->read = string(bases);
        uint8_t quality_char;
        read->quality = (uint8_t*)malloc(sizeof(uint8_t)*nbases);
        for (int j = 0; j < nbases; j++) 
        {
           quality_char = (quality[j] <= 93 ? quality[j] : 93) + 33;
           read->quality[j] = quality_char;
        }
       
        //read->rd = rd;
        read->flowgram = new uint16_t[h.flow_len];
        for(int j=0; j<h.flow_len; j++) {
                read->flowgram[j] = rd.flowgram[j];
                
        }
        
        read->flow_index = (uint8_t*)malloc(sizeof(uint8_t)*nbases);
        for(int j=0; j<nbases; j++) {
                read->flow_index[j] = rd.flow_index[j];
                
        }
        
        read->roche_left_clip = (int) max(1, max(rh.clip_qual_left, rh.clip_adapter_left)) - 1;
        read->roche_right_clip = (int) min( (rh.clip_qual_right    == 0 ? rh.nbases : rh.clip_qual_right   ), (rh.clip_adapter_right == 0 ? rh.nbases : rh.clip_adapter_right) );
        
        reads.push_back(read);
        
        
        string tstr = string(rh.name) + " " + string(itoa(rh.clip_adapter_left,new char[5],10)) +  " " + string(itoa(rh.clip_adapter_right,new char[5],10))+  " " + string(itoa(rh.clip_qual_left,new char[5],10))  +   " " + string(itoa(rh.clip_qual_right,new char[5],10)) + " " + string(itoa(rh.clip_qual_right,new char[5],10)); 
        int t_len = tstr.length();
        
        
        // create read name string 
        int name_length = (int) t_len + 1; // account for NULL termination
        name = (char *) malloc( name_length * sizeof(char) );
        if (!name) {
            fprintf(stderr, "Out of memory! For read name string!\n");
            exit(1);
        }
        memset(name, '\0', (size_t) name_length);
        
        read->readID = (char *) malloc( rh.name_len * sizeof(char) );
        //read->readID = rh.name;
        memcpy( read->readID, rh.name, (size_t) rh.name_len );
        
        //strncpy(name, rh.name, (size_t) rh.name_len);
        strncpy(name, tstr.c_str(), (size_t)t_len);
        
        if ( keep_fastq_orig == true )
            construct_fastq_entry(fastq_fp, name, bases, quality, nbases);
        //printf("%d\n",rh.name_len);
        free(name);
        free(bases);
        free(quality);
        free_sff_read_header(&rh);
        free_sff_read_data(&rd);
        
        
    }
    
    read_manifest(sff_fp);

    //free_sff_common_header(&h);
    if ( keep_fastq_orig == true )
        fclose(fastq_fp);
    
    fclose(sff_fp);
}

void construct_fastq_entry(FILE *fp,
                           char *name,
                           char *bases,
                           uint8_t *quality,
                           int nbases) {
    register int j;
    uint8_t quality_char;

    /* print out the name/sequence blocks */
    fprintf(fp, "@%s\n%s\n+\n", name, bases, name);

    /* print out quality values (as characters)
     * formula taken from http://maq.sourceforge.net/fastq.shtml
     */
    for (j = 0; j < nbases; j++) 
    {
        quality_char = (quality[j] <= 93 ? quality[j] : 93) + 33;
        fprintf(fp, "%c", (char) quality_char );
    }
    fprintf(fp, "\n");
}

void process_fastq_to_sff(char *sff_file) {
    FILE *sff_fp;

    if ( (sff_fp = fopen(sff_file, "w")) == NULL ) {
        fprintf(stderr,
                "[err] Could not open file '%s' for writing.\n", sff_file);
        exit(1);
    }
    
    h.nreads = htobe32(h.nreads);
    
    for(int kk=0; kk<reads.size(); kk++)
    {
        if(reads[kk]->discarded == 1) 
        {
            //h.nreads -= 1;
            discarded_reads++;
        }
    }
    
    sff_common_header ch;
    /* sff files are in big endian notation so adjust appropriately */
    ch.magic        = be32toh(h.magic);//be32toh(779314790);
    //printf("%d\n",h.index_offset );
    ch.index_offset = be64toh(h.index_offset);//be64toh(15730032);
    ch.index_len    = be32toh(h.index_len);//be32toh(61870);
    ch.nreads       = be32toh(h.nreads);//be32toh(reads.size());
    //cout << h.nreads-discarded_reads << endl;
    ch.header_len   = be16toh(h.header_len);//be16toh(1640);
    ch.key_len      = be16toh(h.key_len);//be16toh(4);
    ch.flow =  h.flow;//ff;
    ch.flow_len     = be16toh(h.flow_len);//be16toh(strlen(h.flow));
    ch.flowgram_format = h.flowgram_format;//0x01;
    ch.key = h.key;//"GACT";
    char v[4] = {0x00,0x00,0x00,0x01};
    //ch.version = ;
    memcpy(ch.version,v,4);
    write_sff_common_header(sff_fp, &ch);
    
    int numreads = reads.size();
    for (int i = 0; i < numreads; i++) 
    {
        
        sff_read_header readHeader;
        readHeader.nbases = reads[i]->read.length();
        readHeader.name = (char*)malloc(sizeof(char)*strlen(reads[i]->readID));
        memcpy( readHeader.name, reads[i]->readID, (size_t) strlen(reads[i]->readID) );//This line causes a problem on slarti with sff_extract
        //readHeader.name = reads[i]->readID;
        readHeader.name_len = strlen(reads[i]->readID);
        //printf("%d\n",readHeader.name_len);
        
        //readHeader.name_len = strlen(readHeader.name);
    
        //printf("%d\n",readHeader.name_len);
        /*Working with clip points*/
        if(reads[i]->discarded == 1) 
        {
            readHeader.clip_qual_left = 16;
            readHeader.clip_qual_right = 16;
            readHeader.clip_adapter_left = 16;
            readHeader.clip_adapter_right = 16;
            
            //reads[i]->lclip = reads[i]->rclip = 16;
            //reads[i]->rh.clip_qual_left = 16;
            //reads[i]->rh.clip_qual_right = 16;
            //reads[i]->rh.clip_adapter_left = 16;
            //reads[i]->rh.clip_adapter_right = 16;
            //continue;
        } else 
        {
            //reads[i]->rh.clip_qual_left = reads[i]->lclip;//reads[i]->lucy_lclip;
            //reads[i]->rh.clip_qual_right = reads[i]->rclip;//reads[i]->lucy_rclip;
        
            //reads[i]->rh.clip_adapter_left = reads[i]->lclip;//reads[i]->rlmid.lmid_end;
            //reads[i]->rh.clip_adapter_right = reads[i]->rclip;//reads[i]->rlmid.rmid_start;
            
            readHeader.clip_qual_left = reads[i]->lclip;
            readHeader.clip_qual_right = reads[i]->rclip;
            
            readHeader.clip_adapter_left = reads[i]->lclip;
            readHeader.clip_adapter_right = reads[i]->rclip;
            
        }
        
        readHeader.header_len = sizeof(readHeader.header_len) 
                                    + sizeof(readHeader.name_len)
                                    + sizeof(readHeader.nbases)
                                    + sizeof(readHeader.clip_qual_left)
                                    + sizeof(readHeader.clip_qual_right)
                                    + sizeof(readHeader.clip_adapter_left)
                                    + sizeof(readHeader.clip_adapter_right)
                                    + (sizeof(char) * readHeader.name_len);
    
    //printf("%d\n",readHeader.name_len);
        
        if ( !( readHeader.header_len % 8 == 0) )
        {
                int remainder = 8 - (readHeader.header_len % 8);
                for(int i=0; i< remainder; ++i)
                        readHeader.header_len += 1;
        }
        
        
        
        
        
        sff_read_data readData;
        readData.bases = (char*)reads[i]->read.c_str();
        //printf("%d %d\n", reads[i]->initial_length, reads[i]->read.length());
        readData.flow_index = reads[i]->flow_index;
        readData.flowgram = reads[i]->flowgram;
        
        readData.quality = (uint8_t*)malloc(sizeof(uint8_t)*(readHeader.nbases));
        int j=0;
        for(j=0; j<readHeader.nbases; ++j)
        {
            readData.quality[j] = reads[i]->quality[j] - 33;
            
        }
        //readData.quality[j] = '\0';
        
        
        write_sff_read_header(sff_fp, &(readHeader));
        write_sff_read_data(sff_fp,  &(readData), htobe16(ch.flow_len), htobe32(readHeader.nbases));
        //free(rd.bases);
        free(readData.flow_index);
        free(readData.flowgram);
        free(readData.quality);
        //free(reads[i]->flow_index);
        //free(reads[i]->flowgram);
        
        //printf("%d\n",sff_fp);
        //free(&rd);
        //free(&h);
        //free(&ch);
    }
    
    //printf("%d\n",ftell(sff_fp));
    write_manifest(sff_fp);
    
    fclose(sff_fp);
}
