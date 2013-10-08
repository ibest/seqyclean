/****************************************************************************
*
* Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
* Maryland, U.S.A.  All rights reserved.
*
****************************************************************************/
#include "poly.h"
#include "abi.h"

#define A 0
#define T 3

void PolyAT_Trim(Read* read)
{
    int left, right;
    left = right = 0;
    
    if(qual_trim_flag) {
        left = poly_at_left( (char*)read->read.substr( read->lucy_lclip, read->read.length() - read->lucy_lclip ).c_str(), read->lucy_rclip - read->lucy_lclip); 
        if (left) 
        {
                read->poly_T_clip = read->lucy_lclip + left;
                read->poly_T_found = true;
        }
        right = poly_at_right((char*)read->read.substr( 0, read->lucy_rclip).c_str(), read->lucy_rclip - read->lucy_lclip);
        if (right) 
        {
                read->poly_A_clip = read->lucy_rclip - right;
                read->poly_A_found = true;
        }
    } else {
        left = poly_at_left( (char*)read->read.c_str(), read->read.length()); 
        if (left) 
        {
                read->poly_T_clip = left;
                read->poly_T_found = true;
        }
        right = poly_at_right((char*)read->read.c_str(), read->read.length());
        if (right) 
        {
                read->poly_A_clip = read->read.length()- right;
                read->poly_A_found = true;
        }
    }
 }


int poly_at_left(char *seq, int len)
{
  register int i, err, ttt, pos;

  /* find the first 'cdna' number of connected T's in the first
     'crng' bases of the input sequence */
  for (i=pos=ttt=0; i<crng && i<len; i++/*, seq++*/)
    if (seq[i] == 'T'/*abi_code(*seq)==T*/) {
      ttt++;
      if (ttt>=cdna)
	break;
    } else
      ttt=0;
  
  if (i>=crng || i>=len) /* found nothing within 'crng', return nil */
    return pos;
  
  if (keep) /* keep the poly-T tag for identification purpose */
    return i-ttt+1;

  /* extend span of poly-T within 'cerr' error tolerance */
  for (i++, /*seq++,*/ err=0; i<len; i++/*, seq++*/)
    if (/*abi_code(*seq)==T*/seq[i]=='T') {
      ttt++;
      if (ttt>=cdna)
	err=0;
    } else {
      if (err<=0)
	pos=i;
      err++; ttt=0;
      if (err>c_err)
	return pos;
    }
  return len;
}

int poly_at_right(char *seq, int len)
{
  register int i, err, aaa, pos;
  aaa = 0;
  pos = 0;
  int seq_len = strlen(seq);
  
  /* find the last 'cdna' number of connected A's in the last
     'crng' bases of the input sequence */
  for (i=0; i<crng && i<len; i++/*, seq--*/)
    if (/*abi_code(*seq)*/seq[seq_len-i]=='A') {
      aaa++;
      if (aaa>=cdna)
	break;
    } else
      aaa=0;
  
  if (i>=crng || i>=len) /* found nothing within 'crng', return nil */
    return pos;
  
  if (keep) /* keep the poly-A tag for identification purpose */
    return i-aaa+1;

  /* extend span of poly-A within 'cerr' error tolerance */
  for (i++,/* seq--,*/ err=0; i<len; i++/*, seq--*/) 
    if (/*abi_code(*seq)*/seq[seq_len-i]=='A') {
      aaa++; 
      if (aaa>=cdna)
	err=0;
    } else {
      if (err<=0)
	pos=i;
      err++; aaa=0;
      if (err>c_err)
	return pos;
    }

  return len;
}
