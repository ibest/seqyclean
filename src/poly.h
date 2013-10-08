#ifndef POLY_H
#define	POLY_H

/* 
 * File:   poly.h
 * Author: ilya
 *
 * Created on 5 Ноябрь 2012 г., 13:27
 */
#include "Read.h"


extern unsigned short cdna, c_err, crng, keep;
extern bool qual_trim_flag;

void PolyAT_Trim(Read* read);
int poly_at_left(char *seq, int len);
int poly_at_right(char *seq, int len);

#endif	/* POLY_H */

