#ifndef POLY_H
#define	POLY_H

/* 
 * File:   poly.h
 * Author: ilya
 *
 * Created on 5 Ноябрь 2012 г., 13:27
 */



extern unsigned short cdna, c_err, crng, keep;

int poly_at_left(char *seq, int len);
int poly_at_right(char *seq, int len);

#endif	/* POLY_H */

