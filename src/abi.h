#ifndef ABI_H
#define	ABI_H

/* 
 * File:   abi.h
 * Author: ilya
 *
 * Created on 4 Ноябрь 2012 г., 12:22
 */




#include <stdio.h>
#include <string.h>
#include <stdlib.h>

//using namespace std;

//#include <ctype.h>


//#ifdef	__cplusplus
//extern "C" {
//#endif




int abi_code(int c);

void prepare_abi_mask();
void abi_align(char *a, int width, char *b, int height, int *start, int *end);



//#ifdef	__cplusplus
//}
//#endif





#endif	/* ABI_H */

