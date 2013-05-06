/* 
 * File:   TrimLucy.h
 * Author: ilya
 *
 * Created on 24 Август 2012 г., 12:40
 */

#ifndef TRIMLUCY_H
#define	TRIMLUCY_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include "timer.h"
#include "util.h"
#include <streambuf>
#include <exception>
#include <pthread.h>
#include "Read.h"
using namespace std;

extern vector<Read*> reads;

void Trim(char* filename);


#endif	/* TRIMLUCY_H */

