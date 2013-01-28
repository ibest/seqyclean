/****************************************************************************
*
* Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
* Maryland, U.S.A.  All rights reserved.
*
****************************************************************************/

/**************************************************************************/
/*
* qual_trim.c - Quality trimming algorithm which determines a "clean"
*    range for a sequence based on its phred quality values.
*
* Written by Michael Holmes, 4/13/99.
*
* 1/31/2000, Michael Holmes -- Added default_windows().
* 2/15/2000, Michael Holmes -- Added bracket_clean_range() and 
*    set_bracket().
* 2/18/2000, Michael Holmes -- Minor change to test-mode output only.
*
* Copyright (C) 1999, 2000, The Institute for Genomic Research.  All rights
* reserved.
*
* Note: phred quality values are based on the log (base 10) of the
* probability that the corresponding base call is in error:
*
*    Q = -10 log(P_error)
*
* Thanks to Granger Sutton for valuable suggestions.
*
* The main functions of interest to external code (such as lucy) are
* default_windows() and quality_trim().  The grim() function emulates the
* interface that was used by the old quality trimmimg algorithm (grim),
* but it is not called by lucy.  The main() function is compiled only if
* TEST_THIS_CODE is #defined.
*
* The calling hierarchy is:
*
* main()
*     grim()
*         default_windows()
*         quality_trim()
*             bracket_clean_range()
*             multi_window_trim()
*                 window_trim()
*                 multi_window_trim()    [recursive call]
*                 average_error_trim()
*
*/

#ifndef QUALTRIM_H
#define	QUALTRIM_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <string>
#include "ascii.h"
#include "Read.h"
/* the highest quality value for which we have computed the */
/* corresponding probability of error */
#define MAX_QUALITY 99

/* Note: the max number of windows is 20 */
#define MAX_NUMBER_OF_WINDOWS 20

/* default for maximum average probability of error */
#define DEFAULT_ERROR_THRESHOLD 0.025

/* default for maximum probability of error at each end of clean range */
#define DEFAULT_END_LIMIT 0.02

/* default window size for terminal windows */
#define DEFAULT_BRACKET_WINDOW 10

/* default max average error in terminal windows */
#define DEFAULT_BRACKET_ERROR 0.02

#define TRUE 1
#define FALSE 0

#define max(a,b) (((a) > (b)) ? (a): (b))
#define min(a,b) (((a) < (b)) ? (a): (b))

/* global variables and arrays: */
using namespace std;
/* probability of base call error for phred quality values */
/* from 0 to 99 (indexed by quality) */
double perr[] =
{
	1.000000e+00,
	7.943282e-01,
	6.309573e-01,
	5.011872e-01,
	3.981072e-01,
	3.162278e-01,
	2.511886e-01,
	1.995262e-01,
	1.584893e-01,
	1.258925e-01,
	1.000000e-01,
	7.943282e-02,
	6.309573e-02,
	5.011872e-02,
	3.981072e-02,
	3.162278e-02,
	2.511886e-02,
	1.995262e-02,
	1.584893e-02,
	1.258925e-02,
	1.000000e-02,
	7.943282e-03,
	6.309573e-03,
	5.011872e-03,
	3.981072e-03,
	3.162278e-03,
	2.511886e-03,
	1.995262e-03,
	1.584893e-03,
	1.258925e-03,
	1.000000e-03,
	7.943282e-04,
	6.309573e-04,
	5.011872e-04,
	3.981072e-04,
	3.162278e-04,
	2.511886e-04,
	1.995262e-04,
	1.584893e-04,
	1.258925e-04,
	1.000000e-04,
	7.943282e-05,
	6.309573e-05,
	5.011872e-05,
	3.981072e-05,
	3.162278e-05,
	2.511886e-05,
	1.995262e-05,
	1.584893e-05,
	1.258925e-05,
	1.000000e-05,
	7.943282e-06,
	6.309573e-06,
	5.011872e-06,
	3.981072e-06,
	3.162278e-06,
	2.511886e-06,
	1.995262e-06,
	1.584893e-06,
	1.258925e-06,
	1.000000e-06,
	7.943282e-07,
	6.309573e-07,
	5.011872e-07,
	3.981072e-07,
	3.162278e-07,
	2.511886e-07,
	1.995262e-07,
	1.584893e-07,
	1.258925e-07,
	1.000000e-07,
	7.943282e-08,
	6.309573e-08,
	5.011872e-08,
	3.981072e-08,
	3.162278e-08,
	2.511886e-08,
	1.995262e-08,
	1.584893e-08,
	1.258925e-08,
	1.000000e-08,
	7.943282e-09,
	6.309573e-09,
	5.011872e-09,
	3.981072e-09,
	3.162278e-09,
	2.511886e-09,
	1.995262e-09,
	1.584893e-09,
	1.258925e-09,
	1.000000e-09,
	7.943282e-10,
	6.309573e-10,
	5.011872e-10,
	3.981072e-10,
	3.162278e-10,
	2.511886e-10,
	1.995262e-10,
	1.584893e-10,
	1.258925e-10
};

/* number of windows for window trimming */
int num_windows = 3;

/* the window sizes, and max allowed average probability of error in window */
/* Note: the max number of windows is 20 */
int windows[MAX_NUMBER_OF_WINDOWS];
double err_limits[MAX_NUMBER_OF_WINDOWS];

/* the maximum allowed average probability of error over the entire clean range */
double max_avg_error = DEFAULT_ERROR_THRESHOLD;

/* the maximum allowed probability of error for the 2 bases at the left and right */
/* ends of the clean range */
double end_limit = DEFAULT_END_LIMIT;

/* the terminal window size and allowable error used to do an initial */
/* trim of low-quality base calls from each end of sequence */
double bracket_error = DEFAULT_BRACKET_ERROR;
int bracket_window = DEFAULT_BRACKET_WINDOW;

/* these globals are necessary so that my "grim" function can */
/* simulate the old "grim" function */
int *conf_val_raw; 

extern int minimum_read_length;

void window_trim(
	double *prob_err,	/* array of phred error probabilities */
	int length, 		/* number of probabilities (length of sequence) */
	double err_limit,	/* maximum allowed average probability of error in window */ 
	int window,		/* window size */ 
	int min_frag_length,	/* minimum acceptable fragment size */
	int *num_ranges,	/* number of separate clean ranges found */
	int *range_start,	/* caller-allocated array for storage of range */
				/*   start values (0..length - 1) */
	int *range_end		/* caller-allocated array for storage of range */
				/*   end values (0..length - 1) */
	) ;

/**************************************************************************/
/*
* Finds the largest subsequence of a sequence whose average probability
* of error does not exceed a specified maximum.
*
* Globals used:
*	double max_avg_error,	(maximum allowed average probability of error)
*	double end_limit,	(maximum allowed probability of error for bases)
*				(  at each end of clean range)
*
*/
void average_error_trim(
	double *prob_err,	/* array of phred error probabilities */
	int length, 		/* number of quality values (length of sequence) */
	int min_frag_length,	/* minimum acceptable clean fragment length */
	int *cln_left,		/* base index (0..length-1) of start of clean range */ 
	int *cln_right		/* base index (0..length-1) of end of clean range */
	) ;

/**************************************************************************/
/*
* Determines a single clean range, after recursively considering all of
* the specified window sizes.  Calls average_error_trim to trim the
* results of the final window based on overall average probability of
* error.
*
* Globals used:
*	double max_avg_error,	(maximum allowed average probability of error)
*	double end_limit,	(maximum allowed probability of error for bases)
*				(  at each end of clean range)
*
*/
void multi_window_trim(
        double *prob_err,       /* array of phred error probabilities */
        int length,             /* number of quality values (length of sequence) */
	int num_windows,	/* number of windows */
	int *windows,		/* array of window sizes (largest to smallest) */
	double *err_limits,	/* array of maximum allowed average probability */
				/*   for each window size */
	int min_frag_length,	/* the minimum acceptable length of a fragment */
	int *cln_left,          /* base index (0..length-1) of start of clean range */
	int *cln_right          /* base index (0..length-1) of end of clean range */
	);

/**************************************************************************/
/*
* Sets the window size and average probability of error allowed for the
* terminal windows that "bracket" the candidate clean range.
*
* Globals used:
*	int bracket_window,	(size of terminal window)
*	double bracket_error	(allowable average prob. error in window)
*
*/
void set_bracket(int window_size, double max_error);

void bracket_clean_range(
        double *prob_err,       /* array of phred error probabilities */
        int length,             /* number of quality values (length of sequence) */
	int window_len,		/* size of window */
	double max_err,		/* maximum allowable average error in window */
	int min_frag_length,	/* the minimum acceptable length of a fragment */
	int *cln_left,          /* base index (0..length-1) of start of clean range */
	int *cln_right          /* base index (0..length-1) of end of clean range */
	) ;

void quality_trim(
	int *quality, 		/* array of phred quality values */
	int length, 		/* number of quality values (length of sequence) */
	int min_frag_length,	/* minimum acceptable clean fragment length */
	int *cln_left,		/* base index (0..length-1) of start of clean range */ 
	int *cln_right		/* base index (0..length-1) of end of clean range */
	) ;

void default_windows(void);

void grim(int length, int *left, int *right);

int QualTrim( Read* read, double max_avg_err, double end_lim );


#endif	/* QUALTRIM_H */

