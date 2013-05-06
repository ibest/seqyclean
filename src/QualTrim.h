#ifndef QUALTRIM_H
#define	QUALTRIM_H

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

//#define max(a,b) (((a) > (b)) ? (a): (b))
//#define min(a,b) (((a) < (b)) ? (a): (b))



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
int QualTrimIllumina( Read* read, double max_avg_err, double end_lim );

extern int window0;
extern int window1;

extern int phred_coeff_illumina;
extern bool old_style_illumina_flag;


#endif	/* QUALTRIM_H */

