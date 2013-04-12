#include "QualTrim.h"

/**************************************************************************/
/*
* Determines some number of candidate clean ranges using a particular
* window size.
*
* Globals used:
*	double end_limit,	(maximum allowed probability of error for bases)
*				(  at each end of clean range)
*
*/

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
	)
{
	register double win_limit, win_err;
	register int lidx, ridx, win_idx, end_idx;
	int i;
	int range_idx;

	/* initialize number of clean ranges */
	*num_ranges = 0;

	/* length of sequence must be at least the window size */
	if (length < window || length < min_frag_length)
		return;

	/* calculate max allowable cumulative error for window */
	win_limit = err_limit * (double)window;

	/* calculate cumulative error for first window */
	for (i = 0, win_err = 0.0; i < window; i++)
		win_err += prob_err[i];

	/* initialize index to storage of candidate clean ranges */
	range_idx = 0;

	/* initialize left, right indices for next clean range */
	lidx = ridx = -1;

	/* initialize window beginning and end indices */
	win_idx = 0;
	end_idx = window - 1;

	/* look at each window in turn */
	while (end_idx < length)
	{
		/* is current window's cumulative error probability OK ? */
		if (win_err <= win_limit)
		{
			/* yes, mark the start of a clean range if we're */
			/* not already in one, and if the first base of window */
			/* has acceptable error probability */
			if (lidx < 0 && prob_err[win_idx] <= end_limit)
			{
				lidx = win_idx;
			}

			/* if we're in a clean range and last base of window */
			/* has acceptable error probability, mark it as the */
			/* provisional end of the still-expanding clean range */
			if (lidx >= 0 && prob_err[end_idx] <= end_limit)
			{
				ridx = end_idx;
			}
		}
		else
		{
			/* no, is this the end of a clean range ? */
			if (lidx >= 0 && ridx >= 0 && (ridx - lidx + 1) >= min_frag_length)
			{
				/* yes, save clean range start and end */
				/* as base indices (0..length-1) */
				range_start[range_idx] = lidx;
				range_end[range_idx] = ridx;

				/* increment clean range index */
				range_idx++;
			}

			/* mark the fact that we're not within a clean range */
			lidx = ridx = -1;
		}

		/* advance the window */
		win_err -= prob_err[win_idx];
		win_idx++;
		end_idx++;
		/* (Note: here's why we allocated one extra cell in prob_err array) */
		win_err += prob_err[end_idx];
	}

	/* were we inside a clean range at end of while loop ? */
	if (lidx >= 0 && ridx >= 0 && (ridx - lidx + 1) >= min_frag_length)
	{
		/* yes, save clean range start and end */
		/* as base indices (0..length-1) */
		range_start[range_idx] = lidx;
		range_end[range_idx] = ridx;

		/* increment clean range index */
		range_idx++;
	}

	/* store number of clean ranges found for caller */
	*num_ranges = range_idx;

}  /* window_trim() */


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
	)
{
	int i, j, frag_length;
	int diag_count;
	int done;
	double *diag, this_err;

	/* initialize caller's clean range */
	*cln_left = *cln_right = 0;

	/* if sequence is too short, we're done */
	if (length < min_frag_length)
		return;

	/* allocate space for cumulative error values on the diagonal */
	diag = (double *)malloc((size_t)(length + 1) * sizeof(double));
	if (diag == NULL)
	{
		fprintf(stderr, "Memory allocation failure in find_largest_clean_range.\n");
		exit(1);
	}

	/* calculate the single diagonal value for corner of matrix */
	/* (corresponding to the full-length sequence) */
	for (i = 0, *diag = 0.0; i < length; i++)
	{
		*diag += prob_err[i];
	}

	frag_length = length;
	done = FALSE;
	while (frag_length >= min_frag_length && !done)
	{
		/* calculate cumulative error of last cell on next diagonal */
		diag_count = length - frag_length + 1;
		diag[diag_count] = diag[diag_count - 1] - prob_err[diag_count - 1];

		/* consider each value on this diagonal */
		for (i = 0, j = frag_length - 1;
			j < length; i++, j++)
		{
			/* calculate average error of bases i..j */
			this_err = diag[i] / (double)frag_length;

			/* is it good enough ? */
			if (this_err <= max_avg_error)
			{
				/* yes, this is our clean range */
				*cln_left = i;
				*cln_right = j;
				done = TRUE;
				break;
			}
			else
			{
				/* for next diagonal, subtract out error of last base in range */
				diag[i] -= prob_err[j];
			}
		}
		
		/* decrement fragment length of next diagonal */
		frag_length--;
	}

	/* make sure the ends of the clean range are OK */
	if (*cln_right > 0)
	{
		while (prob_err[*cln_left] > end_limit)
			(*cln_left)++;
		while (prob_err[*cln_right] > end_limit)
			(*cln_right)--;
	}

	/* check fragment length */
	if (*cln_right - *cln_left + 1 < min_frag_length)
	{
		*cln_left = *cln_right = 0;
	}
		
	/* free allocated array */
	free(diag);

}  /* average_error_trim() */


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
	)
{
	int *range_start, *range_end;
	int num_ranges, left, right, i, max_clean;

	/* initialize caller's clean range */
	*cln_left = *cln_right = 0;

	/* allocate space for start and end of each candidate clean range */
	range_start = (int *)malloc((size_t)(length + 1) * sizeof(int));
	if (range_start == NULL)
	{
		fprintf(stderr, "Memory allocation failure in multi_window_trim.\n");
		exit(1);
	}
	range_end = (int *)malloc((size_t)(length + 1) * sizeof(int));
	if (range_end == NULL)
	{
		fprintf(stderr, "Memory allocation failure in multi_window_trim.\n");
		exit(1);
	}

	/* trim sequence to largest window size */
	window_trim(prob_err, length, err_limits[0], windows[0],
		min_frag_length, &num_ranges, range_start, range_end);

	/* any clean ranges found ? */
	if (num_ranges > 0)
	{
#ifdef TEST_THIS_CODE
		if (num_windows == 3)
			fprintf(stderr, "Ranges:  ");
#endif
		/* yes, any smaller windows to be dealt with ? */
		if (num_windows > 1)
		{
			/* yes, trim each candidate clean range by the smaller windows */
			for (i = 0; i < num_ranges; i++)
			{
				/* recursive call */
				multi_window_trim(
					prob_err + range_start[i],
					range_end[i] - range_start[i] + 1,
					num_windows - 1,
					windows + 1, err_limits + 1,
					min_frag_length,
					&left, &right);

				if (right > 0)
				{
					range_end[i] = right + range_start[i];
					range_start[i] = left + range_start[i];
				}
				else
				{
					range_start[i] = range_end[i] = 0;
				}
			}
		}
		else
		{
			/* no, trim each of the bottom-level ranges based on average quality */
			for (i = 0; i < num_ranges; i++)
			{
				average_error_trim(prob_err + range_start[i],
					range_end[i] - range_start[i] + 1,
					min_frag_length,
					&left, &right);
				
				if (right > 0)
				{
					range_end[i] = right + range_start[i];
					range_start[i] = left + range_start[i];
				}
				else
				{
					range_start[i] = range_end[i] = 0;
				}
			}
		}

#ifdef TEST_THIS_CODE
		if (num_windows == 1)
		{
			for (i = 0; i < num_ranges; i++)
				if (range_end[i] > 0)
					fprintf(stderr, "%d  ", range_end[i] - range_start[i] + 1);
		}
		else if (num_windows == 3)
			fprintf(stderr, "\n");
#endif

		/* find the largest clean range */
		max_clean = -1;
		for (i = 0; i < num_ranges; i++)
		{
			if (range_end[i] > 0 && range_end[i] - range_start[i] > max_clean)
			{
				max_clean = range_end[i] - range_start[i];
				*cln_left = range_start[i];
				*cln_right = range_end[i];
			}
		}
	}

	/* free allocated arrays */
	free(range_start);
	free(range_end);

}  /* multi_window_trim() */
	

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
void set_bracket(int window_size, double max_error)
{
	bracket_window = window_size;
	bracket_error = max_error;

}  /* set_bracket() */


/**************************************************************************/
/*
* Finds the leftmost window, and the rightmost window, which meet the
* specified error criterion.  The purpose of this function is to bracket
* the portion of the sequence that is of decent quality -- i.e. to
* eliminate the really bad quality stuff that is often found on either
* end.  The beginning of the first matching window, and the end of the
* last, bracket the sequence range that will be looked at to find the
* final clean range.
*
* This function is called before the multi_window_trim function, and
* its results are used to limit the range of sequence that will be
* considered by that function.
*
*/
void bracket_clean_range(
        double *prob_err,       /* array of phred error probabilities */
        int length,             /* number of quality values (length of sequence) */
	int window_len,		/* size of window */
	double max_err,		/* maximum allowable average error in window */
	int min_frag_length,	/* the minimum acceptable length of a fragment */
	int *cln_left,          /* base index (0..length-1) of start of clean range */
	int *cln_right          /* base index (0..length-1) of end of clean range */
	)
{
	int left, right;
	int i, j;
	double cum_error, max_cum_error;

	/* initialize caller's clean range */
	*cln_left = *cln_right = 0;

	/* initialize my left and right brackets */
	left = right = -1;

	/* calculate cumulative allowable error in window */
	max_cum_error = (double)window_len * max_err;

	/* calculate cumulative error in first window */
	for (i = 0, cum_error = 0.0; i < window_len; i++)
	{
		cum_error += prob_err[i];
	}

	/* test first window */
	if (cum_error <= max_cum_error)
	{
		/* it passes, set left and right brackets */
		left = 0;
		right = window_len - 1;
	}

	/* test the remaining windows */
	/* i is index to beginning of previous window */
	/* j is index to end of current window */
	for (i = 0, j = window_len; j < length; i++, j++)
	{
		/* calculate cumulative error in this window */
		cum_error = cum_error - prob_err[i] + prob_err[j];

		/* test window */
		if (cum_error <= max_cum_error)
		{
			/* window passes */
			/* leave left bracket alone if it's already set */
			if (left < 0)
			{
				/* set left bracket to start of current window */
				left = i + 1;
			}

			/* set right bracket to end of current window */
			right = j;
		}
	}

	/* did we bracket a sequence of sufficient length ? */
	if (right > 0 && (right - left + 1) >= min_frag_length)
	{
		/* yes, set caller's clean range variables */
		*cln_left = left;
		*cln_right = right;
	}

#ifdef TEST_THIS_CODE
	fprintf(stderr, "Bracket:  %d %d\n", *cln_left, *cln_right);
#endif

}  /* bracket_clean_range() */


/**************************************************************************/
/*
* Determines the clean range.  Calls multi_window_trim, which calls both
* window_trim and average_error_trim.
*
* Globals used:
*	int num_windows,	(number of windows)
*	int *windows,		(array of window sizes (largest to smallest))
*	double *err_limits,	(array of maximum allowed average probability)
*				(  for each window size)
*	double max_avg_error,	(maximum allowed average probability of error)
*	double end_limit,	(maximum allowed probability of error for bases)
*				(  at each end of clean range)
*	int bracket_window,	(size of terminal window)
*	double bracket_error	(allowable average prob. error in window)
*
*/
void quality_trim(
	int *quality, 		/* array of phred quality values */
	int length, 		/* number of quality values (length of sequence) */
	int min_frag_length,	/* minimum acceptable clean fragment length */
	int *cln_left,		/* base index (0..length-1) of start of clean range */ 
	int *cln_right		/* base index (0..length-1) of end of clean range */
	)
{
	int i, q;
	int left, right;
	int min_left, max_right;
	double *prob_err;
	//double err, sum;

	/* initialize caller's clean range */
	*cln_left = *cln_right = 0;

	/* if sequence is too short, we're done */
	if (length < min_frag_length)
		return;

	/* allocate space for error probabilities */
	/* (allocate length + 1 locations for computational convenience) */
	prob_err = (double *)malloc((size_t)(length + 1) * sizeof(double));
	if (prob_err == NULL)
	{
		fprintf(stderr, "Memory allocation failure in quality_trim.\n");
		exit(1);
	}

	/* convert qualities to error probabilities */
	for (i = 0; i < length; i++)
	{
		q = quality[i];
                //cout << q << " ";
		if (q > MAX_QUALITY)
			q = MAX_QUALITY;

		prob_err[i] = perr[q];
	}

	/* set the extra value at end of array */
	prob_err[length] = 0.0;

	/* find beginning and end of reasonable-quality sequence */
	left = right = 0;
	bracket_clean_range(prob_err, length, 
		bracket_window, bracket_error,
		min_frag_length,
		&min_left, &max_right);

	if (max_right > 0)
	{
		/* find largest sequence that matches all our window criteria */
		multi_window_trim(prob_err + min_left, 
			max_right - min_left + 1,
			num_windows, windows, err_limits,
			min_frag_length,
			&left, &right);

		if (right > 0)
		{
			left = left + min_left;
			right = right + min_left;
		}
	}

	if (right > 0 && (right - left + 1) >= min_frag_length)
	{
		/* store clean range for caller */
		*cln_left = left;
		*cln_right = right;
	}

#ifdef TEST_THIS_CODE
	if (right > 0)
	{
		/* calculate the probability of error of the final clean range */
		for (i = left, sum = 0.0; i <= right; i++)
		{
			sum += prob_err[i];
		}

		err = sum / (double)(right - left + 1);
		if (err > max_avg_error)
		{
			/* error -- this shouldn't happen */
			fprintf(stderr, "*** Error: %lf %lf\n", err, max_avg_error);
		}
	}
#endif

	/* free allocated array */
	free(prob_err);

}  /* quality_trim() */


/**************************************************************************/
/*
* Sets up 3 windows for quality trimming, with the average error allowed
* in each window calculated from max_avg_error.
*
* Globals used:
*	int num_windows,	(number of windows)
*	int *windows,		(array of window sizes (largest to smallest))
*	double *err_limits,	(array of maximum allowed average probability)
*				(for each window size)
*	double max_avg_error,	(maximum allowed average probability of error)
*
*/
void default_windows(void)
{
	/* 3 fixed window sizes: 100, 30, 5 */
	num_windows = 2;
	windows[0] = window0;//50;//window0;
	windows[1] = window1;//10;//window1;

	/* error criteria for windows are looser than for full-length */
	/* clean range */
	/*
	err_limits[0] = (5.0 + 95.0 * max_avg_error) / 100.0;
	err_limits[1] = (3.0 + 27.0 * max_avg_error) / 30.0;
	err_limits[2] = (1.0 + 4.0 * max_avg_error) / 5.0;
	*/

	err_limits[0] = 0.08;
	err_limits[1] = 0.3;

}  /* default_windows() */


/**************************************************************************/
/*
* Simulates the functionality of the old grim() function, by calling
* quality_trim().
*
* Globals used:
*	int num_windows,	(number of windows)
*	int *windows,		(array of window sizes (largest to smallest))
*	double *err_limits,	(array of maximum allowed average probability)
				(for each window size)
*	double max_avg_error,	(maximum allowed average probability of error)
*	double end_limit,	(maximum allowed probability of error for bases)
*				(  at each end of clean range)
*
*/
void grim(int length, int *left, int *right)
{
	/* set up windows for quality trimming */
	default_windows();

	/* the full-length clean range must have an average probability */
	/* of error no greater than max_avg_error */
	quality_trim(
		conf_val_raw, 		/* the raw quality values */
		length, 		/* length of sequence */
		minimum_read_length,	 		/* minimum acceptable clean range length */
		left, 			/* pointer to caller's start index */
		right			/* pointer to caller's end index */
		);

}  /* grim() */


//#ifdef TEST_THIS_CODE
/**************************************************************************/
//int main(int argc, char **argv)
int QualTrim( Read* read, double max_avg_err, double end_lim )
{       
   // cout << qual_str << endl;
        //FILE * pFile;
	int quality[10000];
	unsigned int qual_count, i;//length, i;
	int left, right;
	//char line_buff[4096] /*seq_name[80]*/, *sptr;

        max_avg_error = max_avg_err;
        end_limit = end_lim;
        
        i=0;
        for( i=0; i< read->read.length(); i++ ) 
        {
            quality[i] = GetNum(read->quality[i]) - 33;
        }
        
        qual_count = i;
        
        /* trim for quality */
	conf_val_raw = quality;
	grim(qual_count, &left, &right);
        
        if (right > 0)
	{
           left++;
	   right++;
	}

	/* display seq name and clean range */
	if (right - left < /*99*/minimum_read_length)
		left = right = 0;
        
        if(left == 0 && right == 0) 
        {
            read->discarded = 1;
            read->discarded_by_read_length = 1;
            read->lucy_lclip = 1;
            read->lucy_rclip = 1;
        } 
        else 
        {
            left == 0 ? read->lucy_lclip = 0 : read->lucy_lclip = left;
            right == 0 ? read->lucy_rclip = 1 : read->lucy_rclip = right;
        }
	

	return 0;

}


int QualTrimIllumina( Read* read, double max_avg_err, double end_lim )
{       
   // cout << qual_str << endl;
        //FILE * pFile;
	int quality[10000];
	unsigned int qual_count, i;//length, i;
        int left, right;
	//char line_buff[4096] /*seq_name[80]*/, *sptr;

        max_avg_error = max_avg_err;
        end_limit = end_lim;
        
        i=0;
        for( i=0; i< read->read.length(); i++ ) 
        {
            quality[i] = GetNum(read->illumina_quality_string[i]) - phred_coeff_illumina;//33 or 64 depending on new or old style;
            
            if (old_style_illumina_flag == true)
            {
                read->illumina_quality_string[i] = read->illumina_quality_string[i] - phred_coeff_illumina + 33;
                
            }
        }
        
        qual_count = i;
        
        /* trim for quality */
	conf_val_raw = quality;
	grim(qual_count, &left, &right);
        /*
        if (right > 0)
	{
           left++; //Lucy always clips one base.
	   right++;
	}**/

	/* display seq name and clean range */
	if (right - left < /*99*/minimum_read_length)
		left = right = 0;
        
        if(left == 0 && right == 0) 
        {
            read->discarded = 1;
            read->discarded_by_read_length = 1;
            read->lucy_lclip = 1;
            read->lucy_rclip = 1;
        } 
        else 
        {
            left == 0 ? read->lucy_lclip = 0 : read->lucy_lclip = left;
            right == 0 ? read->lucy_rclip = 1 : read->lucy_rclip = right;
        }
	

	return 0;

}  /* main() */

/* main() */
//#endif
