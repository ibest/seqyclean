/* 
 * File:   MakeDir.h
 * Author: ilya
 *
 * Created on 14 Январь 2013 г., 18:17
 */

#ifndef MAKEDIR_H
#define	MAKEDIR_H

/*
@(#)File:           $RCSfile: mkpath.c,v $
@(#)Version:        $Revision: 1.13 $
@(#)Last changed:   $Date: 2012/07/15 00:40:37 $
@(#)Purpose:        Create all directories in path
@(#)Author:         J Leffler
@(#)Copyright:      (C) JLSS 1990-91,1997-98,2001,2005,2008,2012
*/

/*TABSTOP=4*/

#include "jlss.h"
#include "emalloc.h"

#include <errno.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif /* HAVE_UNISTD_H */
#include <string.h>
#include "sysstat.h"    /* Fix up for Windows - inc mode_t */

typedef struct stat Stat;

#ifndef lint
/* Prevent over-aggressive optimizers from eliminating ID string */
const char jlss_id_mkpath_c[] = "@(#)$Id: mkpath.c,v 1.13 2012/07/15 00:40:37 jleffler Exp $";
#endif /* lint */

static int do_mkdir(const char *path, mode_t mode);
int mkpath(const char *path, mode_t mode);


#endif	/* MAKEDIR_H */

