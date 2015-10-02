/* ---------------------------------------------------------------- */
/* Copyright (c) Atelier de BioInformatique			    */
/* @file: libaabi.c						    */
/* @desc: Abi general purpose library				    */
/* @+	  error notification / system functions			    */
/*								    */
/* @history:							    */
/* @+       <Gloup> : Jan 96 : PWG version			    */
/* ---------------------------------------------------------------- */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAS_MALLOC_H
#include HAS_MALLOC_H
#endif
#include <errno.h>

#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>

#ifndef errno
extern int 	errno;
#endif

#include "libaabi.h"

static int  sAbortLevel = FATAL_LEVEL,
	    sLastError  = 0;

/* ====================================================	*/
/* Error notifications					*/
/* ====================================================	*/

/* -------------------------------------------- */
/* simple error notifier			*/
/* -------------------------------------------- */

int abi_error(char *filename, int lineno, char *msg, int level)
{
	(void) fprintf(stderr, "*%s* [%d] in file %s at line %d / ",
                       ((level > WARNING_LEVEL) ? "Error" : "Warning"),
	               errno, filename, lineno);
        if (errno != 0)
    	  perror(msg);
	else 
	  (void) fprintf(stderr,"-\n");
	
	if (level >= sAbortLevel) {
	    (void) fprintf(stderr,"*Abort*\n");
	    exit(level);
	}

	return 0;
}

/* -------------------------------------------- */
/* report last error since clear_error		*/
/* -------------------------------------------- */
int abi_last_error()
{
	return sLastError;
}

/* -------------------------------------------- */
/* clear error flag				*/
/* -------------------------------------------- */
void abi_clear_error()
{
	sLastError = 0;
}

/* -------------------------------------------- */
/* set up the current abort level		*/
/* -------------------------------------------- */

void abi_set_abort_level(int level)
{
	sAbortLevel = level;
}

/* -------------------------------------------- */
/* kind of hour-glass :-)			*/
/* -------------------------------------------- */

void abi_play_rotator()
{
	static char rotator[] = "|/-\\";
	static int  rotator_position = 0;

	rotator_position = (rotator_position + 1) % 4;
	(void) fprintf(stderr,"\r%c", rotator[rotator_position]);
}

/* ====================================================	*/
/* Cpu time sys calls					*/
/* ====================================================	*/

/* -------------------------------------------- */
/* Get(/Reset) User Cpu time			*/
/* -------------------------------------------- */
double abi_user_cpu_time(int reset)
{
	static double sLast = 0;
	
	double  now, ust;
     	struct  rusage rusage;

	(void) getrusage(RUSAGE_SELF, &rusage);

     	now =    (double) rusage.ru_utime.tv_sec 
	      + ((double) rusage.ru_utime.tv_usec / 1000000.);
	      
	ust = now - sLast;   

	if (reset) 
	   sLast = now;
	   
	return ust;
}

/* -------------------------------------------- */
/* Get(/Reset) Sys Cpu time			*/
/* -------------------------------------------- */
double abi_sys_cpu_time(int reset)
{
	static  double sLast = 0;
	
	double  now, ust;
     	struct  rusage rusage;
     	
	(void) getrusage(RUSAGE_SELF, &rusage);

     	now =    (double) rusage.ru_stime.tv_sec 
	      + ((double) rusage.ru_stime.tv_usec / 1000000.);
	      
	ust = now - sLast;   

	if (reset) 
	   sLast = now;
	   
	return ust;
}

/* -------------------------------------------- */
/* Get a Cpu Time string			*/
/* -------------------------------------------- */
char *abi_str_cpu_time(int reset)
{
	static char buffer[256];

	double ust, syt, tot;

	ust = abi_user_cpu_time(reset);
	syt = abi_sys_cpu_time(reset);
	tot = ust + syt;

	(void) sprintf(buffer, "cpu time user: %f sys: %f tot: %f",
				(float) ust, (float) syt, (float) tot);
				
	return buffer;
}

#if 0

/* ====================================================	*/
/* Memory state (debug)					*/
/* ====================================================	*/

#define PP(fmt, val) (void) fprintf(stderr, fmt, val)

void abi_memory_info(char *header)
{
	struct mallinfo info;

	info = mallinfo();

	if (header)
	    PP ("--- Memory State at : %s ---\n", header);
	else
	    PP ("--- %s ---\n", "Memory State");
	 
	PP ("total space in arena         : %d\n", info.arena);
	PP ("number of ordinary blocks    : %d\n", info.ordblks);
	PP ("number of small blocks       : %d\n", info.smblks);
	PP ("space in holding block head. : %d\n", info.hblkhd);
	PP ("number of holding blocks     : %d\n", info.hblks);
	PP ("space in small blocks in use : %d\n", info.usmblks);
	PP ("space in free small blocks   : %d\n", info.fsmblks);
	PP ("space in ord. blocks in use  : %d\n", info.uordblks);
	PP ("space in free ord. blocks    : %d\n", info.fordblks);
	PP ("space penalty if keep option : %d\n", info.keepcost);
}

#undef PP

#endif
