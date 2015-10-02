/* ---------------------------------------------------------------- */
/* Copyright (c) Atelier de BioInformatique			    */
/* @file: libaabi.h						    */
/* @desc: Abi general purpose library / include file		    */
/*								    */
/* @history:							    */
/* @+       <Gloup> : Jan 96 : first draft for PWG		    */
/* ---------------------------------------------------------------- */

#ifndef _H_libaabi

#define _H_libaabi

/* ==================================================== */
/* Constants						*/
/* ==================================================== */

#ifndef __FILE__
#define __FILE__ "unknown file"
#endif

#ifndef __LINE__
#define __LINE__ 0
#endif

#define Vrai	1
#define Faux	0

#define TIME_NO_RESET   0
#define TIME_RESET      1

#define WARNING_LEVEL   1
#define FATAL_LEVEL	10
#define NO_ABORT_LEVEL	255

/* ==================================================== */
/* Macros standards					*/
/* ==================================================== */

#ifndef NEW
#define NEW(typ)		(typ*)malloc(sizeof(typ)) 
#define NEWN(typ, dim) 		(typ*)malloc((size_t)(dim) * sizeof(typ))
#define REALLOC(typ, ptr, dim)	(typ*)realloc((void *) (ptr), (size_t)(dim) * sizeof(typ))
#define FREE(ptr)		free(ptr)
#endif

#define Error(message, level)	abi_error(__FILE__, __LINE__, message, level)

#define W_Error(message) Error(message, WARNING_LEVEL)
#define F_Error(message) Error(message, FATAL_LEVEL)

#define MemoryError()	 F_Error("Not enough memory")
#define IOError()	 F_Error("IO Error")

/* ==================================================== */
/* Prototypes of library functions			*/
/* ==================================================== */

					/* libaabi.c	*/

int	abi_error	(char *filename, int lineno, char *msg, int level);
int	abi_last_error	();

void	abi_clear_error	    (),
	abi_set_abort_level (int level),
	abi_play_rotator    ();

double	abi_user_cpu_time   (int reset);
double	abi_sys_cpu_time    (int reset);
char	*abi_str_cpu_time   (int reset);

void	abi_memory_info	    (char *header);

#endif
