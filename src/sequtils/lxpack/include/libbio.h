/* ---------------------------------------------------------------- */
/* Copyright (c) Atelier de BioInformatique			    */
/* @file: libbio.h						    */
/* @desc: bioseq & strings generic library / include file	    */
/*								    */
/* @history:							    */
/* @+       <Gloup> : Jan 96 : first draft for PWG		    */
/* ---------------------------------------------------------------- */

#ifndef _H_libbio

#define _H_libbio

/* ==================================================== */
/* Constants						*/
/* ==================================================== */

#define DNA_ALPHA	"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
#define C_DNA_ALPHA 	"TVGHEFCDIJMLKNOPQYSAABWXRZtvghefcdijmlknopqysaabwxrz"

#ifndef TICKS_PER_SEC
#define TICKS_PER_SEC	    60
#endif

#define PWG_TIME_NO_RESET   TIME_NO_RESET
#define PWG_TIME_RESET	    TIME_RESET

/* ==================================================== */
/* Macros standards					*/
/* ==================================================== */

#ifndef NEW
#define NEW(typ)		(typ*)malloc(sizeof(typ)) 
#define NEWN(typ, dim) 		(typ*)malloc((size_t)(dim) * sizeof(typ))
#define REALLOC(typ, ptr, dim)	(typ*)realloc((void *) (ptr), (size_t)(dim) * sizeof(typ))
#define FREE(ptr)		free((Ptr) ptr)
#endif

#define PWG_STRDUP(s)	pwg_dup_string((s), 0)

/* ==================================================== */
/* Prototypes of library functions			*/
/* ==================================================== */

					/* string.c 	*/

char *str_dup_string		( char *str, int extra			);
char *str_erase_char 		( char *str, int c				);
char *str_replace_char		( char *str, int cfrom, int cto		);
char *str_trim_trailing 	( char *str					);
char *str_trim_leading 		( char *str					);
char *str_pad_right 		( char *str, long size, int padchar		);
char *str_pad_left 		( char *str, long size, int padchar		);
char *str_drop_string 		( char *str, long start, long nchars		);
char *str_insert_string 	( char *dst, char *src, long pos		);
char *str_extract_string 	( char *dst, char *src, long from, long to	);
char *str_extract_to_mark 	( char *dst, char *src, long start, int markchar );
char *str_reverse_string 	( char *str	 			        );
char *str_upper_string		( char *str					);
char *str_lower_string		( char *str					);

					/* bioseq.c		*/

int  bio_base_complement  	( int c	);
char *bio_seq_complement    	( char *str	);
int  bio_codon_translate	( char *codon, int codid);
char *bio_seq_translate		( char *seq, int codid);

#endif
