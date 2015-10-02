/* ---------------------------------------------------------------- */
/* Copyright (c) Atelier de BioInformatique			    */
/* @file: string.c						    */
/* @desc: some strings utility functions			    */
/*								    */
/* @history:							    */
/* @+	    <Hdp>   : Jan 94 : first version from abistr.c	    */
/* @+	    <Gloup> : Feb 94 : cleaned and speedup		    */
/* ---------------------------------------------------------------- */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libaabi.h"
#include "libbio.h"

/* ---------------------------------------------------- */
/* @Function : 	char *str_duplicate_string		*/
/* Purpose : Make a copy of a string buffer	  	*/
/*	     + extra room				*/
/* ---------------------------------------------------- */
char *str_dup_string(char *str, int extra)
{
	char *dst;

	if (! str) 
	    return str;

	if (! (dst = NEWN(char, (strlen(str) + extra + 1)))) {
	    MemoryError();
	    return NULL;
	}
	
	return strcpy(dst, str);
}

/* ---------------------------------------------------- */
/* @Function : 	char *str_erase_char			*/
/* Purpose : Erase a char from a string 		*/
/* ---------------------------------------------------- */
char *str_erase_char(char *str, int c)
{
	char *s, *se;

        if (! str)
	    return str;
	    
	for (se = s = str ; *s ; s++)
	    if ((*se = *s) != c)
		se++;

	*se = '\000';

	return str;
}

/* ---------------------------------------------------- */
/* @Function : 	char *str_replace_char			*/
/* Purpose : replace all occurences of 'cfrom' to 'cto' */
/* ---------------------------------------------------- */
char *str_replace_char(char *str, int cfrom, int cto)
{
	char *s;

        if (! str)
	    return str;
	    
	for (s = str ; *s ; s++)
	    if (*s == cfrom)
		*s = (char) cto;

	return str;
}

/* ---------------------------------------------------- */
/* @Function : 	char *str_trim_trailing			*/
/* Purpose : Trim trailing spaces from 'str'		*/
/* ---------------------------------------------------- */
char *str_trim_trailing(char *str)
{
	char *s;
	
	if (! str)
	   return str;
		
	s = str + strlen(str);

	for (--s ; (s >= str) && isspace(*s) ; s--)
	    /* nop */ ; 

	*++s = '\000';

	return str;	
}

/* ---------------------------------------------------- */
/* @Function : 	char *str_trim_leading			*/
/* Purpose : Trim leading spaces from 'str'		*/
/* ---------------------------------------------------- */
char *str_trim_leading(char *str)
{
	char *sb, *sn;
	
	if (! (str && isspace(*str)))
	    return str;
		
	for (sb = sn = str ; isspace(*sn) ; sn++) 
	    /* nop */ ;

	while (*sn)
	    *sb++ = *sn++;

	*sb = '\000';
	
	return str;
}

/* ---------------------------------------------------- */
/* @Function : 	char *str_pad_right			*/
/* Purpose : Pad 'str' with char 'padchar' on right 	*/
/*	     to built a string of length 'size'		*/
/* ---------------------------------------------------- */
char *str_pad_right(char *str, long size, int padchar)
{
	long len;
	char *s;

	if (! (str && (len = strlen(str)) < size))	
	    return str;

	s = str + len;

	size -= len;
	
	while (size--)
	    *s++ = (char) padchar;
	
	*s = '\000';

	return str;
}

/* ---------------------------------------------------- */
/* @Function : 	char *str_pad_left			*/
/* Purpose : Pad 'str' with char 'padchar' on left	*/
/*	     to built a string of length 'size'		*/
/* ---------------------------------------------------- */
char *str_pad_left(char *str, long size, int padchar)
{
	long len;
	char *s, *t;

	if (! (str && (len = strlen(str)) < size))	
	    return str;

	s = str + len;
	*(t = str + size) = '\000';

	while (len--)
	    *--t = *--s;

	while (t > str)
	    *--t = (char) padchar;
		
	return str;
}

/* ---------------------------------------------------- */
/* @Function : 	char *str_drop_string			*/
/* Purpose : Delete 'nchars' chars from 'str', starting */
/*	     at 'start'					*/
/* ---------------------------------------------------- */
char *str_drop_string(char *str, long start, long nchars)
{
	long len;
	char *sb, *sn;

	if (! (str && (len = strlen(str)) > start))
	    return str;

	if (len < (start + nchars))
	    nchars = len - start;

	sb = str + start;
	sn = sb + nchars;
	
	while (*sn)
	    *sb++ = *sn++;

	*sb = '\000';
	
	return str;
}

/* ---------------------------------------------------- */
/* @Function : 	char *str_insert_string			*/
/* Purpose : Insert 'src' into 'dst', at position 'pos'	*/
/* ---------------------------------------------------- */
char *str_insert_string(char *dst, char *src, long pos)
{
	long srclen, dstlen;
	char *d, *t;

	if (    (! (src && dst))
	     || ((srclen = strlen(src)) == 0)
	     || ((dstlen = strlen(dst)) <= pos))
	     return dst;

	d = dst + dstlen;
	t = d + srclen;
			
	*t = '\000';

	dstlen -= pos;
	
	while (dstlen--)
	    *--t = *--d;

  	t = src;

	while(srclen--)
	    *d++ = *t++;

	return dst;
}

/* ---------------------------------------------------- */
/* @Function : 	char *str_extract_to_mark		*/
/* Purpose : Put into 'dst' the substring from 'src' 	*/
/*	     starting at 'start' up to the char		*/
/*	     'markchar' (or end of string)		*/
/* ---------------------------------------------------- */
char *str_extract_to_mark(char *dst, char *src, long start, int markchar)
{
	long len;
	char *d;

	if (! (src && dst))
	    return dst;

	if ((len = strlen(src)) < start)
	    start = len;

	src += start;

	for (d = dst ; *src && (*src != markchar) ; src++, d++)
	    *d = *src;

	*d = '\000';

	return dst;
}

/* ---------------------------------------------------- */
/* @Function : 	char *str_extract_string		*/
/* Purpose : Extract substring from string 		*/
/*	     'src', starting at 'from' up to 'to'	*/
/*	     included (or end of string)		*/
/* ---------------------------------------------------- */
char *str_extract_string (char *dst, char *src, long from, long to)
{
	long len;
	char *d, *end;

	if (! (src && dst))
	    return dst;

	len = strlen(src);
	
	if (len < from)
	    from = len;
    
	if (len < to)
	    to = len;
	    
	end = src + to;
	
	src += from;

	for (d = dst ; *src && (src <= end) ; src++, d++) 
	    *d = *src;

	*d = '\000';

	return dst;
}

/* ---------------------------------------------------- */
/* @Function : 	char *str_reverse_string		*/
/* Purpose : reverse symbol order in string 'str' 	*/
/* ---------------------------------------------------- */
char *str_reverse_string(char *str)
{
        char *sb, *se, c;

	if (! str)
	    return str;
	    
        sb = str;
        se = str + strlen(str) - 1;

        while(sb <= se) {
           c    = *sb;
          *sb++ = *se;
          *se-- = c;
        }

        return str;
}

/* ---------------------------------------------------- */
/* @Function : 	char *str_upper_string			*/
/* Purpose : uppercase string 'str' 			*/
/* ---------------------------------------------------- */
#define IS_LOWER(c) (((c) >= 'a') && ((c) <= 'z'))
#define TO_UPPER(c) ((c) - 'a' + 'A')

char *str_upper_string(char *str)
{
	char *s;

        if (! str)
	    return str;
    
	for (s = str ; *s ; s++)
	    if (IS_LOWER(*s))
		*s = TO_UPPER(*s);
		
	return str;    
}

#undef IS_LOWER
#undef TO_UPPER

/* ---------------------------------------------------- */
/* @Function : 	char *str_lower_string			*/
/* Purpose : lowercase string 'str' 			*/
/* ---------------------------------------------------- */
#define IS_UPPER(c) (((c) >= 'A') && ((c) <= 'Z'))
#define TO_LOWER(c) ((c) - 'A' + 'a')

char *str_lower_string(char *str)
{
	char *s;

        if (! str)
	    return str;
    
	for (s = str ; *s ; s++)
	    if (IS_UPPER(*s))
		*s = TO_LOWER(*s);
    
	return str;    
}

#undef IS_UPPER
#undef TO_LOWER
