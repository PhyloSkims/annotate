/*
 * sumaclust.h
 *
 *  Created on: april 2, 2012
 *      Author: mercier
 */


#ifndef SUMACLUST_H_
#define SUMACLUST_H_

typedef struct {
    int32_t			next;
    int32_t			threads_number;
    int*			potential_nexts_list;
	fastaSeqPtr*    db;
	int				n;
	int             normalize;
	int 			reference;
	BOOL			lcsmode;
	BOOL			fast;
	double			threshold;
    BOOL 			stop;
    int				sizeForSeqs;
    int16_t**		addresses;
    int16_t**		iseqs1;
    int16_t**		iseqs2;
    int				seeds_counter;
    double			worstscore;
    double			max_ratio;
    int64_t 		elapsedtime;
} thread_control_t;

#endif /* SUMACLUST_H_ */
