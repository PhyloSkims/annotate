/*
 * mtcompare_cumaclust.c
 *
 *  Author: Celine Mercier
 *
 */


#ifdef OMP_SUPPORT
#include <omp.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#include "./sumalibs/libfasta/sequence.h"
#include "./sumalibs/libutils/utilities.h"
#include "./sumalibs/liblcs/upperband.h"
#include "./sumalibs/liblcs/sse_banded_LCS_alignment.h"

#include "sumaclust.h"



static double computeScore(void* c, int32_t seed_number, int32_t i, int thread_number,int32_t maxcount)
{
    thread_control_t 	*control=(thread_control_t*)c;
    fastaSeqPtr*        db   = control->db;
    fastaSeqPtr         db_i = db[i];
    fastaSeqPtr			db_seed_number = db[seed_number];
	int 				LCSmin;
	double	 			score;

	score = control->worstscore;

	if (db_i->count  <= maxcount)
		filters(db_i, db_seed_number,
				control->threshold,
				control->normalize,
				control->reference,
				control->lcsmode,
				&score,
				&LCSmin);

    if (score == -1.0)
		score = alignForSumathings(db_seed_number->sequence, control->iseqs1[thread_number],
								   db_i->sequence,           control->iseqs2[thread_number],
								   db_seed_number->length,   db_i->length,
								   control->normalize,       control->reference,
								   control->lcsmode,         control->addresses[thread_number],
								   control->sizeForSeqs,     LCSmin);

    return score;
}


inline void putSeqInClusterMT(void *c, int32_t center_idx, int32_t seq, double score)
{
	// saves a sequence as belonging to a cluster and its score with the seed

    thread_control_t *control=(thread_control_t*)c;
    fastaSeqPtr*        db   = control->db;
    fastaSeqPtr         pseq = db[seq];


    pseq->center         = db+center_idx;
    pseq->center_index   = center_idx;				// saves cluster
    pseq->score          = score;			    	// saves score with the seed
    pseq->cluster_center = FALSE;
}


int64_t timevaldiff(struct timeval *starttime, struct timeval *finishtime)
{
  int64_t msec;
  msec=(finishtime->tv_sec-starttime->tv_sec)*1000000;
  msec+=(finishtime->tv_usec-starttime->tv_usec)/1000000;
  return msec;
}

void computeOneSeed(void* c)
{
    thread_control_t *control=(thread_control_t*)c;
    BOOL 		  	  found;
    int32_t 		  seed_number;
    int32_t			  nextone=control->n;
	int64_t elapsedtime;
	struct timeval current;
	struct timeval start;

    seed_number = control->next;
    found = FALSE;

    //printf("\n seed = %d, n = %d", seed_number, control->n);

	#ifdef OMP_SUPPORT
    omp_set_num_threads(control->threads_number);
	#endif

	gettimeofday(&start,NULL);

	#ifdef OMP_SUPPORT
	#pragma omp parallel default(none) \
						 firstprivate(found) \
						 firstprivate(seed_number) \
						 firstprivate(control) \
    					 shared(nextone)
	#endif

	{
	    int32_t 		i;
	    double  		score;
	    int32_t 		current_seed;
		#ifdef OMP_SUPPORT
	    int				thread_id=omp_get_thread_num();
		#else
	    int				thread_id=0;
		#endif
	    int     		nseq = control->n;
	    BOOL     		fast = control->fast;
	    BOOL			lcsmode = control->lcsmode;
		int             normalize = control->normalize;
		double			threshold = control->threshold;
	    BOOL            first = TRUE;
	    BOOL			not_already_in_a_cluster;
	    BOOL			threshold_bad;
	    int32_t			priv_nextone=control->n;
	    int32_t         maxcount = (double)(control->db[seed_number]->count) * control->max_ratio;

	    fastaSeqPtr* 	db = control->db;

		#ifdef OMP_SUPPORT
		#pragma omp  for schedule(dynamic,10)
		#endif

    	for (i=seed_number+1;    \
    		 i < nseq;           \
    		 i++)
		{

			current_seed = db[i]->center_index;
			not_already_in_a_cluster = current_seed == i; // At the beginning all the sequences are their own center

    		if ((! fast) || not_already_in_a_cluster)
    		{
				score = computeScore((void*)control, seed_number, i, thread_id,maxcount);		// computes LCS score or 0 if k-mer filter not passed

				if (lcsmode || normalize)
					threshold_bad = (score < threshold);
				else
					threshold_bad = (score > threshold);

    			if (threshold_bad) // similarity under threshold
				{
					if (!found && not_already_in_a_cluster && (i < priv_nextone))
					{
						priv_nextone=i;					// saves potential next seed
//						*potential_nexts_list = i;

						found = TRUE;		        // saves the fact that a next seed
													// has been found for this thread
					}
				}
				else if (not_already_in_a_cluster ||    \
				         ((! fast) &&       \
                         (db[i]->score < score)))
				{ // if seq matching with current seed :
				  // 	clustering with seed if seq doesn't belong to any cluster yet
				  //                         OR in exact mode and the score is better with this seed
					if (! lcsmode && normalize)
						score = 1.0 - score;
					putSeqInClusterMT((void*)control, seed_number, i, score);		// saves new seed for this seq
				}
    		} // if ((! fast) || on_current_seed)

		} // for (i=seed_number+1;...

		#ifdef OMP_SUPPORT
		#pragma omp flush(nextone)
		#endif
    	if (priv_nextone < nextone)
			#ifdef OMP_SUPPORT
			#pragma omp critical
			#endif
    		if (priv_nextone < nextone)
    			nextone=priv_nextone;

	}

	gettimeofday(&current,NULL);
	elapsedtime = timevaldiff(&start,&current);
	control->elapsedtime+=elapsedtime;

	control->next=nextone;

    if (control->next < (control->n)-1)
    	(control->seeds_counter)++;
    else if (control->next == (control->n)-1)
	{
		control->stop = TRUE;
		(control->seeds_counter)++;
	}
	else if (control->next == control->n)
		control->stop = TRUE;
}


void initializeCentersAndScores(void *c)
{
	// Initializing the scores table for each seq :

    thread_control_t *control = (thread_control_t*) c;
    int32_t  i;
    fastaSeqPtr *db_i;
    int scoremax;

	if (control->normalize && control->lcsmode)
		scoremax = 1.0;
	else if (!control->lcsmode)
		scoremax = 0.0;
	else
		scoremax = (*(control->db))->length;

	for (i=0, db_i = control->db;
			i <= control->n-1;
			i++,db_i++)
	{
		(*db_i)->center         = (control->db)+i;
		(*db_i)->center_index   = i;
		(*db_i)->score          = scoremax;
		(*db_i)->cluster_center = TRUE;
	}
}


void freeEverything(void *c)
{
    thread_control_t *control=(thread_control_t*)c;
    int i;

	// free(control->potential_nexts_list);
	if ((control->reference == ALILEN) && (control->normalize || !control->lcsmode))
	{
		for (i=0; i < control->threads_number; i++)
			free(control->addresses[i]);
		free(control->addresses);
	}
	free(control->iseqs1);
	free(control->iseqs2);
}


int mt_compare_sumaclust(fastaSeqPtr* db, int n, BOOL fast, double threshold, BOOL normalize,
		int reference, BOOL lcsmode, int threads_number, double max_ratio)
{
    thread_control_t control;
	int32_t  i;
	int lmax, lmin;

	if (lcsmode || normalize)
		fprintf(stderr,"Clustering sequences when similarity >= %lf\n", threshold);
	else
		fprintf(stderr,"Clustering sequences when distance <= %lf\n", threshold);

	fprintf(stderr,"Aligning and clustering... \n");

	#ifdef OMP_SUPPORT
    control.threads_number = omp_get_max_threads();
	#else
    control.threads_number = 1;
    #endif
    if (threads_number < control.threads_number)
    	control.threads_number = threads_number;

	calculateMaxAndMinLen(db, n, &lmax, &lmin);

    control.addresses  = (int16_t**) malloc(control.threads_number*sizeof(int16_t*));
    control.iseqs1     = (int16_t**) malloc(control.threads_number*sizeof(int16_t*));
    control.iseqs2     = (int16_t**) malloc(control.threads_number*sizeof(int16_t*));

    for (i=0; i < control.threads_number; i++)
    	control.sizeForSeqs = prepareTablesForSumathings(lmax, lmin, threshold, normalize, reference, lcsmode, (control.addresses)+i, (control.iseqs1)+i, (control.iseqs2)+i);

	control.db = db;
	control.next = 0;
    control.normalize = normalize;
    control.reference = reference;
    control.threshold = threshold;
    control.max_ratio = max_ratio;
    control.lcsmode = lcsmode;
    control.stop = FALSE;
    control.fast = fast;
    control.seeds_counter = 1;
//    control.potential_nexts_list = (int*) calloc(control.threads_number, sizeof(int));
    control.n = n;

    if (lcsmode || normalize)
		control.worstscore = 0.0;
	else
		control.worstscore = lmax;

    control.elapsedtime=0;

    fprintf(stderr, "%d threads running\n", control.threads_number);

    // initialize scores table :
    initializeCentersAndScores(&control);

	while (control.stop == FALSE)
	{
		if ((control.next)%100 == 0)
		{
			float p = ((float)(control.next)/(float)n)*100;
			fprintf(stderr,"\rDone : %f %%      ",p);
		}
		computeOneSeed(&control);

	}

	for (i=0; i < control.threads_number; i++)
	{
		free((*((control.iseqs1)+i))-(control.sizeForSeqs)+lmax);
		free((*((control.iseqs2)+i))-(control.sizeForSeqs)+lmax);
	}

	freeEverything(&control);
    fprintf(stderr,"\rDone : 100 %%       %d clusters created.                        \n\n", control.seeds_counter);
    fprintf(stderr,"Pure computation time %f                        \n\n", (double)control.elapsedtime/1000000.);


   return(control.seeds_counter);
}
