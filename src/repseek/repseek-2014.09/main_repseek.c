/*
	Copyright (C) 2006 G achaz, F boyer, E rocha, A viari and E coissac.

	This program is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public License
	as published by the Free Software Foundation; either version 2.1
	of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 
 	for more information, please contact guillaume achaz <achaz@abi.snv.jussieu.fr>
*/
/**
 * @file   main_repseek.c
 * @author Guillaume Achaz <achaz@abi.snv.jussieu.fr>, Eric Coissac <coissac@inrialpes.fr>
 * @date   Feb 26 2004
 * @modif  April 2004
 * @modif  Dec 2004 Adds for option -D : managing of mask for seed detection <EC>
 * @modif  sept 2005 re-write all align_blast2like. Correct Memory bugs <GA>
 * @modif  17 oct 2005 - mode '4' <GA>
 * @modif  dec 2005 - remove mode choice <GA>
 * @modif  jan-fev-march 2006 - stable version (-l, -p or -P options) <GA>
 * @modif  may 2008 - count Xs in the sequence and remove them to compute smin <GA>
 * 
 * @brief  main function of repseek program
 * 
 * 
 */


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if defined(MACOSX) || defined(LINUX) || defined(OSF)  /* for options */
#include <unistd.h>
#endif

#include "repseek_types.h"

#include "sequence.h"
#include "readfst.h"
#include "output.h"
#include "help.h"
#include "families.h"

#include "KMRK.h"
#include "KMRK_Seeds.h"
#include "KMRK_merge_seeds.h"
#include "KMRK_mask.h"

#include "read_seeds.h"

#include "filter.h"
#include "align.h"

#include "memory.h"

#include "smin.h"
#include "lmin.h"

int main(int argc, char **argv){

	
	/*
		Repseek Core Values
	*/ 

	int32_t selected_l=0;             /* the lmin chosen by the user --if any-- */
	
	int32_t lmin=0;                   /* the statistical lmin for seeds detection */

	double selected_s=0;                    /* estimated score minimum if repseek is running in mode '3' */
	double smin=0;                    /* estimated score minimum if repseek is running in mode '3' */

	float pval_lmin=0.0;              /* a pval that gives an Lmin (Karlin & Ost equation) */
	float pval_smin=0.0;              /* a pval that gives an Smin (Watermann & Vingron regression) */


	/*
		Sequence(s)
	*/

	int8_t repseek_nseq;               /* either 1 or 2, 1 being for 1 sequence and 2 for 2 sequences */

	char *sequence;                   /* the first sequence */
	char *sequence2=NULL;             /* the second sequence --only use when repseek_nseq is 2 */
	int32_t size;                     /* the size of the sequence */
	int32_t size2=0;                  /* the size of the second sequence --only use when repseek_nseq is 2*/
	int32_t totalSize;                /* size of (both) sequence(s) */

	int32_t nX=0,                     /* number of X in the first and the second sequence */
	        nX2=0;


	masked_area_table_t *mask=NULL;   /* masking table used during seed detection */
	char *mask_file=NULL;
    
	char opt_shape = 'l';             /* chromosome has to be set at 'l'inear or 'c'ircular */

	

	/*
		Output/Input File(s)
	*/
	char *output_file = NULL;         /* the output file name base */
	char *seed_file  = NULL;          /* the seed file that can be used instead of the seed detection */

	FILE *fout=stdout;                /* the output flux */
	

	/*
		Repseek Options
	*/
		
	int optc;                        /* what U need for options parsing */
	extern int optind;
	extern char *optarg;

	int8_t opt_dir = 1;              /* if 0 user do not want direct repeats */
	int8_t opt_inv = 1;              /* if 0 user do not want inverted repeats */  

	int8_t opt_MergeSeeds=0;         /* if 1, thus user ask for seed merging - remove overlapping seeds */
	int8_t opt_filter=0;             /* if 1, thus user ask for filtering */
	float min=0.0,                   /* if they are non-0, they are used to filter the seeds */
	      Max=0.0;                   /* min is the minimum of n-plication, Max is the maximum */

	float merge_repeats=0.90;        /* merge repeats when x% of both copies is located at the same spot */

	int8_t opt_TableR=0;             /* do we output TableR - 0 = no, default */
	int8_t opt_PrintSeeds=0;         /* do we output Seeds - 0 = no, default */
			               
	int8_t opt_overlap=1;            /* if 0, repeats cannot have their two copies overlapping */
			               


	/*
		repseek storage
	*/
	AllSeeds_type* AllSeeds;         /* A structure that contains all seeds. Defined in KMRK_Seeds.h */
			               
	Repeats AllRepeats;              /* A structure that contains all repeats. */
			               
	int32_t *TableR=NULL;            /* A pointer to the tableR */
	int32_t *TableR2=NULL;           /* A pointer to a second tableR -- If there are 2 sequences */

	
	/*
		Repseek alignment variables
	*/
	float Xg=30.0;                   /* The exploratory value. for BLAST2. same notation than in BLAST2 paper */
	float gap_open=4.0,              /* gap open and gap_ext penality - expressed as a coefficient of the mean score. */
	      gap_ext=1.0;
	char matrix_type='l';            /* matrix is by default 'l'og, can also be set to 'i'dentity */
			             
	

	/*
		Parse options
	*/
	
	while( (optc=getopt(argc, argv, "Hhvl:p:L:P:TSdis:m:M:co:e:X:O:ID:r:BR:")) != -1)

		switch(optc){
		
			/*
				Information
			*/
			case 'h':
				printhelp_repseek(argv);
				exit(1);

			case 'H':
				printmorehelp_repseek(argv);
				exit(1);

			case 'v':
			   	printversion_repseek();
				exit(1);
		
			/*
				Main Inputs
			*/
			case 'l':
				selected_l = (int16_t)atoi(optarg);
				if(selected_l < 1)
					fprintf(stderr,"selected lmin has to be >= 1, bye\n"),exit(1);
				break;
			
			case 'p':
				pval_lmin = atof(optarg);
				if(pval_lmin > 1.0 || pval_lmin < 0.0)fprintf(stderr,"set 0.0>pval_lmin>=1.0, bye\n"),exit(1);
				break;

			case 'L':
				selected_s = atof(optarg);
				if(selected_s < 0.0 )fprintf(stderr,"selected smin has to be >=0.0, bye\n"),exit(1);
				break;
				
			case 'P':
				pval_smin = atof(optarg);
				if(pval_smin > 1.0 || pval_smin < 0.0)fprintf(stderr,"set 0.0>pval_smin>=1.0, bye\n"),exit(1);
				break;
				
			/*
				Output/Input
			*/
			case 'r':
				output_file = optarg;
				break;
				
			case 'T':
				opt_TableR = 1;
				break;

			case 'S':
				opt_PrintSeeds = 1;
				break;

			case 's':
				seed_file = optarg;
				break;

			/*
				Options for detection
			*/
			case 'd':
				if(!opt_dir)
					fprintf(stderr,"options '-d' and '-i' are mutually exclusive, bye\n"),exit(1);
				opt_inv=0;
				break;

			case 'i':
				if(!opt_inv)
					fprintf(stderr,"options '-d' and '-i' are mutually exclusive, bye\n"),exit(1);
				opt_dir=0;
				break;
		
			/*
				Filtering
			*/
			case 'm':
				min = atof(optarg);
				opt_filter = 1;
				break;

			case 'M':
				Max = atof(optarg);
				opt_filter = 1;
				break;

			case 'B':
				opt_MergeSeeds=1;
				break;

			case 'R':
				merge_repeats= atof(optarg);
				if(merge_repeats<0.0 || merge_repeats>1.0)
					fprintf(stderr,"for merging repeats: 0.0 <= R <= 1.0 \n"),exit(1);
				break;


			
			/*
				Alignment options
			*/
			case 'o':
				gap_open = atof(optarg);
				break;

			case 'e':
				gap_ext = atof(optarg);
				break;

			case 'X':
				Xg = ABS(atof(optarg));
				break;

			case 'O':
				opt_overlap = (int8_t )atoi(optarg);
				if(opt_overlap != 0 && opt_overlap != 1)
					fprintf(stderr,"opt_overlap can only be 1-can overlapp or 0- no overlap\n, bye"),exit(1);
				break;

			/*
				Misc options
			*/
			case 'c':
				opt_shape = 'c';
				break;

			case 'I':
				matrix_type='i';
				break;

			case 'D':
				mask_file = optarg;
				break;
            
		}

	/*
		Do we have the minimum to run ? (mode minimum input and an output)
	*/
	repseek_nseq = argc-optind;
		
	if((repseek_nseq < 1) || (repseek_nseq > 2)  )
		printusage_repseek(argv),exit(1);

	if( !selected_l && !pval_lmin && !pval_smin && !selected_s)
		fprintf(stderr,"repseek minium core value is either -l/-p (i.e. set seed size lmin) or -L/-P (i.e. set repeat score smin)\n"),exit(1);
	
	if( selected_l && pval_lmin )
		fprintf(stderr,"Only one between -l, -p should be used (i.e. set seed lmin)\n"),exit(1);
	
	if( selected_s && pval_smin )
		fprintf(stderr,"Only one between -L, -P should be used (i.e. set repeat smin)\n"),exit(1);
	
	
	if( pval_smin && (gap_open != 4.0 || gap_ext != 1.0   ) )
		fprintf(stderr,"Using the Pval on repeats, gaps penalties should be default ones\n"),exit(1);


	if(output_file){
		fout = fopen(output_file,"w");
		if(!fout) fprintf(stderr, "Cannot write into %s file, bye\n", output_file),exit(1);
	}

	fprintf(stderr,"Repseek mode...    '%s' + '%s'\n",
	     pval_lmin?"seed stat":"seed no stat",
	     pval_smin?"repeat stat":"repeat no stat"
	     );


	/*
		Get sequence(s) in memory
	*/

	fprintf(stderr,"get sequence(s)... ");

	sequence = readFastaSeq(argv[optind], &size, &nX); 
	totalSize = size;

	if(repseek_nseq==2){
		sequence2 = readFastaSeq(argv[optind+1], &size2, &nX2);
		totalSize+=size2;
	}
	

	/*
		Output Sequence Length and Xs
	*/
	if(totalSize < 1e3)
	  fprintf(stderr,"total size %d b ",    totalSize  );
	else if(totalSize < 1e6)
	  fprintf(stderr,"total size %.2f kb ", ((float)totalSize)/1e3 );
	else if(totalSize < 1e9)
	  fprintf(stderr,"total size %.2f Mb ", ((float)totalSize)/1e6 );
	else 
	  fprintf(stderr,"total size %.2f Gb ", ((float)totalSize)/1e9 );
	 
	switch( (char)floor( log(nX+nX2)/log(10) ) ){
		case 0: case 1: case 2:
			 fprintf(stderr,"(%d Xs)\n",   (nX+nX2) );
			 break;
		case 3:	case 4: case 5:
			fprintf(stderr,"(%.2f K Xs)\n",   (float)(nX+nX2)/1e3 );
			 break;
		case 6:	case 7: case 8:
			fprintf(stderr,"(%.2f M Xs)\n",   (float)(nX+nX2)/1e6 );
			 break;
		default:
			fprintf(stderr,"(%.2f G Xs)\n",   (float)(nX+nX2)/1e9 );
			break;
	
	}
 


	/*
		handle the seed minimum size
	*/
	fprintf(stderr,"set lmin...        ");
	
	if(selected_l)
		lmin=selected_l; 
	else{
		
		if(pval_lmin){
			if (repseek_nseq==2)
				lmin=set_lmin_2seqs(pval_lmin, sequence,sequence2,size,size2);     /* Lmin is computed ignoring X and N */
			else
				lmin=set_lmin(pval_lmin, sequence);
		}
		else{
			if (repseek_nseq==2)
				lmin=(int32_t)set_lmin_2seqs(0.001, sequence,sequence2,size,size2)*2.0/3.0;
			else
				lmin=(int32_t)set_lmin(0.001, sequence)*2.0/3.0;
		}
		
	}
	fprintf(stderr, "seed >= %d bp ( L5%%: %d ; L1%%: %d ; L0.1%%: %d )\n",
		lmin,
		(repseek_nseq==2)?set_lmin_2seqs(0.05, sequence,sequence2,size,size2):set_lmin(0.05, sequence),
		(repseek_nseq==2)?set_lmin_2seqs(0.01, sequence,sequence2,size,size2):set_lmin(0.01, sequence),
		(repseek_nseq==2)?set_lmin_2seqs(0.001, sequence,sequence2,size,size2):set_lmin(0.001, sequence)
	);


	/*
		Select the repeat minimum score
	*/
	fprintf(stderr,"set smin...        ");
	if(selected_s)
		smin=selected_s;        /* set the repeat minimum score as given in the input */
	else{
		
		if(pval_smin){
			if (repseek_nseq==2)
				smin = compute_smin_2seq(sequence, sequence2, size,  size2, pval_smin, nX, nX2);
			else
				smin = compute_smin(sequence, pval_smin, size, nX);
		}
		else
			smin = 0;
		
	}
		
		
	fprintf(stderr, "repeat >= %.2f score ( S5%%: %.2f ; S1%%: %.2f ; S0.1%%: %.2f )\n",
		smin,
		(repseek_nseq==2)?compute_smin_2seq(sequence, sequence2, size,  size2, 0.05, nX, nX2):compute_smin(sequence, 0.05,  size, nX),
		(repseek_nseq==2)?compute_smin_2seq(sequence, sequence2, size,  size2, 0.01, nX, nX2):compute_smin(sequence, 0.01, size, nX),
		(repseek_nseq==2)?compute_smin_2seq(sequence, sequence2, size,  size2, 0.001, nX, nX2):compute_smin(sequence, 0.001, size, nX)
	);
		
	

	/*
		Begin the seed detection by KMR-K or read file
	*/

	if(seed_file)  	fprintf(stderr,"read seeds...      ");
	else            fprintf(stderr,"detect seeds...    ");
    
	if (repseek_nseq==1)
	{
		if (mask_file)
		{
			if (opt_inv)
				mask = KMRK_ReadMaskArea(mask_file,1,size);  /* 2nd arg: #seq and the 3rd arg: complement */
			else
				mask = KMRK_ReadMaskArea(mask_file,1,0);
		}

		if(seed_file)
			AllSeeds =  read_seeds(seed_file, lmin, opt_dir, opt_inv, size, 0);
		else
			AllSeeds =  KMRK_get_seeds(&sequence, size, lmin, opt_dir, opt_inv, 0, mask);
	
	}else{
	
		if (mask_file)
		{
			if (opt_inv)
				mask = KMRK_ReadMaskArea(mask_file,2,size);
			else
				mask = KMRK_ReadMaskArea(mask_file,2,0);
		}
		
		if(seed_file)
			AllSeeds =  read_seeds(seed_file, lmin, opt_dir, opt_inv, size, size2);
		else
			AllSeeds =  KMRK_get_seeds_2seqs(&sequence, &sequence2, size, size2,  lmin, opt_dir, opt_inv, 0, mask);
	
	}

	
	fprintf(stderr, "%d direct & %d inverted seeds\n", AllSeeds->nDirSeeds, AllSeeds->nInvSeeds);
	
	
	/*
		Optionnal filters
	*/
	if (opt_filter)	
	{
		fprintf(stderr,"filter, keep...    ");	
		if (repseek_nseq==1)
			KMRK_FamilySeeds(AllSeeds, min, Max, opt_dir, opt_inv, size);
		else
			KMRK_FamilySeeds2seq(AllSeeds, min, Max, opt_dir, opt_inv, size, size2);

		fprintf(stderr, "%d direct and %d inverted seeds\n", AllSeeds->nDirSeeds, AllSeeds->nInvSeeds);
	}
	
	if(opt_MergeSeeds){

		fprintf(stderr,"merge, keep...    ");	
		KMRK_MergeSeeds(AllSeeds, opt_dir, opt_inv);
		
		fprintf(stderr,"%d direct and %d inverted seeds\n",  AllSeeds->nDirSeeds, AllSeeds->nInvSeeds);
	}
	
	
	/*
		Final table(s) R of seeds
	*/
	/* filter seeds if they does not fit between min and Max */
		
	if (repseek_nseq==1){
		TableR = KMRK_SeedTableR(AllSeeds, opt_dir, opt_inv, size);
		BuiltFilterFamily_seeds(AllSeeds, TableR, 0, 0, opt_dir, opt_inv);
	}
	else{
		KMRK_SeedTableR2seq(AllSeeds, opt_dir, opt_inv, size,size2, &TableR,&TableR2);
		BuiltFilterFamily_seeds2seq(AllSeeds, TableR, TableR2, 0, 0, opt_dir, opt_inv);      
	}


	/*
		reorder AllSeeds by position
	*/
	if (opt_dir) KMRK_sortSeedsByPos(AllSeeds->dirSeeds,AllSeeds->nDirSeeds);
	if (opt_inv) KMRK_sortSeedsByPos(AllSeeds->invSeeds,AllSeeds->nInvSeeds);


 
 
 	if( opt_PrintSeeds ){

		fprintf(stderr,"ouput seeds...    ");

		if(repseek_nseq == 1)
			WriteSeeds(fout,  AllSeeds, opt_dir, opt_inv);
		else
			WriteSeeds(fout,  AllSeeds, opt_dir, opt_inv);

		if(opt_TableR){
		
			fprintf(fout,"# Seed TableR-1\n");
			WriteTableR(fout, TableR, size);
			if (repseek_nseq == 2){
				fprintf(fout,"# Seed TableR-2\n");
				WriteTableR(fout, TableR2, size2);
			}
		}
		
		fprintf(stderr," done\n");
		
		/*
			Free stuff
		*/
		
		/*  seeds  */
		KMRK_freeSeeds(AllSeeds);

		/* sequence */
		MyFree(sequence, (size+1)*sizeof(char));
		if(repseek_nseq == 2)
			MyFree(sequence2, (size2+1)*sizeof(char));
		
		/* TableR */
		MyFree(TableR, size*sizeof(int32_t));
		if(repseek_nseq==2)MyFree(TableR2, size2*sizeof(int32_t));


		PrintMaxMem();
		exit(0);
	}

	/*
		Free TableR of seeds
	*/
	MyFree(TableR, size*sizeof(int32_t));
	if(repseek_nseq==2)MyFree(TableR2, size2*sizeof(int32_t));

	  
	/*
		Extension : BLAST2-like Alignments
	*/
	if(repseek_nseq == 1)
		AllRepeats = align_seeds(AllSeeds, sequence, Xg, gap_open, gap_ext, matrix_type,
		                         lmin, opt_dir, opt_inv, opt_overlap, merge_repeats);
	else
		AllRepeats = align_seeds_2seq(AllSeeds, sequence, sequence2,  Xg, gap_open, gap_ext, matrix_type,
                                              lmin, opt_dir, opt_inv, merge_repeats);
				
	fprintf(stderr,"align and merge... %d direct & %d inverted repeats (merge with R=%.2f)\n",
	AllRepeats.nDirRep-AllRepeats.nDirBadRep, AllRepeats.nInvRep-AllRepeats.nInvBadRep, merge_repeats);


	/*
		If mode is stat on repeats, remove repeats with a score smaller than smin
	*/
	if( smin ){
		UpdateRepeats(&AllRepeats, smin);
		fprintf(stderr,"After score_min... %d direct & %d inverted repeats\n",
		AllRepeats.nDirRep-AllRepeats.nDirBadRep, AllRepeats.nInvRep-AllRepeats.nInvBadRep);
	}

	/*
		Built the repeats table(s) R.
	*/
	if (repseek_nseq==1)
		TableR = SetRepFamily(AllRepeats, opt_dir, opt_inv, size);
	else
		SetNumberRep_2seqs(AllRepeats, opt_dir, opt_inv, size, size2, &TableR, &TableR2);

	
	/*
		Misc Output
	*/
	
	fprintf(stderr,"output repeats...  ");


       if(repseek_nseq == 1)
               WriteRep(fout,  AllRepeats, opt_dir, opt_inv, opt_shape, size);
       else
               WriteRep_2seqs(fout,  AllRepeats, opt_dir, opt_inv, size, size2);

		
	if(opt_TableR){
		fprintf(fout,"# TableR-1\n");
		WriteTableR(fout, TableR, size);
		if (repseek_nseq == 2){
			fprintf(fout,"# TableR-2\n");
			WriteTableR(fout, TableR2, size2);
		}
	  }


	fprintf(stderr,"done\n");	
	
	fclose(fout);
	
	
	/*
		All free functions
	*/

	/* table(s) R */
	MyFree(TableR, size*sizeof(int32_t));
	if(repseek_nseq == 2)
		MyFree(TableR2, size2*sizeof(int32_t));

	/* seeds and repeats */
	KMRK_freeSeeds(AllSeeds);
	free_Repeats(AllRepeats);

	/* sequence */
	MyFree(sequence, (size+1)*sizeof(char));
	if(repseek_nseq == 2)
		MyFree(sequence2, (size2+1)*sizeof(char));

	/* PrintMem("At the end"); */
	
	PrintMaxMem();
	
	return 0;
}

