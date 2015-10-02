/* ---------------------------------------------------------------- */
/* Copyright (c) Atelier de BioInformatique                         */
/* @file: libkov_proba.c                                            */
/* @desc: prokov probabilities computation                          */
/*                                                                  */
/* @history:                                                        */
/* @+   <Fred Nikitin + Marc Heuveline> : Aug 99 : first version    */
/* @+   <Gloup> : Oct 99 : last revised version                     */
/* ---------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "libkov.h"


#ifdef SGI
#define FLOG10(x) flog10(x)     /* in math library on sgi */
#define FPOW(x,y) powf(x,y)     /*                        */
#else
#define FLOG10(x) (float)log10((double)(x))
#define FPOW(x,y) (float)pow((double)(x), (double)(y)) 
#endif


/* ---------------------------------------------------- */
/* calcul des probas conditionnelles et des probas      */
/* initiales a partir du tableau des occurences         */
/*                                                      */
/* note: on stocke en fait le log10 des probas          */
/* ---------------------------------------------------- */
static void sProbaCond (int kuple, int occ[POW_KMAX],
                float proba[POW_KMAX], float pinit[POW_KMIN])
{
    int i, j, powk, somme, count;

    powk = Pow4(kuple);
    
    /* ------------------------------------------------- */     
    /* calcul somme des occurences                       */
    /* ------------------------------------------------- */     

    somme = 0;

    for (j = 0 ; j < powk ; j++)
        somme += occ[j];

    /* ------------------------------------------------- */     
    /* calcul des probas initiales et des probas         */
    /* conditionnelles                                   */
    /* ------------------------------------------------- */     
    /* les occurence sont rangees dans l'ordre           */
    /* x1...xk-1A, x1...xk-1C, x1...xk-1G, x1...xk-1T    */
    /* donc les frequences des x1...xk-1 correspondent   */
    /* aux sommes 4 par 4 des elements de occ[]          */
    /* ------------------------------------------------- */     

    for (j = 0 ; j < powk ; j += 4) {
                                     /* -1- probas initiales      */
        count = 0;                   /* frequence de x1...xk-1    */
        for (i = 0 ; i < 4 ; i++)
            count += occ[j+i];

        pinit[j/4] = FLOG10((float) count / (float) somme);

                                     /* -2- probas conditionelles */
        for (i = 0 ; i < 4 ; i++)
            proba[j+i] = FLOG10((float) occ[j+i] / (float) count);
    
    }
}

/* ---------------------------------------------------- */
/* calcul des probas conditionnelles et des probas      */
/* initiales a partir du tableau des occurences         */
/*                                                      */
/* API interface call                                   */
/* ---------------------------------------------------- */

void ProbaCond(MarkovMatrix *mat)
{
    int i;
    
    for (i = 0 ; i < 3 ; i++) {
        sProbaCond(mat->kupleC, mat->occpos[i], mat->probpos[i],
                   mat->pinipos[i]);
        sProbaCond(mat->kupleC, mat->occneg[i], mat->probneg[i],
                   mat->pinineg[i]);
    }
    
    sProbaCond(mat->kupleN, mat->occnc, mat->probnc, mat->pininc);
}

/* ---------------------------------------------------- */
/* calcule la probabilite de la chaine seq[0:winmax-1]  */
/* P(window/Codi)                                       */
/*                                                      */
/* ---------------------------------------------------- */
void ProbaMarkov (char *seq, int seqlen, MarkovMatrix *mat,
                  ProbaArray *win)
{
   int i, phase, phi, val; 

   for (phase = 0 ; phase < 3 ; phase++ ) {
                                                /* -------------------- */
       val = Codage(seq, mat->kupleC-1);        /* k-1 uple initial     */
       win->probpos[phase] = mat->pinipos[phase][val];
       win->probneg[phase] = mat->pinineg[phase][val];

       for (i = 0 ; i < seqlen ; i++ ) {        /* autres kuples        */

           val = (  (i == 0) 
                  ? Codage(seq, mat->kupleC)
                  : SuiteCodage(seq + i, val, mat->kupleC));
                    
           phi = (phase+i)%3;
           win->probpos[phase] += mat->probpos[phi][val]; /* probas en log */
           win->probneg[phase] += mat->probneg[phi][val];
       } 
   }

   val = Codage(seq, mat->kupleN-1);
   win->probnc = mat->pininc[val];
   
   for (i = 0 ; i < seqlen ; i++) {          
       val = (  (i == 0)
              ? Codage(seq, mat->kupleN)
              : SuiteCodage(seq + i, val, mat->kupleN));

       win->probnc += mat->probnc[val];
   }
}

/* ---------------------------------------------------- */
/* applique la formule d'inversion de Bayes pour        */
/* calculer P(Codi/wind)                                */
/* ---------------------------------------------------- */
void ProbaBayes(ProbaArray *win, float prior,
                ProbaArray *bay)
{
   int   i, j;
   float somme, pcod, pnc;
   
   pcod = FLOG10(prior/6.);             /* proba codant APRIORI     */
   pnc  = FLOG10(1.0 - prior);          /* proba non-codant APRIORI */

   for (i = 0 ; i < 3 ; i++) {
       somme = 0.;
       for (j = 0 ; j < 3 ; j++) {
           somme += FPOW(10., win->probpos[j] - win->probpos[i]);
           somme += FPOW(10., win->probneg[j] - win->probpos[i]);
       }
       somme += FPOW(10., win->probnc + pnc - win->probpos[i] - pcod);
       bay->probpos[i] = 1. / somme;
   }

   for (i = 0 ; i < 3 ; i++) {
       somme = 0.;
       for (j = 0 ; j < 3 ; j++) {
           somme += FPOW(10., win->probpos[j] - win->probneg[i]);
           somme += FPOW(10., win->probneg[j] - win->probneg[i]);
       }
       somme += FPOW(10., win->probnc + pnc - win->probneg[i] - pcod);
       bay->probneg[i] = 1. / somme;
   }

   somme = 0.;
   for (j = 0 ; j < 3 ; j++) {
           somme += FPOW(10., win->probpos[j] + pcod - win->probnc - pnc);
           somme += FPOW(10., win->probneg[j] + pcod - win->probnc - pnc);
   }
   somme += 1.;
 
   bay->probnc = 1. / somme;
}

