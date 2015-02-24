//
//  wilcox.cpp
//  Mothur
//
//  Created by SarahsWork on 8/6/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//
//  Modified to avoid using Rmath.h 

#include "mothur.h"

/*
 Mathlib : A C Library of Special Functions
 Copyright (C) 1999-2012  The R Core Team
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or (at
 your option) any later version.
 
 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the GNU
 General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, a copy is available at
 http://www.r-project.org/Licenses/
 
 SYNOPSIS
 
 #include <Rmath.h>
 double dwilcox(double x, double m, double n, int give_log)
 double pwilcox(double x, double m, double n, int lower_tail, int log_p)
 double qwilcox(double x, double m, double n, int lower_tail, int log_p);
 double rwilcox(double m, double n)
 
 DESCRIPTION
 
 dwilcox	The density of the Wilcoxon distribution.
 pwilcox	The distribution function of the Wilcoxon distribution.
 qwilcox	The quantile function of the Wilcoxon distribution.
 rwilcox	Random variates from the Wilcoxon distribution.
 
 */

/*
 Note: the checks here for R_CheckInterrupt also do stack checking.
 
 calloc/free are remapped for use in R, so allocation checks are done there.
 freeing is completed by an on.exit action in the R wrappers.
 */
#include "wilcox.h"

#define WILCOX_MAX 50
#define give_log log_p
#define R_D__0	(log_p ? 1 : 0.)		/* 0 */
#define R_D__1	(log_p ? 0. : 1.)			/* 1 */
#define R_DT_0	(lower_tail ? R_D__0 : R_D__1)		/* 0 */
#define R_DT_1	(lower_tail ? R_D__1 : R_D__0)		/* 1 */
#define R_D_val(x)	(log_p	? log(x) : (x))		/*  x  in pF(x,..) */
#define R_D_Clog(p)	(log_p	? log1p(-(p)) : (0.5 - (p) + 0.5)) /* [log](1-p) */
#define R_DT_val(x)	(lower_tail ? R_D_val(x)  : R_D_Clog(x))

/*********************************************************************************************************************************/
//Numerical Recipes pg. 219
double PWilcox::gammln(const double xx) {
    try {
        int j;
        double x,y,tmp,ser;
        static const double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,
            -1.231739572450155,0.120858003e-2,-0.536382e-5};
        
        y=x=xx;
        tmp=x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.0;
        for (j=0;j<6;j++) {
            ser += cof[j]/++y;
        }
        return -tmp+log(2.5066282746310005*ser/x);
    }
    catch(exception& e) {
        mout->errorOut(e, "PWilcox", "gammln");
        exit(1);
    }
}
/*********************************************************************************************************************************/
double PWilcox::choose(double n, double k){
    try {
        n = floor(n + 0.5);
        k = floor(k + 0.5);
        
        double lchoose = gammln(n + 1.0) - gammln(k + 1.0) - gammln(n - k + 1.0);
        
        return (floor(exp(lchoose) + 0.5));
    }
    catch(exception& e) {
        mout->errorOut(e, "PWilcox", "choose");
        exit(1);
    }
}
/*********************************************************************************************************************************/

/* This counts the number of choices with statistic = k */
double PWilcox::cwilcox(int k, int m, int n, double*** w) {
    try {
        int c, u, i, j, l;
        
        if (mout->control_pressed) { return 0; }
        
        u = m * n;
        if (k < 0 || k > u)
            return(0);
        c = (int)(u / 2);
        if (k > c)
            k = u - k; /* hence  k <= floor(u / 2) */
        if (m < n) {
            i = m; j = n;
        } else {
            i = n; j = m;
        } /* hence  i <= j */
        
        if (j == 0) /* and hence i == 0 */
            return (k == 0);
        
        
        /* We can simplify things if k is small.  Consider the Mann-Whitney
         definition, and sort y.  Then if the statistic is k, no more
         than k of the y's can be <= any x[i], and since they are sorted
         these can only be in the first k.  So the count is the same as
         if there were just k y's.
         */
        if (j > 0 && k < j) return cwilcox(k, i, k, w);
        
        if (w[i][j] == 0) {
            w[i][j] = (double *) calloc((size_t) c + 1, sizeof(double));

            for (l = 0; l <= c; l++)
                w[i][j][l] = -1;
        }
        if (w[i][j][k] < 0) {
            if (j == 0) /* and hence i == 0 */
                w[i][j][k] = (k == 0);
            else
                w[i][j][k] = cwilcox(k - j, i - 1, j, w) + cwilcox(k, i, j - 1, w);
            
        }
        return(w[i][j][k]);
    }
    catch(exception& e) {
        mout->errorOut(e, "PWilcox", "cwilcox");
        exit(1);
    }
}
/*********************************************************************************************************************************/

/* args have the same meaning as R function pwilcox */
double PWilcox::pwilcox(double q, double m, double n, bool lower_tail){
    try {
        int i;
        double c, p;
        bool log_p = false;
        double*** w;
        

        if (isnan(m) || isnan(n))
        {return 0;}
        m = floor(m + 0.5);
        n = floor(n + 0.5);
        if (m <= 0 || n <= 0) {return 0;}
        
        q = floor(q + 1e-7);
        
        if (q < 0.0)
            return(R_DT_0);
        if (q >= m * n)
            return(R_DT_1);
        
        int mm = (int) m, nn = (int) n;
        
        if (mout->control_pressed) { return 0; }
        //w_init_maybe(mm, nn);
        /********************************************/
        int thisi;
        if (mm > nn) { thisi = nn; nn = mm; mm = thisi; }
        
        mm = max(mm, 50);
        nn = max(nn, 50);
        w = (double ***) calloc((size_t) mm + 1, sizeof(double **));
            
        for (thisi = 0; thisi <= mm; thisi++) {
            w[thisi] = (double **) calloc((size_t) nn + 1, sizeof(double *));
        }
        allocated_m = m; allocated_n = n;
        /********************************************/    
        
        c = choose(m + n, n);
        p = 0;
        /* Use summation of probs over the shorter range */
        if (q <= (m * n / 2)) {
            for (i = 0; i <= q; i++)
                p += cwilcox(i, m, n, w) / c;
        }
        else {
            q = m * n - q;
            for (i = 0; i < q; i++) {
                p += cwilcox(i, m, n, w) / c;
            }
            lower_tail = !lower_tail; /* p = 1 - p; */
        }
        
        //free w
        /********************************************/
        for (int i = allocated_m; i >= 0; i--) {
            for (int j = allocated_n; j >= 0; j--) {
                if (w[i][j] != 0)
                    free((void *) w[i][j]);
            }
            free((void *) w[i]);
        }
        free((void *) w);
        w = 0; allocated_m = allocated_n = 0;
        /********************************************/
        
        return(R_DT_val(p));
    }
    catch(exception& e) {
        mout->errorOut(e, "PWilcox", "pwilcox");
        exit(1);
    }
} /* pwilcox */

