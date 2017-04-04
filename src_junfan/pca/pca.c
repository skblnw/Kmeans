/*
 * Principal Component Analysis (PCA) of proteins.
 * By Zhiyong "John" Zhang, 12/13/2008.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "pca.h"

extern double cal_rms(int natm, double mass[], rvec x[], rvec xp[]);
extern void read_x(FILE *file, int natm, rvec *x);
extern void reset_x(int natm, rvec x[], double mass[]);
extern void do_fit(int natm, double *mass, rvec *xp, rvec *x);
extern int dsyevr_(char *, char *, char *, long int *, double *,
                   long int *, double *, double *, long int *,
                   long int *, double *, long int *, double *,
                   double *, long int *, long int *, double *,
                   long int *, long int *, long int *, long int *);

int rdpara(int argc, char *argv[], FILE **fp)
{
 int i;

 for (i=1; i<argc; i++)
 {
    if (strcmp(argv[i], "-irf") == 0)
    {
      fp[0] = fopen(argv[++i], "r");
    }
    else if (strcmp(argv[i], "-itj") == 0)
    {
      fp[1] = fopen(argv[++i], "r");
    }
    else if (strcmp(argv[i], "-imw") == 0)
    {
      fp[2] = fopen(argv[++i], "r");
    }
    else if (strcmp(argv[i], "-ipr") == 0)
    {
      fp[3] = fopen(argv[++i], "r");
    }
    else if (strcmp(argv[i], "-ors") == 0)
    {
      fp[4] = fopen(argv[++i], "w");
    }
    else if (strcmp(argv[i], "-oav") == 0)
    {
      fp[5] = fopen(argv[++i], "w");
    }
    else if (strcmp(argv[i], "-ova") == 0)
    {
      fp[6] = fopen(argv[++i], "w");
    }
    else if (strcmp(argv[i], "-ocf") == 0)
    {
      fp[7] = fopen(argv[++i], "w");
    }
    else if (strcmp(argv[i], "-ove") == 0)
    {
      fp[8] = fopen(argv[++i], "w");
    }
 }
 return(1);
}

int main(int argc, char *argv[])
{
 FILE *iref, *itrj, *imas, *ipar, *orms, *oave, *oval, *ocfl, *ovec, *fp[20];
 int i, j, k, l, m, n, d, nn, nfram, natm, nvec;
 long int ndim, il, iu, lm, lwork, liwork, info, *isuppz, *iwork;
 double rms, trace, vl, vu, abstol, cfl, *mass=NULL, *cov, *val, *vec, *work;
 rvec *x, *xref, *xave;
 char str[STRLEN];

 if (argc == 1 || strcmp(argv[1], "-h") == 0)
 {
   fprintf(stderr, "usage: %s -irf REF -itj TRJ -imw mass -ipr control parameters -ors RMSD -oav AVE -ova eigenvalues -ocf cumulative fluctuations -ove eigenvectors\n", argv[0]);
   exit(1);
 }

 if (rdpara(argc, argv, fp) == 0) exit(0);
 iref = fp[0];
 itrj = fp[1];
 imas = fp[2];
 ipar = fp[3];
 orms = fp[4];
 oave = fp[5];
 oval = fp[6];
 ocfl = fp[7];
 ovec = fp[8];

/* control parameters */
 fgets2(str, STRLEN, ipar); // one-line statement of the file

/* number of configurations */
 fgets2(str, STRLEN, ipar); // one-line comment
 fgets2(str, STRLEN, ipar);
 sscanf(str, "%d", &nfram);

/* number of atoms */
 fgets2(str, STRLEN, ipar); // one-line comment
 fgets2(str, STRLEN, ipar);
 sscanf(str, "%d", &natm);

/* number of PCA modes we want to get */
 fgets2(str, STRLEN, ipar); // one-line comment
 fgets2(str, STRLEN, ipar);
 sscanf(str, "%d", &nvec);

 ndim = natm*DIM;
 lwork = 30*ndim;
 liwork = 15*ndim;

/* allocate space */
 snew(xref, natm);
 snew(x, natm);
 snew(xave, natm);
 snew(val, ndim);
 snew(cov, ndim*ndim);
 snew(vec, ndim*nvec);
 snew(mass, natm);
 snew(isuppz, 2*nvec);
 snew(work, lwork);
 snew(iwork, liwork);

/* read the reference structure */
 read_x(iref, natm, xref);

/* translate the c.o.m of reference to origin */
 for (i=0; i<natm; i++)
 {
    fgets2(str, STRLEN, imas);
    sscanf(str, "%d%lf", &j, &mass[i]);
 }

 reset_x(natm, xref, mass);

/* start main loop of all the structures in trajectory */
 for (nn=0; nn<nfram; nn++)
 {

/* read frame. Note: you have to remove PBC by yourself before running this program */
    read_x(itrj, natm, x);

/* translate frame */
    reset_x(natm, x, mass);

/* perform least square fit */
    do_fit(natm, mass, xref, x);

/* calculate RMSD */
    rms = cal_rms(natm, mass, x, xref);
    fprintf(orms, "%8d% 10.4lf\n", nn+1, rms);
    fflush(orms);
    for (i=0; i<natm; i++)
    {

/* add positions to average structure */
       for (m=0; m<DIM; m++)
       {
          xave[i][m] += x[i][m];
       }

/* add componens to covariance matrix */
       for (l=0; l<DIM; l++)
       {
          k = DIM*i+l;
          for (j=i; j<natm; j++)
          {
             for (n=0; n<DIM; n++)
             {
                d = DIM*j+n;
                cov[ndim*d+k] += x[i][l]*x[j][n];
             }
          }
       }
    }
 }

/* divid average by number of frames */
 for (i=0; i<natm; i++)
 {
    for (m=0; m<DIM; m++)
    {
       xave[i][m] /= nfram;
    }
    fprintf(oave, "%14.9lf %14.9lf %14.9lf\n", xave[i][0], xave[i][1], xave[i][2]);
 }

/* obtain the covariance matrix */
 for (i=0; i<natm; i++)
 {
    for (l=0; l<DIM; l++)
    {
       k = DIM*i+l;
       for (j=i; j<natm; j++)
       {
          for (n=0; n<DIM; n++)
          {
             d = DIM*j+n;
             cov[ndim*d+k] = cov[ndim*d+k]/nfram - xave[i][l]*xave[j][n];
          }
       }
    }
 }

/* symmetrize the matrix */
 for (i=0; i<ndim; i++)
 {
    for (j=i; j<ndim; j++)
    {
       cov[ndim*i+j] = cov[ndim*j+i];
    }
 }

/* calculate the trace of the matrix */
 trace = 0.0;
 for (j=0; j<ndim; j++)
 {
    trace += cov[ndim*j+j];
 }

/* call diagonalization routine */
 vl = 0.0;
 vu = 1000.0;
 il = ndim-nvec+1;
 iu = ndim;
 abstol = 0.0;
 dsyevr_("V", "I", "U", &ndim, cov, &ndim, &vl, &vu, &il, &iu,
         &abstol, &lm, val, vec, &ndim, isuppz, work, &lwork, iwork,
         &liwork, &info);

/* calculate the cumulative fluctuations */
 cfl = 0.0;
 for (i=0; i<nvec; i++)
 {
    j = nvec-i-1;
    cfl += val[j];
    fprintf(ocfl, "%8d %10.4lf\n", i+1, cfl/trace);
 }

/* write the eigenvalues and eigenvectors */
 for (i=0; i<nvec; i++)
 {
    fprintf(oval, "%8d %12.4lf\n", i+1, val[nvec-1-i]);
 }

 for (j=1; j<=nvec; j++)
 {
    n = nvec-j;
    for (i=0; i<natm; i++)
    {
       for (m=0; m<DIM; m++)
       {
          x[i][m] = vec[ndim*n+i*DIM+m];
       }
       fprintf(ovec, "%14.9lf %14.9lf %14.9lf\n", x[i][0], x[i][1], x[i][2]);
    }
 }

/* close file */
 fclose(iref);
 fclose(itrj);
 fclose(imas);
 fclose(ipar);
 fclose(orms);
 fclose(oave);
 fclose(oval);
 fclose(ocfl);
 fclose(ovec);

/* the end of main */
 return 0;
}
