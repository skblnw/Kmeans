/*
 * Principal Component Analysis (PCA) of proteins.
 * By Zhiyong "John" Zhang, 12/13/2008.
 */
 
/*
 * PCA Anlysis and print covd.dat representing the geometrical covariance of atoms
 * Modified by Kevin Chun Chan, Apr 2017
 */
 
/*
 * IMPORTANT LOGS
 * 5/4/2017 output of covd added
 * 5/4/2017 output of norm of deviation vector commented out
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
    else if (strcmp(argv[i], "-cov") == 0)
    {
      fp[9] = fopen(argv[++i], "w");
    }
    else if (strcmp(argv[i], "-ocd") == 0)
    {
      fp[10] = fopen(argv[++i], "w");
    }
 }
 return(1);
}

int main(int argc, char *argv[])
{
 FILE *iref, *itrj, *imas, *ipar, *orms, *oave, *oval, *ocfl, *ovec, *fp[20];
// FILE *orco; // Added for covariance matrix of distance r=sqrt(x^2+y^2+z^2)
 FILE *ocov, *ocod;
 int i, j, k, l, m, n, d, nn, nfram, natm, nvec;
 long int ndim, il, iu, lm, lwork, liwork, info, *isuppz, *iwork;
 double rms, trace, vl, vu, abstol, cfl, *mass=NULL, *cov, *val, *vec, *work;
// double *rcov, *rave; // Added for covariance matrix of distance r=sqrt(x^2+y^2+z^2)
 double covdtemp, *covd;
 rvec *x, *xref, *xave;
 char str[STRLEN];

 if (argc == 1 || strcmp(argv[1], "-h") == 0)
 {
   fprintf(stderr, "usage: %s -irf REF -itj TRJ -imw mass -ipr control parameters -ors RMSD -oav AVE -ova eigenvalues -ocf cumulative fluctuations -ove eigenvectors -cov covariance_matrix -ocd cov_atom_distance\n", argv[0]);
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
 ocov = fp[9];
 ocod = fp[10];

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
// snew(rave, natm); // Added for covariance matrix of distance r=sqrt(x^2+y^2+z^2)
 snew(val, ndim);
 snew(cov, ndim*ndim);
// snew(rcov, natm*natm); // Added for covariance matrix of distance r=sqrt(x^2+y^2+z^2)
 snew(covd, natm*natm);
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
 printf("Reading coordinates\n");
 printf("Calculating RMSD and constructing the green covariance matrix\n");
 printf("Looping all structures in trajectory (nframe = %d)\n", nfram);
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
       /* add r to average structure, // Added for covariance matrix of distance r=sqrt(x^2+y^2+z^2) */
//       rave[i] += sqrt(x[i][0]*x[i][0] + x[i][1]*x[i][1] + x[i][2]*x[i][2]); 

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
       
       /* add components to covariance matrix of distance r, // Added for covariance matrix of distance r=sqrt(x^2+y^2+z^2) */
//        for (j=i; j<natm; j++)
//       {
//          rcov[natm*j+i] += sqrt( x[i][0]*x[i][0] + x[i][1]*x[i][1] + x[i][2]*x[i][2] ) * sqrt( x[j][0]*x[j][0] + x[j][1]*x[j][1] + x[j][2]*x[j][2] ); // Added for covariance matrix of distance r=sqrt(x^2+y^2+z^2)
//       }
       
       /* add componenets to covariance matrix of atom separation */
       for (j=i; j<natm; j++)
       {
           covdtemp = 0.0;
           for (l=0; l<DIM; l++)
           {
               covdtemp += x[i][l]*x[j][l];
           }
           covd[natm*i+j] += covdtemp / nfram;
       }
    }
 }

/* divide average by number of frames */
 for (i=0; i<natm; i++)
 {
//    rave[i] /= nfram; // Added for covariance matrix of distance r=sqrt(x^2+y^2+z^2)
    for (m=0; m<DIM; m++)
    {
       xave[i][m] /= nfram;
    }
    fprintf(oave, "%14.9lf %14.9lf %14.9lf\n", xave[i][0], xave[i][1], xave[i][2]);
 }

/* obtain the covariance matrix */
 printf("Finalizing the covariance matrix\n");
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

/* obtain the r covariance matrix, // Added for covariance matrix of distance r=sqrt(x^2+y^2+z^2) */
// printf("Finalizing the r covariance matrix\n");
// for (i=0; i<natm; i++)
// {
//    for (j=i; j<natm; j++)
//    {
//       rcov[natm*j+i] = rcov[natm*j+i]/nfram - rave[i]*rave[j];
//    }
// }

/* symmetrize the matrix */
 for (i=0; i<ndim; i++)
 {
    for (j=i; j<ndim; j++)
    {
       cov[ndim*i+j] = cov[ndim*j+i];
    }
 }

/* write the covariance matrix */
 for (i=0; i<ndim; i++)
 {
    for (j=0; j<ndim; j++)
    {
       fprintf(ocov, "%15.9lf ", cov[ndim*i+j]);
    }
    fprintf(ocov, "\n");
    fflush(ocov);
 }

/* write the r covariance matrix, // Added for covariance matrix of distance r=sqrt(x^2+y^2+z^2) */
// for (j=0; j<natm; j++)
// {
//    for (i=0; i<natm; i++)
//    {
//       fprintf(orco, "%15.9lf ", rcov[natm*j+i]);
//    }
//    fprintf(orco, "\n");
// }
 
/* write the atom separation matrix */
/* Size of covd.dat will be (natm+1)*natm/2 */
for(i=0; i<natm; i++)
{
    for(j=i; j<natm; j++)
    {
        fprintf(ocod, "%15.9lf\n", covd[natm*i+j]);
    }
    fflush(ocod);
}

/* calculate the trace of the matrix */
 trace = 0.0;
 for (j=0; j<ndim; j++)
 {
    trace += cov[ndim*j+j];
 }

/* call diagonalization routine */
 printf("Diagonalizing the covariance matrix. Be patient, almost done :)\n");
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
 fclose(ocov);
 fclose(ocod); // Added for covariance matrix of distance r=sqrt(x^2+y^2+z^2)

/* the end of main */
 return 0;
}
