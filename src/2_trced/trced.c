/*
 * This program is used to calculate the convariance matrix in the
 * essential subspace. By Zhiyong Zhang, 04/01/2008.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "pca.h"

extern void read_x(FILE *file, int natm, rvec *x);

int rdpara(int argc, char *argv[], FILE **fp)
{
 int i;

 for (i=1; i<argc; i++)
 {
    if (strcmp(argv[i], "-iva") == 0)
    {
      fp[0] = fopen(argv[++i], "r");
    }
    else if (strcmp(argv[i], "-ivc") == 0)
    {
      fp[1] = fopen(argv[++i], "r");
    }
    else if (strcmp(argv[i], "-ipr") == 0)
    {
      fp[2] = fopen(argv[++i], "r");
    }
    else if (strcmp(argv[i], "-oce") == 0)
    {
      fp[3] = fopen(argv[++i], "w");
    }
 }
 return(1);
}

int main(int argc, char *argv[])
{
 FILE *ival, *ivec, *ipar, *oced, *fp[10];
 int i, j, m, n, idum, natm, ndim, ned;
 double cov, *val;
 rvec *x, **vec;
 char str[STRLEN];

 if (argc == 1 || strcmp(argv[1], "-h") == 0)
 {
   fprintf(stderr, "usage: %s -iva eigenvalues -ivc eigenvectors -ipr control parameters -oce cov-ed\n", argv[0]);
   exit(1);
 }

 if (rdpara(argc, argv, fp) == 0) exit(0);
 ival = fp[0];
 ivec = fp[1];
 ipar = fp[2];
 oced = fp[3];

/* control parameters */
 fgets2(str, STRLEN, ipar); // one-line statement of the file

/* number of atoms */
 fgets2(str, STRLEN, ipar); // one-line comment
 fgets2(str, STRLEN, ipar);
 sscanf(str, "%d", &natm);

/* number of essential PCA modes */
 fgets2(str, STRLEN, ipar); // one-line comment
 fgets2(str, STRLEN, ipar);
 sscanf(str, "%d", &ned);

 ndim = natm*DIM;

/* allocate space */
 snew(x, natm);
 snew(val, ndim);
 snew(vec, ned);
 for (n=0; n<ned; n++)
 {
    snew(vec[n], natm);
 }

/* read the eigenvalues and eigenvectors in the ED subspace */
 for (n=0; n<ned; n++)
 {
    fgets2(str, STRLEN, ival);
    sscanf(str, "%d%lf", &idum, &val[n]);
    read_x(ivec, natm, vec[n]);  // #natm of eigenvectors were read from ivec for each ED
 }
// By closing this loop, vec is 3 X (natm X ned)
// This means, eigenvector for each ED can be projected to every atom (natm)'s "deviation from average position" in order to visualize the "essential dynamics"

/* obtain the contributions of the essential subspace */
 for (i=0; i<natm; i++)
 {
    for (j=i; j<natm; j++)
    {
       cov = 0.0;
       for (m=0; m<DIM; m++)
       {
          for (n=0; n<ned; n++)
          {
             cov += vec[n][i][m]*val[n]*vec[n][j][m];
          }
       }
       fprintf(oced, "%15.9lf\n", cov);
       fflush(oced);
    }
 }

/* close file */
 fclose(ival);
 fclose(ivec);
 fclose(ipar);
 fclose(oced);

/* the end of main */
 return 0;
}
