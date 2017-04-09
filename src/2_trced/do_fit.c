/*
 * This program, which modified from GROMACS source codes, can be used to
 * do least square fit. By Zhiyong Zhang, 01/29/2004.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "pca.h"

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);

void reset_x(int natm, rvec x[], double mass[])
{
 int i, m;
 double xcm[DIM], tm;

 tm = 0.0;
 for (m=0; m<DIM; m++)
    xcm[m] = 0.0;

 for (i=0; i<natm; i++) {
    tm += mass[i];
    for (m=0; m<DIM; m++)
       xcm[m] += mass[i]*x[i][m];
 }

 for (m=0; m<DIM; m++)
    xcm[m] /= tm;

 for (i=0; i<natm; i++)
    for (m=0; m<DIM; m++)
       x[i][m] -= xcm[m];
}

double cal_rms(int natm, double mass[], rvec x[], rvec xp[])
{
 int i, m;
 double tm, rs;

 tm = 0.0;
 rs = 0.0;
 for (i=0; i<natm; i++) {
    tm += mass[i];
    for (m=0; m<DIM; m++)
       rs += mass[i]*sqr(x[i][m]-xp[i][m]);
 }
 return sqrt(rs/tm);
}

void jacobi(double **a,int n,double d[],double **v,int *nrot)
{
  int j, i;
  int iq, ip;
  double tresh, theta, tau, t, sm, s, h, g, c, *b, *z;

  snew(b,n);
  snew(z,n);
  for (ip=0; ip<n; ip++) {
     for (iq=0; iq<n; iq++) v[ip][iq] = 0.0;
     v[ip][ip] = 1.0;
  }
  for (ip=0; ip<n;ip++) {
     b[ip] = d[ip] = a[ip][ip];
     z[ip] = 0.0;
  }
  *nrot = 0;
  for (i=1; i<=50; i++) {
     sm = 0.0;
     for (ip=0; ip<n-1; ip++) {
        for (iq=ip+1; iq<n; iq++)
           sm += fabs(a[ip][iq]);
     }
     if (sm == 0.0) {
       sfree(z);
       sfree(b);
       return;
     }
     if (i < 4)
       tresh = 0.2*sm/(n*n);
     else
       tresh = 0.0;
     for (ip=0; ip<n-1; ip++) {
        for (iq=ip+1; iq<n; iq++) {
           g = 100.0*fabs(a[ip][iq]);
           if (i > 4 && fabs(d[ip])+g == fabs(d[ip])
               && fabs(d[iq])+g == fabs(d[iq]))
             a[ip][iq] = 0.0;
           else if (fabs(a[ip][iq]) > tresh) {
             h = d[iq]-d[ip];
             if (fabs(h)+g == fabs(h))
               t = (a[ip][iq])/h;
             else {
               theta = 0.5*h/(a[ip][iq]);
               t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
               if (theta < 0.0) t = -t;
             }
             c = 1.0/sqrt(1+t*t);
             s = t*c;
             tau = s/(1.0+c);
             h = t*a[ip][iq];
             z[ip] -= h;
             z[iq] += h;
             d[ip] -= h;
             d[iq] += h;
             a[ip][iq] = 0.0;
             for (j=0; j<ip; j++) {
                ROTATE(a,j,ip,j,iq)
             }
             for (j=ip+1; j<iq; j++) {
                ROTATE(a,ip,j,j,iq)
             }
             for (j=iq+1; j<n; j++) {
                ROTATE(a,ip,j,iq,j)
             }
             for (j=0; j<n; j++) {
                ROTATE(v,j,ip,j,iq)
             }
             ++(*nrot);
           }
        }
     }
     for (ip=0; ip<n; ip++) {
        b[ip] +=  z[ip];
        d[ip] = b[ip];
        z[ip] =  0.0;
     }
  }
}

void calc_fit_R(int natm, double *w_rls, rvec *xp, rvec *x, matrix R)
{
 int c, r, n, j, m, i, irot;
 static double **omega=NULL, **om=NULL;
 double d[2*DIM], xnr, xpc;
 matrix vh, vk, u;
 double mn;
 int index;
 double max_d;

 if (omega == NULL) {
   snew(omega, 2*DIM);
   snew(om, 2*DIM);
   for (i=0; i<2*DIM; i++) {
      snew(omega[i], 2*DIM);
      snew(om[i], 2*DIM);
   }
 }

 for (i=0; i<2*DIM; i++) {
    d[i] = 0.0;
    for (j=0; j<2*DIM; j++) {
       omega[i][j] = 0.0;
       om[i][j] = 0.0;
    }
 }

/* calculate the matrix U */

 clear_mat(u);
 for (n=0; n<natm; n++)
    for (c=0; c<DIM; c++) {
       xpc = xp[n][c];
       for (r=0; r<DIM; r++) {
          xnr = x[n][r];
          u[c][r] += w_rls[n]*xnr*xpc;
       }
    }

/* construct omega, omega is symmetric */

 for (r=0; r<2*DIM; r++)
    for (c=0; c<=r; c++)
       if (r>=DIM && c<DIM) {
         omega[r][c] = u[r-DIM][c];
         omega[c][r] = u[r-DIM][c];
       }
       else {
         omega[r][c] = 0.0;
         omega[c][r] = 0.0;
       }

/* determine h and k */

 jacobi(omega, 2*DIM, d, om, &irot);

 index = 0;

/* copy the first two eigenvectors */

 for (j=0; j<2; j++) {
    max_d = -1000;
    for (i=0; i<2*DIM; i++)
       if (d[i]>max_d) {
         max_d = d[i];
         index = i;
       }
    d[index] = -10000;
    for (i=0; i<DIM; i++) {
       vh[j][i] = M_SQRT2*om[i][index];
       vk[j][i] = M_SQRT2*om[i+DIM][index];
    }
 }

/* Calculate the last eigenvector as the outer-product of the first two.
 * This insures that the conformation is not mirrored and prevents problems
 * with completely flat reference structures.
 */

 oprod(vh[0],vh[1],vh[2]);
 oprod(vk[0],vk[1],vk[2]);

/* determine R */

 for (r=0; r<DIM; r++)
    for (c=0; c<DIM; c++)
       R[r][c] = vk[0][r]*vh[0][c] + vk[1][r]*vh[1][c] + vk[2][r]*vh[2][c];

}

void do_fit(int natm, double *w_rls, rvec *xp, rvec *x)
{
 int j, m, r, c;
 matrix R;
 rvec x_old;

/* calculate the rotation matrix R */

 calc_fit_R(natm, w_rls, xp, x, R);

/* rotate x */

 for (j=0; j<natm; j++) {
    for (m=0; m<DIM; m++)
       x_old[m] = x[j][m];
    for (r=0; r<DIM; r++) {
       x[j][r] = 0.0;
       for (c=0; c<DIM; c++)
          x[j][r] += R[r][c]*x_old[c];
    }
 }
}
