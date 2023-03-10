/*
    Common utility routines for the NewtonLib package.
    
 *  Written by        L. Weimann 
 *  Purpose           Performing certain common tasks of NewtonLib codes,
                      part itlin.
 *  Category          ???. - Utilities
 *  Keywords          Memory allocation, scaled norm, scaled scalarproduct,
                      data output, parameter check.
 *  Version           1.0
 *  Revision          May 2006
 *  Latest Change     May 2006
 *  Library           NewtonLib
 *  Code              C, Double Precision
 *  Environment       Standard C environment on PC's,
                      workstations and hosts.
 *  Copyright     (c) Konrad-Zuse-Zentrum fuer
                      Informationstechnik Berlin (ZIB)
                      Takustrasse 7, D-14195 Berlin-Dahlem
                      phone : + 49/30/84185-0
                      fax   : + 49/30/84185-125
 *  Contact           Lutz Weimann
                      ZIB, Division Scientific Computing, 
                           Department Numerical Analysis and Modelling
                      phone : + 49/30/84185-185
                      fax   : + 49/30/84185-107
                      e-mail: weimann@zib.de
 
   ---------------------------------------------------------------
 
 * Licence
     You may use or modify this code for your own non commercial
     purposes for an unlimited time. 
     In any case you should not deliver this code without a special 
     permission of ZIB.
     In case you intend to use the code commercially, we oblige you
     to sign an according licence agreement with ZIB.
 
 * Warranty 
     This code has been tested up to a certain level. Defects and
     weaknesses, which may be included in the code, do not establish
     any warranties by ZIB. ZIB does not take over any liabilities
     which may follow from acquisition or application of this code.
 
 * Software status 
     This code is under care of ZIB and belongs to ZIB software class 2.
 
      ------------------------------------------------------------
 
    
    This file contains the following routines and functions:
    --------------------------------------------------------
    
    zibnum_fwalloc - allocate memory for a double precision array
    zibnum_iwalloc - allocate memory for an integer array
    zibnum_pfwalloc - allocate memory for an array of pointers to 
                    double precision arrays
    zibnum_scaled_norm2 - compute the scaled Euclidian-norm of a vector
    zibnum_scaled_sprod - compute a scaled scalar product
    zibnum_norm2 - compute the unscaled Euclidian-norm of a vector
    zibnum_scale - compute a scaled vector from an unscaled vector 
    zibnum_descale - compute the unscaled vector from a scaled vector
    zibnum_dataout - write data output
    zibnum_parcheck_and_print - check and print parameter and options settings
      
    The calls of the routines/functions are as follows:
    ---------------------------------------------------
    
    int zibnum_fwalloc(int size, double **ptr, char vname[])
    
    This function allocates memory for a double array via the malloc function. 
    The parameters are:
    int    size    (input)  The number of elements of type double to be allocated.
    double **ptr   (output) The pointer to memory, of type (double *), which has
                            been returned by the malloc function.
    char   vname[] (input)  An identifying character string of the memory portion
                            to be allocated, used for print within a possible error 
                            message. 
    The int return code is 0, if the memory allocation was successful, or -999
    if the allocation failed.
    ---
    int zibnum_iwalloc(int size, int **ptr, char vname[])
    
    This function allocates memory for a int array via the malloc function. 
    The parameters are:
    int    size    (input)  The number of elements of type int to be allocated.
    int    **ptr   (output) The pointer to memory, of type (int *), which has
                            been returned by the malloc function.
    char   vname[] (input)  An identifying character string of the memory portion
                            to be allocated, used for print within a possible error 
                            message. 
    The int return code is 0, if the memory allocation was successful, or -998
    if the allocation failed.
    ---
    int zibnum_pfwalloc(int size, double ***ptr, char vname[])
    
    This function allocates memory for a pointer array via the malloc function.
    The parameters are:
    int    size    (input)  The number of elements of type pointer to be allocated.
    int    ***ptr  (output) The pointer to memory, of type (double **) (i.e.
                            to pointers which are pointing to memory allocated 
                            for double arrays), which has been returned by the
                            malloc function.
    char   vname[] (input)  An identifying character string of the memory portion
                            to be allocated, used for print within a possible error 
                            message. 
    The int return code is 0, if the memory allocation was successful, or -997
    if the allocation failed.
    ---
    Note: In order to activate some debug output on dynamic memory allocation
          set the C preprocessor flag DEBUG, i.e. set the option -DDEBUG if
          you use the GNU C-compiler (gcc).
    ---
    double zibnum_scaled_norm2(int n, double *v, double *scale)
    
    This function computes the scaled norm of the (double) vector v of
    size n, using the (double) vector scale for scaling, as below:
    result := Sqrt ( ( Sum(i=0 to n-1) (v[i]/scale[i])^2 ) / n ) .
    ---
    double zibnum_scaled_sprod(int n, double *v1, double *v2, double *scale)
    
    This function computes the scaled scalar product of the (double) vectors
    v1 and v2 of size n, using the (double) vector scale for scaling, as below:
    result := ( Sum(i=0 to n-1) (v1[i]/scale[i])*(v2[i]/scale[i]) ) / n .
    ---
    double zibnum_norm2(int n, double *v)
    
    This function computes the (ordinary) norm of the (double) vector v of
    size n, as below:
    result := Sqrt ( ( Sum(i=0 to n-1) v[i]^2 ) / n ) .
    ---
    void zibnum_scale(int n, double *v1, double *v2, double *scale)
    
    This routine computes the scaled vector of the vector v1 of size n and
    stores the result to the vector v2, as below:
    v2[i] = v1[i]/scale[i]  for i=0,...,n-1 .
    ---
    void zibnum_descale(int n, double *v1, double *v2, double *scale)
    
    This routine computes the descaled vector of the vector v1 of size n and
    stores the result to the vector v2, as below:
    v2[i] = v1[i]*scale[i]  for i=0,...,n-1 .
    ---
    void itlin_dataout(int k, int n, double *x, struct ITLIN_DATA *data)
    
    This routine is designed to print out computed data within each iteration 
    step. It may be replaced or extended for special purposes.
    The parameters are (all input arguments):

    int    k    The current iteration step number.
    int    n    The size of any double arrays mentioned below.
    double *x   An array which holds the current iterate x^k .

    The fields of the struct ITLIN_DATA are:

    double* data->res      An array which holds the residuum r=b-A*x^k.
    double  data->residuum 
    double  data->normdx   The norm of the latest solution correction.
                           Not used by GMRES.
    double  data->tau      The GBIT parameter tau or the gamma_sum of PCG
                           which is used in the PCG termination criterium.
                           Not used by GMRES.
    double  data->t        The GBIT parameter t or the quantity gamma_i^2 of PCG.
                           Not used by GMRES.
    enum    data->mode     The mode of the current iterate:
                           Initial (=1): This is the first call of itlin_dataout.
                           Intermediate (=2): This is an intermediate call of
                                              itlin_dataout.
                           Solution (=3): This is a final call of itlin_dataout,
                                          and *x holds an approximate solution.
                           Final (=4) : This is a final call of itlin_dataout,
                                        but *x does not hold a solution.
    ---
    int itlin_parcheck_and_print(int n, MATVEC *matvec,
                                 struct ITLIN_OPT *opt, int itlin_code)
                                
    This function checks and prints parameter and options settings.
    The parameters are:
    
    int             n      The number of equations and unknowns.
    struct ITLIN_OPT *opt  A pointer to an options data structure.
    int        itlin_code  The identification number of the calling code.
                           Valid identifications are:
                           0 : GMRES
                           1 : GBIT
                           2 : PCG

*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "itlin.h"

struct ITLIN_IO *itlin_ioctl = NULL;

int zibnum_fwalloc(int size, double **ptr, char vname[])
{  int i;
#ifdef DEBUG
   fprintf(stdout,"\n allocation of %s follows, size=%i\n",vname,size);
#endif
   *ptr = (double *) malloc((int) size*sizeof(double)) ;
#ifdef DEBUG
   fprintf(stdout,"\n allocation of %s done\n",vname);
#endif
   if (*ptr == NULL) 
     { if (ERRORLEVEL>0 && ERROR) 
         fprintf(ERROR,"\n allocation of %s failed!\n",vname);
       return -999;}
   else
     {for(i=0;i<size;i++) (*ptr)[i]=0.0; return 0;};
}

int zibnum_iwalloc(int size, int **ptr, char vname[])
{  int i;
#ifdef DEBUG
   fprintf(stdout,"\n allocation of %s follows, size=%i\n",vname,size);
#endif
   *ptr = (int *) malloc((int) size*sizeof(int)) ;
#ifdef DEBUG
   fprintf(stdout,"\n allocation of %s done\n",vname);
#endif
   if (*ptr == NULL) 
     { if (ERRORLEVEL>0 && ERROR) 
         fprintf(ERROR,"\n allocation of %s failed!\n",vname);
       return -998;}
   else
     {for(i=0;i<size;i++) (*ptr)[i]=0; return 0;};
}

int zibnum_pfwalloc(int size, double ***ptr, char vname[])
{  int i;
#ifdef DEBUG
   fprintf(stdout,"\n allocation of %s follows, size=%i\n",vname,size);
#endif
   *ptr = (double **) malloc((int) size*sizeof(double*)) ;
#ifdef DEBUG
   fprintf(stdout,"\n allocation of %s done\n",vname);
#endif
   if (*ptr == NULL) 
     { if (ERRORLEVEL>0 && ERROR) 
         fprintf(ERROR,"\n allocation of %s failed!\n",vname);
       return -997;}
   else
     {for(i=0;i<size;i++) (*ptr)[i]=NULL; return 0;};
}

double zibnum_scaled_norm2(int n, double *v, double *scale)
{  int i;
   double t, rval = 0.0;
   for (i=0;i<n;i++) {t=v[i]/scale[i]; rval += t*t;};
   return sqrt( rval / (double)n );
}

double zibnum_scaled_sprod(int n, double *v1, double *v2, double *scale)
{  int i;
   double t1, t2, rval = 0.0;
   for (i=0;i<n;i++) 
     {t1=v1[i]/scale[i]; t2=v2[i]/scale[i]; rval += t1*t2;};
   return rval / (double)n ;
}

double zibnum_norm2(int n, double *v)
{  int i;
   double rval = 0.0;
   for (i=0;i<n;i++) rval += v[i]*v[i];
   return sqrt( rval / (double)n );
}

void itlin_noprecon(int n, double *x, double *z)
{ int j;
  for (j=0;j<n;j++) z[j]=x[j];   return;
}

void itlin_dataout(int k, int n, double *x, struct ITLIN_DATA *data)
{  int i;
   double *res = data->res;
   if (FITER)
     { fprintf(FITER,"%5i",k);
       for (i=0;i<n;i++) fprintf(FITER,"  % 14.10e",x[i]);
       fprintf(FITER,"\n");
     };
   if (FRES)
     { fprintf(FRES,"%5i",k);
       for (i=0;i<n;i++) fprintf(FRES,"  % 14.10e",res[i]);
       fprintf(FRES,"\n");
     };
   if (FMISC)
     fprintf(FMISC,"%5i  %1i  % 14.10e  % 14.10e  % 14.10e  % 14.10e\n",
                   k,data->codeid, data->residuum, data->normdx, data->tau,
                   data->t);
   if ( DATALEVEL==0 ) return;
   if ( k==0 )
     { fprintf(DATA,"  Start data:\n  N = %i\n\n",n);
       fprintf(DATA,"  Format: iteration-number, (x(i),i=1,...N), Residuum , Normdx\n\n");
       fprintf(DATA,"    Initial data:\n\n");
     }
   else if ( data->mode==Solution )
      fprintf(DATA,"\n    Solution data:\n\n"); 
   else if ( data->mode==Final )
      fprintf(DATA,"\n    Final data:\n\n"); 
   else if ( k==1 && DATALEVEL>1 )
      fprintf(DATA,"\n    Intermediate data:\n\n"); 
   if ( k==0 || data->mode==Solution || data->mode==Final || DATALEVEL > 1 )
     { fprintf(DATA,"  %4i\n",k);
       for (i=0;i<n-2;i+=3)
         fprintf(DATA,"             % 14.10e  % 14.10e  % 14.10e\n",
                      x[i],x[i+1],x[i+2]);
       if (i<n)
         { fprintf(DATA,"             % 14.10e",x[i]); i++;
           if (i<n) fprintf(DATA,"  % 14.10e",x[i]);
           fprintf(DATA,"\n");
         };
       fprintf(DATA,"             % 14.10e  % 14.10e\n",data->residuum,data->normdx);
     };
}

int itlin_parcheck_and_print(int n, MATVEC *matvec,
                             struct ITLIN_OPT *opt, int itlin_code)
#define TOLMIN 1.0e-15
#define TOLMAX 0.5
{  int i, fail=0;
   double large = 1.0/SMALL, default_scale;
   if (!matvec)
     { if ( ERRORLEVEL>0 )
         fprintf(ERROR,"\n Error - Routine matvec must be supplied\n");
       return -99;
     };
   if ( n<=0 ) 
     { if ( ERRORLEVEL>0 )
         fprintf(ERROR,"\n Error - Number of Equations/unknowns must be >0\n");
       return 20;
     };
   if ( opt->tol <= 0 )
     { if ( ERRORLEVEL>0 )
         fprintf(ERROR,"\n Error - opt->tol must be positive\n");
       return 21;
     }
   else
     { if ( opt->tol > TOLMAX ) 
         {  opt->tol = TOLMAX;
            if ( ERRORLEVEL>1 )
              fprintf(ERROR,
              "\n User prescribed RTOL decreased to reasonable largest value RTOL=%e\n",
              opt->tol);
         }
       else if ( opt->tol < TOLMIN ) 
         { opt->tol = TOLMIN;
           if ( ERRORLEVEL>1 )
              fprintf(ERROR,
              "\n User prescribed RTOL increased to reasonable smallest value RTOL=%e\n",
              opt->tol);
         };
      };
   default_scale = opt->tol;
   if ( opt->scale && itlin_code == 1 )
     { for (i=0;i<n;i++) 
         { if ( (opt->scale)[i] < 0.0 ) 
             { if ( ERRORLEVEL>0 )
                 fprintf(ERROR,
                   "\n Error - negative value (opt->scale)[%i] supplied\n",i);
               return 22;
             }
           else if ( (opt->scale)[i] == 0.0 ) (opt->scale)[i] = default_scale;
           else if ( (opt->scale)[i] < SMALL )
             { if ( ERRORLEVEL>1 )
                 fprintf(ERROR,
                 "\n Warning: (opt->scale)[%i] too small - increased to %e\n",
                 i,SMALL);
               (opt->scale)[i] = SMALL;
             }
           else if ( (opt->scale)[i] > large )
             { if ( ERRORLEVEL>1 )
                 fprintf(ERROR,
                 "\n Warning: (opt->scale)[%i] too large - decreased to %e\n",
                 i,large);
               (opt->scale)[i] = large;
             };
         };
     };
   if ( itlin_code == 0 )
     {
       if ( opt->i_max <= 0 ) 
         { opt->i_max = MIN(10,n);
           if ( ERRORLEVEL > 1 ) 
             fprintf(ERROR,
             " Warning: opt->i_max not set; is reset to i_max=%i\n",opt->i_max);
         }
       else if ( opt->i_max > n ) 
         { opt->i_max = MIN(10,n);
           if ( ERRORLEVEL > 1 ) 
             fprintf(ERROR,
             " Warning: opt->i_max is greater than n; is reset to i_max=%i\n",opt->i_max);
         };
     }
   else if ( itlin_code == 1 )
     {
       if ( opt->i_max <= 0 ) 
         { opt->i_max = 10;
           if ( ERRORLEVEL > 1 ) 
             fprintf(ERROR,
             " Warning: opt->i_max not set; is reset to i_max=%i\n",opt->i_max);
         };
       if ( opt->rho <= 0.0 )
         { opt->rho = 1.0;
           if ( ERRORLEVEL > 1 ) 
             fprintf(ERROR,
             " Warning: opt->rho not set; is reset to rho=%e\n",opt->rho);
         };
     };
   if ( MONITORLEVEL==0 ) return 0;
   fprintf(MONITOR,"\n Problem dimension: n = %i\n",n);
   if ( itlin_code == 0 )
      fprintf(MONITOR,"\n Prescribed residuum threshold: %e\n",opt->tol);
   else if ( itlin_code == 1 )
      fprintf(MONITOR,"\n Prescribed relative precision: %e\n",opt->tol);
   fprintf(MONITOR," The maximum permitted number of iteration steps is: %i\n",
                   opt->maxiter);
   if ( itlin_code == GMRES || itlin_code == GBIT )
     fprintf(MONITOR," The maximum number of iterations before restart is: %i\n",
                     opt->i_max);
   if ( itlin_code == 1 )
     {
       fprintf(MONITOR," The safety factor is rho = %6e\n",opt->rho);                
     }
   else if ( itlin_code == 2 )
     {
       if ( opt->convcheck == Absolute )
         fprintf(MONITOR," The absolute error will be checked for termination\n");                
       else if ( opt->convcheck == Relative )
         fprintf(MONITOR," The relative error will be checked for termination\n");
     };
   return 0;
}
