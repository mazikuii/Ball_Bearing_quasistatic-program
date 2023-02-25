/*
        GBIT - Good Broyden - Error based iterative linear solver

 *  Written by        L. Weimann 
 *  Purpose           Iterative solution of large linear systems
 *  Category          ???. - Linear systems
 *  Keywords          large linear system, iterative solver
 *  Version           1.0
 *  Revision          June 2006
 *  Latest Change     June 2006
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
 
 *    Parameters description
      ======================
 
      The calling interface looks as follows:

      extern void gbit(int n, double *y,
                       MATVEC *matvec, PRECON *precon, double *b,
                       struct ITLIN_OPT *opt, struct ITLIN_INFO *info)

      The structures used within the parameter list are defined
      as follows:
      ---
      struct ITLIN_OPT
      {
         double tol, rho;
         int i_max, maxiter;
         TERM_CHECK termcheck;   /* GMRES only * /
         CONV_CHECK convcheck;   /* PCG only   * /
         LOGICAL rescale;        /* GBIT only  * /
         PRINT_LEVEL errorlevel, monitorlevel, datalevel;
         FILE *errorfile, *monitorfile, *datafile,
              *iterfile, *resfile, *miscfile;
         double *scale;
      };
      
      where the applicable types used within this structure are defined by
      typedef enum {None=0, Minimum=1, Verbose=2, Debug=3} PRINT_LEVEL;
      typedef enum {False=0, True=1} LOGICAL ;
      ---
	  struct ITLIN_INFO
	  {
		 double precision, normdx, residuum;
		 int iter, rcode, subcode, nomatvec, noprecon, noprecl, noprecr;
	  };
      ---
      
      A detailed description of the parameters follows: 
      
      int n :
      The number of equations and unknown variables of the linear system.
      
      double *y :
      A pointer to an array of double values of size n.
      The array must contain on input an initial guess of the linear system
      solution, which is used as the start-vector of the iteration.
      On output, the pointed array contains an approximate solution vector y*,
      which fits the relative error condition epsilon <= rho*||y^i||*opt->tol,
      where ||y||_A = sqrt( sum_1 to n( (y_i/scale_i)^2 ) / n ) denotes a
      scaled Euclidian norm, epsilon is an estimate of the quantity 
      ||y^solution-y^i||, and rho <= 1.0 denotes a safety factor.
      
      void *matvec(int n, double *y, double *z);
      A pointer to the matrix times vector multiplication user routine.
      This routine is required - no default routine is supplied.
      The parameters are:
        int     n     input  Number of vector components.
        double *y     input  Vector of unknowns, of size n .
        double *z     output Vector which holds the matrix-vector product A*y.
 
      void *precon(int n, double *z, double *w);
      A pointer to the preconditioner user routine, which computes w=B*z,
      where B should be an approximation of the inverse of the matrix A.
      If a null pointer is supplied, then a dummy preconditioner routine
      will be used which realizes the preconditioner matrix B=identity.
        int     n     input  Number of vector components.
        double *z     input  Residual vector, of size n .
        double *w     output Vector which holds the matrix-vector product B*z,
                      i.e. the preconditioned residuum.
      
      double *b :
      A pointer to an array of double values of size n.
      The pointed array must hold the right hand side of the linear system 
      to solve.
                                 
      struct ITLIN_OPT *opt:
      A pointer to an options structure. The pointed fields of the structure
      contain input options to gbit.
      
      opt->tol is of type double and must contain the error threshold
      which the relative error of the final iterate y^i must fit.
      
      opt->rho is of type double and must contain a safety factor <= 1.0 by
      which opt->tol will be multiplied.
      
      opt->maxiter is of type int and must contain the maximum number of allowed
      iterations. if a zero or negative value is supplied, then the maximum
      iteration count will be set to 100.
      
      opt->i_max is of type int and must contain the maximum number of 
      iterations before a restart occurs. The main portion of used memory
      by GBIT depends on i_max, such that n*i_max elements of double
      storage will be used. If a nonpositive value is supplied, then i_max
      is set to 10.
      
      opt->scale is a pointer to a double array of size n, which must hold
      nonnegative scaling threshold values. If zero values are supplied,
      the zeros are replaced by the (possibly adjusted) value of opt->tol.
      If a zero pointer is supplied, then all scaling threshold components
      are set to the value of opt->tol. 

      opt->rescale is of type LOGICAL.
      If set to True, the scaling values, which are used in the norm 
      computation, will be adjusted at the start and at each restart
      to the current iterate.
      If set to False, the initial scaling values are used unmodified 
      throughout the whole iteration.
      
      opt->errorlevel is of type PRINT_LEVEL. If it is set to level None,
      then no error message will be printed if an error condition occurs.
      If it is set to level Minimum or any higher level, then error messages
      will be printed, if appropriate.
      
      opt->monitorlevel is of type PRINT_LEVEL. If it is set to level None,
      then no monitor output will be printed.
      If it is set to level Minimum, a few infomation will be printed.
      If set to level Verbose, then some infomation about each iteration
      step, fitting into a single line, will be printed. The higher level Debug
      is reserved for future additional information output.
      
      opt->datalevel is of type PRINT_LEVEL. If it is set to level None,
      then no data output will be printed.
      If it is set to level Minimum, then the values of the initial iteration
      vector x and the final vector x will be printed.
      If set to level Verbose, then the iteration vector x will be printed for
      each step. The higher level Debug is reserved for future additional
      information output.
      
      opt->errorfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer, stderr, stdout,
      or to another file pointer which has been initialized by a fopen call.
      If it is set to NULL, opt->errorfile will be set to stdout. The error 
      messages will be printed to opt->errorfile.
      
      opt->monitorfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer, stderr, stdout,
      or to another file pointer which has been initialized by a fopen call.
      If it is set to NULL, opt->monitorfile will be set to stdout. The monitor 
      output will be printed to opt->monitorfile.
      
      opt->datafile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer, stderr, stdout,
      or to another file pointer which has been initialized by a fopen call.
      If it is set to NULL, a file named "gbit.data" will be opened by a
      fopen call and opt->datafile will be set to the filepointer which the
      fopen returns. The data output will be printed to opt->datafile.
      
      opt->iterfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer or to file pointer
      which has been initialized by a fopen call. The iteration number and
      the iteration vector will be written out to the associated file, for
      each iteration step. If opt->iterfile is set to NULL, no such 
      data will be written out.
      
      opt->miscfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer or to file pointer
      which has been initialized by a fopen call. The iteration number, an
      identification number of the calling code (1 for GBIT), a floating point
      zero number as a dummy placeholder, the scaled norm of the (unrelaxated)
      correction, the parameter tau, and the relaxation factor ti will be
      written out, for each iteration step. If opt->miscfile is set to NULL,
      no such data will be written out.
     
      Note: The output to the files opt->iterfile and opt->miscfile 
            is written as a single for each iteration step.
            Such, the data in these files are suitable as input to the
            graphics utility GNUPLOT.
      
      struct ITLIN_INFO *info:
      A pointer to an info structure. The pointed fields of the structure
      are set output info of GBIT.
      
      info->precision is of type double and is set to the relative error
      estimate epsilon/(rho*||y^i||). For detailed info, refer to the 
      description of the parameter y above.
      
      info->iter is set to number of iteration steps done.
      
      info->nomatvec is set to the number of done calls to the matrix times
      vector multiplication user routine matvec. 
      
      info->noprecon is set to the number of done calls to the preconditioner
      user routine precon or the dummy preconditioner routine, if the user
      didn't supply a preconditioner routine.
      
      info->rcode is set to the return-code of GBIT. A return-code 0
      means that PCG has terminated sucessfully. For the meaning of a
      nonzero return-code, see the error messages list below.


      The following error conditions may occur: (returned via info->rcode)
      --------------------------------------------------------------------
      
      -999 routine zibnum_fwalloc failed to allocate double memory via malloc.
      -997 routine zibnum_pfwalloc failed to allocate double pointer memory 
           via malloc.
      -995 Internal i/o control block could not be allocated via malloc.
      -994 Internally used data structure could not be allocated via malloc.
      -989 Default data-output file could not be opened via fopen call.
       -99 NULL pointer supplied for matvec - the matrix times vector routine
           must be defined!
         2 Maximum number of iterations (as set by opt->maxiter) exceeded.
         3 No convergence, error estimate epsilon exceeded the limit 1.0e20. 
        20 Nonpositive input for dimensional parameter n.
        21 Nonpositive value for opt->tol supplied.
        22 Negative scaling value for some component of vector opt->scale
           supplied.
         
      Summary of changes:
      -------------------
      
      Version  Date        Changes
      1.0      2006/11/30  Initial Release.
      
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "itlin.h"

#define TAU_MIN 1.0e-8
#define TAU_MAX 1.0e2
#define EPSILON_LIMIT 1.0e20

#define MAXITER_DEFAULT 100

extern struct ITLIN_IO *itlin_ioctl;

int gbit_initscale(int n, double **scale, struct ITLIN_OPT *opt);
void gbit_rescale(int n, double *scale, double *y, double *ythresh);

extern void gbit(int n, double *y, MATVEC *matvec, PRECON *precon, double *b,
          struct ITLIN_OPT *opt, struct ITLIN_INFO *info)
{ 
   int i, j, m, iter=0, nomatvec=0, noprecon=0,
       i_max, max_iter = opt->maxiter;
   double **delta, *sigma, *t, *q, *z, *ythresh, *y0,
          *zeta, *delta0, *deltam, *deltamp1,
          *deltai, *deltaip1, *scale = opt->scale;
   double gamma, tau, fac, tm1, ti, epsilon, normyip1, normdiff,
          rho, errtol=opt->tol;
   LOGICAL stop_iter, io_allocated=False, 
           giant_simple=opt->giant_simple;
   struct ITLIN_DATA *data=malloc(sizeof(struct ITLIN_DATA));

   if (!itlin_ioctl) itlin_ioctl=malloc(sizeof(struct ITLIN_IO));
   if (!itlin_ioctl) 
     { fprintf(stderr,"\n could not allocate output controlblock\n");
       RCODE=-995; return; }
   else
     io_allocated = True;
   if (!data)
     { fprintf(stderr,"\n could not allocate struct data\n");
       RCODE=-994; return; };
   data->codeid    = GBIT;
   data->res       = NULL;
   data->residuum  = 0.0;
   data->mode      = Initial;
   ERRORLEVEL   = opt->errorlevel;
   MONITORLEVEL = opt->monitorlevel;
   DATALEVEL    = opt->datalevel;
   ERROR    = opt->errorfile;
   MONITOR  = opt->monitorfile;
   DATA     = opt->datafile;
   FITER    = opt->iterfile;
   FRES     = NULL;
   FMISC    = opt->miscfile;
   if ( !ERROR && ERRORLEVEL>0 )     ERROR   = stdout;
   if ( !MONITOR && MONITORLEVEL>0 ) MONITOR = stdout;
   if ( !DATA && DATALEVEL>0 )
     { DATA=fopen("gbit.data","w");
       if (!DATA && ERRORLEVEL>0)
         { fprintf(ERROR,"\n fopen of file gbit.data failed\n");
           RCODE=-989; return;
         };
     };
   opt->errorfile   = ERROR;
   opt->monitorfile = MONITOR;
   opt->datafile    = DATA;
   if ( max_iter <= 0 ) max_iter = MAXITER_DEFAULT; 
   if ( MONITORLEVEL > 0 ) fprintf(MONITOR,"\n GBIT - Version  1.0\n");
   RCODE = itlin_parcheck_and_print(n,matvec,opt,1);
   if ( RCODE !=0 ) 
     { if (io_allocated) {free(itlin_ioctl); itlin_ioctl=NULL;};
       if (data) free(data);
       return;
     };
   rho   = opt->rho;
   i_max = opt->i_max;
   if (!precon) precon = &itlin_noprecon;
   RCODE = 0;
   RCODE = zibnum_pfwalloc(i_max+2,&delta,"delta");  if ( RCODE !=0 ) return;
   RCODE = zibnum_fwalloc(i_max+2,&sigma,"sigma");   if ( RCODE !=0 ) return;
   RCODE = zibnum_fwalloc(i_max+2,&t,"t");           if ( RCODE !=0 ) return;
   RCODE = zibnum_fwalloc(n,&q,"q");                 if ( RCODE !=0 ) return;
   RCODE = zibnum_fwalloc(n,&z,"z");                 if ( RCODE !=0 ) return;
   RCODE = zibnum_fwalloc(n,&ythresh,"ythresh");     if ( RCODE !=0 ) return;
   RCODE = zibnum_fwalloc(n,&y0,"y0");               if ( RCODE !=0 ) return;
   zeta = z;
   for (j=0;j<n;j++) y0[j]=y[j];
   if ( MONITORLEVEL > 1 )
     fprintf(MONITOR,"\n Iter      sigma   NormIter        tau          t\n\n");
   for (i=0;i<=i_max+1;i++)
     { RCODE = zibnum_fwalloc(n,&delta[i],"delta[i]");    
       if ( RCODE !=0 ) return;
     };
   if (!scale) 
     { RCODE = gbit_initscale(n,&scale,opt); 
       if ( RCODE !=0 ) goto errorexit;
     };  
   for (j=0;j<n;j++) ythresh[j] = scale[j];
   delta0 = delta[0];

   /* Initialization */
   restart:
   if (opt->rescale) gbit_rescale(n,scale,y,ythresh);
   if ( iter > 0 && MONITORLEVEL > 1 ) fprintf(MONITOR," > RESTART\n");
   normyip1 = zibnum_scaled_norm2(n,y,scale);
   matvec(n,y,z); nomatvec++;  for (j=0;j<n;j++) z[j] = b[j]-z[j];
   precon(n,z,delta0);  noprecon++;
   sigma[0] = zibnum_scaled_sprod(n,delta0,delta0,scale);

   /* Iteration loop (restarts excluded) */
   stop_iter = False;
   for (i=0; i<=i_max && !stop_iter && iter <= max_iter; i++)
     { 
       matvec(n,delta[i],q);  nomatvec++;  precon(n,q,zeta);  noprecon++;
       for (m=0;m<=i-1;m++)
         {
           tm1 = 1.0-t[m];  
           fac = zibnum_scaled_sprod(n,delta[m],zeta,scale)/sigma[m];
           deltam = delta[m];  deltamp1 = delta[m+1];
           for (j=0;j<n;j++) zeta[j] += fac*(deltamp1[j]-tm1*deltam[j]);
         };
       gamma = zibnum_scaled_sprod(n,delta[i],z,scale);
       if ( gamma != 0.0 ) tau = sigma[i]/gamma;
       else                tau = 2.0*TAU_MAX;
       if ( tau < TAU_MIN ) 
         { if ( MONITORLEVEL > 1 ) fprintf(MONITOR," > tau = %9e\n",tau);
           if ( i>0 ) goto restart;
           else
             { t[i] = 1.0;
               if ( MONITORLEVEL > 0 )
                 fprintf(MONITOR," > Restart condition ignored, set t[0]=%9e\n",
                                 t[i]);
             };
         }
       else
         if ( tau <= TAU_MAX ) t[i] = tau;  else t[i] = 1.0;
       ti = t[i];  deltai = delta[i];
       if ( MONITORLEVEL > 1 )
         fprintf(MONITOR," %4i % 10.3e % 10.3e % 10.3e % 10.3e\n",
                         iter,sigma[i],normyip1,tau,t[i]);
       data->normdx = sqrt(sigma[i]);  data->tau = tau;  data->t = t[i];
       itlin_dataout(iter,n,y,data);
       data->mode = Intermediate;
       for (j=0;j<n;j++) y[j] += ti*deltai[j];
       deltaip1 = delta[i+1];
       fac = 1.0-ti+tau;
       for (j=0;j<n;j++) deltaip1[j] = fac*deltai[j]-tau*z[j];
       sigma[i+1] = zibnum_scaled_sprod(n,deltaip1,deltaip1,scale);
       if ( i>=1 ) epsilon = 0.5*sqrt(sigma[i-1]+2.0*sigma[i]+sigma[i+1]);
       else        epsilon = sqrt(sigma[1]);
       normyip1 = zibnum_scaled_norm2(n,y,scale);
       if (giant_simple)
         { for (j=0;j<n;j++) q[j]=y[j]-y0[j];
           normdiff = zibnum_scaled_norm2(n,q,scale);
           stop_iter = (epsilon/normdiff <= rho*errtol);
         }
       else 
         stop_iter = epsilon <= rho*normyip1*errtol;
       iter++;
       if ( !stop_iter && epsilon > EPSILON_LIMIT ) 
         { RCODE = 3; goto errorexit; };
     };
   if ( i>=i_max && !stop_iter && iter < max_iter ) goto restart;
   RCODE = ( iter < max_iter ? 0 : 2 );
   
   errorexit:
   if ( MONITORLEVEL > 1 )
     fprintf(MONITOR," %4i % 10.3e % 10.3e % 10.3e % 10.3e\n",
                     iter,sigma[i],normyip1,tau,t[i-1]);
   data->mode = ( RCODE==0 ? Solution : Final );
   data->normdx = sqrt(sigma[i]);  data->tau = tau;  data->t = t[i-1];
   itlin_dataout(iter,n,y,data);
   if ( ERRORLEVEL > 0 && RCODE != 0 )
     {
       switch ( RCODE )
        {
         case     2:
           fprintf(ERROR,"\n Error - Maximum allowed number of iterations exceeded\n");
           break;
         case     3:
           fprintf(ERROR,"\n Error - no convergence, correction too large\n");
           break;
         default   :
           fprintf(ERROR,"\n Error, code=%i\n",RCODE);
        };
     };
   for (i=0;i<=i_max;i++)  free(delta[i]);
   free(delta);  free(sigma); free(t);  free(q);  free(z);  free(y0);
   if (io_allocated) {free(itlin_ioctl); itlin_ioctl=NULL;};
   free(data);   
   info->precision = epsilon/(rho*normyip1);
   info->iter      = iter;
   info->nomatvec  = nomatvec;
   info->noprecon  = noprecon;
}

int gbit_initscale(int n, double **scale, struct ITLIN_OPT *opt)
{  int i, rcode;
   double tol = opt->tol;
   rcode = zibnum_fwalloc(n,scale,"scale");  if ( rcode !=0 ) return rcode;
   for (i=0;i<n;i++) (*scale)[i]=tol;
   return 0;
}

void gbit_rescale(int n, double *scale, double *y, double *ythresh)
{  int i;
   for (i=0;i<n;i++) 
     scale[i] = MAX(ythresh[i],MAX(fabs(y[i]),SMALL));
}
