/*
    GIANT_GBIT - ERRor-based Inexact (damped) Newton method
    
 *  Written by        L. Weimann 
 *  Purpose           Solution of systems of nonlinear equations
 *  Method            ERRor-based Damped Newton algorithm
                      (see reference below)
 *  Category          F2a. - Systems of nonlinear equations
 *  Keywords          Nonlinear equations, Newton methods
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
 
 *    References:
 
      /1/ P. Deuflhard:
          Newton Methods for Nonlinear Problems. -
          Affine Invariance and Adaptive Algorithms.
          Series Computational Mathematics 35, Springer (2004)
 
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

      extern void giant_gbit(struct GIANT_FUN funs, int n, double *x,
                           struct GIANT_OPT *opt, 
                           struct GIANT_INFO *info)

      The structures used within the parameter list are defined
      as follows:
      ---
      struct GIANT_FUN
      {
        GIANT_FFUN *fun;
        GIANT_JFUN *jac;
        MATVEC     *matvec;
        PRECON     *preconr,*preconl;
      };
      
      where the types used within this structure are defined by
      typedef void GIANT_FFUN(int,double*,double*,int,int*); 
      typedef void GIANT_JFUN(int,double*,double*,double*,int,int*); 
      typedef void MATVEC(int, double*, double*);  
      and
      typedef void PRECON(int, double*, double*);
      ---
      struct GIANT_OPT
      {
	   double tol, eta_bar, safetyfactor, rho;
	   int maxiter, lin_maxiter, i_max;
	   LOGICAL restricted; 
	   ADAPT_MODE adaptmode;
	   enum {StandardScale=0,StartValueScale=1} scaleopt;
	   PRINT_LEVEL errorlevel, monitorlevel, datalevel, plotlevel,
				   linmonlevel, lindatalevel;
	   FILE *errorfile, *monitorfile, *datafile, *linmonfile, *lindatafile,
			*iterfile, *resfile, *miscfile;
	   PROBLEM_TYPE nonlin;
	   double *scale;
      };
      
      where the types used within this structure are defined by
      typedef enum {False=0, True=1} LOGICAL ;
      typedef enum {StandardMode=0, QuadraticMode=1} ADAPT_MODE;
      typedef enum {None=0, Minimum=1, Verbose=2, Debug=3} PRINT_LEVEL;
      typedef enum {Mildly_Nonlinear=2, Highly_Nonlinear=3,
                    Extremely_Nonlinear=4} PROBLEM_TYPE ;
      ---
      struct GIANT_INFO
      {
	   double precision, normdx;
	   double *fx, *dx;
	   int iter, rcode, subcode, nofunevals, nojacevals, nomuljac, noprecon,
		   noordlinit, nosimlinit;
      };
      ---
      
      A detailed description of the parameters follows: 
      
      struct GIANT_FUN funs :
      
      The field funs.fun must contain a pointer to the problem function fun -
      The required parameters interface of fun is described in detail below
      
      The field funs.jac must either contain a pointer to the Jacobian function
      jac or a NULL pointer. If a NULL pointer is supplied, then the Jacobian
      will be approximately computed by an internal function of the giant_gbit
      package.
      
      int n :
      The number of equations and unknown variables of the nonlinear system.
      
      double *x :
      A pointer to an array of double values of size n.
      The array must contain on input an initial guess of the problems solution,
      which is used as the start-vector of the damped Newton iteration.
      On output, the pointed array contains an approximate solution vector x^k,
      which fits the error condition
      || x^k - x* || <= opt->tol,
      where x* denotes the exact solution, and ||...|| is a scaled 
      Euclidian norm.
      
      struct GIANT_OPT *opt:
      A pointer to an options structure. The pointed fields of the structure
      contain input options to giant_gbit.
      
      opt->tol is of type double and must contain the error tolerance which
      the final approximate solution x^k must fit.
      
      opt->maxiter is of type int and must contain the maximum number of allowed
      iterations. if a zero or negative value is supplied, then the maximum
      iteration count will be set to 50.
      
      opt->nonlin is of type PROBLEM_TYPE and must classify the problem to
      be solved. The following classifications may be used:
      Mildly_Nonlinear: The problem is considered to be mildly nonlinear and
                        giant_gbit starts up with dampingfactor=1.
      Highly_Nonlinear: The problem is considered to be highly nonlinear and
                        giant_gbit starts up with dampingfactor=1.0e-4.
      Extremely_Nonlinear: The problem is considered to be extremely nonlinear
                        and giant_gbit starts up with dampingfactor=1.0e-6.
                        Moreover, opt->restricted is set automatically to True.

      opt->restricted is of type LOGICAL.
      If set to True, then the restricted monotonicity test will be applied for
      determination whether the next iterate (and the associate damping factor
      lambda) will be accepted. This means, with
      theta = ||dxbar(k+1)|| / ||dx(k)||,
      (where dx(k) denotes the k-th Newton correction, and dxbar(k+1) denotes
       the sucessive simplified Newton correction)
      the condition theta <= 1.0 - lambda/4 must be fit.
      If set to False, then the standard monotonicity test will be applied, i.e.
      the following condition must be fit:
      theta < 1.0.
      
      opt->nleqcalled is of type LOGICAL and only used internally. This field
      should always be set to False.
      
      opt->errorlevel is of type PRINT_LEVEL. If it is set to level None,
      then no error message will be printed if an error condition occurs.
      If it is set to level Minimum or any higher level, then error messages
      will be printed, if appropriate.
      
      opt->monitorlevel is of type PRINT_LEVEL. If it is set to level None,
      then no monitor output will be printed.
      If it is set to level Minimum, a few infomation will be printed.
      If set to level Verbose, then some infomation about each Newton iteration
      step, fitting into a single line, will be printed. The higher level Debug
      is reserved for future additional information output.
      
      opt->datalevel is of type PRINT_LEVEL. If it is set to level None,
      then no data output will be printed.
      If it is set to level Minimum, then the values of the initial iteration
      vector x and the final vector x will be printed.
      If set to level Verbose, then the iteration vector x will be printed for
      each Newton step. The higher level Debug is reserved for future additional
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
      If it is set to NULL, a file named "giant_gbit.data" will be opened by a
      fopen call and opt->datafile will be set to the filepointer which the
      fopen returns. The data output will be printed to opt->datafile.
      
      opt->iterfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer or to file pointer
      which has been initialized by a fopen call. The iteration number and
      the iteration vector will be written out to the associated file, for
      each Newton iteration step. If opt->iterfile is set to NULL, no such 
      data will be written out.
      
      opt->resfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer or to file pointer
      which has been initialized by a fopen call. The iteration number and
      the residuum vector will be written out to the associated file, for
      each Newton iteration step. If opt->resfile is set to NULL, no such 
      data will be written out.
      
      opt->miscfile is of type pointer to FILE, as defined by the <stdio.h>
      header file. It must be set either to a NULL pointer or to file pointer
      which has been initialized by a fopen call. The iteration number, an
      identification number of the calling code (0 for GIANT_GBIT), the norm
      of the residuum, the norm of the Newton correction, the norm of the
      accepted simplified correction, the accepted damping factor, and a
      zero value as a dummy placeholder value will be written out, for
      each Newton iteration step. If opt->miscfile is set to NULL, no such 
      data will be written out. For additional information on the file output,
      refer to the description of this option in the QNERR documentation.
     
      Note: The output to the files opt->iterfile, opt->resfile and
            opt->miscfile is written as a single for each iteration step.
            Such, the data in these files are suitable as input to the
            graphics utility GNUPLOT.
      
      opt->scale is of type pointer to a double array of size n. 
      This array must, if present, contain positive scaling values, which are
      used in computations of scaled norms and Jacobian scaling, as follows:
      || x || = squareroot(sum(1 to n) ( x_i/scale_i )^2)
      The pointer may be initialized with a NULL pointer. In this case, all
      scaling values are internally set to 1.
      
      struct GIANT_INFO *info:
      A pointer to an info structure. The pointed fields of the structure
      are set output info of giant_gbit.
      
      info->precision is of type double and is set to the achieved scaled norm
      of the error of the final iterate.
      
      info->normdx is of type double and is set to the unscaled norm of the
      last Newton correction.
      
      info->fx is a pointer to a double array of size n, which contains the
      final residuum vector.
      
      info->iter is set to number of Newton iteration steps done.
      
      info->nofunevals is set to the number of done calls to the problem
      function funs.fun.
      
      info->nojacevals is set to the number of done calls to the Jacobian
      function funs.jac.
      
      info->rcode is set to the return-code of giant_gbit. A return-code 0
      means that giant_gbit has terminated sucessfully. For the meaning of a
      nonzero return-code, see the error messages list below.
      
      info->subcode is set for certain failure conditions to the error code
      which has been returned by a routine called from giant_gbit.


      Parameter definitions of the required problem routine funs.fun
      and the optional Jacobian routine funs.jac
      --------------------------------------------------------------
      
      void fun(int *n, double *x, double *f, int *fail);
        int    *n     input  Number of vector components.
        double *x     input  Vector of unknowns, of size *n .
        double *f     output Vector of function values.
        int    *fail  output fun evaluation-failure indicator.
                      On input:  undefined.
                      On output: Indicates failure of fun evaluation, 
                                 if having a value <= 2.
                      If <0 or >2: giant_gbit will be terminated with
                                   error code = 82, and *fail will be stored to
                                   info->subcode.
                      If =1: A new trial Newton iterate will
                             computed, with the damping factor
                             reduced to it's half.
                      If =2: A new trial Newton iterate will computed, with the
                             damping factor reduced by a reduction factor, which
                             must be output through f[0] by fun, and it's value
                             must be >0 and <1.
                      Note, that if IFAIL = 1 or 2, additional conditions 
                      concerning the damping factor, e.g. the minimum damping
                      factor may also influence the value of the reduced 
                      damping factor.
 
      void jac(int *n, int *ldjac, double *x, double *dfdx, int *fail);
                        Ext    Jacobian matrix subroutine
        int    *n     input  Number of vector components.
        int    *ldjac input  Leading dimension of Jacobian array, i.e.
                             the total row length for C-style two-dimensional
                             arrays, or the total column length for 
                             Fortran-style two-dimensional arrays.
                             See Note below!
        double *x     input  Vector of unknowns, of size *n .
        double *dfdx  output dfdx[i][k]: partial derivative of i-th component
                             of output parameter *f from fun with respect 
                             to x[k].
        int    *fail  output fun evaluation-failure indicator.
                      On input:  undefined.
                      On output: Indicates failure of jac evaluation
                      and causes termination of giant_gbit, f set to a nonzero
                      value on output.
                                 
      Note: The calling interfaces of the user routines fun and jac has
            been designed to be compatible with routines programmed for
            use with the Fortran codes NLEQ1 and NLEQ2. However, note
            that the Fortran matrix storage mode is columnwise while
            the C matrix storage mode is rowwise. If you intend to link 
            a Jacobian routine, which has been programmed in Fortran for
            use with NLEQ1 or NLEQ2, you must either transpose the Jacobian,
            or you must compile the giant_gbit package for use with Fortran
            matrix storage mode, by setting the C preprocessor flag FMAT,
            i.e. setting the gcc compiler option -DFMAT, when compiling the
            file jacobian_and_linalg.c .


      The following error conditions may occur: (returned via info->rcode)
      --------------------------------------------------------------------
      
      -999 routine zibnum_fwalloc failed to allocate double memory via malloc.
      -995 Internal i/o control block could not be allocated via malloc.
      -994 Internally used data structure could not be allocated via malloc.
      -989 Default data-output file could not be opened via fopen call.
       -99 NULL pointer obtained from funs.fun field - the problem function
           must be defined!
         2 Maximum number of Newton iteration (as set by opt->maxiter) exceeded.
         3 No convergence of Newton iteration, damping factor became too small.
        20 Nonpositive input for dimensional parameter n.
        21 Nonpositive value for opt->tol supplied.
        22 Negative scaling value for some component of vector opt->scale
           supplied.
        81 gbit returned with an error.
           Check info->subcode for the gmres failure code.
        82 The user defined problem function funs.fun returned a nonzero code
           other than 1 or 2. 
           Check info->subcode for the user-function failure code.
        83 The user defined Jacobian function funs.jac returned a nonzero code.
           Check info->subcode for the Jacobian-function failure code.
         
      Summary of changes:
      -------------------
      
      Version  Date        Changes
      1.0      2006/11/30  Initial release.
      
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "giant.h"

int giantgb_initscale(int n, double **scale, struct GIANT_OPT *opt);
void giantgb_rescale(int n, double *scale, double *x, double *xa,
                     double *xthresh);
void giantgb_monitor(int k, int n, double normdx, double normf, double lambda,
                      double tol, double prec, int itlin, char indicator);

#define MAX_ITER_DEFAULT 50
#define LAMBDA_START_DEFAULT 1.0e-2
#define LAMBDA_START_EXTREMELY_DEFAULT 1.0e-4
#define LAMBDA_MIN_DEFAULT 1.0e-4
#define LAMBDA_MIN_EXTREMELY_DEFAULT 1.0e-8
#define DELTA_BAR 1.0e-3

#define NONLIN opt->nonlin

extern struct GIANT_IO *giant_ioctl;

extern void giant_gbit(struct GIANT_FUN fun, int n, double *x,
                  struct GIANT_OPT *opt, struct GIANT_INFO *info)
{  double xtol=opt->tol;
   double lambda, lambda_new, lambda_prev, hk, normfk, normfkp1, normdxkm1, 
          normdxbar, normdiff, reduction_factor, normdx, theta, s, hkpri;
   double lambda_min, rhobar, rhotilde, rhobar_max, rho=opt->rho,
          safetyfact;
   int i, j, k=0, fail=0,  liniter, nordlin=0, nsimlin=0,
            max_iter=opt->maxiter;
   LOGICAL     restricted=opt->restricted,
               io_allocated=False,
               reducted;
   double *dx, *dxbar, *fxk, *xscale, *xthresh, *xkp1, *w;
   int scale_allocated=0, nfcn=0, njac=0, nmuljac=0, nprecon=0;
   ADAPT_MODE adaptmode=opt->adaptmode;
   GIANT_FFUN *f = fun.fun;
   GIANT_JFUN *jac = fun.jac;
   MATVEC     *matvec  = fun.muljac;
   PRECON     *preconr = NULL,
              *preconl = fun.preconl;
   struct GIANT_DATA *data=malloc(sizeof(struct GIANT_DATA));
   struct ITLIN_INFO *infolin = malloc(sizeof(struct ITLIN_INFO));
   struct ITLIN_OPT  *optlin  = malloc(sizeof(struct ITLIN_OPT));
   if (!giant_ioctl) giant_ioctl=malloc(sizeof(struct GIANT_IO));
   if (!giant_ioctl) 
     { fprintf(stderr,"\n could not allocate output controlblock\n");
       RCODE=-995; return; }
   else
     io_allocated = True;
   if (!data)
     { fprintf(stderr,"\n could not allocate struct data\n");
       RCODE=-994; return; };
   data->codeid    = GIANT_GBIT;
   data->theta     = 0.0;
   data->mode      = Initial;
   ERRORLEVEL   = opt->errorlevel;
   MONITORLEVEL = opt->monitorlevel;
   DATALEVEL    = opt->datalevel;
   PLOTLEVEL    = opt->plotlevel;
   ERROR    = opt->errorfile;
   MONITOR  = opt->monitorfile;
   DATA     = opt->datafile;
   FITER    = opt->iterfile;
   FRES     = opt->resfile;
   FMISC    = opt->miscfile;
   if ( !ERROR && ERRORLEVEL>0 )     ERROR   = stdout;
   if ( !MONITOR && MONITORLEVEL>0 ) MONITOR = stdout;
   if ( !DATA && DATALEVEL>0 )
     { DATA=fopen("giant_gbit.data","w");
       if (!DATA && ERRORLEVEL>0)
         { fprintf(ERROR,"\n fopen of file giant_gbit.data failed\n");
           RCODE=-989; return;
         };
     };
   opt->errorfile   = ERROR;
   opt->monitorfile = MONITOR;
   opt->datafile    = DATA;
   if ( MONITORLEVEL > 0 ) fprintf(MONITOR,"\n GIANT_GBIT - Version  1.0\n");
   RCODE = giant_parcheck_and_print(n,opt,fun,0);
   if ( RCODE !=0 ) 
     { if (io_allocated) {free(giant_ioctl); giant_ioctl=NULL;};
       if (data) free(data);
       return;
     };
   if ( max_iter <= 0 ) max_iter = MAX_ITER_DEFAULT;
   if      ( NONLIN==Mildly_Nonlinear ) 
     { lambda = 1.0; lambda_min = LAMBDA_MIN_DEFAULT; }
   else if ( NONLIN==Highly_Nonlinear ) 
     { lambda = LAMBDA_START_DEFAULT; lambda_min = LAMBDA_MIN_DEFAULT; }
   else if ( NONLIN==Extremely_Nonlinear ) 
     { lambda = LAMBDA_START_EXTREMELY_DEFAULT;
       lambda_min = LAMBDA_MIN_EXTREMELY_DEFAULT;
       restricted = True;
     } ;
   RCODE = zibnum_fwalloc(n,&dx,"dx");           if ( RCODE !=0 ) return;
   RCODE = zibnum_fwalloc(n,&dxbar,"dxbar");     if ( RCODE !=0 ) return;
   RCODE = zibnum_fwalloc(n,&xkp1,"xkp1");       if ( RCODE !=0 ) return;
   RCODE = zibnum_fwalloc(n,&fxk,"fxk");         if ( RCODE !=0 ) return;
   RCODE = zibnum_fwalloc(n,&xthresh,"xthresh"); if ( RCODE !=0 ) return;
   RCODE = zibnum_fwalloc(n,&w,"w");             if ( RCODE !=0 ) return;
   optlin->i_max        = opt->i_max;
   if ( opt->i_max <= 0 ) optlin->i_max=10;
   safetyfact=opt->safetyfactor;
   if ( safetyfact <= 0.0 || safetyfact > 1.0 ) safetyfact = 1.0;
   rhotilde = 0.5*rho;  rhobar_max = rhotilde/(1.0+rhotilde);
   optlin->rho          = safetyfact;
   optlin->rescale      = False;
   optlin->maxiter      = opt->lin_maxiter;
   optlin->errorlevel   = opt->linmonlevel;
   optlin->monitorlevel = opt->linmonlevel;
   optlin->datalevel    = opt->lindatalevel;
   optlin->errorfile    = opt->linmonfile;
   optlin->monitorfile  = opt->linmonfile;
   optlin->datafile     = opt->lindatafile;
   optlin->iterfile     = NULL;
   optlin->resfile      = NULL;
   optlin->miscfile     = NULL;
   infolin->iter = 0;
   data->fx = fxk;
   data->dx = dx;
   data->noitersim = 0;
   data->maxiter = max_iter;
   f(n,x,fxk,nfcn,&fail);  nfcn++;
   if (fail != 0) { RCODE=82; goto errorexit ;};
   xscale = opt->scale;
   if ( xscale == NULL ) 
     {RCODE=giantgb_initscale(n,&xscale,opt); if (RCODE !=0) goto errorexit;
      scale_allocated = 1;};
   optlin->scale = xscale;
   if ( MONITORLEVEL > 1 ) 
      fprintf(MONITOR,"\n iter     norm_scl(dx)      norm(fk)  tol_prescr         tol itlin    lambda \n\n");
   normfk = zibnum_norm2(n,fxk);
   normdx = 0.0;  normdxbar = 0.0;
   for (i=0;i<n;i++) { w[i]=x[i]; xthresh[i]=xscale[i]; };
   RCODE = 2;

   do 
     { giantgb_rescale(n,xscale,x,w,xthresh);
       if ( k>0 )
         /* recompute norms after rescaling */
         { normdxkm1 = zibnum_scaled_norm2(n,dx,xscale);
           normdxbar = zibnum_scaled_norm2(n,dxbar,xscale);
         };
       jac(n,x,xscale,fxk,njac,&fail);  njac++;
       if (fail != 0) { RCODE=83; goto errorexit ;};
       for (j=0;j<n;j++) { w[j]=-fxk[j];  dx[j]=dxbar[j]; };
       if ( adaptmode == QuadraticMode )  optlin->tol = 0.25;
       else                               optlin->tol = rho/(2.0*(1.0+rho));
       optlin->giant_simple = False;
       gbit(n,dx,matvec,preconl,w,optlin,infolin);
       fail = infolin->rcode;
       liniter = infolin->iter;
       nmuljac += infolin->nomatvec;  nprecon += infolin->noprecon;
       if ( fail != 0 ) { RCODE=81; goto errorexit; };
       normdx = zibnum_scaled_norm2(n,dx,xscale);
       lambda_prev = lambda;
       if ( k>0 )
         { hkpri = normdx/normdxkm1*hk; 
           lambda = MIN(1.0,1.0/((1.0+rho)*hkpri));
         };
       if ( lambda == 1.0 && k>0 )
         { if ( adaptmode == QuadraticMode )  
             optlin->tol = rho/2.0*hk/(1.0+hk);
           else
             optlin->tol = DELTA_BAR;
           gbit(n,dx,matvec,preconl,w,optlin,infolin);
           fail = infolin->rcode;
           liniter += infolin->iter;
           nmuljac += infolin->nomatvec;  nprecon += infolin->noprecon;
           if ( fail != 0 ) { RCODE=81; goto errorexit; };
           normdx = zibnum_scaled_norm2(n,dx,xscale);
           hkpri = normdx/normdxkm1*hk; 
           lambda = MIN(1.0,1.0/((1.0+rho)*hkpri)); 
         };
       nordlin += liniter;
       data->noiterord = liniter;
       data->tolord    = optlin->tol;
       data->precord   = infolin->precision;
       hk = hkpri;
       if ( MONITORLEVEL > 1 ) 
         giantgb_monitor(k,n,normdx,normfk,lambda_prev,optlin->tol,
                         infolin->precision,liniter,' ');
       if ( normdx <= xtol ) /* solution found, if condition is true */
         { data->normf     = normfk;     data->normdx = normdx;
           data->normdxbar = normdxbar;  data->lambda = lambda;
           giant_dataout(k,n,x,data);
           for (i=0;i<n;i++) x[i] += dx[i];  k++; 
           info->precision = normdx;
           RCODE=0; break; 
         };  
       reducted = False;  data->nocorr = -1;
       checkregularity:   data->nocorr++;  data->monviolated = False;
       if ( lambda < lambda_min ) 
         { RCODE=3; break; };  /* stop, convergence failure! */
       for (i=0;i<n;i++) xkp1[i]=x[i]+lambda*dx[i]; /* new trial iterate */
       f(n,xkp1,fxk,nfcn,&fail);  nfcn++;
       if ( fail<0 || fail>2 ) { RCODE=82; goto errorexit; }
       else if ( fail==1 || fail== 2 )
         { if ( fail==1 ) reduction_factor = 0.5;
           else           reduction_factor = fxk[0];
           if ( reduction_factor <= 0.0 || reduction_factor >= 1.0 )
             { RCODE=82; goto errorexit; };
           if ( MONITORLEVEL>1 )
             fprintf(MONITOR," %4i     FUN could not be evaluated  %7f\n",
                             k,lambda);
           if ( lambda > lambda_min )
             lambda = MAX(lambda*reduction_factor,lambda_min);
           else 
             lambda = lambda*reduction_factor;
           reducted = True;
           goto checkregularity;  
         };
       normfkp1 = zibnum_norm2(n,fxk);
       s = 1.0-lambda;
       for (j=0;j<n;j++) { w[j]=-fxk[j];  dxbar[j]=s*dx[j]; };
       optlin->tol = rhobar_max;
       optlin->giant_simple = True;
       gbit(n,dxbar,matvec,preconl,w,optlin,infolin);
       fail = infolin->rcode;
       nsimlin += infolin->iter;
       data->noitersim = infolin->iter;
       nmuljac += infolin->nomatvec;  nprecon += infolin->noprecon;
       if ( fail != 0 ) { RCODE=81; goto errorexit; };
       normdxbar = zibnum_scaled_norm2(n,dxbar,xscale);
       if ( MONITORLEVEL > 1 ) 
         giantgb_monitor(k,n,normdxbar,normfkp1,lambda,optlin->tol,
                         infolin->precision,infolin->iter,'*');
       theta = normdxbar/normdx;
       s = 1.0-lambda;
       for (i=0;i<n;i++) w[i] = dxbar[i]-s*dx[i];
       normdiff = zibnum_scaled_norm2(n,w,xscale);
       rhobar = (infolin->precision)*safetyfact*normdxbar/normdiff;
       hk = 2.0*(1.0-rhobar)*normdiff/(lambda*lambda*normdx);
       if ( ( !restricted && theta >= 1.0 ) || 
            (  restricted && theta > 1.0-lambda/4.0) )
         { lambda_new = MIN(1.0/(hk*(1.0+rho)),0.5*lambda);
           if ( lambda <= lambda_min ) lambda = lambda_new;
           else                        lambda = MAX(lambda_new,lambda_min);
           reducted = True;
           data->monviolated = True;
           goto checkregularity; };
       lambda_new = MIN(1.0,1.0/(hk*(1.0+rho)));
       if ( lambda==1.0 && lambda_new==1.0 )
         {
           if ( normdxbar <= xtol )
             { for(i=0;i<n;i++) x[i] += dxbar[i];  
               info->precision = normdxbar; RCODE=0;  break; };
         }
       else
         { if( lambda_new >= 4.0*lambda && !reducted ) 
             { lambda=lambda_new; goto checkregularity; };
         };
       data->normf     = normfk;     data->normdx = normdx;
       data->normdxbar = normdxbar;  data->lambda = lambda;
       giant_dataout(k,n,x,data);
       data->mode = Intermediate;
       /* save previous iterate for scaling purposes and accept new iterate */
       for (i=0;i<n;i++) { w[i]=x[i]; x[i]=xkp1[i]; }
       /* next step */
       k++;
       normfk = normfkp1;
     }
   while ( k <= max_iter && RCODE == 2 );

   data->normf     = normfkp1;   data->normdx = normdx;
   data->normdxbar = normdxbar;  data->lambda = lambda;
   data->mode = ( RCODE==0 ? Solution : Final );
   giant_dataout(k,n,x,data);
     
   errorexit:
   if ( ERRORLEVEL > 0 && RCODE != 0 )
     {
       switch ( RCODE )
        {
         case     2:
           fprintf(ERROR,"\n Error - Maximum allowed number of iterations exceeded\n");
           break;
         case     3:
           fprintf(ERROR,"\n Error - no convergence, damping factor became too small\n");
           break;
         case    81:
           fprintf(ERROR,"\n Error return from gbit: fail=%i\n",fail);
           break;
         case    82:
           fprintf(ERROR,"\n Error return from problem function: fail=%i\n",
                         fail);
           break;
         case    83:
           fprintf(ERROR,"\n Error return from Jacobian function: fail=%i\n",
                         fail);
           break;
         default   :
           fprintf(ERROR,"\n Error, code=%i,  subcode=%i\n",RCODE,fail);
        };
     };
   info->subcode = fail;
   if (io_allocated) {free(giant_ioctl); giant_ioctl=NULL;};
   free(data);
   free(dx); free(xkp1);
   free(w);
   if (scale_allocated)  { free(xscale);  opt->scale = NULL; };
   info->normdx    = normdx;
   info->iter       = k; 
   info->noordlinit = nordlin;
   info->nosimlinit = nsimlin;
   info->nofunevals = nfcn;
   info->nojacevals = njac;
   info->nomuljac   = nmuljac;
   info->noprecon   = nprecon;
   info->fx         = fxk;
}

int giantgb_initscale(int n, double **scale, struct GIANT_OPT *opt)
{  int i, rcode;
   double tol = opt->tol;
   rcode = zibnum_fwalloc(n,scale,"scale");  if ( rcode !=0 ) return rcode;
   if ( opt->nonlin <= Mildly_Nonlinear )
     for (i=0;i<n;i++) (*scale)[i]=1.0;
   else
     for (i=0;i<n;i++) (*scale)[i]=tol;
   return 0;
}

void giantgb_rescale(int n, double *scale, double *x, double *xa,
                     double *xthresh)
{  int i;
   for (i=0;i<n;i++) 
     scale[i] = MAX(xthresh[i],MAX(0.5*(fabs(x[i])+fabs(xa[i])),SMALL));
}
void giantgb_monitor(int k, int n, double normdx, double normf, double lambda,
                      double tol, double prec, int itlin, char indicator)
{  fprintf(MONITOR," %4i  %1c  %12e  %12e  %10.4e  %10.4e  %4i  %7f\n",
                   k,indicator,normdx,normf,tol,prec,itlin,lambda);
}
