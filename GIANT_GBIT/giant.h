/*
    Common declarations for programs of the NewtonLib package.
    
 *  Written by        L. Weimann 
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
 
*/
#define _GIANT
#include "itlin.h"

#define RCODE info->rcode
#define MIN(A,B)  ( A < B ? A : B )
#define MAX(A,B)  ( A > B ? A : B )
#define SIGN(A)   ( A > 0 ? 1 : -1 )

#define SMALL  1.0e-150
#define EPMACH 1.0e-17

typedef void GIANT_FFUN(int,double*,double*,int,int*);
typedef void GIANT_JFUN(int,double*,double*,double*,int,int*);
typedef enum {Mildly_Nonlinear=2, Highly_Nonlinear=3, Extremely_Nonlinear=4}
             PROBLEM_TYPE ;
typedef enum {StandardMode=0, QuadraticMode=1} ADAPT_MODE;

struct GIANT_FUN
{
   GIANT_FFUN *fun;
   GIANT_JFUN *jac;
   MATVEC     *muljac;
   PRECON     *preconr, *preconl;
};

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

struct GIANT_INFO
{
   double precision, normdx;
   double *fx, *dx;
   int iter, rcode, subcode, nofunevals, nojacevals, nomuljac, noprecon,
       noordlinit, nosimlinit;
};

struct GIANT_DATA
{
  double *fx, *dx;
  double normf, normdx, normdxbar, lambda, theta, tolord, precord;
  int noiterord, noitersim, maxiter, nocorr; 
  LOGICAL monviolated;
  enum { GIANT_GMRES=0, GIANT_GBIT=1 } codeid;
  DATA_MODE mode;
};

struct GIANT_IO
{
   FILE *errfile, *monfile, *datfile,
        *iterfile, *resfile, *miscfile;
   PRINT_LEVEL errlevel, monlevel, datlevel, plotlevel;
};

#define ERRORLEVEL   giant_ioctl->errlevel
#define MONITORLEVEL giant_ioctl->monlevel
#define DATALEVEL    giant_ioctl->datlevel
#define PLOTLEVEL    giant_ioctl->plotlevel
#define ERROR        giant_ioctl->errfile
#define MONITOR      giant_ioctl->monfile
#define DATA         giant_ioctl->datfile
#define FITER        giant_ioctl->iterfile
#define FRES         giant_ioctl->resfile
#define FMISC        giant_ioctl->miscfile

extern void daxpy_(int *n, double *alpha, double *x, int *incx,
                   double *y, int *incy);

/* routines defined in utils.c */
int    zibnum_fwalloc(int size, double **ptr, char vname[]);
int    zibnum_iwalloc(int size, int **ptr, char vname[]);
int    zibnum_pfwalloc(int size, double ***ptr, char vname[]);
double zibnum_scaled_norm2(int n, double *v, double *scale);
double zibnum_scaled_sprod(int n, double *v1, double *v2, double *scale);
double zibnum_norm2(int n, double *v);
void   giant_dataout(int k, int n, double *x, struct GIANT_DATA *data);
int    giant_parcheck_and_print(int n, struct GIANT_OPT *opt,
                               struct GIANT_FUN fun,
                               int giant_code);

/* main routines ( see file <routine-name>.c ) */
extern void giant_gmres(struct GIANT_FUN fun, int n, double *x,
                  struct GIANT_OPT *opt, struct GIANT_INFO *info);
