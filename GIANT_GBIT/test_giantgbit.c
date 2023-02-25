/*
 
      Testexample for GIANT_GBIT: Artificial test problem.
 
 *  Written by        L. Weimann 
 *  Purpose           Testexample for code GIANT_GBIT
 *  Version           1.0
 *  Revision          June 2006
 *  Latest Change     June 2006
 *  Library           NewtonLib
 *  Code              C, Double Precision
 *  Environment       Standard C environment on PC's,
                      workstations and hosts.
 
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "giant.h"

void fun(int n, double *u, double *fu, int nfcn, int *fail);
void jac(int n, double *u, double *xscale, double *fu, int njac, int *fail);

void giant_gbit(struct GIANT_FUN fun, int n, double *x,
                  struct GIANT_OPT *opt, struct GIANT_INFO *info);

#ifdef  __cplusplus
extern "C" {
#endif

void _stdcall SLV6EQS(double *X,double *F);
void _stdcall X0VCT(double *X);
void _stdcall FDJAC6(double *X,double *FVEC,double *DF);

int __stdcall CMAIN()
{
  int  n=1,j;
  double *u;
  struct GIANT_FUN *problem = malloc( sizeof(struct GIANT_FUN) );
  struct GIANT_OPT  *opt  = malloc( sizeof(struct GIANT_OPT) );
  struct GIANT_INFO *info = malloc( sizeof(struct GIANT_INFO) );

  u=malloc(n*sizeof(double));
  for (j=0;j<n;j++) {u[j]=0.0;};
//	X0VCT(u);
	u[0]=1.0;
	u[1]=1.0;

  problem->fun    = &fun;     
  problem->jac    =NULL; //&jac;  
  problem->muljac =NULL; //&muljac;  
  problem->preconr =NULL; //&precon;
  problem->preconl =NULL; //&precon;
  /* ----------------------------- */
  opt->rho = 0.2;
  opt->tol = 1.0e-12;
  opt->maxiter = 50;
  opt->lin_maxiter = 100;
  opt->i_max =10;
  opt->errorlevel = Verbose;
  opt->monitorlevel = Verbose;
  opt->datalevel = Verbose;
  opt->errorfile = NULL;
  opt->monitorfile = NULL;
  opt->datafile = NULL;
  opt->iterfile = NULL;
  opt->resfile  = NULL;
  opt->linmonfile   = NULL;
  opt->lindatafile  = NULL;
  opt->miscfile = NULL;
  opt->linmonlevel  = Verbose;
  opt->lindatalevel = None;
  opt->nonlin = Mildly_Nonlinear;
  opt->restricted = False;
  opt->scaleopt = StandardScale;
  opt->scale = NULL;
  giant_gbit(*problem, n, u, opt, info);
  if ( info->rcode == 0 )
    fprintf(stdout,"\n precision=%e\n",info->precision);
  else
     fprintf(stdout,"\n GIANT_GBIT failed - no solution found\n");
  fprintf(stdout," iter    = %i\n",info->iter);
  fprintf(stdout," rcode   = %i\n",info->rcode);
  fprintf(stdout," subcode = %i\n",info->subcode);
  fprintf(stdout," nfun    = %i\n",info->nofunevals);
  fprintf(stdout," njac    = %i\n",info->nojacevals);
  fprintf(stdout," nmuljac = %i\n",info->nomuljac);
  fprintf(stdout," nprecon = %i\n",info->noprecon);
  return 0;
}

#ifdef __cplusplus
}
#endif


void fun(int n, double *u, double *f, int nfcn, int *fail)
{
//	SLV6EQS(u,f);
	f[0]=u[0]-2*u[1]+2;
	f[1]=u[0]+u[1]-1;
	*fail=0;
}

void jac(int n, double *u, double *xscale, double *fu, int njac, int *fail)
{
//  double *fvec;
//  fvec = malloc( n*sizeof(double) );
//	SLV6EQS(u,fvec);
//	FDJAC6(u,fvec,fu);
	fu[0]=1;
	fu[1]=1;
	return;
}

