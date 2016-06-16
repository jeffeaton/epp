#include <R.h>
#include <Rinternals.h>

#include <math.h>

#define POP_15TO49_COL  0
#define POP_15ENTER_COL 1
#define POP_50EXIT_COL  2
#define POP_NETMIGR_COL 3
#define POP_A50RATE_COL 4
#define POP_MX_COL      5

#define EPP_RSPLINE 0
#define EPP_RTREND 1

const size_t DS = 8;  // number of disease stages (CD4 stages, including uninfected)
const size_t TS = 4;  // number


SEXP eppC(SEXP s_eppPopTS, SEXP s_projsteps, SEXP s_dt,
	   SEXP s_eppmod,
	   SEXP s_rspline_rvec, SEXP s_iota, SEXP s_relinfectART, SEXP s_tsEpidemicStart,
	   SEXP s_rtrend_beta, SEXP s_rtrend_tstabilize, SEXP s_rtrend_r0,
	   SEXP s_cd4init, SEXP s_cd4prog, SEXP s_cd4artmort,
	   SEXP s_artnumTS, SEXP s_arteligidxTS){

  size_t nsteps = length(s_projsteps);
  double dt = *REAL(s_dt);

  double *projsteps = REAL(s_projsteps);
  double *pop15to49_ts = REAL(VECTOR_ELT(s_eppPopTS, POP_15TO49_COL));
  double *age15enter_ts = REAL(VECTOR_ELT(s_eppPopTS, POP_15ENTER_COL));
  double *netmigr_ts = REAL(VECTOR_ELT(s_eppPopTS, POP_NETMIGR_COL));
  double *age50rate_ts = REAL(VECTOR_ELT(s_eppPopTS, POP_A50RATE_COL));
  double *mx_ts = REAL(VECTOR_ELT(s_eppPopTS, POP_MX_COL));

  int eppmod = *INTEGER(s_eppmod);

  double iota = *REAL(s_iota);
  double relinfectART = *REAL(s_relinfectART);
  double tsEpidemicStart = *REAL(s_tsEpidemicStart);

  double *rspline_rvec;
  double *rtrend_beta, rtrend_tstab, rtrend_r0;
  if(eppmod == EPP_RSPLINE)
    rspline_rvec = REAL(s_rspline_rvec);
  else {
    rtrend_beta = REAL(s_rtrend_beta);
    rtrend_tstab = *REAL(s_rtrend_tstabilize);
    rtrend_r0 = *REAL(s_rtrend_r0);
  }

  double *cd4init = REAL(s_cd4init);
  double *cd4prog = REAL(s_cd4prog);
  double cd4artmort[DS-1][TS];
  for(size_t m = 0; m < DS-1; m++)
    for(size_t u = 0; u < TS; u++)
      cd4artmort[m][u] = REAL(s_cd4artmort)[m + (DS-1)*u];

  // array cd4artmort_alt = {2, INTEGER(getAttrib(s_cd4artmort, R_DimSymbol)), REAL(s_cd4artmort)};

  double *artnum_ts = REAL(s_artnumTS);
  int *arteligidx_ts = INTEGER(s_arteligidxTS);
  
  
  // initialise output
  size_t numOutDates = ceil(dt * nsteps), outIdx = 0;
  SEXP s_Xout, s_Xout_dim, s_rvec; 
  PROTECT(s_Xout = allocVector(REALSXP, numOutDates * DS * TS));
  PROTECT(s_Xout_dim = allocVector(INTSXP, 3));
  INTEGER(s_Xout_dim)[0] = numOutDates;
  INTEGER(s_Xout_dim)[1] = DS;
  INTEGER(s_Xout_dim)[2] = TS;
  setAttrib(s_Xout, R_DimSymbol, s_Xout_dim);

  PROTECT(s_rvec = allocVector(REALSXP, nsteps));
  setAttrib(s_Xout, install("rvec"), s_rvec);

  double *Xout = REAL(s_Xout);
  double *rvec = REAL(s_rvec);


  // initialise population
  double X[DS][TS];
  for(size_t m = 0; m < DS; m++)
    for(size_t u = 0; u < TS; u++)
      X[m][u] = 0.0;
  X[0][0] = pop15to49_ts[0];

  double prevlast = 0.0; // store last prevalence value for r-trend

  // do timesteps
  for(size_t ts = 0; ts < nsteps; ts++){

    // record outputs
    if((projsteps[ts] - floor(projsteps[ts])) == dt*floor(1.0/(dt*2))){
        for(size_t m = 0; m < DS; m++)
    	  for(size_t u = 0; u < TS; u++)
    	    Xout[outIdx + (m + DS*u)*numOutDates] = X[m][u];
    	outIdx++;
    }

    // sum population sizes
    double Xhivp_noart = 0.0, Xonart = 0.0;
    for(size_t m = 1; m < DS; m++){
      Xhivp_noart += X[m][0];
      for(size_t u = 1; u < TS; u++)
	Xonart += X[m][u];
    }
    double Xtot = X[0][0] + Xhivp_noart + Xonart;
    double prevcurr = 1.0 - X[0][0] / Xtot;

    // calculate r(t)
    if(eppmod == EPP_RSPLINE)
      rvec[ts] = rspline_rvec[ts];
    else {
      double t = projsteps[ts];
      if(t > tsEpidemicStart){
	double gamma_t = (t < rtrend_tstab)?0.0:(prevcurr-prevlast) * (t - rtrend_tstab) / (dt * prevlast);
	double logr_diff = rtrend_beta[1]*(rtrend_beta[0] - rvec[ts-1]) + rtrend_beta[2]*prevlast + rtrend_beta[3]*gamma_t;
	rvec[ts] = exp(log(rvec[ts-1]) + logr_diff);
      } else {
	rvec[ts] = rtrend_r0;
      }
    }
    
    // ageing, natural mortality, and net migration
    double grad[DS][TS];
    for(size_t m = 0; m < DS; m++)
      for(size_t u = 0; u < TS; u++)
	grad[m][u] = X[m][u] * (-age50rate_ts[ts] - mx_ts[ts] + netmigr_ts[ts]/(Xtot * dt));

    // new entrants
    grad[0][0] += age15enter_ts[ts] / dt;

    // new HIV infections
    double incrate = rvec[ts] * (Xhivp_noart + relinfectART*Xonart)/Xtot + ((projsteps[ts] == tsEpidemicStart)?iota:0);
    grad[0][0] -= X[0][0] * incrate;
    for(size_t m = 1; m < DS; m++)
      grad[m][0] += X[0][0] * incrate * cd4init[m-1];

    // disease progression and mortality
    for(size_t m = 1; m < (DS-1); m++){  // cd4 stage progression (not on ART)
      grad[m][0] -= cd4prog[m-1] * X[m][0];
      grad[m+1][0] += cd4prog[m-1] * X[m][0];
    }
    for(size_t m = 1; m < DS; m++){
      grad[m][1] -= 2.0 * X[m][1];   // ART stage progression (HARD CODED 6 months mean duration)
      grad[m][2] += 2.0 * X[m][1];
      grad[m][2] -= 2.0 * X[m][2];   
      grad[m][3] += 2.0 * X[m][2];

      for(size_t u = 0; u < TS; u++)
	grad[m][u] -= cd4artmort[m-1][u] * X[m][u];  // HIV mortality
    }
    

    // ART initiation
    if(artnum_ts[ts] > 0){

      // determine number of desired ART initiations
      double artchange = 0.0;
      for(size_t m = 1; m < DS; m++)
	for(size_t u = 1; u < TS; u++)
	  artchange += grad[m][u];

      size_t artelig_idx = arteligidx_ts[ts] - 1; 	// -1 for 0-based indexing vs. 1-based indexing in R
      double art_anninits = (artnum_ts[ts] - Xonart) / dt - artchange;

      // determine number weights for number from each stage (average of expected mortality and eligibility)
      double sum_mortweight = 0.0, artelig = 0.0;
      for(size_t m = artelig_idx; m < DS; m++){
	artelig += X[m][0];
	sum_mortweight += cd4artmort[m-1][0] * X[m][0];
      }
      // determine rate to initiate each stage
      double artinitweight[DS], artstageinit[DS];
      for(size_t m = artelig_idx; m < DS; m++){
	artinitweight[m] = (cd4artmort[m-1][0]/sum_mortweight + 1.0/artelig)/2.0;
	artstageinit[m] = art_anninits * artinitweight[m] * X[m][0];
	artstageinit[m] = (artstageinit[m] > X[m][0]/dt) ? X[m][0]/dt : artstageinit[m]; // check not greater than number in stage
	artstageinit[m] = (artstageinit[m] > (X[m][0]/dt + grad[m][0])) ? (X[m][0]/dt + grad[m][0]) : artstageinit[m]; // check number on stage won't be less than 0;

	grad[m][0] -= artstageinit[m];
	grad[m][1] += artstageinit[m];
      }
    }

    // store prevalence to use in next time step
    prevlast = prevcurr;
      
    // do projection (euler integration)
    for(size_t m = 0; m < DS; m++)
      for(size_t u = 0; u < TS; u++)
	X[m][u] += dt*grad[m][u];
    
  }

  UNPROTECT(3);

  return s_Xout;
}
