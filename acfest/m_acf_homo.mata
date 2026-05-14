mata:mata clear
version 10.1

mata:
// m_acf_homo 1.0.0 MES/CFB 11aug2008
void m_acf_homo(string scalar yname,
	string scalar exname,
	string scalar exnamelag,
	string scalar phi,
	string scalar phi_lag,
	string scalar instname,
	string scalar touse,
	string scalar robust)
{
	real matrix W, X1, Z1, X2, POL, XPOL,QZZ, QZX, omegaht
	real vector cons, p, bols, beta_ols, BOLS, e2
	real scalar K, L, N, S, j, k, q, s

	external y, X, X_lag, Z, e, omega, Wh, robustflag, PHI, PHI_lag, beta_hom

// Robustflag is the robust argument recreated as an external
// MATA scalar
	robustflag=robust

// Views of the data
	st_view(y, ., tokens(yname), touse)
	st_view(X1, ., tokens(exname), touse)
	st_view(X2, ., tokens(exnamelag), touse)
	st_view(Z1, ., tokens(instname), touse)
	st_view(PHI, ., tokens(phi), touse)
	st_view(PHI_lag, ., tokens(phi_lag), touse)

// Related matrices
	
	cons = J(rows(X1),1,1)
	X = (X1, cons)
	X_lag = (X2, cons)
	Z = (Z1, cons)
	K = cols(X)
	L = cols(Z)
	N = rows(y)
	QZZ = 1/N * quadcross(Z, Z)
	QZX = 1/N * quadcross(Z, X)
	
// To pass MATA the initial values from the OLS estimation
	bols = st_matrix("bols")
	q = cols(X1)
	s = cols(bols)
	beta_ols = (bols[.,(1..q)], bols[.,(s)])

	S=optimize_init()
	optimize_init_evaluator(S, &m_acfcrit_homo())
	optimize_init_evaluatortype(S, "d0")
	optimize_init_technique(S, "nm")
	optimize_init_nmsimplexdeltas(S, 0.1)
	optimize_init_which(S, "min")
	optimize_init_params(S,beta_ols)
	p = optimize(S)
	j = optimize_result_value(S)

		
// Weights matrix for the heterokedastic non-linear GMM estimator
	e2 = e:^2
    omegaht = 1/N * quadcross(Z, e2, Z)
    _makesymmetric(omegaht)
	Wh =invsym(omegaht)
		
// Starting values		
	beta_hom = p		
		
// Easiest way to return results to Stata:
// as r-class macros
	
	st_matrix("r(beta)", p)
	st_numscalar("r(j)", j)
	st_numscalar("r(N)", N)
	st_numscalar("r(L)", L)
	st_numscalar("r(K)", K)	
	
}
end

mata: mata mosave m_acf_homo(), dir(PERSONAL) replace

