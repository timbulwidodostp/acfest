mata:mata clear
version 10.1

mata:
// m_acf_het 1.0.0 MES/CFB 11aug2008
void m_acf_het(string scalar yname,
	string scalar exname,
	string scalar exnamelag,
	string scalar phi,
	string scalar phi_lag,
	string scalar instname,
	string scalar touse,
	string scalar robust)
{
	real matrix W, X1, Z1, X2, POL, XPOL, QZZ, QZX
	real vector cons, p, bols, beta_ols
	real scalar K, L, N, S, j, k
	
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
	
// beta_hom passes MATA the initial values from the homokedastic estimation
// Wh are the weights matrix from the homocedastic estimation (see m_acfcrit_het)
	
	S=optimize_init()
	optimize_init_evaluator(S, &m_acfcrit_het())
	optimize_init_evaluatortype(S, "d0")
	optimize_init_technique(S, "nm")
	optimize_init_nmsimplexdeltas(S, 0.1)
	optimize_init_which(S, "min")
	optimize_init_params(S,(beta_hom))
	p = optimize(S)
	j = optimize_result_value(S)
	
// Easiest way to return results to Stata:
// as r-class macros
	
	st_matrix("r(beta)", p)
	st_numscalar("r(j)", j)
	st_numscalar("r(N)", N)
	st_numscalar("r(L)", L)
	st_numscalar("r(K)", K)
}
end

mata: mata mosave m_acf_het(), dir(PERSONAL) replace
