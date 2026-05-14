mata:mata clear
mata:
// GMM evaluator function.
// Handles only d0-type optimization; todo, g and H are just ignored.
// betas is the parameter set over which we optimize, and
// crit is the objective function to minimize.
// m_acfcrit 1.0.0 jamc 25May2015

void m_acfcrit_het(todo,beta,j,g,H) 
{
	external y, X, X_lag, Z , e, omegahm, omegaht, robustflag, PHI, PHI_lag,Wh
	real matrix W, QZZ
	real vector gbar
	real scalar N
		
	cons = J(rows(PHI),1,1)	

	RHS = PHI - X * beta'
	OMEGA_lag = PHI_lag - X_lag*beta'
	OMEGA_lag2 = (OMEGA_lag:^2)
	OMEGA_lag3 = (OMEGA_lag:^3)
	OMEGA_lag_pol = (OMEGA_lag, OMEGA_lag2, OMEGA_lag3,cons)

	N = rows(Z)
	
	g_b = invsym(OMEGA_lag_pol'*OMEGA_lag_pol)*OMEGA_lag_pol'RHS
	g_b = g_b'
	e = RHS - OMEGA_lag_pol * g_b'
	
	gbar = 1/N * quadcross(Z, e)
	j = N * gbar' * Wh * gbar
}
end

mata: mata mosave m_acfcrit_het(), dir(PERSONAL) replace


