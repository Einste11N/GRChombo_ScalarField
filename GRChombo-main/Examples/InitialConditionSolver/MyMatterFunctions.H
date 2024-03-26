/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MYMATTERFUNCTIONS_HPP
#define MYMATTERFUNCTIONS_HPP

#include "SetPsiAndAij.H"

Real my_phi_function(const RealVect loc, const PoissonParameters &a_params,
                     const RealVect &a_dx)
{

    RealVect L = a_params.domainLength;

	Real rr = sqrt(D_TERM((loc[0] - L[0] / 16 * 5) * (loc[0] - L[0] / 16 * 5), +loc[1] * loc[1], +loc[2] * loc[2]));
    	// sqrt(D_TERM((loc[0] - L[0] / 16 * 5) * (loc[0] - L[0] / 16 * 5), +loc[1] * loc[1], +loc[2] * loc[2]));
		// sqrt(D_TERM(loc[0] * loc[0], +loc[1] * loc[1], +loc[2] * loc[2]));

    
    double amp = a_params.phi_amplitude;
    double back = a_params.phi_background;
	
    Real bubble = back + amp * 0.5 * a_params.phi_0 *  (1.0 - tanh((rr - a_params.r_0)/a_params.D_0));
	return  bubble;
}

Real my_Pi_function(const RealVect loc, const PoissonParameters &a_params,
                    const RealVect &a_dx)
{
	return 0;
}

Real my_potential_function(const Real &phi_here,
                           const PoissonParameters &a_params)
{
	Real VF = a_params.pot_VF;
	Real VT = a_params.pot_VT;
	Real a = a_params.pot_a;
	
    return VF * (1.0 + a * pow(phi_here/a_params.phi_0, 2)
			        -(2*a - 4) * pow(phi_here/a_params.phi_0, 3)  
			         + (a + 3) * pow(phi_here/a_params.phi_0, 4) )
            + VT;
}

Real my_rho_function(const IntVect a_iv, const RealVect &a_dx,
                     FArrayBox &a_multigrid_vars_box,
                     const PoissonParameters &a_params)
{

    RealVect loc;
    get_loc(loc, a_iv, a_dx, a_params);

    Real psi_reg = a_multigrid_vars_box(a_iv, c_psi_reg);
    Real psi_bh = set_binary_bh_psi(loc, a_params);
    Real psi_0 = psi_reg + psi_bh;

    Real V_of_phi =
        my_potential_function(a_multigrid_vars_box(a_iv, c_phi_0), a_params);
    Real Pi_0 = a_multigrid_vars_box(a_iv, c_Pi_0);
    Tensor<1, Real, SpaceDim> d1_phi =
        get_d1(a_iv, a_multigrid_vars_box, a_dx, c_phi_0);
    Real d1_phi_squared = 0;
    FOR1(i) { d1_phi_squared += d1_phi[i] * d1_phi[i]; }

    return 0.5 * Pi_0 * Pi_0 + V_of_phi +
           0.5 * d1_phi_squared * pow(psi_0, -4.0);
}

Real my_Si_function(Tensor<1, Real> &Si, const IntVect a_iv,
                    const RealVect &a_dx, FArrayBox &a_multigrid_vars_box,
                    const PoissonParameters &a_params)
{

    RealVect loc;
    get_loc(loc, a_iv, a_dx, a_params);

    Real Pi_0 = a_multigrid_vars_box(a_iv, c_Pi_0);
    Tensor<1, Real, SpaceDim> d1_phi =
        get_d1(a_iv, a_multigrid_vars_box, a_dx, c_phi_0);

    FOR1(i) { Si[i] = -Pi_0 * d1_phi[i]; }
}
#endif /* MYMATTERFUNCTIONS_HPP */
