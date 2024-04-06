/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "Potential.hpp"

// General includes
#include "ArrayTools.hpp"
#include "FilesystemTools.hpp"
#include "VariableType.hpp"
#include "unistd.h" // gives 'access'
#include <algorithm>
#include <string>


class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // read the problem specific params
        read_params(pp);
        check_params();
    }

    void read_params(GRParmParse &pp)
    {
        pp.load("G_Newton", G_Newton, 1.0);

        // Potential data
        pp.load("false_vacuum_potential",
                potential_params.false_vacuum_potential, 0.00001);
        pp.load("true_vacuum_potential", potential_params.true_vacuum_potential,
                0.0);
        pp.load("barrier_amplitude",
                potential_params.barrier_amplitude_parameter, 5.0);
        pp.load("true_vacuum_expectation_value", potential_params.phi_VEV, 0.01);
        pp.load("eta_friction", potential_params.eta, 0.00005);

        // The following are not potential parameters
        pp.load("BH_mass", bh_mass, 0.5);
        pp.load("Grid_Center", my_grid_center);
        pp.load("regrid_th_chi", regrid_th_chi, 0.5);
        // pp.load("regrid_th_phi", regrid_th_phi, 0.1);
        pp.load("max_level", Max_Level);
        
        if (pp.contains("regrid_ths_phi"))
        {
            pp.getarr("regrid_ths_phi", regrid_ths_phi, 0, Max_Level);
            regrid_ths_phi.resize(Max_Level + 1);
            regrid_ths_phi[Max_Level] = regrid_ths_phi[Max_Level - 1];
        }
        else
        {
            double regrid_th_phi;
            pp.load("regrid_th_phi", regrid_th_phi, 0.5);
            regrid_ths_phi = Vector<double>(max_level + 1, regrid_th_phi);
        }
        

#ifdef USE_AHFINDER
        double AH_guess = 0.5 * bh_mass;
        pp.load("AH_initial_guess", AH_initial_guess, AH_guess);
#endif
    }

    void check_params()
    {
        // By now no parameters need to be checked
    }

    // Initial data for matter and potential and BH
    double G_Newton;
    double bh_mass;
    double regrid_th_chi;
    // double regrid_th_phi;
    int Max_Level;
    Vector<double> regrid_ths_phi;

    Potential::params_t potential_params;
    std::array<double, CH_SPACEDIM> my_grid_center;

#ifdef USE_AHFINDER
    double AH_initial_guess;
#endif
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
