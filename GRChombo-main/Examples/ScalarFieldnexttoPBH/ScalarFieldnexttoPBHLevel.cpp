/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "ScalarFieldnexttoPBHLevel.hpp"
#include "BoxLoops.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "SixthOrderDerivatives.hpp"
#include "TraceARemoval.hpp"

// For GW extraction
#include "MatterWeyl4.hpp"
#include "Weyl4.hpp"
#include "WeylExtraction.hpp"

// For RHS update
#include "MatterCCZ4RHS.hpp"

// For constraints calculation
#include "NewMatterConstraints.hpp"

// For tag cells
#include "MyChiAndPhiTaggingCriterion.hpp"
#include "ChiAndPhiTaggingCriterion.hpp"
#include "ChiTaggingCriterion.hpp"
#include "FixedGridsTaggingCriterion.hpp"

// Problem specific includes
#include "ComputePack.hpp"
#include "GammaCalculator.hpp"
#include "Potential.hpp"
#include "ScalarField.hpp"
#include "SetValue.hpp"

// Things to do at each advance step, after the RK4 is calculated
void ScalarFieldnexttoPBHLevel::specificAdvance()
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(),
                          PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(
            NanCheck(m_dx, m_p.center, "NaNCheck in specific Advance"),
            m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}


// Initial data for field and metric variables
void ScalarFieldnexttoPBHLevel::initialData()
{
    CH_TIME("ScalarFieldnexttoPBHLevel::initialData");
    if (m_verbosity)
        pout() << "ScalarFieldnexttoPBHLevel::initialData " << m_level << endl;

    // First set everything to zero then initial conditions for the system
    // here a Schwarzchild Black Hole and a Scalar Field profile
    // BoxLoops::loop(make_compute_pack(SetValue(0.)),
    //    m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // fillAllGhosts();
    // BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
    //               EXCLUDE_GHOST_CELLS);
}

#ifdef CH_USE_HDF5
// Things to do before outputting a checkpoint file
void ScalarFieldnexttoPBHLevel::prePlotLevel()
{
    #ifdef USE_AHFINDER
    // already calculated in 'specificPostTimeStep'
    if (m_bh_amr.m_ah_finder.need_diagnostics(m_dt, m_time))
        return;
    #endif

    fillAllGhosts();
    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);

    BoxLoops::loop(make_compute_pack(
        MatterWeyl4<ScalarFieldWithPotential>(
            scalar_field, m_p.extraction_params.center, m_dx,
    		m_p.formulation, m_p.G_Newton
            ),
        MatterConstraints<ScalarFieldWithPotential>(
            scalar_field, m_dx, m_p.G_Newton, c_Ham, Interval(c_Mom1, c_Mom3)
            )
        ),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    
}
#endif

// Things to do in RHS update, at each RK4 step
void ScalarFieldnexttoPBHLevel::specificEvalRHS(GRLevelData &a_soln,
                                                GRLevelData &a_rhs,
                                       const double a_time)
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(),
                          PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);

    if (m_p.max_spatial_derivative_order == 4)
    {
        MatterCCZ4RHS<ScalarFieldWithPotential, MovingPunctureGauge,
                      FourthOrderDerivatives>
            my_ccz4_matter(scalar_field, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
    else if (m_p.max_spatial_derivative_order == 6)
    {
        MatterCCZ4RHS<ScalarFieldWithPotential, MovingPunctureGauge,
                      SixthOrderDerivatives>
            my_ccz4_matter(scalar_field, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
}

// Things to do at ODE update, after soln + rhs
void ScalarFieldnexttoPBHLevel::specificUpdateODE(GRLevelData &a_soln,
                                         const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void ScalarFieldnexttoPBHLevel::preTagCells()
{
    // In the example of scalar field with no backreaction to the spacetime,
    // The spacetime is coupled with a scalar field and a black hole.
    // 
    // We use chi and phi in the tagging criterion
    // so fill the ghosts for chi and phi
    fillAllGhosts(VariableType::evolution, Interval(c_chi, c_chi));
    fillAllGhosts(VariableType::evolution, Interval(c_phi, c_phi));
}

void ScalarFieldnexttoPBHLevel::computeTaggingCriterion(
    FArrayBox &tagging_criterion, const FArrayBox &current_state,
    const FArrayBox &current_state_diagnostics)
{
    BoxLoops::loop(
        MyChiAndPhiTaggingCriterion(m_dx, m_level, m_p.L, m_p.center,
        	m_p.regrid_th_chi, m_p.regrid_ths_phi),
        	current_state, tagging_criterion);
    // BoxLoops::loop(
    //     FixedGridsTaggingCriterion(m_dx, m_level, 2.0 * m_p.L, m_p.center),
    //     current_state, tagging_criterion);
}

void ScalarFieldnexttoPBHLevel::specificPostTimeStep()
{
    CH_TIME("ScalarFieldnexttoPBHLevel::specificPostTimeStep");
    

    bool first_step =
        (m_time == 0.); // this form is used when 'specificPostTimeStep' was
                        // called during setup at t=0 from Main
    // bool first_step = (m_time == m_dt); // if not called in Main

    if (m_p.activate_extraction == 1)
    {
        int min_level = m_p.extraction_params.min_extraction_level();
        bool calculate_weyl = at_level_timestep_multiple(min_level);
        if (calculate_weyl)
        {
            // Populate the Weyl Scalar values on the grid
            fillAllGhosts();
            
            // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
				Potential potential(m_p.potential_params);
				ScalarFieldWithPotential scalar_field(potential);
           MatterWeyl4<ScalarFieldWithPotential> my_weyl4_matter(
    				scalar_field, m_p.extraction_params.center, m_dx,
    				m_p.formulation, m_p.G_Newton);
				BoxLoops::loop(my_weyl4_matter,
               m_state_new, m_state_diagnostics,
               EXCLUDE_GHOST_CELLS);
               
            // BoxLoops::loop(
            //    MatterWeyl4(m_p.extraction_params.center, m_dx, m_p.formulation),
            //    m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

            // Do the extraction on the min extraction level
            if (m_level == min_level)
            {
                CH_TIME("WeylExtraction");
                // Now refresh the interpolator and do the interpolation
                // fill ghosts manually to minimise communication
                bool fill_ghosts = false;
                m_bh_amr.m_interpolator->refresh(fill_ghosts);
                m_bh_amr.fill_multilevel_ghosts(
                    VariableType::diagnostic, Interval(c_Weyl4_Re, c_Weyl4_Im),
                    min_level);
                WeylExtraction my_extraction(m_p.extraction_params, m_dt,
                                             m_time, first_step,
                                             m_restart_time);
                my_extraction.execute_query(m_bh_amr.m_interpolator);
            }
        }
    }


#ifdef USE_AHFINDER
    // if print is on and there are Diagnostics to write, calculate them!
    if (m_bh_amr.m_ah_finder.need_diagnostics(m_dt, m_time))
    {
        fillAllGhosts();
        BoxLoops::loop(Constraints(m_dx, c_Ham, Interval(c_Mom1, c_Mom3)),
                       m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    }
    if (m_p.AH_activate && m_level == m_p.AH_params.level_to_run)
        m_bh_amr.m_ah_finder.solve(m_dt, m_time, m_restart_time);
#endif



}