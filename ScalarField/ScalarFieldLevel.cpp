/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */


#include "ScalarFieldLevel.hpp"
#include "BoxLoops.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "SixthOrderDerivatives.hpp"
#include "TraceARemoval.hpp"

// Initial data
#include "ComputePack.hpp"
#include "GammaCalculator.hpp"
#include "KerrBH.hpp"

#include "InitialScalarData.hpp"
#include "ScalarField.hpp"

#include "Potential.hpp"
#include "SetValue.hpp"


// For RHS update
#include "MatterCCZ4RHS.hpp"  // #include "CCZ4RHS.hpp"

// For constraints calculation
#include "NewMatterConstraints.hpp"  // #include "NewConstraints.hpp"

// For tag cells
#include "FixedGridsTaggingCriterion.hpp"  // #include "ChiTaggingCriterion.hpp"

void ScalarFieldLevel::specificAdvance()
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(
            NanCheck(m_dx, m_p.center, "NaNCheck in specific Advance"),
            m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}

void ScalarFieldLevel::initialData()
{
    CH_TIME("ScalarFieldLevel::initialData");
    if (m_verbosity)
        pout() << "ScalarFieldLevel::initialData " << m_level << endl;

    // First set everything to zero then calculate initial data  Get the Kerr
    // solution in the variables, then calculate the \tilde\Gamma^i numerically
    // as these are non zero and not calculated in the Kerr ICs
    BoxLoops::loop(
        make_compute_pack(SetValue(0.), KerrBH(m_p.kerr_params, m_dx), InitialScalarData(m_p.initial_params, m_dx)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    fillAllGhosts();
    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);

#ifdef USE_AHFINDER
    // Diagnostics needed for AHFinder
    BoxLoops::loop(Constraints(m_dx, c_Ham, Interval(c_Mom1, c_Mom3)),
                   m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
#endif
}

#ifdef CH_USE_HDF5
void ScalarFieldLevel::prePlotLevel()
{
#ifdef USE_AHFINDER
    // already calculated in 'specificPostTimeStep'
    if (m_bh_amr.m_ah_finder.need_diagnostics(m_dt, m_time))
        return;
#endif

    fillAllGhosts();
    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    BoxLoops::loop(
        MatterConstraints<ScalarFieldWithPotential>(
            scalar_field, m_dx, m_p.G_Newton, c_Ham, Interval(c_Mom1, c_Mom3)),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}
#endif /* CH_USE_HDF5 */

void ScalarFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                  const double a_time)
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(), 
                          PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate CCZ4 right hand side
    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    if (m_p.max_spatial_derivative_order == 4)
    {
        MatterCCZ4RHS<ScalarFieldWithPotential, MovingPunctureGauge,
                      FourthOrderDerivatives> 
            my_ccz4_matter(scalar_field, m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation, m_p.G_Newton);

        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
        // BoxLoops::loop(CCZ4RHS<MovingPunctureGauge, FourthOrderDerivatives>(
        //                    m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation),
        //                a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
    else if (m_p.max_spatial_derivative_order == 6)
    {
        MatterCCZ4RHS<ScalarFieldWithPotential, MovingPunctureGauge,
                      SixthOrderDerivatives>
            my_ccz4_matter(scalar_field, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
        // BoxLoops::loop(CCZ4RHS<MovingPunctureGauge, SixthOrderDerivatives>(
        //                    m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation),
        //                a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
}

void ScalarFieldLevel::specificUpdateODE(GRLevelData &a_soln,
                                    const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void ScalarFieldLevel::preTagCells()
{
    // We only use chi in the tagging criterion so only fill the ghosts for chi
    fillAllGhosts(VariableType::evolution, Interval(c_chi, c_chi));
}

void ScalarFieldLevel::computeTaggingCriterion(
    FArrayBox &tagging_criterion, const FArrayBox &current_state,
    const FArrayBox &current_state_diagnostics)
{
    // BoxLoops::loop(ChiTaggingCriterion(m_dx), current_state, tagging_criterion);
    BoxLoops::loop(
        FixedGridsTaggingCriterion(m_dx, m_level, 2.0 * m_p.L, m_p.center),
        current_state, tagging_criterion);
}

void ScalarFieldLevel::specificPostTimeStep()
{
    CH_TIME("ScalarFieldLevel::specificPostTimeStep");
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
