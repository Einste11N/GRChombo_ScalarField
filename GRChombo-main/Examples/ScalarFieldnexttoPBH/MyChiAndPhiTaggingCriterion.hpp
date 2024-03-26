/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MYCHIANDPHITAGGINGCRITERION_HPP_
#define MYCHIANDPHITAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "ScalarField.hpp"
#include "Tensor.hpp"
#include <vector>

class MyChiAndPhiTaggingCriterion
{
  protected:
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    const double m_L;
    const int m_level;
    const std::array<double, CH_SPACEDIM> m_center;
    const double m_threshold_chi;
    const Vector<double> m_threshold_phi;

    template <class data_t>
    using MatterVars = typename ScalarField<>::template Vars<data_t>;

    /// Vars object for chi
    template <class data_t> struct Vars
    {
        data_t chi; //!< Conformal factor

        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            using namespace VarsTools; // define_enum_mapping is part of
                                       // VarsTools
            define_enum_mapping(mapping_function, c_chi, chi);
        }
    };

  public:
    MyChiAndPhiTaggingCriterion(const double dx, const int a_level, const double a_L,
                              const std::array<double, CH_SPACEDIM> a_center, 
                              const double threshold_chi,
                              const Vector<double> threshold_phi)
        : m_dx(dx), m_deriv(dx), m_level(a_level), 
          m_L(a_L), m_center(a_center),
          m_threshold_chi(threshold_chi),
          m_threshold_phi(threshold_phi){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto d2 = m_deriv.template diff2<MatterVars>(current_cell);
        const auto d2chi = m_deriv.template diff2<Vars>(current_cell);
        
        data_t mod_d2_chi = 0;
        data_t mod_d2_phi = 0;

        FOR(idir, jdir)
        {
            mod_d2_chi += d2chi.chi[idir][jdir] * d2chi.chi[idir][jdir];
            mod_d2_phi += d2.Pi[idir][jdir] * d2.Pi[idir][jdir] +
                          d2.phi[idir][jdir] * d2.phi[idir][jdir];
        }

        // const Coordinates<data_t> coords(current_cell, m_dx, m_center);
        // data_t rr = coords.get_radius();
        // auto regrid_bool = simd_compare_lt(rr, m_L / 8.0);
        // data_t criterion_chi = simd_conditional(regrid_bool, m_dx / m_threshold_chi * sqrt(mod_d2_chi), 0.);
        
        data_t criterion_chi = m_dx / m_threshold_chi * sqrt(mod_d2_chi);
        data_t criterion_phi = m_dx / m_threshold_phi[m_level] * sqrt(mod_d2_phi);        
        
        data_t criterion = simd_max(criterion_chi, criterion_phi);

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* MYCHIANDPHITAGGINGCRITERION_HPP_ */
