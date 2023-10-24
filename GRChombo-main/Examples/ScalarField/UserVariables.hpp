/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "ArrayTools.hpp"
#include "CCZ4UserVariables.hpp"
#include "DiagnosticVariables.hpp"

/// This enum gives the index of every variable stored in the grid
enum
{
    // Note that it is important that the first enum value is set to 1 more than
    // the last CCZ4 var enum
    c_phi = NUM_CCZ4_VARS,  // matter field added
    c_Pi,                   // (minus) conjugate momentum

    NUM_VARS
};

namespace UserVariables
{
    static const std::array<std::string, NUM_VARS - NUM_CCZ4_VARS>
        user_variable_names = {"phi", "Pi"};

    static const std::array<std::string, NUM_VARS> variable_names =
        ArrayTools::concatenate(ccz4_variable_names, user_variable_names);
// static const std::array<std::string, NUM_VARS> variable_names =
//     ccz4_variable_names;
} // namespace UserVariables

#include "UserVariables.inc.hpp"

// uncomment to look for chi instead of expansion
// #define USE_CHI_CONTOURS

#ifdef USE_CHI_CONTOURS
#include "AHFunctions.hpp"
#define AHFunction ChiContourFunction // change default to chi contours
#endif

#endif /* USERVARIABLES_HPP */