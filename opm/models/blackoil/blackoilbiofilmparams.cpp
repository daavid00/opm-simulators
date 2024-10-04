// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#include <config.h>
#include <opm/models/blackoil/blackoilbiofilmparams.hpp>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PermporoTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PefactTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SimpleTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/Biofpara.hpp>
#endif

#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <type_traits>

namespace Opm {

#if HAVE_ECL_INPUT
template<class Scalar>
template<bool enableBiofilm>
void BlackOilBiofilmParams<Scalar>::
initFromState(const EclipseState& eclState)
{
    // Check if module is enabled but BIOFILM is not present
    if constexpr (enableBiofilm) {
        if (!eclState.runspec().biof()) {
            throw std::runtime_error("Biofilm module enabled at compile time, but the deck does not contain "
                                     "BIOFILM!");
        }
    }
    // Check for opposite of the above: module disabled but BIOFILM is in deck
    else {
        if (eclState.runspec().biof()) {
            throw std::runtime_error("Biofilm module disabled at compile time, but deck contains BIOFILM!");
        }
    }

    if (!eclState.runspec().biof()) {
        return; // biofilm is supposed to be disabled
    }

    const auto& tableManager = eclState.getTableManager();
    unsigned numSatRegions = tableManager.getTabdims().getNumSatTables();
    const TableContainer& permporoTables = tableManager.getPermporoTables();
    permporoTable_.resize(numSatRegions);
    for (size_t i = 0; i < permporoTables.size(); ++i) {
        const PermporoTable& permporoTable = permporoTables.getTable<PermporoTable>(i);
        permporoTable_[i].setXYContainers(permporoTable.getPorosityChangeColumn(), permporoTable.getPermeabilityMultiplierColumn());
    }
    const TableContainer& pefactTables = tableManager.getPefactTables();
    if (!pefactTables.empty()) {
        pefactTable_.resize(numSatRegions);
        for (size_t i = 0; i < pefactTables.size(); ++i) {
            const PefactTable& pefactTable = pefactTables.getTable<PefactTable>(i);
            pefactTable_[i].setXYContainers(pefactTable.getPorosityChangeColumn(), pefactTable.getPcMultiplierColumn());
        }
    }
    const auto& biofpara = tableManager.getBiofpara();
    if (!biofpara.empty()) {
        unsigned numSatRegions = tableManager.getTabdims().getNumSatTables();
        biofilmDensity_.resize(numSatRegions);
        maxGrowthRate_.resize(numSatRegions);
        halfVelocityCoeff_.resize(numSatRegions);
        yieldCoeff_.resize(numSatRegions);
        decayCoeff_.resize(numSatRegions);
        for (size_t i = 0; i < biofpara.size(); ++i) {
            biofilmDensity_[i] = biofpara[i].biofilm_density;
            maxGrowthRate_[i] = biofpara[i].max_growth_rate;
            halfVelocityCoeff_[i] = biofpara[i].half_velocity_coefficient;
            yieldCoeff_[i] = biofpara[i].yield_coefficient;
            decayCoeff_[i] = biofpara[i].decay_coefficient;
        }
    }
    else {
        throw std::runtime_error("BIOFPARA must be specified in BIOFILM runs\n");
    }
}
#endif

#define INSTANTIATE_TYPE(T)                                                                             \
    template struct BlackOilBiofilmParams<T>;                                                           \
    template void BlackOilBiofilmParams<T>::initFromState<false>(const EclipseState&); \
    template void BlackOilBiofilmParams<T>::initFromState<true>(const EclipseState&);

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm
