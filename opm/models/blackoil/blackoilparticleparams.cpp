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
#include <opm/models/blackoil/blackoilparticleparams.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/ParticleTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PermfactTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SimpleTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>

#include <algorithm>
#include <stdexcept>
#include <type_traits>

namespace Opm {

template<class Scalar>
template<bool enableParticle>
void BlackOilParticleParams<Scalar>::
initFromState(const EclipseState& eclState)
{
    // some sanity checks:
    // if particles are enabled, the PARTICLE keyword must be present, and
    // if particles are disabled, the keyword must not be present.
    if constexpr (enableParticle) {
        if (!eclState.runspec().particle()) {
            throw std::runtime_error("Non-trivial particle effects requested at compile time, but "
                                     "the deck does not contain the PARTICLE keyword");
        }
    }
    else {
        if (eclState.runspec().particle()) {
            throw std::runtime_error("Particle effects disabled at compile time, but the deck "
                                     "contains the PARTICLE keyword");
        }
    }

    if (!eclState.runspec().particle())
        return; // particles are supposed to be disabled

    const auto& tableManager = eclState.getTableManager();
    unsigned numSatRegions = tableManager.getTabdims().getNumSatTables();

    // initialize the objects which deal with the particle parameters
    const TableContainer& particleTables = tableManager.getParticleTables();
    if (particleTables.empty()) {
        throw std::runtime_error("PARTICLE require the PARTPARA keyword");
    }
    densityParticle_.resize(numSatRegions);
    particleRetainingRate_.resize(numSatRegions);
    particleReleasingRate_.resize(numSatRegions);
    for (unsigned stnRegionIdx = 0; stnRegionIdx < numSatRegions; ++stnRegionIdx) {
        const ParticleTable& particleTable = particleTables.getTable<ParticleTable>(stnRegionIdx);
        densityParticle_[stnRegionIdx] = particleTable.getDensityParticle().front();
        particleRetainingRate_[stnRegionIdx] = particleTable.getParticleAttachmentRate().front();
        particleReleasingRate_[stnRegionIdx] = particleTable.getParticleDetachmentRate().front();
    }

    const TableContainer& permfactTables = tableManager.getPermfactTables();
    if (permfactTables.empty()) {
        throw std::runtime_error("PARTICLE require the PERMFACT keyword");
    }
    permfactTable_.resize(numSatRegions);
    for (std::size_t i = 0; i < permfactTables.size(); ++i) {
        const PermfactTable& permfactTable = permfactTables.getTable<PermfactTable>(i);
        permfactTable_[i].setXYContainers(permfactTable.getPorosityChangeColumn(), permfactTable.getPermeabilityMultiplierColumn());
    }
}

#define INSTANTIATE_TYPE(T)                                                              \
    template struct BlackOilParticleParams<T>;                                          \
    template void BlackOilParticleParams<T>::initFromState<false>(const EclipseState&); \
    template void BlackOilParticleParams<T>::initFromState<true>(const EclipseState&);

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm
