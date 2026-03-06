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
/*!
 * \file
 *
 * \brief Contains the classes required to extend the black-oil model by particles.
 */
#ifndef OPM_BLACK_OIL_PARTICLE_MODULE_HH
#define OPM_BLACK_OIL_PARTICLE_MODULE_HH

#include <dune/common/fvector.hh>

#include <opm/common/utility/gpuDecorators.hpp>

#include <opm/models/blackoil/blackoilparticleparams.hpp>
#include <opm/models/blackoil/blackoilproperties.hh>

#include <opm/models/io/vtkblackoilparticlemodule.hpp>

#include <cmath>
#include <memory>
#include <numeric>
#include <stdexcept>

namespace Opm {
/*!
 * \ingroup BlackOil
 * \brief Contains the high level supplements required to extend the black oil
 *        model by particles.
 */
template <class TypeTag, bool enableParticleV = getPropValue<TypeTag, Properties::EnableParticle>()>
class BlackOilParticleModule
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Model = GetPropType<TypeTag, Properties::Model>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using Toolbox = MathToolbox<Evaluation>;

    using TabulatedFunction = typename BlackOilParticleParams<Scalar>::TabulatedFunction;

    static constexpr unsigned particleConcentrationIdx = Indices::particleConcentrationIdx;
    static constexpr unsigned particleVolumeFractionIdx = Indices::particleVolumeFractionIdx;
    static constexpr unsigned contiSuspendedParticleEqIdx = Indices::contiSuspendedParticleEqIdx;
    static constexpr unsigned contiRetainedParticleEqIdx = Indices::contiRetainedParticleEqIdx;
    static constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;

    static constexpr unsigned enableParticle = enableParticleV;

    static constexpr unsigned numEq = getPropValue<TypeTag, Properties::NumEq>();

public:
    //! \brief Set parameters.
    static void setParams(BlackOilParticleParams<Scalar>&& params)
    {
        params_ = params;
    }

    /*!
     * \brief Register all run-time parameters for the black-oil particles module.
     */
    static void registerParameters()
    {
        if constexpr (enableParticle)
            VtkBlackOilParticleModule<TypeTag>::registerParameters();
    }

    /*!
     * \brief Register all particle specific VTK and ECL output modules.
     */
    static void registerOutputModules(Model& model,
                                      Simulator& simulator)
    {
        if constexpr (enableParticle)
            model.addOutputModule(std::make_unique<VtkBlackOilParticleModule<TypeTag>>(simulator));
    }

    static bool eqApplies(unsigned eqIdx)
    {
        if constexpr (enableParticle)
            return eqIdx == contiSuspendedParticleEqIdx || eqIdx == contiRetainedParticleEqIdx;
        else
            return false;
    }

    static Scalar eqWeight([[maybe_unused]] unsigned eqIdx)
    {
        assert(eqApplies(eqIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    // must be called after water storage is computed
    template <class StorageType>
    OPM_HOST_DEVICE static void addStorage(StorageType& storage,
                                           const IntensiveQuantities& intQuants)
    {
        using LhsEval = typename StorageType::value_type;
        if constexpr (enableParticle) {
            const auto& fs = intQuants.fluidState();
            LhsEval surfaceVolumeWater = Toolbox::template decay<LhsEval>(fs.saturation(waterPhaseIdx)) *
                                         Toolbox::template decay<LhsEval>(fs.invB(waterPhaseIdx)) *
                                         Toolbox::template decay<LhsEval>(intQuants.porosity());
            // avoid singular matrix if no water is present
            surfaceVolumeWater = max(surfaceVolumeWater, 1e-10);
            // suspended particles in water phase
            const LhsEval accumulationSuspendedParticle = surfaceVolumeWater * Toolbox::template decay<LhsEval>(intQuants.particleConcentration());
            storage[contiSuspendedParticleEqIdx] += accumulationSuspendedParticle;
            // retained particles
            const LhsEval accumulationRetainedParticle = Toolbox::template decay<LhsEval>(intQuants.particleVolumeFraction());
            storage[contiRetainedParticleEqIdx] += accumulationRetainedParticle;
        }
    }

    template <class UpEval>
    static void addParticleFluxes_(RateVector& flux,
                                     unsigned phaseIdx,
                                     const Evaluation& volumeFlux,
                                     const IntensiveQuantities& upFs)
    {
        if (phaseIdx == waterPhaseIdx) {
            if constexpr (enableParticle) {
                flux[contiSuspendedParticleEqIdx] =
                    decay<UpEval>(upFs.particleConcentration())
                    * decay<UpEval>(upFs.fluidState().invB(waterPhaseIdx))
                    * volumeFlux;
            }
        }
    }

    static void computeFlux([[maybe_unused]] RateVector& flux,
                            [[maybe_unused]] const ElementContext& elemCtx,
                            [[maybe_unused]] unsigned scvfIdx,
                            [[maybe_unused]] unsigned timeIdx)
    {
        if constexpr (enableParticle) {
            const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
            unsigned focusIdx = elemCtx.focusDofIndex();
            unsigned upIdx = extQuants.upstreamIndex(waterPhaseIdx);
            flux[contiSuspendedParticleEqIdx] = 0.0;
            if (upIdx == focusIdx)
                addParticleFluxes_<Evaluation>(flux, elemCtx, scvfIdx, timeIdx);
            else
                addParticleFluxes_<Scalar>(flux, elemCtx, scvfIdx, timeIdx);
        }
    }

    template <class UpstreamEval>
    static void addParticleFluxes_(RateVector& flux,
                                     const ElementContext& elemCtx,
                                     unsigned scvfIdx,
                                     unsigned timeIdx)
    {
        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
        unsigned upIdx = extQuants.upstreamIndex(waterPhaseIdx);
        const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
        const auto& volFlux = extQuants.volumeFlux(waterPhaseIdx);
        addParticleFluxes_<UpstreamEval>(flux, waterPhaseIdx, volFlux, up);
    }

    static void addSource(RateVector& source,
                          const Problem& problem,
                          const IntensiveQuantities& intQuants,
                          unsigned globalSpaceIdex)
    {
        if constexpr (enableParticle) {
            const auto b = intQuants.fluidState().invB(waterPhaseIdx);
            unsigned satnumIdx = problem.satnumRegionIndex(globalSpaceIdex);
            Scalar rho_p = densityParticle(satnumIdx);
            Scalar k_a = retainmentRate(satnumIdx);
            Scalar k_d = releaseRate(satnumIdx);

            // compute the processes
            source[contiSuspendedParticleEqIdx] += rho_p * intQuants.particleVolumeFraction() * k_d
                                                    - intQuants.particleConcentration() * intQuants.porosity() * b * k_a;

            source[contiRetainedParticleEqIdx] += k_a * intQuants.particleConcentration() * intQuants.porosity() * b / rho_p 
                                                   - intQuants.particleVolumeFraction() * k_d;
        }
    }

    static void addSource([[maybe_unused]] RateVector& source,
                          [[maybe_unused]] const ElementContext& elemCtx,
                          [[maybe_unused]] unsigned dofIdx,
                          [[maybe_unused]] unsigned timeIdx)
    {
    }

    static const Scalar densityParticle(unsigned satnumRegionIdx)
    {
        return params_.densityParticle_[satnumRegionIdx];
    }

    static const Scalar retainmentRate(unsigned satnumRegionIdx)
    {
        return params_.particleRetainingRate_[satnumRegionIdx];
    }

    static const Scalar releaseRate(unsigned satnumRegionIdx)
    {
        return params_.particleReleasingRate_[satnumRegionIdx];
    }

    static const TabulatedFunction& permfactTable(const ElementContext& elemCtx,
                                                  unsigned scvIdx,
                                                  unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.permfactTable_[satnumRegionIdx];
    }

    static const TabulatedFunction& permfactTable(unsigned satnumRegionIdx)
    {
        return params_.permfactTable_[satnumRegionIdx];
    }

private:
    static BlackOilParticleParams<Scalar> params_;

};


template <class TypeTag, bool enableParticleV>
BlackOilParticleParams<typename BlackOilParticleModule<TypeTag, enableParticleV>::Scalar>
BlackOilParticleModule<TypeTag, enableParticleV>::params_;

/*!
 * \ingroup BlackOil
 * \class Opm::BlackOilParticleIntensiveQuantities
 *
 * \brief Provides the volumetric quantities required for the equations needed by the
 *        particle extension of the black-oil model.
 */
template <class TypeTag, bool enableParticleV = getPropValue<TypeTag, Properties::EnableParticle>()>
class BlackOilParticleIntensiveQuantities
{
    using Implementation = GetPropType<TypeTag, Properties::IntensiveQuantities>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    using ParticleModule = BlackOilParticleModule<TypeTag>;

    static constexpr int particleConcentrationIdx = Indices::particleConcentrationIdx;
    static constexpr int particleVolumeFractionIdx = Indices::particleVolumeFractionIdx;
    static constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;

public:

    /*!
     * \brief Update the intensive properties needed to handle particles from the
     *        primary variables
     *
     */
    void particlePropertiesUpdate_(const ElementContext& elemCtx,
                                   unsigned dofIdx,
                                   unsigned timeIdx)
    {
        auto& fs = asImp_().fluidState_;
        const unsigned pvtRegionIdx = asImp_().pvtRegionIndex();
        const auto linearizationType = elemCtx.linearizationType();
        const PrimaryVariables& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        const Scalar referencePorosity_ = elemCtx.problem().referencePorosity(dofIdx, timeIdx);
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, dofIdx, timeIdx);
        const auto& rhoPart = ParticleModule::densityParticle(satnumRegionIdx);

        particleConcentration_ = priVars.makeEvaluation(particleConcentrationIdx, timeIdx, linearizationType);
        particleVolumeFraction_ = priVars.makeEvaluation(particleVolumeFractionIdx, timeIdx, linearizationType);
        particleRetainedMass_ = particleVolumeFraction_ * rhoPart;
        const Evaluation poroFact = min(1.0 - (particleVolumeFraction_) /
                                              (referencePorosity_), 1.0); //phi/phi_0

        const auto& permfactTable = ParticleModule::permfactTable(satnumRegionIdx);
        permFactor_ = permfactTable.eval(poroFact, /*extrapolation=*/true);

        // modify the fluid density using the particle concentration
        const auto& rhoWatRef = FluidSystem::referenceDensity(waterPhaseIdx, pvtRegionIdx);
        Evaluation rhoWatEff = rhoWatRef + particleConcentration_ * (1.0 - rhoWatRef / rhoPart);
        fs.setDensity(waterPhaseIdx, rhoWatEff);
    }

    const Evaluation& particleConcentration() const
    { return particleConcentration_; }

    const Evaluation& particleVolumeFraction() const
    { return particleVolumeFraction_; }

    const Evaluation particleRetainedMass() const
    { return particleRetainedMass_; }

    const Evaluation& permFactor() const
    { return permFactor_; }

protected:
    Evaluation particleConcentration_;
    Evaluation particleVolumeFraction_;
    Evaluation particleRetainedMass_;
    Evaluation permFactor_;

    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

};

template <class TypeTag>
class BlackOilParticleIntensiveQuantities<TypeTag, false>
{
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    void particlePropertiesUpdate_(const ElementContext&,
                                     unsigned,
                                     unsigned)
    {}

    const Evaluation& particleConcentration() const
    { throw std::logic_error("particleConcentration() called but PARTICLE is disabled"); }

    const Evaluation& particleVolumeFraction() const
    { throw std::logic_error("particleVolumeFraction() called but PARTICLE is disabled"); }

    const Evaluation& particleRetainedMass() const
    { throw std::logic_error("particleRetainedMass() called but PARTICLE is disabled"); }

    const Evaluation& permFactor() const
    { throw std::logic_error("permFactor() called but PARTICLE is disabled"); }
};

/*!
 * \ingroup BlackOil
 * \class Opm::BlackOilParticleExtensiveQuantities
 *
 * \brief Provides the particles specific extensive quantities to the generic black-oil
 *        module's extensive quantities.
 */
template <class TypeTag, bool enableParticleV = getPropValue<TypeTag, Properties::EnableParticle>()>
class BlackOilParticleExtensiveQuantities
{
};

template <class TypeTag>
class BlackOilParticleExtensiveQuantities<TypeTag, false>{};

} // namespace Opm

#endif
