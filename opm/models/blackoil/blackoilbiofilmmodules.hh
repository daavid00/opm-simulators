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
 * \brief Contains the classes required to extend the black-oil model by biofilm.
 */
#ifndef EWOMS_BLACK_OIL_BIOFILM_MODULE_HH
#define EWOMS_BLACK_OIL_BIOFILM_MODULE_HH

#include <dune/common/fvector.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/models/blackoil/blackoilbiofilmparams.hpp>
#include <opm/models/blackoil/blackoilproperties.hh>

#include <opm/models/io/vtkblackoilbiofilmmodule.hpp>

#include <opm/models/utils/propertysystem.hh>

#include <cstddef>
#include <stdexcept>

namespace Opm {
/*!
 * \ingroup BlackOil
 * \brief Contains the high level supplements required to extend the black oil
 *        model by biofilm.
 */
template <class TypeTag, bool enableBiofilmV = getPropValue<TypeTag, Properties::EnableBiofilm>()>
class BlackOilBiofilmModule
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using ExtensiveQuantities = GetPropType<TypeTag, Properties::ExtensiveQuantities>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Model = GetPropType<TypeTag, Properties::Model>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using Toolbox = MathToolbox<Evaluation>;

    using TabulatedFunction = typename BlackOilBiofilmParams<Scalar>::TabulatedFunction;

    enum { gasCompIdx = FluidSystem::gasCompIdx };

    static constexpr unsigned conti0EqIdx = Indices::conti0EqIdx;
    static constexpr unsigned biofilmConcentrationIdx = Indices::biofilmConcentrationIdx;
    static constexpr unsigned contiBiofilmEqIdx = Indices::contiBiofilmEqIdx;
    static constexpr unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;

    static constexpr unsigned enableBiofilm = enableBiofilmV;

    static constexpr unsigned numEq = getPropValue<TypeTag, Properties::NumEq>();

public:
    //! \brief Set parameters.
    static void setParams(BlackOilBiofilmParams<Scalar>&& params)
    {
        params_ = params;
    }

    /*!
     * \brief Register all run-time parameters for the black-oil biofilm module.
     */
    static void registerParameters()
    {
        if constexpr (enableBiofilm)
            VtkBlackOilBiofilmModule<TypeTag>::registerParameters();
    }

    /*!
     * \brief Register all biofilm specific VTK and ECL output modules.
     */
    static void registerOutputModules(Model& model,
                                      Simulator& simulator)
    {
        if constexpr (enableBiofilm)
            model.addOutputModule(new VtkBlackOilBiofilmModule<TypeTag>(simulator));
    }

    static bool eqApplies(unsigned eqIdx)
    {
        if constexpr (enableBiofilm)
            return eqIdx == contiBiofilmEqIdx;
        else
            false;
    }

    static Scalar eqWeight([[maybe_unused]] unsigned eqIdx)
    {
        assert(eqApplies(eqIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    // must be called after water storage is computed
    template <class LhsEval>
    static void addStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                           const IntensiveQuantities& intQuants)
    {
        if constexpr (enableBiofilm) {
            const LhsEval solidBiofilm =
                    Toolbox::template decay<LhsEval>(intQuants.referencePorosity())
                    * Toolbox::template decay<LhsEval>(intQuants.biofilmsConcentration());
            storage[contiBiofilmEqIdx] += solidBiofilm;
        }
    }

    template <class UpEval, class Eval, class IntensiveQuantities>
    static void addBiofilmFluxes_(RateVector& flux,
                               const Eval& volumeFlux,
                               const IntensiveQuantities& upFs)
    {
        if constexpr (enableBiofilm) {
            flux[contiBiofilmEqIdx] =
                decay<UpEval>(upFs.biofilmsConcentration())
                * volumeFlux*1e-1;
        }
    }

    static void computeFlux([[maybe_unused]] RateVector& flux,
                            [[maybe_unused]] const ElementContext& elemCtx,
                            [[maybe_unused]] unsigned scvfIdx,
                            [[maybe_unused]] unsigned timeIdx)

    {
        if constexpr (enableBiofilm) {
            const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

            const unsigned upIdx = extQuants.upstreamIndex(waterPhaseIdx);
            unsigned inIdx = extQuants.interiorIndex();
            const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);

            if (upIdx == inIdx) 
                flux[contiBiofilmEqIdx] =
                        extQuants.biofilmVolumeFlux()
                        *up.biofilmsConcentration();
            else
                flux[contiBiofilmEqIdx] =
                        extQuants.biofilmVolumeFlux()
                        *decay<Scalar>(up.biofilmsConcentration());
        }
    }

    static void addSource(RateVector& source,
                          const Problem& problem,
                          const IntensiveQuantities& intQuants,
                          unsigned globalSpaceIdex)
    {
        if constexpr (enableBiofilm) {
            unsigned satnumIdx = problem.satnumRegionIndex(globalSpaceIdex);
            Scalar rhob = densityBiofilm(satnumIdx);
            Scalar mu = maximumGrowthRate(satnumIdx);
            Scalar Kn = halfVelocityOxygen(satnumIdx);
            Scalar Y = yieldGrowthCoefficient(satnumIdx);
            Scalar kd = microbialDeathRate(satnumIdx);

            // Convert Rsw to concentration to use in source term
            const auto& fs = intQuants.fluidState();
            const auto& Rsw = fs.Rsw();
            const auto& rhow = fs.density(waterPhaseIdx);
            unsigned pvtRegionIndex = fs.pvtRegionIndex();

            const auto& xG = RswToMassFraction(pvtRegionIndex, Rsw);

            // Get saturation, porosity, biofilm concentration, and inverse Bg for convenience
            const Evaluation& poroRef = intQuants.referencePorosity();
            const Evaluation& poro = intQuants.porosity();
            const Evaluation& cBiof = intQuants.biofilmsConcentration() * poroRef;
            Scalar rho_gRef = FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, pvtRegionIndex);

            // Calculate biofilm growth rate
            Evaluation kg = mu * (xG * rhow / (xG * rhow + Kn));
            if (xG < 0)
                kg = mu * (xG * rhow / Kn);

            // Compute source terms
            // Biofilm growth and decay rate
            source[contiBiofilmEqIdx] += (kg - kd) * cBiof;

            // Biofilm consumption of dissolved gas is proportional to biofilm growth rate
            unsigned activeGasCompIdx = Indices::canonicalToActiveComponentIndex(gasCompIdx);
            source[activeGasCompIdx] -= cBiof * rhob * kg / (Y * rho_gRef);
        }
    }

    static void addSource([[maybe_unused]] RateVector& source,
                          [[maybe_unused]] const ElementContext& elemCtx,
                          [[maybe_unused]] unsigned dofIdx,
                          [[maybe_unused]] unsigned timeIdx)
    {
        if constexpr (enableBiofilm) {
            const auto& problem = elemCtx.problem();
            const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);
            addSource(source, problem, intQuants, dofIdx);
        }
    }

    static const TabulatedFunction& pcfactTable(unsigned satnumRegionIdx)
    {
        return params_.pcfactTable_[satnumRegionIdx];
    }

    static bool hasPcfactTables()
    {
        if constexpr (enableBiofilm)
            return !params_.pcfactTable_.empty();
        else
            return false;
    }

    static const Scalar densityBiofilm(unsigned satnumRegionIdx)
    {
        return params_.densityBiofilm_[satnumRegionIdx];
    }

    static const Scalar halfVelocityOxygen(unsigned satnumRegionIdx)
    {
        return params_.halfVelocityOxygen_[satnumRegionIdx];
    }

    static const Scalar maximumGrowthRate(unsigned satnumRegionIdx)
    {
        return params_.maximumGrowthRate_[satnumRegionIdx];
    }

    static const Scalar microbialDeathRate(unsigned satnumRegionIdx)
    {
        return params_.microbialDeathRate_[satnumRegionIdx];
    }

    static const Scalar yieldGrowthCoefficient(unsigned satnumRegionIdx)
    {
        return params_.yieldGrowthCoefficient_[satnumRegionIdx];
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
    static BlackOilBiofilmParams<Scalar> params_;

    static Evaluation RswToMassFraction(unsigned regionIdx, const Evaluation& Rsw) {
        Scalar rho_wRef = FluidSystem::referenceDensity(FluidSystem::waterPhaseIdx, regionIdx);
        Scalar rho_gRef = FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, regionIdx);

        const Evaluation rho_oG = Rsw * rho_gRef;

        return rho_oG / (rho_wRef + rho_oG);
    }
};


template <class TypeTag, bool enableBiofilmV>
BlackOilBiofilmParams<typename BlackOilBiofilmModule<TypeTag, enableBiofilmV>::Scalar>
BlackOilBiofilmModule<TypeTag, enableBiofilmV>::params_;

/*!
 * \ingroup BlackOil
 * \class Opm::BlackOilBiofilmIntensiveQuantities
 *
 * \brief Provides the volumetric quantities required for the equations needed by the
 *        Biofilm extension of the black-oil model.
 */
template <class TypeTag, bool enableBiofilmV = getPropValue<TypeTag, Properties::EnableBiofilm>()>
class BlackOilBiofilmIntensiveQuantities
{
    using Implementation = GetPropType<TypeTag, Properties::IntensiveQuantities>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    using BiofilmModule = BlackOilBiofilmModule<TypeTag>;

    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    static constexpr int biofilmConcentrationIdx = Indices::biofilmConcentrationIdx;


public:
    /*!
     * \brief Update the intensive properties needed to handle biofilm from the
     *        primary variables
     *
     */
    void biofilmPropertiesUpdate_(const ElementContext& elemCtx,
                                  unsigned dofIdx,
                                  unsigned timeIdx)
    {
        const auto linearizationType = elemCtx.linearizationType();
        const PrimaryVariables& priVars = elemCtx.primaryVars(dofIdx, timeIdx);

        // update biofilm concentration from primary variables
        biofilmsConcentration_ = priVars.makeEvaluation(biofilmConcentrationIdx, timeIdx, linearizationType);
        const Evaluation porosityFactor  = min(1.0 - biofilmsConcentration_, 1.0); //phi/phi_0
        unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, timeIdx);
        const auto totVolume = elemCtx.simulator().model().dofTotalVolume(globalDofIdx);
        const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, dofIdx, timeIdx);
        biofilmDensity_ =  BiofilmModule::densityBiofilm(satnumRegionIdx);
        biofilmMass_ = biofilmsConcentration_ * intQuants.referencePorosity() * biofilmDensity_;
        const auto& permfactTable = BiofilmModule::permfactTable(satnumRegionIdx);
        permFactor_ = permfactTable.eval(porosityFactor, /*extrapolation=*/true);
        biofilmMobility_ = permFactor_;
        // for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
        //     if (!FluidSystem::phaseIsActive(phaseIdx))
        //         continue;

        //     asImp_().mobility_[phaseIdx] *= permFactor_;
        // }
    }

    const Evaluation& biofilmsConcentration() const
    { return biofilmsConcentration_; }

    const Evaluation& biofilmMasss() const
    { return biofilmMass_; }

    const Evaluation& biofilmDensity() const
    { return biofilmDensity_; }

    const Evaluation& biofilmMobility() const
    { return biofilmMobility_; }

    const Evaluation& permFactor() const
    { return permFactor_; }

protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation biofilmsConcentration_;
    Evaluation biofilmMass_;
    Evaluation biofilmDensity_;
    Evaluation biofilmMobility_;
    Evaluation permFactor_;

};

template <class TypeTag>
class BlackOilBiofilmIntensiveQuantities<TypeTag, false>
{
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    void biofilmPropertiesUpdate_(const ElementContext&,
                                  unsigned,
                                  unsigned)
    { }

    const Evaluation& biofilmsConcentration() const
    { throw std::logic_error("biofilmsConcentration() called but biofilms are disabled"); }

    const Evaluation& biofilmMasss() const
    { throw std::logic_error("biofilmMasss() called but biofilms are disabled"); }

    const Evaluation& biofilmDensity() const
    { throw std::runtime_error("biofilmDensity() called but biofilms are disabled"); }

    const Evaluation& biofilmMobility() const
    { throw std::runtime_error("biofilmMobility() called but biofilms are disabled"); }

    const Evaluation& permFactor() const
    { throw std::logic_error("permFactor() called but biofilms are disabled"); }
};

/*!
 * \ingroup BlackOil
 * \class Opm::BlackOilBiofilmExtensiveQuantities
 *
 * \brief Provides the Biofilm specific extensive quantities to the generic black-oil
 *        module's extensive quantities.
 */
template <class TypeTag, bool enableBiofilmV = getPropValue<TypeTag, Properties::EnableBiofilm>()>
class BlackOilBiofilmExtensiveQuantities
{
    using Implementation = GetPropType<TypeTag, Properties::ExtensiveQuantities>;
};

template <class TypeTag>
class BlackOilBiofilmExtensiveQuantities<TypeTag, false>{};

} // namespace Opm

#endif
