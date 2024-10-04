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

    using Toolbox = MathToolbox<Evaluation>;

    using TabulatedFunction = typename BlackOilBiofilmParams<Scalar>::TabulatedFunction;

    enum { gasCompIdx = FluidSystem::gasCompIdx };

    static constexpr unsigned conti0EqIdx = Indices::conti0EqIdx;
    static constexpr unsigned biofilmsConcentrationIdx = Indices::biofilmsConcentrationIdx;
    static constexpr unsigned contiBiofilmsEqIdx = Indices::contiBiofilmsEqIdx;
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
        if (!enableBiofilm)
            // Biofilm has been disabled at compile time
            return;

        VtkBlackOilBiofilmModule<TypeTag>::registerParameters();
    }

    /*!
     * \brief Register all biofilm specific VTK and ECL output modules.
     */
    static void registerOutputModules(Model& model,
                                      Simulator& simulator)
    {
        if (!enableBiofilm)
            // Biofilm has been disabled at compile time
            return;

        model.addOutputModule(new VtkBlackOilBiofilmModule<TypeTag>(simulator));
    }

    static bool eqApplies(unsigned eqIdx)
    {
        if (!enableBiofilm)
            return false;

        // True if biofilm equation applies
        return eqIdx == contiBiofilmsEqIdx;
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
        if (!enableBiofilm)
            return;

        const LhsEval solidBiofilm =
                Toolbox::template decay<LhsEval>(intQuants.referencePorosity())
                * Toolbox::template decay<LhsEval>(intQuants.biofilmsConcentration());
        storage[contiBiofilmsEqIdx] += solidBiofilm*1e-6;
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

            // if (upIdx == inIdx) min(1.0 - biofilmsConcentration_, 1.0)
            //     flux[contiBiofilmsEqIdx] =
            //             extQuants.volumeFlux(waterPhaseIdx)
            //             *up.biofilmsConcentration()*1e-8;
            // else
            //     flux[contiBiofilmsEqIdx] =
            //             extQuants.volumeFlux(waterPhaseIdx)
            //             *decay<Scalar>(up.biofilmsConcentration())*1e-8;

            if (upIdx == inIdx) 
                flux[contiBiofilmsEqIdx] =
                        extQuants.biofilmVolumeFlux()
                        *up.biofilmsConcentration()*1e-6;
            else
                flux[contiBiofilmsEqIdx] =
                        extQuants.biofilmVolumeFlux()
                        *decay<Scalar>(up.biofilmsConcentration())*1e-6;
}
    }

    static void addSource(RateVector& source,
                            const ElementContext& elemCtx,
                            unsigned dofIdx,
                            unsigned timeIdx)

    {
        if (!enableBiofilm)
            return;

        // Get biofilm parameters
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, dofIdx, timeIdx);
        Scalar rhob = biofilmDensity(satnumRegionIdx);
        Scalar mu = maxGrowthRate(satnumRegionIdx);
        Scalar Kn = halfVelocityCoeff(satnumRegionIdx);
        Scalar Y = yieldCoeff(satnumRegionIdx);
        Scalar kd = decayCoeff(satnumRegionIdx);

        // Convert Rsw to concentration to use in source term
        const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);
        const auto& fs = intQuants.fluidState();
        const auto& Rsw = fs.Rsw();
        const auto& rhow = fs.density(waterPhaseIdx);
        unsigned pvtRegionIndex = fs.pvtRegionIndex();

        const auto& xG = RswToMassFraction(pvtRegionIndex, Rsw);

        // Get saturation, porosity, biofilm concentration, and inverse Bg for convenience
        const Evaluation& poroRef = intQuants.referencePorosity();
        const Evaluation& poro = intQuants.porosity();
        const Evaluation& sw = fs.saturation(waterPhaseIdx);
        const Evaluation& cBiof = intQuants.biofilmsConcentration() * poroRef;
        Scalar rho_gRef = FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, pvtRegionIndex);

        // Calculate biofilm growth rate
        Evaluation kg = mu * (xG * rhow / (xG * rhow + Kn));
        if (xG < 0)
            kg = mu * (xG * rhow / Kn);

        // Compute source terms
        // Biofilm growth and decay rate
        source[contiBiofilmsEqIdx] += (kg - kd) * cBiof * 1e-6;

        // Biofilm consumption of dissolved gas is proportional to biofilm growth rate
        unsigned activeGasCompIdx = Indices::canonicalToActiveComponentIndex(gasCompIdx);
        source[activeGasCompIdx] -= cBiof * rhob * kg / (Y * rho_gRef);
    }

    static const TabulatedFunction& permporoTable(const ElementContext& elemCtx,
                                                  unsigned scvIdx,
                                                  unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.permporoTable_[satnumRegionIdx];
    }

    static const TabulatedFunction& permporoTable(unsigned satnumRegionIdx)
    {
        return params_.permporoTable_[satnumRegionIdx];
    }

    static const TabulatedFunction& pefactTable(unsigned satnumRegionIdx)
    {
        return params_.pefactTable_[satnumRegionIdx];
    }

    static const Scalar biofilmDensity(unsigned satnumRegionIdx)
    {
        return params_.biofilmDensity_[satnumRegionIdx];
    }

    static const Scalar maxGrowthRate(unsigned satnumRegionIdx)
    {
        return params_.maxGrowthRate_[satnumRegionIdx];
    }

    static const Scalar halfVelocityCoeff(unsigned satnumRegionIdx)
    {
        return params_.halfVelocityCoeff_[satnumRegionIdx];
    }

    static const Scalar yieldCoeff(unsigned satnumRegionIdx)
    {
        return params_.yieldCoeff_[satnumRegionIdx];
    }

    static const Scalar decayCoeff(unsigned satnumRegionIdx)
    {
        return params_.decayCoeff_[satnumRegionIdx];
    }

    static bool hasPefactTables()
    {
        if constexpr (enableBiofilm)
            return !params_.pefactTable_.empty();
        else
            return false;
    }

private:
    static BlackOilBiofilmParams<Scalar> params_;

    static Evaluation RswToMassFraction(unsigned regionIdx, const Evaluation& Rsw) {
        Scalar rho_wRef = FluidSystem::referenceDensity(FluidSystem::waterPhaseIdx, regionIdx);
        Scalar rho_gRef = FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, regionIdx);

        const Evaluation rho_oG = Rsw * rho_gRef;

        return rho_oG/(rho_wRef + rho_oG);
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
    static constexpr int biofilmsConcentrationIdx = Indices::biofilmsConcentrationIdx;


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
        biofilmsConcentration_ = priVars.makeEvaluation(biofilmsConcentrationIdx, timeIdx, linearizationType);
        const Evaluation porosityFactor  = min(1.0 - biofilmsConcentration_, 1.0); //phi/phi_0
        unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, timeIdx);
        const auto totVolume = elemCtx.simulator().model().dofTotalVolume(globalDofIdx);
        const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);
        biofilmMass_ = biofilmsConcentration_ * elemCtx.problem().biofilmDensity(dofIdx) * intQuants.referencePorosity();
        if (0 < totVolume*intQuants.referencePorosity()){
            biofilmDensity_ /= (totVolume*intQuants.referencePorosity());
        }

        const auto& permporoTable = BiofilmModule::permporoTable(elemCtx, dofIdx, timeIdx);

        permPoro_ = permporoTable.eval(porosityFactor);
        biofilmMobility_ = permPoro_;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            asImp_().mobility_[phaseIdx] *= permPoro_;
        }
    }

    const Evaluation& biofilmsConcentration() const
    { return biofilmsConcentration_; }

    const Evaluation& biofilmMass() const
    { return biofilmMass_; }

    const Evaluation& biofilmDensity() const
    { return biofilmDensity_; }

    const Evaluation& biofilmMobility() const
    { return biofilmMobility_; }

    const Evaluation& permPoro() const
    { return permPoro_; }


protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation biofilmsConcentration_;
    Evaluation biofilmMass_;
    Evaluation biofilmDensity_;
    Evaluation biofilmMobility_;
    Evaluation permPoro_;

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

    const Evaluation& biofilmMass() const
    { throw std::logic_error("biofilmMass() called but biofilms are disabled"); }

    const Evaluation& biofilmDensity() const
    { throw std::runtime_error("biofilmDensity() called but biofilms are disabled"); }

    const Evaluation& biofilmMobility() const
    { throw std::runtime_error("biofilmMobility() called but biofilms are disabled"); }

     const Evaluation& permPoro() const
    { throw std::logic_error("permPoro() called but biofilms are disabled"); }
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

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using ExtensiveQuantities = GetPropType<TypeTag, Properties::ExtensiveQuantities>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

    using Toolbox = MathToolbox<Evaluation>;

    static constexpr unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr int dimWorld = GridView::dimensionworld;

    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using DimEvalVector = Dune::FieldVector<Evaluation, dimWorld>;

public:
    /*!
     * \brief Method which calculates the volume flux of the biofilm using the
     *        gas pressure potential difference between cells and transmissibilities
     */
    template <class Dummy = bool> // we need to make this method a template to avoid
                                  // compiler errors if it is not instantiated!
    void updateBiofilmFluxTrans(const ElementContext& elemCtx,
                               unsigned scvfIdx,
                               unsigned timeIdx)
    {
        const ExtensiveQuantities& extQuants = asImp_();

        unsigned interiorDofIdx = extQuants.interiorIndex();
        unsigned exteriorDofIdx = extQuants.exteriorIndex();
        assert(interiorDofIdx != exteriorDofIdx);

        const auto& intQuantsIn = elemCtx.intensiveQuantities(interiorDofIdx, timeIdx);
        const auto& intQuantsEx = elemCtx.intensiveQuantities(exteriorDofIdx, timeIdx);

        unsigned I = elemCtx.globalSpaceIndex(interiorDofIdx, timeIdx);
        unsigned J = elemCtx.globalSpaceIndex(exteriorDofIdx, timeIdx);

        Scalar thpres = elemCtx.problem().thresholdPressure(I, J);
        Scalar trans = elemCtx.problem().transmissibility(elemCtx, interiorDofIdx, exteriorDofIdx);
        Scalar g = elemCtx.problem().gravity()[dimWorld - 1];

        Scalar zIn = elemCtx.problem().dofCenterDepth(elemCtx, interiorDofIdx, timeIdx);
        Scalar zEx = elemCtx.problem().dofCenterDepth(elemCtx, exteriorDofIdx, timeIdx);
        Scalar distZ = zIn - zEx;

        const Evaluation& rhoIn = intQuantsIn.biofilmDensity();
        Scalar rhoEx = Toolbox::value(intQuantsEx.biofilmDensity());
        const Evaluation& rhoAvg = rhoIn*0.5 + rhoEx*0.5;

        const Evaluation& pressureInterior = intQuantsIn.fluidState().pressure(gasPhaseIdx);
        Evaluation pressureExterior = Toolbox::value(intQuantsEx.fluidState().pressure(gasPhaseIdx));
        pressureExterior += distZ*g*rhoAvg;

        Evaluation pressureDiffBiofilm = pressureExterior - pressureInterior;
        if (std::abs(scalarValue(pressureDiffBiofilm)) > thpres) {
            if (pressureDiffBiofilm < 0.0)
                pressureDiffBiofilm += thpres;
            else
                pressureDiffBiofilm -= thpres;
        }
        else
            pressureDiffBiofilm = 0.0;

        if (pressureDiffBiofilm > 0.0) {
            biofilmUpstreamDofIdx_ = exteriorDofIdx;
            biofilmDownstreamDofIdx_ = interiorDofIdx;
        }
        else if (pressureDiffBiofilm < 0.0) {
            biofilmUpstreamDofIdx_ = interiorDofIdx;
            biofilmDownstreamDofIdx_ = exteriorDofIdx;
        }
        else {
            // pressure potential gradient is zero; force consistent upstream and
            // downstream indices over the intersection regardless of the side which it
            // is looked at.
            biofilmUpstreamDofIdx_ = std::min(interiorDofIdx, exteriorDofIdx);
            biofilmDownstreamDofIdx_ = std::max(interiorDofIdx, exteriorDofIdx);
            biofilmVolumeFlux_ = 0.0;
            return;
        }

        Scalar faceArea = elemCtx.stencil(timeIdx).interiorFace(scvfIdx).area();
        const IntensiveQuantities& up = elemCtx.intensiveQuantities(biofilmUpstreamDofIdx_, timeIdx);
        if (biofilmUpstreamDofIdx_ == interiorDofIdx)
            biofilmVolumeFlux_ =
                up.biofilmMobility()
                *(-trans/faceArea)
                *pressureDiffBiofilm;
        else
            biofilmVolumeFlux_ =
                scalarValue(up.biofilmMobility())
                *(-trans/faceArea)
                *pressureDiffBiofilm;
    }

    unsigned biofilmUpstreamIndex() const
    { return biofilmUpstreamDofIdx_; }

    unsigned biofilmDownstreamIndex() const
    { return biofilmDownstreamDofIdx_; }

    const Evaluation& biofilmVolumeFlux() const
    { return biofilmVolumeFlux_; }

private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation biofilmVolumeFlux_;
    unsigned biofilmUpstreamDofIdx_;
    unsigned biofilmDownstreamDofIdx_;

};

template <class TypeTag>
class BlackOilBiofilmExtensiveQuantities<TypeTag, false>
{
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;

public:

    void updateBiofilmFluxTrans(const ElementContext&,
                              unsigned,
                              unsigned)
    { }

    unsigned biofilmUpstreamIndex() const
    { throw std::runtime_error("biofilmUpstreamIndex() called but biofilms are disabled"); }

    unsigned biofilmDownstreamIndex() const
    { throw std::runtime_error("biofilmDownstreamIndex() called but biofilms are disabled"); }

    const Evaluation& biofilmVolumeFlux() const
    { throw std::runtime_error("biofilmVolumeFlux() called but biofilms are disabled"); }

};

} // namespace Opm

#endif
