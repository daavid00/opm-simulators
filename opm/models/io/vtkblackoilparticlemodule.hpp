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
 * \copydoc Opm::VtkBlackOilParticleModule
 */
#ifndef OPM_VTK_BLACK_OIL_PARTICLE_MODULE_HPP
#define OPM_VTK_BLACK_OIL_PARTICLE_MODULE_HPP

#include <dune/common/fvector.hh>

#include <opm/material/densead/Math.hpp>

#include <opm/models/blackoil/blackoilproperties.hh>

#include <opm/models/discretization/common/fvbaseparameters.hh>

#include <opm/models/io/baseoutputmodule.hh>
#include <opm/models/io/vtkblackoilparticleparams.hpp>
#include <opm/models/io/vtkmultiwriter.hh>

#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>

namespace Opm {
/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for the particles model's related quantities.
 */
template <class TypeTag>
class VtkBlackOilParticleModule : public BaseOutputModule<TypeTag>
{
    using ParentType = BaseOutputModule<TypeTag>;

    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    static constexpr auto vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkFormat>;

    enum { enableParticle = getPropValue<TypeTag, Properties::EnableParticle>() };

    using BufferType = typename ParentType::BufferType;
    using ScalarBuffer = typename ParentType::ScalarBuffer;

public:
    explicit VtkBlackOilParticleModule(const Simulator& simulator)
        : ParentType(simulator)
    {
        if constexpr (enableParticle) {
            params_.read();
        }
    }

    /*!
     * \brief Register all run-time parameters for the multi-phase VTK output
     * module.
     */
    static void registerParameters()
    {
        if constexpr (enableParticle) {
            VtkBlackOilParticleParams::registerParameters();
        }
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers() override
    {
        if constexpr (enableParticle) {
            if (!Parameters::Get<Parameters::EnableVtkOutput>()) {
                return;
            }

            if (params_.particleConcentrationOutput_) {
                this->resizeScalarBuffer_(particleConcentration_, BufferType::Dof);
            }
            if (params_.particleVolumeFractionOutput_) {
                this->resizeScalarBuffer_(particleVolumeFraction_, BufferType::Dof);
            }
        }
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant for
     *        an element
     */
    void processElement(const ElementContext& elemCtx) override
    {
        if constexpr (enableParticle) {
            if (!Parameters::Get<Parameters::EnableVtkOutput>()) {
                return;
            }

            for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
                const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);
                const unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);

                if (params_.particleConcentrationOutput_) {
                    particleConcentration_[globalDofIdx] =
                        scalarValue(intQuants.particleConcentration());
                }
                if (params_.particleVolumeFractionOutput_) {
                    particleVolumeFraction_[globalDofIdx] =
                        scalarValue(intQuants.particleVolumeFraction());
                }
            }
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(BaseOutputWriter& baseWriter) override
    {
        if constexpr (enableParticle) {
            if (!dynamic_cast<VtkMultiWriter*>(&baseWriter)) {
                return;
            }

            if (params_.particleConcentrationOutput_) {
                this->commitScalarBuffer_(baseWriter, "particle concentration",
                                          particleConcentration_, BufferType::Dof);
            }
            if (params_.particleVolumeFractionOutput_) {
                this->commitScalarBuffer_(baseWriter, "particle volume fraction",
                                          particleVolumeFraction_, BufferType::Dof);
            }
        }
    }

private:
    VtkBlackOilParticleParams params_{};
    ScalarBuffer particleConcentration_{};
    ScalarBuffer particleVolumeFraction_{};
};

} // namespace Opm

#endif // OPM_VTK_BLACK_OIL_PARTICLE_MODULE_HPP
