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
 * \copydoc Opm::VtkBlackOilBiofilmModule
 */
#ifndef OPM_VTK_BLACK_OIL_BIOFILM_MODULE_HPP
#define OPM_VTK_BLACK_OIL_BIOFILM_MODULE_HPP

#include <dune/common/fvector.hh>

#include <opm/material/densead/Math.hpp>

#include <opm/models/blackoil/blackoilproperties.hh>

#include <opm/models/discretization/common/fvbaseparameters.hh>

#include <opm/models/io/baseoutputmodule.hh>
#include <opm/models/io/vtkblackoilbiofilmparams.hpp>
#include <opm/models/io/vtkmultiwriter.hh>

#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>

namespace Opm {
/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for the biofilm model's related quantities.
 */
template <class TypeTag>
class VtkBlackOilBiofilmModule : public BaseOutputModule<TypeTag>
{
    using ParentType = BaseOutputModule<TypeTag>;

    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    static const int vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkFormat>;

    enum { enableBiofilm = getPropValue<TypeTag, Properties::EnableBiofilm>() };

    using ScalarBuffer = typename ParentType::ScalarBuffer;

public:
    VtkBlackOilBiofilmModule(const Simulator& simulator)
        : ParentType(simulator)
    {
        if constexpr (enableBiofilm) {
            params_.read();
        }
    }

    /*!
     * \brief Register all run-time parameters for the multi-phase VTK output
     * module.
     */
    static void registerParameters()
    {
        if constexpr (enableBiofilm) {
            VtkBlackoilBiofilmParams::registerParameters();
        }
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers()
    {
        if constexpr (enableBiofilm) {
            if (!Parameters::Get<Parameters::EnableVtkOutput>()) {
                return;
            }

            if (params_.biofilmsConcentrationOutput_) {
                this->resizeScalarBuffer_(biofilmsConcentration_);
            }
        }
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant for
     *        an element
     */
    void processElement(const ElementContext& elemCtx)
    {
        if constexpr (enableBiofilm) {
            if (!Parameters::Get<Parameters::EnableVtkOutput>()) {
                return;
            }

            for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
                const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);
                unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);

                if (params_.biofilmsConcentrationOutput_) {
                    biofilmsConcentration_[globalDofIdx] =
                        scalarValue(intQuants.biofilmsConcentration());
                }
            }
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(BaseOutputWriter& baseWriter)
    {
        if constexpr (enableBiofilm) {
            VtkMultiWriter* vtkWriter = dynamic_cast<VtkMultiWriter*>(&baseWriter);
            if (!vtkWriter) {
                return;
            }

            if (params_.biofilmsConcentrationOutput_) {
                this->commitScalarBuffer_(baseWriter, "biofilm fraction", biofilmsConcentration_);
            }
        }
    }

private:
    VtkBlackoilBiofilmParams params_{};
    ScalarBuffer biofilmsConcentration_{};
};

} // namespace Opm

#endif // OPM_VTK_BLACKOIL_BIOFILM_MODULE_HPP
