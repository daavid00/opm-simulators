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
#include <opm/simulators/flow/LogOutputHelper.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>

namespace Opm {

template<class Scalar>
LogOutputHelper<Scalar>::LogOutputHelper(const EclipseState& eclState,
                                         const Schedule& schedule,
                                         const SummaryState& summaryState)
    : eclState_(eclState)
    , schedule_(schedule)
    , summaryState_(summaryState)
{}

template<class Scalar>
void LogOutputHelper<Scalar>::
cumulative(const std::size_t reportStepNum,
           std::function<bool(const std::string&)> isDefunct) const
{
    std::vector<Scalar> tmp_values(WellCumDataType::numWCValues, 0.0);
    std::vector<std::string> tmp_names(WellCumDataType::numWCNames, "");
    this->outputCumulativeReport_(tmp_values, tmp_names);

    const auto& st = summaryState_;
    for (const auto& gname : schedule_.groupNames()) {
        auto gName = static_cast<std::string>(gname);
        auto get = [&st, &gName](const std::string& vector)
        {
            const auto key = vector + ':' + gName;
            return st.has(key) ? st.get(key) : 0.0;
        };

        tmp_names[0] = gname;

        if (tmp_names[0] == "FIELD") {
            tmp_values[2] = st.get("FOPT", 0.0); // WellCumDataType::OilProd
            tmp_values[3] = st.get("FWPT", 0.0); // WellCumDataType::WaterProd
            tmp_values[4] = st.get("FGPT", 0.0); // WellCumDataType::GasProd
            tmp_values[5] = st.get("FVPT", 0.0); // WellCumDataType::FluidResVolProd
            tmp_values[6] = st.get("FOIT", 0.0); // WellCumDataType::OilInj
            tmp_values[7] = st.get("FWIT", 0.0); // WellCumDataType::WaterInj
            tmp_values[8] = st.get("FGIT", 0.0); // WellCumDataType::GasInj
            tmp_values[9] = st.get("FVIT", 0.0); // WellCumDataType::FluidResVolInj
        } else {
            tmp_values[2] = get("GOPT"); // WellCumDataType::OilProd
            tmp_values[3] = get("GWPT"); // WellCumDataType::WaterProd
            tmp_values[4] = get("GGPT"); // WellCumDataType::GasProd
            tmp_values[5] = get("GVPT"); // WellCumDataType::FluidResVolProd
            tmp_values[6] = get("GOIT"); // WellCumDataType::OilInj
            tmp_values[7] = get("GWIT"); // WellCumDataType::WaterInj
            tmp_values[8] = get("GGIT"); // WellCumDataType::GasInj
            tmp_values[9] = get("GVIT"); // WellCumDataType::FluidResVolInj
        }

        this->outputCumulativeReport_(tmp_values, tmp_names);
    }

    for (const auto& wname : schedule_.wellNames(reportStepNum)) {
        // don't bother with wells not on this process
        if (isDefunct(wname)) {
            continue;
        }

        const auto& well = schedule_.getWell(wname, reportStepNum);
        tmp_names[0] = wname; // WellCumDataType::WellName
        auto wName = static_cast<std::string>(wname);
        auto get = [&st, &wName](const std::string& vector)
        {
            const auto key = vector + ':' + wName;
            return st.has(key) ? st.get(key) : 0.0;
        };

        if (well.isInjector()) {
            const auto& controls = well.injectionControls(st);
            const auto ctlMode = controls.cmode;
            const auto injType = controls.injector_type;
            using CMode = ::Opm::Well::InjectorCMode;
            using WType = ::Opm::InjectorType;

            auto ftype = [](const auto wtype) -> std::string
            {
                switch (wtype) {
                case WType::OIL:   return "Oil";
                case WType::WATER: return "Wat";
                case WType::GAS:   return "Gas";
                case WType::MULTI: return "Multi";
                default: return "";
                }
            };

            auto fctl = [](const auto wmctl) -> std::string
            {
                switch (wmctl) {
                case CMode::RATE: return "RATE";
                case CMode::RESV: return "RESV";
                case CMode::THP:  return "THP";
                case CMode::BHP:  return "BHP";
                case CMode::GRUP: return "GRUP";
                default: return "";
                }
            };

            tmp_names[1] = "INJ"; // WellCumDataType::WellType
            const auto flowctl = fctl(ctlMode);
            if (flowctl == "RATE") { // WellCumDataType::WellCTRL
                const auto flowtype = ftype(injType);
                if (flowtype == "Oil") {
                    tmp_names[2] = "ORAT";
                } else if (flowtype == "Wat") {
                    tmp_names[2] = "WRAT";
                } else if (flowtype == "Gas") {
                    tmp_names[2] = "GRAT";
                }
            } else {
                tmp_names[2] = flowctl;
            }
        } else if (well.isProducer()) {
            const auto& controls = well.productionControls(st);
            using CMode = ::Opm::Well::ProducerCMode;

            auto fctl = [](const auto wmctl) -> std::string
            {
                switch (wmctl) {
                case CMode::ORAT: return "ORAT";
                case CMode::WRAT: return "WRAT";
                case CMode::GRAT: return "GRAT";
                case CMode::LRAT: return "LRAT";
                case CMode::RESV: return "RESV";
                case CMode::THP:  return "THP";
                case CMode::BHP:  return "BHP";
                case CMode::CRAT: return "CRAT";
                case CMode::GRUP: return "GRUP";
                default: return "none";
                }
            };
            tmp_names[1] = "PROD"; // WellProdDataType::CTRLMode
            tmp_names[2] = fctl(controls.cmode); // WellProdDataType::CTRLMode
        }

        tmp_values[0] = well.getHeadI() + 1; // WellCumDataType::wellLocationi
        tmp_values[1] = well.getHeadJ() + 1; // WellCumDataType::wellLocationj
        tmp_values[2] = get("WOPT"); // WellCumDataType::OilProd
        tmp_values[3] = get("WWPT"); // WellCumDataType::WaterProd
        tmp_values[4] = get("WGPT"); // WellCumDataType::GasProd
        tmp_values[5] = get("WVPT"); // WellCumDataType::FluidResVolProd
        tmp_values[6] = get("WOIT"); // WellCumDataType::OilInj
        tmp_values[7] = get("WWIT"); // WellCumDataType::WaterInj
        tmp_values[8] = get("WGIT"); // WellCumDataType::GasInj
        tmp_values[9] = get("WVIT"); // WellCumDataType::FluidResVolInj

        this->outputCumulativeReport_(tmp_values, tmp_names);
    }
}

template<class Scalar>
void LogOutputHelper<Scalar>::
production(const std::size_t reportStepNum,
           std::function<bool(const std::string&)> isDefunct) const
{
    std::vector<Scalar> tmp_values(WellProdDataType::numWPValues, 0.0);
    std::vector<std::string> tmp_names(WellProdDataType::numWPNames, "");
    this->outputProductionReport_(tmp_values, tmp_names);

    const auto& st = summaryState_;
    for (const auto& gname : schedule_.groupNames()) {
        auto gName = static_cast<std::string>(gname);
        auto get = [&st, &gName](const std::string& vector)
        {
            const auto key = vector + ':' + gName;
            return st.has(key) ? st.get(key) : 0.0;
        };

        tmp_names[0] = gname;
        if (tmp_names[0] == "FIELD") {
            tmp_values[2] = st.get("FOPR", 0.0); // WellProdDataType::OilRate
            tmp_values[3] = st.get("FWPR", 0.0); // WellProdDataType::WaterRate
            tmp_values[4] = st.get("FGPR", 0.0); // WellProdDataType::GasRate
            tmp_values[5] = st.get("FVPR", 0.0); // WellProdDataType::FluidResVol
            tmp_values[6] = st.get("FWCT", 0.0); // WellProdDataType::WaterCut
            tmp_values[7] = st.get("FGOR", 0.0); // WellProdDataType::GasOilRatio
        } else {
            tmp_values[2] = get("GOPR"); // WellProdDataType::OilRate
            tmp_values[3] = get("GWPR"); // WellProdDataType::WaterRate
            tmp_values[4] = get("GGPR"); // WellProdDataType::GasRate
            tmp_values[5] = get("GVPR"); // WellProdDataType::FluidResVol
            tmp_values[6] = get("GWCT"); // WellProdDataType::WaterCut
            tmp_values[7] = get("GGOR"); // WellProdDataType::GasOilRatio
        }

        tmp_values[8] = tmp_values[3] / tmp_values[4]; // WellProdDataType::WaterGasRatio
        if (std::isnan(tmp_values[8])) {
            tmp_values[8] = 0.0;
        }

        this->outputProductionReport_(tmp_values, tmp_names);
    }

    for (const auto& wname : schedule_.wellNames(reportStepNum)) {

        // don't bother with wells not on this process
        if (isDefunct(wname)) {
            continue;
        }

        const auto& well = schedule_.getWell(wname, reportStepNum);

        // Ignore injector wells
        if (well.isInjector()) {
            continue;
        }

        tmp_names[0] = wname; // WellProdDataType::WellName


        auto wName = static_cast<std::string>(wname);
        auto get = [&st, &wName](const std::string& vector)
        {
            const auto key = vector + ':' + wName;
            return st.has(key) ? st.get(key) : 0.0;
        };

        const auto& controls = well.productionControls(st);
        using CMode = Well::ProducerCMode;

        auto fctl = [](const auto wmctl) -> std::string
        {
            switch (wmctl) {
            case CMode::ORAT: return "ORAT";
            case CMode::WRAT: return "WRAT";
            case CMode::GRAT: return "GRAT";
            case CMode::LRAT: return "LRAT";
            case CMode::RESV: return "RESV";
            case CMode::THP:  return "THP";
            case CMode::BHP:  return "BHP";
            case CMode::CRAT: return "CRate";
            case CMode::GRUP: return "GRUP";
            default: return "none";
            }
        };

        tmp_names[1] = fctl(controls.cmode); // WellProdDataType::CTRLMode

        tmp_values[0] = well.getHeadI() + 1; // WellProdDataType::WellLocationi
        tmp_values[1] = well.getHeadJ() + 1; // WellProdDataType::WellLocationj
        tmp_values[2] = get("WOPR"); // WellProdDataType::OilRate
        tmp_values[3] = get("WWPR"); // WellProdDataType::WaterRate
        tmp_values[4] = get("WGPR"); // WellProdDataType::GasRate
        tmp_values[5] = get("WVPR"); // WellProdDataType::FluidResVol
        tmp_values[6] = get("WWCT"); // WellProdDataType::WaterCut
        tmp_values[7] = get("WGOR"); // WellProdDataType::GasOilRatio
        tmp_values[9] = get("WBHP"); // WellProdDataType::BHP
        tmp_values[10] = get("WTHP"); // WellProdDataType::THP
        //tmp_values[11] = 0; //WellProdDataType::SteadyStatePI //

        tmp_values[8] = tmp_values[3] / tmp_values[4]; // WellProdDataType::WaterGasRatio
        if (std::isnan(tmp_values[8])) {
            tmp_values[8] = 0.0;
        }

        this->outputProductionReport_(tmp_values, tmp_names);
    }
}

template<class Scalar>
void LogOutputHelper<Scalar>::
outputCumulativeReport_(const std::vector<Scalar>& wellCum,
                        const std::vector<std::string>& wellCumNames) const
{
    const UnitSystem& units = eclState_.getUnits();
    std::ostringstream ss;
    if (wellCumNames[WellCumDataType::WellName].empty()) {
        ss << "=================================================== CUMULATIVE PRODUCTION/INJECTION REPORT =========================================\n"
           << ":  WELL  :  LOCATION :  WELL  :CTRL:    OIL    :   WATER   :    GAS    :   Prod    :    OIL    :   WATER   :    GAS    :   INJ     :\n"
           << ":  NAME  :  (I,J,K)  :  TYPE  :MODE:    PROD   :   PROD    :    PROD   :  RES.VOL. :    INJ    :   INJ     :    INJ    :  RES.VOL. :\n";
        if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_METRIC) {
            ss << ":        :           :        :    :    MSCM   :   MSCM    :    MMSCM  :   MRCM    :    MSCM   :   MSCM    :    MMSCM  :   MRCM    :\n";
        } else if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_FIELD) {
            ss << ":        :           :        :    :    MSTB   :   MSTB    :    MMSCF  :   MRB     :    MSTB   :   MSTB    :    MMSCF  :   MRB     :\n";
        } else if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_LAB) {
            ss << ":        :           :        :    :     MSCC  :   MSCC    :    MMSCC  :   MRCC    :    MSCC   :   MSCC    :    MMSCC  :   MRCC    :\n";
        }
        ss << "====================================================================================================================================\n";
    } else {
        if (wellCum[WellCumDataType::WellLocationi] < 1) {
            ss << std::right << std::fixed << std::setprecision(0) << ":" << std::setw (8)
               << wellCumNames[WellCumDataType::WellName] << ":"
               << std::setw(11) <<  "" << ":"
               << std::setw(8) << wellCumNames[WellCumDataType::WellType] << ":"
               << std::setw(4) << wellCumNames[WellCumDataType::WellCTRL] << ":"
               << std::setprecision(1) << std::setw(11) << wellCum[WellCumDataType::OilProd] / 1000.0 << ":"
               << std::setw(11) << wellCum[WellCumDataType::WaterProd] / 1000.0 << ":"
               << std::setw(11)<< wellCum[WellCumDataType::GasProd] / 1000.0 << ":"
               << std::setw(11) << wellCum[WellCumDataType::FluidResVolProd] / 1000.0 << ":"
               << std::setw(11) << wellCum[WellCumDataType::OilInj] / 1000.0 << ":"
               << std::setw(11) << wellCum[WellCumDataType::WaterInj] / 1000.0 << ":"
               << std::setw(11) << wellCum[WellCumDataType::GasInj] / 1000.0 << ":"
               << std::setw(11) << wellCum[WellCumDataType::FluidResVolInj] / 1000.0 << ": \n";
        } else {
            ss << std::right << std::fixed << std::setprecision(0) << ":"
               << std::setw (8) << wellCumNames[WellCumDataType::WellName] << ":"
               << std::setw(5) << wellCum[WellCumDataType::WellLocationi] << ","
               << std::setw(5) << wellCum[WellCumDataType::WellLocationj] << ":"
               << std::setw(8) << wellCumNames[WellCumDataType::WellType] << ":"
               << std::setw(4) << wellCumNames[WellCumDataType::WellCTRL] << ":"
               << std::setprecision(1) << std::setw(11) << wellCum[WellCumDataType::OilProd] / 1000.0 << ":"
               << std::setw(11) << wellCum[WellCumDataType::WaterProd] / 1000.0 << ":"
               << std::setw(11)<< wellCum[WellCumDataType::GasProd] / 1000.0 << ":"
               << std::setw(11) << wellCum[WellCumDataType::FluidResVolProd] / 1000.0 << ":"
               << std::setw(11) << wellCum[WellCumDataType::OilInj] / 1000.0 << ":"
               << std::setw(11) << wellCum[WellCumDataType::WaterInj] / 1000.0 << ":"
               << std::setw(11) << wellCum[WellCumDataType::GasInj] / 1000.0 << ":"
               << std::setw(11) << wellCum[WellCumDataType::FluidResVolInj] / 1000.0 << ": \n";
        }
        ss << ":--------:-----------:--------:----:------------:----------:-----------:-----------:------------:----------:-----------:-----------: \n";
    }
    OpmLog::note(ss.str());
}

template<class Scalar>
void LogOutputHelper<Scalar>::
outputProductionReport_(const std::vector<Scalar>& wellProd,
                        const std::vector<std::string>& wellProdNames) const
{
    const UnitSystem& units = eclState_.getUnits();
    std::ostringstream ss;
    if (wellProdNames[WellProdDataType::WellName].empty()) {
        ss << "======================================================= PRODUCTION REPORT =======================================================\n"//=================== \n"
           << ":  WELL  :  LOCATION :CTRL:    OIL    :   WATER   :    GAS    :   FLUID   :   WATER   : GAS/OIL  :  WAT/GAS   : BHP OR : THP OR :\n"// STEADY-ST PI       :\n"
           << ":  NAME  :  (I,J,K)  :MODE:    RATE   :   RATE    :    RATE   :  RES.VOL. :    CUT    :  RATIO   :   RATIO    : CON.PR.: BLK.PR.:\n";// OR POTN OF PREF. PH:\n";
        if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_METRIC) {
            ss << ":        :           :    :  SCM/DAY  :  SCM/DAY  :  SCM/DAY  :  RCM/DAY  :  SCM/SCM  :  SCM/SCM :  SCM/SCM   :  BARSA :  BARSA :\n";//                    :\n";
        } else if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_FIELD) {
            ss << ":        :           :    :  STB/DAY  :  STB/DAY  :  MSCF/DAY :  RB/DAY   :           : MSCF/STB :  STB/MSCF  :  PSIA  :  PSIA  :\n";//                    :\n";
        } else if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_LAB) {
            ss << ":        :           :    :  SCC/HR   :  SCC/HR   :  SCC/HR   :    RCC    :  SCC/SCC  :  SCC/SCC :  SCC/SCC   :  ATMA  :  ATMA  :\n";//                    :\n";
        }
        ss << "=================================================================================================================================\n";//=================== \n";
    } else {
        if (wellProd[WellProdDataType::WellLocationi] < 1) {
            ss << std::right << std::fixed << ":"
               << std::setw(8) << wellProdNames[WellProdDataType::WellName] << ":"
               << std::setprecision(0) << std::setw(11) << "" << ":"
               << std::setw(4) << wellProdNames[WellProdDataType::CTRLMode] << ":"
               << std::setprecision(1) << std::setw(11) << wellProd[WellProdDataType::OilRate] << ":"
               << std::setw(11) << wellProd[WellProdDataType::WaterRate] << ":"
               << std::setw(11)<< wellProd[WellProdDataType::GasRate] << ":"
               << std::setw(11) << wellProd[WellProdDataType::FluidResVol] << std::setprecision(3) << ":"
               << std::setw(11) << wellProd[WellProdDataType::WaterCut] << std::setprecision(2) << ":"
               << std::setw(10) << wellProd[WellProdDataType::GasOilRatio] << std::setprecision(4) << ":"
               << std::setw(12) << wellProd[WellProdDataType::WatGasRatio] << std::setprecision(1) << ":"
               << std::setw(8) << "" << ":" <<  std::setw(8) << "" << ": \n";
        } else {
            ss << std::right << std::fixed << ":"
               << std::setw (8) << wellProdNames[WellProdDataType::WellName] << ":"
               << std::setprecision(0) << std::setw(5) << wellProd[WellProdDataType::WellLocationi] << ","
               << std::setw(5) << wellProd[WellProdDataType::WellLocationj] << ":"
               << std::setw(4) << wellProdNames[WellProdDataType::CTRLMode] << ":"
               << std::setprecision(1) << std::setw(11) << wellProd[WellProdDataType::OilRate] << ":"
               << std::setw(11) << wellProd[WellProdDataType::WaterRate] << ":"
               << std::setw(11)<< wellProd[WellProdDataType::GasRate] << ":"
               << std::setw(11) << wellProd[WellProdDataType::FluidResVol] << std::setprecision(3) << ":"
               << std::setw(11) << wellProd[WellProdDataType::WaterCut] << std::setprecision(2) << ":"
               << std::setw(10) << wellProd[WellProdDataType::GasOilRatio] << std::setprecision(4) << ":"
               << std::setw(12) << wellProd[WellProdDataType::WatGasRatio] << std::setprecision(1) << ":"
               << std::setw(8) << wellProd[WellProdDataType::BHP] << ":"
               << std::setw(8) << wellProd[WellProdDataType::THP] << ": \n";
        }
        ss << ":"<< std::setfill ('-') << std::setw (9) << ":"
           << std::setfill ('-') << std::setw (12) << ":"
           << std::setfill ('-') << std::setw (5) << ":"
           << std::setfill ('-') << std::setw (12) << ":"
           << std::setfill ('-') << std::setw (12) << ":"
           << std::setfill ('-') << std::setw (12) << ":"
           << std::setfill ('-') << std::setw (12) << ":"
           << std::setfill ('-') << std::setw (12) << ":"
           << std::setfill ('-') << std::setw (11) << ":"
           << std::setfill ('-') << std::setw (13) << ":"
           << std::setfill ('-') << std::setw (9) << ":"
           << std::setfill ('-') << std::setw (9) << ":" << "\n";
    }
    OpmLog::note(ss.str());
}

template class LogOutputHelper<double>;

} // namespace Opm
