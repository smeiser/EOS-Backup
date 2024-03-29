/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Méril Reboud
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef EOS_GUARD_EOS_RARE_B_DECAYS_LAMBDA_B_TO_LAMBDA_NU_NU_HH
#define EOS_GUARD_EOS_RARE_B_DECAYS_LAMBDA_B_TO_LAMBDA_NU_NU_HH 1

#include <eos/observable.hh>

#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/reference-name.hh>

namespace eos
{
    /*
     * Decay: L_b -> L nu nu
     */
    class LambdaBToLambdaDineutrino :
        public ParameterUser,
        public PrivateImplementationPattern<LambdaBToLambdaDineutrino>
    {
        public:
            LambdaBToLambdaDineutrino(const Parameters & parameters, const Options & options);
            ~LambdaBToLambdaDineutrino();

            struct AngularCoefficients;

            // Single-differential Observables
            double differential_decay_width(const double & q2) const;
            double differential_branching_ratio(const double & q2) const;
            double differential_longitudinal_polarisation(const double & q2) const;

            // Integrated Observables
            class IntermediateResult;
            const IntermediateResult * prepare(const double & q2_min, const double & q2_max) const;
            double integrated_decay_width(const IntermediateResult * ir) const;
            double integrated_branching_ratio(const IntermediateResult * ir) const;
            double integrated_longitudinal_polarisation(const IntermediateResult * ir) const;

            /*!
             * Descriptions of the process and its kinematics.
             */
            static const std::string description;
            static const std::string kinematics_description_q2;

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;

            /*!
             * Options used in the computation of our observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
    };
}

#endif
