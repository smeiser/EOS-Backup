/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2019 Nico Gubernari
 * Copyright (c) 2021 Méril Reboud
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

#include <test/test.hh>
#include <eos/rare-b-decays/nonlocal-formfactors.hh>

using namespace test;
using namespace eos;

class NonlocalFormFactorGvDV2020Test :
    public TestCase
{
    public:
        NonlocalFormFactorGvDV2020Test() :
            TestCase("nonlocal_formfactor_GvDV2020_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["mass::B_d"]                               = 5.279;
                p["mass::K_d^*"]                             = 0.896;
                p["mass::J/psi"]                             = 3.0969;
                p["mass::psi(2S)"]                           = 3.6860;
                p["mass::D^0"]                               = 1.86723;
                p["b->sccbar::t_0"]                          = 4.0;
                p["b->sccbar::t_s"]                          = -17.4724;
                p["b->sccbar::chiOPE@GvDV2020"]              = 1.81e-4;
                p["B->K^*ccbar::Re{alpha_0^perp}@GvDV2020"]  = 2.0;
                p["B->K^*ccbar::Im{alpha_0^perp}@GvDV2020"]  = 3.0;
                p["B->K^*ccbar::Re{alpha_1^perp}@GvDV2020"]  = 4.0;
                p["B->K^*ccbar::Im{alpha_1^perp}@GvDV2020"]  = 5.0;
                p["B->K^*ccbar::Re{alpha_2^perp}@GvDV2020"]  = 6.0;
                p["B->K^*ccbar::Im{alpha_2^perp}@GvDV2020"]  = 7.0;
                p["B->K^*ccbar::Re{alpha_0^para}@GvDV2020"]  = 8.0;
                p["B->K^*ccbar::Im{alpha_0^para}@GvDV2020"]  = 9.0;
                p["B->K^*ccbar::Re{alpha_1^para}@GvDV2020"]  = 10.0;
                p["B->K^*ccbar::Im{alpha_1^para}@GvDV2020"]  = 11.0;
                p["B->K^*ccbar::Re{alpha_2^para}@GvDV2020"]  = 12.0;
                p["B->K^*ccbar::Im{alpha_2^para}@GvDV2020"]  = 13.0;
                p["B->K^*ccbar::Re{alpha_0^long}@GvDV2020"]  = 14.0;
                p["B->K^*ccbar::Im{alpha_0^long}@GvDV2020"]  = 15.0;
                p["B->K^*ccbar::Re{alpha_1^long}@GvDV2020"]  = 16.0;
                p["B->K^*ccbar::Im{alpha_1^long}@GvDV2020"]  = 17.0;
                p["B->K^*ccbar::Re{alpha_2^long}@GvDV2020"]  = 18.0;
                p["B->K^*ccbar::Im{alpha_2^long}@GvDV2020"]  = 19.0;

                Options o = {
                    { "model", "WET" },
                    { "q", "d" }
                };

                auto nff = NonlocalFormFactor<nff::PToV>::make("B->K^*::GvDV2020", p, o);


                auto diagnostics = nff->diagnostics();

                std::cout << "Diagnostics:" << std::endl;
                for (auto & d : diagnostics)
                {
                    std::cout << d.description << ": " << d.value << std::endl;
                }
                std::cout << "Diagnostics ended" << std::endl;

                static const std::vector<std::pair<double, double>> reference
                {
                    /* outer functions */
                    std::make_pair( 0.0,       eps),        // Re{1/phi_long(q2 = 0.0)}
                    std::make_pair( 0.0,       eps),        // Im{1/phi_long(q2 = 0.0)}

                    std::make_pair(-39.01168,  eps),        // Re{phi_long(q2 = 16.0)}
                    std::make_pair( 10.87513,  eps),        // Im{phi_long(q2 = 16.0)}

                    std::make_pair( 24.64525,  eps),        // Re{phi_perp(q2 = 16.0)}
                    std::make_pair(-18.28392,  eps)         // Im{phi_perp(q2 = 16.0)}
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                TEST_CHECK_RELATIVE_ERROR(real(nff->H_perp(16.0)), -0.761178,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_perp(16.0)), -1.025598,   eps);
                TEST_CHECK_RELATIVE_ERROR(real(nff->H_para(16.0)), -1.467510,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_para(16.0)), -1.877480,   eps);
                TEST_CHECK_RELATIVE_ERROR(real(nff->H_long(16.0)),  2.278766,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_long(16.0)),  1.340658,   eps);

                TEST_CHECK_RELATIVE_ERROR(real(nff->H_perp_residue_jpsi()),  -30.9653,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_perp_residue_jpsi()),  -36.3466,   eps);
                TEST_CHECK_RELATIVE_ERROR(real(nff->H_perp_residue_psi2s()),   3.4865,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_perp_residue_psi2s()),   4.05179,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(nff->H_para_residue_jpsi()),  -63.2528,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_para_residue_jpsi()),  -68.634,    eps);
                TEST_CHECK_RELATIVE_ERROR(real(nff->H_para_residue_psi2s()),   6.87827,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_para_residue_psi2s()),   7.44356,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(nff->H_long_residue_jpsi()),   29.8167,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_long_residue_jpsi()),   31.4962,   eps);
                TEST_CHECK_RELATIVE_ERROR(real(nff->H_long_residue_psi2s()),  -6.0986,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_long_residue_psi2s()),  -6.43428,  eps);

            }
        }
} nonlocal_formfactor_gvdv2020_test;
