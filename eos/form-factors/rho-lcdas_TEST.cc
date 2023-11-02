/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Stefan Meiser
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
#include <eos/form-factors/rho-lcdas.hh>

#include <eos/models/model.hh>

#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;

class RhoLCDAsTest :
    public TestCase
{
    public:
        RhoLCDAsTest() :
            TestCase("rho_lcdas_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            Parameters p = Parameters::Defaults();
            p["QCD::alpha_s(MZ)"] = 0.1176;
            p["mass::d(2GeV)"]    = 0.0048;
            p["mass::u(2GeV)"]    = 0.0032;
            p["rho::a2para@1GeV"] = 0.22;
            p["rho::a4para@1GeV"] = 0.16;
            p["rho::a2perp@1GeV"] = 0.14;
            p["rho::a4perp@1GeV"] = 0.25;
            p["rho::fperp@1GeV"]  = 0.16;

            /* Diagnostics */
            {
                RhoLCDAs rho(p, Options{ });
                Diagnostics diagnostics = rho.diagnostics();
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(+1.00000, 1e-5), // c_rge(mu = 1.0 GeV)
                    std::make_pair(+0.94850, 1e-5), // c_rge(mu = 2.0 GeV)
                    std::make_pair(+0.92874, 1e-5), // c_rge(mu = 3.0 GeV)
                    std::make_pair(+0.91708, 1e-5), // c_rge(mu = 4.0 GeV)
                    std::make_pair(+0.90893, 1e-5), // c_rge(mu = 5.0 GeV)
                };

                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);
            }

            /* Twist 2 */
            {
                RhoLCDAs rho(p, Options{ });

                // coefficients at mu = 1.0 GeV, 2.0 GeV, 3.0 GeV, 4.0 GeV and 5.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.0,       rho.a1para(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,       rho.a1para(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,       rho.a1para(3.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,       rho.a1para(4.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,       rho.a1para(5.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.22,      rho.a2para(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.164004,  rho.a2para(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.145901,  rho.a2para(3.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.136008,  rho.a2para(4.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.129431,  rho.a2para(5.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,       rho.a3para(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,       rho.a3para(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,       rho.a3para(3.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,       rho.a3para(4.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,       rho.a3para(5.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.16,      rho.a4para(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.1043234, rho.a4para(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0879874, rho.a4para(3.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0794371, rho.a4para(4.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0739064, rho.a4para(5.0),   eps);

                TEST_CHECK_NEARLY_EQUAL( 0.0,       rho.a1perp(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,       rho.a1perp(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,       rho.a1perp(3.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,       rho.a1perp(4.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,       rho.a1perp(5.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.14,      rho.a2perp(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.103147,  rho.a2perp(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0913332, rho.a2perp(3.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0849017, rho.a2perp(4.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0806361, rho.a2perp(5.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,       rho.a3perp(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,       rho.a3perp(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,       rho.a3perp(3.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,       rho.a3perp(4.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,       rho.a3perp(5.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.25,      rho.a4perp(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.162241,  rho.a4perp(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.13658,   rho.a4perp(3.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.12317,   rho.a4perp(4.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.11450,   rho.a4perp(5.0),   eps);

                TEST_CHECK_NEARLY_EQUAL( 0.16,      rho.fperp(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.147292,  rho.fperp(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.142517,  rho.fperp(3.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.139726,  rho.fperp(4.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.137788,  rho.fperp(5.0),   eps);

                // scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.0,      rho.phipara(0.0, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.911333, rho.phipara(0.3, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.455,    rho.phipara(0.5, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.911333, rho.phipara(0.7, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,      rho.phipara(1.0, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,      rho.phiperp(0.0, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.792225, rho.phiperp(0.3, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.88813,  rho.phiperp(0.5, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.792225, rho.phiperp(0.7, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,      rho.phiperp(1.0, 1.0), eps);

                // scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.0,      rho.phipara(0.0, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.02489,  rho.phipara(0.3, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.4244,   rho.phipara(0.5, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.02489,  rho.phipara(0.7, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,      rho.phipara(1.0, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,      rho.phiperp(0.0, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.951784, rho.phiperp(0.3, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.72422,  rho.phiperp(0.5, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.951784, rho.phiperp(0.7, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,      rho.phiperp(1.0, 2.0), eps);
            }

        }
} rho_lcdas_test;
