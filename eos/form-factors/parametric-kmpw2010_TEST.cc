/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2013-2016, 2018 Danny van Dyk
 * Copyright (c) 2015 Christoph Bobeth
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
#include <eos/form-factors/parametric-kmpw2010.hh>

using namespace test;
using namespace eos;

class BToKstarKMPW2010FormFactorsTest:
    public TestCase
{
    public:
        BToKstarKMPW2010FormFactorsTest() :
            TestCase("b_to_kstar_kmpw2010_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 5e-6;

            // central values
            {
                Parameters p = Parameters::Defaults();
                p["B->K^*::F^V(0)@KMPW2010"]  = 0.36; p["B->K^*::b^V_1@KMPW2010"]  =  -4.8;
                p["B->K^*::F^A0(0)@KMPW2010"] = 0.29; p["B->K^*::b^A0_1@KMPW2010"] = -18.2;
                p["B->K^*::F^A1(0)@KMPW2010"] = 0.25; p["B->K^*::b^A1_1@KMPW2010"] =  +0.34;
                p["B->K^*::F^A2(0)@KMPW2010"] = 0.23; p["B->K^*::b^A2_1@KMPW2010"] =  -0.85;
                p["B->K^*::F^T1(0)@KMPW2010"] = 0.31; p["B->K^*::b^T1_1@KMPW2010"] =  -4.6;
                p["B->K^*::F^T2(0)@KMPW2010"] = 0.31; p["B->K^*::b^T2_1@KMPW2010"] =  -3.2;
                p["B->K^*::F^T3(0)@KMPW2010"] = 0.22; p["B->K^*::b^T3_1@KMPW2010"] = -10.3;
                std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->K^*::KMPW2010", p, Options{ });
                TEST_CHECK(0 != ff.get());

                // V
                TEST_CHECK_RELATIVE_ERROR(0.421964, ff->v( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.497285, ff->v( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.590298, ff->v( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.707366, ff->v( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.858208, ff->v(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(1.058404, ff->v(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(1.334545, ff->v(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(1.735915, ff->v(18.4), eps);

                // A_0
                TEST_CHECK_RELATIVE_ERROR(0.410875, ff->a_0( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.560285, ff->a_0( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.747894, ff->a_0( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.988041, ff->a_0( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(1.302880, ff->a_0(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(1.728452, ff->a_0(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(2.327350, ff->a_0(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(3.218300, ff->a_0(18.4), eps);

                // A_1
                TEST_CHECK_RELATIVE_ERROR(0.266631, ff->a_1( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.285779, ff->a_1( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.308075, ff->a_1( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.334383, ff->a_1( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.365917, ff->a_1(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.404440, ff->a_1(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.452618, ff->a_1(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(0.514681, ff->a_1(18.4), eps);

                // A_2
                TEST_CHECK_RELATIVE_ERROR(0.250198, ff->a_2( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.273759, ff->a_2( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.301563, ff->a_2( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.334822, ff->a_2( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.375252, ff->a_2(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.425358, ff->a_2(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.488950, ff->a_2(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(0.572108, ff->a_2(18.4), eps);

                // t_1
                TEST_CHECK_RELATIVE_ERROR(0.362235, ff->t_1( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.425698, ff->t_1( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.504029, ff->t_1( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.602573, ff->t_1( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.729487, ff->t_1(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.897853, ff->t_1(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(1.12999 , ff->t_1(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(1.46727 , ff->t_1(18.4), eps);

                // t_2
                TEST_CHECK_RELATIVE_ERROR(0.35026 , ff->t_2( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.397838, ff->t_2( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.45472 , ff->t_2( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.523652, ff->t_2( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.608536, ff->t_2(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.715103, ff->t_2(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.852101, ff->t_2(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(1.03355 , ff->t_2(18.4), eps);

                // t_3
                TEST_CHECK_RELATIVE_ERROR(0.276523, ff->t_3( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.344214, ff->t_3( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.426193, ff->t_3( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.526794, ff->t_3( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.6522  , ff->t_3(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.811519, ff->t_3(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(1.01871 , ff->t_3(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(1.29621 , ff->t_3(18.4), eps);
            }

            // raised values
            {
                Parameters p = Parameters::Defaults();
                p["B->K^*::F^V(0)@KMPW2010"]  = 0.59; p["B->K^*::b^V_1@KMPW2010"]  =  -4.0;
                p["B->K^*::F^A0(0)@KMPW2010"] = 0.39; p["B->K^*::b^A0_1@KMPW2010"] = -16.9;
                p["B->K^*::F^A1(0)@KMPW2010"] = 0.41; p["B->K^*::b^A1_1@KMPW2010"] =  +1.2;
                p["B->K^*::F^A2(0)@KMPW2010"] = 0.42; p["B->K^*::b^A2_1@KMPW2010"] =  +2.03;
                p["B->K^*::F^T1(0)@KMPW2010"] = 0.49; p["B->K^*::b^T1_1@KMPW2010"] =  -4.6 + 0.81;
                p["B->K^*::F^T2(0)@KMPW2010"] = 0.49; p["B->K^*::b^T2_1@KMPW2010"] =  -3.2 + 2.1;
                p["B->K^*::F^T3(0)@KMPW2010"] = 0.39; p["B->K^*::b^T3_1@KMPW2010"] = -10.3 + 2.5;
                std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->K^*::KMPW2010", p, Options{ });
                TEST_CHECK(0 != ff.get());

                // V
                TEST_CHECK_RELATIVE_ERROR(0.683006, ff->v( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.795816, ff->v( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.934831, ff->v( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(1.109443, ff->v( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(1.333989, ff->v(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(1.631445, ff->v(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(2.041006, ff->v(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(2.635294, ff->v(18.4), eps);

                // A_0
                TEST_CHECK_RELATIVE_ERROR(0.543363, ff->a_0( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.732820, ff->a_0( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.970583, ff->a_0( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(1.274770, ff->a_0( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(1.673378, ff->a_0(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(2.211932, ff->a_0(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(2.969510, ff->a_0(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(4.096084, ff->a_0(18.4), eps);

                // A_1
                TEST_CHECK_RELATIVE_ERROR(0.430965, ff->a_1( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.454709, ff->a_1( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.481882, ff->a_1( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.513361, ff->a_1( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.550366, ff->a_1(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.594652, ff->a_1(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.648841, ff->a_1(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(0.717050, ff->a_1(18.4), eps);

                // A_2
                TEST_CHECK_RELATIVE_ERROR(0.435238, ff->a_2( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.451990, ff->a_2( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.470539, ff->a_2( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.491252, ff->a_2( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.514616, ff->a_2(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.541305, ff->a_2(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.572274, ff->a_2(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(0.608954, ff->a_2(18.4), eps);

                // t_1
                TEST_CHECK_RELATIVE_ERROR(0.56538 , ff->t_1( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.656751, ff->t_1( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.769278, ff->t_1( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.910534, ff->t_1( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(1.09208 , ff->t_1(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(1.33244 , ff->t_1(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(1.66321 , ff->t_1(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(2.14292 , ff->t_1(18.4), eps);

                // t_2
                TEST_CHECK_RELATIVE_ERROR(0.535223, ff->t_2( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.588078, ff->t_2( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.650576, ff->t_2( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.725486, ff->t_2( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.816729, ff->t_2(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.930041, ff->t_2(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(1.07414 , ff->t_2(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(1.26297 , ff->t_2(18.4), eps);

                // t_3
                TEST_CHECK_RELATIVE_ERROR(0.472753, ff->t_3( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.571573, ff->t_3( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.690926, ff->t_3( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.837005, ff->t_3( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(1.01864 , ff->t_3(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(1.24883 , ff->t_3(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(1.54748 , ff->t_3(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(1.94657 , ff->t_3(18.4), eps);            }

            // lowered values
            {
                Parameters p = Parameters::Defaults();
                p["B->K^*::F^V(0)@KMPW2010"]  = 0.24; p["B->K^*::b^V_1@KMPW2010"]  =  -5.2;
                p["B->K^*::F^A0(0)@KMPW2010"] = 0.22; p["B->K^*::b^A0_1@KMPW2010"] = -21.2;
                p["B->K^*::F^A1(0)@KMPW2010"] = 0.15; p["B->K^*::b^A1_1@KMPW2010"] =  -0.46;
                p["B->K^*::F^A2(0)@KMPW2010"] = 0.13; p["B->K^*::b^A2_1@KMPW2010"] =  -2.2;
                p["B->K^*::F^T1(0)@KMPW2010"] = 0.21; p["B->K^*::b^T1_1@KMPW2010"] =  -4.6 - 0.41;
                p["B->K^*::F^T2(0)@KMPW2010"] = 0.21; p["B->K^*::b^T2_1@KMPW2010"] =  -3.2 - 2.2;
                p["B->K^*::F^T3(0)@KMPW2010"] = 0.12; p["B->K^*::b^T3_1@KMPW2010"] = -10.3 - 3.1;
                std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->K^*::KMPW2010", p, Options{ });
                TEST_CHECK(0 != ff.get());

                // V
                TEST_CHECK_RELATIVE_ERROR(0.2830474, ff->v( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.3354240, ff->v( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.4001625, ff->v( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.4817170, ff->v( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.5868880, ff->v(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.7265850, ff->v(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.9194250, ff->v(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(1.1999232, ff->v(18.4), eps);

                // A_0
                TEST_CHECK_RELATIVE_ERROR(0.3236650, ff->a_0( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.4519470, ff->a_0( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.6131980, ff->a_0( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.8198090, ff->a_0( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(1.0909360, ff->a_0(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(1.4577380, ff->a_0(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(1.9743500, ff->a_0(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(2.7434500, ff->a_0(18.4), eps);

                // A_1
                TEST_CHECK_RELATIVE_ERROR(0.1621258, ff->a_1( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.1762210, ff->a_1( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.1927954, ff->a_1( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.2125507, ff->a_1( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.2364774, ff->a_1(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.2660210, ff->a_1(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.3033756, ff->a_1(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(0.3520390, ff->a_1(18.4), eps);

                // A_2
                TEST_CHECK_RELATIVE_ERROR(0.1445568, ff->a_2( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.1616854, ff->a_2( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.1820760, ff->a_2( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.2066815, ff->a_2( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.2368550, ff->a_2(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.2745790, ff->a_2(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.3228774, ff->a_2(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(0.3865910, ff->a_2(18.4), eps);

                // t_1
                TEST_CHECK_RELATIVE_ERROR(0.246944, ff->t_1( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.291875, ff->t_1( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.347386, ff->t_1( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.417288, ff->t_1( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.507397, ff->t_1(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.627041, ff->t_1(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.792141, ff->t_1(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(1.03221 , ff->t_1(18.4), eps);

                // t_2
                TEST_CHECK_RELATIVE_ERROR(0.24554 , ff->t_2( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.287805, ff->t_2( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.338645, ff->t_2( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.400627, ff->t_2( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.477404, ff->t_2(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.574348, ff->t_2(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.699678, ff->t_2(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(0.866585, ff->t_2(18.4), eps);

                // t_3
                TEST_CHECK_RELATIVE_ERROR(0.157488, ff->t_3( 2.3), eps);
                TEST_CHECK_RELATIVE_ERROR(0.202489, ff->t_3( 4.6), eps);
                TEST_CHECK_RELATIVE_ERROR(0.257115, ff->t_3( 6.9), eps);
                TEST_CHECK_RELATIVE_ERROR(0.324297, ff->t_3( 9.2), eps);
                TEST_CHECK_RELATIVE_ERROR(0.40822 , ff->t_3(11.5), eps);
                TEST_CHECK_RELATIVE_ERROR(0.515053, ff->t_3(13.8), eps);
                TEST_CHECK_RELATIVE_ERROR(0.654254, ff->t_3(16.1), eps);
                TEST_CHECK_RELATIVE_ERROR(0.841037, ff->t_3(18.4), eps);
            }
        }
} b_to_kstar_kmpw2010_form_factors_test;

class BToKKMPW2010FormFactorsTest :
    public TestCase
{
    public:
        BToKKMPW2010FormFactorsTest() :
            TestCase("b_to_k_kmpw2010_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-10;

            Parameters p = Parameters::Defaults();
            std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("B->K::KMPW2010", p, Options{ });
            TEST_CHECK(0 != ff.get());

            TEST_CHECK_NEARLY_EQUAL(0.3844144474375, ff->f_p( 2.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4381960494587, ff->f_p( 4.6), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5043916755865, ff->f_p( 6.9), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5874824689246, ff->f_p( 9.2), eps);
            TEST_CHECK_NEARLY_EQUAL(0.6943252577904, ff->f_p(11.5), eps);
            TEST_CHECK_NEARLY_EQUAL(0.8359400978914, ff->f_p(13.8), eps);
            TEST_CHECK_NEARLY_EQUAL(1.0311821367700, ff->f_p(16.1), eps);
            TEST_CHECK_NEARLY_EQUAL(1.3150955485260, ff->f_p(18.4), eps);
            TEST_CHECK_NEARLY_EQUAL(1.7607545970800, ff->f_p(20.7), eps);

            TEST_CHECK_NEARLY_EQUAL(0.4417893713152, ff->f_t( 2.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5045402114062, ff->f_t( 4.6), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5818230155226, ff->f_t( 6.9), eps);
            TEST_CHECK_NEARLY_EQUAL(0.6788902047561, ff->f_t( 9.2), eps);
            TEST_CHECK_NEARLY_EQUAL(0.8037805019019, ff->f_t(11.5), eps);
            TEST_CHECK_NEARLY_EQUAL(0.9694161202369, ff->f_t(13.8), eps);
            TEST_CHECK_NEARLY_EQUAL(1.1979110863960, ff->f_t(16.1), eps);
            TEST_CHECK_NEARLY_EQUAL(1.5303740143290, ff->f_t(18.4), eps);
            TEST_CHECK_NEARLY_EQUAL(2.0525389704830, ff->f_t(20.7), eps);

            TEST_CHECK_NEARLY_EQUAL(0.3691335300041, ff->f_0( 2.3), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4001522415666, ff->f_0( 4.6), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4333067770539, ff->f_0( 6.9), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4689025169045, ff->f_0( 9.2), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5073172894003, ff->f_0(11.5), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5490271105427, ff->f_0(13.8), eps);
            TEST_CHECK_NEARLY_EQUAL(0.5946449022903, ff->f_0(16.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.6449812193912, ff->f_0(18.4), eps);
            TEST_CHECK_NEARLY_EQUAL(0.7011445924499, ff->f_0(20.7), eps);
        }
} b_to_k_kmpw2010_form_factors_test;
