/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Christian Wacker
 * Copyright (c) 2022 Viktor Kuschke
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

#include <eos/maths/complex.hh>
#include <eos/maths/power-of.hh>

#include <array>
#include <cmath>
#include <limits>

#include <iostream>

namespace eos
{
    static const int max_iterations = 54;

    namespace dilog_impl
    {
        // for the dilogarithm the coefficients of f1 are : Zeta(2 - m) / m! starting from m = 0 to m = max_iterations
        std::array<complex<double>, max_iterations> series_coefficient_f1
        {{
            +1.6449340668482264365,       0.0,
            -0.25,                       -0.013888888888888888889,
             0.0,                        +0.000069444444444444444444,
             0.0,                        -7.8735197782816830436e-7,
             0.0,                        +1.1482216343327454439e-8,
             0.0,                        -1.8978869988970999072e-10,
             0.0,                        +3.3873013709535212723e-12,
             0.0,                        -6.3726364431831803966e-14,
             0.0,                        +1.2462059912950672305e-15,
             0.0,                        -2.5105444608999545509e-17,
             0.0,                        +5.1782588060906235072e-19,
             0.0,                        -1.0887357368300848844e-20,
             0.0,                        +2.3257441143020872235e-22,
             0.0,                        -5.0351952131473895608e-24,
             0.0,                        +1.1026499294381215333e-25,
             0.0,                        -2.4386585509007344735e-27,
             0.0,                        +5.4401426788562523156e-29,
             0.0,                        -1.2228340131217352117e-30,
             0.0,                        +2.7672634689679505842e-32,
             0.0,                        -6.3000905918320139487e-34,
             0.0,                        +1.4420868388418475211e-35,
             0.0,                        -3.3170939991595428044e-37,
             0.0,                        +7.6639135579206578874e-39,
             0.0,                        -1.7778714733830657873e-40,
             0.0,                        +4.1396058982341373449e-42,
             0.0,                        -9.6715570360811017926e-44,
             0.0,                        +2.2667187016766123705e-45
        }};

        // used for |z| < 0.5 : series expansion of the polylogarithm for n = 2 around the origin
        complex<double> f0(const complex<double> & z)
        {
            static const double eps = std::numeric_limits<double>::epsilon();

            complex<double> result(0.0, 0.0);
            complex<double> x = 1.0;

            for (int i = 1 ; i < max_iterations ; ++i)
            {
                x *= z;

                complex<double> summand = x / power_of<2>(double(i));

                result += summand;

                if (abs(summand / result) < eps)
                    break;
            }

            return result;
        }

        // used for |z| >= 0.5 and |z| <= 2.0
        complex<double> f1(const complex<double> & z)
        {
            complex<double> result(0.0, 0.0);

            const complex<double> lnz = std::log(z);
            complex<double> lnlnz = std::log(-lnz);
            complex<double> x(1.0, 0.0);

            // + 0.0 and - 0.0 are processed differently by the std::log():
            // if x > 0 std::log(- x + 0.0 i) = std::log(x) + i*Pi, that is the main branch of log,
            // whereas std::log(- x - 0.0 i) = std::log(x) - i*Pi.
            // If z > 1 real, lnz > 0 real and -lnz has imaginary part (- 0.0 i),
            // therefore the imaginary part of lnlnz has to be changed from (- i*Pi) to (i*Pi) by complex conjugation.
            if ((lnz.imag() == 0.0) && (lnz.real() > 0.0))
                lnlnz = std::conj(lnlnz);

            for (int i = 0 ; i < max_iterations ; ++i)
            {
                result += series_coefficient_f1[i] * x;
                x *= lnz;
            }
            result += lnz * (1.0 - lnlnz);

            return result;
        }

        // used for |z| > 2.0
        complex<double> g(const complex<double> & z)
        {
            const complex<double> lnz = std::log(z);
            const complex<double> a = -power_of<2>(complex<double>(0.0, 2.0 * M_PI)) / 2.0;
            const complex<double> b = complex<double>(0.0, -0.5 / M_PI) * lnz;
            double theta = 0.0;

            if ((z.imag() < 0.0) || ((z.real() >= 1.0) && (z.imag() == 0.0)))
            {
                theta = 1.0;
            }

            complex<double> result = a * (b * b - b + 1.0 / 6.0);
            result += complex<double>(0.0, -2.0 * M_PI) * theta * lnz;

            return result;
        }
    }

    // Calculation of the dilogarithm based on [C:2006A]
    complex<double> dilog(const complex<double> & z)
    {
        // special cases
        if (z == complex<double>(0.0, 0.0))
            return complex<double>(0.0);

        if (z == complex<double>(1.0, 0.0))
            return M_PI * M_PI / 6.0;

        if (z == complex<double>(-1.0, 0.0))
            return -M_PI * M_PI / 12.0;

        if (std::abs(z) < 0.5)
            return dilog_impl::f0(z); // series expansion around the origin

        if (std::abs(z) > 2.0)
            return dilog_impl::g(z) - dilog_impl::f0(1.0 / z); // reflection formula

        return dilog_impl::f1(z);
    }

    namespace trilog_impl
    {
        // for the trilogarithm the coefficients of f1 are : Zeta(3 - m) / m! starting from m = 0
        std::array<complex<double>, max_iterations> series_coefficient_f1
        {{
            +1.2020569031595943,         +1.6449340668482264,
             0.0,                        -0.083333333333333333,
            -0.0034722222222222222,       0.0,
            +0.000011574074074074074,     0.0,
            -9.8418997228521038e-8,       0.0,
            +1.1482216343327454e-9,       0.0,
            -1.5815724990809166e-11,      0.0,
            +2.4195009792525152e-13,      0.0,
            -3.9828977769894877e-15,      0.0,
            +6.9233666183059291e-17,      0.0,
            -1.2552722304499773e-18,      0.0,
            +2.3537540027684652e-20,      0.0,
            -4.536398903458687e-22,       0.0,
            +8.9451696703926432e-24,      0.0,
            -1.7982840046954963e-25,      0.0,
            +3.6754997647937384e-27,      0.0,
            -7.6208079715647952e-29,      0.0,
            +1.600041964369486e-30,       0.0,
            -3.3967611475603756e-32,      0.0,
            +7.2822722867577647e-34,      0.0,
            -1.5750226479580035e-35,      0.0,
            +3.4335400924805893e-37,      0.0,
            -7.53884999808987e-39,        0.0,
            +1.6660681647653604e-40,      0.0,
            -3.7038989028813871e-42,      0.0,
            +8.2792117964682747e-44,      0.0,
            -1.8599148146309811e-45,      0.0
        }};

        // used for |z| < 0.5 : series expansion of the polylogarithm for n = 3 around the origin
        complex<double> f0(const complex<double> & z)
        {
            static const double eps = std::numeric_limits<double>::epsilon();

            complex<double> result(0.0, 0.0);
            complex<double> x = 1.0;

            for (int i = 1 ; i < max_iterations ; ++i)
            {
                x *= z;

                complex<double> summand = x / power_of<3>(double(i));

                result += summand;

                if (abs(summand / result) < eps)
                    break;
            }

            return result;
        }

        // used for |z| >= 0.5 and |z| <= 2.0
        complex<double> f1(const complex<double> & z)
        {
            complex<double> result(0.0, 0.0);

            const complex<double> lnz = std::log(z);
            complex<double> lnlnz = std::log(-lnz);
            complex<double> x(1.0, 0.0);

            // + 0.0 and - 0.0 are processed differently by the std::log():
            // if x > 0 std::log(- x + 0.0 i) = std::log(x) + i*Pi, that is the main branch of log,
            // whereas std::log(- x - 0.0 i) = std::log(x) - i*Pi.
            // If z > 1 real, lnz > 0 real and -lnz has imaginary part (- 0.0 i),
            // therefore the imaginary part of lnlnz has to be changed from (- i*Pi) to (i*Pi) by complex conjugation.
            if ((lnz.imag() == 0.0) && (lnz.real() > 0.0))
                lnlnz = std::conj(lnlnz);

            for (int i = 0 ; i < max_iterations ; ++i)
            {
                result += series_coefficient_f1[i] * x;
                x *= lnz;
            }
            result += 0.5 * lnz * lnz * (3.0 / 2.0 - lnlnz);

            return result;
        }

        // used for |z| > 2.0
        complex<double> g(const complex<double> & z)
        {
            const complex<double> lnz = std::log(z);
            const complex<double> a = -power_of<3>(complex<double>(0.0, 2.0 * M_PI)) / 6.0;
            const complex<double> b = complex<double>(0.0, -0.5 / M_PI) * lnz;
            double theta = 0.0;

            if ((z.imag() < 0.0) || ((z.real() >= 1.0) && (z.imag() == 0.0)))
            {
                theta = 1.0;
            }

            complex<double> result = a * (b * b * b - 3.0 / 2.0 * b * b + b / 2.0);
            result += complex<double>(0.0, -M_PI) * theta * power_of<2>(lnz);

            return result;
        }
    }

    // Calculation of the trilogarithm based on [C:2006A]
    complex<double> trilog(const complex<double> & z)
    {
        static const double aperys_constant = 1.2020569031595942854;

        // special cases
        if (z == complex<double>(0.0, 0.0))
            return complex<double>(0.0);

        if (z == complex<double>(1.0, 0.0))
            return aperys_constant;

        if (z == complex<double>(-1.0, 0.0))
            return -3.0 / 4.0 * aperys_constant;

        if (std::abs(z) < 0.5)
            return trilog_impl::f0(z); // series expansion around the origin

        if (std::abs(z) > 2.0)
            return trilog_impl::g(z) + trilog_impl::f0(1.0 / z); // reflection formula

        return trilog_impl::f1(z);
    }

    namespace quadlog_impl
    {
        // for the quadlogarithm the coefficients of f1 are : Zeta(4 - m) / m! starting from m = 0
        std::array<complex<double>, max_iterations> series_coefficient_f1
        {{
            +1.082323233711138192,        +1.20205690315959428,
            +0.822467033424113218,         0.0,
            -0.020833333333333333,        -0.000694444444444444,
            0.0,                          +1.65343915343915344e-6,
            0.0,                          -1.09354441365023376e-8,
            0.0,                          +1.04383784939340495e-10,
            0.0,                          -1.21659423006224353e-12,
            0.0,                          +1.61300065283501013e-14,
            0.0,                          -2.34288104528793397e-16,
            0.0,                          +3.64387716752943635e-18,
            0.0,                          -5.97748681166655845e-20,
            0.0,                          +1.02337130555150662e-21,
            0.0,                          -1.81455956138347481e-23,
            0.0,                          +3.31302580384912710e-25,
            0.0,                          -6.20097932653619404e-27,
            0.0,                          +1.18564508541733498e-28,
            0.0,                          -2.30933574895902886e-30,
            0.0,                          +4.57154846962710279e-32,
            0.0,                          -9.18043553394696105e-34,
            0.0,                          +1.86724930429686274e-35,
            0.0,                          -3.84151865355610607e-37,
            0.0,                          +7.98497695925718450e-39,
            0.0,                          -1.67529999957552667e-40,
            0.0,                          +3.54482588247949023e-42,
            0.0,                          -7.55897735281915726e-44,
            0.0,                          +1.62337486205260288e-45,
            0.0,                          -3.50927323515279455e-47
        }};

        // used for |z| < 0.5 : series expansion of the polylogarithm for n = 4 around the origin
        complex<double> f0(const complex<double> & z)
        {
            static const double eps = std::numeric_limits<double>::epsilon();

            complex<double> result(0.0, 0.0);
            complex<double> x = 1.0;

            for (int i = 1 ; i < max_iterations ; ++i)
            {
                x *= z;

                complex<double> summand = x / power_of<4>(double(i));

                result += summand;

                if (abs(summand / result) < eps)
                    break;
            }

            return result;
        }

        // used for |z| >= 0.5 and |z| <= 2.0
        complex<double> f1(const complex<double> & z)
        {
            complex<double> result(0.0, 0.0);

            const complex<double> lnz = std::log(z);
            complex<double> lnlnz = std::log(-lnz);
            complex<double> x(1.0, 0.0);

            // + 0.0 and - 0.0 are processed differently by the std::log():
            // if x > 0 std::log(- x + 0.0 i) = std::log(x) + i*Pi, that is the main branch of log,
            // whereas std::log(- x - 0.0 i) = std::log(x) - i*Pi.
            // If z > 1 real, lnz > 0 real and -lnz has imaginary part (- 0.0 i),
            // therefore the imaginary part of lnlnz has to be changed from (- i*Pi) to (i*Pi) by complex conjugation.
            if ((lnz.imag() == 0.0) && (lnz.real() > 0.0))
                lnlnz = std::conj(lnlnz);

            for (int i = 0 ; i < max_iterations ; ++i)
            {
                result += series_coefficient_f1[i] * x;
                x *= lnz;
            }
            result += (1.0 / 6.0) * lnz * lnz * lnz * (11.0 / 6.0 - lnlnz);

            return result;
        }

        // used for |z| > 2.0
        complex<double> g(const complex<double> & z)
        {
            const complex<double> lnz = std::log(z);
            const complex<double> a = - 2.0 / 3.0 * power_of<4>(M_PI);
            const complex<double> b = complex<double>(0.0, -0.5 / M_PI) * lnz;
            double theta = 0.0;

            if ((z.imag() < 0.0) || ((z.real() >= 1.0) && (z.imag() == 0.0)))
            {
                theta = 1.0;
            }

            complex<double> result = a * (-1.0 / 30.0 + b * b * (1.0 - 2.0 * b * (1.0 - 0.5 * b)))
                + complex<double>(0.0, -M_PI / 3.0) * theta * power_of<3>(lnz);

            return result;
        }
    }

    // Calculation of the quadlogarithm based on [C:2006A]
    complex<double> quadlog(const complex<double> & z)
    {
        // special cases
        if (z == complex<double>(0.0, 0.0))
            return complex<double>(0.0);

        if (z == complex<double>(1.0, 0.0))
            return power_of<4>(M_PI) / 90.0;

        if (z == complex<double>(-1.0, 0.0))
            return - power_of<4>(M_PI) * 7.0 / 720.0;

        if (std::abs(z) < 0.5)
            return quadlog_impl::f0(z); // series around the origin

        if (std::abs(z) > 2.0)
            return - quadlog_impl::f0(1.0 / z) + quadlog_impl::g(z) ; // reflection formula

        return quadlog_impl::f1(z);
    }
}
