/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_NONLOCAL_FORM_FACTORS_LONG_DISTANCE_HH
#define EOS_GUARD_SRC_NONLOCAL_FORM_FACTORS_LONG_DISTANCE_HH 1

#include <complex>

namespace eos
{
    using std::complex;

    struct LongDistance
    {
        static complex<double> g_had_ccbar(const double & s, const double & m_c);
    };
}

#endif
