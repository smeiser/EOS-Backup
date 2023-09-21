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

#include <iostream>
#include <eos/b-decays/bq-to-dq-vec.hh>
#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/k-lcdas.hh>
#include <eos/form-factors/rho-lcdas.hh>
#include <eos/maths/complex.hh>
#include <eos/maths/integrate.hh>
#include <eos/maths/power-of.hh>
#include <eos/maths/polylog.hh>
#include <eos/models/model.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <string>

namespace eos
{
    using std::norm;
    using std::real;
    using std::imag;

    /*
     * Decay: B_q -> D_q V, cf. [BBNS:2000A] (class I only, P = rho^- or Kstar^-)
     */
    template <>
    struct Implementation<BqToDqVector>
    {
        SpecifiedOption opt_model;

        std::shared_ptr<Model> model;

        UsedParameter hbar;

        UsedParameter g_fermi;

        QuarkFlavorOption opt_q;

        UsedParameter m_B;

        UsedParameter f_B;

        UsedParameter tau_B;

        UsedParameter m_D;

        UsedParameter f_D;

        UsedParameter m_V;

        UsedParameter f_V;

        UsedParameter alpha_s;

        std::function<double ()> mu_P;

        SpecifiedOption opt_ff;

        std::shared_ptr<FormFactors<PToP>> ff;

        std::shared_ptr<VectorLCDAs> lcdas;

        SpecifiedOption opt_cp_conjugate;

        bool cp_conjugate;

        UsedParameter mu;

        std::function<complex<double> ()> ckm_factor;
        std::function<WilsonCoefficients<bern::ClassIII> (bool)> wc;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            opt_model(o, options, "model"),
            model(Model::make(opt_model.value(), p, o)),
            opt_q(o, options, "q"),
            hbar(p["QM::hbar"], u),
            g_fermi(p["WET::G_Fermi"], u),
            m_B(p["mass::B_" + opt_q.str()], u),
            f_B(p["decay-constant::B_" + opt_q.str()], u),
            m_D(p["mass::D_" + opt_q.str()], u),
            f_D(p["decay-constant::D_" + opt_q.str()], u),
            m_V(p["mass::" + stringify(opt_q.value() == QuarkFlavor::down ? "K_u^*" : "rho^+")], u),
            f_V(p["decay-constant::" + stringify(opt_q.value() == QuarkFlavor::down ? "K_u^*" : "rho")], u),
            alpha_s(p["QCD::alpha_s(MZ)"], u),
            tau_B(p["life_time::B_" + opt_q.str()], u),
            opt_ff(o, options, "form-factors"),
            opt_cp_conjugate(o, options, "cp-conjugate"),
            cp_conjugate(destringify<bool>(opt_cp_conjugate.value())),
            mu(p[opt_q.str() + "bcu::mu"], u)
        {
            Context ctx("When constructing B_q->D_q V observable");

            // handle the spectator quark flavor q
            switch (opt_q.value())
            {
                case QuarkFlavor::strange:
                    ckm_factor = [this]() { return conj(model->ckm_ud()) * model->ckm_cb(); };
                    wc         = [this](const bool & cp_conjugate) { return model->wet_dbcu(cp_conjugate); };
                    ff         = FormFactorFactory<PToP>::create("B_s->D_s::" + opt_ff.value(), p, o);
                    lcdas      = VectorLCDAs::make("rho", p, o);
                    break;
                case QuarkFlavor::down:
                    ckm_factor = [this]() { return conj(model->ckm_us()) * model->ckm_cb(); };
                    wc         = [this](const bool & cp_conjugate) { return model->wet_sbcu(cp_conjugate); };
                    ff         = FormFactorFactory<PToP>::create("B->D::" + opt_ff.value(), p, o);
                    lcdas      = VectorLCDAs::make("Kstar", p, o);
                    break;
                default:
                    throw InternalError("Invalid quark flavor: " + stringify(opt_q.value()));
            }
            u.uses(*model);
            u.uses(*ff);
        }
        // auxiliary functions for hard-scattering kernels
        complex<double> AVLL(const double & u, const double & z) const
        {
            // evaluate AVLL chosing the correct branch
            if(0.0 < z && z < 1.0)
            {
                return
                    log(u * (1.0 - z * z)) / (1.0 - u * (1.0 - z * z)) - pow(log(u * (1.0 - z * z)), 2) - dilog(1.0 - u * (1.0 - z * z));
            }
            if (z > 1.0)
            {
                return
                    (log(-(u * (1.0 - z * z))) - 1i * M_PI) / (1.0 - u * (1.0 - z * z)) +
                    pow(log(1 - u * (1.0 - (z * z))), 2) / 2.0 - pow(log(-(u * (1.0 - z * z))) - 1i * M_PI, 2) +
                    dilog(1.0 / (1.0 - u * (1.0 - z * z))) +
                    -1.0 / 3.0 * M_PI * M_PI - 1i * M_PI * log(1 - u * (1.0 - z * z));
            }
        };
        complex<double> fVLL(const double & u, const double & z) const
        {
            // evaluate fVLL chosing the correct branch
            if (0.0 < z && z < 1.0)
            {
                return
                    2.0 * (AVLL(u, z) - AVLL(1.0 - u, z)) - (z / (1.0 - u * (1 - z * z))) -
                    (u * (1.0 - z * z) * (z + 3.0 * (1.0 - u * (1.0 - z * z))) * log(u * (1.0 - z * z))) / pow(1.0 - u * (1.0 - z * z), 2);
            }
            if (z > 1.0)
            {
                return
                    2.0 * (AVLL(u, z) - AVLL(1.0 - u, z)) - (z / (1.0 - u * (1 - z * z))) -
                    (u * (1.0 - z * z) * (z + 3.0 * (1.0 - u * (1.0 - z * z))) * (log(-u * (1.0 - z * z)) - 1i * M_PI)) / pow(1.0 - u * (1.0 - z * z), 2);
            }
        };
        complex<double> ASLR(const double & u, const double & z) const
        {
             // evaluate ASLR chosing the correct branch
            if (0.0 < z && z < 1.0)
            {
                return
                    z * z / (pow(1.0 + z, 2) * (1.0 - u * (1.0 - z * z))) +
                    ((-2.0 + u * u * pow(-1.0 + z, 2) * (2.0 + 4.0 * z + 3.0 * z * z)) * log(u * (1.0 - z * z))) / pow(1.0 - u * (1 - z * z), 2) +
                    2.0 * ((2.0 * log(u * (1.0 - z * z))) / (1.0 - u * (1 - z * z)) - pow(log(u * (1.0 - z * z)), 2) - dilog(1.0 - u * (1.0 - z * z)));
            }
            if (z > 1.0)
            {
                return
                    z * z / (pow(1 + z, 2) * (1.0 - u * (1.0 -z * z))) +
                    ((-2.0 + u * u * pow(-1.0 + z, 2) * (2.0 + 4.0 * z + 3.0 * z * z))* (-1i * M_PI + log(-(u * (1 - z * z))))) / pow(1.0 - u * (1 - z * z), 2) +
                    2.0 * (-1.0 / 3.0 * M_PI * M_PI + (2.0 * (-1i * M_PI + log(-(u * (1.0 - z * z))))) / (1.0 - u * (1 - z * z)) -
                    pow(-1i * M_PI + log(-(u * (1 -z * z))), 2) - 1i * M_PI * log(1.0 - u * (1.0 -z * z)) + pow(log(1.0 - u * (1.0 - z * z)), 2) / 2.0 +
                    dilog(1.0 / (1.0 - u * (1.0 - z * z))));
            }
        };
        complex<double> fSLR(const double & u, const double & z) const
        {
            return
                ASLR(u, z) - ASLR(1.0 - u, z);
        };
        complex<double> ATLL(const double & u, const double & z) const
        {
            if (0.0 < z && z < 1.0)
            {
                return
                    ((-1.0 + u * (2.0 - u - 2.0 * z  + (-2.0 + u) * z * z)) * log(u * (1.0 - z * z))) / (1.0 - u * (1.0 - z * z)) +
                    (1.0 - 2.0 * u) * (pow(log(u * (1.0 - z * z)), 2) + dilog(1.0 - u * (1.0 - z * z)));
            }
            if (z > 1.0)
            {
                return
                    ((-1.0 + u * (2.0 - u - 2.0 * z + (-2.0 + u) * z * z)) * (-1i * M_PI + log(-(u * (1.0 - z * z))))) / (1.0 - u * (1.0 - z * z)) +
                    (1.0 - 2.0 * u) * (M_PI * M_PI / 3.0 + pow(-1i * M_PI + log(-(u * (1.0 - z * z))), 2) + 1i * M_PI * log(1.0 - u * (1.0 - z * z)) -
                    pow(log(1.0 - u * (1.0 - z * z)), 2) / 2.0 - dilog(1.0 / (1.0 - u * (1.0 - z * z))));
            }
        };
        complex<double> fTLL(const double & u, const double & z) const
        {
            return
                -((8.0 * (4.0 * u + 3.0)) / (1.0 + z)) + (8.0 * (1.0 - z)) / (1.0 + z) * (ATLL(u, z) + ATLL(1.0 - u, z));
        };
        complex<double> a_1() const
        {
            const WilsonCoefficients<bern::ClassIII> wc = this->wc(cp_conjugate);

            // cf. [BBNS:2000A], converted to the Bern basis
            const double mb = model->m_b_msbar(mu());
            const double mc = model->m_c_msbar(mu());
            const double mu_L = m_V * lcdas->fperp(mu()) / f_V;
            const double z = mc / mb;
            //const double f_3P = lcdas->f3(mu());
            const double a_s_mu = model->alpha_s(mu()) / (4.0 * M_PI);
            const complex<double> a_1_lo =
                -1.0 / 3.0 * (wc.c1() + wc.c1p())  - 4.0 / 9.0 * (wc.c2() + wc.c2p()) -
                16.0 / 3.0 * (wc.c3() + wc.c3p()) - 64.0 / 9.0 * (wc.c4() + wc.c4p()) +
                1.0 / 6.0 * (wc.c5() + wc.c5p()) + 2.0 / 9.0 * (wc.c6() + wc.c6p()) +
                8.0 / 3.0 * (wc.c9() + wc.c9p()) + 32.0 / 9.0 * (wc.c10() + wc.c10p());
            auto a_1_nlo_integrand = [&](const double & u) -> complex<double>
            {
                static const double eps = 1.0e-10;
                complex<double> TVLL ;
                if ((u < eps) || (u > 1.0 - eps))
                {
                    TVLL = 0.0;
                }
                else
                {
                    TVLL = (-18.0 - 6.0 * 2.0 * log(mu() / mb) + fVLL(1.0 - u, 1.0 / z) + fVLL(u, z) + (3.0 + 2.0 * log(u / (1.0 - u))) * log(z * z)) * this->lcdas->phipara(u, mu());
                };
                complex<double> TVLR ;
                if ((u < eps) || (u > 1.0 - eps))
                {
                    TVLR = 0.0;
                }
                else
                {
                    TVLR = (6.0 + 6.0 * 2.0 * log(mu() / mb) - (3.0 + 2.0 * log((1.0 - u) / u)) * log(z * z) - fVLL((1.0 - u), z) - fVLL(u, 1.0 / z)) * this->lcdas->phipara(u, mu());
                };
                complex<double> TSLR ;
                if ((u < eps) || (u > 1.0 - eps))
                {
                    TSLR = 0.0;
                }
                else
                {
                    TSLR = (2.0 * log(u / (1.0 - u)) * log(z * z) - 6.0 + fSLR(u, z) + fSLR((1.0 - u), 1.0 / z)) * lcdas->phiperp(u, mu());
                };
                // TSLL vanishes after integration since LCDA in two-particle limit is just unity and the hard-scattering kernel is antisymetric
                const double TSLL = 0.0;
                complex<double> TTLL ;
                if ((u < eps) || (u > 1.0 - eps))
                {
                    TTLL = 0.0;
                }
                else
                {
                    TTLL = (-48.0 * 2.0 * log(mu() / mb) + 8.0 * (3.0 + ((u - (1.0 - u)) * (1.0 - z)) / (1.0 + z) * log(u / (1.0 - u))) * log(z * z) + fTLL(u, z) + fTLL(1.0 - u, 1.0 / z)) * lcdas->phiperp(u, mu());
                };
                return
                  /*  4.0 / 9.0  * (wc.c1() + wc.c1p()) * (-((2.0 * mu_L * TSLR) / (mb - mc)) - TVLL + 4.0) +
                    2.0 / 27.0 * (wc.c2() + wc.c2p()) * ((2.0 * mu_L * TSLR) / (mb-mc) + TVLL + 14.0) +
                    32.0 /9.0 * (wc.c3() + wc.c3p()) * (-((mu_L * TSLR) / (mb - mc)) - 2.0 * (TVLL + 5.0)) +
                    16.0 /27.0 * (wc.c4() + wc.c4p()) * ((mu_L * TSLR) / (mb - mc) +2.0 * TVLL + 19.0) +
                    1.0 / 18.0 * (wc.c5() + wc.c5p()) * ((mu_L * (TTLL - 4.0 * TSLL)) / (mb - mc) + 4.0 * (TVLR + 6.0)) +
                    1.0 / 108.0 * (wc.c6() + wc.c6p()) * ((mu_L * (4.0 * TSLL - TTLL)) / (mb - mc) - 4.0 * TVLR + 84.0) +
                    1.0 / 9.0 * (wc.c7() + wc.c7p()) * (32.0 - (2.0 * mu_L * (12.0 * TSLL + TTLL)) / (mb - mc))+
                    1.0 / 27.0 * (wc.c8() + wc.c8p()) * ((mu_L * (12.0 * TSLL + TTLL)) / (mb - mc) + 56.0) +
                    32.0 / 9.0 * (wc.c9() + wc.c9p()) * ((2.0 * mu_L * (4.0 * TSLL + TTLL)) / (mb - mc) + TVLR - 40.0) +
                    16.0 / 27.0 * (wc.c10() + wc.c10p()) * (-((2.0 * mu_L * (4.0 * TSLL + TTLL)) / (mb - mc)) - TVLR + 76.0);*/
                    TVLL;
            };
           /* // convoluted 3-particle hard-scattering kernels
            const double TVLL_nlp = 4.0 * (5.0 * lcdas->kappa4(mu()) * m_V * m_V) / (3.0 * (mb * mb - mc * mc));
            const double TTLL_nlp = 4.0 * (3.0 - lcdas->omega3(mu())) / pow(1.0 + z, 2);
            // calculate contributions from three-particle light-meson states
            const complex<double> a_1_nlp =
                    ((wc.c1() - wc.c1p()) * TVLL_nlp) / 3.0 - ((wc.c2() - wc.c2p()) * TVLL_nlp) / 18.0 + (16.0 * (wc.c3() - wc.c3p()) * TVLL_nlp) /3.0 -
                    (8.0 * (wc.c4() - wc.c4p()) * TVLL_nlp) / 9.0 +
                    (wc.c5() - wc.c5p()) * (-1.0 / 6.0 * TVLL_nlp - (f_3P * m_V * m_V * TTLL_nlp) / (24.0 * f_V() * mb * mb * (mb - mc))) +
                    ((wc.c6() - wc.c6p()) * (4.0 * TVLL_nlp + (f_3P * m_V * m_V* TTLL_nlp) / (f_V() * mb * mb * (mb - mc)))) / 144.0 +
                    ((wc.c7() - wc.c7p()) * f_3P * m_V() * m_V()* TTLL_nlp) / (6.0 * f_V() * mb * mb * (mb - mc)) -
                    ((wc.c8() - wc.c8p()) * f_3P * m_V() * m_V() * TTLL_nlp) / (36.0 * f_V() * mb * mb *(mb - mc)) +
                    (8.0 * (wc.c9() - wc.c9p()) * (-TVLL_nlp - (2.0 * f_3P * m_V * m_V * TTLL_nlp) / (f_V() * mb * mb * (mb - mc)))) / 3.0 +
                    (4.0 * (wc.c10() - wc.c10p()) * (TVLL_nlp + (2.0 * f_3P * m_V * m_V * TTLL_nlp) / (f_V()* mb * mb * (mb - mc)))) / 9.0;*/
            auto a_1_nlo_integrand_re = [this, & a_1_nlo_integrand](const double & u) -> const double
            {
                return real(a_1_nlo_integrand(u));
            };
            auto a_1_nlo_integrand_im = [this, & a_1_nlo_integrand](const double & u) -> const double
            {
                return imag(a_1_nlo_integrand(u));
            };
            const double a_1_nlo_re = integrate<GSL::QAGS>(a_1_nlo_integrand_re, 0.0, 1.0);
            const double a_1_nlo_im = integrate<GSL::QAGS>(a_1_nlo_integrand_im, 0.0, 1.0);
            complex<double> a_1_nlo = a_1_nlo_re + a_1_nlo_im * 1i;
            std::cerr << "a1para(mu): " << lcdas->a1para(mu()) << std::endl;
            std::cerr << "a2para(mu): " << lcdas->a2para(mu()) << std::endl;
            std::cerr << "a3para(mu): " << lcdas->a3para(mu()) << std::endl;
            std::cerr << "a4para(mu): " << lcdas->a4para(mu()) << std::endl;
            std::cerr << "a1perp(mu): " << lcdas->a1perp(mu()) << std::endl;
            std::cerr << "a2perp(mu): " << lcdas->a2perp(mu()) << std::endl;
            std::cerr << "a3perp(mu): " << lcdas->a3perp(mu()) << std::endl;
            std::cerr << "a4perp(mu): " << lcdas->a4perp(mu()) << std::endl;
            std::cerr << "mu_V: " << mu_L << std::endl;
            std::cerr << "mb(mu): " << mb << std::endl;
            std::cerr << "mc(mu): " << mc << std::endl;
            std::cerr << "mu: " << mu() << std::endl;
            std::cerr << "TVLR: " << a_1_nlo << std::endl;
            //return a_1_lo + a_s_mu * a_1_nlo;
            return a_1_nlo;
        };

        double decay_width() const
        {
            // cf. [BBNS:2000A], eq. (212), p. 80
            const complex<double> amplitude = g_fermi() / sqrt(2.0) * ckm_factor() * f_V() * ff->f_plus_T(m_V * m_V)
                * sqrt(lambda(m_B * m_B, m_D * m_D, m_V * m_V)) * this->a_1();
            // cf. [BBNS:2000A], eq. (216), p. 80
            const double breakup_momentum = sqrt(lambda(m_B * m_B, m_D * m_D, m_V * m_V)) / (2.0 * m_B);

            // cf. [BBNS:2000A], eq. (221), p. 81
            return norm(amplitude) * breakup_momentum / (8.0 * M_PI * m_B * m_B);
        }

        double branching_ratio() const
        {
            return decay_width() * tau_B / hbar;
        }
    };

    const std::vector<OptionSpecification>
    Implementation<BqToDqVector>::options
    {
        Model::option_specification(),
        FormFactorFactory<PToP>::option_specification(),
        { "cp-conjugate", { "true", "false" },  "false" },
        { "q",            { "s", "d" }                  },
        { "P",            { "Kstar", "rho"}                  }
    };

    BqToDqVector::BqToDqVector(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BqToDqVector>(new Implementation<BqToDqVector>(parameters, options, *this))
    {
    }

    BqToDqVector::~BqToDqVector()
    {
    }

    double
    BqToDqVector::branching_ratio() const
    {
        return _imp->branching_ratio();
    }

    double
    BqToDqVector::decay_width() const
    {
        return _imp->decay_width();
    }

    double
    BqToDqVector::re_a_1() const
    {
        return real(_imp->a_1());
    }

    double
    BqToDqVector::im_a_1() const
    {
        return imag(_imp->a_1());
    }

    const std::set<ReferenceName>
    BqToDqVector::references
    {
        "BBNS:2000A"_rn
    };

    std::vector<OptionSpecification>::const_iterator
    BqToDqVector::begin_options()
    {
        return Implementation<BqToDqVector>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BqToDqVector::end_options()
    {
        return Implementation<BqToDqVector>::options.cend();
    }
}
