/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/rare-b-decays/charm-loops.hh>
#include <src/rare-b-decays/exclusive-b-to-s-dilepton.hh>
#include <src/rare-b-decays/form_factors.hh>
#include <src/utils/concrete_observable.hh>
#include <src/utils/integrate.hh>
#include <src/utils/kinematic.hh>
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/qcd.hh>

#include <cmath>
#include <tr1/functional>
#include <utility>
#include <map>
#include <vector>

#include <gsl/gsl_sf.h>

namespace wf
{
    using std::norm;

    // Large Recoil

    template <>
    struct Implementation<BToKstarDilepton<LargeRecoil>>
    {
        Parameter c1;

        Parameter c2;

        Parameter c3;

        Parameter c4;

        Parameter c5;

        Parameter c6;

        Parameter c7;

        Parameter c7prime;

        Parameter c8;

        Parameter c9;

        Parameter c9prime;

        Parameter c10;

        Parameter c10prime;

        Parameter m_b_MSbar;

        Parameter m_c;

        Parameter m_B;

        Parameter m_Kstar;

        double m_l;

        Parameter mu;

        double lambda_t2;

        complex<double> lambda_u_hat;

        double f_B;

        double f_Kstar_par;

        double f_Kstar_perp;

        double lambda_B_p;

        Parameter a_1_par;

        Parameter a_2_par;

        Parameter a_1_perp;

        Parameter a_2_perp;

        Parameter ckm_A;

        Parameter ckm_lambda;

        Parameter uncertainty_par_left;

        Parameter uncertainty_par_right;

        Parameter uncertainty_perp_left;

        Parameter uncertainty_perp_right;

        Parameter uncertainty_long_left;

        Parameter uncertainty_long_right;

        double e_q;

        std::tr1::shared_ptr<FormFactors<BToKstar>> form_factors;

        Implementation(const Parameters & p, const ObservableOptions & o) :
            c1(p["c1"]),
            c2(p["c2"]),
            c3(p["c3"]),
            c4(p["c4"]),
            c5(p["c5"]),
            c6(p["c6"]),
            c7(p["c7"]),
            c8(p["c8"]),
            c9(p["c9"]),
            c10(p["c10"]),
            c7prime(p["c7prime"]),
            c9prime(p["c9prime"]),
            c10prime(p["c10prime"]),
            m_b_MSbar(p["mass::b(MSbar)"]),
            m_c(p["mass::c"]),
            m_B(p["mass::B0"]),
            m_Kstar(p["mass::K^*0"]),
            mu(p["mu"]),
            f_B(0.180), // +/- 0.03 GeV, cf. [BHP2008], Table 1, p. 8
            f_Kstar_par(0.225), // +/- 0.005 GeV, cf. [BHP2008], Table 1, p. 8
            f_Kstar_perp(0.185), // +/-0.005 GeV, cf. [BHP2008], Table 1, p. 8
            lambda_B_p(0.458),// +/- 0.115 GeV, cf. [BHP2008], Table 1, p. 8
            a_1_par(p["B->K^*::a_1_par"]),
            a_2_par(p["B->K^*::a_1_par"]),
            a_1_perp(p["B->K^*::a_1_perp"]),
            a_2_perp(p["B->K^*::a_2_perp"]),
            ckm_A(p["CKM::A"]),
            ckm_lambda(p["CKM::lambda"]),
            uncertainty_par_left(p["B->K^*ll::A_par^L_uncertainty"]),
            uncertainty_par_right(p["B->K^*ll::A_par^R_uncertainty"]),
            uncertainty_perp_left(p["B->K^*ll::A_perp^L_uncertainty"]),
            uncertainty_perp_right(p["B->K^*ll::A_perp^R_uncertainty"]),
            uncertainty_long_left(p["B->K^*ll::A_0^L_uncertainty"]),
            uncertainty_long_right(p["B->K^*ll::A_0^R_uncertainty"]),
            e_q(-1.0/3.0)
        {
            form_factors = FormFactorFactory<BToKstar>::create(o["form-factors"], p);
            if (! form_factors.get())
                throw std::string("InternalError");

            // TODO: Lepton masses, m_l = m_mu
            m_l = 0.10565836; // (GeV), cf. [PDG2008], p. 13
        }

        double beta_l(const double & s) const
        {
            return std::sqrt(1.0 - 4.0 * m_l * m_l / s);
        }

        double norm(const double & s) const
        {
            static const double alpha_e = 1.0 / 133.0; // cf. [BHP2008]
            static const double g_fermi = 1.16637e-5; // (Gev^-2 (hbar c)^3), cf. [PDG2008], p.5

            double lambda_t = ckm_A * ckm_lambda * ckm_lambda;

            return std::sqrt(g_fermi * g_fermi * alpha_e * alpha_e / 3.0 / 1024 / std::pow(M_PI, 5.0) / m_B
                    * lambda_t * lambda_t * s_hat(s)
                    * std::sqrt(lambda(m_B * m_B, m_Kstar * m_Kstar, s))); // cf. [BHP2008], Eq. (C.6), p. 21
        }

        inline double s_hat(double s) const
        {
            return s / m_B / m_B;
        }

        inline double energy(const double & s) const
        {
            return (m_B * m_B + m_Kstar * m_Kstar - s) / (2.0 * m_B);
        }

        inline double mu_f() const
        {
            static const double Lambda_QCD = 0.5; // (GeV)
            return std::sqrt(mu * Lambda_QCD);
        }

        inline double m_b_PS() const
        {
            // Actually use the PS mass at mu_f = 2 GeV
            return QCD::mb_PS(m_b_MSbar, mu, 2.0);
        }

        /* Effective wilson coefficients */
        // cf. [BFS2001], below Eq. (9), p. 4
        double c7eff() const
        {
            return c7() - 1.0/3.0 * c3() - 4.0/9.0 * c4() - 20.0/3.0 * c5() - 80.0/9.0 * c6();
        }

        // cf. [BFS2001], below Eq. (26), p. 8
        double c8eff() const
        {
            return c8() + c3() - 1.0/6.0 * c4() + 20.0 * c5() - 10.0/3.0 * c6();
        }

        // cf. [BFS2001], Eq. (10), p. 4
        complex<double> Y0(const double & s) const
        {
            double Y_c = 4.0 / 3.0 * c1() + c2() + 6.0 * c3() + 60.0 * c5();
            double Y_b = -0.5 * (7.0 * c3() + 4.0 / 3.0 * c4() + 76.0 * c5() + 64.0 / 3.0 * c6());
            double Y_0 = -0.5 * (c3() + 4.0 / 3.0 * c4() + 16.0 * c5() + 64 / 3.0 * c6());
            double Y = 2.0 / 9.0 * (6.0 * c3() + 32.0 * c5() + 32.0 / 3.0 * c6());

            // Uses b pole mass according to [BFS2001], Sec. 3.1, paragraph Quark Masses
            return Y_c * CharmLoops::h(mu, s, m_c)
                + Y_b * CharmLoops::h(mu, s, m_b_PS())
                + Y_0 * CharmLoops::h(mu, s)
                + Y;
        }

        // cf. [BFS2004], ?
        complex<double> Y0u(const double & s) const
        {
            double a = 4.0 / 3.0 * c1() + c2();

            return a * (CharmLoops::h(mu, s, m_c) - CharmLoops::h(mu, s));
        }

        /* Form factors */
        //  cf. [BHP2008], Eq. (E.4), p. 23
        double xi_perp(const double & s) const
        {
            const double factor = m_B / (m_B + m_Kstar);
            double result = factor * form_factors->v(s_hat(s));

            return result;
        }

        double xi_par(const double & s) const
        {
            const double factor1 = (m_B + m_Kstar) / (2.0 * energy(s));
            const double factor2 = (1.0 - m_Kstar / m_B);
            double result = factor1 * form_factors->a_1(s_hat(s)) - factor2 * form_factors->a_2(s_hat(s));

            return result;
        }

        /* NLO functions */
        // cf. [BFS2001], Eq. (29), p. 8
        complex<double> B0(const double & s, const double & m_q) const
        {
            double z = 4.0 * m_q * m_q / s;
            double rp, ip;

            if (m_q < 0.0)
            {
                throw InternalError("Implementation<BToKstarDilepton<LargeRecoil>>::B0: m_q < 0!");
            };

            if ((0.0 == m_q) && (0.0 == s))
            {
                throw InternalError("Implementation<BToKstarDilepton<LargeRecoil>>::B0: m_q == 0 & s == 0");
            }

            if (0 == s)
                return complex<double>(-2.0, 0.0);

            if (z > 1.0)
            {
                rp = -2.0 * std::sqrt(z - 1.0) * std::atan(1.0 / std::sqrt(z - 1.0));
                ip = 0.0;
            }
            else
            {
                rp = std::sqrt(1.0 - z) * std::log((1.0 - std::sqrt(1 - z)) / (1.0 + std::sqrt(1.0 - z)));
                ip = std::sqrt(1.0 - z) * M_PI;
            }

            return complex<double>(rp, ip);
        }

        // cf. [BFS2001], Eq. (84), p. 30
        double C0(const double & s) const
        {
            static const int points = 40;
            // Integration boundaries of u = (0, u_max]
            static const double dx = (1.0 - 0.0) / points;
            static const double g[2] =
            {
                (1.0 + std::sqrt(3.0 / 5.0)) / 2.0,
                (1.0 - std::sqrt(3.0 / 5.0)) / 2.0
            };

            long double s_hat = s / m_B / m_B;
            long double x1, x2, x;
            long double y[3];
            long double result = 0;

            for (int i = 0; i < points; i++)
            {
                // 0 is lower integration boundary

                x1 = 0.0 + i * dx;
                x2 = 0.0 + (i + 1) * dx;
                for (int j = 0; j < 3; j++)
                {
                    switch (j)
                    {
                        case 0: x = x1 * g[0] + x2 * g[1];
                                break;

                        case 1: x = (x1 + x2) / 2.0;
                                break;

                        case 2: x = x1 * g[1] + x2 * g[0];
                                break;
                    }

                    y[j] = std::log(x * x / (1.0 - s_hat * x * (1.0 - x))) / (1.0 + x * (1.0 - s_hat));
                }

                result += (5.0 * y[0] + 8.0 * y[1] + 5.0 * y[2]) * dx / 18.0;
            };

            return result;
        }

        // cf. [BFS2001], Eqs. (30)-(32), p. 8
        complex<double> I1(const double & s, const double & u, const double & m_q) const
        {
            if (m_q == 0.0)
                return complex<double>(1.0, 0.0);

            int status;
            gsl_sf_result res_re, res_im;

            double ubar = 1.0 - u;
            double m_q2 = m_q * m_q;
            double m_B2 = m_B * m_B;

            double a, a2, sign;
            complex<double> dilogArg, dilog1, dilog2;
            complex<double> LxpLxm, LypLym;

            if (1.0 - 4.0 * m_q2 / (m_B2 - u * (m_B2 - s)) > 0)
            {
                a = (1 - std::sqrt(1.0 - 4.0 * m_q2 / (m_B2 - u * (m_B2 - s))))
                  / (1 + std::sqrt(1.0 - 4.0 * m_q2 / (m_B2 - u * (m_B2 - s))));
                LxpLxm = -M_PI* M_PI / 3.0
                    + std::log(a) * (std::log(a) + complex<double>(0.0, M_PI))
                    + gsl_sf_dilog(-a) + gsl_sf_dilog(-1.0/a);
            }
            else
            {
                a  = sqrt(4.0 * m_q2 / (m_B2 - u * (m_B2 - s)) - 1);
                a2 = a * a;

                if (a2 - 1 > 0)
                    sign = +1.0;
                else
                    sign = -1.0;

                dilogArg = complex<double>((a2 - 1.0)/(a2 + 1.0), -2.0 * a / (a2 + 1.0));
                status = gsl_sf_complex_dilog_e(abs(dilogArg), arg(dilogArg), &res_re, &res_im);
                dilog1 = complex<double>( res_re.val, res_im.val);

                dilogArg = complex<double>((a2 - 1.0)/(a2 + 1.0), +2.0 * a / (a2 + 1.0));
                status = gsl_sf_complex_dilog_e(abs(dilogArg), arg(dilogArg), &res_re, &res_im);
                dilog2 = complex<double>(res_re.val, res_im.val);

                LxpLxm = -1.0 / 3.0 * M_PI * M_PI - std::atan(2.0 * a / (a2 - 1.0)) * (std::atan(2.0 * a / (a2 - 1.0)) - M_PI * sign)
                    + dilog1 + dilog2;
            }

            if (1.0 - 4.0 * m_q2 / s > 0)
            {
                a = (1.0 - std::sqrt(1.0 - 4.0 * m_q2 / s))
                  / (1.0 + std::sqrt(1.0 - 4.0 * m_q2 / s));
                LypLym = -1.0 / 3.0 * M_PI * M_PI + std::log(a) * (std::log(a) + complex<double>(0.0, M_PI))
                    + gsl_sf_dilog(-a) + gsl_sf_dilog(-1./a);
            }
            else
            {
                a  = std::sqrt(4.0 * m_q2 / s - 1.0);
                a2 = a * a;
                if (a2 - 1.0 > 0)
                    sign = +1.0;
                else
                    sign = -1.0;

                dilogArg = complex<double>((a2 - 1.0) / (a2 + 1.0), -2.0 * a / (a2 + 1.0));
                status = gsl_sf_complex_dilog_e(abs(dilogArg), arg(dilogArg), &res_re, &res_im);
                dilog1 = complex<double>(res_re.val, res_im.val);

                dilogArg = complex<double>((a2 - 1.0) / (a2 + 1.0), +2.0 * a / (a2 + 1.0));
                status = gsl_sf_complex_dilog_e(abs(dilogArg), arg(dilogArg), &res_re, &res_im);
                dilog2 = complex<double>(res_re.val, res_im.val);

                LypLym = -1.0 / 3.0 * M_PI * M_PI - std::atan(2.0 * a / (a2 - 1.0)) * (std::atan(2.0 * a / (a2 - 1.0)) - M_PI * sign)
                    + dilog1 + dilog2;
            }

            return 1.0 + 2.0 * m_q2 / (ubar * (m_B2 - s)) * (LxpLxm - LypLym);
        }

        // cf. [BFS2001], Eq. (48)
        double phi_K(const double & u, const double & a_1, const double & a_2) const
        {
            double xi = 2.0 * u - 1.0;

            return 6.0 * u * (1 - u) * (1.0 + a_1 * 3.0 * xi + a_2 * (7.5 * xi * xi - 1.5));
        }

        // cf. [BFS2001], Eq. (18), p. 6, modulo the inverse of lambda_B,-
        double T0_perp_m(const double & s, const double & u) const
        {
            double wilson = c3 + 4.0/3.0 * c4 + 16.0 * c5 + 64.0/3.0 * c6;
            return -1.0 * e_q * (4.0 * m_B / m_b_PS()) * wilson;
        }

        // cf. [BFS2001], Eq. (20), p. 6, module the inverse of lambda_B,+
        double Tf_perp_p(const double & h, const double & s, const double & u) const
        {
            return (c7eff() + h * c7prime()) * (2.0 * m_B / (1.0 - u) / energy(s));
        }

        // cf. [BFS2001], Eq. (21), p. 6
        complex<double> Tf_par_p(const double & s, const double & u) const
        {
            double E = energy(s);

            return ((c7eff() - c7prime) + (s / (2.0 * m_B * m_b_PS())) * Y0(s)) * (2.0 * m_B * m_B / (1 - u) / E / E);
        }

        // cf. [BFS2001], Eq. (27), p. 8
        complex<double> t_perp(const double & s, const double & u, const double & m_q) const
        {
            if (0.0 == s)
                return t_perp_0(u, m_q);

            double ubar = 1.0 - u;
            double E = energy(s);
            double x = ubar * m_B * m_B + u * s;

            complex<double> result = (2.0 * m_B / ubar / E) * I1(s, u, m_q);
            if (m_q > 0.0)
                result = result + (s / ubar / ubar / E / E) * (B0(x, m_q) - B0(s, m_q));

            return result;
        }

        complex<double> t_perp_0(const double & u, const double & m_q) const
        {
            double ubar = 1.0 - u;
            double m_q2 = m_q * m_q, m_B2 = m_B * m_B;
            double a, a2, sign;
            complex<double> dilogArg, dilog1, dilog2;
            complex<double> LxpLxm;
            int status;
            gsl_sf_result res_re, res_im;

            if (m_q > 0)
            { // m != 0
                if (1.0 - 4.0 * m_q2 / (m_B2 - u * m_B2) > 0)
                {
                    a = (1.0 - std::sqrt(1.0 - 4.0 * m_q2 / (m_B2 - u * m_B2)))
                        / (1.0 + std::sqrt(1.0 - 4.0 * m_q2 / (m_B2 - u * m_B2)));
                    LxpLxm = -M_PI * M_PI / 3.0 + std::log(a) * (std::log(a) + complex<double>(0.0, M_PI))
                        + gsl_sf_dilog(-a) + gsl_sf_dilog(-1.0/a);
                }
                else
                {
                    a = std::sqrt(4.0 * m_q2 / (m_B2 - u * m_B2) - 1.0);
                    a2 = a * a;

                    if (a2 - 1 > 0)
                        sign = +1.0;
                    else
                        sign = -1.0;

                    dilogArg = complex<double>((a2 - 1.0) / (a2 + 1.0), -2.0 * a / (a2 + 1.0));
                    status = gsl_sf_complex_dilog_e(abs(dilogArg), arg(dilogArg), &res_re, &res_im);
                    dilog1 = complex<double>(res_re.val, res_im.val);
                    dilogArg = complex<double>((a2 - 1.0) / (a2 + 1.0), +2.0 * a / (a2 + 1.0));
                    status = gsl_sf_complex_dilog_e(abs(dilogArg), arg(dilogArg), &res_re, &res_im);
                    dilog2 = complex<double>(res_re.val, res_im.val);

                    LxpLxm = -M_PI * M_PI / 3.0 - std::atan(2.0 * a / (a2 - 1.0)) * (std::atan(2.0 * a / (a2 - 1.0)) - M_PI * sign)
                        + dilog1 + dilog2;
                }

                return 4.0 / ubar * (1.0 + 2.0 * m_q2 / ubar / m_B2 * LxpLxm);
            }
            else
            {
                return complex<double>(4.0 / ubar, 0.0);
            }
        }

        // cf. [BFS2001], Eq. (28), p. 8
        complex<double> t_par(const double & s, const double & u, const double & m_q) const
        {
            double ubar = 1.0 - u;
            double E = energy(s);
            double x = ubar * m_B * m_B + u * s;

            complex<double> result = (2.0 * m_B / ubar / E) * I1(s, u, m_q);
            if (m_q > 0.0)
                result = result + (x / ubar / ubar / E / E) * (B0(x, m_q) - B0(s, m_q));

            return result;
        }

        // cf. [BFS2001], Eq. (23), p. 7, multiplied by phi_K^*,perp
        complex<double> Tnf_perp_p(const double & s, const double & u) const
        {
            double m_b = m_b_PS();
            double s_hat = s / m_B / m_B;
            double ubar = 1.0 - u;

            double a = (4.0 / 3.0 / (u + ubar * s_hat)) * c8eff();
            complex<double> ba = (+2.0 / 3.0) * (-c1 / 6.0 + c2 + 6.0 * c6) * t_perp(s, u, m_c);
            complex<double> bb = (-1.0 / 3.0)
                * (c3 - c4 / 6.0 + 16.0 * c5 + 10.0/3.0 * c6 - (4.0 * m_b / m_B) * (c3 - c4/6.0 + 4.0 * c5 - 2.0/3.0 * c6))
                * t_perp(s, u, m_b);
            complex<double> bc = (-1.0 / 3.0) * (c3 - c4 / 6.0 + 16.0 * c5 - 8.0/3.0 * c6) * t_perp(s, u, 0.0);

            return (a + (m_B / 2.0 / m_b) * (ba + bb + bc)) * phi_K(u, a_1_perp, a_2_perp);
        }

        // cf. [BFS2001, Eq. (25), p. 7, multiplied by phi_K^*,par
        complex<double> Tnf_par_p(const double & s, const double & u) const
        {
            double m_b = m_b_PS();

            complex<double> a = (+2.0 / 3.0) * (-c1 / 6.0 + c2 + 6.0 * c6) * t_par(s, u, m_c);
            complex<double> b = (-1.0 / 3.0) * (c3 - c1 / 6.0 + 16.0 * c5 + 10.0/3.0 * c6) * t_par(s, u, m_b);
            complex<double> c = (-1.0 / 3.0) * (c3 - c4 / 6.0 + 16.0 * c5 - 8.0/3.0 * c6) * t_par(s, u, 0.0);

            return (m_B / m_b) * (a + b + c) * phi_K(u, a_1_par, a_2_par);
        }

        // cf. [BFS2001], Eq. (26), pp. 7-8, multiplied by phi_K^*,par
        complex<double> Tnf_par_m(const double & s, const double & u) const
        {
            double m_b = m_b_PS();

            double s_hat = s / m_B / m_B;
            double ubar = 1.0 - u;
            double x = ubar * m_B * m_B + u * s;

            double a = (e_q * 8.0 / (ubar + u * s_hat)) * c8eff();
            complex<double> ba = (-c1 / 6.0 + c2 + c4 + 10 * c6) * CharmLoops::h(mu, x, m_c);
            complex<double> bb = (c3 + 5.0/6.0 * c4 + 16.0 * c5 + 22.0/3.0 * c6) * CharmLoops::h(mu, x, m_b);
            complex<double> bc = (c3 + 17.0/6.0 * c4 + 16.0 * c5 + 82.0/3.0 * c6) * CharmLoops::h(mu, x);
            double bd = -8.0 / 27.0 * (-7.5 * c4 + 12.0 * c5 - 32.0 * c6);

            return (6.0 * m_B / m_b) * (ba + bb + bc + bd) * phi_K(u, a_1_par, a_2_par);
        }

        // cf. [BFS2001], Eq. (36), p. 9
        double L(const double & s) const
        {
            double m_b = m_b_PS();
            double m_b2 = m_b * m_b;

            return -1.0 * (m_b2 - s) / s * std::log(1.0 - s / m_b2);
        }

        // cf. [BFS2001], Eq. (54), p. 15
        complex<double> lambda_B_m_inv(const double & s) const
        {
            if (0.0 == s)
                return complex<double>(0.0, 0.0);

            double omega_0 = lambda_B_p;
            double x = s / m_B / omega_0;
            double ei = gsl_sf_expint_Ei(x);

            complex<double> result = complex<double>(-ei, M_PI) * (std::exp(-x) / omega_0);

            return result;
        }

        // cf. [BFS2001], Eq. (82), p. 30
        complex<double> F87(const double & s) const
        {
            double m_b = m_b_PS();
            double s_hat = s / m_b / m_b;
            double s_hat2 = s_hat * s_hat;
            double denom = (1.0 - s_hat);
            double denom2 = denom * denom;

            complex<double> a = complex<double>(-32.0 * std::log(mu / m_b) - 8.0 * s_hat / denom * std::log(s_hat)
                    - 4.0 * (11.0 - 16.0 * s_hat + 8.0 * s_hat2) / denom2,
                    -8.0 * M_PI);
            complex<double> b = (4.0 / denom / denom2)
                * ((9.0 * s_hat - 5.0 * s_hat2 + 2.0 * s_hat * s_hat2) * B0(s, m_b) - (4.0 + 2.0 * s_hat) * C0(s));

            return (1.0 / 9.0) * (a + b);
        }

        // cf. [BFS2001], Eq. (83), p. 30
        complex<double> F89(const double & s) const
        {
            double m_b = m_b_PS();
            double s_hat = s / m_b / m_b;
            double denom = (1.0 - s_hat);
            double denom2 = denom * denom;

            double a = 16.0 * std::log(s_hat) / denom + 8.0 * (5.0 - 2.0 * s_hat) / denom2;
            complex<double> b = (-8.0 * (4.0 - s_hat) / denom / denom2) * ((1.0 + s_hat) * B0(s, m_b) - 2.0 * C0(s));

            return (1.0 / 9.0) * (a + b);
        }

        // cf. [BFS2004], Eq. (51) the integrand of the first term only,
        // or [FM2002], Eq. (15) times a factor of 3.0, respectively
        double Twa_perp(const double & s_hat, const double & u) const
        {
            double ubar = 1.0 - u;

            return phi_K(u, a_1_perp, a_2_perp) / (ubar + u * s_hat);
        }

        // cf. [FM2002], Eq. (17), the integrand only
        double Xperp(const double & s_hat, const double & u) const
        {
            double ubar = 1.0 - u;
            double denom = ubar + u * s_hat;

            return phi_K(u, a_1_perp, a_2_perp) * (1 / denom + 1 / denom / denom) / 3.0;
        }

        // cf. [FM2002], Eq. (22), p. 9
        complex<double> FV(const double & s) const
        {
            return 3.0 / 4.0 * (
                    (-c1 / 6.0 + c2 + c4 + 10.0 * c6) * CharmLoops::h(mu, s, m_c)
                    + (c3 + 5.0/6.0 * c4 + 16.0 * c5 + 22.0/3.0 * c6) * CharmLoops::h(mu, s, m_b_PS())
                    + (c3 + 17.0/6.0 * c4 + 16.0 * c5 + 82.0/3.0 * c6) * CharmLoops::h(mu, s)
                    - 8.0/27.0 * (-7.5 * c4 + 12 * c5 - 32 * c6));
        }

        // cf. [BFS2004], Eq. (52), the integrand of the second term only
        complex<double> Thsa_1_perp(const double & s_hat, const double & u) const
        {
            double ubar = 1.0 - u;
            double x = (ubar + u * s_hat) * m_B * m_B;

            return phi_K(u, a_1_perp, a_2_perp) / (ubar + u * s_hat) * FV(x);
        }

        // cf. [BFS2004], Eq. (52), the integrand of the third term only.
        // the v integration has been executed analytically
        complex<double> Thsa_2_perp(const double & s_hat, const double & u) const
        {
            double ubar = 1.0 - u;
            double x = (ubar + u * s_hat) * m_B * m_B;

            return 3.0 * u * u * (1.0 + a_1_par * (4.0 * u - 3.0) + a_2_par * (15.0 * u * u - 20.0 * u + 6.0))
                * FV(x);
        }

        complex<double> tensor_perp_hsa(const double & shat, const double & u) const
        {
            return 12.0 * c8eff() * m_b_PS() / m_B * f_Kstar_perp * Xperp(shat, u)
                    + 8.0 * f_Kstar_perp * Thsa_1_perp(shat, u)
                    - 4.0 * m_Kstar * f_Kstar_par / ((1.0 - shat) * lambda_B_p) * Thsa_2_perp(shat, u);
        }

        // cf. [BFS2001], Eq. (15) with a = perp, and [BHP2008], Eq. (C.4)
        complex<double> tensor_perp(const double & h, const double & s) const
        {
            // cf. [BFS2004], paragraph below Eq. (42)
            double ff_nlo_factor = QCD::alpha_s(mu) * QCD::casimir_f / 4.0 / M_PI;
            double scatt_nlo_factor = QCD::alpha_s(mu_f()) * QCD::casimir_f / 4.0 / M_PI;

            double scatt_factor = M_PI * M_PI / 3.0 * f_B * f_Kstar_perp / m_B;
            double shat = s_hat(s);

            // Here m_b_PS is used instead of m_b_pole, cf. [BFS2001], comment below Eq. (36), p. 9
            double m_b = m_b_PS();

            /* Form factor corrections */
            complex<double> ff_0 = (c7eff() + h * c7prime()) + s / (2.0 * m_b * m_B) * Y0(s);
            // cf. [BFS2001], Eq. (34), p. 9
            double ff_f = (c7eff() + h * c7prime()) * (8.0 * std::log(m_b / mu) - 4.0 * (1.0 - mu_f() / m_b) - L(s));
            // cf. [BFS2001], Eq. (37), p. 9
            complex<double> ff_nf = (-1.0 / QCD::casimir_f) * (
                    (c2 - c1 / 6.0) * CharmLoops::F27(mu, s, m_b) + c8eff() * F87(s)
                    + (s / (2.0 * m_b * m_B)) * (c2() * CharmLoops::F29(mu, s, m_b) + c1() * CharmLoops::F19(mu, s, m_b) + c8eff() * F89(s)));
            complex<double> ff = ff_0 + ff_nlo_factor * (ff_f + ff_nf);

            /* Specator scattering, folded with phi_K^*,perp */
            // cf. [BFS2001], Eq. (20), p. 6
            double scatt_f_p = (c7eff() + h * c7prime) * 2.0 * m_B / energy(s) * 3.0 * (1.0 + a_1_perp + a_2_perp);
            // cf. [BFS2001], Eq. (23), p. 7
            complex<double> scatt_nf_p = integrate(
                    std::tr1::function<complex<double> (const double &)>(
                        std::tr1::bind(&Implementation<BToKstarDilepton<LargeRecoil>>::Tnf_perp_p, this, s, std::tr1::placeholders::_1)),
                    32, 0.001, 0.999);
            complex<double> scatt_p = (1.0 / lambda_B_p) * scatt_nlo_factor * (scatt_f_p + scatt_nf_p);

            /* Weak annihilation */
            // cf. [BFS2004], Eq. (51), p. 26
            double wa = (e_q * 2.0 * M_PI * M_PI * f_B / (3.0 * m_b * m_B)) * (
                    - f_Kstar_perp * (c3 + 4.0/3.0 * c4 + 4.0 * c5 + 16.0/3.0 * c6) * integrate(
                        std::tr1::function<double (const double &)>(
                            std::tr1::bind(&Implementation<BToKstarDilepton<LargeRecoil>>::Twa_perp, this, shat, std::tr1::placeholders::_1)),
                        32, 0.0, 1.0)
                    + f_Kstar_par * m_Kstar / ((1 - shat) * lambda_B_p));
            /* Hard spectator scattering */
            complex<double> hsa = e_q * scatt_nlo_factor * M_PI * M_PI * f_B / (3.0 * m_b * m_B) * integrate(
                    std::tr1::function<complex<double> (const double &)>(
                        std::tr1::bind(&Implementation<BToKstarDilepton<LargeRecoil>>::tensor_perp_hsa, this, shat, std::tr1::placeholders::_1)),
                    32, 0.0, 1.0);

            complex<double> result = xi_perp(s) * ff + scatt_factor * scatt_p + wa + hsa;

            return result;
        }

        // cf. [BFS2001], Eq. (15) with a = par, and [BHP2008], Eq. (C.4)
        complex<double> tensor_par(const double & s) const
        {
            // cf. [BFS2004], paragraph below Eq. (42)
            double ff_nlo_factor = QCD::alpha_s(mu) * QCD::casimir_f / 4.0 / M_PI;
            double scatt_nlo_factor = QCD::alpha_s(mu_f()) * QCD::casimir_f / 4.0 / M_PI;

            double scatt_factor = M_PI * M_PI / 3.0 * f_B * f_Kstar_perp / m_B * m_Kstar / energy(s);

            // Here m_b_PS is used instead of m_b_pole, cf. [BFS2001], comment below Eq. (36), p. 9
            double m_b = m_b_PS();

            /* Form factor corrections */
            complex<double> ff_0 = -1.0 * (c7eff() - c7prime() + m_B / (2.0 * m_b) * (Y0(s)));
            // cf. [BFS2004], Eq. (44), p. 24
            complex<double> ff_f = -1.0 * (c7eff() - c7prime()) * (8.0 * std::log(m_b / mu) + 2.0 * L(s) - 4.0 * (1.0 - mu_f() / m_b))
                    + (m_B / (2.0 * m_b)) * Y0(s) * (2.0 - 2.0 * L(s));
            // cf. [BFS2001], Eq. (38), p. 9
            complex<double> ff_nf = (+1.0 / QCD::casimir_f) * (
                    (c2 - c1 / 6.0) * CharmLoops::F27(mu, s, m_b) + c8eff() * F87(s)
                    + (m_B / (2.0 * m_b)) * (c2() * CharmLoops::F29(mu, s, m_b) + c1() * CharmLoops::F19(mu, s, m_b) + c8eff() * F89(s)));
            complex<double> ff = ff_0 + ff_nlo_factor * (ff_f + ff_nf);

            /* Spectator scattering */
            // cf. [BFS2001], Eq. (18), p. 6
            double scatt_0_m = -e_q * 4.0 * m_B / m_b * (c3 + 4.0/3.0 * c4 + 16.0 * c5 + 64.0/3.0 * c6);
            // cf. [BFS2001], Eq. (21), p. 6
            complex<double> scatt_f_p = (c7eff() - c7prime + (s / (2.0 * m_b * m_B)) * Y0(s)) * (2.0 * m_B * m_B / energy(s) / energy(s))
                    * 3.0 * (1.0 + a_1_par + a_2_par);
            scatt_f_p = complex<double>(((c7eff() - c7prime) * 4.0 * m_B / energy(s)) * 3.0 * (1.0 + a_1_par + a_2_par), 0.0);
            // cf. [BFS2001], Eq. (25), p. 7
            complex<double> scatt_nf_p = integrate(
                    std::tr1::function<complex<double> (const double &)>(std::tr1::bind(&Implementation<BToKstarDilepton<LargeRecoil>>::Tnf_par_p,
                            this, s, std::tr1::placeholders::_1)),
                    32, 0.001, 0.999);
            // cf. [BFS2001], Eq. (26), pp. 7-8
            complex<double> scatt_nf_m = integrate(
                    std::tr1::function<complex<double> (const double &)>(std::tr1::bind(&Implementation<BToKstarDilepton<LargeRecoil>>::Tnf_par_m,
                            this, s, std::tr1::placeholders::_1)),
                    32, 0.001, 0.999);
            complex<double> scatt_p = (1.0 / lambda_B_p) * scatt_nlo_factor * (scatt_f_p + scatt_nf_p);
            complex<double> scatt_m = lambda_B_m_inv(s)  * (scatt_0_m + scatt_nlo_factor * scatt_nf_m);

            complex<double> result = xi_par(s) * ff + scatt_factor * (scatt_p + scatt_m);

            return result;
        }

        /* Amplitudes */
        // cf. [BHP2008], p. 20
        complex<double> a_long(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            double shat = s_hat(s);
            double mbhat = m_b_PS() / m_B;
            double E = 0.5 * (m_B - s / m_B);
            double mKhat = m_Kstar / m_B;
            double lambdahat = lambda(1.0, mKhat * mKhat, s);

            double uncertainty = (1.0 - h) / 2.0 * uncertainty_long_left + (1.0 + h) / 2.0 * uncertainty_long_right;
            double wilson = (c9() - c9prime()) + h * (c10() - c10prime());
            double prefactor = -1.0 / (2.0 * m_Kstar * std::sqrt(s));

            double a = wilson * ((m_B * m_B - m_Kstar * m_Kstar - s) * 2.0 * energy(s) * xi_perp(s)
                -lambda(m_B * m_B, m_Kstar * m_Kstar, s) * m_B / (m_B * m_B - m_Kstar * m_Kstar) * (xi_perp(s) - xi_par(s)));
            complex<double> b = 2 * m_b_PS() * (((m_B * m_B + 3 * m_Kstar * m_Kstar - s) * 2.0 * energy(s) / m_B
                        - lambda(m_B * m_B, m_Kstar * m_Kstar, s) / (m_B * m_B - m_Kstar * m_Kstar)) * tensor_perp(-1.0, s)
                    - lambda(m_B * m_B, m_Kstar * m_Kstar, s) / (m_B * m_B - m_Kstar * m_Kstar) * tensor_par(s));

            return uncertainty * prefactor * (a + b);
        }

        // cf. [BHP2008], p. 20
        complex<double> a_perp(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            double shat = s_hat(s);
            double mbhat = m_b_PS() / m_B;
            double mKhat = m_Kstar / m_B;

            double uncertainty = (1.0 - h) / 2.0 * uncertainty_perp_left + (1.0 + h) / 2.0 * uncertainty_perp_right;
            double prefactor = +std::sqrt(2.0) * m_B * std::sqrt(lambda(1.0, mKhat * mKhat, shat));
            double wilson = (c9() + c9prime()) + h * (c10() + c10prime());

            return uncertainty * prefactor * (wilson * xi_perp(s) + (2.0 * mbhat / shat) * tensor_perp(+1.0, s));
        }

        // cf. [BHP2008], p. 20
        complex<double> a_par(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            double shat = s_hat(s);
            double mbhat = m_b_PS() / m_B;
            double mKhat = m_Kstar / m_B;

            double uncertainty = (1.0 - h) / 2.0 * uncertainty_par_left + (1.0 + h) / 2.0 * uncertainty_par_right;
            double prefactor = -std::sqrt(2.0) * m_B * (1.0 - shat);
            double wilson = (c9() - c9prime()) + h * (c10() - c10prime());

            return uncertainty * prefactor * (wilson * xi_perp(s) + (2.0 * mbhat / shat) * (1.0 - mKhat * mKhat) * tensor_perp(-1.0, s));
        }
    };

    BToKstarDilepton<LargeRecoil>::BToKstarDilepton(const Parameters & parameters, const ObservableOptions & options) :
        PrivateImplementationPattern<BToKstarDilepton<LargeRecoil>>(new Implementation<BToKstarDilepton<LargeRecoil>>(parameters, options))
    {
    }

    BToKstarDilepton<LargeRecoil>::~BToKstarDilepton()
    {
    }

    complex<double>
    BToKstarDilepton<LargeRecoil>::a_long(const Helicity & h, const double & s) const
    {
        return _imp->a_long(h, s);
    }

    complex<double>
    BToKstarDilepton<LargeRecoil>::a_perp(const Helicity & h, const double & s) const
    {
        return _imp->a_perp(h, s);
    }

    complex<double>
    BToKstarDilepton<LargeRecoil>::a_par(const Helicity & h, const double & s) const
    {
        return _imp->a_par(h, s);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_branching_ratio(const double & s) const
    {
        // cf. [PDG2008] : Gamma = hbar / tau_B, pp. 5, 79
        static const double Gamma(6.58211899e-22 * 1e-3 / 1.53e-12);

        return differential_decay_width(s) / Gamma;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_decay_width(const double & s) const
    {
        return _imp->norm(s) * _imp->norm(s) * (norm(a_long(left_handed, s))
            + norm(a_long(right_handed, s))
            + norm(a_perp(left_handed, s))
            + norm(a_perp(right_handed, s))
            + norm(a_par(left_handed, s))
            + norm(a_par(right_handed, s)));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_forward_backward_asymmetry(const double & s) const
    {
        return 1.5 * _imp->norm(s) * _imp->norm(s) / differential_decay_width(s)
            * (real(a_par(left_handed, s) * conj(a_perp(left_handed, s))) - real(a_par(right_handed, s) * conj(a_perp(right_handed, s))));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_unnormalized_forward_backward_asymmetry(const double & s) const
    {
        return 1.5 * _imp->norm(s) * _imp->norm(s)
            * (real(a_par(left_handed, s) * conj(a_perp(left_handed, s))) - real(a_par(right_handed, s) * conj(a_perp(right_handed, s))));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_2(const double & s) const
    {
        double a = norm(a_perp(left_handed, s)) + norm(a_perp(right_handed, s));
        double b = norm(a_par(left_handed, s)) + norm(a_par(right_handed, s));

        return (a - b) / (a + b);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_3(const double & s) const
    {
        double a = abs(a_long(left_handed, s) * conj(a_par(left_handed, s)) + conj(a_long(right_handed, s)) * a_par(right_handed, s));
        double b = std::sqrt(norm(a_long(left_handed, s)) + norm(a_long(right_handed, s)))
                * std::sqrt(norm(a_perp(left_handed, s)) + norm(a_perp(right_handed, s)));

        return a / b;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_4(const double & s) const
    {
        double a = abs(a_long(left_handed, s) * conj(a_perp(left_handed, s)) - conj(a_long(right_handed, s)) * a_perp(right_handed, s));
        double b = abs(conj(a_long(left_handed, s)) * a_par(left_handed, s) + a_long(right_handed, s) * conj(a_par(right_handed, s)));

        return a / b;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_5(const double & s) const
    {
        double a = abs(a_long(left_handed, s) * conj(a_perp(left_handed, s)) - conj(a_long(right_handed, s)) * a_perp(right_handed, s));
        double b = std::sqrt(norm(a_long(left_handed, s)) + norm(a_long(right_handed, s)))
                * std::sqrt(norm(a_perp(left_handed, s)) + norm(a_perp(right_handed, s)));

        return a / b;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_longitudinal_polarisation(const double & s) const
    {
        return (norm(a_long(left_handed, s)) + norm(a_long(right_handed, s))) * _imp->norm(s) * _imp->norm(s) / differential_decay_width(s);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        std::tr1::function<double (const double &)> f = std::tr1::bind(std::tr1::mem_fn(&BToKstarDilepton<LargeRecoil>::differential_branching_ratio),
                this, std::tr1::placeholders::_1);

        return integrate(f, 128, s_min, s_max);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const
    {
        std::tr1::function<double (const double &)> num = std::tr1::bind(
                std::tr1::mem_fn(&BToKstarDilepton<LargeRecoil>::differential_unnormalized_forward_backward_asymmetry), this, std::tr1::placeholders::_1);

        std::tr1::function<double (const double &)> denom = std::tr1::bind(
                std::tr1::mem_fn(&BToKstarDilepton<LargeRecoil>::differential_decay_width), this, std::tr1::placeholders::_1);

        return integrate(num, 128, s_min, s_max) / integrate(denom, 128, s_min, s_max);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_longitudinal_polarisation(const double & s_min, const double & s_max) const
    {
        std::tr1::function<double (const double &)> num = std::tr1::bind(
                std::tr1::mem_fn(&BToKstarDilepton<LargeRecoil>::differential_longitudinal_polarisation), this, std::tr1::placeholders::_1);

        std::tr1::function<double (const double &)> denom = std::tr1::bind(
                std::tr1::mem_fn(&BToKstarDilepton<LargeRecoil>::differential_decay_width), this, std::tr1::placeholders::_1);

        return integrate(num, 128, s_min, s_max) / integrate(denom, 128, s_min, s_max);
    }

    // Low Recoil

    template <>
    struct Implementation<BToKstarDilepton<LowRecoil>>
    {
        Parameter c1;

        Parameter c2;

        Parameter c3;

        Parameter c4;

        Parameter c5;

        Parameter c6;

        Parameter c7;

        Parameter c7prime;

        Parameter c8;

        Parameter c9;

        Parameter c9prime;

        Parameter c10;

        Parameter c10prime;

        Parameter m_b_MSbar;

        Parameter m_c;

        Parameter m_B;

        Parameter m_Kstar;

        Parameter mu;

        Parameter ckm_A;

        Parameter ckm_lambda;

        Parameter uncertainty_par_left;

        Parameter uncertainty_par_right;

        Parameter uncertainty_perp_left;

        Parameter uncertainty_perp_right;

        Parameter uncertainty_long_left;

        Parameter uncertainty_long_right;

        std::tr1::shared_ptr<FormFactors<BToKstar>> form_factors;

        Implementation(const Parameters & p, const ObservableOptions & o) :
            c1(p["c1"]),
            c2(p["c2"]),
            c3(p["c3"]),
            c4(p["c4"]),
            c5(p["c5"]),
            c6(p["c6"]),
            c7(p["c7"]),
            c8(p["c8"]),
            c9(p["c9"]),
            c10(p["c10"]),
            c7prime(p["c7prime"]),
            c9prime(p["c9prime"]),
            c10prime(p["c10prime"]),
            m_b_MSbar(p["mass::b(MSbar)"]),
            m_c(p["mass::c"]),
            m_B(p["mass::B0"]),
            m_Kstar(p["mass::K^*0"]),
            mu(p["mu"]),
            ckm_A(p["CKM::A"]),
            ckm_lambda(p["CKM::lambda"]),
            uncertainty_par_left(p["B->K^*ll::A_par^L_uncertainty"]),
            uncertainty_par_right(p["B->K^*ll::A_par^R_uncertainty"]),
            uncertainty_perp_left(p["B->K^*ll::A_perp^L_uncertainty"]),
            uncertainty_perp_right(p["B->K^*ll::A_perp^R_uncertainty"]),
            uncertainty_long_left(p["B->K^*ll::A_0^L_uncertainty"]),
            uncertainty_long_right(p["B->K^*ll::A_0^R_uncertainty"])
        {
            form_factors = FormFactorFactory<BToKstar>::create(o["form-factors"], p);
            if (! form_factors.get())
                throw std::string("InternalError");
        }

        // We use the PS mass except for kappa_1
        double m_b_PS() const
        {
            return QCD::mb_PS(m_b_MSbar, mu, 2.0);
        }

        // cf. [BFS2001], Eq. (29), p. 8
        complex<double> B0(const double & s, const double & m_q) const
        {
            double z = 4.0 * m_q * m_q / s;
            double rp, ip;

            if (m_q < 0.0)
            {
                throw InternalError("Implementation<BToKstarDilepton<LargeRecoil>>::B0: m_q < 0!");
            };

            if ((0.0 == m_q) && (0.0 == s))
            {
                throw InternalError("Implementation<BToKstarDilepton<LargeRecoil>>::B0: m_q == 0 & s == 0");
            }

            if (0 == s)
                return complex<double>(-2.0, 0.0);

            if (z > 1.0)
            {
                rp = -2.0 * std::sqrt(z - 1.0) * std::atan(1.0 / std::sqrt(z - 1.0));
                ip = 0.0;
            }
            else
            {
                rp = std::sqrt(1.0 - z) * std::log((1.0 - std::sqrt(1 - z)) / (1.0 + std::sqrt(1.0 - z)));
                ip = std::sqrt(1.0 - z) * M_PI;
            }

            return complex<double>(rp, ip);
        }

        // cf. [BFS2001], Eq. (84), p. 30
        double C0(const double & s) const
        {
            static const int points = 40;
            // Integration boundaries of u = (0, u_max]
            static const double dx = (1.0 - 0.0) / points;
            static const double g[2] =
            {
                (1.0 + std::sqrt(3.0 / 5.0)) / 2.0,
                (1.0 - std::sqrt(3.0 / 5.0)) / 2.0
            };

            long double s_hat = s / m_B / m_B;
            long double x1, x2, x;
            long double y[3];
            long double result = 0;

            for (int i = 0; i < points; i++)
            {
                // 0 is lower integration boundary

                x1 = 0.0 + i * dx;
                x2 = 0.0 + (i + 1) * dx;
                for (int j = 0; j < 3; j++)
                {
                    switch (j)
                    {
                        case 0: x = x1 * g[0] + x2 * g[1];
                                break;

                        case 1: x = (x1 + x2) / 2.0;
                                break;

                        case 2: x = x1 * g[1] + x2 * g[0];
                                break;
                    }

                    y[j] = std::log(x * x / (1.0 - s_hat * x * (1.0 - x))) / (1.0 + x * (1.0 - s_hat));
                }

                result += (5.0 * y[0] + 8.0 * y[1] + 5.0 * y[2]) * dx / 18.0;
            };

            return result;
        }

        // cf. [BFS2001], Eq. (82), p. 30
        complex<double> F87(const double & s) const
        {
            double m_b = m_b_PS();
            double s_hat = s / m_b / m_b;
            double s_hat2 = s_hat * s_hat;
            double denom = (1.0 - s_hat);
            double denom2 = denom * denom;

            complex<double> a = complex<double>(-32.0 * std::log(mu / m_b) - 8.0 * s_hat / denom * std::log(s_hat)
                    - 4.0 * (11.0 - 16.0 * s_hat + 8.0 * s_hat2) / denom2,
                    -8.0 * M_PI);
            complex<double> b = (4.0 / denom / denom2)
                * ((9.0 * s_hat - 5.0 * s_hat2 + 2.0 * s_hat * s_hat2) * B0(s, m_b) - (4.0 + 2.0 * s_hat) * C0(s));

            return (1.0 / 9.0) * (a + b);
        }

        // cf. [BFS2001], Eq. (83), p. 30
        complex<double> F89(const double & s) const
        {
            double m_b = m_b_PS();
            double s_hat = s / m_b / m_b;
            double denom = (1.0 - s_hat);
            double denom2 = denom * denom;

            double a = 16.0 * std::log(s_hat) / denom + 8.0 * (5.0 - 2.0 * s_hat) / denom2;
            complex<double> b = (-8.0 * (4.0 - s_hat) / denom / denom2) * ((1.0 + s_hat) * B0(s, m_b) - 2.0 * C0(s));

            return (1.0 / 9.0) * (a + b);
        }

        // cf. [GP2004], Eq. (56)
        complex<double> c7eff(double s) const
        {
            double m_b = m_b_PS();

            // cf. [BFS2001] Eq. (29), p. 8, and Eqs. (82)-(84), p. 30
            double lo = - 1.0/3.0 * c3 - 4.0/9.0 * c4 + 20.0/3.0 * c5 + 80.0/9.0 * c6;
            complex<double> nlo = -1.0 * (c1() * CharmLoops::F17(mu, s, m_b) + c2() * CharmLoops::F27(mu, s, m_b) + c8() * F87(s));

            return c7() + lo + (QCD::alpha_s(mu) / (4.0 * M_PI)) * nlo;
        }

        // cf. [GP2004], Eq. (55), p. 10
        complex<double> c9eff(const double & s) const
        {
            double m_b = m_b_PS();

            double c_0 = 4.0 / 3.0 * c1() + c2() + 5.5 * c3() - 2.0 / 3.0 * c4() + 52.0 * c5() + 32.0 / 3.0 * c6();
            double c_b = -0.5 * (7.0 * c3() + 4.0 / 3.0 * c4() + 76.0 * c5() + 64.0 / 3.0 * c6());
            double c = 2.0 / 9.0 * (6.0 * c3() + 32.0 * c5() + 32.0 / 3.0 * c6());

            // Uses b pole mass according to [BFS2001], Sec. 3.1, paragraph Quark Masses
            complex<double> lo = c_b * CharmLoops::h(mu, s, m_b) + c_0 * CharmLoops::h(mu, s) + c;
            complex<double> nlo = -1.0 * (c1() * CharmLoops::F19(mu, s, m_b) + c2() * CharmLoops::F29(mu, s, m_b) + c8() * F89(s));

            return c9() + lo + (QCD::alpha_s(mu) / (4.0 * M_PI)) * nlo;
        }

        double kappa_1() const
        {
            // cf. [BHvD2010], Eq. (?), p. ?
            // Use m_b_MSbar(m_b_MSbar) instead m_b_MSbar(mu), as we want kappa_1 up to NLO only.
            return (1.0 - 2.0 * QCD::alpha_s(mu) / (3.0 * M_PI) * std::log(mu / m_b_MSbar));
        }

        double norm(const double & s) const
        {
            static const double alpha_e = 1.0 / 133.0; // cf. [BHP2008]
            static const double g_fermi = 1.16637e-5; // (Gev^-2 (hbar c)^3), cf. [PDG2008], p.5
            double lambda_t = ckm_A * ckm_lambda * ckm_lambda;

            return std::sqrt(g_fermi * g_fermi * alpha_e * alpha_e / 3.0 / 1024 / std::pow(M_PI, 5.0) / m_B
                    * lambda_t * lambda_t * s_hat(s)
                    * std::sqrt(lambda(m_B * m_B, m_Kstar * m_Kstar, s))); // cf. [BHP2008], Eq. (C.6), p. 21
        }

        inline double s_hat(const double & s) const
        {
            return s / m_B / m_B;
        }

        // Amplitudes
        complex<double> a_long(const Helicity & helicity, const double & s) const
        {
            double m_b = m_b_PS();
            double h = helicity;
            double m_Kstar_hat = m_Kstar / m_B;
            complex<double> wilson = (c9eff(s) - c9prime()) + h * (c10() - c10prime()) + kappa_1() * (c7eff(s) + c7prime()) * (2 * m_b * m_B / s);
            complex<double> prefactor = complex<double>(0.0, -0.5 * m_B * m_B * m_B / m_Kstar / std::sqrt(s));
            double uncertainty = (1.0 - h) / 2.0 * uncertainty_long_left + (1.0 + h) / 2.0 * uncertainty_long_right;
            double formfactor = lambda(1.0, m_Kstar_hat * m_Kstar_hat, s_hat(s)) * form_factors->a_1(s_hat(s)) - (1 - s_hat(s)) * form_factors->a_2(s_hat(s));

            return uncertainty * prefactor * wilson * formfactor; // cf. [BHvD2010], Eq. (??)
        }

        complex<double> a_perp(const Helicity & helicity, const double & s) const
        {
            double m_b = m_b_PS();
            double h = helicity;
            double m_Kstar_hat = m_Kstar / m_B;
            complex<double> wilson = (c9eff(s) + c9prime()) + h * (c10() + c10prime()) + kappa_1() * (c7eff(s) - c7prime()) * (2 * m_b * m_B / s);
            complex<double> prefactor = complex<double>(0.0, std::sqrt(2 * lambda(1.0, m_Kstar_hat * m_Kstar_hat, s_hat(s))) * m_B);
            double uncertainty = (1.0 - h) / 2.0 * uncertainty_perp_left + (1.0 + h) / 2.0 * uncertainty_perp_right;

            return uncertainty * prefactor * wilson * form_factors->v(s_hat(s)); // cf. [BHvD2010], Eq. (??)
        }

        complex<double> a_par(const Helicity & helicity, const double & s) const
        {
            double m_b = m_b_PS();
            double h = helicity;
            complex<double> wilson = (c9eff(s) - c9prime()) + h * (c10() - c10prime()) + kappa_1() * (c7eff(s) - c7prime()) * (2 * m_b * m_B / s);
            complex<double> prefactor = complex<double>(0.0, -std::sqrt(2) * m_B);
            double uncertainty = (1.0 - h) / 2.0 * uncertainty_par_left + (1.0 + h) / 2.0 * uncertainty_par_right;

            return prefactor * wilson * form_factors->a_1(s_hat(s)); // cf. [BHvD2010], Eq. (??)
        }
    };

    BToKstarDilepton<LowRecoil>::BToKstarDilepton(const Parameters & parameters, const ObservableOptions & options) :
        PrivateImplementationPattern<BToKstarDilepton<LowRecoil>>(new Implementation<BToKstarDilepton<LowRecoil>>(parameters, options))
    {
    }

    BToKstarDilepton<LowRecoil>::~BToKstarDilepton()
    {
    }

    complex<double>
    BToKstarDilepton<LowRecoil>::a_long(const Helicity & h, const double & s) const
    {
        return _imp->a_long(h, s);
    }

    complex<double>
    BToKstarDilepton<LowRecoil>::a_perp(const Helicity & h, const double & s) const
    {
        return _imp->a_perp(h, s);
    }

    complex<double>
    BToKstarDilepton<LowRecoil>::a_par(const Helicity & h, const double & s) const
    {
        return _imp->a_par(h, s);
    }

    double
    BToKstarDilepton<LowRecoil>::differential_branching_ratio(const double & s) const
    {
        // cf. [PDG2008] : Gamma = hbar / tau_B, pp. 5, 79
        static const double Gamma(6.58211899e-22 * 1e-3 / 1.53e-12);

        return differential_decay_width(s) * (_imp->norm(s) * _imp->norm(s) / Gamma);
    }

    double
    BToKstarDilepton<LowRecoil>::differential_decay_width(const double & s) const
    {
        return (norm(a_long(left_handed, s))
            + norm(a_long(right_handed, s))
            + norm(a_perp(left_handed, s))
            + norm(a_perp(right_handed, s))
            + norm(a_par(left_handed, s))
            + norm(a_par(right_handed, s)));
    }

    double
    BToKstarDilepton<LowRecoil>::differential_forward_backward_asymmetry(const double & s) const
    {
        return 1.5 / differential_decay_width(s) * (real(a_par(left_handed, s) * conj(a_perp(left_handed, s))) - real(a_par(right_handed, s) * conj(a_perp(right_handed, s))));
    }

    double
    BToKstarDilepton<LowRecoil>::differential_unnormalized_forward_backward_asymmetry(const double & s) const
    {
        return 1.5 * (real(a_par(left_handed, s) * conj(a_perp(left_handed, s))) - real(a_par(right_handed, s) * conj(a_perp(right_handed, s))));
    }

    double
    BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_2(const double & s) const
    {
        double a = norm(a_perp(left_handed, s)) + norm(a_perp(right_handed, s));
        double b = norm(a_par(left_handed, s)) + norm(a_par(right_handed, s));

        return (a - b) / (a + b);
    }

    double
    BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_3(const double & s) const
    {
        double a = abs(a_long(left_handed, s) * conj(a_par(left_handed, s)) + conj(a_long(right_handed, s)) * a_par(right_handed, s));
        double b = std::sqrt((norm(a_long(left_handed, s)) + norm(a_long(right_handed, s))))
                * std::sqrt((norm(a_perp(left_handed, s)) + norm(a_perp(right_handed, s))));

        return a / b;
    }

    double
    BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_4(const double & s) const
    {
        double a = abs(a_long(left_handed, s) * conj(a_perp(left_handed, s)) - conj(a_long(right_handed, s)) * a_perp(right_handed, s));
        double b = abs(conj(a_long(left_handed, s)) * a_par(left_handed, s) + a_long(right_handed, s) * conj(a_par(right_handed, s)));

        return a / b;
    }

    double
    BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_5(const double & s) const
    {
        double a = abs(a_long(left_handed, s) * conj(a_perp(left_handed, s)) - conj(a_long(right_handed, s)) * a_perp(right_handed, s));
        double b = std::sqrt((norm(a_long(left_handed, s)) + norm(a_long(right_handed, s))))
                * std::sqrt((norm(a_perp(left_handed, s)) + norm(a_perp(right_handed, s))));

        return a / b;
    }

    double
    BToKstarDilepton<LowRecoil>::differential_longitudinal_polarisation(const double & s) const
    {
        return (norm(a_long(left_handed, s)) + norm(a_long(right_handed, s)))
            / differential_decay_width(s);
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        std::tr1::function<double (const double &)> f = std::tr1::bind(std::tr1::mem_fn(&BToKstarDilepton<LowRecoil>::differential_branching_ratio),
                this, std::tr1::placeholders::_1);

        return integrate(f, 100, s_min, s_max);
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const
    {
        std::tr1::function<double (const double &)> num = std::tr1::bind(
                std::tr1::mem_fn(&BToKstarDilepton<LowRecoil>::differential_unnormalized_forward_backward_asymmetry), this, std::tr1::placeholders::_1);

        std::tr1::function<double (const double &)> denom = std::tr1::bind(
                std::tr1::mem_fn(&BToKstarDilepton<LowRecoil>::differential_decay_width), this, std::tr1::placeholders::_1);

        return integrate(num, 1000, s_min, s_max) / integrate(denom, 1000, s_min, s_max);
    }
}
