#ifndef _screen_H_
#define _screen_H_

#include <AMReX.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>
#include <network_properties.H>
#include <microphysics_math.H>
#include <microphysics_autodiff.H>
#include <fundamental_constants.H>
#include <cmath>
#include <screen_data.H>
#include <extern_parameters.H>

using namespace amrex;
namespace admath = autodiff::math_functions;

#if SCREEN_METHOD == 0
const std::string screen_name = "screen5";
#elif SCREEN_METHOD == 1
const std::string screen_name = "chugunov2007";
#elif SCREEN_METHOD == 2
const std::string screen_name = "chugunov2009";
#elif SCREEN_METHOD == 3
const std::string screen_name = "null";
#elif SCREEN_METHOD == 4
const std::string screen_name = "chabrier1998";
#endif

struct plasma_state_t {

  autodiff::dual qlam0z;

  autodiff::dual taufac;

  autodiff::dual aa;

  autodiff::dual temp;
  Real zbar;
  Real z2bar;
  Real n_e;
  Real gamma_e_fac;
};

inline
std::ostream& operator<< (std::ostream& o, plasma_state_t const& pstate)
{
  o << "qlam0z   = " << autodiff::repr(pstate.qlam0z) << std::endl;
  o << "taufac   = " << autodiff::repr(pstate.taufac) << std::endl;
  o << "aa       = " << autodiff::repr(pstate.aa) << std::endl;
  o << "temp     = " << autodiff::repr(pstate.temp) << std::endl;
  o << "zbar     = " << pstate.zbar << std::endl;
  o << "z2bar    = " << pstate.z2bar << std::endl;
  o << "n_e      = " << pstate.n_e << std::endl;
  o << "gamma_e_fac = " << pstate.gamma_e_fac << std::endl;

  return o;
}

AMREX_FORCE_INLINE
void
screening_init() {

}

AMREX_FORCE_INLINE
void
screening_finalize() {

}


AMREX_FORCE_INLINE
void add_screening_factor(const int i,
                          const Real z1, const Real a1, const Real z2, const Real a2)
{
  using namespace scrn;

  BL_ASSERT(i < NSCREEN);

  scn_facs[i].z1 = z1;
  scn_facs[i].z2 = z2;
  scn_facs[i].a1 = a1;
  scn_facs[i].a2 = a2;

  scn_facs[i].zs13 = std::cbrt(z1 + z2);
  scn_facs[i].zs13inv = 1.0_rt/scn_facs[i].zs13;
  scn_facs[i].zhat = std::pow(z1 + z2, 5.0_rt/3.0_rt) -
                     std::pow(z1, 5.0_rt/3.0_rt) -
                     std::pow(z2, 5.0_rt/3.0_rt);
  scn_facs[i].zhat2 = std::pow(z1 + z2, 5.0_rt/12.0_rt) -
                      std::pow(z1, 5.0_rt/12.0_rt) -
                      std::pow(z2, 5.0_rt/12.0_rt);
  scn_facs[i].lzav = (5.0_rt/3.0_rt) * std::log(z1 * z2 / (z1 + z2));
  scn_facs[i].aznut = std::cbrt(z1 * z1 * z2 * z2 * a1 * a2 / (a1 + a2));
  scn_facs[i].ztilde = 0.5_rt * (std::cbrt(z1) + std::cbrt(z2));

}


template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
fill_plasma_state(plasma_state_t& state, const dual_t temp, const Real dens, Array1D<Real, 1, NumSpec> const& y) {

  Real sum = 0.0_rt;
  for (int n = 1; n <= NumSpec; n++) {
    sum += y(n);
  }
  Real abar = 1.0_rt / sum;
  Real ytot = sum;

  sum = 0.0_rt;
  Real sum2 = 0.0_rt;
  for (int n = 1; n <= NumSpec; n++) {
    sum += zion[n-1]*y(n);
    sum2 += zion[n-1]*zion[n-1]*y(n);
  }

  // Part of Eq.6 in Itoh:1979
  // 4.248719e3 = (27*pi^2*e^4*m_u/(2*k_B*hbar^2))^(1/3)
  // the extra (1/3) to make tau -> tau/3
  const Real co2 = (1.0_rt/3.0_rt) * 4.248719e3_rt;

  Real zbar = sum * abar;
  Real z2bar = sum2 * abar;

  // ntot
  Real rr = dens * ytot;
  dual_t tempi = 1.0_rt / temp;

  // Part of Eq. 19 in Graboske:1973
  // pp = sqrt( \tilde{z}*(rho/u_I/T) )
  dual_t pp = admath::sqrt(rr*tempi*(z2bar + zbar));

  // Part version of Eq. 19 in Graboske:1973
  state.qlam0z = 1.88e8_rt * tempi * pp;

  // Part of Eq.6 in Itoh:1979
  state.taufac = co2 * admath::cbrt(tempi);

  Real xni = admath::cbrt(rr * zbar);

  // Part of Eq.4 in Itoh:1979
  // 2.27493e5 = e^2 / ( (3*m_u/(4pi))^(1/3) *k_B )
  state.aa = 2.27493e5_rt * tempi * xni;

  state.temp = temp;
  state.zbar = zbar;
  state.z2bar = z2bar;

  // Electron number density
  // zbar * ntot works out to sum(z[i] * n[i]), after cancelling terms
  state.n_e = zbar * rr * C::n_A;

  // precomputed part of Gamma_e, from Chugunov 2009 eq. 6
  constexpr Real gamma_e_constants =
      C::q_e*C::q_e/C::k_B * gcem::pow(4.0_rt/3.0_rt*M_PI, 1.0_rt/3.0_rt);
  state.gamma_e_fac = gamma_e_constants * std::cbrt(state.n_e);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
fill_plasma_state(plasma_state_t& state, const autodiff::dual temp, const Real dens, Array1D<Real, 1, NumSpec> const& y) {
    // temp will be automatically upcast from a Real if necessary
    if (jacobian == 1) {
        fill_plasma_state<autodiff::dual>(state, temp, dens, y);
    } else {
        fill_plasma_state<amrex::Real>(state, autodiff::val(temp), dens, y);
    }
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
dual_t actual_screen5 (const plasma_state_t& state,
                       const scrn::screen_factors_t& scn_fac)
{
    // this subroutine calculates screening factors and their derivatives
    // for nuclear reaction rates in the weak, intermediate and strong regimes.
    // based on graboske, dewit, grossman and cooper apj 181 457 1973 for
    // weak screening. based on alastuey and jancovici apj 226 1034 1978,
    // with plasma parameters from itoh et al apj 234 1079 1979, for strong
    // screening.

    // input:
    // state   = plasma state (T, rho, abar, zbar, etc.)
    // scn_fac = screening factors for A and Z


    // fact = 2^(1/3)
    const Real fact    = 1.25992104989487e0_rt;
    const Real gamefx  = 0.3e0_rt;          // lower gamma limit for intermediate screening
    const Real gamefs  = 0.8e0_rt;          // upper gamma limit for intermediate screening
    const Real h12_max = 300.e0_rt;

    // Get the ion data based on the input index
    Real z1 = scn_fac.z1;
    Real z2 = scn_fac.z2;

    // calculate individual screening factors
    Real bb = z1 * z2;
    dual_t gamp = static_cast<dual_t>(state.aa);

    // In Eq.4 in Itoh:1979, this term is 2*Z_1*Z_2/(Z_1^(1/3) + Z_2^(1/3))
    // However here we follow Wallace:1982 Eq. A13, which is Z_1*Z_2*(2/(Z_1+Z_2))^(1/3)
    
    Real qq = fact * bb * scn_fac.zs13inv;

    // Full Equation of Wallace:1982 Eq. A13
    
    dual_t gamef = qq * gamp;

    // Full version of Eq.6 in Itoh:1979 with extra 1/3 factor
    // the extra 1/3 factor is there for convenience.
    // tau12 = Eq.6 / 3
    
    dual_t tau12 = static_cast<dual_t>(state.taufac) * scn_fac.aznut;

    // alph12 = 3*gamma_ij/tau_ij
    
    dual_t alph12 = gamef / tau12;


    // limit alph12 to 1.6 to prevent unphysical behavior.
    // See Introduction in Alastuey:1978

    // this should really be replaced by a pycnonuclear reaction rate formula
    if (alph12 > 1.6_rt) {
        alph12   = 1.6e0_rt;

	// redetermine previous factors if 3*gamma_ij/tau_ij > 1.6
	
        gamef    = 1.6e0_rt * tau12;

        qq = scn_fac.zs13/(fact * bb);
        gamp = gamef * qq;
    }

    // weak screening regime
    // Full version of Eq. 19 in Graboske:1973 by considering weak regime
    // and Wallace:1982 Eq. A14. Here the degeneracy factor is assumed to be 1.
    
    dual_t h12w = bb * static_cast<dual_t>(state.qlam0z);

    dual_t h12 = h12w;

    // intermediate and strong sceening regime
    
    if (gamef > gamefx) {

      // gamma_ij^(1/4)
      
        dual_t gamp14 = admath::pow(gamp, 0.25_rt);
        dual_t rr = 1.0_rt/gamp;

	// Here we follow Eq. A9 in Wallace:1982
	// See Eq. 25 Alastuey:1978, Eq. 16 and 17 in Jancovici:1977 for reference
        dual_t cc = 0.896434e0_rt * gamp * scn_fac.zhat
            - 3.44740e0_rt * gamp14 * scn_fac.zhat2
            - 0.5551e0_rt * (admath::log(gamp) + scn_fac.lzav)
            - 2.996e0_rt;

	// (3gamma_ij/tau_ij)^3
        dual_t a3 = alph12 * alph12 * alph12;

	// Part of Eq. 28 in Alastuey:1978
        dual_t tmp = 0.014e0_rt + 0.0128e0_rt*alph12;

	// Part of Eq. 28 in Alastuey:1978
        rr = (5.0_rt/32.0_rt) - alph12*tmp;

	// Part of Eq. 28 in Alastuey:1978
        dual_t ss = tau12*rr;

	// Part of Eq. 31 in Alastuey:1978
        dual_t tt = -0.0098e0_rt + 0.0048e0_rt*alph12;

	// Part of Eq. 31 in Alastuey:1978
        dual_t uu = 0.0055e0_rt + alph12*tt;

	// Part of Eq. 31 in Alastuey:1978
        dual_t vv = gamef * alph12 * uu;

	// Exponent of Eq. 32 in Alastuey:1978, which uses Eq.28 and Eq.31
	// Strong screening factor
        h12 = cc - a3 * (ss + vv);

	// See conclusion and Eq. 34 in Alastuey:1978
	// This is an extra factor to account for quantum effects
        rr = 1.0_rt - 0.0562e0_rt*a3;

        dual_t xlgfac;

	// In extreme case, rr is 0.77, see conclusion in Alastuey:1978
        if (rr >= 0.77e0_rt) {
            xlgfac = rr;
        } else {
            xlgfac = 0.77e0_rt;
        }

	// Include the extra factor that accounts for quantum effects
        h12 = admath::log(xlgfac) + h12;
        rr = 1.0_rt/xlgfac;

	// If gamma_ij < upper limit of intermediate regime
	// then it is in the intermediate regime, else strong screening.
        if (gamef <= gamefs) {
            Real dgamma = 1.0e0_rt/(gamefs - gamefx);

            rr =  dgamma*(gamefs - gamef);

            ss = dgamma*(gamef - gamefx);

            vv = h12;

	    // Then the screening factor is a combination
	    // of the strong and weak screening factor.
            h12 = h12w*rr + vv*ss;
        }

        // end of intermediate and strong screening
    }

    // machine limit the output
    // further limit to avoid the pycnonuclear regime
    h12 = admath::max(admath::min(h12, h12_max), 0.0_rt);
    return admath::exp(h12);
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
dual_t chugunov2007 (const plasma_state_t& state,
                     const scrn::screen_factors_t& scn_fac)
{
    // Calculates screening factors based on Chugunov et al. 2007, following the
    // the approach in Yakovlev 2006 to extend to a multi-component plasma.
    //
    // References:
    //   Chugunov, DeWitt, and Yakovlev 2007, PhRvD, 76, 025028
    //   Yakovlev, Gasques, Afanasjev, Beard, and Wiescher 2006, PhRvC, 74, 035803
    //   Chugunov and DeWitt 2009, PhRvC, 80, 014611

    // input:
    // state   = plasma state (T, rho, abar, zbar, etc.)
    // scn_fac = screening factors for A and Z

    // Plasma temperature T_p
    // This formula comes from working backwards from zeta_ij (Chugunov 2009 eq. 12)
    // through Chugunov 2007 eq. 3 to Chugunov 2007 eq. 2.
    // Ultimately, the changes from the expression in Chugunov 2007 are:
    //   Z^2 -> Z1 * Z2
    //   n_i -> n_e / ztilde^3, where ztilde = (Z1^(1/3) + Z2^(1/3)) / 2
    //   m_i -> 2 mu12 (reduced mass)
    // This prescription reduces to the expressions from Chugunov 2007 in the case
    // of an OCP, and to Chugunov 2009 in the case of a binary ionic mixture.
    // This also matches Yakovlev et al. 2006, eq. 10.
    //
    // For reference, MESA r21.12.1 does:
    //   Z^2 -> Z1 * Z2
    //   n_i -> n_e / zbar (=ntot)
    //   m_i -> m_u * abar
    // Sam Jones' Fortran implementation (linked in issue #146) does:
    //   Z^2 -> zbar^2
    //   n_i -> ntot
    //   m_i -> m_u * abar
    Real mu12 = scn_fac.a1 * scn_fac.a2 / (scn_fac.a1 + scn_fac.a2);
    Real z_factor = scn_fac.z1 * scn_fac.z2;
    Real n_i = state.n_e / gcem::pow(scn_fac.ztilde, 3);
    Real m_i = 2.0_rt * mu12 / C::n_A;

    constexpr Real T_p_factor = C::hbar/C::k_B*C::q_e*gcem::sqrt(4.0_rt*GCEM_PI);
    Real T_p = T_p_factor * std::sqrt(z_factor * n_i / m_i);

    // Normalized temperature
    Real inv_T_p = 1.0_rt / T_p;
    dual_t T_norm = static_cast<dual_t>(state.temp) * inv_T_p;

    // The fit has only been verified down to T ~ 0.1 T_p, below which the rate
    // should be nearly temperature-independent (in the pycnonuclear regime),
    // and we clip the temperature to 0.1 T_p at small T.
    // start the transition here
    constexpr Real T_norm_fade = 0.2_rt;
    constexpr Real T_norm_min = 0.1_rt;

    if (T_norm < T_norm_min) {
        // clip temperature to the minimum value
        T_norm = T_norm_min;
    } else if (T_norm <= T_norm_fade) {
        // blend using a cosine, after MESA
        constexpr Real delta_T = T_norm_fade - T_norm_min;
        dual_t tmp = M_PI * (T_norm - T_norm_min) / delta_T;
        dual_t f = 0.5_rt * (1.0_rt - admath::cos(tmp));
        T_norm = (1.0_rt - f) * T_norm_min + f * T_norm;
    }

    // Coulomb coupling parameter from Yakovlev 2006 eq. 10
    dual_t Gamma = state.gamma_e_fac*scn_fac.z1*scn_fac.z2 / (scn_fac.ztilde*T_norm*T_p);

    // The fit for Gamma is only applicable up to ~600, so smoothly cap its value
    constexpr Real Gamma_fade = 590;
    constexpr Real Gamma_max = 600;

    if (Gamma > Gamma_max) {
        // clip Gamma to the max value
        Gamma = Gamma_max;
    } else if (Gamma >= Gamma_fade) {
        // blend using a cosine, after MESA
        constexpr Real delta_gamma = Gamma_max - Gamma_fade;
        dual_t tmp = M_PI * (Gamma - Gamma_fade) / delta_gamma;
        dual_t f = 0.5_rt * (1.0_rt - admath::cos(tmp));
        Gamma = (1.0_rt - f) * Gamma + f * Gamma_max;
    }

    // Chugunov 2007 eq. 3
    dual_t zeta = admath::cbrt(4.0_rt / (3.0_rt * GCEM_PI*GCEM_PI * T_norm*T_norm));

    // Gamma tilde from Chugunov 2007 eq. 21
    constexpr Real fit_alpha = 0.022_rt;
    dual_t fit_beta = 0.41_rt - 0.6_rt / Gamma;
    dual_t fit_gamma = 0.06_rt + 2.2_rt / Gamma;
    // Polynomial term in Gamma tilde
    dual_t poly = 1.0_rt + zeta*(fit_alpha + zeta*(fit_beta + fit_gamma*zeta));

    dual_t gamtilde = Gamma / admath::cbrt(poly);

    // fit parameters just after Chugunov 2007 eq. 19
    constexpr Real A1 = 2.7822_rt;
    constexpr Real A2 = 98.34_rt;
    constexpr Real A3 = gcem::sqrt(3.0_rt) - A1 / gcem::sqrt(A2);
    const Real B1 = -1.7476_rt;
    const Real B2 = 66.07_rt;
    const Real B3 = 1.12_rt;
    const Real B4 = 65_rt;
    dual_t gamtilde2 = gamtilde * gamtilde;

    dual_t term1, term2, term3, term4;
    // Chugunov 2007 eq. 19
    term1 = 1.0_rt / admath::sqrt(A2 + gamtilde);

    term2 = 1.0_rt / (1.0_rt + gamtilde);

    dual_t tmp = B2 + gamtilde;
    term3 = gamtilde2 / tmp;

    tmp = B4 + gamtilde2;
    term4 = gamtilde2 / tmp;

    dual_t inner = A1 * term1 + A3 * term2;

    dual_t gamtilde32 = admath::pow(gamtilde, 1.5_rt);
    dual_t h = gamtilde32 * inner + B1 * term3 + B3 * term4;

    // machine limit the output
    constexpr Real h_max = 300.e0_rt;
    h = admath::min(h, h_max);
    return admath::exp(h);
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
dual_t f0 (const dual_t& gamma)
{
    // Calculates the free energy per ion in a OCP, from Chugunov and DeWitt 2009
    // equation 24.

    // fit parameters
    constexpr Real A1 = -0.907_rt;
    constexpr Real A2 = 0.62954_rt;
    constexpr Real A3 = -gcem::sqrt(3.0_rt) / 2.0_rt - A1 / gcem::sqrt(A2);
    constexpr Real B1 = 0.00456_rt;
    constexpr Real B2 = 211.6_rt;
    constexpr Real B3 = -1e-4_rt;
    constexpr Real B4 = 0.00462_rt;
    dual_t gamma_12 = admath::sqrt(gamma);

    dual_t term1 = gamma_12 * admath::sqrt(A2 + gamma);

    dual_t term2 = admath::log(admath::sqrt(gamma / A2) + admath::sqrt(1.0_rt + gamma / A2));

    dual_t term3 = gamma_12 - admath::atan(gamma_12);

    dual_t term4 = admath::log1p(gamma / B2);

    dual_t term5 = admath::log1p(gamma * gamma / B4);

    return A1*(term1 - A2*term2) + 2.0_rt*A3*term3 + B1*(gamma - B2*term4) + 0.5_rt*B3*term5;
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
dual_t chugunov2009 (const plasma_state_t& state,
                     const scrn::screen_factors_t& scn_fac)
{
    // Calculates screening factors based on Chugunov and DeWitt 2009, PhRvC, 80, 014611

    // input:
    // state   = plasma state (T, rho, abar, zbar, etc.)
    // scn_fac = screening factors for A and Z

    Real z1z2 = scn_fac.z1 * scn_fac.z2;
    Real zcomp = scn_fac.z1 + scn_fac.z2;

    // Gamma_e from eq. 6
    dual_t Gamma_e = state.gamma_e_fac / static_cast<dual_t>(state.temp);

    // Coulomb coupling parameters for ions and compound nucleus, eqs. 7 & 9
    dual_t Gamma_1 = Gamma_e * std::pow(scn_fac.z1, 5.0_rt/3.0_rt);
    dual_t Gamma_2 = Gamma_e * std::pow(scn_fac.z2, 5.0_rt/3.0_rt);
    dual_t Gamma_comp = Gamma_e * std::pow(zcomp, 5.0_rt/3.0_rt);

    dual_t Gamma_12 = Gamma_e * z1z2 / scn_fac.ztilde;

    // Coulomb barrier penetrability, eq. 10
    constexpr Real tau_factor = gcem::pow(
            27.0_rt/2.0_rt * gcem::pow(M_PI*C::q_e*C::q_e/C::hbar, 2) / (C::n_A*C::k_B),
            1.0_rt/3.0_rt);
    dual_t tau_12 = tau_factor * scn_fac.aznut / admath::cbrt(static_cast<dual_t>(state.temp));

    // eq. 12
    dual_t zeta = 3.0_rt * Gamma_12 / tau_12;

    // additional fit parameters, eq. 25
    Real y_12 = 4.0_rt * z1z2 / (zcomp * zcomp);
    Real c1 = 0.013_rt * y_12 * y_12;
    Real c2 = 0.406_rt * std::pow(y_12, 0.14_rt);
    dual_t c3 = 0.062_rt * std::pow(y_12, 0.19_rt) + 1.8_rt / Gamma_12;

    dual_t poly = 1.0_rt + zeta*(c1 + zeta*(c2 + c3*zeta));
    dual_t t_12 = admath::cbrt(poly);

    // strong screening enhancement factor, eq. 23, replacing tau_ij with t_ij
    // Using Gamma/tau_ij gives extremely low values, while Gamma/t_ij gives
    // values similar to those from Chugunov 2007.
    dual_t term1 = f0<dual_t>(Gamma_1 / t_12);
    dual_t term2 = f0<dual_t>(Gamma_2 / t_12);
    dual_t term3 = f0<dual_t>(Gamma_comp / t_12);
    dual_t h_fit = term1 + term2 - term3;

    // weak screening correction term, eq. A3
    Real corr_C = 3.0_rt*z1z2 * std::sqrt(state.z2bar/state.zbar) /
                  (std::pow(zcomp, 2.5_rt) - std::pow(scn_fac.z1, 2.5_rt) - std::pow(scn_fac.z2, 2.5_rt));

    // corrected enhancement factor, eq. A4
    dual_t Gamma_12_2 = Gamma_12 * Gamma_12;
    dual_t numer = corr_C + Gamma_12_2;
    dual_t denom = 1.0_rt + Gamma_12_2;
    dual_t h12 = numer / denom * h_fit;

    // machine limit the output
    const Real h12_max = 300.e0_rt;
    h12 = admath::min(h12, h12_max);
    return admath::exp(h12);
}


template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
dual_t chabrier1998 (const plasma_state_t& state,
                     const scrn::screen_factors_t& scn_fac)
{
  // Calculates screening factors based on Chabrier & Potekhin 1998,
  // Calder2007 and partly screen5 routine mentioned in Alastuey 1978.
  
  // This screening is valid for weak screening: Gamma < 0.1
  // and strong screening: 1 <= Gamma <= 160
  // Reference:
  //            Chabrier & Potekhin 1998, PhRvE, 58, 4941
  //            Calder 2007, doi:10.1086/510709
  //            Wallace & Woosley 1982
  //            Alastuey 1978

  // Eq. 2 in Chabrier & Potekhin 1998
  
  dual_t Gamma_e = state.gamma_e_fac / static_cast<dual_t>(state.temp);
  Real zcomp = scn_fac.z1 + scn_fac.z2;

  // See Calder2007 appendix Eq. A9
  
  dual_t Gamma1 = Gamma_e * std::pow(scn_fac.z1, 5.0_rt/3.0_rt);
  dual_t Gamma2 = Gamma_e * std::pow(scn_fac.z2, 5.0_rt/3.0_rt);
  dual_t Gamma12 = Gamma_e * std::pow(zcomp, 5.0_rt/3.0_rt);

  // Fitted parameters, see Chabrier & Potekhin 1998 Sec.IV
  
  const Real A_1 = -0.9052_rt;
  const Real A_2 = 0.6322_rt;
  const Real A_3 = -0.5_rt * std::sqrt(3.0_rt) - A_1 / std::sqrt(A_2);
    
  // Helmholtz free energy, See Chabrier & Potekhin 1998 Eq. 28

  dual_t f1 = A_1 * (admath::sqrt(Gamma1 * (A_2 + Gamma1)) -
		     A_2 * admath::log(admath::sqrt(Gamma1/A_2) + admath::sqrt(1.0_rt + Gamma1/A_2)))
    + 2.0_rt * A_3 * (admath::sqrt(Gamma1) - admath::atan(admath::sqrt(Gamma1)));
  
  dual_t f2 = A_1 * (admath::sqrt(Gamma2 * (A_2 + Gamma2)) -
		     A_2 * admath::log(admath::sqrt(Gamma2/A_2) + admath::sqrt(1.0_rt + Gamma2/A_2)))
    + 2.0_rt * A_3 * (admath::sqrt(Gamma2) - admath::atan(admath::sqrt(Gamma2)));
  
  dual_t f12 = A_1 * (admath::sqrt(Gamma12 * (A_2 + Gamma12)) -
		      A_2 * admath::log(admath::sqrt(Gamma12/A_2) + admath::sqrt(1.0_rt + Gamma12/A_2)))
    + 2.0_rt * A_3 * (admath::sqrt(Gamma12) - admath::atan(admath::sqrt(Gamma12)));

  // Now we add quantum correction terms discussed in Alastuey 1978.
  // Notice in Alastuey 1978, they have a differnt classical term,
  // which is implemented in the strong screening limit of our screen5 routine.

  // See Wallace1982, Eq. A13
  
  dual_t Gamma_eff = std::cbrt(2.0_rt) * scn_fac.z1 * scn_fac.z2 * scn_fac.zs13inv * Gamma_e;
  [[maybe_unused]] Real Gamma_eff_dT;

  // TAU/3, see Wallace1982, Eq. A2
  
  dual_t tau12 = static_cast<dual_t>(state.taufac) * scn_fac.aznut;

  // see Calder 2007 Eq. A8
  
  dual_t b_fac = Gamma_eff / tau12;

  // Quantum correction terms (same as screen5)
  //see Calder 2007 Eq.A8 and Alastuey1978, Eq. 24 and 31
  
  dual_t quantum_corr_1 = -tau12 * (5.0_rt/32.0_rt * admath::pow(b_fac, 3)
				    - 0.014_rt * admath::pow(b_fac, 4)
				    - 0.128_rt * admath::pow(b_fac, 5));
  
  dual_t quantum_corr_2 = -Gamma_eff * (0.0055_rt * admath::pow(b_fac, 4)
				        - 0.0098_rt * admath::pow(b_fac, 5)
				        - 0.0048_rt * admath::pow(b_fac, 6));

  // See Calder2007 Appendix Eq. A8.
  // f1 + f2 - f12 gives the classical terms
  // The difference between this and strong screening of screen5
  // is that we replaced the classical term which is f1 + f2 - f12
  // using results from Chabrier&Potekhin1998.
  
  dual_t h12 = f1 + f2 - f12 + quantum_corr_1 + quantum_corr_2;
  
  Real h12_max = 300.0_rt;
  h12 = admath::min(h12_max, h12);
  
  return admath::exp(h12);
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
dual_t actual_screen(const plasma_state_t& state,
                     const scrn::screen_factors_t& scn_fac)
{
    dual_t scor;
#ifdef AMREX_DEBUG
#endif
#if SCREEN_METHOD == 0
    scor = actual_screen5<dual_t>(state, scn_fac);
#elif SCREEN_METHOD == 1
    scor = chugunov2007<dual_t>(state, scn_fac);
#elif SCREEN_METHOD == 2
    scor = chugunov2009<dual_t>(state, scn_fac);
#elif SCREEN_METHOD == 3
    // null screening
    scor = 1.0_rt;
#elif SCREEN_METHOD == 4
    scor = chabrier1998<dual_t>(state, scn_fac);
#endif
    return scor;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
autodiff::dual actual_screen(const plasma_state_t& state,
                             const scrn::screen_factors_t& scn_fac)
{
    autodiff::dual scor;
    if (jacobian == 1) {
        scor = actual_screen<autodiff::dual>(state, scn_fac);
    } else {
        scor = actual_screen<amrex::Real>(state, scn_fac);
    }
    return scor;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_screen(const plasma_state_t& state,
                   const scrn::screen_factors_t& scn_fac,
                   Real& scor, Real& scordt)
{
    autodiff::dual scor_dual = actual_screen(state, scn_fac);
    scor = autodiff::val(scor_dual);
    scordt = autodiff::derivative<1>(scor_dual);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void screen(const plasma_state_t& state,
            const int jscreen,
            MICROPHYSICS_UNUSED const Real z1_screen, MICROPHYSICS_UNUSED const Real a1_screen,
            MICROPHYSICS_UNUSED const Real z2_screen, MICROPHYSICS_UNUSED const Real a2_screen,
            Real& scor, Real& scordt, Real& /*scordd*/)
{
    using namespace scrn;

#ifndef AMREX_USE_GPU
    AMREX_ASSERT(scn_facs[jscreen].validate_nuclei(z1_screen, a1_screen, z2_screen, a2_screen));
#endif

    actual_screen(state, scn_facs[jscreen], scor, scordt);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void screen(const plasma_state_t& state,
            const scrn::screen_factors_t& scn_fac,
            Real& scor, Real& scordt, Real& /*scordd*/)
{
    using namespace scrn;

    actual_screen(state, scn_fac, scor, scordt);
}

#ifdef AMREX_DEBUG
namespace {
// good for catching bugs with clang-tidy and ccls
[[maybe_unused]] static void __lint_helper(const plasma_state_t& state, const scrn::screen_factors_t& scn_fac)
{
  actual_screen5<autodiff::dual>(state, scn_fac);
  chugunov2007<autodiff::dual>(state, scn_fac);
  chugunov2009<autodiff::dual>(state, scn_fac);
  chabrier1998<autodiff::dual>(state, scn_fac);
  actual_screen5<amrex::Real>(state, scn_fac);
  chugunov2007<amrex::Real>(state, scn_fac);
  chugunov2009<amrex::Real>(state, scn_fac);
  chabrier1998<amrex::Real>(state, scn_fac);
}
}
#endif

#endif
