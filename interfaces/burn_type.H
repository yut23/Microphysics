#ifndef _burn_type_H_
#define _burn_type_H_

#include <iomanip>

#include <AMReX_REAL.H>

#include <network.H>
#include <eos_type.H>
#include <eos_composition.H>
#include <extern_parameters.H>
#include <microphysics_autodiff.H>

#include <ArrayUtilities.H>

using namespace amrex;

// A generic structure holding data necessary to do a nuclear burn.

// Set the number of independent variables -- this should be
// enuc + the number of species which participate
// in the evolution equations.

// For Strang evolution, this will be the number of equations we are
// evolving.  For simplified-SDC, we will need neqs when getting the
// reaction sources from the network's RHS.
const int neqs = 1 + NumSpec;

// Index of the energy variable in the work arrays.

const int net_ienuc = NumSpec + 1;

// this is the data type that is used to get the ydots from the actual
// RHS of the network, regardless of Strang or SDC
typedef amrex::Array1D<Real, 1, neqs> YdotNetArray1D;

#ifdef SIMPLIFIED_SDC

#if defined(SDC_EVOLVE_ENERGY)

// these indicies represent the order that the conserved state comes
// into the ODE integration from the hydro code.
//
// they also represent the order of the advective sources
//
// integrate rho*X, internal energy, total energy
// carry momentum as an unevolved variable

const int SFS   = 0;
const int SEINT = NumSpec;
// the following are not evolved
const int SFX   = SEINT + 1;
const int SRHO  = SFX + NumAux; // this is SFS + NumSpec if NumAux = 0;
const int SMX   = SRHO + 1;
const int SMY   = SRHO + 2;
const int SMZ   = SRHO + 3;
const int SEDEN = SMZ+1;

const int SVAR  = SEDEN+1;
const int SVAR_EVOLVE = SFX;

#elif defined(SDC_EVOLVE_ENTHALPY)

// integrate rho*X (species masses) and rho*h (enthalpy)
// carry pressure for EOS calls in RHS

const int SFS = 0;
const int SENTH = SFS + NumSpec;
const int SVAR  = SENTH+1;
const int SVAR_EVOLVE = SVAR;

#endif

#endif

// this is the data type of the dense Jacobian that the network wants.
// it is not the same size as the Jacobian that VODE cares about when
// we are doing simplified-SDC
typedef ArrayUtil::MathArray2D<1, neqs, 1, neqs> JacNetArray2D;

template <typename dual_t = amrex::Real>
struct burn_t
{

  // on exit, this will be set to the integration time we achieved

  Real time;

  // this first group are the quantities the network RHS uses

  Real rho;
  dual_t T;
  Real e;
  dual_t xn[NumSpec];
#if NAUX_NET > 0
  Real aux[NumAux];
#endif

  // derivatives needed for calling the EOS
  Real dedr;
  Real dedT;
#ifdef EXTRA_THERMO
  Real dedA;
  Real dedZ;
#endif

  // for diagnostics / error reporting

  int i;
  int j;
  int k;
#ifdef NONAKA_PLOT
  int level;
  Real reference_time;
#ifdef STRANG
  int strang_half;
#endif
#endif

  // now come the bits that we need for SDC or Strang evolution

#ifdef SIMPLIFIED_SDC

  // y is the input conserved state.  We will keep this state updated
  // in time as we integrate, such that upon output it will be the
  // final conserved state.
  Real y[SVAR];

  // we need to store a copy of the original state as well so we can
  // handle the non-evolved state evolution
  Real rho_orig;

#if defined(SDC_EVOLVE_ENERGY)
  Real rhoe_orig;
#endif

  // ydot_a are the advective terms that will modify the state y due
  // to hydrodynamics over the timestep.
  Real ydot_a[SVAR];

#if defined(SDC_EVOLVE_ENTHALPY)
  // Pressure in case we wish to use it for EOS calls
  Real p0;
#endif

  int sdc_iter;
  int num_sdc_iters;

#endif

  // for drive_initial_convection, we will fix T during the
  // integration to a passed in value.  We will interpret a positive
  // T_fixed as setting this feature.
  Real T_fixed{-1.0};

  // all coupling types need the specific heats to transform the
  // reaction Jacobian elements from T to e
  Real cv;

  // dx is useful for estimating timescales for equilibriation
  Real dx;

#ifdef TRUE_SDC
  Real E_var;
  Real mom[3];
  Real f_source[NumSpec+2];
#endif

  Real mu;
  Real mu_e;
  Real y_e;
  Real eta;
  dual_t abar;
  dual_t zbar;

#ifdef NSE_NET
  Real mu_p;
  Real mu_n;
#endif
  
  // diagnostics
  int n_rhs, n_jac, n_step;

  // Was the burn successful?
  bool success;
};


template <typename dual_t>
inline
std::ostream& operator<< (std::ostream& o, burn_t<dual_t> const& burn_state)
{
  o << std::setprecision(12);
  o << "rho = " << burn_state.rho << std::endl;
  o << "T =   " << burn_state.T << std::endl;
  o << "xn = ";
  for (int n = 0; n < NumSpec; ++n) {
    o << burn_state.xn[n] << " ";
  }
  o << std::endl;
#if NAUX_NET > 0
  o << "aux = ";
  for (int n = 0; n < NumAux; ++n) {
    o << burn_state.aux[n] << " ";
  }
  o << std::endl;
#endif

#ifdef SIMPLIFIED_SDC
#if defined(SDC_EVOLVE_ENERGY)
  o << "y[SRHO] = " << burn_state.y[SRHO] << std::endl;
#endif
  o << "y[SEINT] = " << burn_state.y[SEINT] << std::endl;

  o << "y[SFS:] = ";
  for (int n = 0; n < NumSpec; n++) {
      o << burn_state.y[SFS+n] << " ";
  }
  o << std::endl;

#if NAUX_NET > 0
  o << "y[SFX:] = ";
  for (int n = 0; n < NumAux; n++) {
      o << burn_state.y[SFX+n] << " ";
  }
  o << std::endl;
#endif

  o << "ydot_a[SRHO] = " << burn_state.ydot_a[SRHO] << std::endl;
  o << "ydot_a[SEINT] = " << burn_state.ydot_a[SEINT] << std::endl;

  o << "ydot_a[SFS:] = ";
  for (int n = 0; n < NumSpec; n++) {
      o << burn_state.ydot_a[SFS+n] << " ";
  }
  o << std::endl;

#if NAUX_NET > 0
  o << "ydot_a[SFX:] = ";
  for (int n = 0; n < NumAux; n++) {
      o << burn_state.ydot_a[SFX+n] << " ";
  }
  o << std::endl;
#endif
#endif
  return o;
}


// Copy data between different template specializations of burn_t
template <typename dual_t, typename dest_t>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void burn_to_burn(const burn_t<dual_t>& state, burn_t<dest_t>& new_state) {
    new_state.time = state.time;

    new_state.rho = state.rho;
    new_state.T = static_cast<dest_t>(state.T);
    new_state.e = state.e;
    for (int i = 0; i < NumSpec; ++i) {
        new_state.xn[i] = static_cast<dest_t>(state.xn[i]);
    }
#if NAUX_NET > 0
    std::copy(std::begin(state.aux), std::end(state.aux), std::begin(new_state.aux));
#endif

    new_state.dedr = state.dedr;
    new_state.dedT = state.dedT;
#ifdef EXTRA_THERMO
    new_state.dedA = state.dedA;
    new_state.dedZ = state.dedZ;
#endif

    new_state.i = state.i;
    new_state.j = state.j;
    new_state.k = state.k;
#ifdef NONAKA_PLOT
    new_state.level = state.level;
    new_state.reference_time = state.reference_time;
#ifdef STRANG
    new_state.strang_half = state.strang_half;
#endif
#endif

#ifdef SIMPLIFIED_SDC

    std::copy(std::begin(state.y), std::end(state.y), std::begin(new_state.y));
    new_state.rho_orig = state.rho_orig;

#if defined(SDC_EVOLVE_ENERGY)
    new_state.rhoe_orig = state.rhoe_orig;
#endif

    std::copy(std::begin(state.ydot_a), std::end(state.ydot_a), std::begin(new_state.ydot_a));

#if defined(SDC_EVOLVE_ENTHALPY)
    new_state.p0 = state.p0;
#endif

    new_state.sdc_iter = state.sdc_iter;
    new_state.num_sdc_iters = state.num_sdc_iters;

#endif

    new_state.T_fixed = state.T_fixed;
    new_state.cv = state.cv;
    new_state.dx = state.dx;

#ifdef TRUE_SDC
    new_state.E_var = state.E_var;
    std::copy(std::begin(state.mom), std::end(state.mom), std::begin(new_state.mom));
    std::copy(std::begin(state.f_source), std::end(state.f_source), std::begin(new_state.f_source));
#endif

    new_state.mu = state.mu;
    new_state.mu_e = state.mu_e;
    new_state.y_e = state.y_e;
    new_state.eta = state.eta;
    new_state.abar = static_cast<dest_t>(state.abar);
    new_state.zbar = static_cast<dest_t>(state.zbar);

#ifdef NSE_NET
    new_state.mu_p = state.mu_p;
    new_state.mu_n = state.mu_n;
#endif

    new_state.n_rhs = state.n_rhs;
    new_state.n_jac = state.n_jac;
    new_state.n_step = state.n_step;
    new_state.success = state.success;
}


// Given an eos type, copy the data relevant to the burn type.

template <typename T, typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eos_to_burn (const T& eos_state, burn_t<dual_t>& burn_state)
{
    burn_state.rho  = eos_state.rho;
    burn_state.T    = eos_state.T;
    burn_state.e    = eos_state.e;
    for (int n = 0; n < NumSpec; ++n) {
        burn_state.xn[n] = eos_state.xn[n];
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; ++n) {
        burn_state.aux[n] = eos_state.aux[n];
    }
#endif

    burn_state.cv   = eos_state.cv;

    burn_state.y_e  = eos_state.y_e;
    burn_state.eta  = eos_state.eta;
    burn_state.abar = eos_state.abar;
    burn_state.zbar = eos_state.zbar;

#ifdef NSE_NET
    burn_state.mu_p = eos_state.mu_p;
    burn_state.mu_n = eos_state.mu_n;
#endif
}



// Given a burn type, copy the data relevant to the eos type. 
// Note that when doing simplified SDC integration, we should
// avoid using this interface because the energy includes a
// contribution from the advection term. However this is useful
// for instantaneous RHS evaluations.

template <typename T, typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void burn_to_eos (const burn_t<dual_t>& burn_state, T& eos_state)
{
    eos_state.rho  = burn_state.rho;
    eos_state.T    = burn_state.T;
    eos_state.e    = burn_state.e;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = burn_state.xn[n];
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; ++n) {
        eos_state.aux[n] = burn_state.aux[n];
    }
#endif

    eos_state.cv   = burn_state.cv;

    eos_state.y_e  = burn_state.y_e;
    eos_state.eta  = burn_state.eta;
    eos_state.abar = burn_state.abar;
    eos_state.zbar = burn_state.zbar;

#ifdef NSE_NET
    eos_state.mu_p = burn_state.mu_p;
    eos_state.mu_n = burn_state.mu_n;
#endif
}


#ifdef SIMPLIFIED_SDC
template <typename T, typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void sdc_burn_to_eos (const burn_t<dual_t>& state, T& eos_state)
{
#if defined(SDC_EVOLVE_ENERGY)
    eos_state.rho = state.y[SRHO];

#elif defined(SDC_EVOLVE_ENTHALPY)
    eos_state.rho = 0.0_rt;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.rho += state.y[SFS+n];
    }
#endif

    eos_state.T = state.T;

    // species

    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state.y[SFS+n] / eos_state.rho;
    }

#if NAUX_NET > 0
    // aux
    for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = state.y[SFX+n] / eos_state.rho;
    }
    //set_nse_aux_from_X(eos_state);
#endif

    // we don't bother filling the other fields, since the EOS call
    // will do that for us
}

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void normalize_abundances_sdc_burn(burn_t<dual_t>& state)
{

    // Constrain the partial densities in burn_t state to sum to the 
    // density.
    //
    // This is meant to be used upon exit, and we assume that
    // vode_to_burn was already called

#ifdef SDC_EVOLVE_ENERGY

    // we have state.y[SRHO] to constrain to

    Real rhoX_sum = 0.0_rt;
    for (int n = 0; n < NumSpec; n++) {
        state.y[SFS+n] = amrex::max(amrex::min(state.y[SFS+n], state.y[SRHO]), 0.0_rt);
        rhoX_sum += state.y[SFS+n];
    }
    rhoX_sum /= state.y[SRHO];

    for (int n = 0; n < NumSpec; n++) {
        state.y[SFS+n] /= rhoX_sum;
    }

#ifdef AUX_THERMO
    // make the aux data consistent with the vode_state X's

    eos_re_t eos_state;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state.y[SFS+n] / state.y[SRHO];
    }

    set_aux_comp_from_X(eos_state);

    // also store it in the burn_t state

    for (int n = 0; n < NumAux; n++) {
        state.y[SFX+n] = state.y[SRHO] * eos_state.aux[n];
    }
#endif

#endif

}
#endif

template <typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void normalize_abundances_burn (burn_t<dual_t>& state)
{

    Real sum = 0.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
        state.xn[n] = admath::max(small_x, admath::min(1.0_rt, state.xn[n]));
        sum += state.xn[n];
    }
    for (int n = 0; n < NumSpec; ++n) {
        state.xn[n] /= sum;
    }

}

#endif
