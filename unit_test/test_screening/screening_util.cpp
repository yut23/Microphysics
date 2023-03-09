#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>

#include <variables.H>
#include <network.H>
#include <eos.H>

#include <screen.H>

#include <cmath>

using namespace amrex;

constexpr scrn::screen_factors_t make_screening_combination(const Species::NetworkSpecies nuc_1, const Species::NetworkSpecies nuc_2)
{
  const Real Z1 = NetworkProperties::zion(nuc_1);
  const Real A1 = NetworkProperties::aion(nuc_1);

  const Real Z2 = NetworkProperties::zion(nuc_2);
  const Real A2 = NetworkProperties::aion(nuc_2);

  return scrn::calculate_screen_factor(Z1, A1, Z2, A2);
}

constexpr scrn::screen_factors_t make_screening_combination(const Real Z1, const Real A1, const Species::NetworkSpecies nuc_2)
{
  const Real Z2 = NetworkProperties::zion(nuc_2);
  const Real A2 = NetworkProperties::aion(nuc_2);

  return scrn::calculate_screen_factor(Z1, A1, Z2, A2);
}

constexpr scrn::screen_factors_t make_screening_combination(const Species::NetworkSpecies nuc_1, const Real Z2, const Real A2)
{
  const Real Z1 = NetworkProperties::zion(nuc_1);
  const Real A1 = NetworkProperties::aion(nuc_1);

  return scrn::calculate_screen_factor(Z1, A1, Z2, A2);
}

void screen_test_C(const Box& bx,
                   const Real dlogrho, const Real dlogT, const Real dmetal,
                   const plot_t vars,
                   Array4<Real> const sp) {
  using namespace Species;

  const int ih1 = network_spec_index("hydrogen-1");
  if (ih1 < 0) amrex::Error("Error: ih1 not found");

  const int ihe3 = network_spec_index("helium-3");
  if (ihe3 < 0) amrex::Error("Error: ihe3 not found");

  const int ihe4 = network_spec_index("helium-4");
  if (ihe4 < 0) amrex::Error("Error: ihe4 not found");

  const int ic12 = network_spec_index("carbon-12");
  if (ic12 < 0) amrex::Error("Error: ic12 not found");

  const int in14 = network_spec_index("nitrogen-14");
  if (in14 < 0) amrex::Error("Error: in14 not found");

  const int io16 = network_spec_index("oxygen-16");
  if (io16 < 0) amrex::Error("Error: io16 not found");

  const int ine20 = network_spec_index("neon-20");
  if (ine20 < 0) amrex::Error("Error: ine20 not found");

  const int img24 = network_spec_index("magnesium-24");
  if (img24 < 0) amrex::Error("Error: img24 not found");

  const int isi28 = network_spec_index("silicon-28");
  if (isi28 < 0) amrex::Error("Error: isi28 not found");

  const int is32 = network_spec_index("sulfur-32");
  if (is32 < 0) amrex::Error("Error: is32 not found");

  const int iar36 = network_spec_index("argon-36");    
  if (iar36 < 0) amrex::Error("Error: iar36 not found");

  const int ica40 = network_spec_index("calcium-40");
  if (ica40 < 0) amrex::Error("Error: ica40 not found");

  const int iti44 = network_spec_index("titanium-44");
  if (iti44 < 0) amrex::Error("Error: iti44 not found");

  const int icr48 = network_spec_index("chromium-48");
  if (icr48 < 0) amrex::Error("Error: icr48 not found");

  const int ife52 = network_spec_index("iron-52");
  if (ife52 < 0) amrex::Error("Error: ife52 not found");

  const int ife54 = network_spec_index("iron-54");
  if (ife54 < 0) amrex::Error("Error: ife54 not found");

  const int ife56 = network_spec_index("iron-56");
  if (ife56 < 0) amrex::Error("Error: ife56 not found");

  std::map<plot_t::index_t, scrn::screen_factors_t> screen_factors{};

  // 3-alpha
  {
    constexpr auto scn_fac = make_screening_combination(He4, He4);
    screen_factors[vars.iscn_he4_he4] = scn_fac;
  }

  {
    constexpr auto scn_fac = make_screening_combination(He4, 4.0_rt, 8.0_rt);
    screen_factors[vars.iscn_he4_be8] = scn_fac;
  }

  // c12(a,g)o16
  {
    constexpr auto scn_fac = make_screening_combination(C12, He4);
    screen_factors[vars.iscn_c12_he4] = scn_fac;
  }

  // c12 + c12
  {
    constexpr auto scn_fac = make_screening_combination(C12, C12);
    screen_factors[vars.iscn_c12_c12] = scn_fac;
  }

  // c12 + o16
  {
    constexpr auto scn_fac = make_screening_combination(C12, O16);
    screen_factors[vars.iscn_c12_o16] = scn_fac;
  }

  // o16 + o16
  {
    constexpr auto scn_fac = make_screening_combination(O16, O16);
    screen_factors[vars.iscn_o16_o16] = scn_fac;
  }

  // o16 + he4
  {
    constexpr auto scn_fac = make_screening_combination(O16, He4);
    screen_factors[vars.iscn_o16_he4] = scn_fac;
  }

  // ne20(a,g)mg24
  {
    constexpr auto scn_fac = make_screening_combination(Ne20, He4);
    screen_factors[vars.iscn_ne20_he4] = scn_fac;
  }

  // mg24(a,g)si28
  {
    constexpr auto scn_fac = make_screening_combination(Mg24, He4);
    screen_factors[vars.iscn_mg24_he4] = scn_fac;
  }

  // al27(p,g)si28
  {
    constexpr auto scn_fac = make_screening_combination(13.0_rt, 27.0_rt, H1);
    screen_factors[vars.iscn_al27_p] = scn_fac;
  }

  // si28 + he4
  {
    constexpr auto scn_fac = make_screening_combination(Si28, He4);
    screen_factors[vars.iscn_si28_he4] = scn_fac;
  }

  // p31(p,g)s32
  {
    constexpr auto scn_fac = make_screening_combination(15.0_rt, 31.0_rt, H1);
    screen_factors[vars.iscn_p31_p] = scn_fac;
  }

  // s32 to ar36
  {
    constexpr auto scn_fac = make_screening_combination(S32, He4);
    screen_factors[vars.iscn_s32_he4] = scn_fac;
  }

  // cl35(p,g)ar36
  {
    constexpr auto scn_fac = make_screening_combination(17.0_rt, 35.0_rt, H1);
    screen_factors[vars.iscn_cl35_p] = scn_fac;
  }

  // ar36 to ca40
  {
    constexpr auto scn_fac = make_screening_combination(Ar36, He4);
    screen_factors[vars.iscn_ar36_he4] = scn_fac;
  }

  // k39(p,g)ca40
  {
    constexpr auto scn_fac = make_screening_combination(19.0_rt, 39.0_rt, H1);
    screen_factors[vars.iscn_k39_p] = scn_fac;
  }

  // ca40 to ti44
  {
    constexpr auto scn_fac = make_screening_combination(Ca40, He4);
    screen_factors[vars.iscn_ca40_he4] = scn_fac;
  }

  // sc43(p,g)ti44
  {
    constexpr auto scn_fac = make_screening_combination(21.0_rt, 43.0_rt, H1);
    screen_factors[vars.iscn_sc43_p] = scn_fac;
  }

  // ti44 to cr48
  {
    constexpr auto scn_fac = make_screening_combination(Ti44, He4);
    screen_factors[vars.iscn_ti44_he4] = scn_fac;
  }

  // v47(p,g)cr48
  {
    constexpr auto scn_fac = make_screening_combination(23.0_rt, 47.0_rt, H1);
    screen_factors[vars.iscn_v47_p] = scn_fac;
  }

  // cr48 to fe52
  {
    constexpr auto scn_fac = make_screening_combination(Cr48, He4);
    screen_factors[vars.iscn_cr48_he4] = scn_fac;
  }

  // mn51(p,g)fe52
  {
    constexpr auto scn_fac = make_screening_combination(25.0_rt, 51.0_rt, H1);
    screen_factors[vars.iscn_mn51_p] = scn_fac;
  }

  // fe to ni
  {
    constexpr auto scn_fac = make_screening_combination(Fe52, He4);
    screen_factors[vars.iscn_fe52_he4] = scn_fac;
  }

  // co55(p,g)ni56
  {
    constexpr auto scn_fac = make_screening_combination(27.0_rt, 55.0_rt, H1);
    screen_factors[vars.iscn_co55_p] = scn_fac;
  }

  // fe54(p,g)co55
  {
    constexpr auto scn_fac = make_screening_combination(Fe54, H1);
    screen_factors[vars.iscn_fe54_p] = scn_fac;
  }

  // fe54(a,p)co57
  {
    constexpr auto scn_fac = make_screening_combination(Fe54, He4);
    screen_factors[vars.iscn_fe54_he4] = scn_fac;
  }

  // fe56(p,g)co57
  {
    constexpr auto scn_fac = make_screening_combination(Fe56, H1);
    screen_factors[vars.iscn_fe56_p] = scn_fac;
  }

  // d + p
  {
    constexpr auto scn_fac = make_screening_combination(1.0_rt, 2.0_rt, H1);
    screen_factors[vars.iscn_d_p] = scn_fac;
  }

  // pp
  {
    constexpr auto scn_fac = make_screening_combination(H1, H1);
    screen_factors[vars.iscn_p_p] = scn_fac;
  }

  // he3 + he3
  {
    constexpr auto scn_fac = make_screening_combination(He3, He3);
    screen_factors[vars.iscn_he3_he3] = scn_fac;
  }

  // he3 + he4
  {
    constexpr auto scn_fac = make_screening_combination(He3, He4);
    screen_factors[vars.iscn_he3_he4] = scn_fac;
  }

  // c12(p,g)n13
  {
    constexpr auto scn_fac = make_screening_combination(C12, H1);
    screen_factors[vars.iscn_c12_p] = scn_fac;
  }

  // n14(p,g)o15
  {
    constexpr auto scn_fac = make_screening_combination(N14, H1);
    screen_factors[vars.iscn_n14_p] = scn_fac;
  }

  // o16(p,g)f17
  {
    constexpr auto scn_fac = make_screening_combination(O16, H1);
    screen_factors[vars.iscn_o16_p] = scn_fac;
  }

  // n14(a,g)f18
  {
    constexpr auto scn_fac = make_screening_combination(N14, He4);
    screen_factors[vars.iscn_n14_he4] = scn_fac;
  }

  using dual_t = autodiff::dual;
  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k)
  {

    // set the composition -- approximately solar
    Real metalicity = 0.0 + static_cast<Real> (k) * dmetal;

    Real xn[NumSpec];

    // for now... the screening using 1-based indexing
    Array1D<dual_t, 1, NumSpec> ymass;

    for (int n = 0; n < NumSpec; n++) {
      xn[n] = metalicity / static_cast<Real>(NumSpec - 2);
    }
    xn[ih1] = 0.75_rt - 0.5_rt * metalicity;
    xn[ihe4] = 0.25_rt - 0.5_rt * metalicity;

    for (int n = 0; n < NumSpec; n++) {
      ymass(n+1) = xn[n] / aion[n];
    }

    dual_t temp_zone = std::pow(10.0, std::log10(temp_min) + static_cast<Real>(j)*dlogT);
    // seed the dual number for temperature before calculating anything with it
    autodiff::seed<1>(temp_zone, 1.0);

    Real dens_zone = std::pow(10.0, std::log10(dens_min) + static_cast<Real>(i)*dlogrho);

    // store default state
    sp(i, j, k, vars.irho) = dens_zone;
    sp(i, j, k, vars.itemp) = static_cast<Real>(temp_zone);
    for (int n = 0; n < NumSpec; n++) {
      sp(i, j, k, vars.ispec+n) = xn[n];
    }

    plasma_state_t<dual_t> pstate;
    fill_plasma_state(pstate, temp_zone, dens_zone, ymass);

    dual_t sc1a;

    for (auto &[var, scn_fac] : screen_factors)
    {
      actual_screen(pstate, scn_fac, sc1a);
      sp(i, j, k, var.value) = std::log(static_cast<Real>(sc1a));
      sp(i, j, k, var.dt) = autodiff::derivative<1>(sc1a) / static_cast<Real>(sc1a);
    }

  });

}
