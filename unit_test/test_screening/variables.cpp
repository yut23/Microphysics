#include <variables.H>
#include <network.H>

plot_t init_variables() {

  plot_t p;

  p.irho = p.next_index(1);
  p.itemp = p.next_index(1);
  p.ispec = p.next_index(NumSpec);

  p.iscn_he4_he4.value = p.next_index(1);
  p.iscn_he4_be8.value = p.next_index(1);
  p.iscn_c12_he4.value = p.next_index(1);
  p.iscn_c12_c12.value = p.next_index(1);
  p.iscn_c12_o16.value = p.next_index(1);
  p.iscn_o16_o16.value = p.next_index(1);
  p.iscn_o16_he4.value = p.next_index(1);
  p.iscn_ne20_he4.value = p.next_index(1);
  p.iscn_mg24_he4.value = p.next_index(1);
  p.iscn_al27_p.value = p.next_index(1);
  p.iscn_si28_he4.value = p.next_index(1);
  p.iscn_p31_p.value = p.next_index(1);
  p.iscn_s32_he4.value = p.next_index(1);
  p.iscn_cl35_p.value = p.next_index(1);
  p.iscn_ar36_he4.value = p.next_index(1);
  p.iscn_k39_p.value = p.next_index(1);
  p.iscn_ca40_he4.value = p.next_index(1);
  p.iscn_sc43_p.value = p.next_index(1);
  p.iscn_ti44_he4.value = p.next_index(1);
  p.iscn_v47_p.value = p.next_index(1);
  p.iscn_cr48_he4.value = p.next_index(1);
  p.iscn_mn51_p.value = p.next_index(1);
  p.iscn_fe52_he4.value = p.next_index(1);
  p.iscn_co55_p.value = p.next_index(1);
  p.iscn_fe54_p.value = p.next_index(1);
  p.iscn_fe54_he4.value = p.next_index(1);
  p.iscn_fe56_p.value = p.next_index(1);
  p.iscn_d_p.value = p.next_index(1);
  p.iscn_p_p.value = p.next_index(1);
  p.iscn_he3_he3.value = p.next_index(1);
  p.iscn_he3_he4.value = p.next_index(1);
  p.iscn_c12_p.value = p.next_index(1);
  p.iscn_n14_p.value = p.next_index(1);
  p.iscn_o16_p.value = p.next_index(1);
  p.iscn_n14_he4.value = p.next_index(1);

  p.iscn_he4_he4.dt = p.next_index(1);
  p.iscn_he4_be8.dt = p.next_index(1);
  p.iscn_c12_he4.dt = p.next_index(1);
  p.iscn_c12_c12.dt = p.next_index(1);
  p.iscn_c12_o16.dt = p.next_index(1);
  p.iscn_o16_o16.dt = p.next_index(1);
  p.iscn_o16_he4.dt = p.next_index(1);
  p.iscn_ne20_he4.dt = p.next_index(1);
  p.iscn_mg24_he4.dt = p.next_index(1);
  p.iscn_al27_p.dt = p.next_index(1);
  p.iscn_si28_he4.dt = p.next_index(1);
  p.iscn_p31_p.dt = p.next_index(1);
  p.iscn_s32_he4.dt = p.next_index(1);
  p.iscn_cl35_p.dt = p.next_index(1);
  p.iscn_ar36_he4.dt = p.next_index(1);
  p.iscn_k39_p.dt = p.next_index(1);
  p.iscn_ca40_he4.dt = p.next_index(1);
  p.iscn_sc43_p.dt = p.next_index(1);
  p.iscn_ti44_he4.dt = p.next_index(1);
  p.iscn_v47_p.dt = p.next_index(1);
  p.iscn_cr48_he4.dt = p.next_index(1);
  p.iscn_mn51_p.dt = p.next_index(1);
  p.iscn_fe52_he4.dt = p.next_index(1);
  p.iscn_co55_p.dt = p.next_index(1);
  p.iscn_fe54_p.dt = p.next_index(1);
  p.iscn_fe54_he4.dt = p.next_index(1);
  p.iscn_fe56_p.dt = p.next_index(1);
  p.iscn_d_p.dt = p.next_index(1);
  p.iscn_p_p.dt = p.next_index(1);
  p.iscn_he3_he3.dt = p.next_index(1);
  p.iscn_he3_he4.dt = p.next_index(1);
  p.iscn_c12_p.dt = p.next_index(1);
  p.iscn_n14_p.dt = p.next_index(1);
  p.iscn_o16_p.dt = p.next_index(1);
  p.iscn_n14_he4.dt = p.next_index(1);


  return p;
}


void get_varnames(const plot_t p, amrex::Vector<std::string>& names) {

  names.resize(p.n_plot_comps);

  names[p.irho] = "density";
  names[p.itemp] = "temperature";
  for (int n = 0; n < NumSpec; n++) {
    names[p.ispec + n] = "X_" + spec_names_cxx[n];
  }

  names[p.iscn_he4_he4.value] = "scn_he4_he4";
  names[p.iscn_he4_be8.value] = "scn_he4_be8";
  names[p.iscn_c12_he4.value] = "scn_c12_he4";
  names[p.iscn_c12_c12.value] = "scn_c12_c12";
  names[p.iscn_c12_o16.value] = "scn_c12_o16";
  names[p.iscn_o16_o16.value] = "scn_o16_o16";
  names[p.iscn_o16_he4.value] = "scn_o16_he4";
  names[p.iscn_ne20_he4.value] = "scn_ne20_he4";
  names[p.iscn_mg24_he4.value] = "scn_mg24_he4";
  names[p.iscn_al27_p.value] = "scn_al27_p";
  names[p.iscn_si28_he4.value] = "scn_si28_he4";
  names[p.iscn_p31_p.value] = "scn_p31_p";
  names[p.iscn_s32_he4.value] = "scn_s32_he4";
  names[p.iscn_cl35_p.value] = "scn_cl35_p";
  names[p.iscn_ar36_he4.value] = "scn_ar36_he4";
  names[p.iscn_k39_p.value] = "scn_k39_p";
  names[p.iscn_ca40_he4.value] = "scn_ca40_he4";
  names[p.iscn_sc43_p.value] = "scn_sc43_p";
  names[p.iscn_ti44_he4.value] = "scn_ti44_he4";
  names[p.iscn_v47_p.value] = "scn_v47_p";
  names[p.iscn_cr48_he4.value] = "scn_cr48_he4";
  names[p.iscn_mn51_p.value] = "scn_mn51_p";
  names[p.iscn_fe52_he4.value] = "scn_fe52_he4";
  names[p.iscn_co55_p.value] = "scn_co55_p";
  names[p.iscn_fe54_p.value] = "scn_fe54_p";
  names[p.iscn_fe54_he4.value] = "scn_fe54_he4";
  names[p.iscn_fe56_p.value] = "scn_fe56_p";
  names[p.iscn_d_p.value] = "scn_d_p";
  names[p.iscn_p_p.value] = "scn_p_p";
  names[p.iscn_he3_he3.value] = "scn_he3_he3";
  names[p.iscn_he3_he4.value] = "scn_he3_he4";
  names[p.iscn_c12_p.value] = "scn_c12_p";
  names[p.iscn_n14_p.value] = "scn_n14_p";
  names[p.iscn_o16_p.value] = "scn_o16_p";
  names[p.iscn_n14_he4.value] = "scn_n14_he4";

  names[p.iscn_he4_he4.dt] = "scn_he4_he4_dt";
  names[p.iscn_he4_be8.dt] = "scn_he4_be8_dt";
  names[p.iscn_c12_he4.dt] = "scn_c12_he4_dt";
  names[p.iscn_c12_c12.dt] = "scn_c12_c12_dt";
  names[p.iscn_c12_o16.dt] = "scn_c12_o16_dt";
  names[p.iscn_o16_o16.dt] = "scn_o16_o16_dt";
  names[p.iscn_o16_he4.dt] = "scn_o16_he4_dt";
  names[p.iscn_ne20_he4.dt] = "scn_ne20_he4_dt";
  names[p.iscn_mg24_he4.dt] = "scn_mg24_he4_dt";
  names[p.iscn_al27_p.dt] = "scn_al27_p_dt";
  names[p.iscn_si28_he4.dt] = "scn_si28_he4_dt";
  names[p.iscn_p31_p.dt] = "scn_p31_p_dt";
  names[p.iscn_s32_he4.dt] = "scn_s32_he4_dt";
  names[p.iscn_cl35_p.dt] = "scn_cl35_p_dt";
  names[p.iscn_ar36_he4.dt] = "scn_ar36_he4_dt";
  names[p.iscn_k39_p.dt] = "scn_k39_p_dt";
  names[p.iscn_ca40_he4.dt] = "scn_ca40_he4_dt";
  names[p.iscn_sc43_p.dt] = "scn_sc43_p_dt";
  names[p.iscn_ti44_he4.dt] = "scn_ti44_he4_dt";
  names[p.iscn_v47_p.dt] = "scn_v47_p_dt";
  names[p.iscn_cr48_he4.dt] = "scn_cr48_he4_dt";
  names[p.iscn_mn51_p.dt] = "scn_mn51_p_dt";
  names[p.iscn_fe52_he4.dt] = "scn_fe52_he4_dt";
  names[p.iscn_co55_p.dt] = "scn_co55_p_dt";
  names[p.iscn_fe54_p.dt] = "scn_fe54_p_dt";
  names[p.iscn_fe54_he4.dt] = "scn_fe54_he4_dt";
  names[p.iscn_fe56_p.dt] = "scn_fe56_p_dt";
  names[p.iscn_d_p.dt] = "scn_d_p_dt";
  names[p.iscn_p_p.dt] = "scn_p_p_dt";
  names[p.iscn_he3_he3.dt] = "scn_he3_he3_dt";
  names[p.iscn_he3_he4.dt] = "scn_he3_he4_dt";
  names[p.iscn_c12_p.dt] = "scn_c12_p_dt";
  names[p.iscn_n14_p.dt] = "scn_n14_p_dt";
  names[p.iscn_o16_p.dt] = "scn_o16_p_dt";
  names[p.iscn_n14_he4.dt] = "scn_n14_he4_dt";

}


