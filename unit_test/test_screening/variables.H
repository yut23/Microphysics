#include <vector>
#include <string>

#ifndef _variables_H_
#define _variables_H_

#include <AMReX_Vector.H>

class plot_t {

public:

  struct index_t {
    int value = -1;
    int dt = -1;

    bool operator==(const index_t& rhs) const { return value == rhs.value; }
    bool operator<(const index_t& rhs)  const { return value < rhs.value;  }
  };

  int irho = -1;
  int itemp = -1;
  int ispec = -1;

  index_t iscn_he4_he4 = {};
  index_t iscn_he4_be8 = {};
  index_t iscn_c12_he4 = {};
  index_t iscn_c12_c12 = {};
  index_t iscn_c12_o16 = {};
  index_t iscn_o16_o16 = {};
  index_t iscn_o16_he4 = {};
  index_t iscn_ne20_he4 = {};
  index_t iscn_mg24_he4 = {};
  index_t iscn_al27_p = {};
  index_t iscn_si28_he4 = {};
  index_t iscn_p31_p = {};
  index_t iscn_s32_he4 = {};
  index_t iscn_cl35_p = {};
  index_t iscn_ar36_he4 = {};
  index_t iscn_k39_p = {};
  index_t iscn_ca40_he4 = {};
  index_t iscn_sc43_p = {};
  index_t iscn_ti44_he4 = {};
  index_t iscn_v47_p = {};
  index_t iscn_cr48_he4 = {};
  index_t iscn_mn51_p = {};
  index_t iscn_fe52_he4 = {};
  index_t iscn_co55_p = {};
  index_t iscn_fe54_p = {};
  index_t iscn_fe54_he4 = {};
  index_t iscn_fe56_p = {};
  index_t iscn_d_p = {};
  index_t iscn_p_p = {};
  index_t iscn_he3_he3 = {};
  index_t iscn_he3_he4 = {};
  index_t iscn_c12_p = {};
  index_t iscn_n14_p = {};
  index_t iscn_o16_p = {};
  index_t iscn_n14_he4 = {};

  int n_plot_comps = 0;

  int next_index(const int num) {
    int next = n_plot_comps;
    n_plot_comps += num;
    return next;
  }

};

plot_t init_variables();

void get_varnames(const plot_t p, amrex::Vector<std::string>& names);


#endif
