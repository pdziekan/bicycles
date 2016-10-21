#pragma once

//TODO: make these functions return arrays

// Grabowski & Smolarkiewicz 1996 eq. (2)
template <class ct_params_t>
void slvr_lgrngn<ct_params_t>::buoyancy(typename parent_t::arr_t &th, typename parent_t::arr_t &rv)
{
  const auto &ijk = this->ijk;

  namespace moist_air = libcloudphxx::common::moist_air;
  const real_t eps = moist_air::R_v<real_t>() / moist_air::R_d<real_t>() - 1.;
  tmp1(ijk).reindex(this->zero) = (libcloudphxx::common::earth::g<setup::real_t>() / si::metres_per_second_squared) * ((th(ijk).reindex(this->zero) - (*params.th_e)(this->vert_idx)) / (*params.th_ref)(this->vert_idx));

  this->smooth(tmp1, F);
}

