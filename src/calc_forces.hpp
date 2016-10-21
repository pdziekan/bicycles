#pragma once
#include "forcings.hpp"

// helper functors
struct calc_c_p
{
  setup::real_t operator()(setup::real_t rv) const
  {return libcloudphxx::common::moist_air::c_p<setup::real_t>(rv) * si::kilograms * si::kelvins / si::joules;}
  BZ_DECLARE_FUNCTOR(calc_c_p)
};

struct calc_T
{
  setup::real_t operator()(setup::real_t th, setup::real_t rhod) const
  {return libcloudphxx::common::theta_dry::T<setup::real_t>(th * si::kelvins, rhod * si::kilograms / si::metres  / si::metres / si::metres) / si::kelvins;}
  BZ_DECLARE_FUNCTOR2(calc_T)
};


template <class ct_params_t>
void slvr_lgrngn<ct_params_t>::th_src(typename parent_t::arr_t &rv)
{
  const auto &ijk = this->ijk;
  if(params.th_src)
  {
    real_t hscale = 25;
    alpha(ijk).reindex(this->zero) = 0.01 * 1. / hscale * exp(- blitz::tensor::k * params.dz / hscale);
  }
  else
    alpha(ijk) = 0.;

  beta(ijk) = 0.;
  // nudging, todo: use some other coeff than vab_coeff
  //alpha(ijk).reindex(this->zero) += (*this->mem->vab_coeff)(ijk).reindex(this->zero) * (*params.th_e)(this->vert_idx);
  //beta(ijk) = - (*this->mem->vab_coeff)(ijk);
}

template <class ct_params_t>
void slvr_lgrngn<ct_params_t>::w_src(typename parent_t::arr_t &th, typename parent_t::arr_t &rv)
{
  const auto &ijk = this->ijk;
  // buoyancy
  buoyancy(th, rv);
  alpha(ijk) = F(ijk);
  // large-scale vertical wind
  //subsidence(ix::w); // TODO: in case 1, w here should be in step n+1, calc it explicitly as w + 0.5 * dt * rhs(w); 
                     //       could also be calculated implicitly, but we would need implicit w^n+1 in other cells;
                     //       also include absorber in w^n+1 estimate...

  //alpha(ijk) += F(ijk);
}


