#pragma once

#include <libmpdata++/solvers/mpdata_rhs_vip_prs.hpp>
#include <libmpdata++/output/hdf5_xdmf.hpp>

using namespace libmpdataxx; // TODO: get rid of it?

template <class ct_params_t>
class slvr_common : public 
  output::hdf5_xdmf<
    solvers::mpdata_rhs_vip_prs<ct_params_t>
  >
{
  using parent_t = output::hdf5_xdmf<solvers::mpdata_rhs_vip_prs<ct_params_t>>;

  public:
  using real_t = typename ct_params_t::real_t;

  protected:

  int spinup; // number of timesteps

  // surface precip stuff
  real_t prec_vol;
  std::ofstream f_prec;
  

  void hook_ante_loop(int nt) 
  {

    parent_t::hook_ante_loop(nt); 

    // open file for output of precitpitation volume
    f_prec.open(this->outdir+"/prec_vol.dat");
    prec_vol = 0.;
  }

  void hook_ante_step()
  {
    parent_t::hook_ante_step(); 
  }

  void hook_post_step()
  {
    parent_t::hook_post_step(); 

    // recording total precipitation volume through the lower boundary
    if(this->rank==0)
    {
      f_prec << this->timestep << " "  << prec_vol << "\n";
      prec_vol = 0.;
    }
  }

  // get shape from a rng_t or an idx_t
  inline int shape(const rng_t &rng) { return rng.length();}
  template<int n_dims>
  blitz::TinyVector<int, n_dims> shape(const idx_t<n_dims> &rng) { return rng.ubound() - rng.lbound() + 1;}

  public:

  struct rt_params_t : parent_t::rt_params_t 
  { 
    int spinup = 0, // number of timesteps during which autoconversion is to be turned off
        nt;         // total number of timesteps
    bool relax_th_rv, rv_src, th_src, uv_src, w_src;
    setup::arr_1D_t *th_e, *rv_e, *th_ref, *pre_ref, *rhod, *w_LS, *hgt_fctr_sclr, *hgt_fctr_vctr;
    typename ct_params_t::real_t dz; // vertical grid size
  };

  // ctor
  slvr_common( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) : 
    parent_t(args, p),
    spinup(p.spinup)
    {}  
};
