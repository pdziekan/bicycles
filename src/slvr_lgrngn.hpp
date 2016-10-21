#pragma once
#include <chrono>
#include "slvr_dim.hpp"
#include "outmom.hpp"

#include <libcloudph++/lgrngn/factory.hpp>

#if defined(STD_FUTURE_WORKS)
#  include <future>
#endif

template <class ct_params_t>
class slvr_lgrngn : public slvr_dim<ct_params_t>
{
  using parent_t = slvr_dim<ct_params_t>; 
  using clock = std::chrono::high_resolution_clock; // TODO: option to disable timing, as it may affect performance a little?

  public:
  using ix = typename ct_params_t::ix;
  using real_t = typename ct_params_t::real_t;
  private:

  // member fields
  std::unique_ptr<libcloudphxx::lgrngn::particles_proto_t<real_t>> prtcls;

  // timing fields
  clock::time_point tbeg, tend, tbeg1, tend1, tbeg_loop;
  std::chrono::milliseconds tdiag, tupdate, tsync, tasync, tasync_wait, tloop, tvip_rhs; 

  // array with index of inversion
  blitz::Array<real_t, parent_t::n_dims-1> k_i;

  // global arrays, shared among threads, TODO: in fact no need to share them?
  typename parent_t::arr_t &tmp1,
                           &tmp2,
                           &r_l,
                           &F,       // forcings helper
                           &alpha,   // 'explicit' rhs part - does not depend on the value at n+1
                           &beta;    // 'implicit' rhs part - coefficient of the value at n+1
  // helper methods
  void diag()
  {
return;
  } 

  libcloudphxx::lgrngn::arrinfo_t<real_t> make_arrinfo(
    typename parent_t::arr_t arr
  ) {
    return libcloudphxx::lgrngn::arrinfo_t<real_t>(
      arr.dataZero(), 
      arr.stride().data()
    );
  }

  std::string aux_name(
    const std::string pfx, 
    const int rng,
    const int mom
  )
  { 
    std::ostringstream tmp;
    tmp << pfx << "_rng" << std::setw(3) << std::setfill('0') << rng << "_mom" << mom;
    return tmp.str();
  }

  protected:

  // deals with nitial supersaturation
  void hook_ante_loop(int nt)
  {
    parent_t::hook_ante_loop(nt); 

    // TODO: barrier?
    if (this->rank == 0) 
    {
      assert(params.backend != -1);
      assert(params.dt != 0); 
    }
    // TODO: barrier?
  }

  void vip_rhs_expl_calc()
  {
    parent_t::vip_rhs_expl_calc();
  }

  void buoyancy(typename parent_t::arr_t &th, typename parent_t::arr_t &rv);
  void radiation(typename parent_t::arr_t &rv);
  void rv_src();
  void th_src(typename parent_t::arr_t &rv);
  void w_src(typename parent_t::arr_t &th, typename parent_t::arr_t &rv);
  void surf_sens();
  void surf_latent();
  void subsidence(const int&);

  void update_rhs(
    arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at 
  )   
  {   
    parent_t::update_rhs(rhs, dt, at);
    this->mem->barrier();
    if(this->rank == 0)
      tbeg = clock::now();

    using ix = typename ct_params_t::ix;

    const auto &ijk = this->ijk;
    auto ix_w = this->vip_ixs[ct_params_t::n_dims - 1];


    // forcing
    switch (at) 
    {   
      // for eulerian integration or used to init trapezoidal integration
      case (0): 
      {   
        // ---- potential temp sources ----
        th_src(this->state(ix::th));
        rhs.at(ix::th)(ijk) += alpha(ijk) + beta(ijk) * this->state(ix::th)(ijk); 

        // vertical velocity sources
        if(params.w_src)
        {
          w_src(this->state(ix::th), this->state(ix::th));
          rhs.at(ix_w)(ijk) += alpha(ijk);
        }
        break;
      }   
      case (1): 
      // trapezoidal rhs^n+1
      {   
        // ---- potential temp sources ----
        tmp2(ijk) = this->state(ix::th)(ijk) + 0.5 * this->dt * rhs.at(ix::th)(ijk);
        // todo: once rv_src beta!=0 (e.g. nudging), rv^n+1 estimate should be implicit here
        th_src(tmp2);
        rhs.at(ix::th)(ijk) += (alpha(ijk) + beta(ijk) * this->state(ix::th)(ijk)) / (1. - 0.5 * this->dt * beta(ijk));
        // TODO: alpha should also take (possibly impolicit) estimate of th^n+1 too
        //       becomes important when nudging is introduced?

        // vertical velocity sources
        if(params.w_src)
        {
          // temporarily use beta to store the th^n+1 estimate
          beta(ijk) = this->state(ix::th)(ijk) + 0.5 * this->dt * rhs.at(ix::th)(ijk);
          // todo: oncethv_src beta!=0 (e.g. nudging), th^n+1 estimate should be implicit here

          // temporarily use F to store the rv^n+1 estimate
          F(ijk) = this->state(ix::th)(ijk) + 0.5 * this->dt * rhs.at(ix::th)(ijk);
          // todo: once rv_src beta!=0 (e.g. nudging), rv^n+1 estimate should be implicit here

          w_src(beta, F);
          rhs.at(ix_w)(ijk) += alpha(ijk);
        }
      }
    }  
    this->mem->barrier();
  }



  void hook_ante_step()
  {
//std::cout << this->timestep << std::endl;
    parent_t::hook_ante_step(); // includes output
  }
  // 
  void hook_post_step()
  {
    parent_t::hook_post_step(); // includes output
    this->mem->barrier();

  }

  public:

  struct rt_params_t : parent_t::rt_params_t 
  { 
    int backend = -1;
    bool async = true;
    libcloudphxx::lgrngn::opts_t<real_t> cloudph_opts;
    libcloudphxx::lgrngn::opts_init_t<real_t> cloudph_opts_init;
    outmom_t<real_t> out_dry, out_wet;
    bool flag_coal; // do we want coal after spinup
  };

  private:

  // per-thread copy of params
  rt_params_t params;

  public:

  // ctor
  slvr_lgrngn( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) : 
    parent_t(args, p),
    params(p),
    tmp1(args.mem->tmp[__FILE__][0][0]),
    tmp2(args.mem->tmp[__FILE__][0][5]),
    r_l(args.mem->tmp[__FILE__][0][2]),
    alpha(args.mem->tmp[__FILE__][0][3]),
    beta(args.mem->tmp[__FILE__][0][4]),
    F(args.mem->tmp[__FILE__][0][1])
  {
    k_i.resize(this->shape(this->hrzntl_domain));
    r_l = 0.;

    // TODO: equip rank() in libmpdata with an assert() checking if not in serial block
  }  

  static void alloc(typename parent_t::mem_t *mem, const int &n_iters)
  {
    parent_t::alloc(mem, n_iters);
    parent_t::alloc_tmp_sclr(mem, __FILE__, 6); // tmp1, tmp2, r_l, alpha, beta, F
  }
};
