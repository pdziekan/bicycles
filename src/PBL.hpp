#pragma once

// like pbl test from limpdataxx
// but with cdrag = 0 and no absorber

#include <iostream>

#include <libcloudph++/common/hydrostatic.hpp>
#include <libcloudph++/common/theta_std.hpp>
#include <libcloudph++/common/theta_dry.hpp>
#include <libcloudph++/common/lognormal.hpp>
#include <libcloudph++/common/unary_function.hpp>

#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>

namespace setup 
{
  const quantity<si::length, real_t> 
    Z    = 50 * 30 * si::metres, // DYCOMS: 1500
    X    = 65 * 50 * si::metres, // DYCOMS: 6400
    Y    = 65 * 50 * si::metres; // DYCOMS: 6400

  const real_t mixed_length = 500; // kg/kg

  const real_t z_abs = 100000; // no absorber

  template <class T, class U>
  void setopts_hlpr(T &params, const U &user_params)
  {
    params.outdir = user_params.outdir;
    params.outfreq = user_params.outfreq;
    params.w_src = user_params.w_src;
    params.uv_src = user_params.uv_src;
    params.th_src = user_params.th_src;
    params.rv_src = user_params.rv_src;
    params.prs_tol=1e-6;
    params.dt = user_params.dt;
    params.nt = user_params.nt;
  }

  // function expecting a libmpdata solver parameters struct as argument
  template <class T, class U>
  void setopts(T &params, int nx, int nz, const U &user_params)
  {
    setopts_hlpr(params, user_params);
    params.di = (X / si::metres) / (nx-1); 
    params.dj = (Z / si::metres) / (nz-1);
    params.dz = params.dj;
  }
  template <class T, class U>
  void setopts(T &params, int nx, int ny, int nz, const U &user_params)
  {
    setopts_hlpr(params, user_params);
    params.di = (X / si::metres) / (nx-1); 
    params.dj = (Y / si::metres) / (ny-1);
    params.dk = (Z / si::metres) / (nz-1);
    params.dz = params.dk;
  }


  template <class concurr_t, class index_t>
  void intcond_hlpr(concurr_t &solver, arr_1D_t &rhod, int rng_seed, index_t index)
  {
    using ix = typename concurr_t::solver_t::ix;
    int nz = solver.advectee().extent(ix::w);  // ix::w is the index of vertical domension both in 2D and 3D
    real_t dz = (Z / si::metres) / (nz-1); 

//    solver.advectee(ix::rv) = r_t()(index * dz); 
    solver.advectee(ix::u)= 0;// setup::u()(index * dz);
    solver.advectee(ix::w) = 0;  
   
    // absorbers
    solver.vab_coefficient() = where(index * dz >= z_abs,  1. / 100 * pow(sin(3.1419 / 2. * (index * dz - z_abs)/ (Z / si::metres - z_abs)), 2), 0);
    solver.vab_relaxed_state(0) = solver.advectee(ix::u);
    solver.vab_relaxed_state(ix::w) = 0; // vertical relaxed state

    // density profile
    solver.g_factor() = 1.;// rhod(index); // copy the 1D profile into 2D/3D array

    // initial potential temperature
//    solver.advectee(ix::th) = th_dry_fctr()(index * dz); 
    // randomly prtrb tht
    std::mt19937 gen(rng_seed);
    std::uniform_real_distribution<> dis(-0.5, 0.5);
    auto rand = std::bind(dis, gen);

    decltype(solver.advectee(ix::th)) prtrb(solver.advectee(ix::th).shape()); // array to store perturbation
    std::generate(prtrb.begin(), prtrb.end(), rand); // fill it, TODO: is it officialy stl compatible?
    prtrb *= where((1 - index * dz / mixed_length) > 0., (1 - index * dz / mixed_length), 0.);

    solver.advectee(ix::th) = 300 * where(index * dz <= mixed_length, 1, 1 + (index * dz - mixed_length) * st) + 1e-3 * prtrb;
    solver.advectee(ix::w) = 0.2 *  prtrb;
  }

  // function enforcing cyclic values in horizontal directions
  // 2D version
  template<int nd, class arr_t>
  void make_cyclic(arr_t arr,
    typename std::enable_if<nd == 2>::type* = 0)
  { arr(arr.extent(0) - 1, blitz::Range::all()) = arr(0, blitz::Range::all()); }

  // 3D version
  template<int nd, class arr_t>
  void make_cyclic(arr_t arr,
    typename std::enable_if<nd == 3>::type* = 0)
  { 
    arr(arr.extent(0) - 1, blitz::Range::all(), blitz::Range::all()) = 
      arr(0, blitz::Range::all(), blitz::Range::all()); 
    arr(blitz::Range::all(), arr.extent(1) - 1, blitz::Range::all()) = 
      arr(blitz::Range::all(), 0, blitz::Range::all());
  }

  // function expecting a libmpdata++ solver as argument
  // 2D version
  template <int nd, class concurr_t>
  void intcond(concurr_t &solver, arr_1D_t &rhod, int rng_seed,
    typename std::enable_if<nd == 2>::type* = 0
  )
  {
    blitz::secondIndex k;
    intcond_hlpr(solver, rhod, rng_seed, k);
    using ix = typename concurr_t::solver_t::ix;
    make_cyclic<2>(solver.advectee(ix::th));
    make_cyclic<2>(solver.advectee(ix::w));
  }

  // 3D version
  template <int nd, class concurr_t>
  void intcond(concurr_t &solver, arr_1D_t &rhod, int rng_seed,
    typename std::enable_if<nd == 3>::type* = 0
  )
  {
    blitz::thirdIndex k;
    intcond_hlpr(solver, rhod, rng_seed, k);
    using ix = typename concurr_t::solver_t::ix;
    make_cyclic<3>(solver.advectee(ix::th));
    make_cyclic<3>(solver.advectee(ix::w));

    solver.advectee(ix::v)= 0;// setup::v()(k * dz);
    solver.vab_relaxed_state(1) = solver.advectee(ix::v);
  }


  // calculate the initial environmental theta and rv profiles
  // alse set w_LS and hgt_fctrs
  // like in Wojtek's BabyEulag
  template<class user_params_t>
  void env_prof(arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &th_ref, arr_1D_t &rhod, arr_1D_t &w_LS, arr_1D_t &hgt_fctr_vctr, arr_1D_t &hgt_fctr_sclr, int nz, const user_params_t &user_params)
  {
    setup::real_t dz = (Z / si::metres) / (nz-1);
    blitz::thirdIndex k;
    real_t st = 1e-5;
    th_e = 300 * where(k * dz <= mixed_length, 1, 1 + (k * dz - mixed_length) * st);
    th_ref = 300;
  }
};
