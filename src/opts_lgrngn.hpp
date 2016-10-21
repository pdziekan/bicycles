/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <boost/assign/ptr_map_inserter.hpp>  // for 'ptr_map_insert()'

#include "opts_common.hpp"
#include "slvr_lgrngn.hpp"
#include "calc_forces.hpp"

// string parsing
#include <boost/spirit/include/qi.hpp>    
#include <boost/fusion/adapted/std_pair.hpp> 
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

// simulation and output parameters for micro=lgrngn
template <class solver_t, class user_params_t>
void setopts_micro(
  typename solver_t::rt_params_t &rt_params, 
  const user_params_t &user_params,
  typename std::enable_if<std::is_same<
    decltype(solver_t::rt_params_t::cloudph_opts),
    libcloudphxx::lgrngn::opts_t<typename solver_t::real_t>
  >::value>::type* = 0
)
{

  // output variables
  rt_params.outvars = {
    // <TODO>: make it common among all three micro?
    {solver_t::ix::th, {"th", "[K]"}},
//    {solver_t::ix::rv, {"rv", "[kg kg-1]"}},
//    {solver_t::ix::u, {"u", "[m/s]"}},
//    {solver_t::ix::v, {"v", "[m/s]"}},
//    {solver_t::ix::w, {"w", "[m/s]"}}
    // </TODO>
  };
}
