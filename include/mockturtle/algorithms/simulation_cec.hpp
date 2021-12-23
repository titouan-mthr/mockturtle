/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2021  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file simulation_cec.hpp
  \brief Simulation-based CEC

  EPFL CS-472 2021 Final Project Option 2
*/

/*
* ==============================================================================================================
* File           :  simulation_cec.hpp
* Description    :  Simulation-based CEC.
* Notes          :  For the project of Design Technologies for Integrated Systems.
*                :  Project option 2 : Combinational Equivalence Checking Using Circuit Simulation
* Author         :  Titouan MATHERET, Grenoble-INP & EPFL (MNIS), titouan.matheret@epfl.ch, Sciper: 339515
* ==============================================================================================================
*/

#pragma once

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>

#include "../utils/node_map.hpp"
#include "miter.hpp"
#include "simulation.hpp"

namespace mockturtle {

/* Statistics to be reported */
struct simulation_cec_stats {

  /*! \brief Split variable (simulation size). */
  uint32_t split_var{ 0 };

  /*! \brief Number of simulation round. */
  uint32_t rounds{ 0 };

};

namespace detail {

/* ---------------------------------------------------------------------------------------------------------
 * ------------------------------------------ CIRCUIT_SIMULATOR --------------------------------------------
 * --------------------------------------------------------------------------------------------------------- */

class circuit_simulator {

public:

  /* constructor of the class */
  circuit_simulator( unsigned int num_variables, unsigned int split_variables, unsigned int round )
      : num_variables( num_variables ), split_variables( split_variables ), round( round )
  {
  }

  /* returns the dynamic truth table or its complement depending on the value of the parameter "select"
   * select = 0 --> return the truth table
   * select = 1 --> return its complement */
  [[nodiscard]] kitty::dynamic_truth_table compute_constant( bool select ) const {
    kitty::dynamic_truth_table dyn_tt( split_variables );
    if ( select ) {
      return ~dyn_tt;
    }
    else {
      return dyn_tt;
    }
  }

  /* returns a small truth table depending on the index "ind" given as parameter */
  [[nodiscard]] kitty::dynamic_truth_table compute_pi( uint64_t ind ) const {
    kitty::dynamic_truth_table dyn_tt( split_variables );
    if ( ind < split_variables ) {
      kitty::create_nth_var( dyn_tt, ind );
    }
    else {
      if ( ( ( round >> ( ind - split_variables ) ) & 1 ) != 0 ) {
        dyn_tt = ~dyn_tt;
      }
    }
    return dyn_tt;
  }

  /* returns the complement of the truth table */
  [[nodiscard]] static kitty::dynamic_truth_table compute_not( kitty::dynamic_truth_table const& dtt ) {
    return ~dtt;
  }

private:

  /* elements of the class */
  unsigned int num_variables;
  unsigned int split_variables;
  unsigned int round;

};

/* ---------------------------------------------------------------------------------------------------------
 * ----------------------------------------- SIMULATION_CEC_IMPL -------------------------------------------
 * --------------------------------------------------------------------------------------------------------- */

template<class Ntk>

class simulation_cec_impl {

public:

  using pattern_t = unordered_node_map<kitty::dynamic_truth_table, Ntk>;
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

public:

  /* constructor of the class */
  explicit simulation_cec_impl( Ntk& ntk, simulation_cec_stats& st )
      : _ntk( ntk ), _st( st )
  {
  }

  /* compute the rounds required, generate the patterns and simulate truth tables
   * check the equivalence of all the small truth tables
   * return FALSE if circuits are NOT equivalent (a difference is spotted)
   * return TRUE if the circuits are equivalent (no difference spotted on any small truth tables) */
  bool run() {
    uint64_t cnt_round;
    compute_split_rounds();
    for ( cnt_round = 0 ; cnt_round < _st.rounds ; ++cnt_round ) {
      circuit_simulator circuit_simulator(_ntk.num_pis(), _st.split_var, cnt_round);
      const std::vector<kitty::dynamic_truth_table> dtt = simulate<kitty::dynamic_truth_table>(_ntk, circuit_simulator);
      for ( auto& po : dtt ) {
        if ( kitty::is_const0(po) == 0 ) {
          return FALSE;
        }
      }
    }
    /* no differences have been observed so the circuits are equivalent */
    return TRUE;
  }

private:

  /* define the variable limit for the truth table
   * in order to limit the memory usage to a maximum of 512 MB */
  void compute_split_rounds() {
    if ( _ntk.num_pis() <= 6 ) {
      _st.split_var = _ntk.num_pis();
    }
    else {
      uint64_t m;
      for ( m = 7 ; m < _ntk.num_pis() && (32 + pow(2, (m-3 + 1))) * _ntk.size() <= pow(2, 29) ; ++m ) {
      }
      _st.split_var = m;
    }
    _st.rounds = pow(2, _ntk.num_pis() - _st.split_var);
  }

private:

  /* elements of the class */
  Ntk& _ntk;
  simulation_cec_stats& _st;

};

} // namespace detail

/* Entry point for users to call */

/*! \brief Simulation-based CEC.
 *
 * This function implements a simulation-based combinational equivalence checker.
 * The implementation creates a miter network and run several rounds of simulation
 * to verify the functional equivalence. For memory and speed reasons this approach
 * is limited up to 40 input networks. It returns an optional which is `nullopt`,
 * if the network has more than 40 inputs.
 */
template<class Ntk>
std::optional<bool> simulation_cec( Ntk const& ntk1, Ntk const& ntk2, simulation_cec_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_num_pis_v<Ntk>, "Ntk does not implement the num_pis method" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_foreach_pi_v<Ntk>, "Ntk does not implement the foreach_pi method" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_is_complemented_v<Ntk>, "Ntk does not implement the is_complemented method" );

  simulation_cec_stats st;

  bool result = false;

  if ( ntk1.num_pis() > 40 )
    return std::nullopt;

  auto ntk_miter = miter<Ntk>( ntk1, ntk2 );

  if ( ntk_miter.has_value() )
  {
    detail::simulation_cec_impl p( *ntk_miter, st );
    result = p.run();
  }

  if ( pst )
    *pst = st;

  return result;
}

} // namespace mockturtle
