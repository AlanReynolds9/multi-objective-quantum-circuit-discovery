// MOQCD2 (Multi-objective quantum circuit discovery - algorithm 2.)
//
// Copyright (c) 2020  Alan Reynolds (alanreynolds9@gmail.com)
//
// This file is part of MOQCD2.
//
// MOQCD2 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//
// MOQCD2 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with MOQCD2.  If not, see
// <http://www.gnu.org/licenses/>.

//----------------------------------------------------------------------------------------------------------------------

#include <cassert>
#include "rng.h"
#include "constants.h"
#include "math.h"
#include "state.h"
#include "gate.h"
#include "circuit.h"
#include "replacements.h"
#include "newton.h"

using std::cerr;
using std::endl;
using std::vector;
using std::unique_ptr;
using std::make_shared;
using std::make_unique;
using constants::pi;
using utils::rand::randInt;
using utils::rand::randIntPair;
using utils::rand::randIntSet;
using utils::rand::randDouble;
using utils::rand::randNormal;  // Used only in mutate().
using circuit::CircuitContext;


int swapValues(int input, int i, int j)
{
  // Implements a function that is the identity on every value except i and j, which get switched.
  // (It might be better to take a reference to 'input' and change the value there. (I just tried using this function as
  // if that was the way it was written.)
  if (input == i)
  {
    return j;
  }
  if (input == j)
  {
    return i;
  }
  return input;
}


bool contains(const vector<int>& big, const vector<int>& small)
{
  // Returns true if the sorted vector of ints 'big' contains all the elements of the sorted vector 'small'.
  auto i = big.begin();
  for (auto bit : small)
  {
    for (; i != big.end() && *i < bit; ++i)
    {
    }
    if (i == big.end() || *i > bit)
    {
      return false;
    }
  }
  return true;
}


GateSequence makeGateSequence(unique_ptr<Gate>&& first)
{
  // Unfortunately, the various Gate::simplification() routines cannot initialized the returned GateSequence using an
  // initializer list such as {clone(), make_unique<PiByEight>(*context_, prev.target(), prev.controls())}, since one
  // cannot move from an initializer list and one cannot copy a unique_ptr. We therefore provide this overloaded
  // function to avoid creating the GateSequence manually in all the simplification routines, and to mimic the
  // initializer list approach. Note that we don't need this when the GateSequence has no gates.
  // (It is likely possible to use a variadic function or template to do this. However, since we only create short
  // GateSequences, we simply manually write each version.)
  GateSequence v(1);
  v[0] = std::move(first);
  return v;
}


GateSequence makeGateSequence(unique_ptr<Gate>&& first, unique_ptr<Gate>&& second)
{
  GateSequence v(2);
  v[0] = std::move(first);
  v[1] = std::move(second);
  return v;
}


GateSequence makeGateSequence(unique_ptr<Gate>&& first, unique_ptr<Gate>&& second, unique_ptr<Gate>&& third)
{
  GateSequence v(3);
  v[0] = std::move(first);
  v[1] = std::move(second);
  v[2] = std::move(third);
  return v;
}


// The following 6 templated functions are provided to avoid code duplication between the more concrete functions such
// as piByEightEquivalent()  and zGateEquivalentAvailable(). They should not be used, except via the more concrete
// versions. (DANGER!)
template <class G>
bool gateEquivalentAvailable(const CircuitContext& context, int target, const Controls& controls)
{
  // Assuming that a gate of type G can also be represented by an ArbitraryPhase of some particular angle, this function
  // returns whether either option is available.
  // ONLY CALL from functions like piByEightEquivalentAvailable(). DO NOT CALL for non phase type gates. (I'd like to
  // make this 'private'.)
  return G::gateAvailable(context, target, controls) || ArbitraryPhase::gateAvailable(context, target, controls);
}


template <class G>
long gateEquivalentCost(const CircuitContext& context, int target, const Controls& controls)
{
  // Assuming that a gate of type G can also be represented by an ArbitraryPhase of some particular angle, this function
  // returns the cost of the cheaper available option. We assume that gateEquivalentAvailable() has already been called
  // and returned 'true'.
  // ONLY CALL from functions like piByEightEquivalentCost(). DO NOT CALL for non phase type gates. (I'd like to make
  // this 'private'.)
  // (If we eliminate redundant gates, then G being available should mean that ArbitraryPhase is either unavailable or
  // more expensive. Hence we could eliminate the nested 'if'. However, this reliance on elimination of redundant gates
  // would make the code 'fragile', i.e. easily broken by future changes.
  if (G::gateAvailable(context, target, controls))
  {
    if (ArbitraryPhase::gateAvailable(context, target, controls) &&
        ArbitraryPhase::gateCost(context, controls.numControls()) < G::gateCost(context, controls.numControls()))
    {
      return ArbitraryPhase::gateCost(context, controls.numControls());
    }
    return G::gateCost(context, controls.numControls());
  }
  return ArbitraryPhase::gateCost(context, controls.numControls());
}


template <class G, int phaseAngle>
unique_ptr<Gate> gateEquivalent(const CircuitContext& context, int target, const Controls& controls)
{
  // Assuming that a gate of type G can also be represented by an ArbitraryPhase of angle 'phaseAngle', this function
  // returns the cheaper available option. No assumptions are made regarding whether a suitable gate is available. If
  // one is not, a nullptr is returned. The phaseAngle template parameter gives the angle for the ArbitraryPhase in
  // units of pi/4 radians.
  // ONLY CALL from functions like piByEightEquivalent(). DO NOT CALL for non phase type gates. MAKE SURE the phaseAngle
  // matches the gate type. (I'd like to make this 'private'.)
  // (If we eliminate redundant gates, then G being available should mean that ArbitraryPhase is either unavailable or
  // more expensive. Hence we could eliminate the nested 'if'. However, this reliance on eliminate of redundant gates
  // would make the code 'fragile', i.e. easily broken by future changes.

  if (G::gateAvailable(context, target, controls))
  {
    if (ArbitraryPhase::gateAvailable(context, target, controls) &&
        ArbitraryPhase::gateCost(context, controls.numControls()) < G::gateCost(context, controls.numControls()))
    {
      return make_unique<ArbitraryPhase>(context, target, controls, phaseAngle * pi / 4);
    }
    return make_unique<G>(context, target, controls);
  }

  if (ArbitraryPhase::gateAvailable(context, target, controls))
  {
    return make_unique<ArbitraryPhase>(context, target, controls, phaseAngle * pi / 4);
  }

  return {};
}


// STRUCTURE: Repeated code. Can we use variadic functions or something like that?
template <class G>
bool gateEquivalentAvailable(const CircuitContext& context, const Controls& controls)
{
  // Assuming that a gate of type G can also be represented by an ArbitraryPhase of some particular angle, this function
  // returns whether either option is available.
  // ONLY CALL from functions like piByEightEquivalentAvailable. DO NOT CALL for non phase type gates. (I'd like to make
  // this 'private'.)
  return G::gateAvailable(context, controls) || ArbitraryPhase::gateAvailable(context, controls);
}


template <class G>
long gateEquivalentCost(const CircuitContext& context, const Controls& controls)
{
  // Assuming that a gate of type G can also be represented by an ArbitraryPhase of some particular angle, this function
  // returns the cost of the cheaper available option. We assume that gateEquivalentAvailable() has already been called
  // and returned 'true'.
  // ONLY CALL from functions like piByEightEquivalentCost. DO NOT CALL for non phase type gates. (I'd like to make this
  // 'private'.)
  // (If we eliminate redundant gates, then G being available should mean that ArbitraryPhase is either unavailable or
  // more expensive. Hence we could eliminate the nested 'if'. However, this reliance on eliminate of redundant gates
  // would make the code 'fragile', i.e. easily broken by future changes.
  if (G::gateAvailable(context, controls))
  {
    if (ArbitraryPhase::gateAvailable(context, controls) &&
       ArbitraryPhase::gateCost(context, controls.numControls() - 1) < G::gateCost(context, controls.numControls() - 1))
    {
      return ArbitraryPhase::gateCost(context, controls.numControls() - 1);
    }
    return G::gateCost(context, controls.numControls() - 1);
  }
  return ArbitraryPhase::gateCost(context, controls.numControls() - 1);  // One of the gates in 'controls' becomes...
}                                                                        // ...the target


template <class G, int phaseAngle>
unique_ptr<Gate> gateEquivalent(const CircuitContext& context, const Controls& controls)
{
  // Assuming that a gate of type G can also be represented by an ArbitraryPhase of angle 'phaseAngle', this function
  // returns the cheaper available option. No assumptions are made regarding whether a suitable gate is available. If
  // one is not, a nullptr is returned. The phaseAngle template parameter gives the angle for the ArbitraryPhase in
  // units of pi/4 radians.
  // ONLY CALL from functions like piByEightEquivalent(). DO NOT CALL for non phase type gates. MAKE SURE the phaseAngle
  // matches the gate type. (I'd like to make this 'private'.)
  // (If we eliminate redundant gates, then G being available should mean that ArbitraryPhase is either unavailable or
  // more expensive. Hence we could eliminate the nested 'if'. However, this reliance on eliminate of redundant gates
  // would make the code 'fragile', i.e. easily broken by future changes.

  if (G::gateAvailable(context, controls))
  {
    if (ArbitraryPhase::gateAvailable(context, controls) &&
       ArbitraryPhase::gateCost(context, controls.numControls() - 1) < G::gateCost(context, controls.numControls() - 1))
    {
      return make_unique<ArbitraryPhase>(context, controls, phaseAngle * pi / 4);
    }
    return make_unique<G>(context, controls);
  }

  if (ArbitraryPhase::gateAvailable(context, controls))
  {
    return make_unique<ArbitraryPhase>(context, controls, phaseAngle * pi / 4);
  }

  return {};
}


bool piByEightEquivalentAvailable(const CircuitContext& context, int target, const Controls& controls)
{
  return gateEquivalentAvailable<PiByEight>(context, target, controls);
}


long piByEightEquivalentCost(const CircuitContext& context, int target, const Controls& controls)
{
  return gateEquivalentCost<PiByEight>(context, target, controls);
}


unique_ptr<Gate> piByEightEquivalent(const CircuitContext& context, int target, const Controls& controls)
{
  return gateEquivalent<PiByEight, 1>(context, target, controls);  // The 1 indicates that this is equivalent to one...
}                                                                  // ...PiByEight.


bool piByEightEquivalentAvailable(const CircuitContext& context, const Controls& controls)
{
  return gateEquivalentAvailable<PiByEight>(context, controls);
}


long piByEightEquivalentCost(const CircuitContext& context, const Controls& controls)
{
  return gateEquivalentCost<PiByEight>(context, controls);
}


unique_ptr<Gate> piByEightEquivalent(const CircuitContext& context, const Controls& controls)
{
  return gateEquivalent<PiByEight, 1>(context, controls);  // The 1 indicates that this is equivalent to one PiByEight.
}


bool piByEightInvEquivalentAvailable(const CircuitContext& context, int target, const Controls& controls)
{
  return gateEquivalentAvailable<PiByEightInv>(context, target, controls);
}


long piByEightInvEquivalentCost(const CircuitContext& context, int target, const Controls& controls)
{
  return gateEquivalentCost<PiByEightInv>(context, target, controls);
}


unique_ptr<Gate> piByEightInvEquivalent(const CircuitContext& context, int target, const Controls& controls)
{
  return gateEquivalent<PiByEightInv, -1>(context, target, controls);  // The -1 indicates that this is equivalent to...
}                                                                      // ...minus one PiByEights.


bool piByEightInvEquivalentAvailable(const CircuitContext& context, const Controls& controls)
{
  return gateEquivalentAvailable<PiByEightInv>(context, controls);
}


long piByEightInvEquivalentCost(const CircuitContext& context, const Controls& controls)
{
  return gateEquivalentCost<PiByEightInv>(context, controls);
}


unique_ptr<Gate> piByEightInvEquivalent(const CircuitContext& context, const Controls& controls)
{
  return gateEquivalent<PiByEightInv, -1>(context, controls);  // The -1 indicates that this is equivalent to minus...
}                                                              // ...one PiByEights.


bool phaseGateEquivalentAvailable(const CircuitContext& context, int target, const Controls& controls)
{
  return gateEquivalentAvailable<PhaseGate>(context, target, controls);
}


long phaseGateEquivalentCost(const CircuitContext& context, int target, const Controls& controls)
{
  return gateEquivalentCost<PhaseGate>(context, target, controls);
}


unique_ptr<Gate> phaseGateEquivalent(const CircuitContext& context, int target, const Controls& controls)
{
  return gateEquivalent<PhaseGate, 2>(context, target, controls);  // The 2 indicates that this is equivalent to two...
}                                                                  // ...PiByEights.


bool phaseGateEquivalentAvailable(const CircuitContext& context, const Controls& controls)
{
  return gateEquivalentAvailable<PhaseGate>(context, controls);
}


long phaseGateEquivalentCost(const CircuitContext& context, const Controls& controls)
{
  return gateEquivalentCost<PhaseGate>(context, controls);
}


unique_ptr<Gate> phaseGateEquivalent(const CircuitContext& context, const Controls& controls)
{
  return gateEquivalent<PhaseGate, 2>(context, controls);  // The 2 indicates that this is equivalent to two PiByEights.
}


bool phaseInvEquivalentAvailable(const CircuitContext& context, int target, const Controls& controls)
{
  return gateEquivalentAvailable<PhaseInv>(context, target, controls);
}


long phaseInvEquivalentCost(const CircuitContext& context, int target, const Controls& controls)
{
  return gateEquivalentCost<PhaseInv>(context, target, controls);
}


unique_ptr<Gate> phaseInvEquivalent(const CircuitContext& context, int target, const Controls& controls)
{
  return gateEquivalent<PhaseInv, -2>(context, target, controls);  // The -2 indicates that this is equivalent to...
}                                                                  // ...minus two PiByEights.


bool phaseInvEquivalentAvailable(const CircuitContext& context, const Controls& controls)
{
  return gateEquivalentAvailable<PhaseInv>(context, controls);
}


long phaseInvEquivalentCost(const CircuitContext& context, const Controls& controls)
{
  return gateEquivalentCost<PhaseInv>(context, controls);
}


unique_ptr<Gate> phaseInvEquivalent(const CircuitContext& context, const Controls& controls)
{
  return gateEquivalent<PhaseInv, -2>(context, controls);  // The -2 indicates that this is equivalent to minus two...
}                                                          // ...PiByEights.


bool zGateEquivalentAvailable(const CircuitContext& context, int target, const Controls& controls)
{
  return gateEquivalentAvailable<ZGate>(context, target, controls);
}


long zGateEquivalentCost(const CircuitContext& context, int target, const Controls& controls)
{
  return gateEquivalentCost<ZGate>(context, target, controls);
}


unique_ptr<Gate> zGateEquivalent(const CircuitContext& context, int target, const Controls& controls)
{
  return gateEquivalent<ZGate, 4>(context, target, controls);  // The 4 indicates that this is equivalent to four...
}                                                              // ...PiByEights.


bool zGateEquivalentAvailable(const CircuitContext& context, const Controls& controls)
{
  return gateEquivalentAvailable<ZGate>(context, controls);
}


long zGateEquivalentCost(const CircuitContext& context, const Controls& controls)
{
  return gateEquivalentCost<ZGate>(context, controls);
}


unique_ptr<Gate> zGateEquivalent(const CircuitContext& context, const Controls& controls)
{
  return gateEquivalent<ZGate, 4>(context, controls);  // The 2 indicates that this is equivalent to two PiByEights.
}


int totalGateCost(const vector<unique_ptr<Gate> >& gateSequence)
{
  int total{0};
  for (const auto& gate : gateSequence)
  {
    total += gate->cost();
  }

  return total;
}


Gate::Gate(const CircuitContext& context) :
context_(&context)
{
}


Gate::~Gate()
{
}


void Gate::createSimulators()
{
  // Default option is to assume that the gate needs no simulators - at least none attached to the gate.
}


void Gate::linkToMarkedState(std::shared_ptr<const int> markedState)
{
  // Most gates require no access to 'circuit data' such as the marked state.
}


void Gate::mutate()
{
  throw std::logic_error("Member function 'mutate()' called for gate with no parameters. This should not happen.");
}


int Gate::numParameters() const
{
  // The number of angle parameters associated with the gate. The default is zero. This will be overridden in
  // RotationGate and SU2Gate derived classes.
  return 0;
}


double Gate::parameter(int paramNum) const
{
  // Get the value of the parameter indicated. The default is to report an error - many gates are parameterless.
  // (Can we restructure so as to avoid 'don't call me' functions in the base class?)
  throw std::logic_error("Function Gate::parameter() cannot be called on a gate with no parameters.");
}


void Gate::setParameter(int paramNum, double value)
{
  // Set the value of the parameter indicated. The default is to report an error - many gates are parameterless.
  // STRUCTURE: Can we restructure so as to avoid 'don't call me' functions in the base class?
  throw std::logic_error("Function Gate::setParameter() cannot be called on a gate with no parameters.");
}


void Gate::randomizeParameters()
{
  // Randomize gate parameters. The default, of course, is to do nothing - many gates are parameterless.
}


long Gate::cost() const
{
  // Calls a function in the CircuitContext that, in turn, calls gateTypeId() (which is virtual) to get the type of
  // gate. Overridden in SingleTargetGate, which passes the correct number of controls.
  return context_->gateCost(*this, 0);
}


bool Gate::operator==(const Gate& rhs) const
{
  if (context_ != rhs.context_)
  {
    // We shouldn't be comparing gates with different 'contexts', i.e. different problems with differing numbers of
    // qbits, different costs, etc. If we do, we return 'false'.
    return false;
  }
  if (typeid(*this) == typeid(rhs))
  {
    return equals_(rhs);
  }
  return false;
}


bool Gate::operator!=(const Gate& rhs) const
{
  return !(*this == rhs);
}


bool Gate::sortBefore(const Gate& rhs) const
{
  if (context_ != rhs.context_)
  {
    // Shouldn't be comparing gates in different contexts, i.e. with different costs, num_qbits, etc.
    throw std::logic_error("Attempting to compare gates in different 'contexts'");
  }
  if (gateTypeId() == rhs.gateTypeId())  // Using gateTypeId() rather than typeid(*this) gives us greater...
  {                                      // ...control - we can put swaps at the end!
    return sort_before(rhs);
  }
  return gateTypeId() < rhs.gateTypeId();
}


void Gate::applyPartTo(State &state, int partIndex) const
{
  throw std::logic_error("Member function 'applyPartTo()' called for gate with no parameters. This should not happen.");
}


void Gate::optimizeAngle(const vector<vector<double> >& errorCoeffs)
{
  throw std::logic_error("Member function 'optimizeAngle()' called for a gate with no angles!");
}


bool Gate::mutatable() const
{
  // Only parameterized gates are mutatable. Default option, therefore, is to return 'false'.
  return false;
}


void Gate::output(std::ostream& out) const
{
  // Default output routine. (At present, only used for the Oracle derived class.
  out << name() << endl;
}


void Gate::extendedOutput(std::ostream &out) const
{
  // Default extended output routine.
  output(out);
}



bool Gate::canReduce() const
{
  // Since only gates with parameters can be 'reduced', i.e. simplified to a cheaper, non-parameterized gate or
  // eliminated, due to the special value of the angle parameter, the default option is to return 'false'.
  return false;
}


GateSequence Gate::reduction() const
{
  // This should never be called, as canReduce should only return 'true' if the gate is a RotationGate or an SU2Gate,
  // in which case the appropriate override runs instead.
  cerr << "Base class version of reduction() should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


bool Gate::canSimplySwap(const Gate& next) const
{
  // This function only runs if 'this' does not point to a SingleTargetGate. Hence 'this' points to an Oracle or a
  // SwapGate (or some potential future gate). Since 'next' could still be any type of gate, we return false - oracles,
  // in general, do not commute with other gates, while SwapGates do not commute in a simply manner - they swap the
  // qbits on the other gate. (The fact that Oracles commute with each other is irrelevant - consecutive oracles
  // cancel.) Whether we wish to detect 'simple' swaps involving a swap gate and some other gate that operates on
  // different qbits (or has both qbits as controls) might be worth thinking about. At present, I think treating swap
  // gates differently is useful - they can be shoved to the end of the circuit and then ignored during the
  // simplification of the rest of the circuit.) Other possible gate swaps are deferred to the overrides of this
  // function in SingleTargetGate and DiagonalGate.
  return false;
}


bool Gate::canSwapSwap(const Gate& next) const
{
  // Base class function for determining whether a Gate can be swapped with the next in the circuit, when one (or both)
  // gates is a SwapGate. Overrides exist for SingleTargetGate and SwapGate. Hence this version of the function only
  // runs if *this is an Oracle (or some future gate), which cannot be swapped with a SwapGate.
  return false;
}


unique_ptr<Gate> Gate::rightSwapMoverChange(const Gate &next) const
{
  // This should never be called, as canSwapSwap() should only return 'true' when *this a SingleTargetGate or a
  // SwapGate, for which we have overrides.
  cerr << "Base class version of rightSwapMoverChange(const Gate& next) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


unique_ptr<Gate> Gate::leftSwapMoverChange(const Gate &next) const
{
  // This should never be called, as canSwapSwap() should only return 'true' when *this a SingleTargetGate or a
  // SwapGate, for which we have overrides.
  cerr << "Base class version of leftSwapMoverChange(const Gate& next) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


bool Gate::canHSwap(const Gate& next) const
{
  // Returns true if one of the gates is a Hadamard, the other is a XGate or ZGate, they share target bits and the gates
  // can be swapped, with X<->Z. Other swaps of such gates, e.g. when the gates involve different sets of qbits, are
  // handled by canSimplySwap(). This function is overridden in the derived classes where these swaps are available.
  return false;
}


unique_ptr<Gate> Gate::rightHMoverChange(const Gate& next) const
{
  // This should never be called, as canHSwap() should only return 'true' when *this a HadamardGate, XGate or ZGate, for
  // which we have overrides.
  cerr << "Base class version of rightHMoverChange(const Gate& next) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


unique_ptr<Gate> Gate::leftHMoverChange(const Gate& next) const
{
  // This should never be called, as canHSwap() should only return 'true' when *this a HadamardGate, XGate or ZGate, for
  // which we have overrides.
  cerr << "Base class version of leftHMoverChange(const Gate& next) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


bool Gate::canRSwap(const Gate& next) const
{
  // Returns true if one of the gates is a Rotation (including ArbitraryPhase), the other is a XGate, YGate or ZGate,
  // they share target bits and the gates can be swapped with the angle being inverted. Other swaps of such gates, e.g.
  // when the gates involve different sets of qbits, are handled by canSimplySwap(). This function is overridden in the
  // derived classes where these swaps are available.
  return false;
}


unique_ptr<Gate> Gate::rightRMoverChange(const Gate& next) const
{
  // This should never be called, as canRSwap() should only return 'true' when *this a Rotation (including
  // ArbitraryPhase), XGate, YGate or ZGate, for which we have overrides.
  cerr << "Base class version of rightRMoverChange(const Gate& next) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


unique_ptr<Gate> Gate::leftRMoverChange(const Gate& next) const
{
  // This should never be called, as canRSwap() should only return 'true' when *this a Rotation (including
  // ArbitraryPhase), XGate, YGate or ZGate, for which we have overrides.
  cerr << "Base class version of leftRMoverChange(const Gate& next) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


bool Gate::canSimplySwap(const SingleTargetGate& prev) const
{
  // Again, since overrides for this function exist when *this is a SingleTargetGate, this leaves only the Oracle, which
  // doesn't commute with SingleTargetGates and the SwapGate. We deal with SwapGates separately, allowing us to shove
  // them out of the way.
  return false;
}


bool Gate::canSwapSwap(const SingleTargetGate& prev) const
{
  // Again, this version only runs if *this is an Oracle, which does not commute with anything.
  return false;
}


unique_ptr<Gate> Gate::rightSwapMoverChange(const SingleTargetGate &prev) const
{
  // This should never be called, as canSwapSwap() should only return 'true' when *this a SingleTargetGate or a
  // SwapGate, for which we have overrides.
  cerr << "Base class version of rightSwapMoverChange(const SingleTargetGate& next) should never run.
          " Aborting." << endl;
  exit(EXIT_FAILURE);
}


unique_ptr<Gate> Gate::leftSwapMoverChange(const SingleTargetGate &prev) const
{
  // This should never be called, as canSwapSwap() should only return 'true' when *this a SingleTargetGate or a
  // SwapGate, for which we have overrides.
  cerr << "Base class version of leftSwapMoverChange(const SingleTargetGate& next) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


bool Gate::canSimplySwap(const DiagonalGate& prev) const
{
  // Again, since overrides for this function exist when *this is a SingleTargetGate, this leaves only the SwapGate. We
  // deal with SwapGates separately, allowing us to shove them out of the way.
  return false;
}


long Gate::canSimplify(const Hadamard& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  // Any simplifications are handled by overrides in derived classes, If no override exists, no simplification takes
  // place and this default implementation returns 0, i.e. no improvement.
  return 0;
}


GateSequence Gate::simplification(const Hadamard& prev) const
{
  // This should never be called, as the fact that an override in one of the derived classes has not been triggered
  // indicates that no simplfication takes place.
  cerr << "Base class version of simplification(const Hadamard& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


bool Gate::canHSwap(const Hadamard& prev) const
{
  return false;
}


unique_ptr<Gate> Gate::rightHMoverChange(const Hadamard& prev) const
{
  // This should never be called, as canHSwap() should only return 'true' when *this a Hadamard, XGate, ZGate, XRotation
  // or ZRotation, for which we have overrides. The error message should catch any future mishap, when, for example, I
  // accidentally call this function directly, instead of as the second dispatch in double-dispatch.
  cerr << "Base class version of rightHMoverChange(const Hadamard& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


unique_ptr<Gate> Gate::leftHMoverChange(const Hadamard& prev) const
{
  // This should never be called, as canHSwap() should only return 'true' when *this a HadamardGate, XGate, ZGate,
  // XRotation or ZRotation, for which we have overrides. The error message should catch any future mishap, when, for
  // example, I accidentally call this function directly, instead of as the second dispatch in double-dispatch.
  cerr << "Base class version of leftHMoverChange(const Hadamard& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


long Gate::canSimplify(const PiByEight& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  // Any simplifications are handled by overrides in derived classes, If no override exists, no simplification takes
  // place and this default implementation returns 0, i.e. no improvement.
  return 0;
}


GateSequence Gate::simplification(const PiByEight& prev) const
{
  // This should never be called, as the fact that an override in one of the derived classes has not been triggered
  // indicates that no simplfication takes place.
  cerr << "Base class version of simplification(const PiByEight& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


long Gate::canSimplify(const PiByEightInv& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  // Any simplifications are handled by overrides in derived classes, If no override exists, no simplification takes
  // place and this default implementation returns 0, i.e. no improvement.
  return 0;
}


GateSequence Gate::simplification(const PiByEightInv& prev) const
{
  // This should never be called, as the fact that an override in one of the derived classes has not been triggered
  // indicates that no simplfication takes place.
  cerr << "Base class version of simplification(const PiByEightInv& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


long Gate::canSimplify(const PhaseGate& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  // Any simplifications are handled by overrides in derived classes, If no override exists, no simplification takes
  // place and this default implementation returns 0, i.e. no improvement.
  return 0;
}


GateSequence Gate::simplification(const PhaseGate& prev) const
{
  // This should never be called, as the fact that an override in one of the derived classes has not been triggered
  // indicates that no simplfication takes place.
  cerr << "Base class version of simplification(const PhaseGate& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


bool Gate::canRSwap(const PhaseGate& prev) const
{
  return false;
}


unique_ptr<Gate> Gate::rightRMoverChange(const PhaseGate& prev) const
{
  // This should never be called, as canRSwap() should only return 'true' when *this a Rotation (including
  // ArbitraryPhase), XGate, YGate, ZGate, PhaseGate or PhaseInv, for which we have overrides. The error message should
  // catch any future mishap, when, for example, I accidentally call this function directly, instead of as the second
  // dispatch in double-dispatch.
  cerr << "Base class version of rightRMoverChange(const PhaseGate& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


unique_ptr<Gate> Gate::leftRMoverChange(const PhaseGate& prev) const
{
  // This should never be called, as canRSwap() should only return 'true' when *this a Rotation (including
  // ArbitraryPhase), XGate, YGate, ZGate, PhaseGate or Phase Inv, for which we have overrides. The error message should
  // catch any future mishap, when, for example, I accidentally call this function directly, instead of as the second
  // dispatch in double-dispatch.
  cerr << "Base class version of leftRMoverChange(const PhaseGate& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


long Gate::canSimplify(const PhaseInv& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  // Any simplifications are handled by overrides in derived classes, If no override exists, no simplification takes
  // place and this default implementation returns 0, i.e. no improvement.
  return 0;
}


GateSequence Gate::simplification(const PhaseInv& prev) const
{
  // This should never be called, as the fact that an override in one of the derived classes has not been triggered
  // indicates that no simplfication takes place.
  cerr << "Base class version of simplification(const PhaseInv& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


bool Gate::canRSwap(const PhaseInv& prev) const
{
  return false;
}


unique_ptr<Gate> Gate::rightRMoverChange(const PhaseInv& prev) const
{
  // This should never be called, as canRSwap() should only return 'true' when *this a Rotation (including
  // ArbitraryPhase), XGate, YGate, ZGate, PhaseGate or PhaseInv, for which we have overrides. The error message should
  // catch any future mishap, when, for example, I accidentally call this function directly, instead of as the second
  // dispatch in double-dispatch.
  cerr << "Base class version of rightRMoverChange(const PhaseInv& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


unique_ptr<Gate> Gate::leftRMoverChange(const PhaseInv& prev) const
{
  // This should never be called, as canRSwap() should only return 'true' when *this a Rotation (including
  // ArbitraryPhase), XGate, YGate, ZGate, PhaseGate or Phase Inv, for which we have overrides. The error message should
  // catch any future mishap, when, for example, I accidentally call this function directly, instead of as the second
  // dispatch in double-dispatch.
  cerr << "Base class version of leftRMoverChange(const PhaseInv& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


long Gate::canSimplify(const XGate& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  // Any simplifications are handled by overrides in derived classes, If no override exists, no simplification takes
  // place and this default implementation returns 0, i.e. no improvement.
  return 0;
}


GateSequence Gate::simplification(const XGate& prev) const
{
  // This should never be called, as the fact that an override in one of the derived classes has not been triggered
  // indicates that no simplfication takes place.
  cerr << "Base class version of simplification(const XGate& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


bool Gate::canHSwap(const XGate& prev) const
{
  return false;
}


unique_ptr<Gate> Gate::rightHMoverChange(const XGate& prev) const
{
  // This should never be called, as canHSwap() should only return 'true' when *this a HadamardGate, XGate, ZGate,
  // XRotation or ZRotation, for which we have overrides. The error message should catch any future mishap, when, for
  // example, I accidentally call this function directly, instead of as the second dispatch in double-dispatch.
  cerr << "Base class version of rightHMoverChange(const XGate& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


unique_ptr<Gate> Gate::leftHMoverChange(const XGate& prev) const
{
  // This should never be called, as canHSwap() should only return 'true' when *this a HadamardGate, XGate, ZGate,
  // XRotation or ZRotation, for which we have overrides. The error message should catch any future mishap, when, for
  // example, I accidentally call this function directly, instead of as the second dispatch in double-dispatch.
  cerr << "Base class version of leftHMoverChange(const XGate& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


bool Gate::canRSwap(const XGate& prev) const
{
  return false;
}


unique_ptr<Gate> Gate::rightRMoverChange(const XGate& prev) const
{
  // This should never be called, as canRSwap() should only return 'true' when *this a Rotation (including
  // ArbitraryPhase), XGate, YGate, ZGate, PhaseGate or PhaseInv, for which we have overrides. The error message should
  // catch any future mishap, when, for example, I accidentally call this function directly, instead of as the second
  // dispatch in double-dispatch.
  cerr << "Base class version of rightRMoverChange(const XGate& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


unique_ptr<Gate> Gate::leftRMoverChange(const XGate& prev) const
{
  // This should never be called, as canRSwap() should only return 'true' when *this a Rotation (including
  // ArbitraryPhase), XGate, YGate, ZGate, PhaseGate or PhaseInv, for which we have overrides. The error message should
  // catch any future mishap, when, for example, I accidentally call this function directly, instead of as the second
  // dispatch in double-dispatch.
  cerr << "Base class version of leftRMoverChange(const XGate& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


long Gate::canSimplify(const YGate& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  // Any simplifications are handled by overrides in derived classes, If no override exists, no simplification takes
  // place and this default implementation returns 0, i.e. no improvement.
  return 0;
}


GateSequence Gate::simplification(const YGate& prev) const
{
  // This should never be called, as the fact that an override in one of the derived classes has not been triggered
  // indicates that no simplfication takes place.
  cerr << "Base class version of simplification(const YGate& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


bool Gate::canHSwap(const YGate& prev) const
{
  return false;
}


unique_ptr<Gate> Gate::rightHMoverChange(const YGate& prev) const
{
  // This should never be called, as canHSwap() should only return 'true' when *this is a HadamardGate and we have an
  // override for that case. The error message should catch any future mishap, when, for example, I accidentally call
  // this function directly, instead of as the second dispatch in double-dispatch.
  cerr << "Base class version of rightHMoverChange(const YGate& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


unique_ptr<Gate> Gate::leftHMoverChange(const YGate& prev) const
{
  // This should never be called, as canHSwap() should only return 'true' when *this is a HadamardGate and we have an
  // override for that case.  The error message should catch any future mishap, when, for example, I accidentally call
  // this function directly, instead of as the second dispatch in double-dispatch.
  cerr << "Base class version of leftHMoverChange(const YGate& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


bool Gate::canRSwap(const YGate& prev) const
{
  return false;
}


unique_ptr<Gate> Gate::rightRMoverChange(const YGate& prev) const
{
  // This should never be called, as canRSwap() should only return 'true' when *this a Rotation (including
  // ArbitraryPhase), XGate, YGate, ZGate, PhaseGate or PhaseInv, for which we have overrides. The error message should
  // catch any future mishap, when, for example, I accidentally call this function directly, instead of as the second
  // dispatch in double-dispatch.
  cerr << "Base class version of rightRMoverChange(const YGate& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


unique_ptr<Gate> Gate::leftRMoverChange(const YGate& prev) const
{
  // This should never be called, as canRSwap() should only return 'true' when *this a Rotation (including
  // ArbitraryPhase), XGate, YGate, ZGate, PhaseGate or PhaseInv, for which we have overrides. The error message should
  // catch any future mishap, when, for example, I accidentally call this function directly, instead of as the second
  // dispatch in double-dispatch.
  cerr << "Base class version of leftRMoverChange(const YGate& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


long Gate::canSimplify(const ZGate& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  // Any simplifications are handled by overrides in derived classes, If no override exists, no simplification takes
  // place and this default implementation returns 0, i.e. no improvement.
  return 0;
}


GateSequence Gate::simplification(const ZGate& prev) const
{
  // This should never be called, as the fact that an override in one of the derived classes has not been triggered
  // indicates that no simplfication takes place.
  cerr << "Base class version of simplification(const ZGate& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


bool Gate::canHSwap(const ZGate& prev) const
{
  return false;
}


unique_ptr<Gate> Gate::rightHMoverChange(const ZGate& prev) const
{
  // This should never be called, as canHSwap() should only return 'true' when *this a HadamardGate, XGate, ZGate,
  // XRotation or ZRotation, for which we have overrides. The error message should catch any future mishap, when, for
  // example, I accidentally call this function directly, instead of as the second dispatch in double-dispatch.
  cerr << "Base class version of rightHMoverChange(const ZGate& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


unique_ptr<Gate> Gate::leftHMoverChange(const ZGate& prev) const
{
  // This should never be called, as canHSwap() should only return 'true' when *this a HadamardGate, XGate, ZGate,
  // XRotation or ZRotation, for which we have overrides. The error message should catch any future mishap, when, for
  // example, I accidentally call this function directly, instead of as the second dispatch in double-dispatch.
  cerr << "Base class version of leftHMoverChange(const ZGate& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


bool Gate::canRSwap(const ZGate& prev) const
{
  return false;
}


unique_ptr<Gate> Gate::rightRMoverChange(const ZGate& prev) const
{
  // This should never be called, as canRSwap() should only return 'true' when *this a Rotation (including
  // ArbitraryPhase), XGate, YGate, ZGate, PhaseGate or PhaseInv, for which we have overrides. The error message should
  // catch any future mishap, when, for example, I accidentally call this function directly, instead of as the second
  // dispatch in double-dispatch.
  cerr << "Base class version of rightRMoverChange(const ZGate& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


unique_ptr<Gate> Gate::leftRMoverChange(const ZGate& prev) const
{
  // This should never be called, as canRSwap() should only return 'true' when *this a Rotation (including
  // ArbitraryPhase), XGate, YGate, ZGate, PhaseGate or PhaseInv, for which we have overrides. The error message should
  // catch any future mishap, when, for example, I accidentally call this function directly, instead of as the second
  // dispatch in double-dispatch.
  cerr << "Base class version of leftRMoverChange(const ZGate& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


long Gate::canSimplify(const XRotation& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  // Any simplifications are handled by overrides in derived classes, If no override exists, no simplification takes
  // place and this default implementation returns 0, i.e. no improvement.
  return 0;
}


GateSequence Gate::simplification(const XRotation& prev) const
{
  // This should never be called, as the fact that an override in one of the derived classes has not been triggered
  // indicates that no simplfication takes place.
  cerr << "Base class version of simplification(const XRotation& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


bool Gate::canHSwap(const XRotation& prev) const
{
  return false;
}


unique_ptr<Gate> Gate::rightHMoverChange(const XRotation& prev) const
{
  // This should never be called, as canHSwap() should only return 'true' when *this a HadamardGate, XGate, ZGate,
  // XRotation or ZRotation, for which we have overrides. The error message should catch any future mishap, when, for
  // example, I accidentally call this function directly, instead of as the second dispatch in double-dispatch.
  cerr << "Base class version of rightHMoverChange(const XRotation& next) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


unique_ptr<Gate> Gate::leftHMoverChange(const XRotation& prev) const
{
  // This should never be called, as canHSwap() should only return 'true' when *this a HadamardGate, XGate, ZGate,
  // XRotation or ZRotation, for which we have overrides. The error message should catch any future mishap, when, for
  // example, I accidentally call this function directly, instead of as the second dispatch in double-dispatch.
  cerr << "Base class version of leftHMoverChange(const XRotation& next) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


bool Gate::canRSwap(const XRotation& prev) const
{
  return false;
}


unique_ptr<Gate> Gate::rightRMoverChange(const XRotation& prev) const
{
  // This should never be called, as canRSwap() should only return 'true' when *this a Rotation (including
  // ArbitraryPhase), XGate, YGate, ZGate, PhaseGate or PhaseInv, for which we have overrides. The error message should
  // catch future mishaps, if, for example, I accidentally call this function directly, instead of as the second
  // dispatch in double-dispatch.
  cerr << "Base class version of rightRMoverChange(const XRotation& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


unique_ptr<Gate> Gate::leftRMoverChange(const XRotation& prev) const
{
  // This should never be called, as canRSwap() should only return 'true' when *this a Rotation (including
  // ArbitraryPhase), XGate, YGate, ZGate, PhaseGate or PhaseInv, for which we have overrides. The error message should
  // catch future mishaps, if, for example, I accidentally call this function directly, instead of as the second
  // dispatch in double-dispatch.
  cerr << "Base class version of leftRMoverChange(const XRotation& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


long Gate::canSimplify(const YRotation& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  // Any simplifications are handled by overrides in derived classes, If no override exists, no simplification takes
  // place and this default implementation returns 0, i.e. no improvement.
  return 0;
}


GateSequence Gate::simplification(const YRotation& prev) const
{
  // This should never be called, as the fact that an override in one of the derived classes has not been triggered
  // indicates that no simplfication takes place.
  cerr << "Base class version of simplification(const YRotation& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


bool Gate::canHSwap(const YRotation& prev) const
{
  return false;
}


unique_ptr<Gate> Gate::rightHMoverChange(const YRotation& prev) const
{
  // This should never be called, as canHSwap() should only return 'true' when *this is a HadamardGate and we have an
  // override for that case. The error message should catch future mishaps, if, for example, I accidentally call this
  // function directly, instead of as the second dispatch in double-dispatch.
  cerr << "Base class version of rightHMoverChange(const YRotation& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


unique_ptr<Gate> Gate::leftHMoverChange(const YRotation& prev) const
{
  // This should never be called, as canHSwap() should only return 'true' when *this is a HadamardGate and we have an
  // override for that case. The error message should catch future mishaps, if, for example, I accidentally call this
  // function directly, instead of as the second dispatch in double-dispatch.
  cerr << "Base class version of leftHMoverChange(const YRotation& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


bool Gate::canRSwap(const YRotation& prev) const
{
  return false;
}


unique_ptr<Gate> Gate::rightRMoverChange(const YRotation& prev) const
{
  // This should never be called, as canRSwap() should only return 'true' when *this a Rotation (including
  // ArbitraryPhase), XGate, YGate, ZGate, PhaseGate or PhaseInv, for which we have overrides. The error message should
  // catch future mishaps, if, for example, I accidentally call this function directly, instead of as the second
  // dispatch in double-dispatch.
  cerr << "Base class version of rightRMoverChange(const YRotation& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


unique_ptr<Gate> Gate::leftRMoverChange(const YRotation& prev) const
{
  // This should never be called, as canRSwap() should only return 'true' when *this a Rotation (including
  // ArbitraryPhase), XGate, YGate, ZGate, PhaseGate or PhaseInv, for which we have overrides. The error message should
  // catch future mishaps, if, for example, I accidentally call this function directly, instead of as the second
  // dispatch in double-dispatch.
  cerr << "Base class version of leftRMoverChange(const YRotation& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


long Gate::canSimplify(const ZRotation& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  // Any simplifications are handled by overrides in derived classes, If no override exists, no simplification takes
  // place and this default implementation returns 0, i.e. no improvement.
  return 0;
}


GateSequence Gate::simplification(const ZRotation& prev) const
{
  // This should never be called, as the fact that an override in one of the derived classes has not been triggered
  // indicates that no simplfication takes place.
  cerr << "Base class version of simplification(const XRotation& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


bool Gate::canHSwap(const ZRotation& prev) const
{
  return false;
}


unique_ptr<Gate> Gate::rightHMoverChange(const ZRotation& prev) const
{
  // This should never be called, as canHSwap() should only return 'true' when *this a HadamardGate, XGate, ZGate,
  // XRotation or ZRotation, for which we have overrides. The error message should catch future mishaps, if, for
  // example, I accidentally call this function directly, instead of as the second dispatch in double-dispatch.
  cerr << "Base class version of rightHMoverChange(const ZRotation& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


unique_ptr<Gate> Gate::leftHMoverChange(const ZRotation& prev) const
{
  // This should never be called, as canHSwap() should only return 'true' when *this a HadamardGate, XGate, ZGate,
  // XRotation or ZRotation, for which we have overrides. The error message should catch future mishaps, if, for
  // example, I accidentally call this function directly, instead of as the second dispatch in double-dispatch.
  cerr << "Base class version of leftHMoverChange(const ZRotation& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


bool Gate::canRSwap(const ZRotation& prev) const
{
  return false;
}


unique_ptr<Gate> Gate::rightRMoverChange(const ZRotation& prev) const
{
  // This should never be called, as canRSwap() should only return 'true' when *this a Rotation (including
  // ArbitraryPhase), XGate, YGate, ZGate, PhaseGate or PhaseInv, for which we have overrides. The error message should
  // catch future mishaps, if, for example, I accidentally call this function directly, instead of as the second
  // dispatch in double-dispatch.
  cerr << "Base class version of rightRMoverChange(const ZRotation& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


unique_ptr<Gate> Gate::leftRMoverChange(const ZRotation& prev) const
{
  // This should never be called, as canRSwap() should only return 'true' when *this a Rotation (including
  // ArbitraryPhase), XGate, YGate, ZGate, PhaseGate or PhaseInv, for which we have overrides. The error message should
  // catch future mishaps, if, for example, I accidentally call this function directly, instead of as the second
  // dispatch in double-dispatch.
  cerr << "Base class version of leftRMoverChange(const ZRotation& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


long Gate::canSimplify(const ArbitraryPhase& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  // Any simplifications are handled by overrides in derived classes, If no override exists, no simplification takes
  // place and this default implementation returns 0, i.e. no improvement.
  return 0;
}


GateSequence Gate::simplification(const ArbitraryPhase& prev) const
{
  // This should never be called, as the fact that an override in one of the derived classes has not been triggered
  // indicates that no simplfication takes place.
  cerr << "Base class version of simplification(const XRotation& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


bool Gate::canRSwap(const ArbitraryPhase& prev) const
{
  return false;
}


unique_ptr<Gate> Gate::rightRMoverChange(const ArbitraryPhase& prev) const
{
  // This should never be called, as canRSwap() should only return 'true' when *this a Rotation (including
  // ArbitraryPhase), XGate, YGate, ZGate, PhaseGate or PhaseInv, for which we have overrides. The error message should
  // catch future mishaps, if, for example, I accidentally call this function directly, instead of as the second
  // dispatch in double-dispatch.
  cerr << "Base class version of rightRMoverChange(const ArbitraryPhase& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


unique_ptr<Gate> Gate::leftRMoverChange(const ArbitraryPhase& prev) const
{
  // This should never be called, as canRSwap() should only return 'true' when *this a Rotation (including
  // ArbitraryPhase), XGate, YGate, ZGate, PhaseGate or PhaseInv, for which we have overrides. The error message should
  // catch future mishaps, if, for example, I accidentally call this function directly, instead of as the second
  // dispatch in double-dispatch.
  cerr << "Base class version of leftRMoverChange(const ArbitraryPhase& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


long Gate::canSimplify(const SwapGate& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  // Any simplifications are handled by overrides in derived classes, If no override exists, no simplification takes
  // place and this default implementation returns 0, i.e. no improvement.
  return 0;
}


GateSequence Gate::simplification(const SwapGate& prev) const
{
  // This should never be called, as the fact that an override in one of the derived classes has not been triggered
  // indicates that no simplfication takes place.
  cerr << "Base class version of simplification(const SwapGate& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


bool Gate::canSwapSwap(const SwapGate& prev) const
{
  // Again, since overrides exist when *this is a SingleTargetGate or a SwapGate, this function only runs when *this is
  // an Oracle, which cannot be swapped with a swap gate.
  return false;
}


unique_ptr<Gate> Gate::rightSwapMoverChange(const SwapGate &prev) const
{
  // This should never be called, as canSwapSwap() should only return 'true' when *this a SingleTargetGate or a
  // SwapGate, for which we have overrides.
  cerr << "Base class version of rightSwapMoverChange(const SwapGate& next) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


unique_ptr<Gate> Gate::leftSwapMoverChange(const SwapGate &prev) const
{
  // This should never be called, as canSwapSwap() should only return 'true' when *this a SingleTargetGate or a
  // SwapGate, for which we have overrides.
  cerr << "Base class version of leftSwapMoverChange(const SwapGate& next) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


bool Gate::canSimplySwap(const Oracle& prev) const
{
  // An override for this function exists when *this is a DiagonalGate. Other gates do not commute with the Oracle.
  // (While the Oracle commutes with the Oracle, there is little point in swapping them!)
  return false;
}


long Gate::canSimplify(const Oracle& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  // Any simplifications are handled by overrides in derived classes, If no override exists, no simplification takes
  // place and this default implementation returns 0, i.e. no improvement.
  return 0;
}


GateSequence Gate::simplification(const Oracle& prev) const
{
  // This should never be called, as the fact that an override in one of the derived classes has not been triggered
  // indicates that no simplfication takes place.
  cerr << "Base class version of simplification(const Oracle& prev) should never run. Aborting." << endl;
  exit(EXIT_FAILURE);
}


int Gate::num_qbits() const
{
  return context_->numQbits();
}


std::ostream& operator<<(std::ostream& out, const Gate& gate)
{
  gate.output(out);
  return out;
}

//----------------------------------------------------------------------------------------------------------------------

SingleTargetGate::SingleTargetGate(const CircuitContext& context) :
Gate(context),
target_(0),  // Ensures that the gate is valid, despite being 'default' constructed.
controls_(num_qbits())
{
}


SingleTargetGate::SingleTargetGate(const CircuitContext& context, int target, const Controls& controls) :
Gate(context),
target_(target),
controls_(controls)
{
}


void SingleTargetGate::createSimulators()
{
  // Create simulators for forward and reverse simulation of this gate and link to them via shared_ptrs. This eradicates
  // the need to repeatedly create the simulator whenever we want to apply the gate (or its inverse) to a state, and
  // should be more efficient.
  //
  // If the shared_ptrs already point somewhere, we assume that the job has already been done. Hence it is important,
  // when changing a SingleTargetGate, to set both forward_simulator and reverse_simulator to nullptr.
  if (!simulator_index_manager)
  {
    vector<int> targets{target()};
    simulator_index_manager = make_shared<SimulatorIndexManager>(targets, controls().asList());
  }

  if (!forward_simulator)
  {
    assert(!reverse_simulator);
    vector<vector<cmplx> > matrix{matrix_()};
    vector<vector<cmplx> > invMatrix{{conj(matrix[0][0]), conj(matrix[1][0])},
                                     {conj(matrix[0][1]), conj(matrix[1][1])}};
    forward_simulator = make_shared<SparseGateSimulator>(sparse(matrix));
    reverse_simulator = make_shared<SparseGateSimulator>(sparse(invMatrix));
  }
}


int SingleTargetGate::numQbitOptions() const
{
  int n = num_qbits();
  int numOptions = 0;
  for (int c = 0; c < n; ++c)  // c is number of controls
  {
    if (context_->gateTypeAvailable(*this, c))
    {
      numOptions += n * binCoeff(n - 1, c);
    }
  }

  return numOptions;
}


void SingleTargetGate::random()
{
  randomize_qbits();

  // Gate has changed - destroy pointers to the old simulators.
  reset_simulator_index_manager();
}


const Controls& SingleTargetGate::controls() const
{
  return controls_;
}


int SingleTargetGate::numControls() const
{
  return controls().numControls();
}


void SingleTargetGate::swapBits(int i, int j)
{
  controls_.swapBits(i, j);
  target_ = swapValues(target_, i, j);

  // Gate has changed - destroy pointers to the old simulators.
  reset_simulator_index_manager();
}


bool SingleTargetGate::available() const
{
  return context_->gateTypeAvailable(*this, numControls());
}


long SingleTargetGate::cost() const
{
  // Calls function in the CircuitContext that, in turn, calls gateTypeId() (which is virtual) to get the type of gate.
  return context_->gateCost(*this, numControls());
}


vector<int> SingleTargetGate::allQbits() const
{
  // Returns a sorted vector of all qbits affected or used by this gate.
  // (Slow - we avoid this function, if possible, whenever efficiency is important.)

  vector<int> qbitSet;
  qbitSet.reserve(num_qbits());
  for (auto qbit = 0; qbit < num_qbits(); ++qbit)
  {
    if (qbit == target_ || controls().isControl(qbit))
    {
      qbitSet.push_back(qbit);
    }
  }
  return qbitSet;
}


bool SingleTargetGate::allAreInvolvedIn(const SingleTargetGate& other) const
{
  // Determines whether the set of all qbits involved in this gate are a subset of those involved in 'other'.
  if (!other.is_involved(target()))
  {
    return false;
  }
  for (auto qbit = 0; qbit < num_qbits(); ++qbit)
  {
    if (controls().isControl(qbit) && !other.is_involved(qbit))
    {
      return false;
    }
  }
  return true;
}


bool SingleTargetGate::allQbitsMatch(const SingleTargetGate& other) const
{
  // Determines whether the sets of qbits involved in the two gates match.
  for (auto qbit = 0; qbit < num_qbits(); ++qbit)
  {
    if (is_involved(qbit) != other.is_involved(qbit))
    {
      return false;
    }
  }
  return true;
}


bool SingleTargetGate::allQbitsMatchControlsOf(const SingleTargetGate& other) const
{
  // Determines whether the set of bits involved in this gate (including the target) match the controls of 'other'.
  for (auto qbit = 0; qbit < num_qbits(); ++qbit)
  {
    if (is_involved(qbit) != other.controls().isControl(qbit))
    {
      return false;
    }
  }
  return true;
}


bool SingleTargetGate::controlsAreControlsOf(const SingleTargetGate& other) const
{
  // (Maybe replace with a member function of Controls called 'contains()'?)
  for (auto qbit = 0; qbit < num_qbits(); ++qbit)
  {
    if (controls().isControl(qbit) && !other.controls().isControl(qbit))
    {
      return false;
    }
  }
  return true;
}


void SingleTargetGate::applyTo(State& state) const
{
  state.transform(*forward_simulator, *simulator_index_manager);
}


void SingleTargetGate::applyInvTo(State& state) const
{
  state.transform(*reverse_simulator, *simulator_index_manager);
}


int SingleTargetGate::target() const
{
  return target_;
}


bool SingleTargetGate::matricesCommute(const SingleTargetGate& other) const
{
  // Default case: Matrices for gates of the same type always commute.
  return typeid(*this) == typeid(other);
}


bool SingleTargetGate::matricesCommute(const DiagonalGate &other) const
{
  // Creating overrides for matricesCommute where one of the parameters is a DiagonalGate is unnecessary, since we only
  // use matricesCommute() to determine whether SingleTargetGates can be swapped and the swap rules for DiagonalGates
  // are overridden without calls to matricesCommute(). (For example, two DiagonalGates may always be swapped.) However,
  // we include these overrides to make the name of the function accurate - if I examine the code in future and spot
  // that I haven't put in that diagonal matrices commute, I may get a bit confused.
  return false;
}


bool SingleTargetGate::matricesCommute(const XTypeGate& other) const
{
  return false;
}


bool SingleTargetGate::matricesCommute(const YTypeGate& other) const
{
  return false;
}


bool SingleTargetGate::matricesAnticommute(const SingleTargetGate& other) const
{
  // Default case: gate matrices don't typically anticommute.
  return false;
}


bool SingleTargetGate::matricesAnticommute(const XGate& other) const
{
  // The XGate matrix does not anticommute with any old gate (unitary) matrix.
  return false;
}


bool SingleTargetGate::matricesAnticommute(const YGate& other) const
{
  return false;
}


bool SingleTargetGate::matricesAnticommute(const ZGate& other) const
{
  return false;
}


void SingleTargetGate::output(std::ostream& out) const
{
  out << name() << " on target qbit " << target_ << controls_ << endl;
}


void SingleTargetGate::extendedOutput(std::ostream& out) const  // We might remove this soon.
{
  output(out);
  out << "Forward simulator:" << endl;
  out << forward_simulator << endl;
  out << "Reverse simulator:" << endl;
  out << reverse_simulator << endl;
}


bool SingleTargetGate::cancelsAtStart() const
{
  // If a control bit is always |0> then the gate never does anything. It may therefore be removed.
  for (int qbit = 0; qbit < num_qbits(); ++qbit)
  {
    if (controls().isControl(qbit) && context_->qbitInputOptions(qbit) == QbitOptions::alwaysZero)
    {
      return true;
    }
  }
  return false;
}


bool SingleTargetGate::canSimplySwap(const Gate& next) const
{
  return next.canSimplySwap(*this);
}


bool SingleTargetGate::canSwapSwap(const Gate& next) const
{
  return next.canSwapSwap(*this);  // Double dispatch.
}


unique_ptr<Gate> SingleTargetGate::rightSwapMoverChange(const Gate &next) const
{
  return next.rightSwapMoverChange(*this);  // Double dispatch.
}


unique_ptr<Gate> SingleTargetGate::leftSwapMoverChange(const Gate &next) const
{
  return next.leftSwapMoverChange(*this);  // Double dispatch
}


bool SingleTargetGate::canSimplySwap(const SingleTargetGate &prev) const
{
  // The general canSwap() function for gates with non-diagonal 2x2 matrices. In general, gates can be swapped if each
  // gate's target is not involved in the other gate. There are two additional cases. If the gates' matrices commute,
  // then the gates can be swapped if both have the same target. If the gates' matrices anticommute, they may be swapped
  // if neither has any controls. We deal with the special cases first.
  if (matricesCommute(prev))
  {
    return !target_is_control(prev);
  }
  if (matricesAnticommute(prev) && numControls() == 0 && prev.numControls() == 0)
  {
    return true;
  }
  return !target_is_involved(prev);
}


bool SingleTargetGate::canSimplySwap(const DiagonalGate& prev) const
{
  // When one of the gates has a diagonal matrix, we find that the gates commute provided the target bit of the
  // non-diagonal gate is not involved in the other gate. If the gates' matrices anticommute, there is an additional
  // chance: the gates can be swapped if they have no controls. (If the matrices anticommute, the diagonal matrix is a
  // multiple of that for the ZGate, while the other gate must have zeroes on the diagonal. If the matrices commute,
  // then the diagonal matrix has to be a multiple of the identity, which does not correspond with any gate.)
  if (matricesAnticommute(prev) && numControls() == 0 && prev.numControls() == 0)
  {
    return true;
  }
  return target() != prev.target() && !prev.controls().isControl(target());
}


bool SingleTargetGate::canSwapSwap(const SwapGate& prev) const
{
  // Can swap the gates, provided we swap the roles of the qbits on the SingleTargetGate appropriately.
  return true;
}


unique_ptr<Gate> SingleTargetGate::rightSwapMoverChange(const SwapGate &prev) const
{
  return {};  // Swap gate is unchanged.
}


unique_ptr<Gate> SingleTargetGate::leftSwapMoverChange(const SwapGate &prev) const
{
  // Swap the bits on the SingleTargetGate as required.
  auto firstGate = clone();
  firstGate->swapBits(prev.bit1(), prev.bit2());
  return firstGate;
}


void SingleTargetGate::reset_simulator_index_manager()
{
  simulator_index_manager.reset();
}


void SingleTargetGate::reset_simulators()
{
  forward_simulator.reset();
  reverse_simulator.reset();
}


bool SingleTargetGate::equivalent_structure(const Gate& rhs) const
{
  // If the code gets here, then rhs must also be of a type derived from SingleTargetGate<Controls>. Indeed both *this
  // and rhs should be of the same type, e.g. Hadamard. Furthermore, both gates should have the same 'context', i.e.
  // num_qbits, cost etc.
  if (typeid(*this) == typeid(rhs))  // Technically unnecessary - types should be the same - but playing safe.
  {
    const SingleTargetGate& other = static_cast<const SingleTargetGate&>(rhs);
    return target_ == other.target_ && controls_ == other.controls_;
  }
  return false;
}


bool SingleTargetGate::equals_(const Gate& rhs) const
{
  return equivalent_structure(rhs);
}


bool SingleTargetGate::sort_before(const Gate& rhs) const
{
  // If the code gets here, then rhs must also be of a type derived from SingleTargetGate<Controls>. Indeed both *this
  // and rhs should be of the same type, e.g. Hadamard. Furthermore, both gates should have the same 'context', i.e.
  // num_qbits, cost, etc.
  const SingleTargetGate& other = static_cast<const SingleTargetGate&>(rhs);
  if (target_ != other.target_)
  {
    return target_ < other.target_;
  }
  return controls_ < other.controls_;
}


void SingleTargetGate::randomize_qbits()
{
  if (numQbitOptions() <= 0)
  {
    throw std::logic_error("SingleTargetGate::randomize_qbits() called for a gate type with no options (or fewer!)"
                           " available.");
  }

  auto n = num_qbits();
  target_ = randInt(0, n);

  // At present, the gate has no control bits. Add a randomly selected set. We first select the 'rank' of the
  // combination of the control bits. Combinations are ranked, starting from 0, with the empty set first, then all the
  // one bit sets and so on.
  auto rank = randInt(0, numQbitOptions() / n);

  // We find the number of control bits, k, by subtracting, from the rank, the number of combinations with zero bits,
  // one bit, etc., until we find that we can no longer do so without the rank becoming negative.
  int numControls = 0;
  while (!context_->gateTypeAvailable(*this, numControls) || binCoeff(n - 1, numControls) <= rank)
  {
    if (context_->gateTypeAvailable(*this, numControls))
    {
      rank -= binCoeff(n - 1, numControls);
    }
    ++numControls;
  }

  // The remaining rank is then interpreted as the standard rank of a k-combination. We use the standard algorithm to
  // unrank.
  vector<int> newControls = combinationUnrank(rank, n - 1, numControls);

  // What we have now is a set of k numbers from 0 to n-2, one of which might be the same as the target. We want to have
  // a set of numbers from 0 to n-2, missing the target bit.
  for (auto& control : newControls)
  {
    if (control >= target_)
    {
      ++control;
    }
  }

  // Transfer to the controls object.
  controls_ = Controls(n, newControls);
}


bool SingleTargetGate::is_involved(int qbit) const
{
  return qbit == target_ || controls_.isControl(qbit);
}


bool SingleTargetGate::target_is_involved(const SingleTargetGate &other) const
{
  return target_ == other.target_ || controls().isControl(other.target_) || other.controls().isControl(target_);
}


bool SingleTargetGate::target_is_control(const SingleTargetGate &other) const
{
  // Any gate commutes with a gate of the same type provided the target bit of one gate is not a control on the other.
  // Alas, we cannot simply write this in SingleTargetGate::canPull(const SingleTargetGate&), since that function is
  // designed for two SingleTargetGates that are not necessarily of the same derived type. We therefore provide this
  // function, to be called by functions such as YGate::canPull(const YGate&). (Could we do something fancy with
  // templates??)
  return controls().isControl(other.target_) || other.controls().isControl(target_);
}

//----------------------------------------------------------------------------------------------------------------------

RotationGate::RotationGate(const CircuitContext& context) :
SingleTargetGate(context),
angle_(0.0)  // Ensures that the gate constructed is valid, despite being 'default' constructed.
{
}


RotationGate::RotationGate(const CircuitContext& context, int target, const Controls& controls, double angle) :
SingleTargetGate(context, target, controls)
{
  // (An alternative would be to allow any angle, but post-process where required. We could also defer all angle fixing
  // to the derived classes.)
  setParameter(0, angle);  // WARNING: Virtual function in constructor will NOT call derived class versions. Angles for
                           //          controlled ArbPhase must be fixed later.
}


int RotationGate::numParameters() const
{
  // RotationGate is the base class for all gates with a single angle parameter.
  return 1;
}


double RotationGate::parameter(int paramNum) const
{
  assert(paramNum == 0);  // Rotation gates have only one parameter - the angle.
  return angle();
}


void RotationGate::setParameter(int paramNum, double value)
{
  // Also ensures that the angle for a gate is in a sensible range. That is usually -2pi to 2pi for a controlled gate
  // and -pi to pi for an uncontrolled gate. (The exception is ArbitraryPhase, which always has a range of -pi to pi.)

  assert(paramNum == 0);  // Rotation gates have only one parameter - the angle.
  if (numControls() == 0)
  {
    angle_ = remainder(value, 2 * pi);  // Numerical optimization might cause the parameter to go 'out of range'. We...
  }                                     // ...fix that here.
  else
  {
    angle_ = remainder(value, 4 * pi);
  }
}


void RotationGate::randomizeParameters()
{
  // Use setParameter(). This will fix the angle for a gate with no controls, if the randomly generated angle is outside
  // the narrower range. (It will also fix the angle for ArbitraryPhase gates too.)

  setParameter(0, randDouble(-2 * pi, 2 * pi));
//  setParameter(0, angle_ + randNormal(0, 0.2));  // Alternative version.
}


void RotationGate::random()
{
  randomize_qbits();
  randomizeParameters();

  // Gate has changed - destroy pointers to the old simulators.
  reset_simulator_index_manager();
  reset_simulators();
}


void RotationGate::mutate()  // Used?
{
  // Mutate the angle parameter for this gate. By using setParameter(), we ensure that the angle is put into the correct
  // range, regardless of gate type or number of controls.
  setParameter(0, randDouble(-2 * pi, 2 * pi));
//  setParameter(0, angle_ + randNormal(0, 0.2));  // Alternative version.

  // Gate has changed - destroy pointers to the old simulators.
  // Simulator manager may remain - it depends only on the identities of the qbits to which the gate is applied.
  reset_simulators();
}


double RotationGate::angle() const
{
  return angle_;
}


bool RotationGate::mutatable() const
{
  return true;
}


void RotationGate::applyPartTo(State &state, int partIndex) const
{
  state.transform(part_simulator(partIndex), *simulator_index_manager);
}


void RotationGate::optimizeAngle(const vector<vector<double> >& errorCoeffs)
{
  // At present, most of the job of optimizing the gate angle is the same regardless of the precise type of rotation
  // gate. Hence this default version. In each case, we need to find the value of theta that minimizes a function of
  // the form
  //          A + B cos(theta) + C sin(theta) + D cos^2(theta) + E^2 sin^2(theta) + F cos(theta) sin(theta)
  // which can be simplified to something of the form
  //          a + b sin(theta + phase1) + c sin(2 theta + phase2).
  // What does differ is the relationship between theta and the Gate's 'angle'. For the X, Y and Z rotations,
  // theta = angle / 2. For the ArbitraryPhase, theta = angle. Hence we have a virtual function to provide this
  // conversion that must be provided by the derived gate classes.

  // We ignore the constant part of the function, i.e. 'a', and refer to the two harmonic parts as the slow and fast
  // parts respectively. (I will want to add a more detailed description of the maths to the code documentation.)
  using constants::pi, constants::tolerance;

  const double& cosCoeff = errorCoeffs[0][1];
  const double& sinCoeff = errorCoeffs[0][2];
  double slowAmp = sqrt(cosCoeff * cosCoeff + sinCoeff * sinCoeff);  // Note that slowAmp is non-negative.
  double slowPhase = atan2(cosCoeff, sinCoeff);

  const double& cosSquaredCoeff = errorCoeffs[1][1];
  const double& sinSquaredCoeff = errorCoeffs[2][2];
  double diff = cosSquaredCoeff - sinSquaredCoeff;  // This subtraction possible causes minor numerical issues.
  const double& sinCosCoeff = errorCoeffs[1][2];
  double fastAmp = sqrt(diff * diff + sinCosCoeff * sinCosCoeff) / 2;  // Note that fastAmp is non-negative.
  double fastPhase = atan2(diff, sinCosCoeff);

  if (slowAmp < tolerance && fastAmp < tolerance)
  {
    // The new gate does nothing to change the error, regardless of angle. Set the angle to zero so that the gate can
    // subsequently be detected and removed by the 'reduction' code.
    setParameter(0, 0.0);
  }
  else if (slowAmp < tolerance)
  {
    // Error is merely a constant plus a fast harmonic compenent. Finding the minimum is easy.
    double theta = -(fastPhase + pi / 2) / 2;
    setParameter(0, thetaToAngle(theta));  // This is the bit that needs to know about gate type.
  }
  else if (fastAmp < tolerance)
  {
    // Error is merely a constant plus a slow harmonic component. Finding the minimum is easy.
    double theta = -(slowPhase + pi / 2);
    setParameter(0, thetaToAngle(theta));  // This is the bit that needs to know about gate type.
  }
  else
  {
    // Error contains both fast and slow harmonic components. We will need numeric optimization.
    // Create the lambdas. (Functions to be passed to the optimizer.)
    // (Currently completely untested! On our current problems, this code is never needed.)
    auto f = [slowAmp, slowPhase, fastAmp, fastPhase](double theta)
    {
      return slowAmp * sin(theta + slowPhase) + fastAmp * sin(2 * theta + fastPhase);
    };
    auto df = [slowAmp, slowPhase, fastAmp, fastPhase](double theta)
    {
      return slowAmp * cos(theta + slowPhase) + 2 * fastAmp * cos(2 * theta + fastPhase);
    };
    auto ddf = [slowAmp, slowPhase, fastAmp, fastPhase](double theta)
    {
      return slowAmp * sin(theta + slowPhase) + 4 * fastAmp * sin(2 * theta + fastPhase);
    };

    // Create our first guess. The minimum must be in the region where the slow component is negative. The fast
    // component has a minimum within this region. We pick a point midway between the minima of the two components.
    // (We could undoubtably improve on this starting point.)
    double slowMin = 3 * pi / 2 - slowPhase;
    double fastMin = (3 * pi / 2 - fastPhase) / 2;  // A minimum, but not necessarily the closest to slowMin.
    while (fastMin < slowMin - pi / 2)
    {
      fastMin += pi;
    }
    while (fastMin > slowMin + pi / 2)
    {
      fastMin -= pi;
    }
    double firstGuess = (slowMin + fastMin) / 2;  // (This is the bit that most likely can be improved.)

    // Create the optimizer.
    Newton newton(f, df, ddf);

    // Optimize.
    newton.minimize(firstGuess, 0.00001);
    double theta = newton.minimum();

    // Set the gate angle to the best value.
    setParameter(0, thetaToAngle(theta));  // This is the bit that needs to know about gate type.
  }

  // The gate's numeric matrix has changed. This means we must change the simulators, since these still reference the
  // old matrix.
  // (This requirement to think about whether the GateSimulator needs updating whenever the Gate's parameter is modified
  // now seems like a symptom of bad structuring. When we added pointers to GateSimulators to the gates and structured
  // things so that the simulators were not created whenever we needed to simulate a Gate, we got faster execution
  // times. However, now most of the stuff has been put into the SimulatorIndexManager. Only the matrix now lives in the
  // simulator. Perhaps we should update the simulator whenever setParameter() is called. Alternatively, we might simply
  // calculate the matrix when simulating the gate.)
  reset_simulators();
  createSimulators();
}


void RotationGate::output(std::ostream& out) const
{
  out << name() << " of angle " << angle_ << " on target qbit " << target_ << controls_ << endl;
}


bool RotationGate::equals_(const Gate& rhs) const
{
  // If the code gets here, then rhs must also be of a type derived from RotationGate<Controls>. Indeed both *this and
  // rhs should be of the same type, e.g. XRotation. Furthermore, both gates should have the same 'context', i.e.
  // num_qbits, costs etc.
  if (typeid(*this) == typeid(rhs))  // Technically unnecessary - types should be the same - but playing safe.
  {
    const RotationGate& other = dynamic_cast<const RotationGate&>(rhs);  // Must use dynamic_cast - using virtual
    return target_ == other.target_ && angle_ == other.angle_ && controls_ == other.controls_;  // ...inheritance.
  }
  return false;
}

bool RotationGate::sort_before(const Gate& rhs) const
{
  // If the code gets here, then rhs must also be of a type derived from RotationGate<Controls>. Indeed both *this and
  // rhs should be of the same type, e.g. XRotation. Furthermore, both gates should have the same 'context', i.e.
  // num_qbits, cost, etc.
  const RotationGate& other = dynamic_cast<const RotationGate&>(rhs);  // Must have dynamic_cast now that we have...
  if (target_ != other.target_)                                        // ...virtual inheritance.
  {
    return target_ < other.target_;
  }
  if (controls_ != other.controls_)
  {
    return controls_ < other.controls_;
  }
  return angle_ < other.angle_;
}


//----------------------------------------------------------------------------------------------------------------------

DiagonalGate::DiagonalGate(const CircuitContext& context) :
SingleTargetGate(context)
{
}


DiagonalGate::DiagonalGate(const CircuitContext& context, int target, const Controls& controls) :
SingleTargetGate(context, target, controls)
{
}


bool DiagonalGate::matricesCommute(const SingleTargetGate &other) const
{
  // Creating overrides for matricesCommute where one of the parameters is a DiagonalGate is unnecessary, since we only
  // use matricesCommute() to determine whether SingleTargetGates can be swapped and the swap rules for DiagonalGates
  // are overridden with calls to matricesCommute(). (For example, two DiagonalGates may always be swapped.) However, we
  // include these overrides to make the name of the function accurate - if I examine the code in future and discover
  // that I haven't put in that diagonal matrices commute, I may get a bit confused.
  return other.matricesCommute(*this);
}


bool DiagonalGate::matricesCommute(const DiagonalGate &other) const
{
  // Creating overrides for matricesCommute where one of the parameters is a DiagonalGate is unnecessary, since we only
  // use matricesCommute() to determine whether SingleTargetGates can be swapped and the swap rules for DiagonalGates
  // are overridden with calls to matricesCommute(). (For example, two DiagonalGates may always be swapped.) However, we
  // include these overrides to make the name of the function accurate - if I examine the code in future and discover
  // that I haven't put in that diagonal matrices commute, I may get a bit confused.
  return true;
}


bool DiagonalGate::cancelsAtStart() const
{
  // Diagonal gates at the beginning of the circuit are redundant whenever SingleTargetGates are. They are also
  // redundant when all the involved qbits have fixed inputs - the gate merely multiplies the state by an irrelevant
  // overall phase.
  if (SingleTargetGate::cancelsAtStart())
  {
    // A control bit must be fixed at |0>.
    return true;
  }
  for (int qbit = 0; qbit < num_qbits(); ++qbit)
  {
    if (controls().isControl(qbit) && context_->qbitInputOptions(qbit) == QbitOptions::varies)
    {
      // No controls fixed at |0> and at least one varies. Hence an input basis state may, or may not be multiplied by a
      // phase. This is physically relevant for non-basis inputs
      return false;
    }
  }
  return context_->qbitInputOptions(target()) != QbitOptions::varies;
}


bool DiagonalGate::canSimplySwap(const Gate& next) const
{
  return next.canSimplySwap(*this);
}


bool DiagonalGate::canSimplySwap(const SingleTargetGate& prev) const
{
  return prev.canSimplySwap(*this);
}


bool DiagonalGate::canSimplySwap(const DiagonalGate &prev) const
{
  return true;
}


bool DiagonalGate::canSimplySwap(const Oracle &prev) const
{
  // Diagonal gates merely mark basis states, by which we mean that basis states are multiplied by a c-number. This is
  // also all that the Oracle does. The order in which we perform the multiplications does not matter.
  return true;
}

//----------------------------------------------------------------------------------------------------------------------

PhaseTypeGate::PhaseTypeGate(const CircuitContext& context) :
DiagonalGate(context)
{
}


PhaseTypeGate::PhaseTypeGate(const CircuitContext& context, int target, const Controls& controls) :
DiagonalGate(context)
{
  // If one of the controls comes before the target qbit, switch them over
  // (Inefficient. The vector of bools in the Controls is converted to a vector of the qbits and then back again.)
  auto qbits = controls.qbits();  // Note that qbits is already sorted.
  if (!qbits.empty() && qbits[0] < target)
  {
    std::swap(qbits[0], target);
  }

  target_ = target;
  controls_ = Controls(num_qbits(), qbits);  // No need to sort 'qbits' - they're stored a vector of bools anyway.
}


PhaseTypeGate::PhaseTypeGate(const CircuitContext& context, const Controls& controls) :  // The 'target' is any one...
DiagonalGate(context)                                                                    // ...of the 'controls'.
{
  // In a controlled PhaseGate, the target bit can be swapped with any control without changing the gate. In a sense, we
  // can consider all bits to be target bits. We pick one to be 'the' target and mark the remainder as controls.

  // Get all the used qbits in a vector.
  auto qbits = controls.qbits();

  // Remove the first qbit from the vector, marking it as the target.
  target_ = qbits.front();
  qbits.front() = qbits.back();
  qbits.pop_back();

  // Set the controls of this gate to be the remaining qbits.
  controls_ = Controls(num_qbits(), qbits);  // No need to sort 'qbits' - they're stored as a vector of bools anyway.
}


int PhaseTypeGate::numQbitOptions() const
{
  auto n = num_qbits();
  int numOptions = 0;
  for (int c = 0; c < num_qbits(); ++c)  // c is number of controls
  {
    if (context_->gateTypeAvailable(*this, c))
    {
      numOptions += binCoeff(n, c + 1);
    }
  }

  return numOptions;
}


void PhaseTypeGate::swapBits(int i, int j)
{
  // While performing the bit swap, we need to ensure that the target bit still comes before the control bits in the
  // result. I.e. if we have a Phase gate on target bit 0, control bit 1, and we swap bits 0 and 2, the result should
  // NOT be target bit 2, control bit 1, but should instead by target bit 1, control bit 2. (One may interchange target
  // and control bits on phase type gates. To avoid duplications, we insist on using a canonical form for such gates,
  // i.e. target comes first.)

  target_ = swapValues(target_, i, j);
  controls_.swapBits(i, j);
  int firstControl = controls_.first();
  if (firstControl < target_)
  {
    controls_.remove(firstControl);
    controls_.add(target_);
    target_ = firstControl;
  }

  // Since the gate has been changed, ensure that the following pointers do not point to the old simulators.
  reset_simulator_index_manager();
}


bool PhaseTypeGate::cancelsAtStart() const
{
  // A Phase type gate at the start of the circuit is redundant whenever a diagonal gate is. It is also redundant when
  // the target qbit is always |0> - the gate does nothing.
  if (DiagonalGate::cancelsAtStart())
  {
    return true;
  }
  return context_->qbitInputOptions(target()) == QbitOptions::alwaysZero;
}



void PhaseTypeGate::randomize_qbits()
{
  // Phase type gates allow the target and a control to be interchanged, without actually changing the gate. To avoid a
  // multiplicity of circuits that look different, but are in fact the same, we force all such gates to have their
  // target bit be the first (least) of the bits involved in the gate. I.e. target_ < all controls.
  
  if (numQbitOptions() <= 0)
  {
    throw std::logic_error("PhaseTypeGate::randomize_qbits() called for a gate type with no options (or fewer!)"
                           " available.");
  }
  auto n = num_qbits();

  // At present, the gate has no target or controls. Add a randomly selected set. We first select the 'rank' of the
  // combination of the involved qbits. These are ranked, starting from 0, with the empty set first, then all the one
  // bit sets and so on.
  auto rank = randInt(0, numQbitOptions());

  // We find the number of involved bits, k, by subtracting, from the rank, the number of combinations with zero bits,
  // one bit, etc., until we find that we can no longer do so without the rank becoming negative.
  int numControls = 0;
  while (!context_->gateTypeAvailable(*this, numControls) || binCoeff(n, numControls + 1) <= rank)
  {
    if (context_->gateTypeAvailable(*this, numControls))
    {
      rank -= binCoeff(n, numControls + 1);
    }
    ++numControls;
  }

  // The remaining rank is then interpreted as the standard rank of a k-combination. We use the standard algorithm to
  // unrank.
  vector<int> involvedQbits = combinationUnrank(rank, n, numControls + 1);

  // Separate the target from the controls.
  target_ = involvedQbits.front();
  involvedQbits.erase(involvedQbits.begin());
  controls_ = Controls(n, involvedQbits);
}

//----------------------------------------------------------------------------------------------------------------------

XTypeGate::XTypeGate(const CircuitContext& context) :
SingleTargetGate(context)
{
}


XTypeGate::XTypeGate(const CircuitContext& context, int target, const Controls& controls) :
SingleTargetGate(context, target, controls)
{
}


bool XTypeGate::matricesCommute(const SingleTargetGate& other) const
{
  return other.matricesCommute(*this);
}


bool XTypeGate::matricesCommute(const XTypeGate& other) const
{
  return true;
}

//----------------------------------------------------------------------------------------------------------------------

YTypeGate::YTypeGate(const CircuitContext& context) :
SingleTargetGate(context)
{
}


YTypeGate::YTypeGate(const CircuitContext& context, int target, const Controls& controls) :
SingleTargetGate(context, target, controls)
{
}


bool YTypeGate::matricesCommute(const SingleTargetGate& other) const
{
  return other.matricesCommute(*this);
}


bool YTypeGate::matricesCommute(const YTypeGate& other) const
{
  return true;
}

//----------------------------------------------------------------------------------------------------------------------

// (While each gate doesn't require that much, this still seems like 'code repetition' to me. Could we make an array of
// name/matrix pairs and then make a class templated on the array index, so that we can reduce the repetition?)

Hadamard::Hadamard(const CircuitContext& context) :
SingleTargetGate(context)
{
//  // Gate has not been initialized, so no need to call 'create_simulators()' yet.
}


Hadamard::Hadamard(const CircuitContext& context, int target, const Controls& controls) :
SingleTargetGate(context, target, controls)
{
}


int Hadamard::gateTypeId() const
{
  return gate_type_id;
}


unique_ptr<Gate> Hadamard::clone() const
{
  return make_unique<Hadamard>(*this);  // Will get converted to unique_ptr<Gate>.
}


Hadamard* Hadamard::rawClone() const
{
  return new Hadamard(*this);
}


std::string Hadamard::name() const
{
  return "Hadamard";
}


bool Hadamard::gateAvailable(const CircuitContext& context, int target, const Controls& controls)  // Static
{
  // Static function for determining the availability of a Gate without needing to create one. (The corresponding
  // non-static function is useful when we have a Gate but don't necessarily know the precise gate type.)
  // (We don't need to provide the target bit at present - just the number of controls. However, in future we may have
  // more complicated rules on availability.)
  return context.gateTypeAvailable(gate_type_id, controls.numControls());
}


long Hadamard::gateCost(const CircuitContext& context, int numControls)  // Static
{
  return context.gateCost(gate_type_id, numControls);
}


long Hadamard::canSimplify(const Gate& next) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  return next.canSimplify(*this);  // Double dispatch.
}


GateSequence Hadamard::simplification(const Gate& next) const
{
  return next.simplification(*this);  // Double dispatch.
}


bool Hadamard::canHSwap(const Gate& next) const
{
  return next.canHSwap(*this);  // Double dispatch.
}


unique_ptr<Gate> Hadamard::rightHMoverChange(const Gate& next) const
{
  return next.rightHMoverChange(*this);  // Double dispatch.
}


unique_ptr<Gate> Hadamard::leftHMoverChange(const Gate& next) const
{
  return next.leftHMoverChange(*this);  // Double dispatch.
}


long Hadamard::canSimplify(const Hadamard& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  if (target() == prev.target() && controls() == prev.controls())
  {
    return cost() + prev.cost();
  }
  return 0;
}


GateSequence Hadamard::simplification(const Hadamard& prev) const
{
  return {};  // Gates cancel.
}


long Hadamard::canSimplify(const XGate& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  // Note that we don't actually check that the replacement is cheaper here - we simply return the improvement. The
  // calling function will filter out negative improvements.
  if (target() == prev.target() && controlsAreControlsOf(prev) &&
      zGateEquivalentAvailable(*context_, prev.target(), prev.controls()))
  {
    return prev.cost() - zGateEquivalentCost(*context_, prev.target(), prev.controls());
  }
  return 0;
}


GateSequence Hadamard::simplification(const XGate& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  return makeGateSequence(clone(), zGateEquivalent(*context_, prev.target(), prev.controls()));
//  return {clone(), zGateEquivalent(*context_, prev.target(), prev.controls())};
}


bool Hadamard::canHSwap(const XGate& prev) const
{
  // The case considered here is where the XGate and the Hadamard have the same target qbit. Cases where the target qbit
  // differ are handled in canSimplySwap(). For the considered case, the transformation is possible whenever the
  // Hadamard's control bits are also control bits of the XGate, and when a ZGate is available. However, it is only
  // considered to be a HSwap if costs remain the same, i.e. when the ZGate is available and has the same cost as the
  // XGate.
  // (Should we allow the use of an ArbitraryPhase instead of a ZGate. Note that, at present, we would not be able to
  // undo such a swap.)
  if (target() == prev.target() && controlsAreControlsOf(prev) &&
      ZGate::gateAvailable(*context_, prev.target(), prev.controls()))
  {
    return ZGate::gateCost(*context_, numControls()) == prev.cost();
  }
  return false;
}


unique_ptr<Gate> Hadamard::rightHMoverChange(const XGate& prev) const
{
  // This function is only called if canHSwap returns true. Hence there is no need to check whether the gates can be
  // swapped, or whether the swap changes the XGate. We already know that it does.
  return make_unique<ZGate>(*context_, prev.target(), prev.controls());
}


unique_ptr<Gate> Hadamard::leftHMoverChange(const XGate& prev) const
{
  return {};  // Hadamard gate is unchanged
}


bool Hadamard::canHSwap(const YGate& prev) const
{
  // Swapping a YGate with a Hadamard results in a change of phase. Provided neither gate has controls, this just means
  // that the overall state is multiplied by -1. We can ignore this. (If, in future, we were to add a MinusYGate, then
  // we could also include cases with controls.)
  return target() == prev.target() && numControls() == 0 && prev.numControls() == 0;
}


unique_ptr<Gate> Hadamard::rightHMoverChange(const YGate& prev) const
{
  // Swapping uncontrolled Hadamards and YGates results in no change to the gates. (As such, we did consider adding
  // these scenarios into the basic canSwap() function. However, with the addition of HSwaps for YRotations too, where
  // the YRotation does change, it felt better to have these as HSwaps.
  return {};
}


unique_ptr<Gate> Hadamard::leftHMoverChange(const YGate& prev) const
{
  return {};  // Hadamard gate is unchanged
}


long Hadamard::canSimplify(const ZGate& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  // (Note that the availability check for the XGate does not use the correct target and controls. However, the number
  // of controls is correct and this is all that matters, at present, for the check.)
  if (allAreInvolvedIn(prev) && XGate::gateAvailable(*context_, prev.target(), prev.controls()))
  {
    // We don't check for cost improvement here. Calling function will deal with negative improvments.
    return prev.cost() - XGate::gateCost(*context_, prev.numControls());
  }
  return 0;
}


GateSequence Hadamard::simplification(const ZGate& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  Controls xControls(prev.controls());  // New gate's controls will be the set of all bits involved in the old Z,...
  xControls.add(prev.target());         // ...minus the new gate's target bit, which happens to be the target of the...
  xControls.remove(target());           // ...Hadamard.

  return makeGateSequence(clone(), make_unique<XGate>(*context_, target(), xControls));
//  return {clone(), make_unique<XGate>(*context_, target(), xControls)};
}


bool Hadamard::canHSwap(const ZGate& prev) const
{
  // The case considered here is where all the qbits involved in the Hadamard are also involved in the ZGate. Then the
  // target of the ZGate can be changed to match that of the Hadamard before converting it into an XGate. It is also
  // necessary for the XGate and ZGate to have the same cost - the case where X is cheaper is found in the
  // 'simplifications'. The case where the Hadamard's target is not involved in the ZGate is found in canSimplySwap().
  // (Note that the check on the availability of the XGate hasn't quite got the right target and controls - the target
  // might need to be switched with one of the controls. However, at present only the number of controls matters when
  // checking availability, so this does not matter.)
  return allAreInvolvedIn(prev) && XGate::gateAvailable(*context_, prev.target(), prev.controls()) &&
         XGate::gateCost(*context_, prev.numControls()) == prev.cost();
}


unique_ptr<Gate> Hadamard::rightHMoverChange(const ZGate& prev) const
{
  // This function is only called if canHSwap returns true. Hence there is no need to check whether the gates can be
  // swapped, or whether the swap changes the ZGate. We already know that it does.
  Controls xControls(prev.controls());  // New gate's controls will be the set of all qbits involved in the old Z,...
  xControls.add(prev.target());         // ...minus the new gate's target qbit, which is the target of the Hadamard.
  xControls.remove(target());

  return make_unique<XGate>(*context_, target(), xControls);
}


unique_ptr<Gate> Hadamard::leftHMoverChange(const ZGate& prev) const
{
  return {};  // Hadamard gate is unchanged
}


long Hadamard::canSimplify(const XRotation& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  if (target() == prev.target() && controlsAreControlsOf(prev) &&
      ZRotation::gateAvailable(*context_, prev.target(), prev.controls()))
  {
    // We don't check that a cost improvement occurs here. Calling function will have to deal with negative improvments.
    return prev.cost() - ZRotation::gateCost(*context_, prev.numControls());
  }
  return 0;
}


GateSequence Hadamard::simplification(const XRotation& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  return makeGateSequence(clone(), make_unique<ZRotation>(*context_, prev.target(), prev.controls(), prev.angle()));
//  return {clone(), make_unique<ZRotation>(*context_, prev.target(), prev.controls(), prev.angle())};
}


bool Hadamard::canHSwap(const XRotation& prev) const
{
  return target() == prev.target() && controlsAreControlsOf(prev) &&
         ZRotation::gateAvailable(*context_, prev.target(), prev.controls()) &&
         ZRotation::gateCost(*context_, prev.numControls()) == prev.cost();
}


unique_ptr<Gate> Hadamard::rightHMoverChange(const XRotation& prev) const
{
  // This function is only called if canHSwap returns true. Hence there is no need to check whether the gates can be
  // swapped, or whether the swap changes the XRotation. We already know that it does.
  return make_unique<ZRotation>(*context_, prev.target(), prev.controls(), prev.angle());
}


unique_ptr<Gate> Hadamard::leftHMoverChange(const XRotation& prev) const
{
  return {};  // Hadamard gate is unchanged
}


bool Hadamard::canHSwap(const YRotation& prev) const
{
  return target() == prev.target() && controlsAreControlsOf(prev);
}


unique_ptr<Gate> Hadamard::rightHMoverChange(const YRotation& prev) const
{
  // This function is only called if canHSwap returns true. Hence there is no need to check whether the gates can be
  // swapped, or whether the swap changes the XRotation. We already know that it does. The YRotation gets its angle
  // flipped.
  return make_unique<YRotation>(*context_, prev.target(), prev.controls(), -prev.angle());
}


unique_ptr<Gate> Hadamard::leftHMoverChange(const YRotation& prev) const
{
  return {};  // Hadamard gate is unchanged
}


long Hadamard::canSimplify(const ZRotation& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  if (target() == prev.target() && controlsAreControlsOf(prev) &&
      XRotation::gateAvailable(*context_, prev.target(), prev.controls()))
  {
    // We don't check that a cost improvement occurs here. Calling function will have to deal with negative improvments.
    return prev.cost() - XRotation::gateCost(*context_, prev.numControls());
  }
  return 0;
}


GateSequence Hadamard::simplification(const ZRotation& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  return makeGateSequence(clone(), make_unique<XRotation>(*context_, prev.target(), prev.controls(), prev.angle()));
//  return {clone(), make_unique<XRotation>(*context_, prev.target(), prev.controls(), prev.angle())};
}


bool Hadamard::canHSwap(const ZRotation& prev) const
{
  return target() == prev.target() && controlsAreControlsOf(prev) &&
         XRotation::gateAvailable(*context_, prev.target(), prev.controls()) &&
         XRotation::gateCost(*context_, prev.numControls()) == prev.cost();
}


unique_ptr<Gate> Hadamard::rightHMoverChange(const ZRotation& prev) const
{
  // This function is only called if canHSwap returns true. Hence there is no need to check whether the gates can be
  // swapped, or whether the swap changes the ZRotation. We already know that it does.
  return make_unique<XRotation>(*context_, prev.target(), prev.controls(), prev.angle());
}


unique_ptr<Gate> Hadamard::leftHMoverChange(const ZRotation& prev) const
{
  return {};  // Hadamard gate is unchanged
}


vector<vector<cmplx> > Hadamard::matrix_() const
{
  using std::sqrt;
  return {{1 / sqrt(2), 1 / sqrt(2)}, {1 / sqrt(2), -1 / sqrt(2)}};
}

//----------------------------------------------------------------------------------------------------------------------

PiByEight::PiByEight(const CircuitContext& context) :
SingleTargetGate(context),
//DiagonalGate(context),
PhaseTypeGate(context)
{
}


PiByEight::PiByEight(const CircuitContext& context, int target, const Controls& controls) :
SingleTargetGate(context, target, controls),
//DiagonalGate(context, target, controls),
PhaseTypeGate(context, target, controls)
{
}


PiByEight::PiByEight(const CircuitContext& context, const Controls& controls) :  // Target is any one of the controls.
SingleTargetGate(context),
//DiagonalGate(context),
PhaseTypeGate(context, controls)
{
}


int PiByEight::gateTypeId() const
{
  return gate_type_id;
}


unique_ptr<Gate> PiByEight::clone() const
{
  return make_unique<PiByEight>(*this);  // Will get converted to unique_ptr<Gate>.
}


PiByEight* PiByEight::rawClone() const
{
  return new PiByEight(*this);
}


std::string PiByEight::name() const
{
  return "PiByEight";
}


bool PiByEight::gateAvailable(const CircuitContext& context, int target, const Controls& controls)  // Static
{
  // Static function for determining the availability of a Gate without needing to create one. (The corresponding
  // non-static function is useful when we have a Gate but don't necessarily know the precise gate type.)
  // (We don't need to provide the target bit at present - just the number of controls. However, in future we may have
  // more complicated rules on availability.)
  return context.gateTypeAvailable(gate_type_id, controls.numControls());
}


bool PiByEight::gateAvailable(const CircuitContext& context, const Controls& controls)  // Static
{
  // Static function for determining the availability of a Gate without needing to create one. (The corresponding
  // non-static function is useful when we have a Gate but don't necessarily know the precise gate type.)
  // Parameter 'controls' is actually the controls AND the target. Hence the number of actual controls for the desired
  // gate is numControls() - 1.
  assert(controls.numControls() > 0);
  return context.gateTypeAvailable(gate_type_id, controls.numControls() - 1);
}


long PiByEight::gateCost(const CircuitContext& context, int numControls)  // Static
{
  return context.gateCost(gate_type_id, numControls);
}


long PiByEight::canSimplify(const Gate& next) const
{
  return next.canSimplify(*this);  // Double dispatch.
}


GateSequence PiByEight::simplification(const Gate& next) const
{
  return next.simplification(*this);  // Double dispatch.
}


long PiByEight::canSimplify(const PiByEight& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  if (target() == prev.target() && controls() == prev.controls() &&
      phaseGateEquivalentAvailable(*context_, target(), controls()))
  {
    return cost() + prev.cost() - phaseGateEquivalentCost(*context_, target(), controls());
  }
  return 0;
}


GateSequence PiByEight::simplification(const PiByEight& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  // Gates become a single PhaseGate or an equivalent ArbitraryPhase.
  return makeGateSequence(phaseGateEquivalent(*context_, target(), controls()));
//  return {phaseGateEquivalent(*context_, target(), controls())};
}


long PiByEight::canSimplify(const PiByEightInv& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  if (target() == prev.target() && controls() == prev.controls())
  {
    return cost() + prev.cost();
  }
  return 0;
}


GateSequence PiByEight::simplification(const PiByEightInv& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  return {};  // Gates cancel.
}


long PiByEight::canSimplify(const PhaseGate& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  if (target() == prev.target() && controls() == prev.controls() &&
      ArbitraryPhase::gateAvailable(*context_, target(), controls()))
  {
    return cost() + prev.cost() - ArbitraryPhase::gateCost(*context_, numControls());
  }
  return 0;
}


GateSequence PiByEight::simplification(const PhaseGate& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  // Gates become a single ArbitraryPhase.
  return makeGateSequence(make_unique<ArbitraryPhase>(*context_, target(), controls(), 3 * pi / 4));
//  return {make_unique<ArbitraryPhase>(*context_, target(), controls(), 3 * pi / 4)};
}


long PiByEight::canSimplify(const PhaseInv& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  if (target() == prev.target() && controls() == prev.controls() &&
      piByEightInvEquivalentAvailable(*context_, target(), controls()))
  {
    return cost() + prev.cost() - piByEightInvEquivalentCost(*context_, target(), controls());
  }
  return 0;
}


GateSequence PiByEight::simplification(const PhaseInv& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  // Gates become a PhaseGate or an equivalent ArbitraryPhase.
  return makeGateSequence(piByEightInvEquivalent(*context_, target(), controls()));
//  return {piByEightInvEquivalent(*context_, target(), controls())};
}


long PiByEight::canSimplify(const ZGate& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  if (target() == prev.target() && controls() == prev.controls() &&
      ArbitraryPhase::gateAvailable(*context_, target(), controls()))
  {
    return cost() + prev.cost() - ArbitraryPhase::gateCost(*context_, numControls());
  }
  return 0;
}


GateSequence PiByEight::simplification(const ZGate& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  // Gates become a single ArbitraryPhase.
  return makeGateSequence(make_unique<ArbitraryPhase>(*context_, target(), controls(), -3 * pi / 4));
//  return {make_unique<ArbitraryPhase>(*context_, target(), controls(), -3 * pi / 4)};
}


long PiByEight::canSimplify(const ZRotation& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  //
  // If both gates have no controls then they simplify to a single ZRotation, up to an overall phase which we can
  // ignore. If both have controls and the controls are the same, then the phase difference becomes important because it
  // only applies when the control bits are set. Hence, in this case, the gates simplify to a ZRotation (with the same
  // controls) and an ArbitraryPhase on the control bits.
  //
  // Note that while one may freely swap the target qbits with a control qbits on a PiByEight, with no effect on its
  // functions, this is not the case for the ZRotation.
  if (allQbitsMatch(prev))
  {
    if (numControls() == 0)
    {
      return cost();
    }
    else if (ArbitraryPhase::gateAvailable(*context_, prev.controls()))
    {
      return cost() - ArbitraryPhase::gateCost(*context_, prev.numControls() - 1);
    }
  }
  return 0;
}


GateSequence PiByEight::simplification(const ZRotation& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  //
  // If both gates have no controls then they simplify to a single ZRotation, up to an overall phase which we can
  // ignore. If both have controls and the controls are the same, then the phase difference becomes important because it
  // only applies when the control bits are set. Hence, in this case, the gates simplify to a ZRotation (with the same
  // controls) and an ArbitraryPhase on the control bits.
  auto firstGate = make_unique<ZRotation>(*context_, prev.target(), prev.controls(), prev.angle() - pi / 4);
  if (numControls() == 0)
  {
    return makeGateSequence(std::move(firstGate));
//    return {std::move(firstGate)};
  }
  return makeGateSequence(std::move(firstGate), make_unique<ArbitraryPhase>(*context_, prev.controls(), pi / 8));
//  return {std::move(firstGate), make_unique<ArbitraryPhase>(*context_, prev.controls(), pi / 8)};
}


long PiByEight::canSimplify(const ArbitraryPhase& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  if (target() == prev.target() && controls() == prev.controls())
  {
    return cost();
  }
  return 0;
}


GateSequence PiByEight::simplification(const ArbitraryPhase& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  return makeGateSequence(make_unique<ArbitraryPhase>(*context_, target(), controls(), prev.angle() + pi / 4));
//  return {make_unique<ArbitraryPhase>(*context_, target(), controls(), prev.angle() + pi / 4)};
}


vector<vector<cmplx> > PiByEight::matrix_() const
{
  using std::exp, constants::i;
  return {{1, 0}, {0, exp(i * pi / 4.0)}};
}

//----------------------------------------------------------------------------------------------------------------------

PiByEightInv::PiByEightInv(const CircuitContext& context) :
SingleTargetGate(context),
//DiagonalGate(context),
PhaseTypeGate(context)
{
}


PiByEightInv::PiByEightInv(const CircuitContext& context, int target, const Controls& controls) :
SingleTargetGate(context, target, controls),
//DiagonalGate(context, target, controls),
PhaseTypeGate(context, target, controls)
{
}


PiByEightInv::PiByEightInv(const CircuitContext& context, const Controls& controls) :  // The target is any one of...
SingleTargetGate(context),                                                             // ...the controls.
//DiagonalGate(context),
PhaseTypeGate(context, controls)
{
}


int PiByEightInv::gateTypeId() const
{
  return gate_type_id;
}


unique_ptr<Gate> PiByEightInv::clone() const
{
  return make_unique<PiByEightInv>(*this);  // Will get converted to unique_ptr<Gate>.
}


PiByEightInv* PiByEightInv::rawClone() const
{
  return new PiByEightInv(*this);
}


std::string PiByEightInv::name() const
{
  return "PiByEightInv";
}


bool PiByEightInv::gateAvailable(const CircuitContext& context, int target, const Controls& controls)  // Static
{
  // Static function for determining the availability of a Gate without needing to create one. (The corresponding
  // non-static function is useful when we have a Gate but don't necessarily know the precise gate type.)
  // (We don't need to provide the target bit at present - just the number of controls. However, in future we may have
  // more complicated rules on availability.)
  return context.gateTypeAvailable(gate_type_id, controls.numControls());
}


bool PiByEightInv::gateAvailable(const CircuitContext& context, const Controls& controls)  // Static
{
  // Static function for determining the availability of a Gate without needing to create one. (The corresponding
  // non-static function is useful when we have a Gate but don't necessarily know the precise gate type.)
  // Parameter 'controls' is actually the controls AND the target. Hence the number of actual controls for the desired
  // gate is numControls() - 1.
  assert(controls.numControls() > 0);
  return context.gateTypeAvailable(gate_type_id, controls.numControls() - 1);
}


long PiByEightInv::gateCost(const CircuitContext& context, int numControls)  // Static
{
  return context.gateCost(gate_type_id, numControls);
}


long PiByEightInv::canSimplify(const Gate& next) const
{
  return next.canSimplify(*this);  // Double dispatch.
}


GateSequence PiByEightInv::simplification(const Gate& next) const
{
  return next.simplification(*this);  // Double dispatch.
}


long PiByEightInv::canSimplify(const PiByEight& prev) const
{
  return prev.canSimplify(*this);  // Simplification is exactly the same as if the gates were switched.
}


GateSequence PiByEightInv::simplification(const PiByEight& prev) const
{
//  return prev.simplification(*this);  // This would work too.
  return {};  // Gates cancel.
}


long PiByEightInv::canSimplify(const PiByEightInv& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  if (target() == prev.target() && controls() == prev.controls() &&
      phaseInvEquivalentAvailable(*context_, target(), controls()))
  {
    return cost() + prev.cost() - phaseInvEquivalentCost(*context_, target(), controls());
  }
  return 0;
}


GateSequence PiByEightInv::simplification(const PiByEightInv& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  // Gates become a single PhaseInv or an equivalent ArbitraryPhase.
  return makeGateSequence(phaseInvEquivalent(*context_, target(), controls()));
//  return {phaseInvEquivalent(*context_, target(), controls())};
}


long PiByEightInv::canSimplify(const PhaseGate& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  if (target() == prev.target() && controls() == prev.controls() &&
      piByEightEquivalentAvailable(*context_, target(), controls()))
  {
    return cost() + prev.cost() - piByEightEquivalentCost(*context_, target(), controls());
  }
  return 0;
}


GateSequence PiByEightInv::simplification(const PhaseGate& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  // Gates becomes a single PiByEight or an equivalant ArbitraryPhase.
  return makeGateSequence(piByEightEquivalent(*context_, target(), controls()));
//  return {piByEightEquivalent(*context_, target(), controls())};
}


long PiByEightInv::canSimplify(const PhaseInv& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  if (target() == prev.target() && controls() == prev.controls() &&
      ArbitraryPhase::gateAvailable(*context_, target(), controls()))
  {
    return cost() + prev.cost() - ArbitraryPhase::gateCost(*context_, numControls());
  }
  return 0;
}


GateSequence PiByEightInv::simplification(const PhaseInv& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  // Gates become a single ArbitraryPhase.
  return makeGateSequence(make_unique<ArbitraryPhase>(*context_, target(), controls(), -3 * pi / 4));
//  return {make_unique<ArbitraryPhase>(*context_, target(), controls(), -3 * pi / 4)};
}


long PiByEightInv::canSimplify(const ZGate& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  if (target() == prev.target() && controls() == prev.controls() &&
      ArbitraryPhase::gateAvailable(*context_, target(), controls()))
  {
    return cost() + prev.cost() - ArbitraryPhase::gateCost(*context_, numControls());
  }
  return 0;
}


GateSequence PiByEightInv::simplification(const ZGate& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  // Gates become a single ArbitraryPhase.
  return makeGateSequence(make_unique<ArbitraryPhase>(*context_, target(), controls(), 3 * pi / 4));
//  return {make_unique<ArbitraryPhase>(*context_, target(), controls(), 3 * pi / 4)};
}


long PiByEightInv::canSimplify(const ZRotation& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  //
  // If both gates have no controls then they simplify to a single ZRotation, up to an overall phase which we can
  // ignore. If both have controls and the controls are the same, then the phase difference becomes important because it
  // only applies when the control bits are set. Hence, in this case, the gates simplify to a ZRotation (with the same
  // controls) and an ArbitraryPhase on the control bits.
  //
  // Note that while one may freely swap the target qbits with a control qbits on a PiByEight, with no effect on its
  // functions, this is not the case for the ZRotation.
  if (allQbitsMatch(prev))
  {
    if (numControls() == 0)
    {
      return cost();
    }
    else if (ArbitraryPhase::gateAvailable(*context_, prev.controls()))
    {
      return cost() - ArbitraryPhase::gateCost(*context_, prev.numControls() - 1);
    }
  }
  return 0;
}


GateSequence PiByEightInv::simplification(const ZRotation& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  //
  // If both gates have no controls then they simplify to a single ZRotation, up to an overall phase which we can
  // ignore. If both have controls and the controls are the same, then the phase difference becomes important because it
  // only applies when the control bits are set. Hence, in this case, the gates simplify to a ZRotation (with the same
  // controls) and an ArbitraryPhase on the control bits.
  auto firstGate = make_unique<ZRotation>(*context_, prev.target(), prev.controls(), prev.angle() + pi / 4);
  if (numControls() == 0)
  {
    return makeGateSequence(std::move(firstGate));
//    return {std::move(firstGate)};
  }
  return makeGateSequence(std::move(firstGate), make_unique<ArbitraryPhase>(*context_, prev.controls(), -pi / 8));
//  return {std::move(firstGate), make_unique<ArbitraryPhase>(*context_, prev.controls(), -pi / 8)};
}


long PiByEightInv::canSimplify(const ArbitraryPhase& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  if (target() == prev.target() && controls() == prev.controls())
  {
    return cost();
  }
  return 0;
}


GateSequence PiByEightInv::simplification(const ArbitraryPhase& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  return makeGateSequence(make_unique<ArbitraryPhase>(*context_, target(), controls(), prev.angle() - pi / 4));
//  return {make_unique<ArbitraryPhase>(*context_, target(), controls(), prev.angle() - pi / 4)};
}


vector<vector<cmplx> > PiByEightInv::matrix_() const
{
  using std::exp, constants::i;
  return {{1, 0}, {0, exp(-i * pi / 4.0)}};
}

//----------------------------------------------------------------------------------------------------------------------

PhaseGate::PhaseGate(const CircuitContext& context) :
SingleTargetGate(context),
//DiagonalGate(context),
PhaseTypeGate(context)
{
}


PhaseGate::PhaseGate(const CircuitContext& context, int target, const Controls& controls) :
SingleTargetGate(context, target, controls),
//DiagonalGate(context, target, controls),
PhaseTypeGate(context, target, controls)
{
}


PhaseGate::PhaseGate(const CircuitContext& context, const Controls& controls) :  // Target is any one of the controls.
SingleTargetGate(context),
//DiagonalGate(context),
PhaseTypeGate(context, controls)
{
}


int PhaseGate::gateTypeId() const
{
  return gate_type_id;
}


unique_ptr<Gate> PhaseGate::clone() const
{
  return make_unique<PhaseGate>(*this);  // Will get converted to unique_ptr<Gate>.
}


PhaseGate* PhaseGate::rawClone() const
{
  return new PhaseGate(*this);
}


std::string PhaseGate::name() const
{
  return "Phase";
}


bool PhaseGate::gateAvailable(const CircuitContext& context, int target, const Controls& controls)  // Static
{
  // Static function for determining the availability of a Gate without needing to create one. (The corresponding
  // non-static function is useful when we have a Gate but don't necessarily know the precise gate type.)
  // (We don't need to provide the target bit at present - just the number of controls. However, in future we may have
  // more complicated rules on availability.)
  return context.gateTypeAvailable(gate_type_id, controls.numControls());
}


bool PhaseGate::gateAvailable(const CircuitContext& context, const Controls& controls)  // Static
{
  // Static function for determining the availability of a Gate without needing to create one. (The corresponding
  // non-static function is useful when we have a Gate but don't necessarily know the precise gate type.)
  // Parameter 'controls' is actually the controls AND the target. Hence the number of actual controls for the desired
  // gate is numControls() - 1.
  assert(controls.numControls() > 0);
  return context.gateTypeAvailable(gate_type_id, controls.numControls() - 1);
}


long PhaseGate::gateCost(const CircuitContext& context, int numControls)  // Static
{
  return context.gateCost(gate_type_id, numControls);
}


long PhaseGate::canSimplify(const Gate& next) const
{
  return next.canSimplify(*this);  // Double dispatch.
}


GateSequence PhaseGate::simplification(const Gate& next) const
{
  return next.simplification(*this);  // Double dispatch.
}


bool PhaseGate::canRSwap(const Gate& next) const
{
  return next.canRSwap(*this);  // Double dispatch.
}


unique_ptr<Gate> PhaseGate::rightRMoverChange(const Gate& next) const
{
  return next.rightRMoverChange(*this);  // Double dispatch.
}


unique_ptr<Gate> PhaseGate::leftRMoverChange(const Gate& next) const
{
  return next.leftRMoverChange(*this);  // Double dispatch.
}


long PhaseGate::canSimplify(const PiByEight& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.canSimplify(*this);
}


GateSequence PhaseGate::simplification(const PiByEight& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.simplification(*this);
}


long PhaseGate::canSimplify(const PiByEightInv& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.canSimplify(*this);
}


GateSequence PhaseGate::simplification(const PiByEightInv& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.simplification(*this);
}


long PhaseGate::canSimplify(const PhaseGate& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  if (target() == prev.target() && controls() == prev.controls() &&
      zGateEquivalentAvailable(*context_, target(), controls()))
  {
    return cost() + prev.cost() - zGateEquivalentCost(*context_, target(), controls());
  }
  return 0;
}


GateSequence PhaseGate::simplification(const PhaseGate& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  // Gates become a single ZGate or an equivalent ArbitraryPhase.
  return makeGateSequence(zGateEquivalent(*context_, target(), controls()));
//  return {zGateEquivalent(*context_, target(), controls())};
}


long PhaseGate::canSimplify(const PhaseInv& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  if (target() == prev.target() && controls() == prev.controls())
  {
    return cost() + prev.cost();
  }
  return 0;
}


GateSequence PhaseGate::simplification(const PhaseInv& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  return {};  // Gates cancel.
}


long PhaseGate::canSimplify(const ZGate& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  if (target() == prev.target() && controls() == prev.controls() &&
      phaseInvEquivalentAvailable(*context_, target(), controls()))
  {
    return cost() + prev.cost() - phaseInvEquivalentCost(*context_, target(), controls());
  }
  return 0;
}


GateSequence PhaseGate::simplification(const ZGate& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  // Gates become a single PhaseInv or an equivalent ArbitraryPhase.
  return makeGateSequence(phaseInvEquivalent(*context_, target(), controls()));
//  return {phaseInvEquivalent(*context_, target(), controls())};
}


long PhaseGate::canSimplify(const XRotation& prev) const
{
  // If a YRotation is cheaper than an XRotation (which seems unlikely), then we can swap the gates over, turning the
  // XRotation into a YRotation of the same angle, provided targets match and all controls of the PhaseGate are also
  // controls of the rotation.
  if (YRotation::gateCost(*context_, prev.numControls()) < prev.cost() && target() == prev.target() &&
      controlsAreControlsOf(prev))
  {
    return prev.cost() - YRotation::gateCost(*context_, prev.numControls());
  }
  return 0;
}


GateSequence PhaseGate::simplification(const XRotation& prev) const
{
  return makeGateSequence(clone(), make_unique<YRotation>(*context_, prev.target(), prev.controls(), prev.angle()));
//  return {clone(), make_unique<YRotation>(*context_, prev.target(), prev.controls(), prev.angle())};
}


bool PhaseGate::canRSwap(const XRotation& prev) const
{
  // Gates can be swapped, with the XRotation becoming a YRotation of the same angle, provided targets match and all
  // controls of the PhaseGate are also controls of the rotation. We must also check that the costs of XRotations and
  // YRotations match.
  return target() == prev.target() && controlsAreControlsOf(prev) &&
         YRotation::gateCost(*context_, prev.numControls()) == prev.cost();
}


unique_ptr<Gate> PhaseGate::rightRMoverChange(const XRotation& prev) const
{
  // XRotation turns into a YRotation of the same angle.
  return make_unique<YRotation>(*context_, prev.target(), prev.controls(), prev.angle());
}


unique_ptr<Gate> PhaseGate::leftRMoverChange(const XRotation& prev) const
{
  return {};  // PhaseGate is unchanged.
}


long PhaseGate::canSimplify(const YRotation& prev) const
{
  // If an XRotation is cheaper than a YRotation (which seems unlikely), then we can swap the gates over, turning the
  // YRotation into an XRotation of opposite angle, provided targets match and all controls of the PhaseGate are also
  // controls of the rotation.
  if (XRotation::gateCost(*context_, prev.numControls()) < prev.cost() && target() == prev.target() &&
      controlsAreControlsOf(prev))
  {
    return prev.cost() - XRotation::gateCost(*context_, prev.numControls());
  }
  return 0;
}


GateSequence PhaseGate::simplification(const YRotation& prev) const
{
  return makeGateSequence(clone(), make_unique<XRotation>(*context_, prev.target(), prev.controls(), -prev.angle()));
//  return {clone(), make_unique<XRotation>(*context_, prev.target(), prev.controls(), -prev.angle())};
}


bool PhaseGate::canRSwap(const YRotation& prev) const
{
  // Gates can be swapped, with the YRotation becoming a XRotation of opposite angle, provided targets match and all
  // controls of the PhaseGate are also controls of the rotation. We must also check that the costs of XRotations and
  // YRotations match.
  return target() == prev.target() && controlsAreControlsOf(prev) &&
         XRotation::gateCost(*context_, prev.numControls()) == prev.cost();;
}


unique_ptr<Gate> PhaseGate::rightRMoverChange(const YRotation& prev) const
{
  // YRotation turns into a XRotation of opposite angle.
  return make_unique<XRotation>(*context_, prev.target(), prev.controls(), -prev.angle());
}


unique_ptr<Gate> PhaseGate::leftRMoverChange(const YRotation& prev) const
{
  return {};  // PhaseGate is unchanged.
}


long PhaseGate::canSimplify(const ZRotation& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  //
  // If both gates have no controls then they simplify to a single ZRotation, up to an overall phase which we can
  // ignore. If both have controls and the controls are the same, then the phase difference becomes important because it
  // only applies when the control bits are set. Hence, in this case, the gates simplify to a ZRotation (with the same
  // controls) and a PiByEight (or equivalent) on the control bits.
  if (allQbitsMatch(prev))
  {
    if (numControls() == 0)
    {
      return cost();
    }
    if (piByEightEquivalentAvailable(*context_, prev.controls()))
    {
      return cost() - piByEightEquivalentCost(*context_, prev.controls());
    }
  }
  return 0;
}


GateSequence PhaseGate::simplification(const ZRotation& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  //
  // If both gates have no controls then they simplify to a single ZRotation, up to an overall phase which we can
  // ignore. If both have controls and the controls are the same, then the phase difference becomes important because
  // it only applies when the control bits are set. Hence, in this case, the gates simplify to a ZRotation (with the
  // same controls) and a PiByEight (or equivalent) on the control bits.
  auto firstGate = make_unique<ZRotation>(*context_, prev.target(), prev.controls(), prev.angle() - pi / 2);
  if (numControls() == 0)
  {
    return makeGateSequence(std::move(firstGate));
//    return {std::move(firstGate)};
  }
  return makeGateSequence(std::move(firstGate), piByEightEquivalent(*context_, prev.controls()));
//  return {std::move(firstGate), piByEightEquivalent(*context_, prev.controls())};
}


long PhaseGate::canSimplify(const ArbitraryPhase& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  if (target() == prev.target() && controls() == prev.controls())
  {
    return cost();
  }
  return 0;
}


GateSequence PhaseGate::simplification(const ArbitraryPhase& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  return makeGateSequence(make_unique<ArbitraryPhase>(*context_, target(), controls(), prev.angle() + pi / 2));
//  return {make_unique<ArbitraryPhase>(*context_, target(), controls(), prev.angle() + pi / 2)};
}


vector<vector<cmplx> > PhaseGate::matrix_() const
{
  using constants::i;
  return {{1, 0}, {0, i}};
}

//----------------------------------------------------------------------------------------------------------------------

PhaseInv::PhaseInv(const CircuitContext& context) :
SingleTargetGate(context),
//DiagonalGate(context),
PhaseTypeGate(context)
{
}


PhaseInv::PhaseInv(const CircuitContext& context, int target, const Controls& controls) :
SingleTargetGate(context, target, controls),
//DiagonalGate(context, target, controls),
PhaseTypeGate(context, target, controls)
{
}


PhaseInv::PhaseInv(const CircuitContext& context, const Controls& controls) :  // The target is any one of the controls.
SingleTargetGate(context),
//DiagonalGate(context),
PhaseTypeGate(context, controls)
{
}


int PhaseInv::gateTypeId() const
{
  return gate_type_id;
}


unique_ptr<Gate> PhaseInv::clone() const
{
  return make_unique<PhaseInv>(*this);  // Will get converted to unique_ptr<Gate>.
}


PhaseInv* PhaseInv::rawClone() const
{
  return new PhaseInv(*this);
}


std::string PhaseInv::name() const
{
  return "PhaseInv";
}


bool PhaseInv::gateAvailable(const CircuitContext& context, int target, const Controls& controls)  // Static
{
  // Static function for determining the availability of a Gate without needing to create one. (The corresponding
  // non-static function is useful when we have a Gate but don't necessarily know the precise gate type.)
  // (We don't need to provide the target bit at present - just the number of controls. However, in future we may have
  // more complicated rules on availability.)
  return context.gateTypeAvailable(gate_type_id, controls.numControls());
}


bool PhaseInv::gateAvailable(const CircuitContext& context, const Controls& controls)  // Static
{
  // Static function for determining the availability of a Gate without needing to create one. (The corresponding
  // non-static function is useful when we have a Gate but don't necessarily know the precise gate type.)
  // Parameter 'controls' is actually the controls AND the target. Hence the number of actual controls for the desired
  // gate is numControls() - 1.
  assert(controls.numControls() > 0);
  return context.gateTypeAvailable(gate_type_id, controls.numControls() - 1);
}


long PhaseInv::gateCost(const CircuitContext& context, int numControls)  // Static
{
  return context.gateCost(gate_type_id, numControls);
}


long PhaseInv::canSimplify(const Gate& next) const
{
  return next.canSimplify(*this);  // Double dispatch.
}


GateSequence PhaseInv::simplification(const Gate& next) const
{
  return next.simplification(*this);  // Double dispatch.
}


bool PhaseInv::canRSwap(const Gate& next) const
{
  return next.canRSwap(*this);  // Double dispatch.
}


unique_ptr<Gate> PhaseInv::rightRMoverChange(const Gate& next) const
{
  return next.rightRMoverChange(*this);  // Double dispatch.
}


unique_ptr<Gate> PhaseInv::leftRMoverChange(const Gate& next) const
{
  return next.leftRMoverChange(*this);  // Double dispatch.
}


long PhaseInv::canSimplify(const PiByEight& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.canSimplify(*this);
}


GateSequence PhaseInv::simplification(const PiByEight& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.simplification(*this);
}


long PhaseInv::canSimplify(const PiByEightInv& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.canSimplify(*this);
}


GateSequence PhaseInv::simplification(const PiByEightInv& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.simplification(*this);
}


long PhaseInv::canSimplify(const PhaseGate& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.canSimplify(*this);
}


GateSequence PhaseInv::simplification(const PhaseGate& prev) const
{
  // Gate order is unimportant - they commute.
  return {};  // Not much code duplicated! :-)
}


long PhaseInv::canSimplify(const PhaseInv& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  if (target() == prev.target() && controls() == prev.controls() &&
      zGateEquivalentAvailable(*context_, target(), controls()))
  {
    return cost() + prev.cost() - zGateEquivalentCost(*context_, target(), controls());
  }
  return 0;
}


GateSequence PhaseInv::simplification(const PhaseInv& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  // Gates becomes a single ZGate or an equivalent ArbitraryPhase.
  return makeGateSequence(zGateEquivalent(*context_, target(), controls()));
//  return {zGateEquivalent(*context_, target(), controls())};
}


long PhaseInv::canSimplify(const ZGate& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  if (target() == prev.target() && controls() == prev.controls() &&
      phaseGateEquivalentAvailable(*context_, target(), controls()))
  {
    return cost() + prev.cost() - phaseGateEquivalentCost(*context_, target(), controls());
  }
  return 0;
}


GateSequence PhaseInv::simplification(const ZGate& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  // Gates become a single PhaseGate or an equivalent ArbitraryPhase.
  return makeGateSequence(phaseGateEquivalent(*context_, target(), controls()));
//  return {phaseGateEquivalent(*context_, target(), controls())};
}


long PhaseInv::canSimplify(const XRotation& prev) const
{
  // If a YRotation is cheaper than an XRotation (which seems unlikely), then we can swap the gates over, turning the
  // XRotation into a YRotation of opposite angle, provided targets match and all controls of the PhaseInv are also
  // controls of the rotation.
  if (YRotation::gateCost(*context_, prev.numControls()) < prev.cost() && target() == prev.target() &&
      controlsAreControlsOf(prev))
  {
    return prev.cost() - YRotation::gateCost(*context_, prev.numControls());
  }
  return 0;
}


GateSequence PhaseInv::simplification(const XRotation& prev) const
{
  return makeGateSequence(clone(), make_unique<YRotation>(*context_, prev.target(), prev.controls(), -prev.angle()));
//  return {clone(), make_unique<YRotation>(*context_, prev.target(), prev.controls(), prev.angle())};
}


bool PhaseInv::canRSwap(const XRotation& prev) const
{
  // Gates can be swapped, with the XRotation becoming a YRotation of opposite angle, provided targets match and all
  // controls of the PhaseInv are also controls of the rotation. We must also check that the cost of XRotations and
  // YRotations match.
  return target() == prev.target() && controlsAreControlsOf(prev) &&
         YRotation::gateCost(*context_, prev.numControls()) == prev.cost();;
}


unique_ptr<Gate> PhaseInv::rightRMoverChange(const XRotation& prev) const
{
  // XRotation turns into a YRotation of opposite angle.
  return make_unique<YRotation>(*context_, prev.target(), prev.controls(), -prev.angle());
}


unique_ptr<Gate> PhaseInv::leftRMoverChange(const XRotation& prev) const
{
  return {};  // PhaseInv is unchanged.
}


long PhaseInv::canSimplify(const YRotation& prev) const
{
  // If an XRotation is cheaper than a YRotation (which seems unlikely), then we can swap the gates over, turning the
  // YRotation into an XRotation of the same angle, provided targets match and all controls of the PhaseInv are also
  // controls of the rotation.
  if (XRotation::gateCost(*context_, prev.numControls()) < prev.cost() && target() == prev.target() &&
      controlsAreControlsOf(prev))
  {
    return prev.cost() - XRotation::gateCost(*context_, prev.numControls());
  }
  return 0;
}


GateSequence PhaseInv::simplification(const YRotation& prev) const
{
  return makeGateSequence(clone(), make_unique<XRotation>(*context_, prev.target(), prev.controls(), prev.angle()));
//  return {clone(), make_unique<XRotation>(*context_, prev.target(), prev.controls(), prev.angle())};
}


bool PhaseInv::canRSwap(const YRotation& prev) const
{
  // Gates can be swapped, with the YRotation becoming a XRotation of the same angle, provided targets match and all
  // controls of the PhaseInv are also controls of the rotation. We must also check that the costs of XRotations and
  // YRotations match.
  return target() == prev.target() && controlsAreControlsOf(prev) &&
         XRotation::gateCost(*context_, prev.numControls()) == prev.cost();
}


unique_ptr<Gate> PhaseInv::rightRMoverChange(const YRotation& prev) const
{
  // YRotation turns into a XRotation of the same angle.
  return make_unique<XRotation>(*context_, prev.target(), prev.controls(), prev.angle());
}


unique_ptr<Gate> PhaseInv::leftRMoverChange(const YRotation& prev) const
{
  return {};  // PhaseGate is unchanged.
}


long PhaseInv::canSimplify(const ZRotation& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  //
  // If both gates have no controls then they simplify to a single ZRotation, up to an overall phase which we can
  // ignore. If both have controls and the controls are the same, then the phase difference becomes important because it
  // only applies when the control bits are set. Hence, in this case, the gates simplify to a ZRotation (with the same
  // controls) and a PiByEightInv (or equivalent) on the control bits.
  //
  // Note that while one may freely swap the target qbits with a control qbits on a PiByEightInv, with no effect on its
  // function, this is not the case for the ZRotation.
  if (allQbitsMatch(prev))
  {
    if (numControls() == 0)
    {
      return cost();
    }
    if (piByEightInvEquivalentAvailable(*context_, prev.controls()))
    {
      return cost() - piByEightInvEquivalentCost(*context_, prev.controls());
    }
  }
  return 0;
}


GateSequence PhaseInv::simplification(const ZRotation& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  //
  // If both gates have no controls then they simplify to a single ZRotation, up to an overall phase which we can
  // ignore. If both have controls and the controls are the same, then the phase difference becomes important because it
  // only applies when the control bits are set. Hence, in this case, the gates simplify to a ZRotation (with the same
  // controls) and a PiByEight (or equivalent) on the control bits.
  auto firstGate = make_unique<ZRotation>(*context_, prev.target(), prev.controls(), prev.angle() + pi / 2);
  if (numControls() == 0)
  {
    return makeGateSequence(std::move(firstGate));
//    return {std::move(firstGate)};
  }
  return makeGateSequence(std::move(firstGate), piByEightInvEquivalent(*context_, prev.controls()));
//  return {std::move(firstGate), piByEightInvEquivalent(*context_, prev.controls())};
}


long PhaseInv::canSimplify(const ArbitraryPhase& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  if (target() == prev.target() && controls() == prev.controls())
  {
    return cost();
  }
  return 0;
}


GateSequence PhaseInv::simplification(const ArbitraryPhase& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  return makeGateSequence(make_unique<ArbitraryPhase>(*context_, target(), controls(), prev.angle() - pi / 2));
//  return {make_unique<ArbitraryPhase>(*context_, target(), controls(), prev.angle() - pi / 2)};
}


vector<vector<cmplx> > PhaseInv::matrix_() const
{
  using constants::i;
  return {{1, 0}, {0, -i}};
}

//----------------------------------------------------------------------------------------------------------------------

XGate::XGate(const CircuitContext& context) :
SingleTargetGate(context),
XTypeGate(context)
{
}


XGate::XGate(const CircuitContext& context, int target, const Controls& controls) :
SingleTargetGate(context, target, controls),
XTypeGate(context, target, controls)
{
}


int XGate::gateTypeId() const
{
  return gate_type_id;
}


unique_ptr<Gate> XGate::clone() const
{
  return make_unique<XGate>(*this);  // Will get converted to unique_ptr<Gate>.
}


XGate* XGate::rawClone() const
{
  return new XGate(*this);
}


std::string XGate::name() const
{
  return "X";
}


bool XGate::gateAvailable(const CircuitContext& context, int target, const Controls& controls)  // Static
{
  // Static function for determining the availability of a Gate without needing to create one. (The corresponding
  // non-static function is useful when we have a Gate but don't necessarily know the precise gate type.)
  // (We don't need to provide the target bit at present - just the number of controls. However, in future we may have
  // more complicated rules on availability.)
  return context.gateTypeAvailable(gate_type_id, controls.numControls());
}


long XGate::gateCost(const CircuitContext& context, int numControls)  // Static
{
  return context.gateCost(gate_type_id, numControls);
}


bool XGate::matricesAnticommute(const SingleTargetGate& other) const
{
  return other.matricesAnticommute(*this);  // Double dispatch.
}


bool XGate::matricesAnticommute(const YGate& other) const
{
  return true;
}


bool XGate::matricesAnticommute(const ZGate& other) const
{
  return true;
}


long XGate::canSimplify(const Gate& next) const
{
  return next.canSimplify(*this);  // Double dispatch.
}


GateSequence XGate::simplification(const Gate& next) const
{
  return next.simplification(*this);  // Double dispatch.
}


bool XGate::canHSwap(const Gate& next) const
{
  return next.canHSwap(*this);  // Double dispatch.
}


unique_ptr<Gate> XGate::rightHMoverChange(const Gate& next) const
{
  return next.rightHMoverChange(*this);  // Double dispatch.
}


unique_ptr<Gate> XGate::leftHMoverChange(const Gate& next) const
{
  return next.leftHMoverChange(*this);  // Double dispatch.
}


bool XGate::canRSwap(const Gate& next) const
{
  return next.canRSwap(*this);  // Double dispatch.
}


unique_ptr<Gate> XGate::rightRMoverChange(const Gate& next) const
{
  return next.rightRMoverChange(*this);  // Double dispatch.
}


unique_ptr<Gate> XGate::leftRMoverChange(const Gate& next) const
{
  return next.leftRMoverChange(*this);  // Double dispatch.
}


long XGate::canSimplify(const Hadamard& prev) const
{
  // The conditions under which we can simplify HX and the resulting cost change are the same as those for XH.
  return prev.canSimplify(*this);
}


GateSequence XGate::simplification(const Hadamard& prev) const
{
  // The result of a simplification of HX is NOT the same as for XH. We get ZH not HZ.
  // This function is only called if a simplification exists. Hence we need not check again.
  return makeGateSequence(zGateEquivalent(*context_, target(), controls()), prev.clone());
//  return {zGateEquivalent(*context_, target(), controls()), prev.clone()};
}


bool XGate::canHSwap(const Hadamard& prev) const
{
  // We can do HX->ZH with no change in cost iff we can do XH->HZ with no change in cost. Requirements on qbits and gate
  // availability are the same.
  return prev.canHSwap(*this);
}


unique_ptr<Gate> XGate::rightHMoverChange(const Hadamard& prev) const
{
  return {};  // The Hadamard remains unchanged.
}


unique_ptr<Gate> XGate::leftHMoverChange(const Hadamard& prev) const
{
  // The gate that X changes to in HX->ZH swap is the same (i.e. a ZGate) as in the XH->HZ swap.
  return prev.rightHMoverChange(*this);
}


long XGate::canSimplify(const XGate& prev) const
{
  if (target() == prev.target() && controls() == prev.controls())
  {
    return cost() + prev.cost();
  }
  return 0;
}


GateSequence XGate::simplification(const XGate& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  return {};
}


long XGate::canSimplify(const YGate& prev) const
{
  // If both gates have the same target and no controls, they simplify to a single ZGate (or equivalent), up to an
  // overall phase which we can ignore. If they have controls and target and controls agree, then the phase shift
  // becomes important because it only applies when the control bits are set. Hence, in this case, the gates 'simplify'
  // to a ZGate (with the same controls) and a PhaseGate (because XY = iZ) on the control bits (or equivalents).
  // Naturally, if targets or controls disagree, no simplification takes place.
  if (target() == prev.target() && controls() == prev.controls() &&
      zGateEquivalentAvailable(*context_, target(), controls()))
  {
    auto costImprovement = cost() + prev.cost() - zGateEquivalentCost(*context_, target(), controls());
    if (numControls() == 0)
    {
      return costImprovement;
    }

    // In the presence of control qbits, we need an additional PhaseGate.
    if (phaseGateEquivalentAvailable(*context_, controls()))
    {
      return costImprovement - phaseGateEquivalentCost(*context_, controls());
    }
  }
  return 0;
}


GateSequence XGate::simplification(const YGate& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  // Note that we could place the PhaseGate before the ZGate, but since swaps will be applied to get the circuit into
  // canonical form, the order of these gates is unimportant. (Would placing PhaseGate first be more efficient,
  // eliminating the need for a swap?)
  if (numControls() == 0)
  {
    return makeGateSequence(zGateEquivalent(*context_, target(), controls()));
//    return {zGateEquivalent(*context_, target(), controls())};
  }
  return makeGateSequence(zGateEquivalent(*context_, target(), controls()), phaseGateEquivalent(*context_, controls()));
//  return {zGateEquivalent(*context_, target(), controls()), phaseGateEquivalent(*context_, controls())};
}


long XGate::canSimplify(const ZGate& prev) const
{
  // If both gates have the same target and no controls, they simplify to a single YGate, up to an overall phase which
  // we ignore. Otherwise, if we can arrange for the target and control bits to agree (i.e. if the same set of qbits is
  // involved in each gate), then the phase shift becomes important because it only applies when the control bits are
  // set. Hence, in this case, the gates 'simplify' to a YGate and a PhaseInv (because XZ = -iY) on the control bits
  // (or equivalent). Naturally, if targets or controls disagree, no simplification takes place.
  //
  // Note that while we may swap the target for a control on the ZGate, we cannot do the same on the XGate. Hence the
  // new YGate must have the same target and controls as the old XGate.
  if (allQbitsMatch(prev) && YGate::gateAvailable(*context_, target(), controls()))
  {
    auto costImprovement = cost() + prev.cost() - YGate::gateCost(*context_, numControls());
    if (numControls() == 0)
    {
      return costImprovement;
    }

    // In the presence of control qbits, we need an additional PhaseInv.
    if (phaseInvEquivalentAvailable(*context_, controls()))
    {
      return costImprovement - phaseInvEquivalentCost(*context_, controls());
    }
  }
  return 0;
}


GateSequence XGate::simplification(const ZGate& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  // Note that we could place the PhaseInv before the YGate, but since swaps will be applied to get the circuit into
  // canonical form, the order of these gates is unimportant. (Would placing PhaseInv first be more efficient,
  // eliminating the need for a swap?)
  if (numControls() == 0)
  {
    return makeGateSequence(make_unique<YGate>(*context_, target(), controls()));
//    return {make_unique<YGate>(*context_, target(), controls())};
  }
  return makeGateSequence(make_unique<YGate>(*context_, target(), controls()),
                          phaseInvEquivalent(*context_, controls()));
//  return {make_unique<YGate>(*context_, target(), controls()), phaseInvEquivalent(*context_, controls())};
}


long XGate::canSimplify(const XRotation& prev) const
{
  // If both gates have the same target and no controls, they simplify to a single XRotation, up to an overall phase
  // which we ignore. Otherwise, if target and controls agree, the phase shift becomes important because it only applies
  // when the control bits are set. In this case the gates 'simplify' to an XRotation (with the same controls) and
  // either a PhaseGate or a PhaseInv on the control bits (or the equivalent ArbitraryPhase). Naturally, if targets or
  // controls disagree, no simplification takes place.
  //
  // The angle of the new XRotation depends on whether we use a PhaseGate or a PhaseInv. It is possible that one option
  // may lead to better subsequent simplification. For example, if the new gate occurs next to a PiByEight, a PhaseInv
  // would produce a PiByEightInv, while a PhaseGate would produce a more expensive ArbitraryPhase, if available.
  // However, since the calling function only expects one simplification to be returned, we must choose. If one is
  // cheaper than the other, then we choose that one. If both cost the same, we choose the PhaseGate. If only
  // ArbitraryPhase is available, we choose one equivalent to the PhaseGate.
  if (target() == prev.target() && controls() == prev.controls())
  {
    if (numControls() == 0)
    {
      return cost();
    }
    if (phaseGateEquivalentAvailable(*context_, controls()))
    {
      if (phaseInvEquivalentAvailable(*context_, controls()) &&
          phaseInvEquivalentCost(*context_, controls()) < phaseGateEquivalentCost(*context_, controls()))
      {
        return cost() - phaseInvEquivalentCost(*context_, controls());
      }
      return cost() - phaseGateEquivalentCost(*context_, controls());
    }
    if (phaseInvEquivalentAvailable(*context_, controls()))
    {
      return cost() - phaseInvEquivalentCost(*context_, controls());
    }
  }
  return 0;
}


GateSequence XGate::simplification(const XRotation& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  if (numControls() == 0)
  {
    return makeGateSequence(make_unique<XRotation>(*context_, target(), controls(), prev.angle() - pi));
//    return {make_unique<XRotation>(*context_, target(), controls(), prev.angle() - pi)};
  }
  if (phaseGateEquivalentAvailable(*context_, controls()))
  {
    if (phaseInvEquivalentAvailable(*context_, controls()) &&
        phaseInvEquivalentCost(*context_, controls()) < phaseGateEquivalentCost(*context_, controls()))
    {
      return makeGateSequence(make_unique<XRotation>(*context_, target(), controls(), prev.angle() + pi),
                              phaseInvEquivalent(*context_, controls()));
//      return {make_unique<XRotation>(*context_, target(), controls(), prev.angle() + pi),
//              phaseInvEquivalent(*context_, controls())};
    }
    return makeGateSequence(make_unique<XRotation>(*context_, target(), controls(), prev.angle() - pi),
                            phaseGateEquivalent(*context_, controls()));
//    return {make_unique<XRotation>(*context_, target(), controls(), prev.angle() - pi),
//            phaseGateEquivalent(*context_, controls())};
  }

  // PhaseGate (or equivalent) unavailable, which leaves PhaseInv. (Could just create a PhaseInv rather than
  // phaseInvEquivalent.)
  return makeGateSequence(make_unique<XRotation>(*context_, target(), controls(), prev.angle() + pi),
                          phaseInvEquivalent(*context_, controls()));
//  return {make_unique<XRotation>(*context_, target(), controls(), prev.angle() + pi),
//          phaseInvEquivalent(*context_, controls())};
}


bool XGate::canRSwap(const YRotation& prev) const
{
  // Gates can be swapped, with the rotation angle being inverted, if targets match and all controls of the XGate are
  // also controls of the rotation.
  return target() == prev.target() && controlsAreControlsOf(prev);
}


unique_ptr<Gate> XGate::rightRMoverChange(const YRotation& prev) const
{
  // YRotation's angle gets inverted.
  return make_unique<YRotation>(*context_, prev.target(), prev.controls(), -prev.angle());
}


unique_ptr<Gate> XGate::leftRMoverChange(const YRotation& prev) const
{
  return {};  // XGate is unchanged.
}


bool XGate::canRSwap(const ZRotation& prev) const
{
  // Gates can be swapped, with the rotation angle being inverted, if targets match and all controls of the XGate are
  // also controls of the rotation.
  return target() == prev.target() && controlsAreControlsOf(prev);
}


unique_ptr<Gate> XGate::rightRMoverChange(const ZRotation& prev) const
{
  // ZRotation's angle gets inverted.
  return make_unique<ZRotation>(*context_, prev.target(), prev.controls(), -prev.angle());
}


unique_ptr<Gate> XGate::leftRMoverChange(const ZRotation& prev) const
{
  return {};  // XGate is unchanged.
}


bool XGate::canRSwap(const ArbitraryPhase& prev) const
{
  // Gates can be swapped, with the rotation angle being inverted, if targets match and there are no controls.
  return target() == prev.target() && numControls() == 0 && prev.numControls() == 0;
}


unique_ptr<Gate> XGate::rightRMoverChange(const ArbitraryPhase& prev) const
{
  // ArbitraryPhase's angle gets inverted.
  return make_unique<ArbitraryPhase>(*context_, prev.target(), prev.controls(), -prev.angle());
}


unique_ptr<Gate> XGate::leftRMoverChange(const ArbitraryPhase& prev) const
{
  return {};  // XGate is unchanged.
}


vector<vector<cmplx> > XGate::matrix_() const
{
  return {{0, 1}, {1, 0}};
}

//----------------------------------------------------------------------------------------------------------------------

YGate::YGate(const CircuitContext& context) :
SingleTargetGate(context),
YTypeGate(context)
{
}


YGate::YGate(const CircuitContext& context, int target, const Controls& controls) :
SingleTargetGate(context, target, controls),
YTypeGate(context, target, controls)
{
}


int YGate::gateTypeId() const
{
  return gate_type_id;
}


unique_ptr<Gate> YGate::clone() const
{
  return make_unique<YGate>(*this);  // Will get converted to unique_ptr<Gate>.
}


YGate* YGate::rawClone() const
{
  return new YGate(*this);
}


std::string YGate::name() const
{
  return "Y";
}


bool YGate::gateAvailable(const CircuitContext& context, int target, const Controls& controls)  // Static
{
  // Static function for determining the availability of a Gate without needing to create one. (The corresponding
  // non-static function is useful when we have a Gate but don't necessarily know the precise gate type.)
  // (We don't need to provide the target bit at present - just the number of controls. However, in future we may have
  // more complicated rules on availability.)
  return context.gateTypeAvailable(gate_type_id, controls.numControls());
}


long YGate::gateCost(const CircuitContext& context, int numControls)  // Static
{
  return context.gateCost(gate_type_id, numControls);
}


bool YGate::matricesAnticommute(const SingleTargetGate& other) const
{
  return other.matricesAnticommute(*this);  // Double dispatch.
}


bool YGate::matricesAnticommute(const XGate& other) const
{
  return true;
}


bool YGate::matricesAnticommute(const ZGate& other) const
{
  return true;
}


long YGate::canSimplify(const Gate& next) const
{
  return next.canSimplify(*this);  // Double dispatch.
}


GateSequence YGate::simplification(const Gate& next) const
{
  return next.simplification(*this);  // Double dispatch.
}


bool YGate::canHSwap(const Gate& next) const
{
  return next.canHSwap(*this);  // Double dispatch.
}


unique_ptr<Gate> YGate::rightHMoverChange(const Gate& next) const
{
  return next.rightHMoverChange(*this);  // Double dispatch.
}


unique_ptr<Gate> YGate::leftHMoverChange(const Gate& next) const
{
  return next.leftHMoverChange(*this);  // Double dispatch.
}


bool YGate::canRSwap(const Gate& next) const
{
  return next.canRSwap(*this);  // Double dispatch.
}


unique_ptr<Gate> YGate::rightRMoverChange(const Gate& next) const
{
  return next.rightRMoverChange(*this);  // Double dispatch.
}


unique_ptr<Gate> YGate::leftRMoverChange(const Gate& next) const
{
  return next.leftRMoverChange(*this);  // Double dispatch.
}


bool YGate::canHSwap(const Hadamard& prev) const
{
  // We can do HY->YH with no change in cost iff we can do YH->HY with no change in cost. Requirements on qbits and gate
  // availability are the same.
  return prev.canHSwap(*this);
}


unique_ptr<Gate> YGate::rightHMoverChange(const Hadamard& prev) const
{
  return {};  // The Hadamard remains unchanged.
}


unique_ptr<Gate> YGate::leftHMoverChange(const Hadamard& prev) const
{
  // The gate that Y changes to in HY->YH swap is the same (i.e. just a YGate) as in the YH->HY swap.
  return prev.rightHMoverChange(*this);
}


long YGate::canSimplify(const XGate& prev) const
{
  // If both gates have the same target and no controls, they simplify to a single ZGate (or equivalent), up to an
  // overall phase which we can ignore. If they have controls and target and controls agree, then the phase shift
  // becomes important because it only applies when the control bits are set. Hence, in this case, the gates 'simplify'
  // to a ZGate (with the same controls) and a PhaseInv (because YX = -iZ) on the control bits (or equivalents).
  // Naturally, if targets or controls disagree, no simplification takes place.
  if (target() == prev.target() && controls() == prev.controls() &&
      zGateEquivalentAvailable(*context_, target(), controls()))
  {
    auto costImprovement = cost() + prev.cost() - zGateEquivalentCost(*context_, target(), controls());
    if (numControls() == 0)
    {
      return costImprovement;
    }

    // In the presence of control qbits, we need an additional PhaseInv, or equivalent.
    if (phaseInvEquivalentAvailable(*context_, controls()))
    {
      return costImprovement - phaseInvEquivalentCost(*context_, controls());
    }
  }
  return 0;
}


GateSequence YGate::simplification(const XGate& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  // Note that we could place the PhaseInv before the ZGate, but since swaps will be applied to get the circuit into
  // canonical form, the order of these gates is unimportant. (Would placing PhaseInv first be more efficient,
  // eliminating the need for a swap?)
  if (numControls() == 0)
  {
    return makeGateSequence(zGateEquivalent(*context_, target(), controls()));
//    return {zGateEquivalent(*context_, target(), controls())};
  }
  return makeGateSequence(zGateEquivalent(*context_, target(), controls()), phaseInvEquivalent(*context_, controls()));
//  return {zGateEquivalent(*context_, target(), controls()), phaseInvEquivalent(*context_, controls())};
}


long YGate::canSimplify(const YGate& prev) const
{
  if (target() == prev.target() && controls() == prev.controls())
  {
    return cost() + prev.cost();
  }
  return 0;
}


GateSequence YGate::simplification(const YGate& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  return {};
}


long YGate::canSimplify(const ZGate& prev) const
{
  // If both gates have the same target and no controls, they simplify to a single XGate, up to an overall phase which
  // we ignore. Otherwise, if we can arrange for the target and control bits to agree (i.e. if the same set of qbits is
  // involved in each gate), then the phase shift becomes important because it only applies when the control bits are
  // set. Hence, in this case, the gates 'simplify' to an XGate and a PhaseGate (because YZ = iX) on the control bits.
  // Naturally, if targets or controls disagree, no simplification takes place.
  //
  // Note that while we may swap the target for a control on the ZGate, we cannot do the same on the YGate. Hence the
  // new XGate must have the same target and controls as the old YGate.
  if (allQbitsMatch(prev) && XGate::gateAvailable(*context_, target(), controls()))
  {
    auto costImprovement = cost() + prev.cost() - XGate::gateCost(*context_, numControls());
    if (numControls() == 0)
    {
      return costImprovement;
    }

    // In the presence of control qbits, we need an additional PhaseGate.
    if (phaseGateEquivalentAvailable(*context_, controls()))
    {
      return costImprovement - phaseGateEquivalentCost(*context_, controls());
    }
  }
  return 0;
}


GateSequence YGate::simplification(const ZGate& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  // Note that we could place the PhaseGate before the XGate, but since swaps will be applied to get the circuit into
  // canonical form, the order of these gates is unimportant. (Would placing PhaseGate first be more efficient,
  // eliminating the need for a swap?)
  if (numControls() == 0)
  {
    return makeGateSequence(make_unique<XGate>(*context_, target(), controls()));
//    return {make_unique<XGate>(*context_, target(), controls())};
  }
  return makeGateSequence(make_unique<XGate>(*context_, target(), controls()),
                          phaseGateEquivalent(*context_, controls()));
//  return {make_unique<XGate>(*context_, target(), controls()), phaseGateEquivalent(*context_, controls())};
}


bool YGate::canRSwap(const XRotation& prev) const
{
  // Gates can be swapped, with the rotation angle being inverted, if targets match and all controls of the YGate are
  // also controls of the rotation.
  return target() == prev.target() && controlsAreControlsOf(prev);
}


unique_ptr<Gate> YGate::rightRMoverChange(const XRotation& prev) const
{
  // XRotation's angle gets inverted.
  return make_unique<XRotation>(*context_, prev.target(), prev.controls(), -prev.angle());
}


unique_ptr<Gate> YGate::leftRMoverChange(const XRotation& prev) const
{
  return {};  // YGate is unchanged.
}


long YGate::canSimplify(const YRotation& prev) const  // (Can we avoid code duplication with a template function??)
{
  // If both gates have the same target and no controls, they simplify to a single YRotation, up to an overall phase
  // which we ignore. Otherwise, if target and controls agree, the phase shift becomes important because it only applies
  // when the control bits are set. In this case the gates 'simplify' to an YRotation (with the same controls) and
  // either a PhaseGate or a PhaseInv on the control bits (or the equivalent ArbitraryPhase). Naturally, if targets or
  // controls disagree, no simplification takes place.
  //
  // The angle of the new YRotation depends on whether we use a PhaseGate or a PhaseInv. It is possible that one option
  // may lead to better subsequent simplification. For example, if the new gate occurs next to a PiByEight, a PhaseInv
  // would produce a PiByEightInv, while a PhaseGate would produce a more expensive ArbitraryPhase, if available.
  // However, since the calling function only expects one simplification to be returned, we must choose. If one is
  // cheaper than the other, then we choose that one. If both cost the same, we choose the PhaseGate. If only
  // ArbitraryPhase is available, we choose one equivalent to the PhaseGate.
  if (target() == prev.target() && controls() == prev.controls())
  {
    if (numControls() == 0)
    {
      return cost();
    }
    if (phaseGateEquivalentAvailable(*context_, controls()))
    {
      if (phaseInvEquivalentAvailable(*context_, controls()) &&
          phaseInvEquivalentCost(*context_, controls()) < phaseGateEquivalentCost(*context_, controls()))
      {
        return cost() - phaseInvEquivalentCost(*context_, controls());
      }
      return cost() - phaseGateEquivalentCost(*context_, controls());
    }
    if (phaseInvEquivalentAvailable(*context_, controls()))
    {
      return cost() - phaseInvEquivalentCost(*context_, controls());
    }
  }
  return 0;
}


GateSequence YGate::simplification(const YRotation& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  if (numControls() == 0)
  {
    return makeGateSequence(make_unique<YRotation>(*context_, target(), controls(), prev.angle() - pi));
//    return {make_unique<YRotation>(*context_, target(), controls(), prev.angle() - pi)};
  }
  if (phaseGateEquivalentAvailable(*context_, controls()))
  {
    if (phaseInvEquivalentAvailable(*context_, controls()) &&
        phaseInvEquivalentCost(*context_, controls()) < phaseGateEquivalentCost(*context_, controls()))
    {
      return makeGateSequence(make_unique<YRotation>(*context_, target(), controls(), prev.angle() + pi),
                              phaseInvEquivalent(*context_, controls()));
//      return {make_unique<YRotation>(*context_, target(), controls(), prev.angle() + pi),
//              phaseInvEquivalent(*context_, controls())};
    }
    return makeGateSequence(make_unique<YRotation>(*context_, target(), controls(), prev.angle() - pi),
                            phaseGateEquivalent(*context_, controls()));
//    return {make_unique<YRotation>(*context_, target(), controls(), prev.angle() - pi),
//            phaseGateEquivalent(*context_, controls())};
  }

  // PhaseGate (or equivalent) unavailable, which leaves PhaseInv. (Could just create a PhaseInv rather than
  // phaseInvEquivalent.)
  return makeGateSequence(make_unique<YRotation>(*context_, target(), controls(), prev.angle() + pi),
                          phaseInvEquivalent(*context_, controls()));
//  return {make_unique<YRotation>(*context_, target(), controls(), prev.angle() + pi),
//          phaseInvEquivalent(*context_, controls())};
}


bool YGate::canRSwap(const ZRotation& prev) const
{
  // Gates can be swapped, with the rotation angle being inverted, if targets match and all controls of the YGate are
  // also controls of the rotation.
  return target() == prev.target() && controlsAreControlsOf(prev);
}


unique_ptr<Gate> YGate::rightRMoverChange(const ZRotation& prev) const
{
  // ZRotation's angle gets inverted.
  return make_unique<ZRotation>(*context_, prev.target(), prev.controls(), -prev.angle());
}


unique_ptr<Gate> YGate::leftRMoverChange(const ZRotation& prev) const
{
  return {};  // YGate is unchanged.
}


bool YGate::canRSwap(const ArbitraryPhase& prev) const
{
  // Gates can be swapped, with the rotation angle being inverted, if targets match and neither gate has controls.
  return target() == prev.target() && numControls() == 0 && prev.numControls() == 0;
}


unique_ptr<Gate> YGate::rightRMoverChange(const ArbitraryPhase& prev) const
{
  // ArbitraryPhase's angle gets inverted.
  return make_unique<ArbitraryPhase>(*context_, prev.target(), prev.controls(), -prev.angle());
}


unique_ptr<Gate> YGate::leftRMoverChange(const ArbitraryPhase& prev) const
{
  return {};  // YGate is unchanged.
}


vector<vector<cmplx> > YGate::matrix_() const
{
  using constants::i;
  return {{0, -i}, {i, 0}};
}

//----------------------------------------------------------------------------------------------------------------------

ZGate::ZGate(const CircuitContext& context) :
SingleTargetGate(context),
//DiagonalGate(context),
PhaseTypeGate(context)
{
}


ZGate::ZGate(const CircuitContext& context, int target, const Controls& controls) :
SingleTargetGate(context, target, controls),
//DiagonalGate(context, target, controls),
PhaseTypeGate(context, target, controls)
{
}


ZGate::ZGate(const CircuitContext& context, const Controls& controls) :
SingleTargetGate(context),
//DiagonalGate(context, target, controls),
PhaseTypeGate(context, controls)
{
}


int ZGate::gateTypeId() const
{
  return gate_type_id;
}


unique_ptr<Gate> ZGate::clone() const
{
  return make_unique<ZGate>(*this);  // Will get converted to unique_ptr<Gate>.
}


ZGate* ZGate::rawClone() const
{
  return new ZGate(*this);
}


std::string ZGate::name() const
{
  return "Z";
}


bool ZGate::gateAvailable(const CircuitContext& context, int target, const Controls& controls)  // Static
{
  // Static function for determining the availability of a Gate without needing to create one. (The corresponding
  // non-static function is useful when we have a Gate but don't necessarily know the precise gate type.)
  // (We don't need to provide the target bit at present - just the number of controls. However, in future we may have
  // more complicated rules on availability.)
  return context.gateTypeAvailable(gate_type_id, controls.numControls());
}


bool ZGate::gateAvailable(const CircuitContext& context, const Controls& controls)  // Static
{
  // Static function for determining the availability of a Gate without needing to create one. (The corresponding
  // non-static function is useful when we have a Gate but don't necessarily know the precise gate type.)
  // Parameter 'controls' is actually the controls AND the target. Hence the number of actual controls for the desired
  // gate is numControls() - 1.
  assert(controls.numControls() > 0);
  return context.gateTypeAvailable(gate_type_id, controls.numControls() - 1);
}


long ZGate::gateCost(const CircuitContext& context, int numControls)  // Static
{
  return context.gateCost(gate_type_id, numControls);
}


bool ZGate::matricesAnticommute(const SingleTargetGate& other) const
{
  return other.matricesAnticommute(*this);  // Double dispatch.
}


bool ZGate::matricesAnticommute(const XGate& other) const
{
  return true;
}


bool ZGate::matricesAnticommute(const YGate& other) const
{
  return true;
}


long ZGate::canSimplify(const Gate& next) const
{
  return next.canSimplify(*this);  // Double dispatch.
}


GateSequence ZGate::simplification(const Gate& next) const
{
  return next.simplification(*this);  // Double dispatch.
}


bool ZGate::canHSwap(const Gate& next) const
{
  return next.canHSwap(*this);  // Double dispatch.
}


unique_ptr<Gate> ZGate::rightHMoverChange(const Gate& next) const
{
  return next.rightHMoverChange(*this);  // Double dispatch.
}


unique_ptr<Gate> ZGate::leftHMoverChange(const Gate& next) const
{
  return next.leftHMoverChange(*this);  // Double dispatch.
}


bool ZGate::canRSwap(const Gate& next) const
{
  return next.canRSwap(*this);  // Double dispatch.
}


unique_ptr<Gate> ZGate::rightRMoverChange(const Gate& next) const
{
  return next.rightRMoverChange(*this);  // Double dispatch.
}


unique_ptr<Gate> ZGate::leftRMoverChange(const Gate& next) const
{
  return next.leftRMoverChange(*this);  // Double dispatch.
}


long ZGate::canSimplify(const Hadamard& prev) const
{
  // The conditions under which we can simplify HZ and the resulting cost change are the same as those for ZH.
  return prev.canSimplify(*this);
}


GateSequence ZGate::simplification(const Hadamard& prev) const
{
  // The result of a simplification of HZ is NOT the same as for ZH. We get XH not HX.
  // This function is only called if a simplification exists. Hence we need not check again.
  Controls xControls(controls());  // New gate's controls will be the set of all bits involved in the old Z, minus...
  xControls.add(target());         // ...the new gate's target bit, which happens to be the target of the Hadamard.
  xControls.remove(prev.target());

  return makeGateSequence(make_unique<XGate>(*context_, prev.target(), xControls), prev.clone());
//  return {make_unique<XGate>(*context_, prev.target(), xControls), prev.clone()};
}


bool ZGate::canHSwap(const Hadamard& prev) const
{
  // Can do HZ->XH whenever we can do ZH->HX, and vice versa.
  return prev.canHSwap(*this);
}


unique_ptr<Gate> ZGate::rightHMoverChange(const Hadamard& prev) const
{
  return {};  // The Hadamard remains unchanged.
}


unique_ptr<Gate> ZGate::leftHMoverChange(const Hadamard& prev) const
{
  // The gate that Z changes to in HZ->XH swap is the same (i.e. an XGate) as in the ZH->HX swap.
  return prev.rightHMoverChange(*this);
}


long ZGate::canSimplify(const PiByEight& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.canSimplify(*this);
}


GateSequence ZGate::simplification(const PiByEight& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.simplification(*this);
}


long ZGate::canSimplify(const PiByEightInv& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.canSimplify(*this);
}


GateSequence ZGate::simplification(const PiByEightInv& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.simplification(*this);
}


long ZGate::canSimplify(const PhaseGate& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.canSimplify(*this);
}


GateSequence ZGate::simplification(const PhaseGate& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.simplification(*this);
}


long ZGate::canSimplify(const PhaseInv& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.canSimplify(*this);
}


GateSequence ZGate::simplification(const PhaseInv& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.simplification(*this);
}


long ZGate::canSimplify(const XGate& prev) const
{
  // If both gates have the same target and no controls, they simplify to a single YGate, up to an overall phase which
  // we ignore. Otherwise, if we can arrange for the target and control bits to agree (i.e. if the same set of qbits is
  // involved in each gate), then the phase shift becomes important because it only applies when the control bits are
  // set. Hence, in this case, the gates 'simplify' to a YGate and a PhaseGate (because ZX = iY) on the control bits (or
  // equivalent). Naturally, if targets or controls disagree, no simplification takes place.
  //
  // Note that while we may swap the target for a control on the ZGate, we cannot do the same on the XGate. Hence the
  // new YGate must have the same target and controls as the old XGate.
  if (allQbitsMatch(prev) && YGate::gateAvailable(*context_, prev.target(), prev.controls()))
  {
    auto costImprovement = cost() + prev.cost() - YGate::gateCost(*context_, prev.numControls());
    if (numControls() == 0)
    {
      return costImprovement;
    }

    // In the presence of control qbits, we need an additional PhaseGate.
   if (phaseGateEquivalentAvailable(*context_, prev.controls()))
    {
      return costImprovement - phaseGateEquivalentCost(*context_, prev.controls());
    }
  }
  return 0;
}


GateSequence ZGate::simplification(const XGate& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  // Note that we could place the PhaseGate before the YGate, but since swaps will be applied to get the circuit into
  // canonical form, the order of these gates is unimportant. (Would placing PhaseGate first be more efficient,
  // eliminating the need for a swap?)
  if (numControls() == 0)
  {
    return makeGateSequence(make_unique<YGate>(*context_, prev.target(), prev.controls()));
//    return {make_unique<YGate>(*context_, prev.target(), prev.controls())};
  }
  return makeGateSequence(make_unique<YGate>(*context_, prev.target(), prev.controls()),
                          phaseGateEquivalent(*context_, prev.controls()));
//  return {make_unique<YGate>(*context_, prev.target(), prev.controls()),
//          phaseGateEquivalent(*context_, prev.controls())};
}


long ZGate::canSimplify(const YGate& prev) const
{
  // If both gates have the same target and no controls, they simplify to a single XGate, up to an overall phase which
  // we ignore. Otherwise, if we can arrange for the target and control bits to agree (i.e. if the same set of qbits is
  // involved in each gate), then the phase shift becomes important because it only applies when the control bits are
  // set. Hence, in this case, the gates 'simplify' to an XGate and a PhaseInv (because ZY = -iX) on the control bits
  // (or equivalent). Naturally, if targets or controls disagree, no simplification takes place.
  //
  // Note that while we may swap the target for a control on the ZGate, we cannot do the same on the YGate. Hence the
  // new XGate must have the same target and controls as the old YGate.
  if (allQbitsMatch(prev) && XGate::gateAvailable(*context_, prev.target(), prev.controls()))
  {
    auto costImprovement = cost() + prev.cost() - XGate::gateCost(*context_, prev.numControls());
    if (numControls() == 0)
    {
      return costImprovement;
    }

    // In the presence of control qbits, we need an additional PhaseInv.
    if (phaseInvEquivalentAvailable(*context_, prev.controls()))
    {
      return costImprovement - phaseInvEquivalentCost(*context_, prev.controls());
    }
  }
  return 0;
}


GateSequence ZGate::simplification(const YGate& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  // Note that we could place the PhaseInv before the XGate, but since swaps will be applied to get the circuit into
  // canonical form, the order of these gates is unimportant. (Would placing PhaseInv first be more efficient,
  // eliminating the need for a swap?)
  if (numControls() == 0)
  {
    return makeGateSequence(make_unique<XGate>(*context_, prev.target(), prev.controls()));
//    return {make_unique<XGate>(*context_, prev.target(), prev.controls())};
  }
  return makeGateSequence(make_unique<XGate>(*context_, prev.target(), prev.controls()),
                          phaseInvEquivalent(*context_, prev.controls()));
//  return {make_unique<XGate>(*context_, prev.target(), prev.controls()),
//          phaseInvEquivalent(*context_, prev.controls())};
}


long ZGate::canSimplify(const ZGate& prev) const
{
  if (target() == prev.target() && controls() == prev.controls())
  {
    return cost() + prev.cost();
  }
  return 0;
}


GateSequence ZGate::simplification(const ZGate& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  return {};
}


long ZGate::canSimplify(const XRotation& prev) const
{
  // If the ZGate operates solely on the control bits of the rotation gate, this simplifies to a single rotation with
  // matrix equal to -1 times the original. For XRotations, this means adding 2pi to the angle.
  if (allQbitsMatchControlsOf(prev))
  {
    return cost();
  }
  return 0;
}


GateSequence ZGate::simplification(const XRotation& prev) const
{
  // Only called if a simplification exists.
  return makeGateSequence(make_unique<XRotation>(*context_, prev.target(), prev.controls(), prev.angle() + 2 * pi));
//  return {make_unique<XRotation>(*context_, prev.target(), prev.controls(), prev.angle() + 2 * pi)};
}


bool ZGate::canRSwap(const XRotation& prev) const
{
  // Gates can be swapped, with the rotation angle being inverted, if we can arrange for the targets to match and all
  // controls of the ZGate to also be controls of the rotation. (Recall that the target of a ZGate can be swapped with a
  // control, without changing the gate.)
  return allAreInvolvedIn(prev) && is_involved(prev.target());
}


unique_ptr<Gate> ZGate::rightRMoverChange(const XRotation& prev) const
{
  // XRotation's angle gets inverted.
  return make_unique<XRotation>(*context_, prev.target(), prev.controls(), -prev.angle());
}


unique_ptr<Gate> ZGate::leftRMoverChange(const XRotation& prev) const
{
  return {};  // ZGate is unchanged.
}


long ZGate::canSimplify(const YRotation& prev) const
{
  // If the ZGate operates solely on the control bits of the rotation gate, this simplifies to a single rotation with
  // matrix equal to -1 times the original. For YRotations, this means adding 2pi to the angle.
  if (allQbitsMatchControlsOf(prev))
  {
    return cost();
  }
  return 0;
}


GateSequence ZGate::simplification(const YRotation& prev) const
{
  // Only called if a simplification exists.
  return makeGateSequence(make_unique<YRotation>(*context_, prev.target(), prev.controls(), prev.angle() + 2 * pi));
//  return {make_unique<YRotation>(*context_, prev.target(), prev.controls(), prev.angle() + 2 * pi)};
}


bool ZGate::canRSwap(const YRotation& prev) const
{
  // Gates can be swapped, with the rotation angle being inverted, if we can arrange for the targets to match and all
  // controls of the ZGate to also be controls of the rotation. (Recall that the target of a ZGate can be swapped with a
  // control, without changing the gate.)
  return allAreInvolvedIn(prev) && is_involved(prev.target());
}


unique_ptr<Gate> ZGate::rightRMoverChange(const YRotation& prev) const
{
  // YRotation's angle gets inverted.
  return make_unique<YRotation>(*context_, prev.target(), prev.controls(), -prev.angle());
}


unique_ptr<Gate> ZGate::leftRMoverChange(const YRotation& prev) const
{
  return {};  // ZGate is unchanged.
}


long ZGate::canSimplify(const ZRotation& prev) const
{
  // There are two cases where a simplification may occur:
  // 1. Both gates work with the same set of bits. This is similar to the situation that occurs with an XGate and
  //    XRotation, or with a YGate and YRotation, with minor differences due to the fact that the target bit of a ZGate
  //    may be swapped for a control bit without changing the gate. With no controls, the gates simplify to a ZRotation,
  //    with angle increased (or decreased) by pi. If the gates have control bits, we must also add a PhaseGate (or
  //    PhaseInv) or equivalent to the control bits.
  // 2. The ZGate works on the control bits of the ZRotation. In this case, the old rotation gate is multiplied by -1,
  //    which is achieved by adding 2pi to the rotation angle.

  // Case 1.
  if (allQbitsMatch(prev))
  {
    if (numControls() == 0)
    {
      return cost();
    }
    if (phaseGateEquivalentAvailable(*context_, prev.controls()))
    {
      if (phaseInvEquivalentAvailable(*context_, prev.controls()) &&
          phaseInvEquivalentCost(*context_, prev.controls()) < phaseGateEquivalentCost(*context_, prev.controls()))
      {
        return cost() - phaseInvEquivalentCost(*context_, prev.controls());
      }
      return cost() - phaseGateEquivalentCost(*context_, prev.controls());
    }
    if (phaseInvEquivalentAvailable(*context_, prev.controls()))
    {
      return cost() - phaseInvEquivalentCost(*context_, prev.controls());
    }
  }

  // Case 2.
  if (allQbitsMatchControlsOf(prev))
  {
    return cost();
  }

  return 0;
}


GateSequence ZGate::simplification(const ZRotation& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.

  // Case 1.
  if (allQbitsMatch(prev))
  {
    if (numControls() == 0)
    {
      return makeGateSequence(make_unique<ZRotation>(*context_, prev.target(), prev.controls(), prev.angle() - pi));
//      return {make_unique<ZRotation>(*context_, prev.target(), prev.controls(), prev.angle() - pi)};
    }
    if (phaseGateEquivalentAvailable(*context_, prev.controls()))
    {
      if (phaseInvEquivalentAvailable(*context_, prev.controls()) &&
          phaseInvEquivalentCost(*context_, prev.controls()) < phaseGateEquivalentCost(*context_, prev.controls()))
      {
        return makeGateSequence(make_unique<ZRotation>(*context_, prev.target(), prev.controls(), prev.angle() + pi),
                                phaseInvEquivalent(*context_, prev.controls()));
//        return {make_unique<ZRotation>(*context_, prev.target(), prev.controls(), prev.angle() + pi),
//                phaseInvEquivalent(*context_, prev.controls())};
      }
      return makeGateSequence(make_unique<ZRotation>(*context_, prev.target(), prev.controls(), prev.angle() - pi),
                              phaseGateEquivalent(*context_, prev.controls()));
//      return {make_unique<ZRotation>(*context_, prev.target(), prev.controls(), prev.angle() - pi),
//              phaseGateEquivalent(*context_, prev.controls())};
    }
    // PhaseGate (or equivalent) unavailable, which leaves PhaseInv. (Could just create a PhaseInv rather than
    // phaseInvEquivalent.)
    return makeGateSequence(make_unique<ZRotation>(*context_, prev.target(), prev.controls(), prev.angle() + pi),
                            phaseInvEquivalent(*context_, prev.controls()));
//    return {make_unique<ZRotation>(*context_, prev.target(), prev.controls(), prev.angle() + pi),
//            phaseInvEquivalent(*context_, prev.controls())};
  }

  // Must be case 2.
  return makeGateSequence(make_unique<ZRotation>(*context_, prev.target(), prev.controls(), prev.angle() + 2 * pi));
//  return {make_unique<ZRotation>(*context_, prev.target(), prev.controls(), prev.angle() + 2 * pi)};
}


long ZGate::canSimplify(const ArbitraryPhase& prev) const
{
  // Returns the cost improvement of the best simplification of the two gates.
  if (target() == prev.target() && controls() == prev.controls())
  {
    return cost();
  }
  return 0;
}


GateSequence ZGate::simplification(const ArbitraryPhase& prev) const
{
  // This function is only called if a simplification exists. Hence we need not check again.
  return makeGateSequence(make_unique<ArbitraryPhase>(*context_, target(), controls(), prev.angle() + pi));
//  return {make_unique<ArbitraryPhase>(*context_, target(), controls(), prev.angle() + pi)};
}


vector<vector<cmplx> > ZGate::matrix_() const
{
  return {{1, 0}, {0, -1}};
}


//----------------------------------------------------------------------------------------------------------------------

// Set up XRotation class simulators for just the constant part, the cosine part and the sine part of the rotation
// matrix. Simulator for the constant part must handle the case where a control bit is not set. The other two must
// ignore this case - hence the 'false'.
const SparseGateSimulator XRotation::const_sim(sparse({{0, 0}, {0, 0}}));
const SparseGateSimulator XRotation::cos_sim(sparse({{1, 0}, {0, 1}}), false);
const SparseGateSimulator XRotation::sin_sim(sparse({{0, constants::i}, {constants::i, 0}}), false);


XRotation::XRotation(const CircuitContext& context) :
SingleTargetGate(context),
RotationGate(context),
XTypeGate(context)
{
}


XRotation::XRotation(const CircuitContext& context, int target, const Controls& controls, double angle) :
SingleTargetGate(context, target, controls),
RotationGate(context, target, controls, angle),
XTypeGate(context, target, controls)
{
}


int XRotation::gateTypeId() const
{
  return gate_type_id;
}


unique_ptr<Gate> XRotation::clone() const
{
  return make_unique<XRotation>(*this);  // Will get converted to unique_ptr<Gate>.
}


XRotation* XRotation::rawClone() const
{
  return new XRotation(*this);
}


std::string XRotation::name() const
{
  return "XRotation";
}


bool XRotation::gateAvailable(const CircuitContext& context, int target, const Controls& controls)  // Static
{
  // Static function for determining the availability of a Gate without needing to create one. (The corresponding
  // non-static function is useful when we have a Gate but don't necessarily know the precise gate type.)
  // (We don't need to provide the target bit at present - just the number of controls. However, in future we may have
  // more complicated rules on availability.)
  return context.gateTypeAvailable(gate_type_id, controls.numControls());
}


long XRotation::gateCost(const CircuitContext& context, int numControls)  // Static
{
  return context.gateCost(gate_type_id, numControls);
}


bool XRotation::canReduce() const
{
  double angleTolerance = context_->angleTolerance();

  if (abs(angle()) < angleTolerance)  // Angle close to zero.
  {
    // Gate essentially does nothing, so we can eliminate it.
    return true;
  }

  if (abs(angle()) > 2 * pi - angleTolerance)  // Angle close to +-2pi. (Angle should NEVER be more than 2pi or less...
  {                                            // ...than -2pi.)
    // Gate (which must have controls) is equivalent to a ZGate on the control bits. If this (or an equivalent) is
    // available and cheaper than the XRotation, we take it.
    assert(numControls() > 0);  // The size of angles for uncontrolled gates should not exceed pi.
    return zGateEquivalentAvailable(*context_, controls()) && zGateEquivalentCost(*context_, controls()) < cost();
  }

  if (abs(angle() - pi) < angleTolerance)  // Angle close to pi.
  {
    // If uncontrolled, the gate is in effect just an XGate. If there are controls, then a PhaseGate is also required on
    // the control bits.
    // (The likelihood that the sum of the costs of an XGate and a PhaseGate exceed that of the original XRotation is
    // motivation for the future consideration of 'reduction' functions that consider pairs of gates. Then we might
    // detect that the cost increase is only temporary, eliminated by the cancellation (or simplification) of the
    // introduced PhaseGate with a neighbouring gate.)
    if (numControls() == 0)
    {
      return XGate::gateAvailable(*context_, target(), controls()) &&
             XGate::gateCost(*context_, numControls()) < cost();
    }
    return XGate::gateAvailable(*context_, target(), controls()) &&
           phaseGateEquivalentAvailable(*context_, controls()) &&
           XGate::gateCost(*context_, numControls()) + phaseGateEquivalentCost(*context_, controls()) < cost();
  }

  if (abs(angle() + pi) < angleTolerance)  // Angle close to -pi.
  {
    // If uncontrolled, the gate is in effect just an XGate. If there are controls, then a PhaseInv is also required on
    // the control bits.
    if (numControls() == 0)
    {
      return XGate::gateAvailable(*context_, target(), controls()) &&
             XGate::gateCost(*context_, numControls()) < cost();
    }
    return XGate::gateAvailable(*context_, target(), controls()) &&
           phaseInvEquivalentAvailable(*context_, controls()) &&
           XGate::gateCost(*context_, numControls()) + phaseInvEquivalentCost(*context_, controls()) < cost();
  }

  return false;  // Angle does not take a special value.
}


GateSequence XRotation::reduction() const
{
  // This should only be called if a reduction has been found to exist. That allows for some simplifications.
  double angleTolerance = context_->angleTolerance();

  if (abs(angle()) < angleTolerance)  // Angle close to zero.
  {
    return {};
  }

  if (abs(angle()) > 2 * pi - angleTolerance)  // Angle close to +-2pi. (Angle should NEVER be more than 2pi or less...
  {                                            // ...than -2pi.)
    return makeGateSequence(zGateEquivalent(*context_, controls()));
  }

  if (abs(angle() - pi) < angleTolerance)  // Angle close to pi.
  {
    auto firstGate = make_unique<XGate>(*context_, target(), controls());
    if (numControls() == 0)
    {
      return makeGateSequence(std::move(firstGate));
    }
    return makeGateSequence(std::move(firstGate), phaseGateEquivalent(*context_, controls()));
  }

  assert(abs(angle() + pi) < angleTolerance);  // Angle must be close to -pi.
  {                                            // Added braces mainly just to match the format of the previous 'if's!
    auto firstGate = make_unique<XGate>(*context_, target(), controls());
    if (numControls() == 0)
    {
      return makeGateSequence(std::move(firstGate));
    }
    return makeGateSequence(std::move(firstGate), phaseInvEquivalent(*context_, controls()));
  }
}


long XRotation::canSimplify(const Gate& next) const
{
  return next.canSimplify(*this);  // Double dispatch.
}


GateSequence XRotation::simplification(const Gate& next) const
{
  return next.simplification(*this);  // Double dispatch.
}


bool XRotation::canHSwap(const Gate& next) const
{
  return next.canHSwap(*this);  // Double dispatch.
}


unique_ptr<Gate> XRotation::rightHMoverChange(const Gate& next) const
{
  return next.rightHMoverChange(*this);  // Double dispatch.
}


unique_ptr<Gate> XRotation::leftHMoverChange(const Gate& next) const
{
  return next.leftHMoverChange(*this);  // Double dispatch.
}


bool XRotation::canRSwap(const Gate& next) const
{
  return next.canRSwap(*this);  // Double dispatch.
}


unique_ptr<Gate> XRotation::rightRMoverChange(const Gate& next) const
{
  return next.rightRMoverChange(*this);  // Double dispatch.
}


unique_ptr<Gate> XRotation::leftRMoverChange(const Gate& next) const
{
  return next.leftRMoverChange(*this);  // Double dispatch.
}


long XRotation::canSimplify(const Hadamard& prev) const
{
  // Situations where H-XRot simplify are precisely the same as those were XRot-H simplify.
  return prev.canSimplify(*this);
}


GateSequence XRotation::simplification(const Hadamard& prev) const
{
  // Simplification of H-XRot does not produce the same sequence as simplification of XRot-H.
  return makeGateSequence(make_unique<ZRotation>(*context_, target(), controls(), angle()), prev.clone());
//  return {make_unique<ZRotation>(*context_, target(), controls(), angle()), prev.clone()};
}


bool XRotation::canHSwap(const Hadamard& prev) const
{
  return prev.canHSwap(*this);
}


unique_ptr<Gate> XRotation::rightHMoverChange(const Hadamard& prev) const
{
  return {};  // The Hadamard remains unchanged.
}


unique_ptr<Gate> XRotation::leftHMoverChange(const Hadamard& prev) const
{
  // The gate that X changes to in HX->ZH swap is the same (i.e. a ZRotation) as in the XH->HZ swap.
  return prev.rightHMoverChange(*this);
}


long XRotation::canSimplify(const PhaseGate& prev) const
{
  // Situations where S-XRot simplify are precisely the same as those where XRot-S simplify.
  return prev.canSimplify(*this);
}


GateSequence XRotation::simplification(const PhaseGate& prev) const
{
  // Simplification of S-XRot does not produce the same sequence as simplification of XRot-S. Here the XRotation turns
  // into a YRotation of opposite angle.
  return makeGateSequence(make_unique<YRotation>(*context_, target(), controls(), -angle()), prev.clone());
  // return {make_unique<YRotation>(*context_, target(), controls, -angle()), prev.clone()};
}


bool XRotation::canRSwap(const PhaseGate& prev) const
{
  return prev.canRSwap(*this);
}


unique_ptr<Gate> XRotation::rightRMoverChange(const PhaseGate& prev) const
{
  return {};  // PhaseGate is unchanged.
}


unique_ptr<Gate> XRotation::leftRMoverChange(const PhaseGate& prev) const
{
  // XRotation turns into a YRotation of opposite angle. Note that this is NOT the same as when the XRotation precedes
  // the PhaseGate.
  return make_unique<YRotation>(*context_, prev.target(), prev.controls(), -angle());
}


long XRotation::canSimplify(const PhaseInv& prev) const
{
  // Situations where S-XRot simplify are precisely the same as those where XRot-S simplify.
  return prev.canSimplify(*this);
}


GateSequence XRotation::simplification(const PhaseInv& prev) const
{
  // Simplification of PhaseInv-XRot does not produce the same sequence as simplification of XRot-PhaseInv. Here the
  // XRotation turns into a YRotation of the same angle.
  return makeGateSequence(make_unique<YRotation>(*context_, target(), controls(), angle()), prev.clone());
  // return {make_unique<YRotation>(*context_, target(), controls, angle()), prev.clone()};
}


bool XRotation::canRSwap(const PhaseInv& prev) const
{
  return prev.canRSwap(*this);
}


unique_ptr<Gate> XRotation::rightRMoverChange(const PhaseInv& prev) const
{
  return {};  // PhaseInv is unchanged.
}


unique_ptr<Gate> XRotation::leftRMoverChange(const PhaseInv& prev) const
{
  // XRotation turns into a YRotation of the same angle. Note that this is NOT the same as when the XRotation precedes
  // the PhaseInv.
  return make_unique<YRotation>(*context_, prev.target(), prev.controls(), angle());
}


long XRotation::canSimplify(const XGate& prev) const
{
  // Gate order is unimportant - they commute whenever a simplification exists. To avoid code duplication, call the
  // function with the gates reversed.
  return prev.canSimplify(*this);
}


GateSequence XRotation::simplification(const XGate& prev) const
{
  // Gate order is unimportant - they commute whenever a simplification exists. To avoid code duplication, call the
  // function with the gates reversed.
  return prev.simplification(*this);
}


bool XRotation::canRSwap(const YGate& prev) const
{
  return prev.canRSwap(*this);
}


unique_ptr<Gate> XRotation::rightRMoverChange(const YGate& prev) const
{
  return prev.leftRMoverChange(*this);
}


unique_ptr<Gate> XRotation::leftRMoverChange(const YGate& prev) const
{
  return prev.rightRMoverChange(*this);
}


long XRotation::canSimplify(const ZGate& prev) const
{
  // Gate order is unimportant - they commute whenever a simplification exists. To avoid code duplication, call the
  // function with the gates reversed.
  return prev.canSimplify(*this);
}


GateSequence XRotation::simplification(const ZGate& prev) const
{
  // Gate order is unimportant - they commute whenever a simplification exists. To avoid code duplication, call the
  // function with the gates reversed.
  return prev.simplification(*this);
}


bool XRotation::canRSwap(const ZGate& prev) const
{
  return prev.canRSwap(*this);
}


unique_ptr<Gate> XRotation::rightRMoverChange(const ZGate& prev) const
{
  return prev.leftRMoverChange(*this);
}


unique_ptr<Gate> XRotation::leftRMoverChange(const ZGate& prev) const
{
  return prev.rightRMoverChange(*this);
}


long XRotation::canSimplify(const XRotation& prev) const
{
  // Provided both target and control bits agree, this gives another XRotation - just add the angles modulo 2pi.
  if (target() == prev.target() && controls() == prev.controls())
  {
    return cost();
  }
  return 0;
}


GateSequence XRotation::simplification(const XRotation& prev) const
{
  return makeGateSequence(make_unique<XRotation>(*context_, target(), controls(), angle() + prev.angle()));
//  return {make_unique<XRotation>(*context_, target(), controls(), angle() + prev.angle())};
}


const GateSimulator& XRotation::part_simulator(int partIndex) const
{
  switch (partIndex)
  {
    case 0:
      return const_sim;
    case 1:
      return cos_sim;
    case 2:
      return sin_sim;
    default:
      throw std::invalid_argument("Invalid partIndex passed to XRotation::part_simulator().");
  }
}


double XRotation::thetaToAngle(double theta) const
{
  return 2 * theta;
}


vector<vector<cmplx> > XRotation::matrix_() const
{
  using std::cos, std::sin, constants::i;
  double angle = this->angle_;  // Yes, the 'this' is necessary! (Apparently not on Visual Studio, but that's a...
  return {{cos(angle / 2), i * sin(angle / 2)}, {i * sin(angle / 2), cos(angle / 2)}};  // ... 'non-conformance'.)
}

//----------------------------------------------------------------------------------------------------------------------

// Set up YRotation class simulators for just the constant part, the cosine part and the sine part of the rotation
// matrix. Simulator for the constant part must handle the case where a control bit is not set. The other two must
// ignore this case - hence the 'false'.
const SparseGateSimulator YRotation::const_sim(sparse({{0, 0}, {0, 0}}));
const SparseGateSimulator YRotation::cos_sim(sparse({{1, 0}, {0, 1}}), false);
const SparseGateSimulator YRotation::sin_sim(sparse({{0, 1}, {-1, 0}}), false);


YRotation::YRotation(const CircuitContext& context) :
SingleTargetGate(context),
RotationGate(context),
YTypeGate(context)
{
}


YRotation::YRotation(const CircuitContext& context, int target, const Controls& controls, double angle) :
SingleTargetGate(context, target, controls),
RotationGate(context, target, controls, angle),
YTypeGate(context, target, controls)
{
}


int YRotation::gateTypeId() const
{
  return gate_type_id;
}


unique_ptr<Gate> YRotation::clone() const
{
  return make_unique<YRotation>(*this);  // Will get converted to unique_ptr<Gate>.
}


YRotation* YRotation::rawClone() const
{
  return new YRotation(*this);
}


std::string YRotation::name() const
{
  return "YRotation";
}


bool YRotation::gateAvailable(const CircuitContext& context, int target, const Controls& controls)  // Static
{
  // Static function for determining the availability of a Gate without needing to create one. (The corresponding
  // non-static function is useful when we have a Gate but don't necessarily know the precise gate type.)
  // (We don't need to provide the target bit at present - just the number of controls. However, in future we may have
  // more complicated rules on availability.)
  // NAME: Called 'gateAvailable' just to avoid clashing with 'available' - the non-static version
  return context.gateTypeAvailable(gate_type_id, controls.numControls());
}


long YRotation::gateCost(const CircuitContext& context, int numControls)  // Static
{
  return context.gateCost(gate_type_id, numControls);
}


bool YRotation::canReduce() const
{
  double angleTolerance = context_->angleTolerance();

  if (abs(angle()) < angleTolerance)  // Angle close to zero.
  {
    // Gate essentially does nothing, so we can eliminate it.
    return true;
  }

  if (abs(angle()) > 2 * pi - angleTolerance)  // Angle close to +-2pi. (Angle should NEVER be more than 2pi or less...
  {                                            // ...than -2pi.)
    // Gate (which must have controls) is equivalent to a ZGate on the control bits. If this (or an equivalent) is
    // available and cheaper than the YRotation, we take it.
    assert(numControls() > 0);  // The size of angles for uncontrolled gates should not exceed pi.
    return zGateEquivalentAvailable(*context_, controls()) && zGateEquivalentCost(*context_, controls()) < cost();
  }

  if (abs(angle() - pi) < angleTolerance)  // Angle close to pi.
  {
    // If uncontrolled, the gate is in effect just an YGate. If there are controls, then a PhaseGate is also required on
    // the control bits.
    // (The likelihood that the sum of the costs of an YGate and a PhaseGate exceed that of the original YRotation a
    // motivation for the future consideration of 'reduction' functions that consider pairs of gates. Then we might
    // detect that the cost increase is only temporary, eliminated by the cancellation (or simplification) of the
    // introduced PhaseGate with a neighbouring gate.)
    if (numControls() == 0)
    {
      return YGate::gateAvailable(*context_, target(), controls()) &&
             YGate::gateCost(*context_, numControls()) < cost();
    }
    return YGate::gateAvailable(*context_, target(), controls()) &&
           phaseGateEquivalentAvailable(*context_, controls()) &&
           YGate::gateCost(*context_, numControls()) + phaseGateEquivalentCost(*context_, controls()) < cost();
  }

  if (abs(angle() + pi) < angleTolerance)  // Angle close to -pi.
  {
    // If uncontrolled, the gate is in effect just a YGate. If there are controls, then a PhaseInv is also required on
    // the control bits.
    if (numControls() == 0)
    {
      return YGate::gateAvailable(*context_, target(), controls()) &&
             YGate::gateCost(*context_, numControls()) < cost();
    }
    return YGate::gateAvailable(*context_, target(), controls()) &&
           phaseInvEquivalentAvailable(*context_, controls()) &&
           YGate::gateCost(*context_, numControls()) + phaseInvEquivalentCost(*context_, controls()) < cost();
  }

  if (abs(abs(angle()) - pi / 2) < angleTolerance)  // Angle close to pi/2 or -pi/2. (Different reductions, but with...
  {                                                 // ...the same conditions.)
    // Gate is equivalent to an XGate followed by a Hadamard (or vice versa), or a Hadamard followed by a ZGate (or vice
    // versa), regardless of the existence of controls. (Note that if ZGate is unavailable, we can always use an
    // equivalent ArbitraryPhase.)
    if (Hadamard::gateAvailable(*context_, target(), controls()))
    {
      if (zGateEquivalentAvailable(*context_, target(), controls()) &&
          Hadamard::gateCost(*context_, numControls()) + zGateEquivalentCost(*context_, target(), controls()) < cost())
      {
        return true;
      }
      if (XGate::gateAvailable(*context_, target(), controls()) &&
          Hadamard::gateCost(*context_, numControls()) + XGate::gateCost(*context_, numControls()) < cost())
      {
        return true;
      }
    }
    return false;  // H, X or Z are unavailable, or both XH and HZ are too expensive.
  }

  if (abs(abs(angle()) - 3 * pi / 2) < angleTolerance)  // Angle close to 3pi/2 or -3pi/2. (Different reductions.)
  {
    // Gate must have controls. Gate is equivalent to a controlled Hadamard, a controlled XGate or ZGate (or
    // equivalent), and a ZGate (or equivalent) on the control bits to handle the phase difference.
    assert(numControls() > 0);  // The size of angles for uncontrolled gates should not exceed pi.
    if (Hadamard::gateAvailable(*context_, target(), controls()) && zGateEquivalentAvailable(*context_, controls()))
    {
      if (zGateEquivalentAvailable(*context_, target(), controls()) &&  // Yes, we need this check. (No. of controls.)
          Hadamard::gateCost(*context_, numControls()) + zGateEquivalentCost(*context_, target(), controls()) +  // HZ
          zGateEquivalentCost(*context_, controls()) < cost())  // Z on controls.
      {
        return true;
      }
      if (XGate::gateAvailable(*context_, target(), controls()) &&
          Hadamard::gateCost(*context_, numControls()) + XGate::gateCost(*context_, numControls()) +  // XH
          zGateEquivalentCost(*context_, controls()) < cost())  // Z on controls.
      {
        return true;
      }
    }
  }

  return false;  // Angle does not take a special value.
}


GateSequence YRotation::reduction() const
{
  // This should only be called if a reduction has been found to exist. That allows for some simplifications.
  double angleTolerance = context_->angleTolerance();

  if (abs(angle()) < angleTolerance)  // Angle close to zero.
  {
    // Gate can be eliminated.
    return {};
  }

  if (abs(angle()) > 2 * pi - angleTolerance)  // Angle close to +-2pi. (Angle should NEVER be more than 2pi or less...
  {                                            // ...than -2pi.)
    // Gate can be replaced by a ZGate on the control bits.
    return makeGateSequence(zGateEquivalent(*context_, controls()));
  }

  if (abs(angle() - pi) < angleTolerance)  // Angle close to pi.
  {
    // Gate can be replaced with a YGate and a PhaseGate (or equivalent) on the control bits, if there are any.
    auto firstGate = make_unique<YGate>(*context_, target(), controls());
    if (numControls() == 0)
    {
      return makeGateSequence(std::move(firstGate));
    }
    return makeGateSequence(std::move(firstGate), phaseGateEquivalent(*context_, controls()));
  }

  if (abs(angle() + pi) < angleTolerance)   // Angle close to -pi.
  {
    // Gate can be replaced with a YGate and a PhaseInv (or equivalent) on the control bits, if there are any.
    auto firstGate = make_unique<YGate>(*context_, target(), controls());
    if (numControls() == 0)
    {
      return makeGateSequence(std::move(firstGate));
    }
    return makeGateSequence(std::move(firstGate), phaseInvEquivalent(*context_, controls()));
  }

  if (abs(angle() - pi / 2) < angleTolerance)  // Angle close to pi/2.
  {
    // Gate can be replaced by an XGate followed by a Hadamard, or a Hadamard followed by a ZGate, whichever is cheaper,
    // regardless of the existence of controls, with the new gates having the same controls as the old. (This is because
    // of a fortuitous cancellation of phases.) If XGate and ZGate have the same cost, prefer Hadamard followed by
    // ZGate. If the ZGate is actually an equivalent ArbitraryPhase, still prefer Hadamard followed by the
    // ArbitraryPhase. (The extra degree of freedom might be useful if we do yet more numerical optimization. Yes, we
    // might miss simplification with a neighbouring XGate (or other).)
    auto hadamard = make_unique<Hadamard>(*context_, target(), controls());
    if (zGateEquivalentAvailable(*context_, target(), controls()))
    {
      if (XGate::gateAvailable(*context_, target(), controls()))
      {
        if (XGate::gateCost(*context_, numControls()) < zGateEquivalentCost(*context_, target(), controls()))
        {
          // XH cheaper than HZ.
          return makeGateSequence(make_unique<XGate>(*context_, target(), controls()), std::move(hadamard));
        }
      }
      // XH as expensive as HZ or unavailable.
      return makeGateSequence(std::move(hadamard), zGateEquivalent(*context_, target(), controls()));
    }
    // HZ not available. (XH must be.)
    return makeGateSequence(make_unique<XGate>(*context_, target(), controls()), std::move(hadamard));
  }

  if (abs(angle() + pi / 2) < angleTolerance)  // Angle close to -pi/2.
  {
    // Gate can be replaced by a ZGate followed by a Hadamard, or a Hadamard followed by an XGate, whichever is cheaper,
    // regardless of the existence of controls, with the new gates having the same controls as the old. (This is because
    // of a fortuitous cancellation of phases.) If XGate and ZGate have the same cost, prefer ZGate followed by
    // Hadamard. If the ZGate is actually an equivalent ArbitraryPhase, still prefer ArbitraryPhase followed by
    // Hadamard. (The extra degree of freedom might be useful if we do yet more numerical optimization. Yes, we might
    // miss simplification with a neighbouring XGate (or other).)
    auto hadamard = make_unique<Hadamard>(*context_, target(), controls());
    if (zGateEquivalentAvailable(*context_, target(), controls()))
    {
      if (XGate::gateAvailable(*context_, target(), controls()))
      {
        if (XGate::gateCost(*context_, numControls()) < zGateEquivalentCost(*context_, target(), controls()))
        {
          // HX cheaper than ZH.
          return makeGateSequence(std::move(hadamard), make_unique<XGate>(*context_, target(), controls()));
        }
      }
      // HX as expensive as ZH or unavailable.
      return makeGateSequence(zGateEquivalent(*context_, target(), controls()), std::move(hadamard));
    }
    // ZH not available. (HX must be.)
    return makeGateSequence(std::move(hadamard), make_unique<XGate>(*context_, target(), controls()));
  }

  if (abs(angle() - 3 * pi / 2) < angleTolerance)  // Angle close to 3pi/2.
  {
    // Gate must have controls, as the angle exceeds pi. This case is the same as the previous (-pi/2) case, except that
    // we also need a ZGate on the control bits to account for the difference in phase.
    auto hadamard = make_unique<Hadamard>(*context_, target(), controls());
    auto fixPhase = zGateEquivalent(*context_, controls());
    if (zGateEquivalentAvailable(*context_, target(), controls()))
    {
      if (XGate::gateAvailable(*context_, target(), controls()))
      {
        if (XGate::gateCost(*context_, numControls()) < zGateEquivalentCost(*context_, target(), controls()))
        {
          return makeGateSequence(std::move(hadamard), make_unique<XGate>(*context_, target(), controls()),
                                  std::move(fixPhase));
        }
      }
      return makeGateSequence(zGateEquivalent(*context_, target(), controls()), std::move(hadamard),
                              std::move(fixPhase));
    }
    return makeGateSequence(std::move(hadamard), make_unique<XGate>(*context_, target(), controls()),
                            std::move(fixPhase));
  }

  assert(abs(angle() + 3 * pi / 2) < angleTolerance);  // Angle must be close to -3pi/2.
  {                                              // Added braces mainly just to match the format of the previous 'if's!
    // Gate must have controls, as the angle exceeds pi. This case is the same as the pi/2 case, except that we also
    // need a ZGate on the control bits to account for the difference in phase.
    auto hadamard = make_unique<Hadamard>(*context_, target(), controls());
    auto fixPhase = zGateEquivalent(*context_, controls());
    if (zGateEquivalentAvailable(*context_, target(), controls()))
    {
      if (XGate::gateAvailable(*context_, target(), controls()))
      {
        if (XGate::gateCost(*context_, numControls()) < zGateEquivalentCost(*context_, target(), controls()))
        {
          return makeGateSequence(make_unique<XGate>(*context_, target(), controls()), std::move(hadamard),
                                  std::move(fixPhase));
        }
      }
      return makeGateSequence(std::move(hadamard), zGateEquivalent(*context_, target(), controls()),
                              std::move(fixPhase));
    }
    return makeGateSequence(make_unique<XGate>(*context_, target(), controls()), std::move(hadamard),
                            std::move(fixPhase));
  }
}


long YRotation::canSimplify(const Gate& next) const
{
  return next.canSimplify(*this);  // Double dispatch.
}


GateSequence YRotation::simplification(const Gate& next) const
{
  return next.simplification(*this);  // Double dispatch.
}


bool YRotation::canHSwap(const Gate& next) const
{
  return next.canHSwap(*this);  // Double dispatch.
}


unique_ptr<Gate> YRotation::rightHMoverChange(const Gate& next) const
{
  return next.rightHMoverChange(*this);  // Double dispatch.
}


unique_ptr<Gate> YRotation::leftHMoverChange(const Gate& next) const
{
  return next.leftHMoverChange(*this);  // Double dispatch.
}


bool YRotation::canRSwap(const Gate& next) const
{
  return next.canRSwap(*this);  // Double dispatch.
}


unique_ptr<Gate> YRotation::rightRMoverChange(const Gate& next) const
{
  return next.rightRMoverChange(*this);  // Double dispatch.
}


unique_ptr<Gate> YRotation::leftRMoverChange(const Gate& next) const
{
  return next.leftRMoverChange(*this);  // Double dispatch.
}


bool YRotation::canHSwap(const Hadamard& prev) const
{
  // We can do HY->YH with no change in cost iff we can do YH->HY with no change in cost. Requirements on qbits and gate
  // availability are the same.
  return prev.canHSwap(*this);
}


unique_ptr<Gate> YRotation::rightHMoverChange(const Hadamard& prev) const
{
  return {};  // The Hadamard remains unchanged.
}


unique_ptr<Gate> YRotation::leftHMoverChange(const Hadamard& prev) const
{
  // The gate that X changes to in HY->YH swap is the same (i.e. a YRotation with opposite angle) as in the YH->HY swap.
  return prev.rightHMoverChange(*this);
}


long YRotation::canSimplify(const PhaseGate& prev) const
{
  // Situations where S-YRot simplify are precisely the same as those where YRot-S simplify.
  return prev.canSimplify(*this);
}


GateSequence YRotation::simplification(const PhaseGate& prev) const
{
  // Simplification of S-YRot does not produce the same sequence as simplification of YRot-S. Here the YRotation turns
  // into an XRotation of the same angle.
  return makeGateSequence(make_unique<XRotation>(*context_, target(), controls(), angle()), prev.clone());
  // return {make_unique<XRotation>(*context_, target(), controls, angle()), prev.clone()};
}


bool YRotation::canRSwap(const PhaseGate& prev) const
{
  return prev.canRSwap(*this);
}


unique_ptr<Gate> YRotation::rightRMoverChange(const PhaseGate& prev) const
{
  return prev.leftRMoverChange(*this);
}


unique_ptr<Gate> YRotation::leftRMoverChange(const PhaseGate& prev) const
{
  return prev.rightRMoverChange(*this);
}


long YRotation::canSimplify(const PhaseInv& prev) const
{
  // Situations where (S-1)-YRot simplify are precisely the same as those where YRot-(S-1) simplify.
  return prev.canSimplify(*this);
}


GateSequence YRotation::simplification(const PhaseInv& prev) const
{
  // Simplification of (S-1)-YRot does not produce the same sequence as simplification of YRot-(S-1). Here the YRotation
  // turns into an XRotation of opposite angle.
  return makeGateSequence(make_unique<XRotation>(*context_, target(), controls(), -angle()), prev.clone());
  // return {make_unique<XRotation>(*context_, target(), controls, -angle()), prev.clone()};
}


bool YRotation::canRSwap(const PhaseInv& prev) const
{
  return prev.canRSwap(*this);
}


unique_ptr<Gate> YRotation::rightRMoverChange(const PhaseInv& prev) const
{
  return {};  // PhaseGate is unchanged.
}


unique_ptr<Gate> YRotation::leftRMoverChange(const PhaseInv& prev) const
{
  // YRotation turns into a XRotation of the opposite angle. Note that this is NOT the same as when the YRotation
  // precedes the PhaseInv.
  return make_unique<XRotation>(*context_, prev.target(), prev.controls(), -angle());
}


bool YRotation::canRSwap(const XGate& prev) const
{
  return prev.canRSwap(*this);
}


unique_ptr<Gate> YRotation::rightRMoverChange(const XGate& prev) const
{
  return prev.leftRMoverChange(*this);
}


unique_ptr<Gate> YRotation::leftRMoverChange(const XGate& prev) const
{
  return prev.rightRMoverChange(*this);
}


long YRotation::canSimplify(const YGate& prev) const
{
  // Gate order is unimportant - they commute whenever a simplification exists. To avoid code duplication, call the
  // function with the gates reversed.
  return prev.canSimplify(*this);
}


GateSequence YRotation::simplification(const YGate& prev) const
{
  // Gate order is unimportant - they commute whenever a simplification exists. To avoid code duplication, call the
  // function with the gates reversed.
  return prev.simplification(*this);
}


long YRotation::canSimplify(const ZGate& prev) const
{
  // Gate order is unimportant - they commute whenever a simplification exists. To avoid code duplication, call the
  // function with the gates reversed.
  return prev.canSimplify(*this);
}


GateSequence YRotation::simplification(const ZGate& prev) const
{
  // Gate order is unimportant - they commute whenever a simplification exists. To avoid code duplication, call the
  // function with the gates reversed.
  return prev.simplification(*this);
}


bool YRotation::canRSwap(const ZGate& prev) const
{
  return prev.canRSwap(*this);
}


unique_ptr<Gate> YRotation::rightRMoverChange(const ZGate& prev) const
{
  return prev.leftRMoverChange(*this);
}


unique_ptr<Gate> YRotation::leftRMoverChange(const ZGate& prev) const
{
  return prev.rightRMoverChange(*this);
}


long YRotation::canSimplify(const YRotation& prev) const
{
  // Provided both target and control bits agree, this gives another XRotation - just add the angles modulo 2pi.
  if (target() == prev.target() && controls() == prev.controls())
  {
    return cost();
  }
  return 0;
}


GateSequence YRotation::simplification(const YRotation& prev) const
{
  return makeGateSequence(make_unique<YRotation>(*context_, target(), controls(), angle() + prev.angle()));
//  return {make_unique<YRotation>(*context_, target(), controls(), angle() + prev.angle())};
}


const GateSimulator& YRotation::part_simulator(int partIndex) const
{
  switch (partIndex)
  {
    case 0:
      return const_sim;
    case 1:
      return cos_sim;
    case 2:
      return sin_sim;
    default:
      throw std::invalid_argument("Invalid partIndex passed to XRotation::part_simulator().");
  }
}


double YRotation::thetaToAngle(double theta) const
{
  return 2 * theta;
}


vector<vector<cmplx> > YRotation::matrix_() const
{
  using std::cos, std::sin;
  double angle = this->angle_;  // Yes, the 'this' is necessary! (Apparently not on Visual Studio, but that's a...
  return {{cos(angle / 2), sin(angle / 2)}, {-sin(angle / 2), cos(angle / 2)}};  // ...'non-conformance'.)
}

//----------------------------------------------------------------------------------------------------------------------

// Set up ZRotation class simulators for just the constant part, the cosine part and the sine part of the rotation
// matrix. Simulator for the constant part must handle the case where a control bit is not set. The other two must
// ignore this case - hence the 'false'.
const SparseGateSimulator ZRotation::const_sim(sparse({{0, 0}, {0, 0}}));
const SparseGateSimulator ZRotation::cos_sim(sparse({{1, 0}, {0, 1}}), false);
const SparseGateSimulator ZRotation::sin_sim(sparse({{constants::i, 0}, {0, -constants::i}}), false);


ZRotation::ZRotation(const CircuitContext& context) :
SingleTargetGate(context),
DiagonalGate(context),
RotationGate(context)
{
}


ZRotation::ZRotation(const CircuitContext& context, int target, const Controls& controls, double angle) :
SingleTargetGate(context, target, controls),
DiagonalGate(context, target, controls),
RotationGate(context, target, controls, angle)
{
}


int ZRotation::gateTypeId() const
{
  return gate_type_id;
}


unique_ptr<Gate> ZRotation::clone() const
{
  return make_unique<ZRotation>(*this);  // Will get converted to unique_ptr<Gate>.
}


ZRotation* ZRotation::rawClone() const
{
  return new ZRotation(*this);
}


std::string ZRotation::name() const
{
  return "ZRotation";
}


bool ZRotation::gateAvailable(const CircuitContext& context, int target, const Controls& controls)  // Static
{
  // Static function for determining the availability of a Gate without needing to create one. (The corresponding
  // non-static function is useful when we have a Gate but don't necessarily know the precise gate type.)
  // (We don't need to provide the target bit at present - just the number of controls. However, in future we may have
  // more complicated rules on availability.)
  return context.gateTypeAvailable(gate_type_id, controls.numControls());
}


long ZRotation::gateCost(const CircuitContext& context, int numControls)  // Static
{
  return context.gateCost(gate_type_id, numControls);
}


bool ZRotation::canReduce() const
{
  double angleTolerance = context_->angleTolerance();

  if (abs(angle()) < angleTolerance)  // Angle close to zero.
  {
    // Gate essentially does nothing, so we can eliminate it.
    return true;
  }

  if (abs(angle()) > 2 * pi - angleTolerance)  // Angle close to +-2pi. (Angle should NEVER be more than 2pi or less...
  {                                            // ...than -2pi.)
    // Gate (which must have controls) is equivalent to a ZGate on the control bits. If this (or an equivalent) is
    // available and cheaper than the XRotation, we take it.
    assert(numControls() > 0);  // The size of angles for uncontrolled gates should not exceed pi.
    return zGateEquivalentAvailable(*context_, controls()) && zGateEquivalentCost(*context_, controls()) < cost();
  }

  if (abs(angle() - pi) < angleTolerance)  // Angle close to pi.
  {
    // If uncontrolled, the gate is in effect just an ZGate. If there are controls, then a PhaseGate is also required on
    // the control bits.
    // (The likelihood that the sum of the costs of an ZGate and a PhaseGate exceed that of the original ZRotation is
    // motivation for the future consideration of 'reduction' functions that consider pairs of gates. Then we might
    // detect that the cost increase is only temporary, eliminated by the cancellation (or simplification) of the
    // introduced PhaseGate with a neighbouring gate.)
    if (numControls() == 0)
    {
      return zGateEquivalentAvailable(*context_, target(), controls()) &&
             zGateEquivalentCost(*context_, target(), controls()) < cost();
    }
    return zGateEquivalentAvailable(*context_, target(), controls()) &&
           phaseGateEquivalentAvailable(*context_, controls()) &&
           zGateEquivalentCost(*context_, target(), controls()) +
                   phaseGateEquivalentCost(*context_, controls()) < cost();
  }

  if (abs(angle() + pi) < angleTolerance)  // Angle close to -pi.
  {
    // If uncontrolled, the gate is in effect just a ZGate. If there are controls, then a PhaseInv is also required on
    // the control bits.
    if (numControls() == 0)
    {
      return zGateEquivalentAvailable(*context_, target(), controls()) &&
             zGateEquivalentCost(*context_, target(), controls()) < cost();
    }
    return zGateEquivalentAvailable(*context_, target(), controls()) &&
           phaseInvEquivalentAvailable(*context_, controls()) &&
           zGateEquivalentCost(*context_, target(), controls()) +
                   phaseInvEquivalentCost(*context_, controls()) < cost();
  }

  if (abs(angle() - pi / 2) < angleTolerance)  // Angle close to pi/2.
  {
    // Gate is equivalent to an PhaseInv with a PiByEight on the control bits, if there are any.
    if (phaseInvEquivalentAvailable(*context_, target(), controls()))
    {
      if (numControls() == 0)
      {
        return phaseInvEquivalentCost(*context_, target(), controls()) < cost();
      }
      return piByEightEquivalentAvailable(*context_, controls()) &&
             phaseInvEquivalentCost(*context_, target(), controls()) +
                     piByEightEquivalentCost(*context_, controls()) < cost();
    }
    return false;  // PhaseInv (or equivalent) not available.
  }

  if (abs(angle() + pi / 2) < angleTolerance)  // Angle close to -pi/2
  {
    // Gate is equivalent to a PhaseGate with a PiByEightInv on the control bits, if there are any.
    if (phaseGateEquivalentAvailable(*context_, target(), controls()))
    {
      if (numControls() == 0)
      {
        return phaseGateEquivalentCost(*context_, target(), controls()) < cost();
      }
      return piByEightInvEquivalentAvailable(*context_, controls()) &&
             phaseGateEquivalentCost(*context_, target(), controls()) +
                     piByEightInvEquivalentCost(*context_, controls()) < cost();
    }
  }

  if (abs(angle() - 3 * pi / 2) < angleTolerance)  // Angle close to 3pi/2.
  {
    // Gate must have controls. Gate is equivalent to a PhaseGate with an ArbitraryPhase on the control bits.
    // (The ArbitraryPhase is equivalent to a PiByEightInv and a ZGate, or a PiByEight and a PhaseGate. If we were to
    // consider those gates in the group generated by the PiByEight as atomic gates, then we would have a simple
    // unparameterized gate to use instead of ArbitraryPhase. We don't, at present, consider using a pair of phase style
    // gates, though perhaps we should.)
    assert(numControls() > 0);  // The size of angles for uncontrolled gates should not exceed pi.
    return phaseGateEquivalentAvailable(*context_, target(), controls()) &&
           ArbitraryPhase::gateAvailable(*context_, controls()) &&
           phaseGateEquivalentCost(*context_, target(), controls()) +
                   ArbitraryPhase::gateCost(*context_, numControls() - 1) < cost();
  }

  if (abs(angle() + 3 * pi / 2) < angleTolerance)  // Angle close to -3pi/2.
  {
    // Gate must have controls. Gate is equivalent to a PhaseInv with an ArbitraryPhase on the control bits.
    // (The ArbitraryPhase is equivalent to a PiByEight and a ZGate, or a PiByEightInv and a PhaseInv. If we were to
    // consider those gates in the group generated by the PiByEight as atomic gates, then we would have a simple
    // unparameterized gate to use instead of ArbitraryPhase. We don't, at present, consider using a pair of phase style
    // gates, though perhaps we should.)
    assert(numControls() > 0);  // The size of angles for uncontrolled gates should not exceed pi.
    return phaseInvEquivalentAvailable(*context_, target(), controls()) &&
           ArbitraryPhase::gateAvailable(*context_, controls()) &&
           phaseInvEquivalentCost(*context_, target(), controls()) +
                  ArbitraryPhase::gateCost(*context_, numControls() - 1) < cost();
  }

  if (abs(angle() - pi / 4) < angleTolerance)  // Angle close to pi/4.
  {
    // Gate is equivalent to a PiByEightInv with an ArbitraryPhase of angle pi/8 on the control bits, if there are any.
    // Note that if PiByEightInv is not available, but ArbitraryPhase is, then the reduction being considered is
    // actually the reduction that could be applied to ANY ZRotation, i.e. 'reduction' to an Arbitrary phase on target
    // and controls, and an ArbitraryPhase on controls only. If this is really a cheaper option then (at present at
    // least) the ZRotation is a redundant gate which perhaps ought to be removed, or at least a warning could be
    // provided to the user.
    if (piByEightInvEquivalentAvailable(*context_, target(), controls()))
    {
      if (numControls() == 0)
      {
        return piByEightInvEquivalentCost(*context_, target(), controls()) < cost();
      }
      return ArbitraryPhase::gateAvailable(*context_, controls()) &&
             piByEightInvEquivalentCost(*context_, target(), controls()) +
                     ArbitraryPhase::gateCost(*context_, numControls() - 1) < cost();
    }
  }

  if (abs(angle() + pi / 4) < angleTolerance)  // Angle close to -pi/4.
  {
    // Gate is equivalent to a PiByEight with an ArbitraryPhase of angle -pi/8 on the control bits, if there are any.
    // Note that if PiByEight is not available, but ArbitraryPhase is, then the reduction being considered is actually
    // the reduction that could be applied to ANY ZRotation, i.e. 'reduction' to an Arbitrary phase on target and
    // controls, and an ArbitraryPhase on controls only. If this is really a cheaper option then (at present at least)
    // the ZRotation is a redundant gate which perhaps ought to be removed, or at least a warning could be provided to
    // the user.
    if (piByEightEquivalentAvailable(*context_, target(), controls()))
    {
      if (numControls() == 0)
      {
        return piByEightEquivalentCost(*context_, target(), controls()) < cost();
      }
      return ArbitraryPhase::gateAvailable(*context_, controls()) &&
             piByEightEquivalentCost(*context_, target(), controls()) +
                     ArbitraryPhase::gateCost(*context_, numControls() - 1) < cost();
    }
  }

  if (abs(angle() - 3 * pi / 4) < angleTolerance)  // Angle close to 3pi/4.
  {
    // Gate is equivalent to a PhaseInv and PiByEightInv (or ZGate and PiByEight), with an ArbitraryPhase of angle 3pi/8
    // on the control bits. At the moment, we don't consider reductions involving sequences such as PhaseInv-
    // PiByEightInv (see control bits in the 3pi/2 case). (We might wish to include gates in the group generated by
    // PiByEight in future?) Furthermore, if using an ArbitraryPhase instead results in a cheaper circuit, then the
    // ZGate (with the number of controls as here) is actually redundant - it can ALWAYS be replaced by ArbitraryPhases
    // in this way.
  }

  if (abs(angle() + 3 * pi / 4) < angleTolerance)  // Angle close to -3pi/4.
  {
    // Gate is equivalent to a PhaseGate and PiByEight (or ZGate and PiByEightInv), with an ArbitraryPhase of angle
    // -3pi/8 on the control bits. At the moment, we don't consider reductions involving sequences such as PhaseInv-
    // PiByEightInv (see control bits in the 3pi/2 case). (We might wish to include gates in the group generated by
    // PiByEight in future?) Furthermore, if using an ArbitraryPhase instead results in a cheaper circuit, then the
    // ZGate (with the number of controls as here) is actually redundant - it can ALWAYS be replaced by ArbitraryPhases
    // in this way.
  }

  if (abs(angle() - 5 * pi / 4) < angleTolerance)  // Angle close to 5pi/4.
  {
    // Gate is equivalent to a PhaseGate and PiByEight (or ZGate and PiByEightInv), with an ArbitraryPhase of angle
    // -5pi/8 on the control bits. At the moment, we don't consider reductions involving sequences such as PhaseInv-
    // PiByEightInv (see control bits in the 3pi/2 case). (We might wish to include gates in the group generated by
    // PiByEight in future?) Furthermore, if using an ArbitraryPhase instead results in a cheaper circuit, then the
    // ZGate (with the number of controls as here) is actually redundant - it can ALWAYS be replaced by ArbitraryPhases
    // in this way.
    assert(numControls() > 0);  // The size of angles for uncontrolled gates should not exceed pi.
  }

  if (abs(angle() + 5 * pi / 4) < angleTolerance)  // Angle close to -5pi/4.
  {
    // Gate is equivalent to a PhaseInv and PiByEightInv, with an ArbitraryPhase of angle -5pi/8 on the control bits. At
    // the moment, we don't consider reductions involving sequences such as PhaseInv PiByEightInv (see control bits in
    // the 3pi/2 case). (We might wish to include gates in the group generated by PiByEight in future?) Furthermore, if
    // using an ArbitraryPhase instead results in a cheaper circuit, then the ZGate (with the number of controls as
    // here) is actually redundant - it can ALWAYS be replaced by ArbitraryPhases in this way.
    assert(numControls() > 0);  // The size of angles for uncontrolled gates should not exceed pi.
  }

  if (abs(angle() - 7 * pi / 4) < angleTolerance)  // Angle close to 7pi/4.
  {
    // Gate must have controls. Gate is equivalent to a PiByEight with an ArbitraryPhase of angle 7pi/8 on the control
    // bits, if there are any. Note that if PiByEight is not available, but ArbitraryPhase is, then the reduction being
    // considered is actually the reduction that could be applied to ANY ZRotation, i.e. 'reduction' to an Arbitrary
    // phase on target and controls, and an ArbitraryPhase on controls only. If this is really a cheaper option then
    // (at present at least) the ZRotation is a redundant gate which perhaps ought to be removed, or at least a warning
    // could be provided to the user.
    assert(numControls() > 0);  // The size of angles for uncontrolled gates should not exceed pi.
    return piByEightEquivalentAvailable(*context_, target(), controls()) &&
           ArbitraryPhase::gateAvailable(*context_, controls()) &&
           piByEightEquivalentCost(*context_, target(), controls()) +
                   ArbitraryPhase::gateCost(*context_, numControls() - 1) < cost();
  }

  if (abs(angle() + 7 * pi / 4) < angleTolerance)  // Angle close to -7pi/4.
  {
    // Gate must have controls. Gate is equivalent to a PiByEightInv with an ArbitraryPhase of angle -7pi/8 on the
    // control bits, if there are any. Note that if PiByEightInv is not available, but ArbitraryPhase is, then the
    // reduction being considered is actually the reduction that could be applied to ANY ZRotation, i.e. 'reduction' to
    // an Arbitrary phase on target and controls, and an ArbitraryPhase on controls only. If this is really a cheaper
    // option then (at present at least) the ZRotation is a redundant gate which perhaps ought to be removed, or at
    // least a warning could be provided to the user.
    assert(numControls() > 0);  // The size of angles for uncontrolled gates should not exceed pi.
    return piByEightInvEquivalentAvailable(*context_, target(), controls()) &&
           ArbitraryPhase::gateAvailable(*context_, controls()) &&
           piByEightInvEquivalentCost(*context_, target(), controls()) +
                   ArbitraryPhase::gateCost(*context_, numControls() - 1) < cost();
  }

  // Angle does not take a special value. While a general ZRotation can be 'reduced' to one or two ArbitraryPhases, we
  // do not consider this option, since then the ZRotation (with this number of controls) is actually redundant.
  return false;
}


GateSequence ZRotation::reduction() const
{
  // This should only be called if a reduction has been found to exist. That allows for some simplifications.
  double angleTolerance = context_->angleTolerance();

  if (abs(angle()) < angleTolerance)  // Angle close to zero.
  {
    // Gate can be eliminated.
    return {};
  }

  if (abs(angle()) > 2 * pi - angleTolerance)  // Angle close to +-2pi. (Angle should NEVER be more than 2pi or less...
  {                                            // ...than -2pi.)
    // Gate can be replaced by a ZGate on the control bits.
    return makeGateSequence(zGateEquivalent(*context_, controls()));
  }

  if (abs(angle() - pi) < angleTolerance)  // Angle close to pi.
  {
    // Gate can be replaced with a ZGate (or equivalent) and a PhaseGate (or equivalent) on the control bits, if there
    // are any.
    auto firstGate = zGateEquivalent(*context_, target(), controls());
    if (numControls() == 0)
    {
      return makeGateSequence(std::move(firstGate));
    }
    return makeGateSequence(std::move(firstGate), phaseGateEquivalent(*context_, controls()));
  }

  if (abs(angle() + pi) < angleTolerance)   // Angle close to -pi.
  {
    // Gate can be replaced with a ZGate (or equivalent) and a PhaseInv (or equivalent) on the control bits, if there
    // are any.
    auto firstGate = zGateEquivalent(*context_, target(), controls());
    if (numControls() == 0)
    {
      return makeGateSequence(std::move(firstGate));
    }
    return makeGateSequence(std::move(firstGate), phaseInvEquivalent(*context_, controls()));
  }

  if (abs(angle() - pi / 2) < angleTolerance)  // Angle close to pi/2.
  {
    // Gate can be replaced by a PhaseInv (or equivalent) with a PiByEight (or equivalent) on the control bits.
    auto firstGate = phaseInvEquivalent(*context_, target(), controls());
    if (numControls() == 0)
    {
      return makeGateSequence(std::move(firstGate));
    }
    return makeGateSequence(std::move(firstGate), piByEightEquivalent(*context_, controls()));
  }

  if (abs(angle() + pi / 2) < angleTolerance)  // Angle close to -pi/2.
  {
    // Gate can be replaced by a PhaseGate (or equivalent) with a PiByEightInv (or equivalent) on the control bits, if
    // there are any.
    auto firstGate = phaseGateEquivalent(*context_, target(), controls());
    if (numControls() == 0)
    {
      return makeGateSequence(std::move(firstGate));
    }
    return makeGateSequence(std::move(firstGate), piByEightInvEquivalent(*context_, controls()));
  }

  if (abs(angle() - 3 * pi / 2) < angleTolerance)  // Angle close to 3pi/2.
  {
    // Gate must have controls, as the angle exceeds pi. Gate can be replaced by a PhaseGate (or equivalent) with an
    // ArbitraryPhase of angle 3pi/4 on the control bits.
    return makeGateSequence(phaseGateEquivalent(*context_, target(), controls()),
                            make_unique<ArbitraryPhase>(*context_, controls(), 3 * pi / 4));
  }

  if (abs(angle() + 3 * pi / 2) < angleTolerance)  // Angle must be close to -3pi/2.
  {
    // Gate must have controls, as the angle exceeds pi. Gate can be replaced by a PhaseInv (or equivalent) with an
    // ArbitraryPhase of angle -3pi/4 on the control bits.
    return makeGateSequence(phaseInvEquivalent(*context_, target(), controls()),
                            make_unique<ArbitraryPhase>(*context_, controls(), -3 * pi / 4));
  }

  if (abs(angle() - pi / 4) < angleTolerance)  // Angle close to pi/4.
  {
    // Gate can be replaced by a PiByEightInv (or equivalent) and an ArbitraryPhase of angle pi/8 on the control bits,
    // if there are any.
    auto firstGate = piByEightInvEquivalent(*context_, target(), controls());
    if (numControls() == 0)
    {
      return makeGateSequence(std::move(firstGate));
    }
    return makeGateSequence(std::move(firstGate), make_unique<ArbitraryPhase>(*context_, controls(), pi / 8));
  }

  if (abs(angle() + pi / 4) < angleTolerance)  // Angle close to -pi/4.
  {
    // Gate can be replaced by a PiByEight (or equivalent) and an ArbitraryPhase of angle -pi/8 on the control bits, if
    // there are any.
    auto firstGate = piByEightEquivalent(*context_, target(), controls());
    if (numControls() == 0)
    {
      return makeGateSequence(std::move(firstGate));
    }
    return makeGateSequence(std::move(firstGate), make_unique<ArbitraryPhase>(*context_, controls(), -pi / 8));
  }

  if (abs(angle() - 7 * pi /4) < angleTolerance)  // Angle close to 7pi/4.
  {
    // Gate must have controls, as the angle exceeds pi. Gate can be replaced by a PiByEight (or equivalent) and an
    // ArbitraryPhase of angle 7pi/8 on the control bits.
    return makeGateSequence(piByEightEquivalent(*context_, target(), controls()),
                            make_unique<ArbitraryPhase>(*context_, controls(), 7 * pi / 8));
  }

  assert(abs(angle() + 7 * pi / 4) < angleTolerance);  // Angle close to -7pi/4.
  {                                     // Added braces mainly just to match the format of the previous 'if's!
    // Gate must have controls, as the size of the angle exceeds pi. Gate can be replaced by a PiByEightInv (or
    // equivalent) and an ArbitraryPhase of angle -7pi/8 on the control bits.
    return makeGateSequence(piByEightInvEquivalent(*context_, target(), controls()),
                            make_unique<ArbitraryPhase>(*context_, controls(), -7 * pi / 8));
  }
}


long ZRotation::canSimplify(const Gate& next) const
{
  return next.canSimplify(*this);  // Double dispatch.
}


GateSequence ZRotation::simplification(const Gate& next) const
{
  return next.simplification(*this);  // Double dispatch.
}


bool ZRotation::canHSwap(const Gate& next) const
{
  return next.canHSwap(*this);  // Double dispatch.
}


unique_ptr<Gate> ZRotation::rightHMoverChange(const Gate& next) const
{
  return next.rightHMoverChange(*this);
}


unique_ptr<Gate> ZRotation::leftHMoverChange(const Gate& next) const
{
  return next.leftHMoverChange(*this);
}


bool ZRotation::canRSwap(const Gate& next) const
{
  return next.canRSwap(*this);  // Double dispatch.
}


unique_ptr<Gate> ZRotation::rightRMoverChange(const Gate& next) const
{
  return next.rightRMoverChange(*this);  // Double dispatch.
}


unique_ptr<Gate> ZRotation::leftRMoverChange(const Gate& next) const
{
  return next.leftRMoverChange(*this);  // Double dispatch.
}


long ZRotation::canSimplify(const Hadamard& prev) const
{
  // Situations where we can simplify H-ZRot are exactly the same as those where we can simplify ZRot-H.
  return prev.canSimplify(*this);
}


GateSequence ZRotation::simplification(const Hadamard& prev) const
{
  // H-ZRot does not simplify to the same sequence as ZRot-H
  return makeGateSequence(make_unique<XRotation>(*context_, target(), controls(), angle()), prev.clone());
//  return {make_unique<XRotation>(*context_, target(), controls(), angle()), prev.clone()};
}


bool ZRotation::canHSwap(const Hadamard& prev) const
{
  return prev.canHSwap(*this);
}


unique_ptr<Gate> ZRotation::rightHMoverChange(const Hadamard& prev) const
{
  return {};  // The Hadamard remains unchanged.
}


unique_ptr<Gate> ZRotation::leftHMoverChange(const Hadamard& prev) const
{
  // The gate that Z changes to in HZ->XH swap is the same (i.e. an XRotation) as in the ZH->HX swap.
  return prev.rightHMoverChange(*this);
}


long ZRotation::canSimplify(const PiByEight& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.canSimplify(*this);
}


GateSequence ZRotation::simplification(const PiByEight& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.simplification(*this);
}


long ZRotation::canSimplify(const PiByEightInv& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.canSimplify(*this);
}


GateSequence ZRotation::simplification(const PiByEightInv& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.simplification(*this);
}


long ZRotation::canSimplify(const PhaseGate& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.canSimplify(*this);
}


GateSequence ZRotation::simplification(const PhaseGate& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.simplification(*this);
}


long ZRotation::canSimplify(const PhaseInv& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.canSimplify(*this);
}


GateSequence ZRotation::simplification(const PhaseInv& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.simplification(*this);
}


bool ZRotation::canRSwap(const XGate& prev) const
{
  return prev.canRSwap(*this);
}


unique_ptr<Gate> ZRotation::rightRMoverChange(const XGate& prev) const
{
  return prev.leftRMoverChange(*this);
}


unique_ptr<Gate> ZRotation::leftRMoverChange(const XGate& prev) const
{
  return prev.rightRMoverChange(*this);
}


bool ZRotation::canRSwap(const YGate& prev) const
{
  return prev.canRSwap(*this);
}


unique_ptr<Gate> ZRotation::rightRMoverChange(const YGate& prev) const
{
  return prev.leftRMoverChange(*this);
}


unique_ptr<Gate> ZRotation::leftRMoverChange(const YGate& prev) const
{
  return prev.rightRMoverChange(*this);
}


long ZRotation::canSimplify(const ZGate& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.canSimplify(*this);
}


GateSequence ZRotation::simplification(const ZGate& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.simplification(*this);
}


long ZRotation::canSimplify(const ZRotation& prev) const
{
  // Provided both target and control bits agree, this gives another XRotation - just add the angles modulo 2pi.
  if (target() == prev.target() && controls() == prev.controls())
  {
    return cost();
  }
  return 0;
}


GateSequence ZRotation::simplification(const ZRotation& prev) const
{
  return makeGateSequence(make_unique<ZRotation>(*context_, target(), controls(), angle() + prev.angle()));
//  return {make_unique<ZRotation>(*context_, target(), controls(), angle() + prev.angle())};
}


long ZRotation::canSimplify(const ArbitraryPhase& prev) const
{
  // These gates will simplify if the same set of qbits is involved in each gate. (Recall that one may swap the target
  // of an ArbitraryPhase with a control, without changing the gate's effect.) If both gates have no controls then they
  // simplify to a single ZRotation, or a single ArbitraryPhase, up to an overall phase which we can ignore. If both
  // have controls then the phase difference becomes important because it only applies when the control bits are set. In
  // this case, the gates simplify to a ZRotation or an ArbitraryPhase as before, along with an ArbitraryPhase on the
  // control bits. (Note that we do need to check the availability of the ArbitraryPhase gate across the control qbits,
  // as it is possible that this has no controls in a situation where only ArbitraryPhase gates with at least one
  // control are permitted.)
  //
  // The choice for the main replacement gate currently depends on the relative costs of the gate types. Whichever is
  // cheaper is chosen. If both have the same cost, this function picks the ArbitraryPhase. Note that when a ZRotation
  // is selected, its target must match that of the original ZRotation. The auxiliary ArbitraryPhase across the control
  // bits must be across the controls of the original ZRotation in either case.
  //
  // (If there is only one control bit on the original gates, then the additional ArbitraryPhase, having no control
  // bits, could be replaced with a ZRotation (of opposite angle). Indeed, this is the case in many of our
  // simplification routines. Of course, if there is more than one control bit then the notion of replacing an
  // ArbitraryPhase with a ZRotation becomes more complicated.)
  if (allQbitsMatch(prev))
  {
    auto mainGateCost = std::min(ZRotation::gateCost(*context_, numControls()),
                                 ArbitraryPhase::gateCost(*context_, numControls()));
    if (numControls() == 0)
    {
      return cost() + prev.cost() - mainGateCost;
    }
    if (ArbitraryPhase::gateAvailable(*context_, controls()))
    {
      return cost() + prev.cost() - mainGateCost - ArbitraryPhase::gateCost(*context_, numControls() - 1);
    }
  }
  return 0;
}


GateSequence ZRotation::simplification(const ArbitraryPhase& prev) const
{
  // This is only called if a simplification exists.
  unique_ptr<Gate> mainGate;
  bool zCheaper = (ZRotation::gateCost(*context_, numControls()) < ArbitraryPhase::gateCost(*context_, numControls()));
  if (zCheaper)
  {
    mainGate = make_unique<ZRotation>(*context_, target(), controls(), angle() - prev.angle());
  }
  else
  {
    mainGate = make_unique<ArbitraryPhase>(*context_, target(), controls(), prev.angle() - angle());
  }

  if (numControls() == 0)
  {
    return makeGateSequence(std::move(mainGate));
//    return {std::move(mainGate)};
  }

  unique_ptr<Gate> secondGate;
  if (zCheaper)
  {
    secondGate = make_unique<ArbitraryPhase>(*context_, controls(), prev.angle() / 2);
  }
  else
  {
    secondGate = make_unique<ArbitraryPhase>(*context_, controls(), angle() / 2);
  }
  return makeGateSequence(std::move(mainGate), std::move(secondGate));
//  return {std::move(mainGate), std::move(secondGate)};
}


const GateSimulator& ZRotation::part_simulator(int partIndex) const
{
  switch (partIndex)
  {
    case 0:
      return const_sim;
    case 1:
      return cos_sim;
    case 2:
      return sin_sim;
    default:
      throw std::invalid_argument("Invalid partIndex passed to XRotation::part_simulator().");
  }
}


double ZRotation::thetaToAngle(double theta) const
{
  return 2 * theta;
}


vector<vector<cmplx> > ZRotation::matrix_() const
{
  using std::exp, constants::i;
  double angle = this->angle_;  // Yes, the 'this' is necessary! (Apparently not on Visual Studio, but that's a...
  return {{exp(i * angle / 2.0), 0}, {0, exp(-i * angle / 2.0)}};  // ..'non-conformance'.)
}

//----------------------------------------------------------------------------------------------------------------------

// Set up ArbitraryPhase class simulators for just the constant part, the cosine part and the sine part of the rotation
// matrix. Simulator for  the constant part must handle the case where a control bit is not set. The other two must
// ignore this case - hence the 'false'.
// Warning: while the XRotation, YRotation and ZRotation matrices contain cosines and sines of angle/2, here we have
// sine and cosine of just 'angle'.
const SparseGateSimulator ArbitraryPhase::const_sim(sparse({{1, 0}, {0, 0}}));
const SparseGateSimulator ArbitraryPhase::cos_sim(sparse({{0, 0}, {0, 1}}), false);
const SparseGateSimulator ArbitraryPhase::sin_sim(sparse({{0, 0}, {0, constants::i}}), false);


ArbitraryPhase::ArbitraryPhase(const CircuitContext& context) :
SingleTargetGate(context),
PhaseTypeGate(context),
RotationGate(context)
{
  //  // Gate has not been initialized, so no need to call 'create_simulators()' yet. Pointers to simulators left null.
}


ArbitraryPhase::ArbitraryPhase(const CircuitContext& context, int target, const Controls& controls, double angle) :
SingleTargetGate(context, target, controls),
PhaseTypeGate(context, target, controls),
RotationGate(context, target, controls, angle)
{
  // Controlled ArbPhase has a narrower range (-pi to pi) than other (default) rotation gates (-2pi to 2pi).
  angle_ = remainder(angle, 2 * pi);
}


ArbitraryPhase::ArbitraryPhase(const CircuitContext& context, const Controls& controls, double angle) :  // Target is...
SingleTargetGate(context),                                                                // ...any one of the controls.
PhaseTypeGate(context, controls),
RotationGate(context)
{
  angle_ = remainder(angle, 2 * pi);
}


int ArbitraryPhase::gateTypeId() const
{
  return gate_type_id;
}


unique_ptr<Gate> ArbitraryPhase::clone() const
{
  return make_unique<ArbitraryPhase>(*this);  // Will get converted to unique_ptr<Gate>.
}


ArbitraryPhase* ArbitraryPhase::rawClone() const
{
  return new ArbitraryPhase(*this);
}


void ArbitraryPhase::setParameter(int paramNum, double value)
{
  // The (controlled) ArbitraryPhase has a narrower range for the angle parameter than the typical RotationGate, so we
  // provide this override.

  assert(paramNum == 0);  // ArbitraryPhase should have only one parameter - the angle.
  angle_ = remainder(value, 2 * pi);  // Gives us an angle between -pi and pi.
}


std::string ArbitraryPhase::name() const
{
  return "ArbPhase";  // Was just "Phase" - the same as the S gate. However, this caused ambiguity on output of...
}                     // ...CircuitContexts.


bool ArbitraryPhase::gateAvailable(const CircuitContext& context, int target, const Controls& controls)  // Static
{
  // Static function for determining the availability of a Gate without needing to create one. (The corresponding
  // non-static function is useful when we have a Gate but don't necessarily know the precise gate type.)
  // (We don't need to provide the target bit at present - just the number of controls. However, in future we may have
  // more complicated rules on availability.)
  // NAME: Called 'gateAvailable' just to avoid clashing with 'available' - the non-static version
  return context.gateTypeAvailable(gate_type_id, controls.numControls());
}


bool ArbitraryPhase::gateAvailable(const CircuitContext& context, const Controls& controls)  // Static
{
  // Static function for determining the availability of a Gate without needing to create one. (The corresponding
  // non-static function is useful when we have a Gate but don't necessarily know the precise gate type.)
  // Parameter 'controls' is actually the controls AND the target. Hence the number of actual controls for the desired
  // gate is numControls() - 1.
  assert(controls.numControls() > 0);
  return context.gateTypeAvailable(gate_type_id, controls.numControls() - 1);
}


long ArbitraryPhase::gateCost(const CircuitContext& context, int numControls)  // Static
{
  return context.gateCost(gate_type_id, numControls);
}


bool ArbitraryPhase::canReduce() const
{
  // (Note that when we seek to replace an ArbitraryPhase gate of appropriate angle with, for example, a ZGate, we
  // check the availability of the ZGate only, not for a zGateEquivalent. The only gate equivalent to the ZGate would be
  // the ArbitraryPhase we are aiming to reduce!)
  double angleTolerance = context_->angleTolerance();

  if (abs(angle()) < angleTolerance)  // Angle close to zero.
  {
    // Gate essentially does nothing, so we can eliminate it.
    return true;
  }

  if (abs(abs(angle()) - pi) < angleTolerance)  // Angle close to +-pi.
  {
    // The gate is in effect just an ZGate.
    // (If we were to reintroduce the redundancy check, making redundant a ZGate that is more expensive than an
    // ArbitraryPhase, we should still keep the cost check here. (Defensive.))
    return ZGate::gateAvailable(*context_, target(), controls()) && ZGate::gateCost(*context_, numControls()) < cost();
  }

  if (abs(angle() - pi / 2) < angleTolerance)  // Angle close to pi/2.
  {
    // Gate is in effect just a PhaseGate.
    return PhaseGate::gateAvailable(*context_, target(), controls()) &&
           PhaseGate::gateCost(*context_, numControls()) < cost();
  }

  if (abs(angle() + pi / 2) < angleTolerance)  // Angle close to -pi/2
  {
    // Gate is in effect just a PhaseInv.
    return PhaseInv::gateAvailable(*context_, target(), controls()) &&
           PhaseInv::gateCost(*context_, numControls()) < cost();
  }

  if (abs(angle() - pi / 4) < angleTolerance)  // Angle close to pi/4.
  {
    // Gate is in effect just a PiByEight.
    return PiByEight::gateAvailable(*context_, target(), controls()) &&
           PiByEight::gateCost(*context_, numControls()) < cost();
  }

  if (abs(angle() + pi / 4) < angleTolerance)  // Angle close to -pi/4.
  {
    // Gate is in effect just a PiByEightInv.
    return PiByEightInv::gateAvailable(*context_, target(), controls()) &&
           PiByEightInv::gateCost(*context_, numControls()) < cost();
  }

  if (abs(angle() - 3 * pi / 4) < angleTolerance)  // Angle close to 3pi/4.
  {
    // Gate is equivalent to a PhaseGate and PiByEight combined (or a ZGate and a PiByEightInv). As this is simpler than
    // similar but disregarded cases for the ZRotation, we decide to include these possibilities.
    if (PhaseGate::gateAvailable(*context_, target(), controls()) &&
        PiByEight::gateAvailable(*context_, target(), controls()) &&
        PhaseGate::gateCost(*context_, numControls()) + PiByEight::gateCost(*context_, numControls()) < cost())
    {
      // PhaseGate-PiByEight available and cheaper than original ArbitraryPhase.
      return true;
    }
    if (ZGate::gateAvailable(*context_, target(), controls()) &&
        PiByEightInv::gateAvailable(*context_, target(), controls()) &&
        ZGate::gateCost(*context_, numControls()) + PiByEightInv::gateCost(*context_, numControls()) < cost())
    {
      // ZGate-PiByEightInv available and cheaper than original ArbitraryPhase.
      return true;
    }
    return false;  // Neither option available and cheap enough.
  }

  if (abs(angle() + 3 * pi / 4) < angleTolerance)  // Angle close to -3pi/4.
  {
    // Gate is equivalent to a PhaseInv and PiByEightInv combined (or ZGate and a PiByEight). As this is simpler than
    // similar but disregarded cases for the ZRotation, we decide to include these possibilities.
    if (PhaseInv::gateAvailable(*context_, target(), controls()) &&
        PiByEightInv::gateAvailable(*context_, target(), controls()) &&
        PhaseInv::gateCost(*context_, numControls()) + PiByEightInv::gateCost(*context_, numControls()) < cost())
    {
      // PhaseGate-PiByEight available and cheaper than original ArbitraryPhase.
      return true;
    }
    if (ZGate::gateAvailable(*context_, target(), controls()) &&
        PiByEight::gateAvailable(*context_, target(), controls()) &&
        ZGate::gateCost(*context_, numControls()) + PiByEight::gateCost(*context_, numControls()) < cost())
    {
      // ZGate-PiByEightInv available and cheaper than original ArbitraryPhase.
      return true;
    }
    return false;  // Neither option available and cheap enough. (Unnecessary, but ensures that we only get to final
  }            // return if the angle is not special - nice and clean - and it gives us somewhere to hang this comment.)

  return false;  // Not a special angle - no reduction.
}


GateSequence ArbitraryPhase::reduction() const
{
  // This should only be called if a reduction has been found to exist. That allows for some simplifications.
  // Recall that we aim to replace the ArbitraryPhase with specific gates, e.g. ZGate, not with 'equivalents', since the
  // only gate currently 'equivalent' to the ZGate is the ArbitraryPhase we are trying to replace!)
  double angleTolerance = context_->angleTolerance();

  if (abs(angle()) < angleTolerance)  // Angle close to zero.
  {
    // Gate can be eliminated.
    return {};
  }

  if (abs(abs(angle()) - pi) < angleTolerance)  // Angle close to +-pi.
  {
    // Gate can be replaced with a ZGate.
    // (In future, we could consider options like two PhaseGates if ZGate is unavailable??)
    return makeGateSequence(make_unique<ZGate>(*context_, target(), controls()));
  }

  if (abs(angle() - pi / 2) < angleTolerance)  // Angle close to pi/2.
  {
    // Gate can be replaced by a PhaseGate.
    return makeGateSequence(make_unique<PhaseGate>(*context_, target(), controls()));
  }

  if (abs(angle() + pi / 2) < angleTolerance)  // Angle close to -pi/2.
  {
    // Gate can be replaced by a PhaseInv.
    return makeGateSequence(make_unique<PhaseInv>(*context_, target(), controls()));
  }

  if (abs(angle() - pi / 4) < angleTolerance)  // Angle close to pi/4.
  {
    // Gate can be replaced by a PiByEight.
    return makeGateSequence(make_unique<PiByEight>(*context_, target(), controls()));
  }

  if (abs(angle() + pi / 4) < angleTolerance)  // Angle close to -pi/4.
  {
    // Gate can be replaced by a PiByEightInv.
    return makeGateSequence(make_unique<PiByEightInv>(*context_, target(), controls()));
  }

  if (abs(angle() - 3 * pi /4) < angleTolerance)  // Angle close to 3pi/4.
  {
    // Gate can be replaced by a PhaseGate and a PiByEight, or by a ZGate and a PiByEightInv. Of course, we choose the
    // cheaper of the two options. If both have the same case, we choose the PhaseGate-PiByEight combination (for no
    // particular reason).
    if (PhaseGate::gateAvailable(*context_, target(), controls()) &&
        PiByEight::gateAvailable(*context_, target(), controls()))
    {
      if (ZGate::gateAvailable(*context_, target(), controls()) &&
          PiByEightInv::gateAvailable(*context_, target(), controls()) &&
          ZGate::gateCost(*context_, numControls()) + PiByEightInv::gateCost(*context_, numControls()) <
                  PhaseGate::gateCost(*context_, numControls()) + PiByEight::gateCost(*context_, numControls()))
      {
        // ZGate-PiByEightInv option avaialble and cheaper.
        return makeGateSequence(make_unique<ZGate>(*context_, target(), controls()),
                                make_unique<PiByEightInv>(*context_, target(), controls()));
      }
      // PhaseGate-PiByEight option available. ZGate-PiByEightInv either unavailable or not cheaper.
      return makeGateSequence(make_unique<PhaseGate>(*context_, target(), controls()),
                              make_unique<PiByEight>(*context_, target(), controls()));
    }
    // PhaseGate-PiByEight option unavailable, so ZGate-PiByEightInv must be available.
    return makeGateSequence(make_unique<ZGate>(*context_, target(), controls()),
                            make_unique<PiByEightInv>(*context_, target(), controls()));
  }

  assert(abs(angle() + 3 * pi / 4) < angleTolerance);  // Angle close to -3pi/4.
  {           // Added braces mainly just to match the format of the previous 'if's!
    // Gate can be replaced by a PhaseInv and a PiByEightInv, or by a ZGate and a PiByEight. Of course, we choose the
    // cheaper of the two options. If both have the same case, we choose the PhaseInv-PiByEightInv combination (for no
    // particular reason).
    if (PhaseInv::gateAvailable(*context_, target(), controls()) &&
        PiByEightInv::gateAvailable(*context_, target(), controls()))
    {
      if (ZGate::gateAvailable(*context_, target(), controls()) &&
          PiByEight::gateAvailable(*context_, target(), controls()) &&
          ZGate::gateCost(*context_, numControls()) + PiByEight::gateCost(*context_, numControls()) <
          PhaseInv::gateCost(*context_, numControls()) + PiByEightInv::gateCost(*context_, numControls()))
      {
        // ZGate-PiByEightInv option avaialble and cheaper.
        return makeGateSequence(make_unique<ZGate>(*context_, target(), controls()),
                                make_unique<PiByEight>(*context_, target(), controls()));
      }
      // PhaseGate-PiByEight option available. ZGate-PiByEightInv either unavailable or not cheaper.
      return makeGateSequence(make_unique<PhaseInv>(*context_, target(), controls()),
                              make_unique<PiByEightInv>(*context_, target(), controls()));
    }
    // PhaseGate-PiByEight option unavailable, so ZGate-PiByEightInv must be available.
    return makeGateSequence(make_unique<ZGate>(*context_, target(), controls()),
                            make_unique<PiByEight>(*context_, target(), controls()));
  }
}


long ArbitraryPhase::canSimplify(const Gate& next) const
{
  return next.canSimplify(*this);  // Double dispatch.
}


GateSequence ArbitraryPhase::simplification(const Gate& next) const
{
  return next.simplification(*this);  // Double dispatch.
}


bool ArbitraryPhase::canRSwap(const Gate& next) const
{
  return next.canRSwap(*this);  // Double dispatch.
}


unique_ptr<Gate> ArbitraryPhase::rightRMoverChange(const Gate& next) const
{
  return next.rightRMoverChange(*this);  // Double dispatch.
}


unique_ptr<Gate> ArbitraryPhase::leftRMoverChange(const Gate& next) const
{
  return next.leftRMoverChange(*this);  // Double dispatch.
}


long ArbitraryPhase::canSimplify(const PiByEight& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.canSimplify(*this);
}


GateSequence ArbitraryPhase::simplification(const PiByEight& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.simplification(*this);
}


long ArbitraryPhase::canSimplify(const PiByEightInv& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.canSimplify(*this);
}


GateSequence ArbitraryPhase::simplification(const PiByEightInv& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.simplification(*this);
}


long ArbitraryPhase::canSimplify(const PhaseGate& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.canSimplify(*this);
}


GateSequence ArbitraryPhase::simplification(const PhaseGate& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.simplification(*this);
}


long ArbitraryPhase::canSimplify(const PhaseInv& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.canSimplify(*this);
}


GateSequence ArbitraryPhase::simplification(const PhaseInv& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.simplification(*this);
}


bool ArbitraryPhase::canRSwap(const XGate& prev) const
{
  return prev.canRSwap(*this);
}


unique_ptr<Gate> ArbitraryPhase::rightRMoverChange(const XGate& prev) const
{
  return prev.leftRMoverChange(*this);
}


unique_ptr<Gate> ArbitraryPhase::leftRMoverChange(const XGate& prev) const
{
  return prev.rightRMoverChange(*this);
}


bool ArbitraryPhase::canRSwap(const YGate& prev) const
{
  return prev.canRSwap(*this);
}


unique_ptr<Gate> ArbitraryPhase::rightRMoverChange(const YGate& prev) const
{
  return prev.leftRMoverChange(*this);
}


unique_ptr<Gate> ArbitraryPhase::leftRMoverChange(const YGate& prev) const
{
  return prev.rightRMoverChange(*this);
}


long ArbitraryPhase::canSimplify(const ZGate& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.canSimplify(*this);
}


GateSequence ArbitraryPhase::simplification(const ZGate& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.simplification(*this);
}


long ArbitraryPhase::canSimplify(const ZRotation& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.canSimplify(*this);
}


GateSequence ArbitraryPhase::simplification(const ZRotation& prev) const
{
  // Gate order is unimportant - they commute. To avoid code duplication, call the function with the gates switched.
  return prev.simplification(*this);
}


long ArbitraryPhase::canSimplify(const ArbitraryPhase& prev) const
{
  // Provided the set of all bits involved in one gate match those involved in the other, two ArbitraryPhases give a new
  // ArbitraryPhase with angle equal to the sum of the original angles.
  if (target() == prev.target() && controls() == prev.controls())
  {
    return cost();
  }
  return 0;
}


GateSequence ArbitraryPhase::simplification(const ArbitraryPhase& prev) const
{
  return makeGateSequence(make_unique<ArbitraryPhase>(*context_, target(), controls(), angle() + prev.angle()));
//  return {make_unique<ArbitraryPhase>(*context_, target(), controls(), angle() + prev.angle())};
}


const GateSimulator& ArbitraryPhase::part_simulator(int partIndex) const
{
  switch (partIndex)
  {
    case 0:
      return const_sim;
    case 1:
      return cos_sim;
    case 2:
      return sin_sim;
    default:
      throw std::invalid_argument("Invalid partIndex passed to XRotation::part_simulator().");
  }
}


double ArbitraryPhase::thetaToAngle(double theta) const
{
  return theta;
}


vector<vector<cmplx> > ArbitraryPhase::matrix_() const
{
  using std::exp, constants::i;
  double angle = this->angle_;  // Yes, the 'this' is necessary! (Apparently not on Visual Studio, but that's a...
  return {{1, 0}, {0, exp(i * angle)}};  // ...'non-conformance'.)
}

//----------------------------------------------------------------------------------------------------------------------

SwapGate::SwapGate(const CircuitContext& context) :
Gate(context),
swap_bit1(0),
swap_bit2(1)  // We ensure that the gate is valid, despite being 'default' constructed.
{
}


SwapGate::SwapGate(const CircuitContext& context, int bit1, int bit2) :
Gate(context),
swap_bit1(bit1),
swap_bit2(bit2)
{
  assert(0 <= bit1 && bit1 < context.numQbits());
  assert(0 <= bit2 && bit2 < context.numQbits());
  assert(bit1 != bit2);

  fix_order();  // Ensure that swap_bit1 < swap_bit2.
}


int SwapGate::gateTypeId() const
{
  return gate_type_id;
}


int SwapGate::numQbitOptions() const
{
  // This is only called from the CircuitContext constructor, on the default SwapGate created there, i.e. one that swaps
  // bits 0 and 1. If the permitted controls selected by the user are something like 'one' or 'atLeastOne', then this
  // SwapGate and indeed any other is unavailable. In this case we throw a useful exception.
  if (!available())
  {
    throw std::logic_error("In SwapGate::numQbitOptions(), we find that SwapGate is 'available' but not without"
                           " controls, which is the only option.");
  }
  return num_qbits() * (num_qbits() - 1) / 2;
}


unique_ptr<Gate> SwapGate::clone() const
{
  return make_unique<SwapGate>(*this);  // Will get converted to unique_ptr<Gate>.
}


SwapGate* SwapGate::rawClone() const
{
  return new SwapGate(*this);
}


void SwapGate::random()
{
  // Inefficient, but clear.
  assert(num_qbits() > 1);
  swap_bit1 = randInt(0, num_qbits());
  do
  {
    swap_bit2 = randInt(0, num_qbits());
  }
  while (swap_bit1 == swap_bit2);

  fix_order();  // Ensure that bit1 < bit2.
}


void SwapGate::swapBits(int i, int j)
{
  swap_bit1 = swapValues(swap_bit1, i, j);
  swap_bit2 = swapValues(swap_bit2, i, j);
  fix_order();
}


int SwapGate::bit1() const
{
  return swap_bit1;
}


int SwapGate::bit2() const
{
  return swap_bit2;
}


bool SwapGate::available() const
{
  return context_->gateTypeAvailable(*this, 0);
}


void SwapGate::applyTo(State& state) const
{
  state.swapQbits(swap_bit1, swap_bit2);
}


void SwapGate::applyInvTo(State& state) const
{
  applyTo(state);  // Swap gate is its own inverse.
}


int SwapGate::target() const
{
  throw std::logic_error("Function SwapGate::target() called, but SwapGates do not have a single target.");
}


std::string SwapGate::name() const
{
  return "Swap";
}


bool SwapGate::gateAvailable(const CircuitContext& context, int bit1, int bit2)  // Static
{
  // Static function for determining the availability of a Gate without needing to create one. (The corresponding
  // non-static function is useful when we have a Gate but don't necessarily know the precise gate type.)
  // (We don't need to provide the target bits yet. However, in future we may have more complicated rules on
  // availability.)
  return context.gateTypeAvailable(gate_type_id, 0);
}


long SwapGate::gateCost(const CircuitContext& context)  // Static
{
  return context.gateCost(gate_type_id, 0);
}


void SwapGate::output(std::ostream& out) const
{
  out << "Swap on target qbits " << swap_bit1 << ", " << swap_bit2 << endl;
}


bool SwapGate::cancelsAtStart() const
{
  // A swap gate at the beginning of the circuit does nothing if both inputs are always |0>, or if both inputs are
  // always |1>.
  auto bit1Options = context_->qbitInputOptions(bit1());
  auto bit2Options = context_->qbitInputOptions(bit2());
  return bit1Options == bit2Options && bit1Options != QbitOptions::varies;
}


long SwapGate::canSimplify(const Gate& next) const
{
  return next.canSimplify(*this);  // Double dispatch.
}


GateSequence SwapGate::simplification(const Gate& next) const
{
  return next.simplification(*this);  // Double dispatch.
}


bool SwapGate::canSwapSwap(const Gate& next) const
{
  return next.canSwapSwap(*this);  // Double dispatch.
}


unique_ptr<Gate> SwapGate::rightSwapMoverChange(const Gate &next) const
{
  return next.rightSwapMoverChange(*this);  // Double dispatch.
}


unique_ptr<Gate> SwapGate::leftSwapMoverChange(const Gate &next) const
{
  return next.leftSwapMoverChange(*this);  // Double dispatch.
}


bool SwapGate::canSwapSwap(const SingleTargetGate& prev) const
{
  // Can swap the gates, provided we swap the roles of the qbits on the SingleTargetGate appropriately.
  return true;
}


unique_ptr<Gate> SwapGate::rightSwapMoverChange(const SingleTargetGate &prev) const
{
  // It doesn't matter in which direction the SingleTargetGate moves past the SwapGate, it's bits get swapped just the
  // same.
  return prev.leftSwapMoverChange(*this);
}


unique_ptr<Gate> SwapGate::leftSwapMoverChange(const SingleTargetGate &prev) const
{
  return {};  // Swap gate is unchanged.
}


long SwapGate::canSimplify(const SwapGate& prev) const
{
  // Two identical SwapGates, i.e. gates with the same qbits, cancel.
  if (*this == prev)
  {
    return cost() + prev.cost();
  }
  return 0;
}


GateSequence SwapGate::simplification(const SwapGate& prev) const
{
  return {};  // Gates cancel.
}


bool SwapGate::canSwapSwap(const SwapGate& prev) const
{
  // Can swap the gates, provided we swap the roles of the qbits on one of the SwapGates appropriately.
  return true;
}


unique_ptr<Gate> SwapGate::rightSwapMoverChange(const SwapGate &prev) const
{
  // Here we consider only the one gate swap that leads to the 'least' result, i.e. the pair of swap gates that would be
  // 'sorted before' the other possible pair. We can do this as, currently, the code only uses that result. (Note that
  // we could extend this to canSwapSwap() - that is, we could say that a swap is not possible if, in fact, in can only
  // result in a 'increase' according to sortBefore().

  // We are looking for those scenarios where changing the left mover (i.e. *this) would result in an increase,
  // according to sortOrder(). This happens when the first qbits are equal and the others are not - changing the left
  // mover would result in the first qbit being increased to the value of prev's second qbit. It also happens when the
  // second qbit of the left mover equals the first qbit of the right mover - changing the left mover would result in
  // its second qbit increasing to that of the right mover.
  if (swap_bit1 == prev.swap_bit1)
  {
    if (swap_bit2 != prev.swap_bit2)
    {
      // If we change the left mover (i.e. *this), then it's swap_bit1 will increase to prev.swap_bit2, i.e. the gate
      // increases according to sort order. This is not what we want, so we keep left mover fixed and change right
      // mover.
      return make_unique<SwapGate>(*context_, swap_bit2, prev.swap_bit2);  // Constructor will sort the qbits.
    }
  }
  else if (swap_bit2 == prev.swap_bit1)
  {
    // If we change the left mover (i.e. *this), then it's swap_bit2 will increase to prev.swap_bit2. Hence we change
    // the right mover.
    return make_unique<SwapGate>(*context_, swap_bit1, prev.swap_bit2);
  }
  return {};
}


unique_ptr<Gate> SwapGate::leftSwapMoverChange(const SwapGate &prev) const
{
  // Here we consider only the one gate swap that leads to the 'least' result, i.e. the pair of swap gates that would be
  // 'sorted before' the other possible pair. We can do this as, currently, the code only uses that result. Note that we
  // could extend this to canSwapSwap() - that is, we could say that a swap is not possible if, in fact, in can only
  // result in a 'increase' according to sortBefore().

  // The scenarios we seek here those where the left mover changes, i.e. those not covered in rightSwapMoverChange() but
  // where the two gates share one qbit.
  if (swap_bit2 == prev.swap_bit2)
  {
    if (swap_bit1 != prev.swap_bit1)
    {
      // Changing the left mover (i.e. *this) results in swap_bit2 decreasing to prev.swap_bit1, which is good.
      return make_unique<SwapGate>(*context_, swap_bit1, prev.swap_bit1);  // Constructor will sort the qbits.
    }
  }
  else if (swap_bit1 == prev.swap_bit2)
  {
    // Changing the left mover (i.e. *this) results in swap_bit1 decreasing to prev.swap_bit1, which is good.
    return make_unique<SwapGate>(*context_, prev.swap_bit1, swap_bit2);
  }
  return {};
}


bool SwapGate::equivalent_structure(const Gate& rhs) const
{
  // If the code gets here, then rhs must also be of SwapGate type and the gates should have the same 'context', i.e.
  // costs, numQbits, etc.
  const SwapGate& other = static_cast<const SwapGate&>(rhs);

  assert(swap_bit1 < swap_bit2);
  assert(other.swap_bit1 < other.swap_bit2);

  return swap_bit1 == other.swap_bit1 && swap_bit2 == other.swap_bit2;
}



bool SwapGate::equals_(const Gate& rhs) const
{
  return equivalent_structure(rhs);
}


bool SwapGate::sort_before(const Gate& rhs) const
{
  // If the code gets here, then rhs must also be of SwapGate type and both gates should have the same 'context', i.e.
  // num_qbits, cost etc.
  const SwapGate& other = static_cast<const SwapGate&>(rhs);

  assert(swap_bit1 < swap_bit2);
  assert(other.swap_bit1 < other.swap_bit2);

  if (swap_bit1 != other.swap_bit1)
  {
    return swap_bit1 < other.swap_bit1;
  }
  return swap_bit2 < other.swap_bit2;
}


void SwapGate::fix_order()
{
  assert(swap_bit1 != swap_bit2);
  if (swap_bit2 < swap_bit1)
  {
    std::swap(swap_bit1, swap_bit2);
  }
}

//----------------------------------------------------------------------------------------------------------------------

Oracle::Oracle(const CircuitContext& context) :
Gate(context)
{
}


void Oracle::linkToMarkedState(std::shared_ptr<const int> markedState)
{
  // At present, the Oracle is the only gate that requires access to 'circuit data', namely the identity of the 'marked
  // state'.
  marked_state = markedState;
}


int Oracle::gateTypeId() const
{
  return gate_type_id;
}


int Oracle::numQbitOptions() const
{
  // This is only called from the CircuitContext constructor, on the default Oracle created there. If the permitted
  // controls selected by the user are something like 'one' or 'atLeastOne', then this Oracle and any other is
  // unavailable, so we throw a useful exception.
  if (!available())
  {
    throw std::logic_error("In Oracle::numQbitOptions(), we find that Oracle is 'available' but not without controls,"
                           " which is the only option.");
  }
  return 1;
}


unique_ptr<Gate> Oracle::clone() const
{
  return make_unique<Oracle>(*this);  // Will get converted to unique_ptr<Gate>.
}


Oracle* Oracle::rawClone() const
{
  return new Oracle(*this);
}


void Oracle::random()
{
  // All oracle gates are the same!
}


void Oracle::swapBits(int i, int j)
{
  cerr << "Function swapBits() called for an Oracle, for which this function is invalid." << endl;
  assert(false);
}


bool Oracle::available() const
{
  return context_->gateTypeAvailable(*this, 0);
}


void Oracle::applyTo(State& state) const
{
  state.mark(*marked_state);
}


void Oracle::applyInvTo(State& state) const
{
  applyTo(state);  // Oracle is its own inverse.
}


int Oracle::target() const
{
  throw std::logic_error("Function Oracle::target() called, but an Oracle does not have a single target qbit.");
}


std::string Oracle::name() const
{
  return "Oracle";
}


bool Oracle::gateAvailable(const CircuitContext& context)  // Static
{
  // Static function for determining the availability of a Gate without needing to create one. (The corresponding
  // non-static function is useful when we have a Gate but don't necessarily know the precise gate type.)
  return context.gateTypeAvailable(gate_type_id, 0);
}


long Oracle::gateCost(const CircuitContext& context)  // Static
{
  return context.gateCost(gate_type_id, 0);
}


bool Oracle::cancelsAtStart() const
{
  // We assume for now that, if Oracles are present, the evaluation of the circuit involves simulating the circuit for
  // each possible Oracle setting, i.e. for each value of the special state. An error value is produced for each
  // setting. This error value should not change if all output states (for the chosen Oracle setting) are multiplied by
  // an overall phase, since this phase is unphysical and unmeasurable. The overall error is produced by simply summing
  // (or averaging) the error over the different Oracle settings.
  // If ALL of the inputs to the circuit are fixed, then an Oracle at the start of the circuit either does nothing or,
  // if its special state is set to match the input, multiplies the state by an overall phase of -1. Both of these will
  // have no effect on the error, so the gate is redundant. If just one of the inputs to the circuit varies, an Oracle
  // at the start of the circuit is NOT redundant, since one can choose a special state such that the Oracle may, or may
  // not mark the input, depending on the value of the varying qbit.
  // (Note that an Oracle that occurs second in the circuit, after an XGate, is redundant if the input to the circuit is
  // fixed. At present we do not catch this sort of simplification - we do not trace the 'fixedness' of the state beyond
  // the input to the circuit. (We may have similar simplications available for other gates too, in particular the
  // PhaseTypeGates.))
  for (int qbit = 0; qbit < num_qbits(); ++qbit)
  {
    if (context_->qbitInputOptions(qbit) == QbitOptions::varies)
    {
      return false;
    }
  }
  return true;
}


bool Oracle::canSimplySwap(const Gate& next) const
{
  return next.canSimplySwap(*this);  // Double dispatch.
}


long Oracle::canSimplify(const Gate& next) const
{
  return next.canSimplify(*this);  // Double dispatch.
}


GateSequence Oracle::simplification(const Gate& next) const
{
  return next.simplification(*this);  // Double dispatch.
}


bool Oracle::canSimplySwap(const DiagonalGate& prev) const
{
  return true;
}


long Oracle::canSimplify(const Oracle& prev) const
{
  // Oracles just cancel.
  return cost() + prev.cost();
}


GateSequence Oracle::simplification(const Oracle& prev) const
{
  return {};  // Oracles cancel.
}


bool Oracle::equivalent_structure(const Gate& rhs) const
{
  // If the code gets here, then rhs must also be an Oracle and both gates have the same 'context'. All Oracle gates
  // have the same structure.
  return true;
}


bool Oracle::equals_(const Gate& rhs) const
{
  // If the code gets here, then rhs must also be an Oracle and both gates have the same 'context'. All Oracle gates
  // are the same.
  return true;
}


bool Oracle::sort_before(const Gate& rhs) const
{
  // If the code gets here, then rhs must also be an Oracle and both gates should have the same 'context'. All Oracle
  // gates are the same.
  return false;
}

//----------------------------------------------------------------------------------------------------------------------

// A master list of GateCreators for each (concrete) Gate type, indexed by gate name. The idea is that it will be
// possible for a Problem to read in the names of the Gate types available to it and then pick the types out of this
// list to make its own vector of GateCreators for use by Solutions.
// (Unpleasant.)
const std::map<std::string, GateCreator> masterGateCreatorIndex =
{
  {"Hadamard", &createGate<Hadamard>}, {"PiByEight", &createGate<PiByEight>},
  {"PiByEightInv", &createGate<PiByEightInv>}, {"Phase", &createGate<PhaseGate>}, {"PhaseInv", &createGate<PhaseInv>},
  {"X", &createGate<XGate>}, {"Y", &createGate<YGate>}, {"Z", &createGate<ZGate>}, {"XRot", &createGate<XRotation>},
  {"YRot", &createGate<YRotation>}, {"ZRot", &createGate<ZRotation>}, {"ArbPhase", &createGate<ArbitraryPhase>},
  {"Swap", &createGate<SwapGate>}, {"Oracle", &createGate<Oracle>}
};

//----------------------------------------------------------------------------------------------------------------------

// A helper function for GateSequences.
long gateCost(const GateSequence& gates)
{
  return std::accumulate(gates.begin(), gates.end(), 0,
                         [](int total, const unique_ptr<Gate>& gate) {return total + gate->cost();} );
}
