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

#ifndef GATE_H
#define GATE_H

#include <string>
#include <map>
#include <iostream>
#include <memory>
#include "rng.h"
#include "constants.h"
#include "simulator.h"
#include "controls.h"

namespace circuit
{
  class CircuitContext;
  class Problem;
  class Circuit;
}
using circuit::CircuitContext;
using circuit::Problem;
using circuit::Circuit;

class State;
class Replacement;

// Using double dispatch for simplification, the base class needs to know about the existence of the base classes.
class Gate;
class SingleTargetGate;
class RotationGate;
class DiagonalGate;
class PhaseTypeGate;
class XTypeGate;
class YTypeGate;
class Hadamard;
class PiByEight;
class PiByEightInv;
class PhaseGate;
class PhaseInv;
class XGate;
class YGate;
class ZGate;
class XRotation;
class YRotation;
class ZRotation;
class ArbitraryPhase;
class SwapGate;
class Oracle;

using GateSequence = std::vector<std::unique_ptr<Gate> >;

// Note on double dispatch (needs updating?)
// The routines that allow for gate pair simplification and swapping make extensive use of double-dispatch. While I used
// to use routines in pairs, with a function called somthing like 'simplifies(const Gate& rhs)' calling something like
// 'rhs.isSimplifiedBy(*this)', at present, to simplify naming, we simply have 'canSimplify(const Gate& next)' calling
// 'next.canSimplify(*this)'. We note below that this creates some dangers.
//
// We want functions of the form SomeGate::swapsTo(const Gate& next) to be public, as these are the 1st dispatch. Some
// of the remaining functions can easily be hidden. For example, SingleTargetGate::swapsTo(const SwapGate& prev) can
// simply be marked as private, as it is an override that is never explicitly called. (We call Gate::swapsTo(const
// SwapGate& prev) instead, which 'redirects'.) However, functions of the form Gate::swapsTo(const DerivedGate& prev)
// are trickier. Making them protected doesn't work, as then DerivedGate::swapsTo(const Gate& next) (1st dispatch) is
// denied access. (DerivedGate gets access to protected members of its own base Gate object, but not to those in
// 'next'.) We cannot add 'friend someType DerivedGate::swapsTo(const Gate& next) const' to class Gate, since at that
// point the DerivedGate is still an incomplete type. I am yet to find a solution.
//
// As for the 'issues' that we are trying to avoid, they seem to be of two types. The first is illustrated by swapsTo().
// The function that is called in Circuit::simplify() determines what the current gate (*this) and the next become upon
// applying the 'swap'. However, it calls functions of the same name that determine what the current gate and the
// PREVIOUS one become. I have endeavoured to reduce possible confusion in this case by naming the parameter 'prev' or
// 'next' depending on how the particular function should be called. Then, if for some unforeseen reason, I call one of
// the functions meant for 2nd dispatch directly, then provided I checked that parameter, I should get the correct
// results. The second issue is illustrated best by canSwap(). Consider what happens if I call Gate::canSwap(const
// SingleTargetGate& other) directly. Provided 'other' is just a SingleTargetGate, then there is no problem. However, if
// 'other' is a DiagonalGate - a class deriving from SingleTargetGate, then I have a problem. If *this is also a
// DiagonalGate, then I wish to call DiagonalGate::canSwap(const DiagonalGate& other) which always returns true.
// However, DiagonalGate::canSwap(const SingleTargetGate& other) is called instead, potentially producing the wrong
// result.
//
// (I suspect that this latter problem could also be an issue when the functions come in pairs, since one could call
// Gate::isSimplifiedBy(const SingleTargetGate& lhs) when lhs is a DiagonalGate. This is not a problem at present, since
// only the leaf classes override the 'simplification' routines.)
//
// At present, all such functions are public, including those I could easily make private.


int totalGateCost(const GateSequence& gateSequence);


class Gate
{
  // (This base Gate class is slowly filling up with functions such as mutate() and applyCosinePartTo(), and even
  // target(), which only apply to derived types. (Rotation gates for the first two, and SingleTargetGates for
  // target().) Is there a neat and tidy way of avoiding this class pollution?)
public:
  explicit Gate(const CircuitContext& context);
  virtual ~Gate();

  // Function that creates and attaches simulators, if the gate needs them and doesn't yet have them.
  virtual void createSimulators();

  // Function that ensures that the gate gets access to the correct 'circuit data'. (Think of this as the gate getting
  // to know its 'parent circuit.) The only circuit data at present is the state marked by the Oracle gate.
  virtual void linkToMarkedState(std::shared_ptr<const int> markedState);

  // (Neither of the next two functions actually refer to the gate referenced by 'this'. I would prefer not to have
  // 'if type == ' style code, which is what having a gateTypeId enables, but I haven't worked out a better way yet. The
  // function numQbitOptions() might belong better in the Context object - certainly if I ever create individual Context
  // classes for each gate type, then this would be moved.)
  virtual int gateTypeId() const = 0;  // A number (from 0 to 16) indicating gate type.
  virtual int numQbitOptions() const = 0;  // Number of options for target/control bits, given the gate's 'context'.

  virtual std::unique_ptr<Gate> clone() const = 0;  // Virtual copy constructor.
  virtual Gate* rawClone() const = 0;  // Use only where unique or shared_ptrs are unsuitable, likely due to efficiency.

  virtual void random() = 0;
  virtual void mutate();  // Works only on parameterized gates.
  virtual void swapBits(int i, int j) = 0;  // Should not be called on Oracle gates!

  virtual int numParameters() const;  // Number of angles associated with the gate. Default is zero.
  virtual double parameter(int paramNum) const;           // (Could better code structure avoid the need for these...
  virtual void setParameter(int paramNum, double value);  // ...rather clunky looking accessor functions.)
  virtual void randomizeParameters();

  // Functions for getting the availability and costs of the different gate types. These will access data stored in the
  // CircuitContext.
  virtual bool available() const = 0;
  virtual long cost() const;

  // Operator== first checks that the types of gate are the same. If so, calls equals() polymorphically.
  bool operator==(const Gate& rhs) const;
  bool operator!=(const Gate& rhs) const;

  // sortBefore() is only used by the MO framework, as far as I can tell, to calculate the entropy of a population.
  bool sortBefore(const Gate& rhs) const;

  virtual void applyTo(State& state) const = 0;
  virtual void applyInvTo(State& state) const = 0;

  // Version of 'applyTo' used for symbolic application of a gate to a state, i.e. leaving any rotation angle
  // unspecified. Each part of the matrix, e.g. the constant part, the cosine part and the sine part, are applied
  // separately - hence the name.
  virtual void applyPartTo(State& state, int partIndex) const;  // Should be called ONLY if the gate is mutatable.

  // Function to optimize gate angle, given the symbolic form of the error
  virtual void optimizeAngle(const std::vector<std::vector<double> >& errorCoeffs);  // Should be called ONLY if the...
                                                                                     // ...gate is mutatable.
  virtual int target() const = 0;  // Should be called ONLY when the gate has a single target.
  virtual bool mutatable() const;  // Alas, can't make this static and virtual.

  virtual std::string name() const = 0;
  virtual void output(std::ostream& out) const;
  virtual void extendedOutput(std::ostream& out) const;  // Available for testing purposes.


  // The remaining public member functions are for circuit simplification and getting the circuit into canonical form,
  // indicating whether a pair of gates can be simplified or swapped, and the outcome of performing such a
  // simplification or swap. The routines for swapping gates are grouped into three sets. Simple swaps are (most of)
  // those that require no change to the gates. HSwaps are those that involve a Hadamard that shares its target with the
  // other (X, Y, Z) gate. The non-Hadamard will usually change, except for the YGate. RSwaps involve a rotation gate
  // and an X, Y or Z. The rotation gate will usually change. SwapSwaps are swaps involving SwapGates, which may (or may
  // not) require a change to a gate.
  // We also include in this section, routines for gate 'reduction' - i.e. routines that take a parameterized gate and
  // determine whether it can be removed or replaced by a cheaper non-parameterized gate. These are useful after
  // numerical optimization, when gates are likely to be reducible.

  virtual bool canReduce() const;
  virtual GateSequence reduction() const;
  // (There is also an argument for having reduction functions that consider neighbouring gates. In particular, the
  // reduction of a controlled XRotation of angle pi to a controlled XGate and a PhaseGate on the controls is unlikely
  // to be taken at present, since the new gates would likely cost more than the old. However, we could detect the
  // presence of, for example, a neighbouring PhaseInv that would cancel the introduced PhaseGate.)

  virtual bool cancelsAtStart() const = 0;
  virtual long canSimplify(const Gate& next) const = 0;  // Returns cost improvement.
  virtual GateSequence simplification(const Gate& next) const = 0;
  virtual bool canSimplySwap(const Gate& next) const;
  virtual bool canSwapSwap(const Gate& next) const;  // All swaps involving SwapGates.
  virtual std::unique_ptr<Gate> rightSwapMoverChange(const Gate& next) const;  // Returns nullptr if there is no change.
  virtual std::unique_ptr<Gate> leftSwapMoverChange(const Gate& next) const;  // Returns nullptr if there is no change.
  virtual bool canHSwap(const Gate& next) const;  // Swaps involving Hadamards sharing targets with X, Y, Z.
  virtual std::unique_ptr<Gate> rightHMoverChange(const Gate& next) const;  // Returns nullptr if there is no change.
  virtual std::unique_ptr<Gate> leftHMoverChange(const Gate& next) const;  // Returns nullptr if there is no change.
  virtual bool canRSwap(const Gate& next) const;  // Swaps involving rotation gates with X, Y, Z, S, S-1.
  virtual std::unique_ptr<Gate> rightRMoverChange(const Gate& next) const;  // Returns nullptr if there is no change.
  virtual std::unique_ptr<Gate> leftRMoverChange(const Gate& next) const;  // Returns nullptr if there is no change.

  // The above simplification routines are those that should be called from outside the gate classes. However, we use
  // 'double-dispatch'. Unfortunately, this requires extra virtual functions whenever we add a new derived class, or at
  // least any derived gate class that permits simplification.

  virtual bool canSimplySwap(const SingleTargetGate& prev) const;
  virtual bool canSwapSwap(const SingleTargetGate& prev) const;  // Swaps involving SwapGates.
  virtual std::unique_ptr<Gate> rightSwapMoverChange(const SingleTargetGate& prev) const;  // These return nullptr if...
  virtual std::unique_ptr<Gate> leftSwapMoverChange(const SingleTargetGate& prev) const;  // ...there is no change.

  virtual bool canSimplySwap(const DiagonalGate& prev) const;

  virtual long canSimplify(const Hadamard& prev) const;  // Returns cost improvement.
  virtual GateSequence simplification(const Hadamard& prev) const;
  virtual bool canHSwap(const Hadamard& prev) const;
  virtual std::unique_ptr<Gate> rightHMoverChange(const Hadamard& prev) const; // Returns nullptr if there is no change.
  virtual std::unique_ptr<Gate> leftHMoverChange(const Hadamard& prev) const;  // Returns nullptr if there is no change.

  virtual long canSimplify(const PiByEight& prev) const;  // Returns cost improvement.
  virtual GateSequence simplification(const PiByEight& prev) const;

  virtual long canSimplify(const PiByEightInv& prev) const;  // Returns cost improvement.
  virtual GateSequence simplification(const PiByEightInv& prev) const;

  virtual long canSimplify(const PhaseGate& prev) const;
  virtual GateSequence simplification(const PhaseGate& prev) const;
  virtual bool canRSwap(const PhaseGate& prev) const;
  virtual std::unique_ptr<Gate> rightRMoverChange(const PhaseGate& prev) const;
  virtual std::unique_ptr<Gate> leftRMoverChange(const PhaseGate& prev) const;

  virtual long canSimplify(const PhaseInv& prev) const;
  virtual GateSequence simplification(const PhaseInv& prev) const;
  virtual bool canRSwap(const PhaseInv& prev) const;
  virtual std::unique_ptr<Gate> rightRMoverChange(const PhaseInv& prev) const;
  virtual std::unique_ptr<Gate> leftRMoverChange(const PhaseInv& prev) const;

  virtual long canSimplify(const XGate& prev) const;
  virtual GateSequence simplification(const XGate& prev) const;
  virtual bool canHSwap(const XGate& prev) const;
  virtual std::unique_ptr<Gate> rightHMoverChange(const XGate& prev) const;
  virtual std::unique_ptr<Gate> leftHMoverChange(const XGate& prev) const;
  virtual bool canRSwap(const XGate& prev) const;
  virtual std::unique_ptr<Gate> rightRMoverChange(const XGate& prev) const;
  virtual std::unique_ptr<Gate> leftRMoverChange(const XGate& prev) const;

  virtual long canSimplify(const YGate& prev) const;
  virtual GateSequence simplification(const YGate& prev) const;
  virtual bool canHSwap(const YGate& prev) const;
  virtual std::unique_ptr<Gate> rightHMoverChange(const YGate& prev) const;  // These always return nullptr which...
  virtual std::unique_ptr<Gate> leftHMoverChange(const YGate& prev) const;   // ...indicates no change.
  virtual bool canRSwap(const YGate& prev) const;
  virtual std::unique_ptr<Gate> rightRMoverChange(const YGate& prev) const;  // These always return nullptr which...
  virtual std::unique_ptr<Gate> leftRMoverChange(const YGate& prev) const;   // ...indicates no change.

  virtual long canSimplify(const ZGate& prev) const;
  virtual GateSequence simplification(const ZGate& prev) const;
  virtual bool canHSwap(const ZGate& prev) const;
  virtual std::unique_ptr<Gate> rightHMoverChange(const ZGate& prev) const;
  virtual std::unique_ptr<Gate> leftHMoverChange(const ZGate& prev) const;
  virtual bool canRSwap(const ZGate& prev) const;
  virtual std::unique_ptr<Gate> rightRMoverChange(const ZGate& prev) const;
  virtual std::unique_ptr<Gate> leftRMoverChange(const ZGate& prev) const;

  virtual long canSimplify(const XRotation& prev) const;
  virtual GateSequence simplification(const XRotation& prev) const;
  virtual bool canHSwap(const XRotation& prev) const;
  virtual std::unique_ptr<Gate> rightHMoverChange(const XRotation& prev) const;
  virtual std::unique_ptr<Gate> leftHMoverChange(const XRotation& prev) const;
  virtual bool canRSwap(const XRotation& prev) const;
  virtual std::unique_ptr<Gate> rightRMoverChange(const XRotation& prev) const;
  virtual std::unique_ptr<Gate> leftRMoverChange(const XRotation& prev) const;

  virtual long canSimplify(const YRotation& prev) const;
  virtual GateSequence simplification(const YRotation& prev) const;
  virtual bool canHSwap(const YRotation& prev) const;
  virtual std::unique_ptr<Gate> rightHMoverChange(const YRotation& prev) const;
  virtual std::unique_ptr<Gate> leftHMoverChange(const YRotation& prev) const;
  virtual bool canRSwap(const YRotation& prev) const;
  virtual std::unique_ptr<Gate> rightRMoverChange(const YRotation& prev) const;
  virtual std::unique_ptr<Gate> leftRMoverChange(const YRotation& prev) const;

  virtual long canSimplify(const ZRotation& prev) const;
  virtual GateSequence simplification(const ZRotation& prev) const;
  virtual bool canHSwap(const ZRotation& prev) const;
  virtual std::unique_ptr<Gate> rightHMoverChange(const ZRotation& prev) const;
  virtual std::unique_ptr<Gate> leftHMoverChange(const ZRotation& prev) const;
  virtual bool canRSwap(const ZRotation& prev) const;
  virtual std::unique_ptr<Gate> rightRMoverChange(const ZRotation& prev) const;
  virtual std::unique_ptr<Gate> leftRMoverChange(const ZRotation& prev) const;

  virtual long canSimplify(const ArbitraryPhase& prev) const;
  virtual GateSequence simplification(const ArbitraryPhase& prev) const;
  virtual bool canRSwap(const ArbitraryPhase& prev) const;
  virtual std::unique_ptr<Gate> rightRMoverChange(const ArbitraryPhase& prev) const;
  virtual std::unique_ptr<Gate> leftRMoverChange(const ArbitraryPhase& prev) const;

  virtual long canSimplify(const SwapGate& prev) const;
  virtual GateSequence simplification(const SwapGate& prev) const;
  virtual bool canSwapSwap(const SwapGate& prev) const;  // Swaps involving SwapGates.
  virtual std::unique_ptr<Gate> rightSwapMoverChange(const SwapGate& prev) const;
  virtual std::unique_ptr<Gate> leftSwapMoverChange(const SwapGate& prev) const;

  virtual bool canSimplySwap(const Oracle& prev) const;
  virtual long canSimplify(const Oracle& prev) const;  // Returns cost improvement.
  virtual GateSequence simplification(const Oracle& prev) const;

protected:
  int num_qbits() const;

private:
  virtual bool equivalent_structure(const Gate& rhs) const = 0;  // Called ONLY if rhs is of the same type.
  virtual bool equals_(const Gate& rhs) const = 0;  // Called ONLY if rhs is of the same type.
  virtual bool sort_before(const Gate& rhs) const = 0;  // Called ONLY if rhs is of the same type.

protected:
  const CircuitContext* context_;  // (Might it be possible to point to a smaller object - a GateContext for each...
};                                 // ...type of gate?)

std::ostream& operator<<(std::ostream& out, const Gate& gate);


class SingleTargetGate : public Gate
{
  // Class used to represent gates that work on a single target qbit. Such gates are permitted to have multiple control
  // qbits, but only one qbit is changed.
  //
  // This class is designed to be manipulated by the optimization algorithm and hence doesn't hold the transformation
  // matrix directly, but holds shared_ptrs to GateSimulators (for both forward and reverse simulation) that are
  // constructed whenever a new Gate is created. (Copies of the same gate have shared_ptrs to the same simulators.)
  //
  // Where default versions of member functions such as random() are provided, it is assumed that the typical derived
  // class will be something like a Hadamard or an X gate, i.e. a parameter (angle) free gate. For rotation gates, these
  // default implementations will be overridden to account for the need to set or mutate the angle parameter.
public:
  explicit SingleTargetGate(const CircuitContext& context);
  SingleTargetGate(const CircuitContext& context, int target, const Controls& controls);

  // Here I include the default copy operations, in order to disable the default move operations. This is because this
  // class will be a virtual base class and certain gates - ZRotation, ArbitraryPhase - will inherit from it twice. The
  // default move operations would then move from the SingleTargetGate twice - very bad! We might, in future, want to
  // consider putting move operations back in, for efficiency, but this will require care.
  SingleTargetGate(const SingleTargetGate&) = default;
  SingleTargetGate& operator=(const SingleTargetGate&) = default;

  // Function for creating the forward and reverse simulators for the gate.
  virtual void createSimulators() override;

  virtual int numQbitOptions() const override;

  void random() override;

  const Controls& controls() const;  // Added for use by simplification code.
  int numControls() const;

  virtual void swapBits(int i, int j) override;

  bool available() const override;
  long cost() const override;

  std::vector<int> allQbits() const;  // Returns a sorted vector of all qbits used or affected by the gate. Slow!

  // The following three functions are used to avoid the inefficient allQbits().
  // (I should come up with better names. In particular, allQbitsMatch() might give the impression that the targets
  // match and the controls match. This is NOT the case. Instead, the SET of all qbits involved in one gate match the
  // SET of all qbits involved in the other, i.e. allQbits() produces the same results for the two gates.)
  // (Should these be protected rather than public?)
  bool allAreInvolvedIn(const SingleTargetGate& other) const;  // All qbits of 'this' are involved in 'other'. (Subset.)
  bool allQbitsMatch(const SingleTargetGate& other) const;  // Qbits of 'this' are those of 'other'. (Equality.)
  bool allQbitsMatchControlsOf(const SingleTargetGate& other) const;  // Qbits of 'this' are the controls of 'other'.
                                                                      // (Equality.)
  bool controlsAreControlsOf(const SingleTargetGate& other) const; // Controls of 'this' are controls 'other'. (Subset.)

  void applyTo(State& state) const override;
  void applyInvTo(State& state) const override;

  int target() const override;

  virtual bool matricesCommute(const SingleTargetGate& other) const;
  virtual bool matricesCommute(const DiagonalGate& other) const;  // Unnecessary, but makes function name more accurate!
  virtual bool matricesCommute(const XTypeGate& other) const;
  virtual bool matricesCommute(const YTypeGate& other) const;
  virtual bool matricesAnticommute(const SingleTargetGate& other) const;  // Note that if 2x2 matrices commute up to...
  virtual bool matricesAnticommute(const XGate& other) const;             // ...a phase factor, then that factor must...
  virtual bool matricesAnticommute(const YGate& other) const;             // ...be +-1, i.e. the matrices commute or...
  virtual bool matricesAnticommute(const ZGate& other) const;             // ...anticommute.

  void output(std::ostream& out) const override;
  void extendedOutput(std::ostream& out) const override;  // Adds the simulators to the output.


  // Remaining public member functions are overrides for the simplification routines - starting with the 'first
  // dispatch' versions and then the 'second dispatch' versions grouped by gate type.
  virtual bool cancelsAtStart() const override;  // Should be able to refer to CircuitContext instead?
  virtual bool canSimplySwap(const Gate& next) const override;
  virtual bool canSwapSwap(const Gate& next) const override;
  virtual std::unique_ptr<Gate> rightSwapMoverChange(const Gate& next) const override;  // These return nullptr if...
  virtual std::unique_ptr<Gate> leftSwapMoverChange(const Gate& next) const override;   // ...there is no change.

  virtual bool canSimplySwap(const SingleTargetGate& prev) const override;

  virtual bool canSimplySwap(const DiagonalGate& prev) const override;

  virtual bool canSwapSwap(const SwapGate& prev) const override;
  virtual std::unique_ptr<Gate> rightSwapMoverChange(const SwapGate& prev) const override;
  virtual std::unique_ptr<Gate> leftSwapMoverChange(const SwapGate& prev) const override;

protected:
  void reset_simulator_index_manager();
  void reset_simulators();

  bool equivalent_structure(const Gate& rhs) const override;  // (Could these
  bool equals_(const Gate& rhs) const override;               // three functions
  bool sort_before(const Gate& rhs) const override;           // be private?)
  virtual void randomize_qbits();
  bool is_involved(int qbit) const;
  bool target_is_involved(const SingleTargetGate& other) const;  // True if target of 'this' is the target or a...
                                                                 // ...control of 'other'.
  bool target_is_control(const SingleTargetGate& other) const;  // True if target of 'this' is a control of 'other'.

private:
  // The matrix elements depend on the details of the derived Gate class, e.g. whether the gate is a NOT gate or a
  // Hadamard.
  virtual std::vector<std::vector<cmplx> > matrix_() const = 0;

protected:
  int target_;
  Controls controls_;
  std::shared_ptr<SimulatorIndexManager> simulator_index_manager;
  std::shared_ptr<GateSimulator> forward_simulator;
  std::shared_ptr<GateSimulator> reverse_simulator;
};


class RotationGate : public virtual SingleTargetGate
{
  // Base class for Pauli rotations. Typically, angle is restricted to be between -2pi and 2pi - if a mutation takes it
  // out of this range, then a suitable multiple of 4pi is added/subtracted. For uncontrolled gates or ArbitraryPhase,
  // we only need a range of -pi to pi. In this case, a suitable multiple of 2pi is added/subtracted to get into range.
public:
  explicit RotationGate(const CircuitContext& context);
  RotationGate(const CircuitContext& context, int target, const Controls& controls, double angle);

  int numParameters() const override;
  double parameter(int paramNum) const override;  // Note that we already have an accessor function: angle().
  void setParameter(int paramNum, double value) override;
  void randomizeParameters() override;

  void random() override;
  void mutate() override;

  double angle() const;

  bool mutatable() const override;

  void applyPartTo(State& state, int partIndex) const override;
  void optimizeAngle(const std::vector<std::vector<double> >& errorCoeffs) override;

  void output(std::ostream& out) const override;

private:
  bool equals_(const Gate& rhs) const override;
  bool sort_before(const Gate& rhs) const override;

  // Simulators for the various parts of the gate matrix, e.g. the constant part, cosine part and sine part. (One of
  // these simulators - typically the one for the constant part - also handles copying the coefficient when a control is
  // unset, while the others should ignore these coefficients.) Used for 'symbolic simulation', i.e. simulation using
  // the general rotation matrix, leaving the rotation angle unspecified. This allows us to quickly find the best angle
  // for a gate being inserted into a circuit.
  virtual const GateSimulator& part_simulator(int partIndex) const = 0;

  // While the optimizeAngle() function does most of the work, it doesn't know that the angle, theta, being optimized is
  // not necessarily the gate angle. Hence we have a virtual function to provide this final conversion.
  virtual double thetaToAngle(double theta) const = 0;  // While angle is usually theta/2, making this the default...
                                                        // ...feels wrong.
protected:
  double angle_;
};


class DiagonalGate : public virtual SingleTargetGate
{
  // A class to handle all the gates that have diagonal matrices, i.e. PhaseGate, PhaseInv, PiByEight, PiByEightInv,
  // ZGate, ZRotation and ArbitraryPhase. This is provided merely to simplify the simplification code, where such code
  // considers the commutation of gates.
  // Note that this is NOT the set of gates where target and control bits may be interchanged - the ZRotation gate does
  // not allow this.
public:
  explicit DiagonalGate(const CircuitContext& context);
  DiagonalGate(const CircuitContext& context, int target, const Controls& controls);  // (Do we want a constructor...
                                                                                      // ...taking a Controls&& too?)
  // The following two overrides are technically unnecessary - matricesCommute() is never called for DiagonalGates.
  // However, if they were, then the default (SingleTargetGate) versions would give the wrong result. We therefore
  // include them for accuracy.
  bool matricesCommute(const SingleTargetGate& other) const override;
  bool matricesCommute(const DiagonalGate& other) const override;

  // Overrides for the simple swap routines - the reason for this class.
  bool cancelsAtStart() const override;
  bool canSimplySwap(const Gate& next) const override;
  bool canSimplySwap(const SingleTargetGate& prev) const override;
  bool canSimplySwap(const DiagonalGate& prev) const override;
  bool canSimplySwap(const Oracle& prev) const override;
};


class PhaseTypeGate : public DiagonalGate
{
  // A class to handle all the diagonal gates that have 1 as their first matrix element. The special treatment of these
  // gates is solely due to the fact that, for such gates, a target can be switched with a control without leading to
  // any real change. To prevent circuits being considered different just because two Phase gates, affecting the same
  // bits, have chosen to label a different one as the target, we ensure that such gates are always placed in 'canonical
  // form', i.e. out of those qbits involved in the gate, the first one is always the target.
public:
  explicit PhaseTypeGate(const CircuitContext& context);
  PhaseTypeGate(const CircuitContext& context, int target, const Controls& controls);  // Will switch target and...
                                                                                       // ...first control if necessary.
  PhaseTypeGate(const CircuitContext& context, const Controls& controls);  // Target is the first control.

  int numQbitOptions() const override;

  void swapBits(int i, int j) override;

  bool cancelsAtStart() const override;

protected:
  void randomize_qbits() override;
};


class XTypeGate : public virtual SingleTargetGate
{
  // Handles all gates that handle X Rotations, up to a phase. This class has only been created in an attempt to reduce
  // the number of functions involved in implementing simplification (specifically swapping) of gate pairs via double
  // dispatch. The reason that these gate types have been collected together is that their 2x2 matrices commute, meaning
  // that such gates may be swapped more often than standard SingleTargetGates.
  // (There must be a neater way of simplifying all this 'double dispatch'!)
public:
  explicit XTypeGate(const CircuitContext& context);
  XTypeGate(const CircuitContext& context, int target, const Controls& controls);

  bool matricesCommute(const SingleTargetGate& other) const override;
  bool matricesCommute(const XTypeGate& other) const override;
};


class YTypeGate : public virtual SingleTargetGate
{
  // Handles all gates that handle Y Rotations, up to a phase. This class has only been created in an attempt to reduce
  // the number of functions involved in implementing simplification (specifically swapping) of gate pairs via double
  // dispatch. The reason that these gate types have been collected together is that their 2x2 matrices commute, meaning
  // that such gates may be swapped more often standard SingleTargetGates.
  // (There must be a neater way of simplifying all this 'double dispatch'!)
public:
  explicit YTypeGate(const CircuitContext& context);
  YTypeGate(const CircuitContext& context, int target, const Controls& controls);

  bool matricesCommute(const SingleTargetGate& other) const override;
  bool matricesCommute(const YTypeGate& other) const override;
};



class Hadamard final : public SingleTargetGate
{
public:
  explicit Hadamard(const CircuitContext& context);
  Hadamard(const CircuitContext& context, int target, const Controls& controls);

  int gateTypeId() const override;  // (I prefer to avoid 'if type == ' code, but am yet to find a better way.)

  std::unique_ptr<Gate> clone() const override;
  Hadamard* rawClone() const override;  // Used only in the simplification routines to avoid inefficient shared_ptrs.

  std::string name() const override;

  static bool gateAvailable(const CircuitContext& context, int target, const Controls& controls);
  static long gateCost(const CircuitContext& context, int numControls);


  // The remaining public member functions are those for simplification, first the 'first-dispatch' versions, then the
  // 'second-dispatch' versions grouped by gate type.
  long canSimplify(const Gate& next) const override;  // Returns cost improvement.
  GateSequence simplification(const Gate& next) const override;
  bool canHSwap(const Gate& next) const override;
  std::unique_ptr<Gate> rightHMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.
  std::unique_ptr<Gate> leftHMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.

  long canSimplify(const Hadamard& prev) const override;
  GateSequence simplification(const Hadamard& prev) const override;

  long canSimplify(const XGate& prev) const override;
  GateSequence simplification(const XGate& prev) const override;  // HX->ZH is a 'simplification' if the cost changes.
  bool canHSwap(const XGate& prev) const override;  // HX->ZH is a HSwap if X and Z have the same cost.
  std::unique_ptr<Gate> rightHMoverChange(const XGate& prev) const override;
  std::unique_ptr<Gate> leftHMoverChange(const XGate& prev) const override;

  bool canHSwap(const YGate& prev) const override;
  std::unique_ptr<Gate> rightHMoverChange(const YGate& prev) const override;  // These always return nullptr, which...
  std::unique_ptr<Gate> leftHMoverChange(const YGate& prev) const override;   // ...indicates no change.

  long canSimplify(const ZGate& prev) const override;
  GateSequence simplification(const ZGate& prev) const override;
  bool canHSwap(const ZGate& prev) const override;
  std::unique_ptr<Gate> rightHMoverChange(const ZGate& prev) const override;
  std::unique_ptr<Gate> leftHMoverChange(const ZGate& prev) const override;

  long canSimplify(const XRotation& prev) const override;
  GateSequence simplification(const XRotation& prev) const override;
  bool canHSwap(const XRotation& prev) const override;
  std::unique_ptr<Gate> rightHMoverChange(const XRotation& prev) const override;
  std::unique_ptr<Gate> leftHMoverChange(const XRotation& prev) const override;

  bool canHSwap(const YRotation& prev) const override;
  std::unique_ptr<Gate> rightHMoverChange(const YRotation& prev) const override;
  std::unique_ptr<Gate> leftHMoverChange(const YRotation& prev) const override;

  long canSimplify(const ZRotation& prev) const override;
  GateSequence simplification(const ZRotation& prev) const override;
  bool canHSwap(const ZRotation& prev) const override;
  std::unique_ptr<Gate> rightHMoverChange(const ZRotation& prev) const override;
  std::unique_ptr<Gate> leftHMoverChange(const ZRotation& prev) const override;

private:
  static const int gate_type_id = 0;
  std::vector<std::vector<cmplx> > matrix_() const override;
};


class PiByEight final : public PhaseTypeGate
{
  // Otherwise called T
public:
  explicit PiByEight(const CircuitContext& context);
  PiByEight(const CircuitContext& context, int target, const Controls& controls);
  PiByEight(const CircuitContext& context, const Controls& controls);  // Target and controls are interchangeable.
                                                                       // Select one to be 'the' target.
  int gateTypeId() const override;  // I try to avoid 'if type == ' code, but have yet to find a better way.

  std::unique_ptr<Gate> clone() const override;
  PiByEight* rawClone() const override;  // Used only in the simplification routines to avoid inefficient shared_ptrs.

  std::string name() const override;

  static bool gateAvailable(const CircuitContext& context, int target, const Controls& controls);
  static bool gateAvailable(const CircuitContext& context, const Controls& controls);  // Target is any one of the...
  static long gateCost(const CircuitContext& context, int numControls);                // ...'controls'.


  // The remaining public member functions are overrides for the simplification routines. First are those for 'first
  // dispatch', then come those for 'second dispatch' grouped by gate type.
  long canSimplify(const Gate& next) const override;  // Returns cost improvement.
  GateSequence simplification(const Gate& next) const override;

  long canSimplify(const PiByEight& prev) const override;
  GateSequence simplification(const PiByEight& prev) const override;

  long canSimplify(const PiByEightInv& prev) const override;
  GateSequence simplification(const PiByEightInv& prev) const override;

  long canSimplify(const PhaseGate& prev) const override;
  GateSequence simplification(const PhaseGate& prev) const override;

  long canSimplify(const PhaseInv& prev) const override;
  GateSequence simplification(const PhaseInv& prev) const override;

  long canSimplify(const ZGate& prev) const override;
  GateSequence simplification(const ZGate& prev) const override;

  long canSimplify(const ZRotation& prev) const override;
  GateSequence simplification(const ZRotation& prev) const override;

  long canSimplify(const ArbitraryPhase& prev) const override;
  GateSequence simplification(const ArbitraryPhase& prev) const override;

private:
  static const int gate_type_id = 1;  
  std::vector<std::vector<cmplx> > matrix_() const override;
};


class PiByEightInv final : public PhaseTypeGate
{
  // Otherwise called TInv
public:
  explicit PiByEightInv(const CircuitContext& context);
  PiByEightInv(const CircuitContext& context, int target, const Controls& controls);
  PiByEightInv(const CircuitContext& context, const Controls& controls);  // Target and controls are interchangeable.
                                                                          // Select one to be 'the' target.
  int gateTypeId() const override;  // I prefer to avoid 'if type == ' code, but I haven't found a better way yet.

  std::unique_ptr<Gate> clone() const override;
  PiByEightInv* rawClone() const override;  // Used only in simplification routines to avoid inefficient shared_ptrs.

  std::string name() const override;

  static bool gateAvailable(const CircuitContext& context, int target, const Controls& controls);
  static bool gateAvailable(const CircuitContext& context, const Controls& controls);  // Target is any one of the...
  static long gateCost(const CircuitContext& context, int numControls);                // ...'controls'.


  // The remaining public member functions are overrides for the simplification routines. First are those for 'first
  // dispatch', then come those for 'second dispatch' grouped by gate type.
  long canSimplify(const Gate& next) const override;  // Returns cost improvement.
  GateSequence simplification(const Gate& next) const override;

  long canSimplify(const PiByEight& prev) const override;
  GateSequence simplification(const PiByEight& prev) const override;

  long canSimplify(const PiByEightInv& prev) const override;
  GateSequence simplification(const PiByEightInv& prev) const override;

  long canSimplify(const PhaseGate& prev) const override;
  GateSequence simplification(const PhaseGate& prev) const override;

  long canSimplify(const PhaseInv& prev) const override;
  GateSequence simplification(const PhaseInv& prev) const override;

  long canSimplify(const ZGate& prev) const override;
  GateSequence simplification(const ZGate& prev) const override;

  long canSimplify(const ZRotation& prev) const override;
  GateSequence simplification(const ZRotation& prev) const override;

  long canSimplify(const ArbitraryPhase& prev) const override;
  GateSequence simplification(const ArbitraryPhase& prev) const override;

private:
  static const int gate_type_id = 2;  
  std::vector<std::vector<cmplx> > matrix_() const override;
};


class PhaseGate final : public PhaseTypeGate
{
  // Otherwise called S
public:
  explicit PhaseGate(const CircuitContext& context);
  PhaseGate(const CircuitContext& context, int target, const Controls& controls);
  PhaseGate(const CircuitContext& context, const Controls& controls);  // Target and controls are interchangeable.
                                                                       // Select one to be 'the' target.
  int gateTypeId() const override;  // I prefer to avoid 'if type == ' code, but I haven't found a better way yet.

  std::unique_ptr<Gate> clone() const override;
  PhaseGate* rawClone() const override;  // Used only in the simplification routines to avoid inefficient shared_ptrs.

  std::string name() const override;

  static bool gateAvailable(const CircuitContext& context, int target, const Controls& controls);
  static bool gateAvailable(const CircuitContext& context, const Controls& controls);  // Target is any one of the...
  static long gateCost(const CircuitContext& context, int numControls);                // ...'controls'.


  // The remaining public member functions are overrides for the simplification routines. First are those for 'first
  // dispatch', then come those for 'second dispatch' grouped by gate type.
  long canSimplify(const Gate& next) const override;  // Returns cost improvement.
  GateSequence simplification(const Gate& next) const override;
  bool canRSwap(const Gate& next) const override;
  std::unique_ptr<Gate> rightRMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.
  std::unique_ptr<Gate> leftRMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.

  long canSimplify(const PiByEight& prev) const override;
  GateSequence simplification(const PiByEight& prev) const override;

  long canSimplify(const PiByEightInv& prev) const override;
  GateSequence simplification(const PiByEightInv& prev) const override;

  long canSimplify(const PhaseGate& prev) const override;
  GateSequence simplification(const PhaseGate& prev) const override;

  long canSimplify(const PhaseInv& prev) const override;
  GateSequence simplification(const PhaseInv& prev) const override;

  long canSimplify(const ZGate& prev) const override;
  GateSequence simplification(const ZGate& prev) const override;

  long canSimplify(const XRotation& prev) const override;
  GateSequence simplification(const XRotation& prev) const override;
  bool canRSwap(const XRotation& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const XRotation& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const XRotation& prev) const override;

  long canSimplify(const YRotation& prev) const override;
  GateSequence simplification(const YRotation& prev) const override;
  bool canRSwap(const YRotation& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const YRotation& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const YRotation& prev) const override;

  long canSimplify(const ZRotation& prev) const override;
  GateSequence simplification(const ZRotation& prev) const override;

  long canSimplify(const ArbitraryPhase& prev) const override;
  GateSequence simplification(const ArbitraryPhase& prev) const override;

private:
  static const int gate_type_id = 3;  
  std::vector<std::vector<cmplx> > matrix_() const override;
};


class PhaseInv final : public PhaseTypeGate
{
  // Otherwise called SInv
public:
  explicit PhaseInv(const CircuitContext& context);
  PhaseInv(const CircuitContext& context, int target, const Controls& controls);
  PhaseInv(const CircuitContext& context, const Controls& controls);  // Target and controls are interchangeable.
                                                                      // Select one to be 'the' target.
  int gateTypeId() const override;  // I prefer to avoid 'if type == ' code, but I haven't found a better way yet.

  std::unique_ptr<Gate> clone() const override;
  PhaseInv* rawClone() const override;  // Used only in the simplification routines to avoid inefficient shared_ptrs.

  std::string name() const override;

  static bool gateAvailable(const CircuitContext& context, int target, const Controls& controls);
  static bool gateAvailable(const CircuitContext& context, const Controls& controls);  // Target is any one of the...
  static long gateCost(const CircuitContext& context, int numControls);                // ...'controls'.


  // The remaining public member functions are overrides for the simplification routines. First are those for 'first
  // dispatch', then come those for 'second dispatch' grouped by gate type.
  long canSimplify(const Gate& next) const override;  // Returns cost improvement.
  GateSequence simplification(const Gate& next) const override;
  bool canRSwap(const Gate& next) const override;  // Calls next.canRSwap(*this). (Double dispatch.)
  std::unique_ptr<Gate> rightRMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.
  std::unique_ptr<Gate> leftRMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.

  long canSimplify(const PiByEight& prev) const override;
  GateSequence simplification(const PiByEight& prev) const override;

  long canSimplify(const PiByEightInv& prev) const override;
  GateSequence simplification(const PiByEightInv& prev) const override;

  long canSimplify(const PhaseGate& prev) const override;
  GateSequence simplification(const PhaseGate& prev) const override;

  long canSimplify(const PhaseInv& prev) const override;
  GateSequence simplification(const PhaseInv& prev) const override;

  long canSimplify(const ZGate& prev) const override;
  GateSequence simplification(const ZGate& prev) const override;

  long canSimplify(const XRotation& prev) const override;
  GateSequence simplification(const XRotation& prev) const override;
  bool canRSwap(const XRotation& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const XRotation& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const XRotation& prev) const override;

  long canSimplify(const YRotation& prev) const override;
  GateSequence simplification(const YRotation& prev) const override;
  bool canRSwap(const YRotation& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const YRotation& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const YRotation& prev) const override;

  long canSimplify(const ZRotation& prev) const override;
  GateSequence simplification(const ZRotation& prev) const override;

  long canSimplify(const ArbitraryPhase& prev) const override;
  GateSequence simplification(const ArbitraryPhase& prev) const override;

private:
  static const int gate_type_id = 4;  
  std::vector<std::vector<cmplx> > matrix_() const override;
};


class XGate final : public XTypeGate
{
public:
  explicit XGate(const CircuitContext& context);
  XGate(const CircuitContext& context, int target, const Controls& controls);

  int gateTypeId() const override;  // I prefer to avoid 'if type == ' code, but I haven't found a better way yet.

  std::unique_ptr<Gate> clone() const override;
  XGate* rawClone() const override;  // Used only in the simplification routines to avoid inefficient shared_ptrs.

  std::string name() const override;

  static bool gateAvailable(const CircuitContext& context, int target, const Controls& controls);
  static long gateCost(const CircuitContext& context, int numControls);

  virtual bool matricesAnticommute(const SingleTargetGate& other) const override;
  virtual bool matricesAnticommute(const YGate& other) const override;
  virtual bool matricesAnticommute(const ZGate& other) const override;


  // The remaining public member functions are overrides for the simplification routines. First are those for 'first
  // dispatch', then come those for 'second dispatch' grouped by gate type.
  long canSimplify(const Gate& next) const override;  // Returns cost improvement.
  GateSequence simplification(const Gate& next) const override;
  bool canHSwap(const Gate& next) const override;
  std::unique_ptr<Gate> rightHMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.
  std::unique_ptr<Gate> leftHMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.
  bool canRSwap(const Gate& next) const override;
  std::unique_ptr<Gate> rightRMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.
  std::unique_ptr<Gate> leftRMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.

  long canSimplify(const Hadamard& prev) const override;
  GateSequence simplification(const Hadamard& prev) const override;
  bool canHSwap(const Hadamard& prev) const override;
  std::unique_ptr<Gate> rightHMoverChange(const Hadamard& prev) const override;
  std::unique_ptr<Gate> leftHMoverChange(const Hadamard& prev) const override;

  long canSimplify(const XGate& prev) const override;
  GateSequence simplification(const XGate& prev) const override;

  long canSimplify(const YGate& prev) const override;
  GateSequence simplification(const YGate& prev) const override;

  long canSimplify(const ZGate& prev) const override;
  GateSequence simplification(const ZGate& prev) const override;

  long canSimplify(const XRotation& prev) const override;
  GateSequence simplification(const XRotation& prev) const override;

  bool canRSwap(const YRotation& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const YRotation& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const YRotation& prev) const override;

  bool canRSwap(const ZRotation& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const ZRotation& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const ZRotation& prev) const override;

  bool canRSwap(const ArbitraryPhase& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const ArbitraryPhase& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const ArbitraryPhase& prev) const override;

private:
  static const int gate_type_id = 5;  
  std::vector<std::vector<cmplx> > matrix_() const override;
};


class YGate final : public YTypeGate
{
public:
  explicit YGate(const CircuitContext& context);
  YGate(const CircuitContext& context, int target, const Controls& controls);

  int gateTypeId() const override;  // I prefer to avoid 'if type == ' code, but I haven't found a better way yet.

  std::unique_ptr<Gate> clone() const override;
  YGate* rawClone() const override;  // Used only in the simplification routines to avoid inefficient shared_ptrs.

  std::string name() const override;

  static bool gateAvailable(const CircuitContext& context, int target, const Controls& controls);
  static long gateCost(const CircuitContext& context, int numControls);

  virtual bool matricesAnticommute(const SingleTargetGate& other) const override;
  virtual bool matricesAnticommute(const XGate& other) const override;
  virtual bool matricesAnticommute(const ZGate& other) const override;


  // The remaining public member functions are overrides for the simplification routines. First are those for 'first
  // dispatch', then come those for 'second dispatch' grouped by gate type.
  long canSimplify(const Gate& next) const override;  // Returns cost improvement.
  GateSequence simplification(const Gate& next) const override;
  bool canHSwap(const Gate& next) const override;  // Calls next.canHSwap(*this). (Double dispatch.)
  std::unique_ptr<Gate> rightHMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.
  std::unique_ptr<Gate> leftHMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.
  bool canRSwap(const Gate& next) const override;  // Calls next.canRSwap(*this). (Double dispatch.)
  std::unique_ptr<Gate> rightRMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.
  std::unique_ptr<Gate> leftRMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.

  bool canHSwap(const Hadamard& prev) const override;
  std::unique_ptr<Gate> rightHMoverChange(const Hadamard& prev) const override;  // These always returns nullptr,...
  std::unique_ptr<Gate> leftHMoverChange(const Hadamard& prev) const override;   // ...indicating no change.

  long canSimplify(const XGate& prev) const override;
  GateSequence simplification(const XGate& prev) const override;

  long canSimplify(const YGate& prev) const override;
  GateSequence simplification(const YGate& prev) const override;

  long canSimplify(const ZGate& prev) const override;
  GateSequence simplification(const ZGate& prev) const override;

  bool canRSwap(const XRotation& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const XRotation& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const XRotation& prev) const override;

  long canSimplify(const YRotation& prev) const override;
  GateSequence simplification(const YRotation& prev) const override;

  bool canRSwap(const ZRotation& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const ZRotation& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const ZRotation& prev) const override;

  bool canRSwap(const ArbitraryPhase& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const ArbitraryPhase& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const ArbitraryPhase& prev) const override;

private:
  static const int gate_type_id = 6;
  std::vector<std::vector<cmplx> > matrix_() const override;
};


class ZGate final : public PhaseTypeGate
{
public:
  explicit ZGate(const CircuitContext& context);
  ZGate(const CircuitContext& context, int target, const Controls& controls);
  ZGate(const CircuitContext& context, const Controls& controls);  // Target and controls are interchangeable. Select...
                                                                   // ...one to be 'the' target.
  int gateTypeId() const override;  // I prefer to avoid 'if type == ' code, but I haven't found a better way yet.

  std::unique_ptr<Gate> clone() const override;
  ZGate* rawClone() const override;  // Used only in the simplification routines to avoid inefficient shared_ptrs.

  std::string name() const override;

  static bool gateAvailable(const CircuitContext& context, int target, const Controls& controls);
  static bool gateAvailable(const CircuitContext& context, const Controls& controls);  // Target is any one of the...
  static long gateCost(const CircuitContext& context, int numControls);                // ...'controls'.

  virtual bool matricesAnticommute(const SingleTargetGate& other) const override;
  virtual bool matricesAnticommute(const XGate& other) const override;
  virtual bool matricesAnticommute(const YGate& other) const override;


  // The remaining public member functions are overrides for the simplification routines. First are those for 'first
  // dispatch', then come those for 'second dispatch' grouped by gate type.
  long canSimplify(const Gate& next) const override;  // Returns cost improvement.
  GateSequence simplification(const Gate& next) const override;
  bool canHSwap(const Gate& next) const override;  // Calls next.canHSwap(*this). (Double dispatch.)
  std::unique_ptr<Gate> rightHMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.
  std::unique_ptr<Gate> leftHMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.
  bool canRSwap(const Gate& next) const override;  // Calls next.canRSwap(*this). (Double dispatch.)
  std::unique_ptr<Gate> rightRMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.
  std::unique_ptr<Gate> leftRMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.


  long canSimplify(const Hadamard& prev) const override;
  GateSequence simplification(const Hadamard& prev) const override;
  bool canHSwap(const Hadamard& prev) const override;
  std::unique_ptr<Gate> rightHMoverChange(const Hadamard& prev) const override;
  std::unique_ptr<Gate> leftHMoverChange(const Hadamard& prev) const override;

  long canSimplify(const PiByEight& prev) const override;
  GateSequence simplification(const PiByEight& prev) const override;

  long canSimplify(const PiByEightInv& prev) const override;
  GateSequence simplification(const PiByEightInv& prev) const override;

  long canSimplify(const PhaseGate& prev) const override;
  GateSequence simplification(const PhaseGate& prev) const override;

  long canSimplify(const PhaseInv& prev) const override;
  GateSequence simplification(const PhaseInv& prev) const override;

  long canSimplify(const XGate& prev) const override;
  GateSequence simplification(const XGate& prev) const override;

  long canSimplify(const YGate& prev) const override;
  GateSequence simplification(const YGate& prev) const override;

  long canSimplify(const ZGate& prev) const override;
  GateSequence simplification(const ZGate& prev) const override;

  long canSimplify(const XRotation& prev) const override;
  GateSequence simplification(const XRotation& prev) const override;
  bool canRSwap(const XRotation& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const XRotation& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const XRotation& prev) const override;

  long canSimplify(const YRotation& prev) const override;
  GateSequence simplification(const YRotation& prev) const override;
  bool canRSwap(const YRotation& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const YRotation& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const YRotation& prev) const override;

  long canSimplify(const ZRotation& prev) const override;
  GateSequence simplification(const ZRotation& prev) const override;

  long canSimplify(const ArbitraryPhase& prev) const override;
  GateSequence simplification(const ArbitraryPhase& prev) const override;

private:
  static const int gate_type_id = 7;  
  std::vector<std::vector<cmplx> > matrix_() const override;
};


class XRotation final : public RotationGate, public XTypeGate
{
public:
  explicit XRotation(const CircuitContext& context);
  XRotation(const CircuitContext& context, int target, const Controls& controls, double angle);

  int gateTypeId() const override;  // I prefer to avoid 'if type == ' code, but I haven't found a better way yet.

  std::unique_ptr<Gate> clone() const override;
  XRotation* rawClone() const override;  // Used only in the simplification routines to avoid inefficient shared_ptrs.

  std::string name() const override;

  static bool gateAvailable(const CircuitContext& context, int target, const Controls& controls);
  static long gateCost(const CircuitContext& context, int numControls);


  // The remaining public member functions are overrides for the simplification routines. First are those for gate
  // 'reduction', i.e. the detection of special values of the angle parameter that allow the gate to be replaced by a
  // cheaper non-parameterized gate. Then come the 'first dispatch' simplification and swap functions. Finally, we have
  // the 'second dispatch' functions grouped by gate type.

  bool canReduce() const override;
  GateSequence reduction() const override;

  long canSimplify(const Gate& next) const override;  // Returns cost improvement.
  GateSequence simplification(const Gate& next) const override;
  bool canHSwap(const Gate& next) const override;  // Calls next.canHSwap(*this). (Double dispatch.)
  std::unique_ptr<Gate> rightHMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.
  std::unique_ptr<Gate> leftHMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.
  bool canRSwap(const Gate& next) const override;  // Calls next.canRSwap(*this). (Double dispatch.)
  std::unique_ptr<Gate> rightRMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.
  std::unique_ptr<Gate> leftRMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.

  long canSimplify(const Hadamard& prev) const override;
  GateSequence simplification(const Hadamard& prev) const override;
  bool canHSwap(const Hadamard& prev) const override;
  std::unique_ptr<Gate> rightHMoverChange(const Hadamard& prev) const override;
  std::unique_ptr<Gate> leftHMoverChange(const Hadamard& prev) const override;

  long canSimplify(const PhaseGate& prev) const override;
  GateSequence simplification(const PhaseGate& prev) const override;
  bool canRSwap(const PhaseGate& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const PhaseGate& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const PhaseGate& prev) const override;

  long canSimplify(const PhaseInv& prev) const override;
  GateSequence simplification(const PhaseInv& prev) const override;
  bool canRSwap(const PhaseInv& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const PhaseInv& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const PhaseInv& prev) const override;

  long canSimplify(const XGate& prev) const override;
  GateSequence simplification(const XGate& prev) const override;

  bool canRSwap(const YGate& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const YGate& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const YGate& prev) const override;

  long canSimplify(const ZGate& prev) const override;
  GateSequence simplification(const ZGate& prev) const override;
  bool canRSwap(const ZGate& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const ZGate& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const ZGate& prev) const override;

  long canSimplify(const XRotation& prev) const override;
  GateSequence simplification(const XRotation& prev) const override;

private:
  const GateSimulator& part_simulator(int partIndex) const override;
  double thetaToAngle(double theta) const override;

private:
  static const int gate_type_id = 8;
  std::vector<std::vector<cmplx> > matrix_() const override;

  static const SparseGateSimulator const_sim;
  static const SparseGateSimulator cos_sim;
  static const SparseGateSimulator sin_sim;
};


class YRotation final : public RotationGate, public YTypeGate
{
public:
  explicit YRotation(const CircuitContext& context);
  YRotation(const CircuitContext& context, int target, const Controls& controls, double angle);

  int gateTypeId() const override;  // I prefer to avoid 'if type == ' code, but I haven't found a better way yet.

  std::unique_ptr<Gate> clone() const override;
  YRotation* rawClone() const override;  // Used only in the simplification routines to avoid inefficient shared_ptrs.

  std::string name() const override;

  static bool gateAvailable(const CircuitContext& context, int target, const Controls& controls);
  static long gateCost(const CircuitContext& context, int numControls);


  // The remaining public member functions are overrides for the simplification routines. First are those for gate
  //'reduction', i.e. the detection of special values of the angle parameter that allow the gate to be replaced by a
  // cheaper non-parameterized gate. Then come the 'first dispatch' simplification and swap functions. Finally, we have
  // the 'second dispatch' functions grouped by gate type.

  bool canReduce() const override;
  GateSequence reduction() const override;

  long canSimplify(const Gate& next) const override;  // Returns cost improvement.
  GateSequence simplification(const Gate& next) const override;
  bool canHSwap(const Gate& next) const override;  // Calls next.canHSwap(*this). (Double dispatch.)
  std::unique_ptr<Gate> rightHMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.
  std::unique_ptr<Gate> leftHMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.
  bool canRSwap(const Gate& next) const override;  // Calls next.canRSwap(*this). (Double dispatch.)
  std::unique_ptr<Gate> rightRMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.
  std::unique_ptr<Gate> leftRMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.

  bool canHSwap(const Hadamard& prev) const override;
  std::unique_ptr<Gate> rightHMoverChange(const Hadamard& prev) const override;
  std::unique_ptr<Gate> leftHMoverChange(const Hadamard& prev) const override;

  long canSimplify(const PhaseGate& prev) const override;
  GateSequence simplification(const PhaseGate& prev) const override;
  bool canRSwap(const PhaseGate& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const PhaseGate& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const PhaseGate& prev) const override;

  long canSimplify(const PhaseInv& prev) const override;
  GateSequence simplification(const PhaseInv& prev) const override;
  bool canRSwap(const PhaseInv& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const PhaseInv& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const PhaseInv& prev) const override;

  bool canRSwap(const XGate& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const XGate& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const XGate& prev) const override;

  long canSimplify(const YGate& prev) const override;
  GateSequence simplification(const YGate& prev) const override;

  long canSimplify(const ZGate& prev) const override;
  GateSequence simplification(const ZGate& prev) const override;
  bool canRSwap(const ZGate& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const ZGate& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const ZGate& prev) const override;

  long canSimplify(const YRotation& prev) const override;
  GateSequence simplification(const YRotation& prev) const override;

private:
  const GateSimulator& part_simulator(int partIndex) const override;
  double thetaToAngle(double theta) const override;

private:
  static const int gate_type_id = 9;
  std::vector<std::vector<cmplx> > matrix_() const override;

  static const SparseGateSimulator const_sim;
  static const SparseGateSimulator cos_sim;
  static const SparseGateSimulator sin_sim;
};


class ZRotation final : public DiagonalGate, public RotationGate
{
public:
  explicit ZRotation(const CircuitContext& context);
  ZRotation(const CircuitContext& context, int target, const Controls& controls, double angle);

  int gateTypeId() const override;  // I prefer to avoid 'if type == ' code, but I haven't found a better way yet.

  std::unique_ptr<Gate> clone() const override;
  ZRotation* rawClone() const override;  // Used only in the simplification routines to avoid inefficient shared_ptrs.

  std::string name() const override;

  static bool gateAvailable(const CircuitContext& context, int target, const Controls& controls);
  static long gateCost(const CircuitContext& context, int numControls);


  // The remaining public member functions are overrides for the simplification routines. First are those for gate
  // 'reduction', i.e. the detection of special values of the angle parameter that allow the gate to be replaced by a
  // cheaper non-parameterized gate. Then come the 'first dispatch' simplification and swap functions. Finally, we have
  // the 'second dispatch' functions grouped by gate type.

  bool canReduce() const override;
  GateSequence reduction() const override;

  long canSimplify(const Gate& next) const override;  // Returns cost improvement.
  GateSequence simplification(const Gate& next) const override;
  bool canHSwap(const Gate& next) const override;  // Calls next.canHSwap(*this). (Double dispatch.)
  std::unique_ptr<Gate> rightHMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.
  std::unique_ptr<Gate> leftHMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.
  bool canRSwap(const Gate& next) const override;  // Calls next.canRSwap(*this). (Double dispatch.)
  std::unique_ptr<Gate> rightRMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.
  std::unique_ptr<Gate> leftRMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.

  long canSimplify(const Hadamard& prev) const override;
  GateSequence simplification(const Hadamard& prev) const override;
  bool canHSwap(const Hadamard& prev) const override;
  std::unique_ptr<Gate> rightHMoverChange(const Hadamard& prev) const override;
  std::unique_ptr<Gate> leftHMoverChange(const Hadamard& prev) const override;

  long canSimplify(const PiByEight& prev) const override;
  GateSequence simplification(const PiByEight& prev) const override;

  long canSimplify(const PiByEightInv& prev) const override;
  GateSequence simplification(const PiByEightInv& prev) const override;

  long canSimplify(const PhaseGate& prev) const override;
  GateSequence simplification(const PhaseGate& prev) const override;

  long canSimplify(const PhaseInv& prev) const override;
  GateSequence simplification(const PhaseInv& prev) const override;

  bool canRSwap(const XGate& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const XGate& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const XGate& prev) const override;

  bool canRSwap(const YGate& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const YGate& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const YGate& prev) const override;

  long canSimplify(const ZGate& prev) const override;
  GateSequence simplification(const ZGate& prev) const override;

  long canSimplify(const ZRotation& prev) const override;
  GateSequence simplification(const ZRotation& prev) const override;

  long canSimplify(const ArbitraryPhase& prev) const override;
  GateSequence simplification(const ArbitraryPhase& prev) const override;

private:
  const GateSimulator& part_simulator(int partIndex) const override;
  double thetaToAngle(double theta) const override;

private:
  static const int gate_type_id = 10;
  std::vector<std::vector<cmplx> > matrix_() const override;

  static const SparseGateSimulator const_sim;
  static const SparseGateSimulator cos_sim;
  static const SparseGateSimulator sin_sim;
};


class ArbitraryPhase final : public PhaseTypeGate, public RotationGate
{
public:
  explicit ArbitraryPhase(const CircuitContext& context);
  ArbitraryPhase(const CircuitContext& context, int target, const Controls& controls, double angle);
  ArbitraryPhase(const CircuitContext& context, const Controls& controls, double angle);  // Target and controls are...
                                                                  // ...interchangeable. Select one to be 'the' target.
  int gateTypeId() const override;  // I prefer to avoid 'if type == ' code, but I haven't found a better way yet.

  std::unique_ptr<Gate> clone() const override;
  ArbitraryPhase* rawClone() const override;  // Used only in simplification routines to avoid inefficient shared_ptrs.

  void setParameter(int paramNum, double value) override;  // Controlled ArbitraryPhase has a narrower angle range...
                                                           // ...than other controlled RotationGates.
  std::string name() const override;

  static bool gateAvailable(const CircuitContext& context, int target, const Controls& controls);
  static bool gateAvailable(const CircuitContext& context, const Controls& controls);  // Target is any one of the...
  static long gateCost(const CircuitContext& context, int numControls);                // ...'controls'.


  // The remaining public member functions are overrides for the simplification routines. First are those for gate
  // 'reduction', i.e. the detection of special values of the angle parameter that allow the gate to be replaced by a
  // cheaper non-parameterized gate. Then come the 'first dispatch' simplification and swap functions. Finally, we have
  // the 'second dispatch' functions grouped by gate type.

  bool canReduce() const override;
  GateSequence reduction() const override;

  long canSimplify(const Gate& next) const override;  // Returns cost improvement.
  GateSequence simplification(const Gate& next) const override;
  bool canRSwap(const Gate& next) const override;  // Calls next.canRSwap(*this). (Double dispatch.)
  std::unique_ptr<Gate> rightRMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.
  std::unique_ptr<Gate> leftRMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.

  long canSimplify(const PiByEight& prev) const override;
  GateSequence simplification(const PiByEight& prev) const override;

  long canSimplify(const PiByEightInv& prev) const override;
  GateSequence simplification(const PiByEightInv& prev) const override;

  long canSimplify(const PhaseGate& prev) const override;
  GateSequence simplification(const PhaseGate& prev) const override;

  long canSimplify(const PhaseInv& prev) const override;
  GateSequence simplification(const PhaseInv& prev) const override;

  bool canRSwap(const XGate& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const XGate& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const XGate& prev) const override;

  bool canRSwap(const YGate& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const YGate& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const YGate& prev) const override;

  long canSimplify(const ZGate& prev) const override;
  GateSequence simplification(const ZGate& prev) const override;

  long canSimplify(const ZRotation& prev) const override;
  GateSequence simplification(const ZRotation& prev) const override;

  long canSimplify(const ArbitraryPhase& prev) const override;
  GateSequence simplification(const ArbitraryPhase& prev) const override;

private:
  const GateSimulator& part_simulator(int partIndex) const override;
  double thetaToAngle(double theta) const override;

private:
  static const int gate_type_id = 11;
  std::vector<std::vector<cmplx> > matrix_() const override;

  static const SparseGateSimulator const_sim;
  static const SparseGateSimulator cos_sim;
  static const SparseGateSimulator sin_sim;
};


class SwapGate final : public Gate
{
  // Class for the simple swap gate. (I would like to work out how best to implement a controlled swap gate too.)
public:
  explicit SwapGate(const CircuitContext& context);
  SwapGate(const CircuitContext& context, int bit1, int bit2);  // Provided for testing.

  int gateTypeId() const override;  // I prefer to avoid 'if type == ' code, but I haven't found a better way yet.
  int numQbitOptions() const override;

  std::unique_ptr<Gate> clone() const override;
  SwapGate* rawClone() const override;  // Used only in the simplification routines to avoid inefficient shared_ptrs.

  void random() override;
  void swapBits(int i, int j) override;

  int bit1() const;  // Access required when commuting this gate with others.
  int bit2() const;  // Access required when commuting this gate with others.

  bool available() const override;

//  State applyTo(const State& state) const override;
//  State applyInvTo(const State& state) const override;
  void applyTo(State& state) const override;
  void applyInvTo(State& state) const override;

  int target() const override;  // Shouldn't be called - target() is only called when setting control bits and swap...
                                // ...gates have none (at present).
  std::string name() const override;

  static bool gateAvailable(const CircuitContext& context, int bit1, int bit2);
  static long gateCost(const CircuitContext& context);

  void output(std::ostream& out) const override;


  // The remaining public member functions are overrides for the simplification routines. First are those for 'first
  // dispatch', then come those for 'second dispatch' grouped by gate type.
  bool cancelsAtStart() const override;
  long canSimplify(const Gate& next) const override;  // Returns cost improvement.
  GateSequence simplification(const Gate& next) const override;
  bool canSwapSwap(const Gate& next) const override;
  std::unique_ptr<Gate> rightSwapMoverChange(const Gate& next) const override; // Returns nullptr if there is no change.
  std::unique_ptr<Gate> leftSwapMoverChange(const Gate& next) const override;  // Returns nullptr if there is no change.

  bool canSwapSwap(const SingleTargetGate& prev) const override;
  std::unique_ptr<Gate> rightSwapMoverChange(const SingleTargetGate& prev) const override;
  std::unique_ptr<Gate> leftSwapMoverChange(const SingleTargetGate& prev) const override;

  long canSimplify(const SwapGate& prev) const override;
  GateSequence simplification(const SwapGate& prev) const override;
  bool canSwapSwap(const SwapGate& prev) const override;
  std::unique_ptr<Gate> rightSwapMoverChange(const SwapGate& prev) const override;
  std::unique_ptr<Gate> leftSwapMoverChange(const SwapGate& prev) const override;

private:
  bool equivalent_structure(const Gate& rhs) const override;
  bool equals_(const Gate& rhs) const override;
  bool sort_before(const Gate& rhs) const override;
  void fix_order();

private:
  static const int gate_type_id = 12;  

  // Keep swap_bit1 < swap_bit2.
  int swap_bit1;
  int swap_bit2;
  // (Consider adding a shared_ptr to the SwapSimulator here, and applying the gate to a state is a way similar to that
  // for a SingleTargetGate.)
};


class Oracle final : public Gate
{
  // Class for the Oracle, inheriting from Gate as for all other Gate classes.
  // (Was orginally in grover.h, but the Gate base class needs to know about the Oracle gate anyway, in order to perform
  // the double dispatch involved in circuit simplification.)
public:
  explicit Oracle(const CircuitContext& context);

  void linkToMarkedState(std::shared_ptr<const int> markedState) override;

  int gateTypeId() const override;  // I prefer to avoid 'if type == ' code, but I haven't found a better way yet.
  int numQbitOptions() const override;

  std::unique_ptr<Gate> clone() const override;
  Oracle* rawClone() const override;  // Used only in the simplification routines to avoid inefficient shared_ptrs.

  void random() override;
  void swapBits(int i, int j) override;

  bool available() const override;

//  State applyTo(const State& state) const override;
//  State applyInvTo(const State& state) const override;
  void applyTo(State& state) const override;
  void applyInvTo(State& state) const override;

  int target() const override;  // Shouldn't be called - target() is only called when setting control bits and Oracle...
                                // gates can have no controls.
  std::string name() const override;

  static bool gateAvailable(const CircuitContext& context);
  static long gateCost(const CircuitContext& context);


  // The remaining public member functions are overrides for the simplification routines. First are those for 'first
  // dispatch', then come those for 'second dispatch' grouped by gate type.
  bool cancelsAtStart() const override;
  bool canSimplySwap(const Gate& next) const override;
  long canSimplify(const Gate& next) const override;  // Returns cost improvement.
  GateSequence simplification(const Gate& next) const override;

  bool canSimplySwap(const DiagonalGate& prev) const override;
  long canSimplify(const Oracle& prev) const override;
  GateSequence simplification(const Oracle& prev) const override;

private:
  bool equivalent_structure(const Gate& rhs) const override;
  bool equals_(const Gate& rhs) const override;
  bool sort_before(const Gate& rhs) const override;

public:
  static const int gate_type_id = 13;  // Public to allow 'isOracle' checks.

private:
  std::shared_ptr<const int> marked_state = nullptr;  // (Unwieldy?)
};

//----------------------------------------------------------------------------------------------------------------------

// We create a master map from Gate names to 'GateCreators'. Given a list of Gate names, the Problem class can then
// create its own list of GateCreators, allowing the creation of a random gate of arbitrary type, selected from this
// list.

// Template for converting all the Gate constructors into GateCreators.
template <class DerivedGate>
std::unique_ptr<Gate> createGate(const CircuitContext& context)
{
  return std::make_unique<DerivedGate>(context);
}

// A function that takes a circuit context and produces a unique_ptr to a Gate.
using GateCreator = std::unique_ptr<Gate> (*)(const CircuitContext& context);

// The master list of gate creators
extern const std::map<std::string, GateCreator> masterGateCreatorIndex;

//----------------------------------------------------------------------------------------------------------------------

// A helper function for GateSequences.
long gateCost(const GateSequence& gates);

#endif  // GATE_H
