// MOQCD1 (Multi-objective quantum circuit discovery - algorithm 1.)
//
// Copyright (c) 2020  Alan Reynolds (alanreynolds9@gmail.com)
//
// This file is part of MOQCD1.
//
// MOQCD1 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//
// MOQCD1 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with MOQCD1.  If not, see
// <http://www.gnu.org/licenses/>.

//----------------------------------------------------------------------------------------------------------------------

#ifndef GATE_H
#define GATE_H

#define QICLIB_DONT_USE_NLOPT  // We haven't installed NLopt and don't want to use NLopt dependent features
#define ARMA_DONT_USE_WRAPPER

#include <string>
#include <map>
#include <iostream>
#include <memory>
#include <armadillo>
#include "rng.h"
#include "dictionary.h"
#include "constants.h"
#include "transformation.h"
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
class SU2Gate;
class SwapGate;
class Oracle;

using GateSequence = std::vector<std::unique_ptr<Gate> >;

// NOTE ON DOUBLE-DISPATCH: (Requires updating!)
// The routines that allow for gate pair simplification and swapping make extensive use of double-dispatch. While I used
// to use routines in pairs, with a function called somthing like 'simplifies(const Gate& rhs)' calling something like
// 'rhs.isSimplifiedBy(*this)', at present, to  simplify naming, we simply have 'canSimplify(const Gate& next)' calling
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
// results. The second issue is illustrated best by canSwap(). Consider what happens if I call
// Gate::canSwap(const SingleTargetGate& other) directly. Provided 'other' is just a SingleTargetGate, then there is no
// problem. However, if 'other' is a DiagonalGate - a class deriving from SingleTargetGate, then I have a problem. If
// *this is also a DiagonalGate, then I wish to call DiagonalGate::canSwap(const DiagonalGate& other) which always
// returns true. However, DiagonalGate::canSwap(const SingleTargetGate& other) is called instead, potentially producing
// the wrong result.
//
// (I suspect that this latter problem could also be an issue when the functions come in pairs, since one could call
// Gate::isSimplifiedBy(const SingleTargetGate& lhs) when lhs is a DiagonalGate. This is not a problem at present, since
// only the leaf classes override the 'simplification' routines.)
//
// At present, all such functions are public, including those I could easily make private.


int totalGateCost(const GateSequence& gateSequence);


class Gate
{
  // While the Transformation class is how a Gate appears to the quantum state, the Gate class is designed to be
  // manipulated by the optimization algorithm. When applied to a State, the Gate creates a Transformation to be passed
  // along.
public:
  explicit Gate(const CircuitContext& context);
  virtual ~Gate();

  // Function that ensures that the gate gets access to the correct 'circuit data'. (Think of this as the gate getting
  // to know its 'parent circuit.) The only circuit data at present is the state marked by the Oracle gate.
  virtual void imprint(const Circuit& circuit);

  // (Neither of the next two functions actually refer to the gate referenced by 'this'. I would prefer not to have
  // 'if type == ' style code, which is what having a gateTypeId enables, but I haven't worked out a better way yet. The
  // function numQbitOptions() might belong better in the Context object - certainly if I ever create individual Context
  // classes for each gate type, then this would be moved.)
  virtual int gateTypeId() const = 0;  // Returns a number (from 0 to 16) indicating gate type.
  virtual int numQbitOptions() const = 0;  // Number of options for target/control bits, given the gate's 'context'.

  virtual size_t gateOptionId() const;  // Returns a unique number for each combination of gate type, target qbit and...
                                        // ...control qbits. For hashing.
  virtual std::unique_ptr<Gate> clone() const = 0;  // Virtual copy constructor.
  virtual Gate* rawClone() const = 0;  // Use only where unique or shared_ptrs are unsuitable, likely due to efficiency.

  virtual void random() = 0;
  virtual void mutate();  // NEW: Added again to allow 'adjust parameter and reoptimize' functionality.
  virtual void swapBits(int i, int j) = 0;  // Should not be called on Oracle gates!

  virtual int numParameters() const;  // The number of associated continuous parameters (angles). Default is zero.
  virtual double parameter(int paramNum) const;
  virtual void setParameter(int paramNum, double value);
  virtual void randomizeParameters();

  // Functions for getting the availability and costs of the different gate types. These will access data stored in the
  // CircuitContext.
  virtual bool available() const = 0;
  virtual long cost() const;

  // Operator== first checks that the types of gate are the same. If so, calls equals() polymorphically.
  bool equivalentStructure(const Gate& rhs) const; // 'True' <=> same types, targets and controls. Angles may differ.
  bool operator==(const Gate& rhs) const;
  bool operator!=(const Gate& rhs) const;

  // Function sortBefore is only used by the MO framework, as far as I can tell, to calculate population entropy.
  bool sortBefore(const Gate& rhs) const;

  virtual State applyTo(const State& state) const = 0;
  virtual State applyInvTo(const State& state) const = 0;
  virtual State applyGradTo(const State& state, int paramNum) const;  // For derivatives.

  virtual int target() const = 0;  // Should be called ONLY when the gate has a single target.
  virtual bool mutatable() const;  // Alas, can't make this static and virtual.

  virtual std::string name() const = 0;
  virtual void output(std::ostream& out) const;
  virtual void extendedOutput(std::ostream& out) const;  // For testing.


  // The remaining public member functions are for circuit simplification and getting the circuit into canonical form,
  // indicating whether a pair of gates can be simplified or swapped, and the outcome of performing such a
  // simplification or swap. The routines for swapping gates are grouped into three sets. Simple swaps are (most of)
  // those that require no change to the gates. HSwaps are those that involve a Hadamard that shares its target with the
  // other (X, Y, Z, Rotation) gate. The non-Hadamard will usually change, except for the YGate. RSwaps involve a
  // rotation gate and an X, Y or Z. The rotation gate will usually change. SwapSwaps are swaps involving SwapGates,
  // which may (or may not) require a change to a gate.
  // We also include in this section, routines for gate 'reduction' - i.e. routines that take a parameterized gate and
  // determine whether it can be removed or replaced by a cheaper non-parameterized gate. These are useful after
  // numerical optimization, when gates are likely to be reducible.

  virtual bool canReduce() const;
  virtual GateSequence reduction() const;
  // (There is also an argument for having reduction functions that consider neighbouring gates. In particular, the
  // reduction of a controlled XRotation of angle pi to a controlled XGate and a PhaseGate on the controls is unlikely
  // to be taken at present, since the new gates would likely cost more than the old. However, we could detect the
  // presence of, for example, a neighbouring PhaseInv that would cancel the introduced PhaseGate.)

  virtual bool cancelsAtStart() const = 0;  // Should be able to refer to CircuitContext instead?
  virtual std::pair<long, bool> canSimplify(const Gate& next) const = 0;  // Returns cost improvement and whether...
  virtual GateSequence simplification(const Gate& next) const = 0;        // ...#dof increases.
  virtual bool canSimplySwap(const Gate& next) const;
  virtual bool canSwapSwap(const Gate& next) const;  // All swaps involving SwapGates.
  virtual std::unique_ptr<Gate> rightSwapMoverChange(const Gate& next) const;  // These return a null pointer if...
  virtual std::unique_ptr<Gate> leftSwapMoverChange(const Gate& next) const;   // ...there is no change.
  virtual bool canHSwap(const Gate& next) const;  // Swaps involving Hadamards sharing targets with X, Y, Z.
  virtual std::unique_ptr<Gate> rightHMoverChange(const Gate& next) const;  // These return a null pointer if there...
  virtual std::unique_ptr<Gate> leftHMoverChange(const Gate& next) const;  // ...is no change.
  virtual bool canRSwap(const Gate& next) const;  // Swaps involving rotation gates with X, Y, Z, S, S-1.
  virtual std::unique_ptr<Gate> rightRMoverChange(const Gate& next) const;  // These return a null pointer if there...
  virtual std::unique_ptr<Gate> leftRMoverChange(const Gate& next) const;  // ...is no change.

  // The above simplification routines are those that should be called from outside the gate classes. However, we use
  // 'double-dispatch'. Unfortunately, this requires extra virtual functions whenever we add a new derived class, or at
  // least any derived gate class that permits simplification.

  virtual bool canSimplySwap(const SingleTargetGate& prev) const;
  virtual bool canSwapSwap(const SingleTargetGate& prev) const;  // Swaps involving SwapGates.
  virtual std::unique_ptr<Gate> rightSwapMoverChange(const SingleTargetGate& prev) const;  // nullptr <=> no change.
  virtual std::unique_ptr<Gate> leftSwapMoverChange(const SingleTargetGate& prev) const;  // nullptr <=> no change.

  virtual bool canSimplySwap(const DiagonalGate& prev) const;

  virtual std::pair<long, bool> canSimplify(const Hadamard& prev) const;  // Returns cost improvement and whether...
  virtual GateSequence simplification(const Hadamard& prev) const;        // ...#dof increases.
  virtual bool canHSwap(const Hadamard& prev) const;
  virtual std::unique_ptr<Gate> rightHMoverChange(const Hadamard& prev) const;  // nullptr <=> no change.
  virtual std::unique_ptr<Gate> leftHMoverChange(const Hadamard& prev) const;  // nullptr <=> no change.

  virtual std::pair<long, bool> canSimplify(const PiByEight& prev) const;
  virtual GateSequence simplification(const PiByEight& prev) const;

  virtual std::pair<long, bool> canSimplify(const PiByEightInv& prev) const;
  virtual GateSequence simplification(const PiByEightInv& prev) const;

  virtual std::pair<long, bool> canSimplify(const PhaseGate& prev) const;
  virtual GateSequence simplification(const PhaseGate& prev) const;
  virtual bool canRSwap(const PhaseGate& prev) const;
  virtual std::unique_ptr<Gate> rightRMoverChange(const PhaseGate& prev) const;
  virtual std::unique_ptr<Gate> leftRMoverChange(const PhaseGate& prev) const;

  virtual std::pair<long, bool> canSimplify(const PhaseInv& prev) const;
  virtual GateSequence simplification(const PhaseInv& prev) const;
  virtual bool canRSwap(const PhaseInv& prev) const;
  virtual std::unique_ptr<Gate> rightRMoverChange(const PhaseInv& prev) const;
  virtual std::unique_ptr<Gate> leftRMoverChange(const PhaseInv& prev) const;

  virtual std::pair<long, bool> canSimplify(const XGate& prev) const;
  virtual GateSequence simplification(const XGate& prev) const;
  virtual bool canHSwap(const XGate& prev) const;
  virtual std::unique_ptr<Gate> rightHMoverChange(const XGate& prev) const;
  virtual std::unique_ptr<Gate> leftHMoverChange(const XGate& prev) const;
  virtual bool canRSwap(const XGate& prev) const;
  virtual std::unique_ptr<Gate> rightRMoverChange(const XGate& prev) const;
  virtual std::unique_ptr<Gate> leftRMoverChange(const XGate& prev) const;

  virtual std::pair<long, bool> canSimplify(const YGate& prev) const;
  virtual GateSequence simplification(const YGate& prev) const;
  virtual bool canHSwap(const YGate& prev) const;
  virtual std::unique_ptr<Gate> rightHMoverChange(const YGate& prev) const;  // These return a null pointer if there...
  virtual std::unique_ptr<Gate> leftHMoverChange(const YGate& prev) const;   // ...is no change. (Always happens!)
  virtual bool canRSwap(const YGate& prev) const;
  virtual std::unique_ptr<Gate> rightRMoverChange(const YGate& prev) const;  // These always return a null pointer...
  virtual std::unique_ptr<Gate> leftRMoverChange(const YGate& prev) const;   // ...indicating no change.

  virtual std::pair<long, bool> canSimplify(const ZGate& prev) const;
  virtual GateSequence simplification(const ZGate& prev) const;
  virtual bool canHSwap(const ZGate& prev) const;
  virtual std::unique_ptr<Gate> rightHMoverChange(const ZGate& prev) const;
  virtual std::unique_ptr<Gate> leftHMoverChange(const ZGate& prev) const;
  virtual bool canRSwap(const ZGate& prev) const;
  virtual std::unique_ptr<Gate> rightRMoverChange(const ZGate& prev) const;
  virtual std::unique_ptr<Gate> leftRMoverChange(const ZGate& prev) const;

  virtual std::pair<long, bool> canSimplify(const XRotation& prev) const;
  virtual GateSequence simplification(const XRotation& prev) const;
  virtual bool canHSwap(const XRotation& prev) const;
  virtual std::unique_ptr<Gate> rightHMoverChange(const XRotation& prev) const;
  virtual std::unique_ptr<Gate> leftHMoverChange(const XRotation& prev) const;
  virtual bool canRSwap(const XRotation& prev) const;
  virtual std::unique_ptr<Gate> rightRMoverChange(const XRotation& prev) const;
  virtual std::unique_ptr<Gate> leftRMoverChange(const XRotation& prev) const;

  virtual std::pair<long, bool> canSimplify(const YRotation& prev) const;
  virtual GateSequence simplification(const YRotation& prev) const;
  virtual bool canHSwap(const YRotation& prev) const;
  virtual std::unique_ptr<Gate> rightHMoverChange(const YRotation& prev) const;
  virtual std::unique_ptr<Gate> leftHMoverChange(const YRotation& prev) const;
  virtual bool canRSwap(const YRotation& prev) const;
  virtual std::unique_ptr<Gate> rightRMoverChange(const YRotation& prev) const;
  virtual std::unique_ptr<Gate> leftRMoverChange(const YRotation& prev) const;

  virtual std::pair<long, bool> canSimplify(const ZRotation& prev) const;
  virtual GateSequence simplification(const ZRotation& prev) const;
  virtual bool canHSwap(const ZRotation& prev) const;
  virtual std::unique_ptr<Gate> rightHMoverChange(const ZRotation& prev) const;
  virtual std::unique_ptr<Gate> leftHMoverChange(const ZRotation& prev) const;
  virtual bool canRSwap(const ZRotation& prev) const;
  virtual std::unique_ptr<Gate> rightRMoverChange(const ZRotation& prev) const;
  virtual std::unique_ptr<Gate> leftRMoverChange(const ZRotation& prev) const;

  virtual std::pair<long, bool> canSimplify(const ArbitraryPhase& prev) const;
  virtual GateSequence simplification(const ArbitraryPhase& prev) const;
  virtual bool canRSwap(const ArbitraryPhase& prev) const;
  virtual std::unique_ptr<Gate> rightRMoverChange(const ArbitraryPhase& prev) const;
  virtual std::unique_ptr<Gate> leftRMoverChange(const ArbitraryPhase& prev) const;

  virtual std::pair<long, bool> canSimplify(const SU2Gate& prev) const;
  virtual GateSequence simplification(const SU2Gate& prev) const;

  virtual std::pair<long, bool> canSimplify(const SwapGate& prev) const;
  virtual GateSequence simplification(const SwapGate& prev) const;
  virtual bool canSwapSwap(const SwapGate& prev) const;  // Swaps involving SwapGates.
  virtual std::unique_ptr<Gate> rightSwapMoverChange(const SwapGate& prev) const;
  virtual std::unique_ptr<Gate> leftSwapMoverChange(const SwapGate& prev) const;

  virtual bool canSimplySwap(const Oracle& prev) const;
  virtual std::pair<long, bool> canSimplify(const Oracle& prev) const;
  virtual GateSequence simplification(const Oracle& prev) const;

protected:
  int num_qbits() const;

private:
  virtual bool equivalent_structure(const Gate& rhs) const = 0;  // Called ONLY if rhs is of the same type.
  virtual bool equals_(const Gate& rhs) const = 0;  // Called ONLY if rhs is of the same type.
  virtual bool sort_before(const Gate& rhs) const = 0;  // Called ONLY if rhs is of the same type.
  virtual void calculate_option_id() = 0;

protected:
  const CircuitContext* context_;  // (Might it be possible to point to a smaller object - a GateContext for each...
  size_t option_id;                // ...type of gate?)
};

std::ostream& operator<<(std::ostream& out, const Gate& gate);


class SingleTargetGate : public Gate
{
  // Class used to represent gates that work on a single target qbit. Such gates are permitted to have multiple control
  // qbits, but only one qbit is changed. This class is designed to be manipulated by the optimization algorithm and
  // hence doesn't hold the transformation matrix directly, but gate type (indirectly) and angle parameters. When
  // applied to a state, a SingleTargetGate creates a 2x2 matrix and passes it, along with target and control gate
  // indices, to the State object.
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

  virtual int numQbitOptions() const override;

  void random() override;

  const Controls& controls() const;  // Added for use by simplification code.
  int numControls() const;

  virtual void swapBits(int i, int j) override;

  bool available() const override;
  long cost() const override;

  std::vector<int> allQbits() const;  // Returns a sorted vector of all qbits used or affected by the gate. (Slow!)

  // The following three functions are used to avoid the inefficient allQbits().
  // (I should come up with better names. In particular, allQbitsMatch() might give the impression that the targets
  // match and the controls match. This is NOT the case. Instead, the SET of all qbits involved in one gate match the
  // SET of all qbits involved in the other, i.e. allQbits() produces the same results for the two gates.)
  // (Should these be protected rather than public?)
  bool allAreInvolvedIn(const SingleTargetGate& other) const;  // All qbits of 'this' are involved in 'other'. (Subset.)
  bool allQbitsMatch(const SingleTargetGate& other) const;  // The qbits of 'this' are those of 'other'. (Equality.)
  bool allQbitsMatchControlsOf(const SingleTargetGate& other) const;  // The qbits of 'this' are the controls of...
                                                                      // ...other. (Equality.)
  bool controlsAreControlsOf(const SingleTargetGate& other) const;  // The controls of 'this' are controls of 'other'.
                                                                    // (Subset.)

  State applyTo(const State& state) const override;
  State applyInvTo(const State& state) const override;

  int target() const override;

  virtual bool matricesCommute(const SingleTargetGate& other) const;
  virtual bool matricesCommute(const DiagonalGate& other) const;  // Unnecessary, but makes function name more accurate!
  virtual bool matricesCommute(const XTypeGate& other) const;
  virtual bool matricesCommute(const YTypeGate& other) const;
  virtual bool matricesAnticommute(const SingleTargetGate& other) const;  // Note that if 2x2 matrices commute up to...
  virtual bool matricesAnticommute(const XGate& other) const;       // ...a phase factor, then that factor must be...
  virtual bool matricesAnticommute(const YGate& other) const;       // ...+-1, i.e. the matrices commute or anticommute.
  virtual bool matricesAnticommute(const ZGate& other) const;

  void output(std::ostream& out) const override;
  void extendedOutput(std::ostream& out) const override;  // Adds transformation_() to the output.


  // Remaining public member functions are overrides for the simplification routines - starting with the 'first
  // dispatch' versions and then the 'second dispatch' versions grouped by gate type.
  virtual bool cancelsAtStart() const override;  // Should be able to refer to CircuitContext instead?
  virtual bool canSimplySwap(const Gate& next) const override;
  virtual bool canSwapSwap(const Gate& next) const override;
  virtual std::unique_ptr<Gate> rightSwapMoverChange(const Gate& next) const override; // These return a null pointer..
  virtual std::unique_ptr<Gate> leftSwapMoverChange(const Gate& next) const override;  // ...if there is no change.

  virtual bool canSimplySwap(const SingleTargetGate& prev) const override;

  virtual bool canSimplySwap(const DiagonalGate& prev) const override;

  virtual bool canSwapSwap(const SwapGate& prev) const override;
  virtual std::unique_ptr<Gate> rightSwapMoverChange(const SwapGate& prev) const override; // These returns a null...
  virtual std::unique_ptr<Gate> leftSwapMoverChange(const SwapGate& prev) const override;  // ...pointer if there is...
                                                                                           // ...no change.

protected:
  bool equivalent_structure(const Gate& rhs) const override;  // (Could these
  bool equals_(const Gate& rhs) const override;               //  three functions
  bool sort_before(const Gate& rhs) const override;           //  be private?)
  virtual void calculate_option_id() override;
  virtual void randomize_qbits();

  bool is_involved(int qbit) const;
  bool target_is_involved(const SingleTargetGate& other) const;  // Returns true if the target of 'this' is also the...
                                                                 // ...target or a control of 'other'.
  bool target_is_control(const SingleTargetGate& other) const;  // Returns true if the target of 'this' is a control...
                                                                // of other.
private:
  Transformation transformation_() const;
  Transformation inv_transformation() const;

  // The matrix elements depend on the details of the derived Gate class.
  virtual std::tuple<cmplx, cmplx, cmplx, cmplx> matrix_elements() const = 0;

protected:
  int target_;
  Controls controls_;
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
  double parameter(int paramNum) const override;  // (Note that we also have an accessor function: angle().)
  void setParameter(int paramNum, double value) override;
  void randomizeParameters() override;

  void random() override;
  void mutate() override;

  double angle() const;

  bool mutatable() const override;

  State applyGradTo(const State& state, int paramNum) const override;

  void output(std::ostream& out) const override;

private:
  bool equals_(const Gate& rhs) const override;
  bool sort_before(const Gate& rhs) const override;

  virtual Transformation grad_transformation() const;  // Transformation using the gradient matrix

  virtual std::tuple<cmplx, cmplx, cmplx, cmplx> grad_matrix_elements() const = 0;  // Derivative of the...
                                                                                    // ...transformation matrix.
protected:
  double angle_;
};


class DiagonalGate : public virtual SingleTargetGate
{
  // A class to handle all the gates that have diagonal matrices, i.e. PhaseGate, PhaseInv, PiByEight, PiByEightInv,
  // ZGate, ZRotation and ArbitraryPhase. This is provided merely to simplify the simplification code, where such code
  // considers the commutation of gates. Note that this is NOT the set of gates where target and control bits may be
  // interchanged - the ZRotation gate does not allow this.
public:
  explicit DiagonalGate(const CircuitContext& context);
  DiagonalGate(const CircuitContext& context, int target, const Controls& controls);

  bool matricesCommute(const SingleTargetGate& other) const override;  // Unnecessary, but included for completeness.
  bool matricesCommute(const DiagonalGate& other) const override;  // Unnecessary, but included for completeness.

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
  PhaseTypeGate(const CircuitContext& context, int target, const Controls& controls);  // Switches target and first...
                                                                                       // ...control if necessary.
  PhaseTypeGate(const CircuitContext& context, const Controls& controls);  // First 'control' becomes the target.

  int numQbitOptions() const override;

  void swapBits(int i, int j) override;

  bool cancelsAtStart() const override;

protected:
  void calculate_option_id() override;
  void randomize_qbits() override;
};


class XTypeGate : public virtual SingleTargetGate
{
  // Handles all gates that handle X Rotations, up to a phase. This class has only been created in an attempt to reduce
  // the number of functions involved in implementing simplification (specifically swapping) of gate pairs via double
  // dispatch. The reason that these gate types have been collected together is that their 2x2 matrices commute, meaning
  // that such gates may be swapped more often standard SingleTargetGates.
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



class Hadamard : public SingleTargetGate
{
public:
  explicit Hadamard(const CircuitContext& context);
  Hadamard(const CircuitContext& context, int target, const Controls& controls);

  int gateTypeId() const override;  // (I prefer not to have 'if type == ' code, but haven't found a better way yet.

  std::unique_ptr<Gate> clone() const override;
  Hadamard* rawClone() const override;  // Used only in the simplification routines to avoid inefficient shared_ptrs.

  std::string name() const override;

  static bool gateAvailable(const CircuitContext& context, int target,  // Can't use 'available'.
                            const Controls& controls);
  static long gateCost(const CircuitContext& context, int numControls);


  // The remaining public member functions are those for simplification, first the 'first-dispatch' versions, then the
  // 'second-dispatch' versions grouped by gate type.
  std::pair<long, bool> canSimplify(const Gate& next) const override;  // Returns cost improvement and whether #dof...
  GateSequence simplification(const Gate& next) const override;        // ...increases.
  bool canHSwap(const Gate& next) const override;
  std::unique_ptr<Gate> rightHMoverChange(const Gate& next) const override;  // nullptr <=> no change.
  std::unique_ptr<Gate> leftHMoverChange(const Gate& next) const override;

  std::pair<long, bool> canSimplify(const Hadamard& prev) const override;
  GateSequence simplification(const Hadamard& prev) const override;

  std::pair<long, bool> canSimplify(const XGate& prev) const override;
  GateSequence simplification(const XGate& prev) const override;  // HX->ZH is a 'simplification' if cost changes.
  bool canHSwap(const XGate& prev) const override;  // HX->ZH is considered a HSwap if X and Z have the same cost.
  std::unique_ptr<Gate> rightHMoverChange(const XGate& prev) const override;
  std::unique_ptr<Gate> leftHMoverChange(const XGate& prev) const override;

  bool canHSwap(const YGate& prev) const override;
  std::unique_ptr<Gate> rightHMoverChange(const YGate& prev) const override;  // These always return nullptr...
  std::unique_ptr<Gate> leftHMoverChange(const YGate& prev) const override;   // ...indicating no change.

  std::pair<long, bool> canSimplify(const ZGate& prev) const override;
  GateSequence simplification(const ZGate& prev) const override;
  bool canHSwap(const ZGate& prev) const override;
  std::unique_ptr<Gate> rightHMoverChange(const ZGate& prev) const override;
  std::unique_ptr<Gate> leftHMoverChange(const ZGate& prev) const override;

  std::pair<long, bool> canSimplify(const XRotation& prev) const override;
  GateSequence simplification(const XRotation& prev) const override;
  bool canHSwap(const XRotation& prev) const override;
  std::unique_ptr<Gate> rightHMoverChange(const XRotation& prev) const override;
  std::unique_ptr<Gate> leftHMoverChange(const XRotation& prev) const override;

  bool canHSwap(const YRotation& prev) const override;
  std::unique_ptr<Gate> rightHMoverChange(const YRotation& prev) const override;
  std::unique_ptr<Gate> leftHMoverChange(const YRotation& prev) const override;

  std::pair<long, bool> canSimplify(const ZRotation& prev) const override;
  GateSequence simplification(const ZRotation& prev) const override;
  bool canHSwap(const ZRotation& prev) const override;
  std::unique_ptr<Gate> rightHMoverChange(const ZRotation& prev) const override;
  std::unique_ptr<Gate> leftHMoverChange(const ZRotation& prev) const override;

private:
  static const int gate_type_id = 0;
  std::tuple<cmplx, cmplx, cmplx, cmplx> matrix_elements() const override;
};


class PiByEight : public PhaseTypeGate
{
  // Otherwise called T
public:
  explicit PiByEight(const CircuitContext& context);
  PiByEight(const CircuitContext& context, int target, const Controls& controls);
  PiByEight(const CircuitContext& context, const Controls& controls);  // Target and controls are interchangeable.
                                                                       // First 'control' is selected to be target.
  int gateTypeId() const override;  // I prefer not to have 'if type == ' code, but I haven't found a better way yet.

  std::unique_ptr<Gate> clone() const override;
  PiByEight* rawClone() const override;  // Used only in the simplification routines to avoid inefficient shared_ptrs.

  std::string name() const override;

  static bool gateAvailable(const CircuitContext& context, int target, const Controls& controls);
  static bool gateAvailable(const CircuitContext& context, const Controls& controls);  // Target is any one of the...
  static long gateCost(const CircuitContext& context, int numControls);                // ...'controls'.


  // The remaining public member functions are overrides for the simplification routines. First are those for 'first
  // dispatch', then come those for 'second dispatch' grouped by gate type.
  std::pair<long, bool> canSimplify(const Gate& next) const override;  // Returns cost improvement and whether #dof...
  GateSequence simplification(const Gate& next) const override;        // ...increases.

  std::pair<long, bool> canSimplify(const PiByEight& prev) const override;
  GateSequence simplification(const PiByEight& prev) const override;

  std::pair<long, bool> canSimplify(const PiByEightInv& prev) const override;
  GateSequence simplification(const PiByEightInv& prev) const override;

  std::pair<long, bool> canSimplify(const PhaseGate& prev) const override;
  GateSequence simplification(const PhaseGate& prev) const override;

  std::pair<long, bool> canSimplify(const PhaseInv& prev) const override;
  GateSequence simplification(const PhaseInv& prev) const override;

  std::pair<long, bool> canSimplify(const ZGate& prev) const override;
  GateSequence simplification(const ZGate& prev) const override;

  std::pair<long, bool> canSimplify(const ZRotation& prev) const override;
  GateSequence simplification(const ZRotation& prev) const override;

  std::pair<long, bool> canSimplify(const ArbitraryPhase& prev) const override;
  GateSequence simplification(const ArbitraryPhase& prev) const override;

  std::pair<long, bool> canSimplify(const SU2Gate& prev) const override;
  GateSequence simplification(const SU2Gate& prev) const override;

private:
  static const int gate_type_id = 1;
  std::tuple<cmplx, cmplx, cmplx, cmplx> matrix_elements() const override;
};


class PiByEightInv : public PhaseTypeGate
{
  // Otherwise called TInv
public:
  explicit PiByEightInv(const CircuitContext& context);
  PiByEightInv(const CircuitContext& context, int target, const Controls& controls);
  PiByEightInv(const CircuitContext& context, const Controls& controls);  // Target and controls are interchangeable.
                                                                          // Select one to be 'the' target.
  int gateTypeId() const override;  // (I prefer not to have 'if type == ' code, but I haven't found a better way yet.)

  std::unique_ptr<Gate> clone() const override;
  PiByEightInv* rawClone() const override; // Used only in the simplification routines to avoid inefficient shared_ptrs.

  std::string name() const override;

  static bool gateAvailable(const CircuitContext& context, int target, const Controls& controls);
  static bool gateAvailable(const CircuitContext& context, const Controls& controls);  // Target is any one of the...
  static long gateCost(const CircuitContext& context, int numControls);                // ...'controls'.


  // The remaining public member functions are overrides for the simplification routines. First are those for 'first
  // dispatch', then come those for 'second dispatch' grouped by gate type.
  std::pair<long, bool> canSimplify(const Gate& next) const override;  // Returns cost improvement and whether #dof...
  GateSequence simplification(const Gate& next) const override;        // ...increases.

  std::pair<long, bool> canSimplify(const PiByEight& prev) const override;
  GateSequence simplification(const PiByEight& prev) const override;

  std::pair<long, bool> canSimplify(const PiByEightInv& prev) const override;
  GateSequence simplification(const PiByEightInv& prev) const override;

  std::pair<long, bool> canSimplify(const PhaseGate& prev) const override;
  GateSequence simplification(const PhaseGate& prev) const override;

  std::pair<long, bool> canSimplify(const PhaseInv& prev) const override;
  GateSequence simplification(const PhaseInv& prev) const override;

  std::pair<long, bool> canSimplify(const ZGate& prev) const override;
  GateSequence simplification(const ZGate& prev) const override;

  std::pair<long, bool> canSimplify(const ZRotation& prev) const override;
  GateSequence simplification(const ZRotation& prev) const override;

  std::pair<long, bool> canSimplify(const ArbitraryPhase& prev) const override;
  GateSequence simplification(const ArbitraryPhase& prev) const override;

  std::pair<long, bool> canSimplify(const SU2Gate& prev) const override;
  GateSequence simplification(const SU2Gate& prev) const override;

private:
  static const int gate_type_id = 2;
  std::tuple<cmplx, cmplx, cmplx, cmplx> matrix_elements() const override;
};


class PhaseGate : public PhaseTypeGate
{
  // Otherwise called S
public:
  explicit PhaseGate(const CircuitContext& context);
  PhaseGate(const CircuitContext& context, int target, const Controls& controls);
  PhaseGate(const CircuitContext& context, const Controls& controls);  // Target and controls are interchangeable.
                                                                       //Select one to be 'the' target.
  int gateTypeId() const override;  // (I prefer not to have 'if type == ' code, but haven't found a better way yet.)

  std::unique_ptr<Gate> clone() const override;
  PhaseGate* rawClone() const override;  // Used only in the simplification routines to avoid inefficient shared_ptrs.

  std::string name() const override;

  static bool gateAvailable(const CircuitContext& context, int target, const Controls& controls);
  static bool gateAvailable(const CircuitContext& context, const Controls& controls);  // Target is any one of the...
  static long gateCost(const CircuitContext& context, int numControls);                // ...'controls'.


  // The remaining public member functions are overrides for the simplification routines. First are those for 'first
  // dispatch', then come those for 'second dispatch' grouped by gate type.
  std::pair<long, bool> canSimplify(const Gate& next) const override;  // Returns cost improvement and whether #dof...
  GateSequence simplification(const Gate& next) const override;        // ...increases.
  bool canRSwap(const Gate& next) const override;
  std::unique_ptr<Gate> rightRMoverChange(const Gate& next) const override;  // These return a null pointer if there...
  std::unique_ptr<Gate> leftRMoverChange(const Gate& next) const override;   // ...is no change.

  std::pair<long, bool> canSimplify(const PiByEight& prev) const override;
  GateSequence simplification(const PiByEight& prev) const override;

  std::pair<long, bool> canSimplify(const PiByEightInv& prev) const override;
  GateSequence simplification(const PiByEightInv& prev) const override;

  std::pair<long, bool> canSimplify(const PhaseGate& prev) const override;
  GateSequence simplification(const PhaseGate& prev) const override;

  std::pair<long, bool> canSimplify(const PhaseInv& prev) const override;
  GateSequence simplification(const PhaseInv& prev) const override;

  std::pair<long, bool> canSimplify(const ZGate& prev) const override;
  GateSequence simplification(const ZGate& prev) const override;

  std::pair<long, bool> canSimplify(const XRotation& prev) const override;
  GateSequence simplification(const XRotation& prev) const override;
  bool canRSwap(const XRotation& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const XRotation& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const XRotation& prev) const override;

  std::pair<long, bool> canSimplify(const YRotation& prev) const override;
  GateSequence simplification(const YRotation& prev) const override;
  bool canRSwap(const YRotation& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const YRotation& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const YRotation& prev) const override;

  std::pair<long, bool> canSimplify(const ZRotation& prev) const override;
  GateSequence simplification(const ZRotation& prev) const override;

  std::pair<long, bool> canSimplify(const ArbitraryPhase& prev) const override;
  GateSequence simplification(const ArbitraryPhase& prev) const override;

  std::pair<long, bool> canSimplify(const SU2Gate& prev) const override;
  GateSequence simplification(const SU2Gate& prev) const override;

private:
  static const int gate_type_id = 3;
  std::tuple<cmplx, cmplx, cmplx, cmplx> matrix_elements() const override;
};


class PhaseInv : public PhaseTypeGate
{
  // Otherwise called SInv
public:
  explicit PhaseInv(const CircuitContext& context);
  PhaseInv(const CircuitContext& context, int target, const Controls& controls);
  PhaseInv(const CircuitContext& context, const Controls& controls);  // Target and controls are interchangeable.
                                                                      //Select one to be 'the' target.
  int gateTypeId() const override;  // (I prefer not to have 'if type == ' code, but haven't found a better way yet.).

  std::unique_ptr<Gate> clone() const override;
  PhaseInv* rawClone() const override;  // Used only in the simplification routines to avoid inefficient shared_ptrs.

  std::string name() const override;

  static bool gateAvailable(const CircuitContext& context, int target, const Controls& controls);
  static bool gateAvailable(const CircuitContext& context, const Controls& controls);  // Target is any one of the...
  static long gateCost(const CircuitContext& context, int numControls);                // ...controls.


  // The remaining public member functions are overrides for the simplification routines. First are those for 'first
  // dispatch', then come those for 'second dispatch' grouped by gate type.
  std::pair<long, bool> canSimplify(const Gate& next) const override;  // Returns cost improvement and whether #dof...
  GateSequence simplification(const Gate& next) const override;        // ...increases.
  bool canRSwap(const Gate& next) const override;  // Calls next.canRSwap(*this). (Double dispatch.)
  std::unique_ptr<Gate> rightRMoverChange(const Gate& next) const override;  // These returns a null pointer if there...
  std::unique_ptr<Gate> leftRMoverChange(const Gate& next) const override;  // ...is no change.

  std::pair<long, bool> canSimplify(const PiByEight& prev) const override;
  GateSequence simplification(const PiByEight& prev) const override;

  std::pair<long, bool> canSimplify(const PiByEightInv& prev) const override;
  GateSequence simplification(const PiByEightInv& prev) const override;

  std::pair<long, bool> canSimplify(const PhaseGate& prev) const override;
  GateSequence simplification(const PhaseGate& prev) const override;

  std::pair<long, bool> canSimplify(const PhaseInv& prev) const override;
  GateSequence simplification(const PhaseInv& prev) const override;

  std::pair<long, bool> canSimplify(const ZGate& prev) const override;
  GateSequence simplification(const ZGate& prev) const override;

  std::pair<long, bool> canSimplify(const XRotation& prev) const override;
  GateSequence simplification(const XRotation& prev) const override;
  bool canRSwap(const XRotation& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const XRotation& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const XRotation& prev) const override;

  std::pair<long, bool> canSimplify(const YRotation& prev) const override;
  GateSequence simplification(const YRotation& prev) const override;
  bool canRSwap(const YRotation& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const YRotation& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const YRotation& prev) const override;

  std::pair<long, bool> canSimplify(const ZRotation& prev) const override;
  GateSequence simplification(const ZRotation& prev) const override;

  std::pair<long, bool> canSimplify(const ArbitraryPhase& prev) const override;
  GateSequence simplification(const ArbitraryPhase& prev) const override;

  std::pair<long, bool> canSimplify(const SU2Gate& prev) const override;
  GateSequence simplification(const SU2Gate& prev) const override;

private:
  static const int gate_type_id = 4;
  std::tuple<cmplx, cmplx, cmplx, cmplx> matrix_elements() const override;
};


class XGate : public XTypeGate
{
public:
  explicit XGate(const CircuitContext& context);
  XGate(const CircuitContext& context, int target, const Controls& controls);

  int gateTypeId() const override;  // (I prefer not to have 'if type == ' code, but haven't found a better way yet.).

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
  std::pair<long, bool> canSimplify(const Gate& next) const override;  // Returns cost improvement and whether #dof...
  GateSequence simplification(const Gate& next) const override;        // ...increases
  bool canHSwap(const Gate& next) const override;
  std::unique_ptr<Gate> rightHMoverChange(const Gate& next) const override;  // These return a null pointer if there...
  std::unique_ptr<Gate> leftHMoverChange(const Gate& next) const override;   // ...is no change.
  bool canRSwap(const Gate& next) const override;
  std::unique_ptr<Gate> rightRMoverChange(const Gate& next) const override;
  std::unique_ptr<Gate> leftRMoverChange(const Gate& next) const override;

  std::pair<long, bool> canSimplify(const Hadamard& prev) const override;
  GateSequence simplification(const Hadamard& prev) const override;
  bool canHSwap(const Hadamard& prev) const override;
  std::unique_ptr<Gate> rightHMoverChange(const Hadamard& prev) const override;
  std::unique_ptr<Gate> leftHMoverChange(const Hadamard& prev) const override;

  std::pair<long, bool> canSimplify(const XGate& prev) const override;
  GateSequence simplification(const XGate& prev) const override;

  std::pair<long, bool> canSimplify(const YGate& prev) const override;
  GateSequence simplification(const YGate& prev) const override;

  std::pair<long, bool> canSimplify(const ZGate& prev) const override;
  GateSequence simplification(const ZGate& prev) const override;

  std::pair<long, bool> canSimplify(const XRotation& prev) const override;
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

  std::pair<long, bool> canSimplify(const SU2Gate& prev) const override;
  GateSequence simplification(const SU2Gate& prev) const override;

private:
  static const int gate_type_id = 5;
  std::tuple<cmplx, cmplx, cmplx, cmplx> matrix_elements() const override;
};


class YGate : public YTypeGate
{
public:
  explicit YGate(const CircuitContext& context);
  YGate(const CircuitContext& context, int target, const Controls& controls);

  int gateTypeId() const override;  // (I prefer not to have 'if type == ' code, but haven't found a better way yet.).

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
  std::pair<long, bool> canSimplify(const Gate& next) const override;  // Returns cost improvement and whether #dof...
  GateSequence simplification(const Gate& next) const override;        // ...increases.
  bool canHSwap(const Gate& next) const override;  // Calls next.canHSwap(*this). (Double dispatch.)
  std::unique_ptr<Gate> rightHMoverChange(const Gate& next) const override;  // These return a null pointer if there...
  std::unique_ptr<Gate> leftHMoverChange(const Gate& next) const override;  // ...is no change.
  bool canRSwap(const Gate& next) const override;  // Calls next.canRSwap(*this). (Double dispatch.)
  std::unique_ptr<Gate> rightRMoverChange(const Gate& next) const override;
  std::unique_ptr<Gate> leftRMoverChange(const Gate& next) const override;

  bool canHSwap(const Hadamard& prev) const override;
  std::unique_ptr<Gate> rightHMoverChange(const Hadamard& prev) const override;  // These always return a null...
  std::unique_ptr<Gate> leftHMoverChange(const Hadamard& prev) const override;  // ...pointer indicating no change.

  std::pair<long, bool> canSimplify(const XGate& prev) const override;
  GateSequence simplification(const XGate& prev) const override;

  std::pair<long, bool> canSimplify(const YGate& prev) const override;
  GateSequence simplification(const YGate& prev) const override;

  std::pair<long, bool> canSimplify(const ZGate& prev) const override;
  GateSequence simplification(const ZGate& prev) const override;

  bool canRSwap(const XRotation& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const XRotation& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const XRotation& prev) const override;

  std::pair<long, bool> canSimplify(const YRotation& prev) const override;
  GateSequence simplification(const YRotation& prev) const override;

  bool canRSwap(const ZRotation& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const ZRotation& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const ZRotation& prev) const override;

  bool canRSwap(const ArbitraryPhase& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const ArbitraryPhase& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const ArbitraryPhase& prev) const override;

  std::pair<long, bool> canSimplify(const SU2Gate& prev) const override;
  GateSequence simplification(const SU2Gate& prev) const override;

private:
  static const int gate_type_id = 6;
  std::tuple<cmplx, cmplx, cmplx, cmplx> matrix_elements() const override;
};


class ZGate : public PhaseTypeGate
{
public:
  explicit ZGate(const CircuitContext& context);
  ZGate(const CircuitContext& context, int target, const Controls& controls);
  ZGate(const CircuitContext& context, const Controls& controls);  // Target and controls are interchangeable.
                                                                   // Select one to be 'the' target.
  int gateTypeId() const override;  // (I prefer not to have 'if type == ' code, but haven't found a better way yet.).

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
  std::pair<long, bool> canSimplify(const Gate& next) const override;  // Returns cost improvement and whether #dof...
  GateSequence simplification(const Gate& next) const override;        // ...increases.
  bool canHSwap(const Gate& next) const override;  // Calls next.canHSwap(*this). (Double dispatch.)
  std::unique_ptr<Gate> rightHMoverChange(const Gate& next) const override;  // These return a null pointer if there...
  std::unique_ptr<Gate> leftHMoverChange(const Gate& next) const override;   // ...is no change.
  bool canRSwap(const Gate& next) const override;  // Calls next.canRSwap(*this). (Double dispatch.)
  std::unique_ptr<Gate> rightRMoverChange(const Gate& next) const override;
  std::unique_ptr<Gate> leftRMoverChange(const Gate& next) const override;

  std::pair<long, bool> canSimplify(const Hadamard& prev) const override;
  GateSequence simplification(const Hadamard& prev) const override;
  bool canHSwap(const Hadamard& prev) const override;
  std::unique_ptr<Gate> rightHMoverChange(const Hadamard& prev) const override;
  std::unique_ptr<Gate> leftHMoverChange(const Hadamard& prev) const override;

  std::pair<long, bool> canSimplify(const PiByEight& prev) const override;
  GateSequence simplification(const PiByEight& prev) const override;

  std::pair<long, bool> canSimplify(const PiByEightInv& prev) const override;
  GateSequence simplification(const PiByEightInv& prev) const override;

  std::pair<long, bool> canSimplify(const PhaseGate& prev) const override;
  GateSequence simplification(const PhaseGate& prev) const override;

  std::pair<long, bool> canSimplify(const PhaseInv& prev) const override;
  GateSequence simplification(const PhaseInv& prev) const override;

  std::pair<long, bool> canSimplify(const XGate& prev) const override;
  GateSequence simplification(const XGate& prev) const override;

  std::pair<long, bool> canSimplify(const YGate& prev) const override;
  GateSequence simplification(const YGate& prev) const override;

  std::pair<long, bool> canSimplify(const ZGate& prev) const override;
  GateSequence simplification(const ZGate& prev) const override;

  std::pair<long, bool> canSimplify(const XRotation& prev) const override;
  GateSequence simplification(const XRotation& prev) const override;
  bool canRSwap(const XRotation& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const XRotation& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const XRotation& prev) const override;

  std::pair<long, bool> canSimplify(const YRotation& prev) const override;
  GateSequence simplification(const YRotation& prev) const override;
  bool canRSwap(const YRotation& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const YRotation& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const YRotation& prev) const override;

  std::pair<long, bool> canSimplify(const ZRotation& prev) const override;
  GateSequence simplification(const ZRotation& prev) const override;

  std::pair<long, bool> canSimplify(const ArbitraryPhase& prev) const override;
  GateSequence simplification(const ArbitraryPhase& prev) const override;

  std::pair<long, bool> canSimplify(const SU2Gate& prev) const override;
  GateSequence simplification(const SU2Gate& prev) const override;

private:
  static const int gate_type_id = 7;
  std::tuple<cmplx, cmplx, cmplx, cmplx> matrix_elements() const override;
};


class XRotation : public RotationGate, public XTypeGate
{
public:
  explicit XRotation(const CircuitContext& context);
  XRotation(const CircuitContext& context, int target, const Controls& controls, double angle);

  int gateTypeId() const override;  // (I prefer not to have 'if type == ' code, but haven't found a better way yet.).

  std::unique_ptr<Gate> clone() const override;  // Can't return std::unique_ptr<XRotation>. (Would fail to override.)
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

  std::pair<long, bool> canSimplify(const Gate& next) const override;  // Returns cost improvement and whether #dof...
  GateSequence simplification(const Gate& next) const override;        // ...increases.
  bool canHSwap(const Gate& next) const override;  // Calls next.canHSwap(*this). (Double dispatch.)
  std::unique_ptr<Gate> rightHMoverChange(const Gate& next) const override;  // These returns a null pointer if there...
  std::unique_ptr<Gate> leftHMoverChange(const Gate& next) const override;  // ...is no change.
  bool canRSwap(const Gate& next) const override;  // Calls next.canRSwap(*this). (Double dispatch.)
  std::unique_ptr<Gate> rightRMoverChange(const Gate& next) const override;
  std::unique_ptr<Gate> leftRMoverChange(const Gate& next) const override;

  std::pair<long, bool> canSimplify(const Hadamard& prev) const override;
  GateSequence simplification(const Hadamard& prev) const override;
  bool canHSwap(const Hadamard& prev) const override;
  std::unique_ptr<Gate> rightHMoverChange(const Hadamard& prev) const override;
  std::unique_ptr<Gate> leftHMoverChange(const Hadamard& prev) const override;

  std::pair<long, bool> canSimplify(const PhaseGate& prev) const override;
  GateSequence simplification(const PhaseGate& prev) const override;
  bool canRSwap(const PhaseGate& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const PhaseGate& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const PhaseGate& prev) const override;

  std::pair<long, bool> canSimplify(const PhaseInv& prev) const override;
  GateSequence simplification(const PhaseInv& prev) const override;
  bool canRSwap(const PhaseInv& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const PhaseInv& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const PhaseInv& prev) const override;

  std::pair<long, bool> canSimplify(const XGate& prev) const override;
  GateSequence simplification(const XGate& prev) const override;

  bool canRSwap(const YGate& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const YGate& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const YGate& prev) const override;

  std::pair<long, bool> canSimplify(const ZGate& prev) const override;
  GateSequence simplification(const ZGate& prev) const override;
  bool canRSwap(const ZGate& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const ZGate& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const ZGate& prev) const override;

  std::pair<long, bool> canSimplify(const XRotation& prev) const override;
  GateSequence simplification(const XRotation& prev) const override;

  std::pair<long, bool> canSimplify(const SU2Gate& prev) const override;
  GateSequence simplification(const SU2Gate& prev) const override;

private:
  static const int gate_type_id = 8;
  std::tuple<cmplx, cmplx, cmplx, cmplx> matrix_elements() const override;
  std::tuple<cmplx, cmplx, cmplx, cmplx> grad_matrix_elements() const override;
};


class YRotation : public RotationGate, public YTypeGate
{
public:
  explicit YRotation(const CircuitContext& context);
  YRotation(const CircuitContext& context, int target, const Controls& controls, double angle);

  int gateTypeId() const override;  // (I prefer not to have 'if type == ' code, but haven't found a better way yet.).

  std::unique_ptr<Gate> clone() const override;  // Can't return std::unique_ptr<YRotation>. (Would fail to override.)
  YRotation* rawClone() const override;  // Used only in the simplification routines to avoid inefficient shared_ptrs.

  std::string name() const override;

  static bool gateAvailable(const CircuitContext& context, int target, const Controls& controls);
  static long gateCost(const CircuitContext& context, int numControls);


  // The remaining public member functions are overrides for the simplification routines. First are those for gate
  // 'reduction', i.e. the detection of special values of the angle parameter that allow the gate to be replaced by a
  // cheaper non-parameterized gate. Then come the 'first dispatch' simplification and swap functions. Finally, we have
  // the 'second dispatch' functions grouped by gate type.

  bool canReduce() const override;
  GateSequence reduction() const override;

  std::pair<long, bool> canSimplify(const Gate& next) const override;  // Returns cost improvement and whether #dof...
  GateSequence simplification(const Gate& next) const override;        // ...increases.
  bool canHSwap(const Gate& next) const override;  // Calls next.canHSwap(*this). (Double dispatch.)
  std::unique_ptr<Gate> rightHMoverChange(const Gate& next) const override;  // These returns a null pointer if there...
  std::unique_ptr<Gate> leftHMoverChange(const Gate& next) const override;  // ...is no change.
  bool canRSwap(const Gate& next) const override;  // Calls next.canRSwap(*this). (Double dispatch.)
  std::unique_ptr<Gate> rightRMoverChange(const Gate& next) const override;
  std::unique_ptr<Gate> leftRMoverChange(const Gate& next) const override;

  bool canHSwap(const Hadamard& prev) const override;
  std::unique_ptr<Gate> rightHMoverChange(const Hadamard& prev) const override;
  std::unique_ptr<Gate> leftHMoverChange(const Hadamard& prev) const override;

  std::pair<long, bool> canSimplify(const PhaseGate& prev) const override;
  GateSequence simplification(const PhaseGate& prev) const override;
  bool canRSwap(const PhaseGate& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const PhaseGate& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const PhaseGate& prev) const override;

  std::pair<long, bool> canSimplify(const PhaseInv& prev) const override;
  GateSequence simplification(const PhaseInv& prev) const override;
  bool canRSwap(const PhaseInv& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const PhaseInv& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const PhaseInv& prev) const override;

  bool canRSwap(const XGate& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const XGate& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const XGate& prev) const override;

  std::pair<long, bool> canSimplify(const YGate& prev) const override;
  GateSequence simplification(const YGate& prev) const override;

  std::pair<long, bool> canSimplify(const ZGate& prev) const override;
  GateSequence simplification(const ZGate& prev) const override;
  bool canRSwap(const ZGate& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const ZGate& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const ZGate& prev) const override;

  std::pair<long, bool> canSimplify(const YRotation& prev) const override;
  GateSequence simplification(const YRotation& prev) const override;

  std::pair<long, bool> canSimplify(const SU2Gate& prev) const override;
  GateSequence simplification(const SU2Gate& prev) const override;

private:
  static const int gate_type_id = 9;
  std::tuple<cmplx, cmplx, cmplx, cmplx> matrix_elements() const override;
  std::tuple<cmplx, cmplx, cmplx, cmplx> grad_matrix_elements() const override;
};


class ZRotation : public DiagonalGate, public RotationGate
{
public:
  explicit ZRotation(const CircuitContext& context);
  ZRotation(const CircuitContext& context, int target, const Controls& controls, double angle);

  int gateTypeId() const override;  // (I prefer not to have 'if type == ' code, but haven't found a better way yet.).

  std::unique_ptr<Gate> clone() const override;  // Can't return std::unique_ptr<ZRotation>. (Would fail to override.)
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

  std::pair<long, bool> canSimplify(const Gate& next) const override;  // Returns cost improvement and whether #dof...
  GateSequence simplification(const Gate& next) const override;        // ...increases.
  bool canHSwap(const Gate& next) const override;  // Calls next.canHSwap(*this). (Double dispatch.)
  std::unique_ptr<Gate> rightHMoverChange(const Gate& next) const override;  // These return a null pointer if there...
  std::unique_ptr<Gate> leftHMoverChange(const Gate& next) const override;  // ...is no change.
  bool canRSwap(const Gate& next) const override;  // Calls next.canRSwap(*this). (Double dispatch.)
  std::unique_ptr<Gate> rightRMoverChange(const Gate& next) const override;
  std::unique_ptr<Gate> leftRMoverChange(const Gate& next) const override;

  std::pair<long, bool> canSimplify(const Hadamard& prev) const override;
  GateSequence simplification(const Hadamard& prev) const override;
  bool canHSwap(const Hadamard& prev) const override;
  std::unique_ptr<Gate> rightHMoverChange(const Hadamard& prev) const override;
  std::unique_ptr<Gate> leftHMoverChange(const Hadamard& prev) const override;

  std::pair<long, bool> canSimplify(const PiByEight& prev) const override;
  GateSequence simplification(const PiByEight& prev) const override;

  std::pair<long, bool> canSimplify(const PiByEightInv& prev) const override;
  GateSequence simplification(const PiByEightInv& prev) const override;

  std::pair<long, bool> canSimplify(const PhaseGate& prev) const override;
  GateSequence simplification(const PhaseGate& prev) const override;

  std::pair<long, bool> canSimplify(const PhaseInv& prev) const override;
  GateSequence simplification(const PhaseInv& prev) const override;

  bool canRSwap(const XGate& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const XGate& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const XGate& prev) const override;

  bool canRSwap(const YGate& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const YGate& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const YGate& prev) const override;

  std::pair<long, bool> canSimplify(const ZGate& prev) const override;
  GateSequence simplification(const ZGate& prev) const override;

  std::pair<long, bool> canSimplify(const ZRotation& prev) const override;
  GateSequence simplification(const ZRotation& prev) const override;

  std::pair<long, bool> canSimplify(const ArbitraryPhase& prev) const override;
  GateSequence simplification(const ArbitraryPhase& prev) const override;

  std::pair<long, bool> canSimplify(const SU2Gate& prev) const override;
  GateSequence simplification(const SU2Gate& prev) const override;

private:
  static const int gate_type_id = 10;
  std::tuple<cmplx, cmplx, cmplx, cmplx> matrix_elements() const override;
  std::tuple<cmplx, cmplx, cmplx, cmplx> grad_matrix_elements() const override;
};


class ArbitraryPhase : public PhaseTypeGate, public RotationGate
{
public:
  explicit ArbitraryPhase(const CircuitContext& context);
  ArbitraryPhase(const CircuitContext& context, int target, const Controls& controls, double angle);
  ArbitraryPhase(const CircuitContext& context, const Controls& controls, double angle);  // Target and controls are...
                                                                   // ...interchangeable. Select one to be 'the' target.
  int gateTypeId() const override;  // (I prefer not to have 'if type == ' code, but haven't found a better way yet.).

  std::unique_ptr<Gate> clone() const override;  // Can't return std::unique_ptr<ArbitraryPhase>.
  ArbitraryPhase* rawClone() const override;  // Used only in the simplification routines to avoid inefficient...
                                              // ...shared_ptrs.
  void setParameter(int paramNum, double value) override;  // Controlled ArbitraryPhase has a narrower angle range...
                                                           // ...than other RotationGates.
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

  std::pair<long, bool> canSimplify(const Gate& next) const override;  // Returns cost improvement and whether #dof...
  GateSequence simplification(const Gate& next) const override;        // ...increases.
  bool canRSwap(const Gate& next) const override;  // Calls next.canRSwap(*this). (Double dispatch.)
  std::unique_ptr<Gate> rightRMoverChange(const Gate& next) const override;  // These return a null pointer if there...
  std::unique_ptr<Gate> leftRMoverChange(const Gate& next) const override;  // ...is no change.

  std::pair<long, bool> canSimplify(const PiByEight& prev) const override;
  GateSequence simplification(const PiByEight& prev) const override;

  std::pair<long, bool> canSimplify(const PiByEightInv& prev) const override;
  GateSequence simplification(const PiByEightInv& prev) const override;

  std::pair<long, bool> canSimplify(const PhaseGate& prev) const override;
  GateSequence simplification(const PhaseGate& prev) const override;

  std::pair<long, bool> canSimplify(const PhaseInv& prev) const override;
  GateSequence simplification(const PhaseInv& prev) const override;

  bool canRSwap(const XGate& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const XGate& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const XGate& prev) const override;

  bool canRSwap(const YGate& prev) const override;
  std::unique_ptr<Gate> rightRMoverChange(const YGate& prev) const override;
  std::unique_ptr<Gate> leftRMoverChange(const YGate& prev) const override;

  std::pair<long, bool> canSimplify(const ZGate& prev) const override;
  GateSequence simplification(const ZGate& prev) const override;

  std::pair<long, bool> canSimplify(const ZRotation& prev) const override;
  GateSequence simplification(const ZRotation& prev) const override;

  std::pair<long, bool> canSimplify(const ArbitraryPhase& prev) const override;
  GateSequence simplification(const ArbitraryPhase& prev) const override;

  std::pair<long, bool> canSimplify(const SU2Gate& prev) const override;
  GateSequence simplification(const SU2Gate& prev) const override;

private:
  static const int gate_type_id = 11;
  std::tuple<cmplx, cmplx, cmplx, cmplx> matrix_elements() const override;
  std::tuple<cmplx, cmplx, cmplx, cmplx> grad_matrix_elements() const override;
};


// SU2Gates are currently not fully implemented - do not use!

class SU2Gate : public SingleTargetGate
{
  // Class for general rotations in SU2, requiring three rotation angles. (Each between -2pi and 2pi?)
public:
  explicit SU2Gate(const CircuitContext& context);
  SU2Gate(const CircuitContext& context, int target, const Controls& controls, double xAngle, double yAngle,
          double zAngle);

  int gateTypeId() const override;  // (I prefer not to have 'if type == ' code, but haven't found a better way yet.).

  std::unique_ptr<Gate> clone() const override;    // Can't return std::unique_ptr<SU2Gate>. (Would fail to override.)
  SU2Gate* rawClone() const override;  // Used only in the simplification routines to avoid inefficient shared_ptrs.

  int numParameters() const override;
  double parameter(int paramNum) const override;  // (Note that we already have an accessor functions: xAngle(), etc.)
  void setParameter(int paramNum, double value) override;
  void randomizeParameters() override;

  void random() override;
  void mutate() override;

  double xAngle() const;
  double yAngle() const;
  double zAngle() const;

  bool mutatable() const override;

  State applyGradTo(const State& state, int paramNum) const override;

  std::string name() const override;

  static bool gateAvailable(const CircuitContext& context, int target, const Controls& controls);
  static long gateCost(const CircuitContext& context, int numControls);

  void output(std::ostream& out) const override;


  // The remaining public member functions are overrides for the simplification routines. First are those for gate
  // 'reduction', i.e. the detection of special values of the angle parameter that allow the gate to be replaced by a
  // cheaper non-parameterized gate. Then come the 'first dispatch' simplification and swap functions. Finally, we have
  // the 'second dispatch' functions grouped by gate type.

  bool canReduce() const override;
  GateSequence reduction() const override;

  std::pair<long, bool> canSimplify(const Gate& next) const override;  // Returns cost improvement and whether #dof...
  GateSequence simplification(const Gate& next) const override;        // ...increases.

  std::pair<long, bool> canSimplify(const PiByEight& prev) const override;
  GateSequence simplification(const PiByEight& prev) const override;

  std::pair<long, bool> canSimplify(const PiByEightInv& prev) const override;
  GateSequence simplification(const PiByEightInv& prev) const override;

  std::pair<long, bool> canSimplify(const PhaseGate& prev) const override;
  GateSequence simplification(const PhaseGate& prev) const override;

  std::pair<long, bool> canSimplify(const PhaseInv& prev) const override;
  GateSequence simplification(const PhaseInv& prev) const override;

  std::pair<long, bool> canSimplify(const XGate& prev) const override;
  GateSequence simplification(const XGate& prev) const override;

  std::pair<long, bool> canSimplify(const YGate& prev) const override;
  GateSequence simplification(const YGate& prev) const override;

  std::pair<long, bool> canSimplify(const ZGate& prev) const override;
  GateSequence simplification(const ZGate& prev) const override;

  std::pair<long, bool> canSimplify(const XRotation& prev) const override;
  GateSequence simplification(const XRotation& prev) const override;

  std::pair<long, bool> canSimplify(const YRotation& prev) const override;
  GateSequence simplification(const YRotation& prev) const override;

  std::pair<long, bool> canSimplify(const ZRotation& prev) const override;
  GateSequence simplification(const ZRotation& prev) const override;

  std::pair<long, bool> canSimplify(const ArbitraryPhase& prev) const override;
  GateSequence simplification(const ArbitraryPhase& prev) const override;

  std::pair<long, bool> canSimplify(const SU2Gate& prev) const override;
  GateSequence simplification(const SU2Gate& prev) const override;

private:
  bool equals_(const Gate& rhs) const override;
  bool sort_before(const Gate& rhs) const override;

  Transformation grad_transformation(int paramNum) const;  // Transformation using the gradient matrix

  // While there is no need to separate out these functions, since we are not going to be inheriting from SU2Gate, I'll
  // keep them for consistency - and just in case!
  std::tuple<cmplx, cmplx, cmplx, cmplx> matrix_elements() const override;
  std::tuple<cmplx, cmplx, cmplx, cmplx> grad_matrix_elements(int paramNum) const;

protected:
  static const int gate_type_id = 12;
  double x_angle;
  double y_angle;
  double z_angle;
};


class SwapGate : public Gate
{
  // Class for the simple swap gate.
public:
  explicit SwapGate(const CircuitContext& context);
  SwapGate(const CircuitContext& context, int bit1, int bit2);  // Provided for testing.

  int gateTypeId() const override;  // (I prefer not to have 'if type == ' code, but haven't found a better way yet.).
  int numQbitOptions() const override;

  std::unique_ptr<Gate> clone() const override;
  SwapGate* rawClone() const override;  // Used only in the simplification routines to avoid inefficient shared_ptrs.

  void random() override;
  void swapBits(int i, int j) override;

  int bit1() const;  // Access required when commuting this gate with others.
  int bit2() const;  // Access required when commuting this gate with others.

  bool available() const override;

  State applyTo(const State& state) const override;
  State applyInvTo(const State& state) const override;

  int target() const override;  // Shouldn't be called - target() is only called when setting control bits and swap...
                                // ...gates have none (at present).
  std::string name() const override;

  static bool gateAvailable(const CircuitContext& context, int bit1, int bit2);
  static long gateCost(const CircuitContext& context);

  void output(std::ostream& out) const override;


  // The remaining public member functions are overrides for the simplification routines. First are those for 'first
  // dispatch', then come those for 'second dispatch' grouped by gate type.
  bool cancelsAtStart() const override;
  std::pair<long, bool> canSimplify(const Gate& next) const override;  // Returns cost improvement and whether #dof...
  GateSequence simplification(const Gate& next) const override;        // ...increases.
  bool canSwapSwap(const Gate& next) const override;
  std::unique_ptr<Gate> rightSwapMoverChange(const Gate& next) const override;  // These return a null pointer if...
  std::unique_ptr<Gate> leftSwapMoverChange(const Gate& next) const override;   // ...there is no change.

  bool canSwapSwap(const SingleTargetGate& prev) const override;
  std::unique_ptr<Gate> rightSwapMoverChange(const SingleTargetGate& prev) const override;
  std::unique_ptr<Gate> leftSwapMoverChange(const SingleTargetGate& prev) const override;

  std::pair<long, bool> canSimplify(const SwapGate& prev) const override;
  GateSequence simplification(const SwapGate& prev) const override;
  bool canSwapSwap(const SwapGate& prev) const override;
  std::unique_ptr<Gate> rightSwapMoverChange(const SwapGate& prev) const override;
  std::unique_ptr<Gate> leftSwapMoverChange(const SwapGate& prev) const override;

private:
  bool equivalent_structure(const Gate& rhs) const override;
  bool equals_(const Gate& rhs) const override;
  bool sort_before(const Gate& rhs) const override;
  void calculate_option_id() override;
  void fix_order();

private:
  static const int gate_type_id = 13;

  // Keep swap_bit1 < swap_bit2.
  int swap_bit1;
  int swap_bit2;
};


class Oracle : public Gate
{
  // Class for the Oracle, inheriting from Gate as for all other Gate classes.
  // (Was orginally in grover.h, but the Gate base class needs to know about the Oracle gate anyway, in order to perform
  // the double dispatch involved in circuit simplification.)
public:
  explicit Oracle(const CircuitContext& context);

  void imprint(const Circuit& circuit) override;

  int gateTypeId() const override;  // (I prefer not to have 'if type == ' code, but haven't found a better way yet.).
  int numQbitOptions() const override;

  std::unique_ptr<Gate> clone() const override;
  Oracle* rawClone() const override;  // Used only in the simplification routines to avoid inefficient shared_ptrs.

  void random() override;
  void swapBits(int i, int j) override;

  bool available() const override;

  State applyTo(const State& state) const override;
  State applyInvTo(const State& state) const override;

  int target() const override;  // Shouldn't be called - target() is only called when setting control bits and Oracle
                                // gate can have no controls.
  std::string name() const override;

  static bool gateAvailable(const CircuitContext& context);
  static long gateCost(const CircuitContext& context);  // (I would like to get rid of the 'context' parameter.)


  // The remaining public member functions are overrides for the simplification routines. First are those for 'first
  // dispatch', then come those for 'second dispatch' grouped by gate type.
  bool cancelsAtStart() const override;
  bool canSimplySwap(const Gate& next) const override;
  std::pair<long, bool> canSimplify(const Gate& next) const override;  // Returns cost improvement and whether #dof...
  GateSequence simplification(const Gate& next) const override;        // ...increases.

  bool canSimplySwap(const DiagonalGate& prev) const override;
  std::pair<long, bool> canSimplify(const Oracle& prev) const override;
  GateSequence simplification(const Oracle& prev) const override;

private:
  bool equivalent_structure(const Gate& rhs) const override;
  bool equals_(const Gate& rhs) const override;
  bool sort_before(const Gate& rhs) const override;
  void calculate_option_id() override;

private:
  static const int gate_type_id = 14;
  const int* marked_state;
};

//----------------------------------------------------------------------------------------------------------------------

// We create a master vector of 'GateCreators' - functions that create a Gate and return a unique pointer to it - and a
// map from Gate names to the location of the Gate type in the master vector. Given a list of Gate names, the Problem
// class can then create its own list of  GateCreators, allowing the creation of a random gate of arbitrary type,
// selected from this list. (I have left in some comments showing the code that used bare pointers (until wrapping in a
// unique pointer in create_selected_gate()), just in case I have problems compiling this code on other machines.)

// Template for converting all the Gate constructors into GateCreators.
template <class DerivedGate>
std::unique_ptr<Gate> createGate(const CircuitContext& context)
{
  return std::make_unique<DerivedGate>(context);
}

using GateCreator = std::unique_ptr<Gate> (*)(const CircuitContext& context);  // A function that takes a circuit...
                                                                      // ...context and produces a unique_ptr to a Gate.
// The master list of gate creators
extern const std::map<std::string, GateCreator> masterGateCreatorIndex;

//----------------------------------------------------------------------------------------------------------------------

// A helper function for GateSequences.
long gateCost(const GateSequence& gates);

#endif  // GATE_H
