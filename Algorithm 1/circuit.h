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

#ifndef CIRCUIT_H
#define CIRCUIT_H

#include <vector>
#include "state.h"
#include "controls.h"
#include "gate.h"
#include "front.h"

// At present, the 'problem' is split into the CircuitContext and the Problem classes, while the 'solution' is split
// into the Circuit and Solution classes. (In future, we may wish to consider the 'workspace' too.) The most obvious
// difference between these classes is that the Problem and Solution classes are designed to be derived from, e.g. by
// the Fourier and Grover versions. As such, they contain virtual functions. In contrast, the CircuitContext and Circuit
// classes are not derived from and do not contain virtual functions. The Circuit contains a pointer to the
// CircuitContext object, but no reference to the Problem object. We can consider the CircuitContext and Circuit classes
// as containing everything that is needed to manipulate and simplify circuits. Everything to do with evaluation is
// problem specific and hence is placed in the Problem and Solution classes. Things like the cache and the front, which
// perhaps ought to go in a separate Workspace class, are in the Problem class with references to them in the Solution
// class.

// Up to now (01/08/19) this separation has not caused any issues. However, I am now at a point where I am discovering
// that some aspects of circuit simplification are problem specific too. In particular, it is now useful to know
// something about the inputs on which the circuits are evaluated. For Grover's problem, the input qbits are always |0>,
// meaning that diagonal gates at the beginning of the circuit essentially do nothing and can be removed. For the
// Fourier problem this is not the case, but one might consider adding auxiliary qbits, in which case just some of the
// input qbits will be always |0>. Adding a function to the Problem classes that returns information about qbit inputs
// is straightforward - for Grover's we simply add a function that returns 'AlwaysOne' regardless of the qbit. However,
// since the Circuit class is where simplification takes place, we would need to add a reference to the Problem class -
// this breaks the neat separation of affairs described above. The alternative is to have the CircuitContext class store
// a vector of input values - alwaysOne, alwaysZero, varies - which is initialized by the Problem class. This just seems
// less pleasant than the simple (virtual) function in the Problem class. This issue is making me question whether this
// separation is useful.

// Circuit::simplify() will need to know about these classes from cache.h.
class CachedStats;
class CircuitCache;

enum class QbitOptions
{
  alwaysZero,
  alwaysOne,
  varies
};


class GateCost
{
  // Gate costs have a polynomial dependency on the number of controls. Polynomial coefficients are integers, to ensure
  // an integral overall circuit cost. (GateCost objects for gate types with no controls should have only a single
  // coefficient.)
public:
  GateCost();
  GateCost(const std::vector<long>& coeff);
  GateCost(std::vector<long>&& coeff);
  GateCost(std::initializer_list<long> coeff);

  long cost(int numControls) const;

  void output(std::ostream& out) const;
  
private:
  std::vector<long> coeff_;  // Starts with the constant term. Stops whenever the polynomial stops.
};

std::ostream& operator<<(std::ostream& out, const GateCost& gateCost);


namespace circuit
{
  void outputEvals(std::ostream& out);

  class CircuitContext
  {
    // The context in which Circuits are created. Enforces the number of qbits, the types of permitted gates, etc.
  public:
    CircuitContext(const std::vector<GateCreator>& gateCreator, const std::vector<PermittedControls>& permittedControls,
                   const std::vector<GateCost>& gateCosts, int numQbits, double meanRandomLength);

    int numPermittedGateOptions() const;

    // Two versions of functions for checking gate availability and cost. The first version, in each case, is provided
    // so that we don't need to create a gate, just to check availability or cost.
    bool gateTypeAvailable(int gateTypeId, int numControls) const;
    bool gateTypeAvailable(const Gate& gate, int numControls) const;
    long gateCost(int gateTypeId, int numControls) const;
    long gateCost(const Gate& gate, int numControls) const;
    size_t gateOptionBaseId(const Gate& gate) const;

    int numQbits() const;
    void setQbitInputOptions(int qbit, QbitOptions options);
    QbitOptions qbitInputOptions(int qbit) const;
    double meanRandomLength() const;

    std::unique_ptr<Gate> randomGate() const;

    void output(std::ostream& out) const;

  private:
    void make_available_with_many_controls(int gateTypeId);
    void remove_redundancies();
    std::unique_ptr<Gate> create_selected_gate(int gateTypeNum) const;

  private:
    int num_permitted_gate_types;  // Was 'const', but we now detect and remove redundant gate types.
    int num_permitted_gate_options;  // An 'option' is a combination of gate type, target bit and control bits.

    // The size of the GateCreator vector is just the number of permitted gate types...
    std::vector<GateCreator> gate_creator;  // Was 'const', but we now detect and remove redundant gate types.

    // ...while the following vectors contain an element for EVERY type of gate.
    std::vector<std::vector<bool> > gate_type_available;  // Indices are type_id and numControls.
    std::vector<GateCost> gate_cost;
    std::vector<size_t> gate_option_base_id;  // Each gate/qbit combo gets a unique id. This array hold the smallest...
                                              // ...id for each gate type.
    const int num_qbits;
    std::vector<QbitOptions> qbit_input_options;  // Indicates whether an input qbit is alwaysZero, alwaysOne or varies.
    const double mean_random_length;  // This feels like it belongs elsewhere.
  };

  std::ostream& operator<<(std::ostream& out, const CircuitContext& circuitContext);


  class Problem
  {
    // Acts as the base class for all quantum circuit Problem classes. (As it happens, there is unlikely to be much, if
    // any, difference between Problem classes.)
  public:
    Problem(const std::vector<GateCreator>& gateCreator, const std::vector<PermittedControls>& permittedControls,
            const std::vector<GateCost>& gateCosts, int numQbits, double meanRandomLength, double gateCostRange,
            int maxLength, double targetError, long expectedReduction, bool randomizeParameters, int iterationQuota,
            bool useCache);

    const CircuitContext& circuitContext() const;
    int numQbits() const;
    virtual QbitOptions qbitInputOptions(int qbit) const = 0;
    void transferQbitInputOptionsToContext();  // Must be called AFTER constructor, not within. (Calls virtual...
                                               // ...function qbitInputOptions().)
    // The following 4 functions are required by the MOMH framework.
    virtual int numObj() const = 0;
    virtual bool objMaximized(int objNum) const = 0;
    virtual double objLowerBound(int objNum) const = 0;
    virtual double objUpperBound(int objNum) const = 0;

    int maxLength() const;
    double targetError() const;
    long expectedReduction() const;
    bool randomizeParameters() const;
    int iterationQuota() const;

    std::shared_ptr<CircuitCache> cache() const;
    Front& front() const;

    virtual void output(std::ostream& out) const;  // This default version works fine with the problems I've...
                                                   // implemented thus far.
  protected:
    CircuitContext circuit_context;
    double gate_cost_range;  // A value used to scale circuit cost when calculating distances in objective space.
                             // (Artificial value.)
    // The following members are really algorithm parameters. We should think about calving off a class for handling the
    // (non-GA)algorithmic side of things.
    int max_length;
    double target_error;  // Value for overall error considered 'good enough'.
    long expected_reduction;  // Used to determine whether to evaluate a circuit when a 'good enough' circuit already...
                              // ...exists at lesser cost.
    bool randomize_parameters;  // Indicates whether gate parameters are randomized before numerical optimization.
    int iteration_quota;  // Numerical optimization proceeds only while expected to improve the front within this...
                          // ...quota of (additional) iterations.

    // The following members should perhaps be in a separate 'Workspace' class. One possibility would be to adapt the GA
    // code to require the user to provide a Workspace class in addition to a Problem and Solution class. Then those
    // items that, in a sense, have Problem scope, are not really part of the Problem, but are involved in enhancing the
    // evaluation of Solutions, could go in there. (Then we wouldn't need the pointer, but could just store the cache,
    // and we wouldn't need 'mutable'.)
    mutable std::shared_ptr<CircuitCache> cache_;
    mutable Front front_;
  };

  std::ostream& operator<<(std::ostream& out, const Problem& problem);


  class Circuit
  {
    friend void crossover(Circuit& lhs, Circuit& rhs);

  public:
    using CircuitContext = circuit::CircuitContext;  // Needed?

    explicit Circuit(const CircuitContext& circuitContext);
    Circuit(const CircuitContext& circuitContext, GateSequence&& gates);  // Created for testing only
    Circuit(const Circuit& rhs);
    Circuit(Circuit&& rhs);

    Circuit& operator=(const Circuit& rhs);  // Enables us to copy the unsimplified circuit back in to a Solution, if...
    Circuit& operator=(Circuit&& rhs);       // ...desired.

    bool equivalentStructure(const Circuit& rhs) const;
    bool operator==(const Circuit& rhs) const;
    bool operator!=(const Circuit& rhs) const;
    bool sortBefore(const Circuit& rhs) const;  // Enable sorting, efficient counting of duplicates, etc. (For...
                                                // ...population entropy calculations.)
    const int& markedState() const;
    void setMarkedState(int markedState);

    int length() const;
    int numParameters() const;
    const Gate& gate(int pos) const;
    std::vector<double> parameters() const;
    void setParameters(const std::vector<double>& values);

    void random();
    void randomizeParameters();

    // Various mutation operators. At present, this list contains only 'small' mutations. Each returns true if the
    // circuit is changed. Each function name is prefaced with 'mutate' to emphasize that these are mutation operators,
    // applied by the GA, that change the circuit. This contrasts them with operations such as those used by the
    // simplification routines that, for example, move a gate within the circuit via a sequence of swap operations that
    // leave the operation of the circuit unchanged.
    bool mutateReplaceGate();
    bool mutateInsertGate();  // How about sequences of gates?
    bool mutateRemoveGate();  // How about sequences of gates?
    bool mutateSwapGates();
    bool mutateMoveGate();
    bool mutateMutateGate();  // Included again to get 'adjust and reoptimize' functionality.

    State simulate(const State& startState) const;
    std::vector<State> simulateAndLog(const State& startState) const;  // Make private?
    std::vector<State> simulateInverseAndLog(const State& startState) const;  // Make private?
    cmplx simulatedOverlap(const State& startState, const State& targetState) const;
    std::pair<cmplx, std::vector<cmplx> > simulatedOverlapAndGrad(const State& startState,
                                                                  const State& targetState) const;

    // Primary simplification and reduction routines. Function makeCanonicalForm() uses the swaps permitted by the
    // various gate classes to arrange gates to be 'as sorted as possible' according to function sortBefore().
    // (Essentially sorts the gates according to sortBefore under the constraint that many gates cannot move past each
    // other.) Function simplifyStructure() simplifies the circuit structure, ignoring gate parameters and making using
    // of the cache to shortcut the process when possible. Function simplifyCircuit() simplies the circuit, taking care
    // to ensure that gate parameters match - shortcutting using the cache is not permitted. Both of these routines add
    // visited solutions to the cache and return a pointer to the (possibly new) statistics object that holds simplified
    // circuit, objective values, numerical optimization stats etc. Function reduce() is applied after numerical
    // optimization, detecting gates that, due to the special value of the angle parameter, can be replaced with cheaper
    // options or simply removed.
    void makeCanonicalForm(int startPos = 1);
    void makeUncanonical();
    std::shared_ptr<CachedStats> simplifyCircuit(std::shared_ptr<CircuitCache> cache,
                                                 std::shared_ptr<CachedStats> emptyStats);
    std::shared_ptr<CachedStats> simplifyStructure(std::shared_ptr<CircuitCache> cache,
                                                   std::shared_ptr<CachedStats> emptyStats);
    bool reduce();

    size_t hash() const;

    void output(std::ostream& out) const;

  private:
    void imprint_gates();  // Ensure gates know their 'parent circuit'.
    std::unique_ptr<Gate> random_gate() const;  // Wrapper for CircuitContext::randomGate(), that ensures the new...
                                                // ...gate knows its 'parent' circuit.
  private:
    // In the following functions, used in simplification, it is assumed that all SwapGates have been moved to the end
    // of the circuit, as should be the case in canonical form. HX->ZH is considered to be a 'swap' only if the XGate
    // and ZGate have the same cost.
    int leftmost_shift(int pos, bool involveSwapGates) const;
    int rightmost_shift(int pos, bool involveSwapGates) const;
    int leftmost_shift(int pos, std::vector<const Gate*>& movedGate, std::vector<const Gate*>& pointers,
                       bool involveSwapGates) const;
    int rightmost_shift(int pos, std::vector<const Gate*>& movedGate, std::vector<const Gate*>& pointers,
                        bool involveSwapGates) const;
    void swap_gate_right(int moverPos, int finalPos,
                         bool involveSwapGates);  // Move the gate by swapping with each intervening gate.
    void swap_gate_left(int moverPos, int finalPos,
                        bool involveSwapGates);  // Move the gate by swapping with each intervening gate.
    Replacement best_replacement() const;
    void apply_replacement(int leftPos, int rightPos, int meetPos);

    bool gates_available() const;  // Checks availability of all the gates - for testing only.

  private:
    const CircuitContext* context_;  // Owned by Problem class. Pointer since we need to be able to copy circuits.
    GateSequence gates_;
    int marked_state;  // Was a static member of Oracle. Moved here for parallelization.
  };

  std::ostream& operator<<(std::ostream& out, const Circuit& circuit);


  class Solution
  {
    // Acts as the base class for all quantum circuit Solution classes.
    // We plan to use the derived classes directly, i.e. to wrap fourier::Solution in an NSGA2Solution when we apply
    // NSGA2 to the Fourier problem. In other words, we do not plan on using polymorphism here - we have merely created
    // this base class to reduce duplication of code associated with circuit manipulation and cost.
    friend void crossover(Solution& lhs, Solution& rhs);
    friend class OverallErrorAndGradient;

  public:
    explicit Solution(const Problem& problem);
    // Default destructor and copy constructor are fine. No assignment operator - what would it do with problem_!?

    bool operator==(const Solution& rhs) const;  // To ensure only unique solutions get into the store.
    bool operator!=(const Solution& rhs) const;  // Might as well! (Probably not needed.)
    bool sortBefore(const Solution& rhs) const;  // To enable sorting, efficient counting of duplicates, etc.

    void random();
    void mutate(double mutateProb);
    int numParameters() const;
    std::vector<double> parameters() const;
    void setParameters(const std::vector<double>& values);

    void evaluateObjectives();  // Required by MOO code - includes circuit simplification, parameter optimization, etc.

    double objective(int objNum) const;  // Access functions required by the MOO code.
    long cost() const;  // Access cost, by name. (In other words, not by calling objective(2).)
    void output(std::ostream& out) const;
    void outputGradient(std::ostream& out) const;  // Gradient for overall (primary) error only.

  private:
    void evaluate_cost();  // Just evaluate the total gate cost of the circuit - the easy bit!
    virtual void assign_worst_error() = 0;  // E.g. for solutions we don't bother evaluating
    double primary_error() const;

    virtual void evaluate_error() = 0;  // The hard bit, involving multiple simulations of the circuit.
    virtual void evaluate_error_and_gradient() = 0;
    void optimize_parameters();
    void optimize_and_cache();  // If there are no gate parameters, simply evaluates and caches.

  protected:
    const Problem& problem_;
    Circuit circuit_;

    long total_cost = -1;
    std::vector<double> error_;
    std::vector<double> gradient_;  // Gradient of the primary error, i.e. of error_[0].
    bool evaluated_ = false;  // When, e.g., mutate() modifies the solution, we must be able to mark it unevaluated.

    std::shared_ptr<CircuitCache> cache_;  // (Doesn't 'own' the cache. Should this be a weak ptr? An ordinary ptr?)
    std::shared_ptr<CachedStats> stats_;
    Front& front_;
  };

  bool sortBefore(const Solution& lhs, const Solution& rhs);  // Enables sorting, efficient counting of duplicates, etc.
  std::ostream& operator<<(std::ostream& out, const Solution& solution);

}  // namespace circuit


#endif // CIRCUIT_H
