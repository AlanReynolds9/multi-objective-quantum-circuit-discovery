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

#ifndef CIRCUIT_H
#define CIRCUIT_H

#include <vector>
#include "state.h"
#include "controls.h"
#include "gate.h"
#include "population.h"

// At present, the 'problem' is split into the CircuitContext and the Problem classes, while the 'solution' is split
// into the Circuit and Solution classes. (In future, we may wish to consider the 'workspace' too.) The most obvious
// difference between these classes is that the  Problem and Solution classes are designed to be derived from, e.g. by
// the fourier and grover versions. As such, they contain virtual functions. In contrast, the CircuitContext and Circuit
// classes are not derived from and do not contain virtual functions. The Circuit contains a pointer to the
// CircuitContext object, but no reference to the Problem object. We can consider the CircuitContext and Circuit classes
// as containing everything that is needed to manipulate and simplify circuits. Everything to do with evaluation is
// problem specific and hence is placed in the Problem and Solution classes.
// Up to now (01/08/19) this separation has not caused any issues. However, I am now at a point where I am discovering
// that some aspects of circuit simplification are problem specific too. In particular, it is now useful to now
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

enum class QbitOptions
{
  alwaysZero,
  alwaysOne,
  varies
};


class GateCost
{
  // Gate costs have a polynomial dependency on the number of controls. Polynomial coefficients are integers, to ensure
  // an integral overall circuit cost.
  // (GateCost objects for gate types with no controls should have only a single coefficient.)
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
  class CircuitContext
  {
    // The context in which Circuits are created. Enforces the number of qbits, the types of permitted gates, etc.
  public:
    CircuitContext(const std::vector<GateCreator>& gateCreator, const std::vector<PermittedControls>& permittedControls,
                   const std::vector<GateCost>& gateCosts, int numQbits, double meanRandomLength,
                   double angleTolerance);

    bool gateTypeAvailable(int gateTypeId, int numControls) const;
    bool gateTypeAvailable(const Gate& gate, int numControls) const;
    long gateCost(int gateTypeId, int numControls) const;
    long gateCost(const Gate& gate, int numControls) const;

    int numQbits() const;
    void setQbitInputOptions(int qbit, QbitOptions options);
    QbitOptions qbitInputOptions(int qbit) const;
    double meanRandomLength() const;
    double angleTolerance() const;

    std::unique_ptr<Gate> randomGate() const;

    void output(std::ostream& out) const;

  private:
    void make_available_with_many_controls(int gateTypeId);
    std::unique_ptr<Gate> create_selected_gate(int gateTypeNum) const;

  private:
    int num_permitted_gate_types;  // Was 'const', but we now detect and remove redundant gate types.

    // The size of the GateCreator vector is just the number of permitted gate types...
    std::vector<GateCreator> gate_creator;  // Leave non-const for if we reintroduce removal of redundant gate types.

    // ...while the following vectors contain an element for EVERY type of gate.
    std::vector<std::vector<bool> > gate_type_available;  // Indices are type_id and numControls.
    std::vector<GateCost> gate_cost;

    const int num_qbits;
    std::vector<QbitOptions> qbit_input_options; // Indicates whether a qbit's input is alwaysZero, alwaysOne or varies.
    const double mean_random_length;  // (This feels like it belongs elsewhere.)
    const double angle_tolerance;
  };

  std::ostream& operator<<(std::ostream& out, const CircuitContext& circuitContext);


  class Circuit
  {
//    friend void crossover(Circuit& lhs, Circuit& rhs);

  public:
    using CircuitContext = circuit::CircuitContext;  // Needed?

    Circuit(const CircuitContext& circuitContext);
    Circuit(const CircuitContext& circuitContext, GateSequence&& gates);
    Circuit(const Circuit& rhs);
    Circuit(Circuit&& rhs);

    Circuit& operator=(const Circuit& rhs);  // Enable copying of the unsimplified circuit back into a Solution, if...
    Circuit& operator=(Circuit&& rhs);       // ...desired.

    void createGateSimulators();  // If any gate needs a simulator but doesn't have one, create and attach the...
                                  // ...simulator. Does not overwrite.
    bool operator==(const Circuit& rhs) const;
    bool operator!=(const Circuit& rhs) const;
    bool sortBefore(const Circuit& rhs) const;  // Enable sorting, efficient counting of duplicates, etc.

    void linkToMarkedState(std::shared_ptr<const int> markedState);
    std::shared_ptr<const int> markedState() const;

    int length() const;
    int numMutatable() const;  // Number of 'mutatable' gates.
    const Gate& gate(int pos) const;

    void random();

    State simulate(const State& startState) const;
    std::vector<State> simulateAndLog(const State& startState) const;
    std::vector<State> reverseSimulateAndLog(const State& targetState) const;

    // Primary simplification routines. Function makeCanonicalForm() uses the swaps permitted by the various gate
    // classes to arrange gates to be 'as sorted as possible' according to function sortBefore(). (Essentially sorts the
    // gates according to sortBefore under the constraint that many gates cannot move past each other.)
    void makeCanonicalForm(int startPos = 1);
    void makeUncanonical();
    void simplify();

    void output(std::ostream& out) const;

  private:
    void link_gates_to_marked_state();  // Ensure gates that need it get access to the marked state.
    std::unique_ptr<Gate> random_gate() const;  // Wrapper for CircuitContext::randomGate(), ensuring the new gate...
                                                // ...knows about the marked state if necessary.
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
    void swap_gate_right(int moverPos, int finalPos, bool involveSwapGates);  // These move the gate by swapping with...
    void swap_gate_left(int moverPos, int finalPos, bool involveSwapGates);   // ...each intervening gate.
    Replacement best_replacement() const;
    void apply_replacement(int leftPos, int rightPos, int meetPos);

    bool gates_available() const;  // Checks availability of all the gates - for testing only.

  private:
    const CircuitContext* context_;  // Owned by Problem class. Pointer since we need to be able to copy circuits.
    GateSequence gates_;
    std::shared_ptr<const int> marked_state;
  };

  std::ostream& operator<<(std::ostream& out, const Circuit& circuit);


  class Problem;


  class Insertion
  {
    // A gate sequence for inserting into solutions. If the sequence is long enough to make it worthwhile, it is
    // pre-simulated, by which we mean that the matrix for the gate sequence is calculated in advance.
    //
    // To keep things simple, we only permit one Oracle gate, at most, in the Insertion. Moreover, if an Oracle is
    // include, it should be the last gate. This allows us to apply the Insertion to a State with at most one matrix
    // multiplication and one application of an Oracle. Oracles in the middle of the Insertion, or multiple Oracles,
    // would require multiple matrices for the gate subsequences between Oracles.
  public:
    // Two constructors. The first moves data from a random externally created sequence. The second copies from a
    // subsequence of a solution.
    Insertion(const Problem& problem, GateSequence&& sequence);
    Insertion(const Problem& problem, GateSequence::const_iterator begin, GateSequence::const_iterator end);

    void createGateSimulators();  // Creates and attaches simulators for any gate that needs them.
    void prepare();  // Creates the simulator for the Insertion by taking each basis state and applying each gate in...
                     // ...turn. Called just once.
    int length() const;
    const Gate& gate(int pos) const;

    // Function that applies the Insertion (i.e. the matrix and optional Oracle) to an arbitrary State. Called often.
    void applyTo(State& state, std::shared_ptr<int> marked_state) const;

    void output(std::ostream& out) const;

  private:
    void handle_oracle();
    void simulate_to_oracle(State& state);

  private:
    const Problem& problem_;
    GateSequence sequence_;  // Includes the Oracle, if there is one.
    std::unique_ptr<SimulatorIndexManager> simulator_index_manager;
    std::unique_ptr<SparseGateSimulator> simulator_;
    bool includes_oracle = false;
  };

  std::ostream& operator<<(std::ostream& out, const Insertion& insertion);


  class Problem
  {
    // Acts as the base class for all quantum circuit Problem classes. (As it happens, there is unlikely to be much, if
    // any, difference between Problem classes.)
    // (This class is becoming polluted with things that are really algorithm parameters, such as the bad error cutoff
    // and the gate optimization probability. Consider refactoring?)
  public:
    Problem(const std::vector<GateCreator>& gateCreator, const std::vector<PermittedControls>& permittedControls,
            const std::vector<GateCost>& gateCosts, double gateCostRange, int numQbits, int maxLength,
            double meanRandomLength, double badErrorCutoff, double gateOptProb, int numInsertions,
            int minInsertionLength, int maxInsertionLength, double angleTolerance);

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
    double badErrorCutoff() const;
    double gateOptProb() const;
    int numInsertions() const;
    const Insertion& insertion(int i) const;
    void prepareForStep();

    virtual void output(std::ostream& out) const;  // This default version works fine with the problems I've...
                                                   // ...implemented thus far.
  protected:
    CircuitContext circuit_context;
    double gate_cost_range;  // Used to scale circuit cost when calculating distances in objective space. (Artificial.)
    int max_length;

    double bad_error_cutoff;  // We do not simplify, or even construct, a circuit if the overall error is found to be...
    double gate_opt_prob;     // ...worse than this.

    int num_insertions;
    int min_insertion_length;
    int max_insertion_length;
    std::vector<Insertion> insertions_;
  };

  std::ostream& operator<<(std::ostream& out, const Problem& problem);


  class Solution
  {
    // Acts as the base class for all quantum circuit Solution classes.
    // We plan to use the derived classes directly, i.e. to wrap fourier::Solution in an NSGA2Solution when we apply
    // NSGA2 to the Fourier problem. In other words, we do not plan on using polymorphism here - we have merely created
    // this base class to reduce duplication of code associated with circuit manipulation and cost.
    friend void crossover(Solution& lhs, Solution& rhs);

  public:
    explicit Solution(const Problem& problem);
    // Default destructor and copy constructor are fine. No assignment operator - what would it do with problem_!?
 
    bool operator==(const Solution& rhs) const;  // To ensure only unique solutions get into the store.
    bool operator!=(const Solution& rhs) const;  // Might as well! (Probably not needed.)
    bool sortBefore(const Solution& rhs) const;  // To enable sorting, efficient counting of duplicates, etc.

    void random();
    void mutate();

    void evaluateObjectives();  // Required by MOO code. This function determines for itself whether a quick or a...
                                // ...full evaluation is required.
    int length() const;
    bool bad() const;
    double objective(int objNum) const;  // Access functions required by the MOO code.
    long cost() const;  // Access cost, by name. (In other words, not by calling objective(0).)
    void output(std::ostream& out) const;

  protected:
    void clone_from(const Solution& parent);  // Copy what is necessary to become a clone of the parent.
    bool mutation_valid(int mutationType) const;

  private:
    void tweak_gate();
    void single_point_crossover();
    void cut_circuit();
    void remove_gate();
    void remove_sequence();
    void replicate_sequence();
    void insert_gate();
    void insert_sequence();

    void construct_circuit();

    void evaluate_cost();  // Just evaluate the total gate cost of the circuit - the easy bit!
    void assign_worst_cost();
    virtual void assign_worst_error() = 0;  // E.g. for solutions we don't bother evaluating.

    // Wrappers for the same functions in the Circuit, but store the intermediate states in this solution. (We assume
    // that Solution evaluation always consists of the same sequence of start states, target states and Oracle
    // settings.)
    void simulate_and_log(const State& startState);
    void reverse_simulate_and_log(const State& targetState);

    // Evaluation
    void record_intermediate_states();
    void calculate_overlaps();  // Used after intermediate states are been found, only on initial solutions.
    void quickly_evaluate();  // Also calculates the overlaps (naturally).
    void optimize_and_evaluate();

    // Some functions used by the evaluation routines, that will need to be defined in derived Solution classes.
    virtual int num_simulations() const = 0;  // The number of different input conditions requiring simulation.
                                              // (Typically the number of basis states.)
    virtual void prepare_(int sim) = 0;  // Prepare for the i'th input. (E.g. by setting the Oracles to mark the...
                                         // ...correct basis state.)
    virtual State start_state(int sim) const = 0;
    virtual State target_state(int sim) const = 0;
    virtual void calculate_error() = 0;
    virtual void calculate_symbolic() = 0;  // Used for gate parameter optimization
    

  protected:
    const Problem& problem_;

    // Initially, the decription of a child solution merely contains the information required to construct the solution
    // from parent solutions and a possible 'insertion'. This is sufficient to allow 'quick' evaluation.
    // (Weak pointers might be preferable to bare, but this requires access to a shared pointer.)
    // (Consider recording the cut points, rather than the number of gates from each parent.)

    // First data regarding the primary parent from which the solution is cloned, which provides the left of the circuit
    // if crossover is applied. This includes a pointer to the parent, a pointer to the circuit (in case the child
    // circuit is only constructed from the parents' circuits after the parents are dead) and the cut point.
    const Solution* left_parent = nullptr;
    std::shared_ptr<Circuit> left_circuit;
    int left_cut = 0;  // This solution will claim the gates to the left of the cut.

    // The same data for the other parent.
    const Solution* right_parent = nullptr;
    std::shared_ptr<Circuit> right_circuit;
    int right_cut = 0;  // This solution will claim the gates to the right of the cut.

    // Any new gate or gate sequence that might be inserted.
    std::unique_ptr<Gate> new_gate = nullptr;  // Single gate.
    int insertion_num = -1;  // Gate sequence. This is an index into the Problem's list of Insertions. Default value...
                             // of -1 indicates no sequence insertion.

    // A pointer to an integer specifying the marked state for the Oracle. We must have one per solution, to enable
    // parallel evaluation to work effectively. (Otherwise a thread would have to lock the marked state before each
    // simulation.) The circuit and contained Oracles will also have shared pointers to this.
    std::shared_ptr<int> marked_state;

    // The circuit. Solutions in the initial population start from this circuit. For child solutions, the circuit is
    // constructed after quick evaluation, provided it achieves reasonable error rates.
    std::shared_ptr<Circuit> circuit_;

    // Vectors of the intermediate states, used for the quick evaluation of child solutions. These are filled during a
    // 'full' evaluation. First index is the simulation number. Second gives the position in the circuit reached.
    std::vector<std::vector<State> > forward_intermediate_states;
    std::vector<std::vector<State> > reverse_intermediate_states;

    // Overlaps - both numeric for when we are not about to optimize a gate parameter and symbolic for when we are.
    // A 'symbolic overlap' is given as a vector of the coefficients for the constant, cosine and sine terms.
    std::vector<cmplx> overlaps_;
    std::vector<std::vector<cmplx> > symbolic_overlaps;

    bool bad_ = false;
    long total_cost = -1;
    std::vector<double> error_;
    bool evaluated_ = false;
    bool intermediate_states_recorded = false;

    std::vector<std::vector<double> > primary_error_coeffs;  // All the coefficients of the symbolic error, leaving a...
  };                                                         // ...single gate angle unfixed.

  bool sortBefore(const Solution& lhs, const Solution& rhs);  // Enables sorting, efficient counting of duplicates, etc.
  std::ostream& operator<<(std::ostream& out, const Solution& solution);

}  // namespace circuit


#endif // CIRCUIT_H
