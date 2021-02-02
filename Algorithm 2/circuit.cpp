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

#include <queue>
#include <unordered_set>
#include "circuit.h"
#include "replacements.h"  // Not needed at present - included in gate.h, which is included in circuit.h.
#include "simulator.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ostream;
using std::swap;
using std::make_pair;
using std::unique_ptr;
using std::make_unique;
using std::shared_ptr;
using std::make_shared;
using std::make_move_iterator;
using std::vector;
using utils::rand::randInt;
using utils::rand::rand01;
using utils::rand::randGeometric;


void deletePointers(vector<const Gate*>& pointers)
{
  for (auto pointer : pointers)
  {
    delete pointer;
  }
}


// GateCost class

GateCost::GateCost()
{
}


GateCost::GateCost(const vector<long>& coeff) :
coeff_(coeff)
{
}


GateCost::GateCost(vector<long>&& coeff) :
coeff_(std::move(coeff))
{
}


GateCost::GateCost(std::initializer_list<long> coeff) :
coeff_(coeff)
{
}


long GateCost::cost(int numControls) const
{
  long power{1};
  long total{0};
  for (auto c : coeff_)
  {
    total += c * power;
    power *= numControls;
  }

  return total;
}


void GateCost::output(std::ostream& out) const
{
  if (coeff_.empty())
  {
    out << "uninitialized" << endl;
    return;
  }
  out << coeff_[0];
  for (auto i = 1; i < coeff_.size(); ++i)
  {
    out << " + " << coeff_[i] << "c";
    if (i > 1)
    {
      out << "^i";
    }
  }
}


ostream& operator<<(ostream& out, const GateCost& gateCost)
{
  gateCost.output(out);
  return out;
}


namespace circuit
{
  CircuitContext::CircuitContext(const vector<GateCreator>& gateCreator,
                                 const vector<PermittedControls>& permittedControls, const vector<GateCost>& gateCosts,
                                 int numQbits, double meanRandomLength, double angleTolerance) :
  num_permitted_gate_types(static_cast<int>(gateCreator.size())),
  gate_creator(gateCreator),
  gate_type_available(constants::totalNumGateTypes, vector<bool>(numQbits, false)),
  gate_cost(constants::totalNumGateTypes),
  num_qbits(numQbits),
  qbit_input_options(numQbits, QbitOptions::varies),
  mean_random_length(meanRandomLength),
  angle_tolerance(angleTolerance)
  {
    assert(permittedControls.size() == num_permitted_gate_types);
    assert(gateCosts.size() == num_permitted_gate_types);

    for (auto i = 0; i < num_permitted_gate_types; ++i)
    {
      // Get a gate of the permitted type and the gate type id.
      auto gate = gate_creator[i](*this);  // Unpleasant. Gate created (with broken option ID) before CircuitContext...
      int id = gate->gateTypeId();         // ...is completely constructed.

      // Set gate costs
      gate_cost[id] = gateCosts[i];

      // Set gate availability. Index of gate_type_available[] is 0, 1 or 2 but means, none, one and many respectively.
      // Gates that don't have controls will look only at gate_type_available[0].
      switch (permittedControls[i])
      {
        case PermittedControls::none:
          gate_type_available[id][0] = true;
          break;

        case PermittedControls::one:
          gate_type_available[id][1] = true;
          break;

        case PermittedControls::notMany:
          gate_type_available[id][0] = gate_type_available[id][1] = true;
          break;

        case PermittedControls::many:
          make_available_with_many_controls(id);
          break;

        case PermittedControls::notOne:
          gate_type_available[id][0] = true;
          make_available_with_many_controls(id);
          break;

        case PermittedControls::atLeastOne:
          gate_type_available[id][1] = true;
          make_available_with_many_controls(id);
         break;

        case PermittedControls::any:
          gate_type_available[id][0] = gate_type_available[id][1] = true;
          make_available_with_many_controls(id);
          break;

        default:
          throw std::domain_error("Parameter 'permittedControls' in CircuitContext::CircuitContext contains an invalid"
                                  " element.");
      }
    }
  }


  bool CircuitContext::gateTypeAvailable(int gateTypeId, int numControls) const
  {
    // For use when we want to avoid creating a gate just to check availability. The Gate classes contain static
    // 'available' routines that call this function, with uncontrolled gate classes such as SwapGate calling it with
    // numControls set to zero.
    assert(numControls >= 0);
    assert(numControls < num_qbits);
    return gate_type_available[gateTypeId][numControls];
  }


  bool CircuitContext::gateTypeAvailable(const Gate& gate, int numControls) const
  {
    // The Gate classes contain their own 'available()' routines that call this function. Uncontrolled gate classes,
    // such as SwapGate, just call this function with numControls set to zero.
    return gate_type_available[gate.gateTypeId()][numControls];
  }


  long CircuitContext::gateCost(int gateTypeId, int numControls) const
  {
    // For use when we want to avoid creating a gate just to check availability.
    return gate_cost[gateTypeId].cost(numControls);
  }


  long CircuitContext::gateCost(const Gate& gate, int numControls) const
  {
    return gate_cost[gate.gateTypeId()].cost(numControls);
  }


  int CircuitContext::numQbits() const
  {
    return num_qbits;
  }


  void CircuitContext::setQbitInputOptions(int qbit, QbitOptions options)
  {
    qbit_input_options[qbit] = options;
  }


  QbitOptions CircuitContext::qbitInputOptions(int qbit) const
  {
    return qbit_input_options[qbit];
  }


  double CircuitContext::meanRandomLength() const
  {
    return mean_random_length;
  }


  double CircuitContext::angleTolerance() const
  {
    return angle_tolerance;
  }


  unique_ptr<Gate> CircuitContext::randomGate() const
  {
    // Creates a random Gate, permitted by the CircuitContext, i.e. of a type in the set selected by the user.

    // Get a gate of random type.
    auto gateType = randInt(0, num_permitted_gate_types);
    auto gate = create_selected_gate(gateType);  // We will create the required GateSimulators in 'random()'.

    // Make it a random gate of that type.
    gate->random();

    // Don't use std::move(gate) - compiler should use 'return value optimization' or 'copy elision'.
    return gate;
  }


  void CircuitContext::output(ostream& out) const
  {
    out << "Number of qbits = " << numQbits() << endl;
    out << "Mean length of a random solution = " << meanRandomLength() << endl;
    out << "Number of available gate types = " << num_permitted_gate_types << endl;
    out << "Gate types available:" << endl;
    for (auto i = 0; i < num_permitted_gate_types; ++i)
    {
      // Get the gate type id and name.
      auto gate = gate_creator[i](*this);
      int id = gate->gateTypeId();
      out << "  " << gate->name();

      // Output which types of controls are available.
      bool allMany = true;
      bool anyMany = false;
      for (int numControls = 2; numControls < numQbits(); ++numControls)
      {
        if (gate_type_available[id][numControls])
        {
          anyMany = true;
        }
        else
        {
          allMany = false;
        }
      }
      if (!anyMany)
      {
        if (gate_type_available[id][0])
        {
          if (gate_type_available[id][1])
          {
            out << " with zero or one control.";
          }
          else
          {
            out << " with no controls.";
          }
        }
        else
        {
          if (gate_type_available[id][1])
          {
            out << " with one control.";
          }
          else
          {
            out << ", but there is no permitted option for controls!";
          }
        }
      }
      else if (allMany)
      {
        if (gate_type_available[id][0])
        {
          if (gate_type_available[id][1])
          {
            out << " with any number of controls.";
          }
          else
          {
            out << " with no controls, or at least two!";
          }
        }
        else
        {
          if (gate_type_available[id][1])
          {
            out << " with at least one control.";
          }
          else
          {
            out << " with at least two controls!";
          }
        }
      }
      else
      {
        out << " with number of controls selected from the set \{ ";
        bool first = true;
        for (int numControls = 0; numControls < numQbits(); ++numControls)
        {
          if (gate_type_available[id][numControls])
          {
            if (first)
            {
              out << numControls;
              first = false;
            }
            else
            {
              out << ", " << numControls;
            }
          }
        }
      }

      out << " Cost = " << gate_cost[id] << "." << endl;
    }
  }


  void CircuitContext::make_available_with_many_controls(int gateTypeId)
  {
    for (int numControls = 2; numControls < num_qbits; ++numControls)
    {
      gate_type_available[gateTypeId][numControls] = true;
    }
  }


  unique_ptr<Gate> CircuitContext::create_selected_gate(int gateTypeNum) const
  {
    return gate_creator[gateTypeNum](*this);
  }


  ostream& operator<<(ostream& out, const CircuitContext& circuitContext)
  {
    circuitContext.output(out);
    return out;
  }

  //--------------------------------------------------------------------------------------------------------------------

  void Insertion::handle_oracle()
  {
    // We restrict Insertions to contain at most one Oracle and this Oracle must be at the end of the gate sequence. The
    // reason for doing this is that, while we can create a big matrix that applies a sequence of non-Oracle gates to a
    // State in one go, the Oracle complicates matters. Having an Oracle in the middle of the sequence requires us to
    // create two big matrices. More than one Oracle makes things even worse. Since an Insertion is designed to be a
    // short sequence of gates, this restriction is not overly restrictive.

    // If the last gate is an Oracle note that the Insertion ends with an Oracle.
    if (sequence_.back()->gateTypeId() == Oracle::gate_type_id)
    {
      includes_oracle = true;
    }

    // Check that the remaining gate sequence, excluding the last gate, contains no Oracle. (Defensive.)
    for (auto i = sequence_.cbegin(); i < sequence_.cend() - 1; ++i)
    {
      const auto& gate = *i;
      assert(gate->gateTypeId() != Oracle::gate_type_id);
    }
  }


  Insertion::Insertion(const Problem& problem, GateSequence&& sequence) :
  problem_(problem),
  sequence_(std::move(sequence))
  {
    // Creates an Insertion from an externally created GateSequence. Used for random Insertions.
    // (Should we simplify Insertions?)

    handle_oracle();
  }


  Insertion::Insertion(const Problem& problem, GateSequence::const_iterator begin, GateSequence::const_iterator end) :
  problem_(problem)
  {
    // Creates an Insertion by copying a sequence from a Solution.
    // (Should we simplify Insertions?)
    for (auto i = begin; i != end; ++i)
    {
      sequence_.push_back((*i)->clone());
    }

    handle_oracle();
  }


  void Insertion::createGateSimulators()
  {
    for (auto& gate : sequence_)
    {
      gate->createSimulators();
    }
  }


  void Insertion::simulate_to_oracle(State& state)
  {
    // Simply apply each gate in turn to the provided state.
    for (auto i = sequence_.cbegin(); i < sequence_.cend() - 1; ++i)  // All but the last gate.
    {
      const auto& gate = *i;
      gate->applyTo(state);
    }
    if (!includes_oracle)
    {
      // Apply the last gate too, if it is not an Oracle.
      sequence_.back()->applyTo(state);
    }
  }


  void Insertion::prepare()
  {
    // Pre-calculate the matrix representing the gate sequence.

    // Allocate memory for the matrix
    auto numQbits = problem_.numQbits();
    vector<vector<std::pair<int, cmplx> > > matrix(dim(numQbits));
    for (auto& row : matrix)
    {
      row.reserve(dim(numQbits));
    }

    // Simulate the gates on each basis state and insert the results into the matrix
    for (auto col = 0; col < dim(numQbits); ++col)
    {
      State state{numQbits, col};
      simulate_to_oracle(state);
      for (auto row = 0; row < dim(numQbits); ++row)
      {
        if (abs(state[row]) > constants::tolerance)
        {
          matrix[row].push_back({col, state[row]});
        }
      }
    }

    vector<int> sys(numQbits);  // (This vector will always contain the same elements.)
    std::iota(sys.begin(), sys.end(), 0);
    vector<int> controls;  // Empty
    simulator_index_manager = make_unique<SimulatorIndexManager>(sys, controls);
    simulator_ = make_unique<SparseGateSimulator>(std::move(matrix));
  }


  int Insertion::length() const
  {
    return static_cast<int>(sequence_.size());
  }


  const Gate& Insertion::gate(int pos) const
  {
    return *sequence_[pos];
  }

  
  void Insertion::applyTo(State& state, std::shared_ptr<int> markedState) const
  {
    state.transform(*simulator_, *simulator_index_manager);
    if (includes_oracle)
    {
      Oracle oracle(problem_.circuitContext());
      oracle.linkToMarkedState(markedState);
      oracle.applyTo(state);
    }
  }


  void Insertion::output(ostream& out) const
  {
    out << "Gate sequence:" << endl;
    for (const auto& gate : sequence_)
    {
      out << *gate;
    }
    out << endl;
    out << "Simulator index manager:" << endl;
    out << *simulator_index_manager << endl;
    out << "Simulator:" << endl;
    out << *simulator_ << endl;
  }


  ostream& operator<<(ostream& out, const Insertion& insertion)
  {
    insertion.output(out);
    return out;
  }

  //--------------------------------------------------------------------------------------------------------------------

  Problem::Problem(const vector<GateCreator>& gateCreator, const vector<PermittedControls>& permittedControls,
                   const vector<GateCost>& gateCosts, double gateCostRange, int numQbits, int maxLength,
                   double meanRandomLength, double badErrorCutoff, double gateOptProb, int numInsertions,
                   int minInsertionLength, int maxInsertionLength, double angleTolerance) :
  circuit_context(gateCreator, permittedControls, gateCosts, numQbits, meanRandomLength, angleTolerance),
  gate_cost_range(gateCostRange),
  max_length(maxLength),
  bad_error_cutoff(badErrorCutoff),
  gate_opt_prob(gateOptProb),
  num_insertions(numInsertions),
  min_insertion_length(minInsertionLength),
  max_insertion_length(maxInsertionLength)
  {
  }


  const CircuitContext& Problem::circuitContext() const
  {
    return circuit_context;
  }


  int Problem::numQbits() const
  {
    return circuit_context.numQbits();
  }


  void Problem::transferQbitInputOptionsToContext()
  {
    for (int qbit = 0; qbit < numQbits(); ++qbit)
    {
      circuit_context.setQbitInputOptions(qbit, qbitInputOptions(qbit));
    }
  }


  int Problem::maxLength() const
  {
    return max_length;
  }


  double Problem::badErrorCutoff() const
  {
    return bad_error_cutoff;
  }


  double Problem::gateOptProb() const
  {
    return gate_opt_prob;
  }


  int Problem::numInsertions() const
  {
    return num_insertions;
  }


  const Insertion& Problem::insertion(int i) const
  {
    return insertions_[i];
  }


  void Problem::prepareForStep()
  {
    insertions_.clear();
//    insertions_.resize(0);   // Function clear() in Visual C++ deallocates memory, or at least used to.

    for (int i = 0; i < num_insertions; ++i)
    {
      int length = randInt(min_insertion_length, max_insertion_length + 1);
      int numNonOracle = 0;
      bool oracleIncluded = false;
      unique_ptr<Gate> oracle;
      GateSequence sequence;
      while (numNonOracle + oracleIncluded < length)
      {
        auto gate = circuitContext().randomGate();
        if (gate->gateTypeId() == Oracle::gate_type_id)
        {
          oracle = std::move(gate);
          oracleIncluded = true;
        }
        else
        {
          sequence.push_back(std::move(gate));
          ++numNonOracle;
        }
      }
      if (oracleIncluded)
      {
        sequence.push_back(std::move(oracle));
      }

      // Create the insertion from the gate sequence, calculate the associated big matrix and add to the insertions
      // list.
      Insertion insertion(*this, std::move(sequence));
      insertion.createGateSimulators(); // Ensure that all SingleTargetGates have their simulators created and attached.
      insertion.prepare();
      insertions_.push_back(std::move(insertion));
    }
  }


  void Problem::output(ostream& out) const
  {
    out << circuit_context;
    out << "Gate cost range = " << gate_cost_range << endl;
  }


  ostream& operator<<(ostream& out, const Problem& problem)
  {
    problem.output(out);
    return out;
  }

  //--------------------------------------------------------------------------------------------------------------------

  Circuit::Circuit(const CircuitContext& context) :
  context_(&context),
  marked_state(nullptr)
  {
  }


  Circuit::Circuit(const CircuitContext& context, GateSequence&& gates) :
  context_(&context),
  gates_(std::move(gates)),
  marked_state(nullptr)
  {
    link_gates_to_marked_state();  // Ensures that the marked_state in each gate is also reset to nullptr. (Defensive.)
  }


  Circuit::Circuit(const Circuit& rhs) :
  context_(rhs.context_),
  marked_state(nullptr)  // DO NOT copy marked_state. Circuits require their own marked state, shared with the...
  {                      // ...Solution and contained Gates.
    // Perform a deep copy of the gates_ vector. (We can't have the new Circuit having pointers into the old one!)

    gates_.reserve(rhs.length());
    for (const auto& gate : rhs.gates_)
    {
      gates_.push_back(gate->clone());
    }

    link_gates_to_marked_state();  // Ensures that the marked_state in each gate is also reset to nullptr. (Defensive.)
  }


  Circuit::Circuit(Circuit&& rhs) :
  context_(std::exchange(rhs.context_, nullptr)),  // Don't have to exchange with nullptr - the (default) destructor...
  gates_(std::move(rhs.gates_)),                   // ...won't delete it anyway.
  marked_state(rhs.marked_state)
  {
  }


  Circuit& Circuit::operator=(const Circuit& rhs)
  {
    // Perform a deep copy of the 'gates_' vector. (We can't have two Circuits having pointers to the same gates!)

    // Check for self assignment.
    if (&rhs != this)
    {
      // Context.
      context_ = rhs.context_;
      marked_state = nullptr;  // DO NOT copy marked_state. Each circuit requires their own.

      // Deep copy of gates.
      gates_.clear();
      gates_.reserve(rhs.length());
      for (const auto& gate : rhs.gates_)  //  STYLE: Do we need const here?
      {
        gates_.push_back(gate->clone());
      }
      link_gates_to_marked_state();  // Ensures that the marked_state in each gate is reset to nullptr. (Defensive.)
    }

    return *this;
  }


  Circuit& Circuit::operator=(Circuit&& rhs)
  {
    // Check for self assignment
    if (&rhs != this)
    {
      context_ = std::exchange(rhs.context_, nullptr);  // Don't have to exchange with nullptr - the (default)...
      marked_state = rhs.marked_state;                  // ...destructor won't delete it anyway.
      gates_ = std::move(rhs.gates_);
    }

    return *this;
  }


  void Circuit::createGateSimulators()
  {
    for (auto& gate : gates_)
    {
      gate->createSimulators();
    }
  }


  bool Circuit::operator==(const Circuit& rhs) const
  {
    if (length() != rhs.length())
    {
      return false;
    }
    for (auto i = 0; i < length(); ++i)
    {
      if (*gates_[i] != *rhs.gates_[i])
      {
        return false;
      }
    }
    return true;
  }


  bool Circuit::operator!=(const Circuit& rhs) const
  {
    return !operator==(rhs);
  }


  bool Circuit::sortBefore(const Circuit& rhs) const
  {
    // To enable sorting, efficient counting of duplicates, etc. Must have a == b, sortBefore(a, b) or sortBefore(b, a).
    if (length() != rhs.length())
    {
      return length() < rhs.length();
    }
    for (auto i = 0; i < length(); ++i)
    {
      if (*gates_[i] != *rhs.gates_[i])
      {
        return gates_[i]->sortBefore(*rhs.gates_[i]);
      }
    }
    return false;
  }


  void Circuit::linkToMarkedState(std::shared_ptr<const int> markedState)
  {
    // Point the circuit to where the marked state is stored. Also ensure that all the gates point to it too.
    marked_state = markedState;
    link_gates_to_marked_state();
  }


  shared_ptr<const int> Circuit::markedState() const
  {
    return marked_state;
  }


  int Circuit::length() const
  {
    return static_cast<int>(gates_.size());
  }


  int Circuit::numMutatable() const
  {
    return static_cast<int>(std::count_if(gates_.begin(), gates_.end(),
                                          [](const auto& gate){return gate->mutatable();}));
  }


  const Gate& Circuit::gate(int pos) const
  {
    return *gates_[pos];
  }


  void Circuit::random()
  {
    // Clear any existing gates.
    gates_.clear();

    // Determine how long the solution should be.
    auto length = randGeometric(context_->meanRandomLength());
    gates_.reserve(length);

    // Create a vector of random gates.
    for (auto i = 0; i < length; ++i)
    {
      // Get a random gate of random type and add it to the circuit.
      gates_.push_back(random_gate());  // random_gate() automatically links the gate to the marked state if necessary.
    }
  }


  State Circuit::simulate(const State& startState) const
  {
    // Simply apply each gate in turn to the provided startState, returning the result.
    State state(startState);
    for (auto& gate : gates_)
    {
      gate->applyTo(state);
    }

    return state;
  }


  vector<State> Circuit::simulateAndLog(const State& startState) const
  {
    // Like simulate(), but returns intermediate states too. These are used to make later calculations more efficient.
    vector<State> visitedStates;
    visitedStates.reserve(length() + 1);

    State state(startState);
    visitedStates.push_back(state);
    for (auto& gate : gates_)
    {
      gate->applyTo(state);
      visitedStates.push_back(state);  // (Note: Result of applyTo is not merely copied into its parameter here.)
    }

    return visitedStates;
  }


  vector<State> Circuit::reverseSimulateAndLog(const State& targetState) const
  {
    // Simulates the application of the inverse gates, in reverse order, applied to the target state. Used to speed up
    // subsequent calculations.
    //
    // Given target state t, start state s and gates ABCDEF, which we simulate, we may find that we subsequently wish to
    // simulate ABCGEF. We can do this by taking the stored result for CBAs, multiplying it by G to get GCBAs and then
    // take the dot product with the stored result for E^F^t. Hence, by doing both a forward and reverse simulation of
    // the circuit and storing the results, we can significantly reduce the effort required to simulate minor variations
    // of the circuit.
    //
    // (Unlike in algorithm 1, here the vector index indicates the position in the circuit, counting from the left. This
    // comes from the addition of the call to reverse().)
    vector<State> visitedStates;
    visitedStates.reserve(length() + 1);

    State state(targetState);
    visitedStates.push_back(state);
    for (auto pos = static_cast<int>(length()) - 1; pos >= 0; --pos)
    {
      gates_[pos]->applyInvTo(state);
      visitedStates.push_back(state);
    }

    reverse(begin(visitedStates), end(visitedStates));

    return visitedStates;
  }


  void Circuit::makeCanonicalForm(int startPos)
  {
    // To avoid dealing with very many versions of the same circuit that merely have gates reordered, we convert to a
    // canonical form. Gates are reordered, as allowed given the rules encoded in swaps(), so that gates that are
    // 'sortedBefore' other gates come earlier where possible. Implemented as a simple modification to insertion sort.
    //
    // Parameter startPos allows us to save time - it indicates that gates before startPos are already in the correct
    // order, i.e. that insertion sort up to this position will simply leave the gates in place.
    //
    // (Is this even in this version of the algorithm, given that we do not use a solution cache? Moving the SwapGates
    // to the end and sorting them makes sense, since this allows simplification of the SwapGates. Using a canonical
    // form in the Store makes sense too, to prevent it filling with myriad copies of essentially the same circuit.
    // However, for the running of the algorithm, just soring out the SwapGate should be sufficient.)

    // Create a vector of Gates, of the same length as the circuit, to act as storage for replacement gates, i.e. what
    // gates will become after some other gate is moved over them. (Note that, at present, gates only change in this way
    // if the moving gate is a SwapGate, and that a SwapGate will only (at present) be moved over other SwapGates.
    // However, we may, in future, consider HX -> ZH to be a valid 'swap'.)
    GateSequence changedGates(length());

    // Consider each element in turn.
    for (int movingPos = startPos; movingPos < length(); ++movingPos)
    {
      auto movingGate = std::move(gates_[movingPos]);  // We std::move it back into the circuit's gate sequence later.
      auto bestPos = movingPos;             // Current best position to insert the gate and...
      auto bestGate = movingGate->clone();  // ...what the inserted gate will be.

      // Consider swapping with previous gates until this is no longer permitted.
      auto candidatePos = movingPos - 1;
      while (candidatePos >= 0)
      {
        if (gates_[candidatePos]->canSimplySwap(*movingGate))
        {
          changedGates[candidatePos] = {};  // Indicates that the gate does not change.
        }
        else if (gates_[candidatePos]->canSwapSwap(*movingGate))
        {
          changedGates[candidatePos] = gates_[candidatePos]->rightSwapMoverChange(*movingGate);
          auto changedMover = gates_[candidatePos]->leftSwapMoverChange(*movingGate);
          if (changedMover)
          {
            movingGate = std::move(changedMover);
          }
        }
        else if (gates_[candidatePos]->canHSwap(*movingGate))
        {
          changedGates[candidatePos] = gates_[candidatePos]->rightHMoverChange(*movingGate);
          auto changedMover = gates_[candidatePos]->leftHMoverChange(*movingGate);
          if (changedMover)
          {
            movingGate = std::move(changedMover);
          }
        }
        else if (gates_[candidatePos]->canRSwap(*movingGate))
        {
          changedGates[candidatePos] = gates_[candidatePos]->rightRMoverChange(*movingGate);
          auto changedMover = gates_[candidatePos]->leftRMoverChange(*movingGate);
          if (changedMover)
          {
            movingGate = std::move(changedMover);
          }
        }
       else
        {
          break;  // Cannot swap.
        }

        // We want to move the gate to the earliest position that currently contains a gate that is 'greater' than it,
        // under the constraint that it must be possible to swap the gate to that position.
        if (movingGate->sortBefore(*gates_[candidatePos]))
        {
          // Moving the gate to here results in a new circuit that 'sorts before' the current circuit (and previous
          // 'new' circuits.)
          bestPos = candidatePos;
          bestGate = movingGate->clone();
        }

        --candidatePos;
      }

      // Insert the gate into its new position, shifting the jumped (and possibly modified) gates up.
      int newMovingPos = movingPos;  // Needed if changes are made to jumped gates. Remember that 1 will be added...
      for (auto changePos = movingPos; changePos > bestPos; --changePos)  // ...to movingPos at loop end.
      {
        if (changedGates[changePos - 1])
        {
          gates_[changePos] = std::move(changedGates[changePos - 1]);

          // We may need to go back over bits that were previously sorted if any of these gates is changed.
          // (If we were to make Hadamard::gate_type_id large (i.e. one less than the maximum (which is taken by swap),
          // then the following lines would be unnecessary - Hadamards would only jump back over Swaps, which would be
          // unchanged. (Check))
          if (changePos - 1 == bestPos)
          {
            // The changed gate, being just after the insertion position, is correctly placed. The subsequent gate might
            // now wish to move back.
            newMovingPos = changePos;
          }
          else
          {
            // The changed gate might now wish to move back.
            newMovingPos = changePos - 1;
          }
        }
        else
        {
          gates_[changePos] = std::move(gates_[changePos - 1]);
        }
      }
      gates_[bestPos] = std::move(bestGate);
      movingPos = newMovingPos;
    }
  }


  void Circuit::makeUncanonical()
  {
    // While converting circuits to a canonical form has advantages for simplification and the logistics of cacheing, it
    // may result in a lack of diversity in the population for the genetic operators. Scrambling the solutions, using
    // legal swaps, should increase the number of possible circuits that could result from application of these
    // operators.

    vector<char> moved(length(), false);  // Array of bools indicating which gates have been moved.
                                          // (DO NOT USE vector<bool> here!!)
    // Consider each gate in turn
    int movingPos = 0;
    while (movingPos < length())
    {
      // Get range of permissible new places for the gate.
      int leftmost = leftmost_shift(movingPos, true);
      int rightmost = rightmost_shift(movingPos, true);

      // Select at random from this range
      int newPos = randInt(leftmost, rightmost + 1);

      // Perform the move, also rearranging the array indicating which positions have been considered.
      if (newPos > movingPos)
      {
        swap_gate_right(movingPos, newPos, true);
        std::rotate(moved.begin() + movingPos, moved.begin() + movingPos + 1, moved.begin() + newPos + 1);
      }
      else if (newPos < movingPos)
      {
        swap_gate_left(movingPos, newPos, true);
        std::rotate(moved.begin() + newPos, moved.begin() + movingPos, moved.begin() + movingPos + 1);
      }
      moved[newPos] = true;

      // Move to the next unmoved element.
      while (moved[movingPos] && movingPos < length())
      {
        ++movingPos;
      }
    }
  }


  void Circuit::simplify()
  {
    // My fifth circuit simplification routine. Perform steepest descent search until no 'simplification' can be found
    // that reduces cost.

    makeCanonicalForm();

    auto bestReplacement = best_replacement();
    while (bestReplacement.valid())
    {
      // Make the simplification, making sure the circuit ends up in canonical form.
      apply_replacement(bestReplacement.leftPos(), bestReplacement.rightPos(), bestReplacement.meetPos());

      // Look for further simplifications
      bestReplacement = best_replacement();
    }
  }


  void Circuit::link_gates_to_marked_state()
  {
    // Make sure each gate knows the identity of its 'parent' circuit. Used when circuits are copied (requiring the
    // copied gates to refer to the new parent) and created.
    // (At present, this is called after circuit construction, assignment and crossover. Furthermore,
    // Gate::linkToMarkedState is also called in Circuit::random_gate(). We could, instead, delay linking, calling only
    // this function before circuit evaluation and return to using CircuitContext::randomGate() directly.)
    for (auto& gate : gates_)
    {
      gate->linkToMarkedState(marked_state);
    }
  }


  std::unique_ptr<Gate> Circuit::random_gate() const
  {
    auto newGate = context_->randomGate();
    newGate->linkToMarkedState(marked_state);
    return newGate;
  }


  int Circuit::leftmost_shift(int pos, bool involveSwapGates) const
  {
    // Determine the leftmost position the gate can be pulled to using legal gate swaps. Parameter 'involveSwapGates'
    // indicates whether SwapSwaps are considered.
    // Currently only used for making a circuit 'uncanonical'. Hence efficiency is not a concern and 'involveSwapGates'
    // is always 'true'.

    auto pulledPos = pos;
    auto movingGate = gates_[pos]->clone();
    while (pulledPos > 0)
    {
      if (gates_[pulledPos - 1]->canSimplySwap(*movingGate))
      {
        // Simple swap.
        --pulledPos;
      }
      else
      {
        // HSwaps and SwapSwaps may change the movingGate.
        unique_ptr<Gate> possibleChange;
        if (gates_[pulledPos - 1]->canHSwap(*movingGate))
        {
          // Can HSwap. Note possible change to the moving gate.
          possibleChange = gates_[pulledPos - 1]->leftHMoverChange(*movingGate);
        }
        else if (gates_[pulledPos - 1]->canRSwap(*movingGate))
        {
          // Can RSwap. Note possible change to the moving gate.
          possibleChange = gates_[pulledPos - 1]->leftRMoverChange(*movingGate);
        }
        else if (involveSwapGates && gates_[pulledPos - 1]->canSwapSwap(*movingGate))
        {
          // Can SwapSwap and we permit such operations. Note possible change to the moving gate.
          possibleChange = gates_[pulledPos - 1]->leftSwapMoverChange(*movingGate);
        }
        else
        {
          // No swap available.
          break;
        }

        // Can HSwap or SwapSwap. Change movingGate if necessary and record gate identity in the movedGate array.
        if (possibleChange)  // If possibleChange is null, gate is unchanged.
        {
          movingGate = std::move(possibleChange);
        }
        --pulledPos;
      }
    }

    return pulledPos;
  }


  int Circuit::leftmost_shift(int pos, vector<const Gate*>& movedGate, vector<const Gate*>& pointers,
                              bool involveSwapGates) const
  {
    // Determine the leftmost position the gate can be pulled to using legal gate swaps. Also note how the moving gate
    // changes by recording the form of the gate in the movedGate array. Parameter 'involveSwapGates' indicates whether
    // SwapSwaps are considered.
    // Only called, at present, when simplifying a circuit that is in canonical form. Hence 'involveSwapGates' is always
    // false. Efficiency may be a concern when the circuits being considered are parameterless. We therefore keep gate
    // creation to a minimum. These efficiency worries also explain the (risky) use of bare pointers, rather than
    // shared_ptrs, which were discovered to be too slow. The 'pointers' parameter gets a copy of all the pointers that
    // need to be deleted by the caller.

    auto pulledPos = pos;
    const Gate* movingGate = gates_[pos].get();
    movedGate[pos] = movingGate;
    while (pulledPos > 0)
    {
      if (gates_[pulledPos - 1]->canSimplySwap(*movingGate))
      {
        // Simple swap.
        --pulledPos;
        movedGate[pulledPos] = movingGate;
      }
      else
      {
        // HSwaps and SwapSwaps may change the movingGate.
        unique_ptr<Gate> possibleChange;
        if (gates_[pulledPos - 1]->canHSwap(*movingGate))
        {
          // Can HSwap. Note possible change to the moving gate.
          possibleChange = gates_[pulledPos - 1]->leftHMoverChange(*movingGate);
        }
        else if (gates_[pulledPos - 1]->canRSwap(*movingGate))
        {
          // Can RSwap. Note possible change to the moving gate.
          possibleChange = gates_[pulledPos - 1]->leftRMoverChange(*movingGate);
        }
        else if (involveSwapGates && gates_[pulledPos - 1]->canSwapSwap(*movingGate))
        {
          // Can SwapSwap and we permit such operations. Note possible change to the moving gate.
          possibleChange = gates_[pulledPos - 1]->leftSwapMoverChange(*movingGate);
        }
        else
        {
          // No swap available.
          break;
        }

        // Can HSwap or SwapSwap. Change movingGate if necessary and record gate identity in the movedGate array.
        if (possibleChange)  // If possibleChange is null, gate is unchanged.
        {
          movingGate = possibleChange.release();  // Claims ownership. Dangerous - must delete somewhere, so we...
          pointers.push_back(movingGate);         // ...make a note of the pointer for future deletion.
        }
        --pulledPos;
        movedGate[pulledPos] = movingGate;
      }
    }

    return pulledPos;
  }


  int Circuit::rightmost_shift(int pos, bool involveSwapGates) const
  {
    // Determine the leftmost position the gate can be pulled to using legal gate swaps. Parameter 'involveSwapGates'
    // indicates whether SwapSwaps are considered.
    // Currently only used for making a circuit 'uncanonical'. Hence efficiency is not a concern and 'involveSwapGates'
    // is always 'true'.

    auto pushedPos = pos;
    auto movingGate = gates_[pos]->clone();
    while (pushedPos < length() - 1)
    {
      if (movingGate->canSimplySwap(*gates_[pushedPos + 1]))
      {
        // Simple swap.
        ++pushedPos;
      }
      else
      {
        // HSwaps and SwapSwaps may change the movingGate.
        unique_ptr<Gate> possibleChange;
        if (movingGate->canHSwap(*gates_[pushedPos + 1]))
        {
          // Can HSwap. Note possible change to the moving gate.
          possibleChange = movingGate->rightHMoverChange(*gates_[pushedPos + 1]);
        }
        else if (movingGate->canRSwap(*gates_[pushedPos + 1]))
        {
          // Can RSwap. Note possible change to the moving gate.
          possibleChange = movingGate->rightRMoverChange(*gates_[pushedPos + 1]);
        }
        else if (involveSwapGates && movingGate->canSwapSwap(*gates_[pushedPos + 1]))
        {
          // Can SwapSwap and we permit such operations. Note possible change to the moving gate.
          possibleChange = movingGate->rightSwapMoverChange(*gates_[pushedPos + 1]);
        }
        else
        {
          // No swap available.
          break;
        }

        // Can HSwap or SwapSwap. Change movingGate if necessary and record gate identity in the movedGate array.
        if (possibleChange)  // If possibleChange is null, gate is unchanged.
        {
          movingGate = std::move(possibleChange);
        }
        ++pushedPos;
      }
    }

    return pushedPos;
  }


  int Circuit::rightmost_shift(int pos, vector<const Gate*>& movedGate, vector<const Gate*>& pointers,
                               bool involveSwapGates) const
  {
    // Determine the rightmost position the gate can be pushed to using legal gate swaps. Also note how the moving gate
    // changes by recording the form of the gate in the movedGate array. Parameter 'involveSwapGates' indicates whether
    // SwapSwaps are considered.
    // Only called, at present, when simplifying a circuit that is in canonical form. Hence 'involveSwapGates' is always
    // false. Efficiency may be a concern when the circuits being considered are parameterless. We therefore keep gate
    // creation to a minimum. These efficiency worries also explain the (risky) use of bare pointers, rather than
    // shared_ptrs, which were discovered to be too slow. The 'pointers' parameter gets a copy of all the pointers that
    // need to be deleted by the caller.

    auto pushedPos = pos;
    const Gate* movingGate = gates_[pos].get();
    movedGate[pos] = movingGate;
    while (pushedPos < length() - 1)
    {
      if (movingGate->canSimplySwap(*gates_[pushedPos + 1]))
      {
        // Simple swap.
        ++pushedPos;
        movedGate[pushedPos] = movingGate;
      }
      else
      {
        // HSwaps and SwapSwaps may change the movingGate.
        unique_ptr<Gate> possibleChange;
        if (movingGate->canHSwap(*gates_[pushedPos + 1]))
        {
          // Can HSwap. Note possible change to the moving gate.
          possibleChange = movingGate->rightHMoverChange(*gates_[pushedPos + 1]);
        }
        else if (movingGate->canRSwap(*gates_[pushedPos + 1]))
        {
          // Can RSwap. Note possible change to the moving gate.
          possibleChange = movingGate->rightRMoverChange(*gates_[pushedPos + 1]);
        }
        else if (involveSwapGates && movingGate->canSwapSwap(*gates_[pushedPos + 1]))
        {
          // Can SwapSwap and we permit such operations. Note possible change to the moving gate.
          possibleChange = movingGate->rightSwapMoverChange(*gates_[pushedPos + 1]);
        }
        else
        {
          // No swap available.
          break;
        }

        // Can HSwap or SwapSwap. Change movingGate if necessary and record gate identity in the movedGate array.
        if (possibleChange)  // If possibleChange is null, gate is unchanged.
        {
          movingGate = possibleChange.release();  // Claims ownership. Dangerous - must delete somewhere, so we...
          pointers.push_back(movingGate);         // ...make a note of the pointer for future deletion
        }
        ++pushedPos;
        movedGate[pushedPos] = movingGate;
      }
    }

    return pushedPos;
  }


  Replacement Circuit::best_replacement() const
  {
    // Find the best replacement resulting from considering simplifications of pairs of gates, provided it improves cost
    // or keeps cost fixed but increases the number of useful degrees of freedom.
    // (We assume that the circuit is in canonical form. In particular, this means that all SwapGates should be at the
    // end (with the possible exception of Oracles).

    // First make a note of how far each gate may be shifted to the left and right and what the gate would look like in
    // each attainable position.
    vector<int> leftmost(length());
    vector<int> rightmost(length());
    vector<vector<const Gate*> > movedGate(length(), vector<const Gate*>(length(), nullptr));
    vector<const Gate*> pointers;
    pointers.reserve(2 * length());
    for (int pos = 0; pos < length(); ++pos)
    {
      leftmost[pos] = leftmost_shift(pos, movedGate[pos] , pointers, false);
      rightmost[pos] = rightmost_shift(pos, movedGate[pos] , pointers, false);
    }

    // Search through all possible replacements.
    long bestImprovement{0};
    Replacement bestReplacement;

    // Start by checking for any gates that can be moved to the start of the circuit and then discarded as redundant.
     for (int rightPos = 0; rightPos < length(); ++rightPos)
    {
      if (leftmost[rightPos] == 0)
      {
        // Gate can be shifted to the beginning.
        auto& newRightGate = movedGate[rightPos][0];
        if (newRightGate->cancelsAtStart())
        {
          auto improvement = newRightGate->cost();
          if (improvement > bestImprovement)
          {
            bestImprovement = improvement;
            bestReplacement = {-1, rightPos, -1};
          }
        }
      }
    }

    // Then look through all possible simplifications of gate pairs
    for (int leftPos = 0; leftPos < length() - 1; ++leftPos)
    {
      for (int rightPos = leftPos + 1; rightPos < length(); ++rightPos)
      {
        if (rightmost[leftPos] >= leftmost[rightPos] - 1)
        {
          // It is possible to move the two gates, using legal gate swaps, so that they meet, giving us the chance to
          // simplify.
          int meetPos = std::min(rightmost[leftPos], rightPos - 1);  // Choosing the rightmost position reduces work...
          auto& newLeftGate = movedGate[leftPos][meetPos];           // ...makeCanonicalForm().
          auto& newRightGate = movedGate[rightPos][meetPos + 1];
          auto improvement = newLeftGate->canSimplify(*newRightGate);
          if (improvement > bestImprovement)
          {
            bestImprovement = improvement;
            bestReplacement = {leftPos, rightPos, meetPos};
          }
        }
      }
    }
    deletePointers(pointers);
    return bestReplacement;
  }


  void Circuit::swap_gate_right(int moverPos, int finalPos, bool involveSwapGates)
  {
    // Perform gate swaps such that the gate at 'moverPos' moves to the right until it reaches 'finalPos'. In the
    // process we adjust the moving gate and any gates jumped, if required. Parameter 'involveSwapGates' indicates
    // whether SwapSwaps are considered. (When simplifying a circuit that is in canonical form, they need not be
    // considered since all the SwapGates are out of the way. When 'un-canonicalizing' a circuit, SwapGates may be in
    // the way.)
    // This function assumes that it has already been determined that it is possible to move the gate into the desired
    // position. If not, it will throw an exception.

    assert(finalPos >= moverPos);
    auto mover = std::move(gates_[moverPos]);
    for (int pos = moverPos; pos < finalPos; ++pos)
    {
      if (mover->canSimplySwap(*gates_[pos + 1]))
      {
        // Just a simple swap
        gates_[pos] = std::move(gates_[pos + 1]);
      }
      else
      {
        // HSwaps and SwapSwaps may change the identify of one of the gates. (This code can handle both gates...
        // ...changing, if necessary.)
        unique_ptr<Gate> newMover;
        unique_ptr<Gate> newJumpedGate;
        if (mover->canHSwap(*gates_[pos + 1]))
        {
          // We can HSwap. Make a note of gate changes.
          newMover = mover->rightHMoverChange(*gates_[pos + 1]);
          newJumpedGate = mover->leftHMoverChange(*gates_[pos + 1]);
        }
        else if (mover->canRSwap(*gates_[pos + 1]))
        {
          // We can RSwap. Make a note of gate changes.
          newMover = mover->rightRMoverChange(*gates_[pos + 1]);
          newJumpedGate = mover->leftRMoverChange(*gates_[pos + 1]);
        }
        else if (involveSwapGates && mover->canSwapSwap(*gates_[pos + 1]))
        {
          // We can SwapSwap, and permit such operations. Make a note of gate changes.
          newMover = mover->rightSwapMoverChange(*gates_[pos + 1]);
          newJumpedGate = mover->leftSwapMoverChange(*gates_[pos + 1]);
        }
        else
        {
          std::cerr << "Moving gate = " << *mover << std::endl;
          std::cerr << "Gate being jumped = " << *gates_[pos + 1] << std::endl;
          std::cerr << "Circuit:" << std::endl << *this;
          throw std::logic_error("Impossible to move gate to the desired position using the permitted swap operations");
        }

        // Perform the swap.
        if (newMover)
        {
          mover = std::move(newMover);
        }
        if (newJumpedGate)
        {
          gates_[pos] = std::move(newJumpedGate);
        }
        else
        {
          gates_[pos] = std::move(gates_[pos + 1]);
        }
      }
    }

    // Insert the moving gate into its final position.
    gates_[finalPos] = std::move(mover);
  }


  void Circuit::swap_gate_left(int moverPos, int finalPos, bool involveSwapGates)
  {
    // Perform gate swaps such that the gate at 'moverPos' moves to the left until it reaches 'finalPos'. In the process
    // we adjust the moving gate and any gates jumped, if required. Parameter 'involveSwapGates' indicates whether
    // SwapSwaps are considered. (When simplifying a circuit that is in canonical form, they need not be considered
    // since all the SwapGates are out of the way. When 'un-canonicalizing' a circuit, SwapGates may be in the way.)
    // This function assumes that it has already been determined that it is possible to move the gate into the desired
    // position. If not, it will throw an exception.

    assert(finalPos <= moverPos);
    auto mover = std::move(gates_[moverPos]);
    for (int pos = moverPos; pos > finalPos; --pos)
    {
      if (gates_[pos - 1]->canSimplySwap(*mover))
      {
        // Just a simple swap
        gates_[pos] = std::move(gates_[pos - 1]);
      }
      else
      {
        // HSwaps and SwapSwaps may change the identity of one of the gates. (This code can handle both gates...
        // ...changing, if necessary.)
        unique_ptr<Gate> newMover;
        unique_ptr<Gate> newJumpedGate;
        if (gates_[pos - 1]->canHSwap(*mover))
        {
          // We can HSwap. Make a note of gate changes.
          newMover = gates_[pos - 1]->leftHMoverChange(*mover);
          newJumpedGate = gates_[pos - 1]->rightHMoverChange(*mover);
        }
        else if (gates_[pos - 1]->canRSwap(*mover))
        {
          // We can RSwap. Make a note of gate changes.
          newMover = gates_[pos - 1]->leftRMoverChange(*mover);
          newJumpedGate = gates_[pos - 1]->rightRMoverChange(*mover);
        }
        else if (involveSwapGates && gates_[pos - 1]->canSwapSwap(*mover))
        {
          // We can SwapSwap, and permit such operations. Make a note of the gate changes.
          newMover = gates_[pos - 1]->leftSwapMoverChange(*mover);
          newJumpedGate = gates_[pos - 1]->rightSwapMoverChange(*mover);
        }
        else
        {
          throw std::logic_error("Impossible to move gate to the desired position using the permitted swap operations");
        }

        // Perform the swap.
        if (newMover)
        {
          mover = std::move(newMover);
        }
        if (newJumpedGate)
        {
          gates_[pos] = std::move(newJumpedGate);
        }
        else
        {
          gates_[pos] = std::move(gates_[pos - 1]);
        }
      }
    }

    // Insert the moving gate into its final position.
    gates_[finalPos] = std::move(mover);
  }


  void Circuit::apply_replacement(int leftPos, int rightPos, int meetPos)
  {
    // The gates at leftPos and rightPos can be moved so as to meet at meetPos and a simplication/replacement can be
    // performed. Gates will change in the process of moving if HSwaps are used.
    // If leftPos (and meetPos) are -1, this means that the gate at rightPos is moved to the start and then eliminated
    // instead.

    if (leftPos == -1)
    {
      // Move the right gate to the start, adjusting any gates jumped if an HSwap is used, then remove the gate.
      swap_gate_left(rightPos, 0, false);
      gates_.erase(gates_.begin());
    }
    else
    {
      // Move the left gate to the right, adjusting it and any gates jumped if an HSwap is used.
      swap_gate_right(leftPos, meetPos, false);

      // Move the right gate to the left, adjusting it and any gates jumped if a HSwap is used.
      swap_gate_left(rightPos, meetPos + 1, false);

      // Now get the replacement sequence.
      GateSequence newSequence = gates_[meetPos]->simplification(*gates_[meetPos + 1]);
      auto replacementSize = static_cast<int>(newSequence.size());
      assert(replacementSize < 3);  // Circuit shouldn't grow. This means we can avoid expensive insertions into the...
                                    // circuit's gate sequence.
      // Replace the old gates with the new ones.
      gates_.erase(gates_.begin() + meetPos + replacementSize, gates_.begin() + meetPos + 2);
      std::copy(make_move_iterator(newSequence.begin()), make_move_iterator(newSequence.end()),
                gates_.begin() + meetPos);
    }

    // Put the new circuit into canonical form.
    makeCanonicalForm();
  }


  bool Circuit::gates_available() const
  {
    // A simple check on the availability of all the gates in the circuit.
    for (auto& gate : gates_)
    {
      if (!gate->available())
      {
        cerr << "Error: " << *gate << "is unavailable, yet appears in the circuit after simplification." << endl;
        return false;
      }
    }
    return true;
  }


  void Circuit::output(ostream& out) const
  {
    for (auto& gate : gates_)  // Must have & as 'gate' is now a unique_ptr.
    {
      if (gate)
      {
        out << *gate;
      }
      else
      {
        out << "Null gate!" << endl;
      }
    }
  }


  ostream& operator<<(ostream& out, const Circuit& circuit)
  {
    circuit.output(out);
    return out;
  }

  //--------------------------------------------------------------------------------------------------------------------

  // A base class for all quantum circuit Solution classes. Doesn't need to implement everything required by the genetic
  // algorithms, as classes such as NSGA2Solution will be templated on the derived Solution classes, not this base
  // class. (Indeed, the genetic algorithm code has no need for a Solution hierarchy, being template based.)

  Solution::Solution(const Problem& problem) :
  problem_(problem),
  marked_state(make_shared<int>(-1)),
  error_(problem.numObj() - 1, -1.0)  // Member 'evaluated_' set to false by in-class initializer.
  {
  }


  bool Solution::operator==(const Solution& rhs) const
  {
    // Equality operator just considers the solution itself, i.e. the sequence of gates. This means that we can use this
    // to detect that an unevaluated solution is, in fact, the same as an evaluated one.

    if (!circuit_ || !rhs.circuit_)
    {
      // One of the circuits hasn't even been constructed yet. While it is possible that the two circuits may end up
      // being the same, we will return false. Note that this only occurs, at present, when an overlarge circuit is
      // considered for entry into the store. Such a solution is highly unlikely to get into the store - it is assigned
      // the worst possible error, as well as being huge - and even if it does, it should soon get dominated by another
      // solution. (Unpleasant!)
      return false;
    }
    return *circuit_ == *rhs.circuit_;
  }


  bool Solution::operator!=(const Solution& rhs) const
  {
    return !operator==(rhs);
  }


  bool Solution::sortBefore(const Solution& rhs) const
  {
    // To enable sorting, efficient counting of duplicates, etc. Must have a == b, sortBefore(a, b) or sortBefore(b, a).

    assert(circuit_ && rhs.circuit_);
    return circuit_->sortBefore(*rhs.circuit_);
  }


  void Solution::random()
  {
    // Create a random circuit and mark as unevaluated and put it into canonical form.

    // Set various members to indicate that the solution has been created with no parents, no insertions. Also indicate
    // that it has not been evaluated.
    left_parent = nullptr;
    left_circuit = nullptr;
    right_parent = nullptr;
    right_circuit = nullptr;
    new_gate = nullptr;
    insertion_num = -1;
    forward_intermediate_states.resize(0);
    reverse_intermediate_states.resize(0);
    overlaps_.resize(0);
    symbolic_overlaps.resize(0);
    total_cost = -1;
    evaluated_ = false;
    intermediate_states_recorded = false;

    // Construct the circuit.
    do
    {
      circuit_ = make_shared<Circuit>(problem_.circuitContext());
      circuit_->random();
    }
    while (circuit_->length() > problem_.maxLength());

    // Ensure that the circuit and all the contained gates link to the same marked state.
    circuit_->linkToMarkedState(marked_state);
  }


  void Solution::tweak_gate()
  {
    // Merely cuts out a single, mutatable gate, by defining the cut points on either side of it, then reinserts it at
    // the same point! While this would seem to do nothing, during evaluation the gate is either optimized or mutated,
    // changing the circuit, but not the circuit structure. Use of this function implies that we are not using
    // crossover.
    //
    // (At present, there is a decent chance that we might attempt to optimize a gate that has just been optimized.)

    // Make parents the same, since we are not performing crossover.
    right_parent = left_parent;
    right_circuit = left_circuit;

    // Choose a mutatable gate at random.
    assert(left_circuit->numMutatable() > 0);
    auto r = randInt(0, left_circuit->numMutatable());
    int pos{0};
    while (r >= 0)
    {
      if (left_circuit->gate(pos).mutatable())
      {
        --r;
      }
      ++pos;
    }
    --pos;

    // Set cut points that remove the gate.
    left_cut = pos;
    right_cut = pos + 1;

    // Make the gate to be inserted to be that which has just been removed.
    new_gate = left_circuit->gate(pos).clone();
  }


  void Solution::single_point_crossover()
  {
    // Randomly select the cut points in the parent solutions. A parent can provide all of its gates. However, we insist
    // that each parent must provide at least one gate. Hence the cut point in 'left_parent', which provides the gates
    // to the left of the cut, must take a position from one to length, inclusive. The cut point in 'right_parent',
    // which provides the gates to the right of the cut, must take a position from zero to length - 1 inclusive.
    left_cut = randInt(0, left_circuit->length()) + 1;
    right_cut = randInt(0, right_circuit->length());
  }


  void Solution::cut_circuit()
  {
    // Creates a cut point in the circuit to prepare for the insertion of a Gate or GateSequence at the cut point. This
    // divides the circuit into a left section and a right section, without removing anything. Any subsequent insertion
    // will happen at the cut point. Use of this function implies that we are not performing crossover. (This is like
    // remove_sequence() for an empty sequence.)

    // Make parents the same, since we are not performing crossover.
    right_parent = left_parent;
    right_circuit = left_circuit;

    // Choose a cut point at random. For a circuit of length 7, cut points go from 0 (before the first gate) to 7 (after
    // the last).
    left_cut = right_cut = randInt(0, left_circuit->length() + 1);
  }


  void Solution::remove_gate()
  {
    // Make parents the same, since we are not performing crossover.
    right_parent = left_parent;
    right_circuit = left_circuit;
    int length = left_circuit->length();
    assert(length > 1);

    // Choose a single gate at random.
    left_cut = randInt(0, length);
    right_cut = left_cut + 1;
  }


  void Solution::remove_sequence()
  {
    // Make parents the same, since we are not performing crossover.
    right_parent = left_parent;
    right_circuit = left_circuit;
    int length = left_circuit->length();
    assert(length > 1);

    // Choose two cut points at random, ensuring that they do not both refer to the same point. We also ensure that we
    // do not removed the entire circuit! Given a circuit of length 7, cut points go from 0 (before the first gate) to 7
    // (after the last gate). (Recall that randInt(x, y) gets integer r such that x <= r < y.)
    left_cut = randInt(0, length + 1);
    if (left_cut == 0)
    {
      right_cut = randInt(1, length);
    }
    else if (left_cut == length)
    {
      right_cut = left_cut;
      left_cut = randInt(1, length);
    }
    else
    {
      right_cut = randInt(0, length);
      if (right_cut >= left_cut)
      {
        ++right_cut;
      }
      else
      {
        swap(left_cut, right_cut);
      }
    }
  }


  void Solution::replicate_sequence()
  {
    // Make parents the same, since we are not performing crossover.
    right_parent = left_parent;
    right_circuit = left_circuit;
    int length = left_circuit->length();

    // Choose two cut points at random, ensuring that they do not both refer to the same point. Given a circuit of
    // length 7, cut points go from 0 (before the first gate) to 7 (after the last gate). (Recall that randInt(x, y)
    // gets integer r such that x <= r < y.)
    // Note that, confusingly, left_cut should be to the right of right_cut. Recall that we will take all gates to the
    // left of left_cut and then add all the gates to the right of right_cut. In this case, by ensuring that
    // left_cut > right_cut, we replicate the gates between the cuts.
    left_cut = randInt(0, length + 1);
    right_cut = randInt(0, length);
    if (right_cut >= left_cut)
    {
      ++right_cut;
      swap(left_cut, right_cut);
    }
  }


  void Solution::insert_gate()
  {
    // Inserts a Gate into the crossover/removal/cut point created by the above three functions.

    new_gate = problem_.circuitContext().randomGate();
    new_gate->createSimulators();  // Ensure, if this is a SingleTargetGate, that the gate has its simulators prepared.
    new_gate->linkToMarkedState(marked_state);  // Some gates (Oracles) require access to the marked state.
  }


  void Solution::insert_sequence()
  {
    // Inserts a GateSequence into the crossover/removal/cut point created by functions above.

    insertion_num = randInt(0, problem_.numInsertions());
  }


  bool Solution::mutation_valid(int mutationType) const
  {
    // Used to prevent the application of unsuitable genetic operators, e.g. crossover when one of the parents is empty,
    // or anything involving gate removal when the parent circuit has only one gate.

    // If the problem includes no sequences to insert, we cannot insert sequences. (Basically the user has switched off
    // this option.)
    if (problem_.numInsertions() == 0 && mutationType >= 10)
    {
      return false;
    }

    // If the primary (left) parent has no mutatable gates, we cannot 'tweak' a gate.
    if (left_circuit->numMutatable() == 0 && mutationType == 0)
    {
      return false;
    }

    // If the primary (left) parent is empty, the only permittable mutations are those that insert a gate or sequence.
    if (left_circuit->length() == 0)
    {
      return mutationType == 5 || mutationType == 10;
    }

    // If the donor (right) parent is empty, we do not permit crossover, as this parent cannot contribute anything.
    if (right_circuit->length() == 0)
    {
      if (mutationType == 1 || mutationType == 6 || mutationType == 11)
      {
        return false;
      }
    }

    // If the primary (left) parent has only one gate, we do not permit removal or replacement.
    if (left_circuit->length() == 1)
    {
      if (mutationType == 2 || mutationType == 3 || mutationType == 7 || mutationType == 8 || mutationType == 12 ||
          mutationType == 13)
      {
        return false;
      }
    }
    return true;
  }

  void Solution::mutate()
  {
    int r;
    do
    {
      r = randInt(0, 15);
    }
    while (!mutation_valid(r));

    switch (r)
    {
      case 0:
        // 'Tweak' a mutatable gate. This may be either gate angle optimization or mutation.
        tweak_gate();
        break;

      case 1:
        // Single point crossover. Quick evaluation of the resulting solution should be O(n^2), where n = 2^q and q is
        // the number of qbits.
        single_point_crossover();
        break;

      case 2:
        // Remove a single gate. Quick evaluation should be O(n^2). We have added this option, despite the fact that
        // remove_sequence() will sometimes remove a single gate, since we may wish to encourage small, circuit
        // reducing, genetic operators.
        remove_gate();
        break;

      case 3:
        // Remove gate sequence. Quick evaluation should be O(n^2).
        remove_sequence();
        break;

      case 4:
        // Replicate gate sequence. Quick evaluation should be O(n^2).
        replicate_sequence();
        break;

      case 5:
        // Insert gate. O(n^3) in most cases.
        cut_circuit();
        insert_gate();
        break;

      case 6:
        // Crossover and insert gate. O(n^3) in most cases.
        single_point_crossover();
        insert_gate();
        break;

      case 7:
        // Replace gate. O(n^3) in most cases.
        remove_gate();
        insert_gate();
        break;

      case 8:
        // Replace sequence with gate. O(n^3) in most cases.
        remove_sequence();
        insert_gate();
        break;

      case 9:
        // Replicate gate sequence and insert a gate between the copies. O(n^3) in most cases.
        replicate_sequence();
        insert_gate();
        break;

      case 10:
        // Insert gate sequence. O(n^3). Sequence should be pre-simulated.
        cut_circuit();
        insert_sequence();
        break;

      case 11:
        // Crossover and insert gate sequence at crossover point. O(n^3). Sequence should be pre-simulated.
        single_point_crossover();
        insert_sequence();
        break;

      case 12:
        // Replace gate with sequence. O(n^3). Sequence should be pre-simulated.
        remove_gate();
        insert_sequence();
        break;

      case 13:
        // Replace gate sequence. O(n^3). Sequence should be pre-simulated.
        remove_sequence();
        insert_sequence();
        break;

      case 14:
        // Replicate gate sequence and insert a gate sequence between the copies.
        replicate_sequence();
        insert_sequence();
        break;
    }
  }


  void Solution::construct_circuit()
  {
    if (circuit_)
    {
      return;
    }

    if (!left_parent || !right_parent)
    {
      throw std::logic_error("Attempting to construct child circuit with at least one absent parent.");
    }

    bool canReduce = false;
    GateSequence newGateReduction;
    if (new_gate && new_gate->canReduce())
    {
      canReduce = true;
      newGateReduction = new_gate->reduction();
    }

    // New circuit is obtained by taking the right parent and replacing those gates that come before the cut with those
    // from the left parent that come before the cut.
    auto length = right_circuit->length() - right_cut + left_cut;

    // We may also insert additional gates at the join.
    if (canReduce)
    {
      length += newGateReduction.size();
    }
    else if (new_gate)
    {
      ++length;
    }
    if (insertion_num >= 0)
    {
      length += problem_.insertion(insertion_num).length();
    }

    // Create the gate sequence.
    GateSequence sequence;
    sequence.reserve(length);

    // Left chunk, from left parent.
    for (int i = 0; i < left_cut; ++i)
    {
      sequence.push_back(left_circuit->gate(i).clone());
    }

    // Any insertions at the join.
    if (canReduce)
    {
      sequence.insert(sequence.end(), std::make_move_iterator(newGateReduction.begin()),
                      std::make_move_iterator(newGateReduction.end()));
    }
    else if (new_gate)
    {
      sequence.push_back(std::move(new_gate));  // Note that this means that new_gate is emptied, so we cannot use it...
    }                                           // ...after circuit construction.
    if (insertion_num >= 0)
    {
      const auto& insertion = problem_.insertion(insertion_num);
      for (int i = 0; i < insertion.length(); ++i)
      {
        sequence.push_back(insertion.gate(i).clone());
      }
    }

    // Right chunk, from right parent.
    for (int i = right_cut; i < right_circuit->length(); ++i)
    {
      sequence.push_back(right_circuit->gate(i).clone());
    }

    // Create the Circuit.
    circuit_ = make_shared<Circuit>(problem_.circuitContext(), std::move(sequence));

    // Some unnecessary cleanup. We have both evaluated and constructed the circuit, so we don't need to keep the
    // pointers to the parents.
    left_parent = right_parent = nullptr;
    left_circuit = right_circuit = nullptr;
  }


  void Solution::simulate_and_log(const State& startState)
  {
    forward_intermediate_states.push_back(circuit_->simulateAndLog(startState));
  }


  void Solution::reverse_simulate_and_log(const State& targetState)
  {
    reverse_intermediate_states.push_back(circuit_->reverseSimulateAndLog(targetState));
  }


  void Solution::record_intermediate_states()
  {
    for (auto i = 0; i < num_simulations(); ++i)
    {
      // Prepare the circuit, using the function provided by the derived class, e.g grover::Solution. This might, for
      // example, set the marked state for the Oracle.
      prepare_(i);

      // Simulate in both directions, using the start and target states provided by the derived class.
      simulate_and_log(start_state(i));
      reverse_simulate_and_log(target_state(i));
    }
    intermediate_states_recorded = true;
  }


  void Solution::calculate_overlaps()
  {
    // Used only if intermediate states have been recorded ('full' evaluation) and only on solutions in the initial
    // population.
    assert(intermediate_states_recorded);
    overlaps_.resize(num_simulations());
    for (auto i = 0; i < num_simulations(); ++i)
    {
      // We could choose any point in the circuit to calculate the overlap - for simplicity we choose the start. As a
      // result, we actually calculate the overlap of the start state with the result of simulating backwards from the
      // target state.
      overlaps_[i] = stateOverlap(forward_intermediate_states[i][0], reverse_intermediate_states[i][0]);
    }
  }


  void Solution::quickly_evaluate()
  {
    overlaps_.resize(num_simulations());
    for (auto i = 0; i < num_simulations(); ++i)
    {
      // Prepare the circuit, using the function provided by the derived class, e.g grover::Solution. This might, for
      // example, set the marked state for the Oracle.
      prepare_(i);

      State state = left_parent->forward_intermediate_states[i][left_cut];
      if (new_gate)
      {
        new_gate->linkToMarkedState(marked_state);
        new_gate->applyTo(state);
      }
      if (insertion_num >= 0)
      {
        problem_.insertion(insertion_num).applyTo(state, marked_state);
      }
      overlaps_[i] = stateOverlap(state, right_parent->reverse_intermediate_states[i][right_cut]);
    }

    // Calculate the error based objectives, again using the function provided by the derived class.
    calculate_error();

    evaluated_ = true;
  }


  void Solution::optimize_and_evaluate()
  {
    // Circuit has an insertion of a single gate with a parameter that can be optimized. We perform symbolic simulation,
    // obtaining a simple expression to be minimized or maximized. The optimal parameter value is then found either
    // algebraically (if possible) or numerically (if not).
    symbolic_overlaps.assign(num_simulations(), vector<cmplx>(3, 0.0));
    for (auto i = 0; i < num_simulations(); ++i)
    {
      // Prepare the circuit, using the function provided by the derived class, e.g grover::Solution. This might, for
      // example, set the marked state for the Oracle.
      prepare_(i);

      new_gate->linkToMarkedState(marked_state);
      vector<State> stateParts(3, left_parent->forward_intermediate_states[i][left_cut]);
      for (int j = 0; j < 3; ++j)
      {
        State state = left_parent->forward_intermediate_states[i][left_cut];
        new_gate->applyPartTo(state, j);
        symbolic_overlaps[i][j] = stateOverlap(state, right_parent->reverse_intermediate_states[i][right_cut]);
      }
    }

    // Calculate symbolically a quantity that, when minimized, also minimizes the overall error. (In the case of
    // Grovers, this is just the overall error minus one. For Fourier and Tofolli, it is -(1 - E)^2, where E is the
    // overall error.) The coefficients that we get are stored in 'primary_error_coeffs'.
    calculate_symbolic();

    // Get the newly inserted gate to calculate the best parameter for itself, given the symbolic primary error
    new_gate->optimizeAngle(primary_error_coeffs);

    // Calculate the error based objectives.
    quickly_evaluate();
  }


  void Solution::evaluateObjectives()
  {
    // Function required by the MOMH library. This function performs two roles. (This is not ideal, but fixing this
    // would require adding a member function to Population, which would require some thought about the MOGA framework.)
    // The first is the evaluation of objectives. This is performed as quickly as possible, using intermediate states
    // recorded in the parents, if they exist, to speed up the process. The second is the recording of these
    // intermediate states. To do this, full simulations, both forward and in reverse, are applied to all of the input
    // and target states. The expense of full simulations means that we only record intermediate states for solutions
    // that are good enough to enter the adult population. Moreover, we will evaluate many child solutions for each
    // new adult, taking full advantage of the quick evaluations.
    //
    // A quick evaluation may be followed by the creation of the actual circuit and circuit simplification, but only if
    // the error rate is good enough to make such work worthwhile.

    int l = length();
    if (l > problem_.maxLength())
    {
      bad_ = true;
      assign_worst_cost();
      assign_worst_error();  // Overkill, now that we can simply label solutions as 'bad'?
      evaluated_ = true;  // (Need to ensure that this solution is also marked as unworthy!!!)
      return;
    }

    if (evaluated_)
    {
      if (intermediate_states_recorded)
      {
        // Circuit is already fully evaluated. There is no point in repeating the task!
        return;
      }

      // Solution is a child that has been (quickly) evaluated already and has just been transferred to the adult
      // population. We need to evaluate the intermediate states to allow for the quick evaluation of children.
      if (error_[0] >= problem_.badErrorCutoff())  // (Replace with use of a boolean flag called 'simplified_'.)
      {
        // Solution is of sufficiently low quality that it was not simplified. However, since it has reached the adult
        // population, we simplify now and reevaluate its cost.
        circuit_->simplify();
        circuit_->createGateSimulators();  // Create GateSimulators for any SingleTargetGate that hasn't got one.
        evaluate_cost();
      }

      assert(length() <= problem_.maxLength());
      circuit_->linkToMarkedState(marked_state);  // Ensure that all gates in the circuit link to the marked state.
      record_intermediate_states();
    }
    else if (!left_parent)
    {
      // Solution has no parents, so must be an initial solution with complete (and simplified) circuit. Count up the
      // gate cost, and find the intermediate states obtained in forward and reverse simulations. Calculate overlaps,
      // from which the derived Solution class can produce error based objectives, using calculate_error().
      assert(circuit_);
      circuit_->linkToMarkedState(marked_state);  // Ensure that all gates in the circuit link to the marked state.
      circuit_->simplify();   // May result in new Gates appearing without GateSimulators
      circuit_->createGateSimulators();  // Create GateSimulators for any SingleTargetGate that hasn't got one.
      evaluate_cost();

      record_intermediate_states();
      calculate_overlaps();
      calculate_error();

      evaluated_ = true;
    }
    else
    {
      // Solution is an unevaluated child. Perform a quick evaluation, using the intermediate states stored in the
      // parents. If the solution is not awful, follow up by constructing the circuit and performing circuit
      // simplification.
      // (Gates in both the parents solutions and any insertion should already have simulators.)
      if (new_gate && new_gate->mutatable())
      {
        if (rand01() < problem_.gateOptProb())
        {
          optimize_and_evaluate();  // Optimize parameter of the inserted gate, exploiting quick evaluation methods.
        }
        else
        {
          // (The next two lines are only necessary if we have just performed a 'tweak gate' move. If we are inserting
          // a truly new gate then it already has a random angle.)
          new_gate->mutate();            // Also resets simulators...
          new_gate->createSimulators();  // ...which must be recreated before evaluation.
          quickly_evaluate();
        }
      }
      else
      {
        quickly_evaluate();  // Quick evaluation
      }
      construct_circuit();   // (We would like to omit this if the circuit error is worse than the cutoff. However...
                             // ...this causes issues with the Store, as it tests solutions for equality and this...
                             // ...requires the circuit to be constructed.
      if (error_[0] < problem_.badErrorCutoff())
      {
        // (Not simplifying bad circuits runs the risk that unsimplified circuits reach the adult population (and the
        // Store). Then the section of a parent solution chosen for crossover might simplify to nothing, which isn't
        // really crossover. However, we hope to gain speedups that make such concerns irrelevant.
        circuit_->simplify();
        circuit_->createGateSimulators();  // Create GateSimulators for any SingleTargetGate that hasn't got one.
      }
      evaluate_cost();
    }
  }


  int Solution::length() const
  {
    if (circuit_)
    {
      return circuit_->length();
    }

    assert(left_parent && right_parent);
    int l = right_parent->length() - right_cut + left_cut;
    if (new_gate)
    {
      ++l;
    }
    if (insertion_num >= 0)
    {
      l += problem_.insertion(insertion_num).length();
    }

    return l;
  }


  bool Solution::bad() const
  {
    return bad_;
  }


  double Solution::objective(int objNum) const
  {
    assert(evaluated_);
    assert(0 <= objNum && objNum < problem_.numObj());
    if (objNum == 0)
    {
      return total_cost;  // Casts total_cost to double.
    }
    return error_[objNum - 1] + 0.00001 * total_cost;  // Prevents search heading towards crazy large solutions, just...
  }                                                    // ...because their error is 1e-7 rather than 2e-7! Bit of a...
                                                       // ...cludge. Affects order of Store output.

  long Solution::cost() const
  {
    return total_cost;
  }


  void Solution::output(ostream& out) const
  {
    if (circuit_)
    {
      out << "Circuit:" << endl << *circuit_;
    }
    else
    {
      if (left_circuit == right_circuit)
      {
        out << "Parent circuit:" << endl << *left_circuit << endl;
      }
      else
      {
        out << "Left parent circuit:" << endl << *left_circuit << endl;
        out << "Right parent circuit:" << endl << *right_circuit << endl;
      }
      out << "Cut points are at " << left_cut << " and " << right_cut << "." << endl;
      if (new_gate)
      {
        out << "New gate: " << *new_gate << endl;
      }
      if (insertion_num >= 0)
      {
        out << "Insertion:" << endl << problem_.insertion(insertion_num) << endl;
      }
    }

    if (evaluated_)
    {
      out << "Error: " << error_[0];
      for (auto i = 1; i < problem_.numObj() - 1; ++i)
      {
        out << ", " << error_[i];
      }
      out << endl;
      out << "Gate cost: " << total_cost << endl;
    }
    else
    {
      out << "Solution unevaluated." << endl;
    }
  }


  void Solution::clone_from(const Solution& parent)
  {
    // We don't actually copy the circuit from the parent at all, here! Instead, we just give the child solution a
    // pointer to the parent, which will be the source of the genetic material for the left hand part of the circuit, if
    // crossover is performed. This will be used later to efficiently evaluate the circuit and construct it if
    // necessary.
    left_parent = &parent;
    left_circuit = left_parent->circuit_;
    right_parent = nullptr;
    right_circuit = nullptr;
    left_cut = left_circuit->length();
    right_cut = 0;
  }


  void Solution::evaluate_cost()
  {
    // The easy bit of evaluating a circuit - adding up the cost of the gates.

    if (circuit_)
    {
      total_cost = 0;
      for (auto i = 0; i < circuit_->length(); ++i)
      {
        total_cost += circuit_->gate(i).cost();
      }
    }
    else
    {
      // We may need to evaluate the cost of an 'unconstructed' circuit, if the error exceeds the bad_error_cutoff.
      total_cost = 0;

      // Cost of gates from left parent.
      for (int i = 0; i < left_cut; ++i)
      {
        total_cost += left_circuit->gate(i).cost();
      }

      // Cost of insertion. (Since the circuit has not been constructed, any insertion should still be there.
      if (new_gate)
      {
        total_cost += new_gate->cost();
      }
      if (insertion_num >= 0)
      {
        const auto& insertion = problem_.insertion(insertion_num);
        for (int i = 0; i < insertion.length(); ++i)
        {
          total_cost += insertion.gate(i).cost();
        }
      }

      // Cost of gates from right parent.
      for (int i = right_cut; i < right_circuit->length(); ++i)
      {
        total_cost += right_circuit->gate(i).cost();
      }
    }
  }


  void Solution::assign_worst_cost()
  {
    total_cost = std::numeric_limits<long>::max();
  }


  void crossover(Solution& lhs, Solution& rhs)
  {
    // This merely notifies each solution of its other parent. Deciding on where the crossover takes place happens in
    // the mutation function.
    lhs.right_parent = rhs.left_parent;
    lhs.right_circuit = rhs.left_circuit;
    rhs.right_parent = lhs.left_parent;
    rhs.right_circuit = lhs.left_circuit;
  }


  ostream& operator<<(ostream& out, const Solution& solution)
  {
    solution.output(out);
    return out;
  }


}
